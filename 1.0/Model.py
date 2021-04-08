import os
import os.path as op
import pandas as pd
import json
import numpy as np
import scipy as sp
import scipy.stats as sps
import functools
import gc
import matplotlib.pyplot as plt
from enum import Enum

from Material import Element, Mixture
from Cone import *
from Parameters import *
from cosmicflux import getIntegratedFlux


class ModelState(Enum):
	INITIAL = 0
	WARMUP = 1
	MAINRUN = 2


class Model:
	def __init__(self, modeldir):
		self.datadir = op.join(modeldir, 'Data')
		self.matdir = op.join(modeldir, 'Materials')
		self.accepted = 0
		self.steps = 0
		self.parameters = 0
		self.detectors = {}
		self.materials = {}

		self.nwarmupsteps = 0
		self.nmainsteps = 0

		self.parameter_names = []
		self.jucova = np.empty(())

		self.cones = {}
		self.logprob = 0

		self.conlprobs = pd.Series()

		self.verbosity = 0

		self.initialized = False

		self.parameters = GlobalParameterGroup()
		self.cs_error_parameters = LocalParameterGroup() # Cross-section errors
		self.parameters.declare_child(self.cs_error_parameters)

		self.state = ModelState.INITIAL

		self.initialise_cs_errors()

		self.read_materials()
		self.read_data()

	def __del__(self):
		print('Model destroyed')

	def initialize(self):
		if not self.initialized:
			self.calculate_synthetic_data(self.cs_error_parameters.parameter_df.iloc[-1])
			self.logprob = self.calculate_model_logprob()

			self.parameters.set_jucova_scaling()
			self.initialized = True

	def get_detectors(self):
		return self.detectors

	def read_materials(self):
		matlist = os.listdir(self.matdir)

		for mat in matlist:
			material_directory = op.join(self.matdir, mat)
			infoloc = op.join(material_directory, 'info.json')
			compoloc = op.join(material_directory, 'composition.json')
			denloc = op.join(material_directory, 'density.json')

			with open(infoloc) as file:
				info = json.load(file)

			with open(compoloc) as file:
				compo = json.load(file)

			with open(denloc) as file:
				density = json.load(file)

			if info['composition'] == 'Groom element':
				thisele = Element(mat, compo, self, matpath=material_directory,
								  denmode=info['density'], den_dict=density)
				self.materials.update({mat: thisele})
			elif info['composition'] == 'Groom mixture':
				thismix = Mixture(mat, compo, self, matpath=material_directory,
				                  denmode=info['density'], den_dict=density)
				self.materials.update({mat: thismix})
			elif info['composition'] == 'Oxides':
				pass

	def read_data(self):
		data_files = os.listdir(self.datadir)
		data_files.remove('data.json')

		dataloc = op.join(self.datadir, 'data.json')

		for det in data_files:
			detpath = op.join(self.datadir, det)
			detdf = pd.read_json(detpath)
			self.detectors[detdf.at['ID', 'value']] = detdf

		self.data = pd.read_json(dataloc)
		self.data['Detector'] = self.data['ID'].str.split('_').apply(lambda row: row[0])

		# Organise data in cones
		for index, row in self.data.iterrows():
			thiscone = Cone(self, row['ID'], row)
			# self.cones.append(thiscone)
			self.cones[row['ID']] = thiscone

	def calculate_synthetic_data(self, cs_err_df):
		for name, con in self.cones.items():
			con.calculate_flux(cs_err_df)

	def calculate_model_logprob(self):
		lprob = 0
		for matkey in self.materials:
			matlprob = self.materials[matkey].calculate_logprob()
			lprob += matlprob

		for name, con in self.cones.items():
			conelprob = con.calculate_logprob()
			lprob += conelprob

		errlprob = self.cs_error_parameters.calculate_logprob()
		lprob += errlprob

		return lprob

	def get_material(self, matname):
		return self.materials[matname]

	def propose_model(self):
		# Change parameter dataframes and parameter objects
		self.parameters.propose_model()

		for key in self.materials:
			self.materials[key].update()

		for name, con in self.cones.items():
			con.update()

		# print(self.parameters.parameter_df)
		# for name, con in self.cones.items():
		# 	for key in con.parameters.parameters:
		# 		print(con.parameters.parameters[key].value)
		#
		# for mat in self.materials:
		# 	for key in self.materials[mat].parameters.parameters:
		# 		print(self.materials[mat].parameters.parameters[key].value)
			# print(self.materials[mat].parameters.parameter_df)

	def MH_step(self):
		if not self.initialized:
			self.initialize()
			self.initialized = True

		if self.state == ModelState.WARMUP:
			self.nwarmupsteps += 1
		elif self.state == ModelState.MAINRUN:
			self.nmainsteps += 1

		self.steps += 1
		# Propose a new model
		self.propose_model()
		if self.verbosity > 0:
			print('---------------------------')
			print('Old/Proposed model:')
			print()
			print(self.parameters.parameter_df.tail(2))
			print()

		# Check if parameters are within bounds (otherwise a calculation is not necessary)
		oob = self.parameters.check_oob()

		if not oob:
			# Perform forward calculation
			self.calculate_synthetic_data(self.cs_error_parameters.parameter_df.iloc[-1])

			# Calculate new log probability
			lprob = self.calculate_model_logprob()

			if self.verbosity > 0:
				print('old: ', self.logprob, ' || new: ', lprob)

			# Calculate transition probability
			transprob = np.exp(lprob - self.logprob)

			if transprob > 1: # Accept upward step
				if self.verbosity > 0:
					print('-- accepted (up) --')
				self.logprob = lprob
				self.accepted += 1
			else:
				alpha = sps.uniform.rvs()
				# print('alpha: ', alpha)
				if alpha < transprob: # Accept downward step
					if self.verbosity > 0:
						print('-- accepted (down) --')
						print('alpha: ', alpha)
					self.logprob = lprob
					self.accepted += 1
				else: # Reject downward step
					if self.verbosity > 0:
						print('-- rejected --')
					self.parameters.parameter_df.iloc[-1, :] = self.parameters.parameter_df.iloc[-2, :]

			if self.verbosity > 0:
				print('Log-Probability: ', lprob)
				print('Transition probability: ', transprob)

		else: # Outright rejection
			if self.verbosity > 0:
				print('-- rejected (oob) --')
			self.parameters.parameter_df.iloc[-1, :] = self.parameters.parameter_df.iloc[-2, :]

	def show_models(self):
		print(self.parameters.parameter_df)

	def set_verbosity(self, level):
		self.verbosity = level

	def adapt_proposal_dist(self, n_sims=None):
		self.parameters.adapt_cov_matrix(n_sims)

	def acc_prob(self):
		return float(self.accepted)/float(self.steps)

	def plot(self):
		axes = self.parameters.parameter_df.plot(subplots=True)
		axes[1].set_yscale('log')
		axes[1].set_xlabel('# of steps')
		plt.show()

	def get_parameternames(self):
		return list(self.parameters.parameter_df.columns)

	def get_cone_length(self, cone):
		conename = '_'.join(cone.split('_')[:-1])
		conelength = self.data.loc[self.data['ID'] == conename]['d_topo'].values[0]

		return conelength

	def get_data_info(self):
		return self.data

	def get_data(self, parameter_name=None):
		if parameter_name:
			return self.parameters.parameter_df[parameter_name]
		else:
			return self.parameters.parameter_df

	def get_mr_data(self, parameter_name=None):
		if self.state == ModelState.MAINRUN:
			if parameter_name:
				return self.parameters.parameter_df[parameter_name].iloc[-self.nmainsteps:]
			else:
				return self.parameters.parameter_df.iloc[-self.nmainsteps:]
		else:
			return None

	def save_to_file(self, savepath):
		self.parameters.save_to_file(savepath)

	def show_parametrs(self):
		self.parameters.show()
		self.parameters.show_children()

	def import_simulations(self, filepath, resauxdf):
		import_df = pd.read_json(filepath)

		import_header = import_df.columns
		existing_header = self.parameters.parameter_df.columns

		difference = import_header.difference(existing_header)

		if difference.empty:
			self.parameters.replace_df(import_df)
			self.parameters.complete_sync_to_children()

			self.nwarmupsteps = resauxdf['n_warmup_sim']
			self.nmainsteps = resauxdf['n_main_sim']

		if self.nwarmupsteps > 0:
			if self.nmainsteps > 0:
				self.change_state(ModelState.MAINRUN)
			else:
				self.change_state(ModelState.WARMUP)
		else:
			self.change_state(ModelState.INITIAL)

		# print(self.nwarmupsteps)
		# print(self.nmainsteps)

	def change_state(self, newModelState):
		self.state = newModelState

	def initialise_cs_errors(self):
		ion_err = Parameter(chr(0x03C3) + '_ion', 'cs_error')
		ion_err.set_pdf(sps.lognorm, loc=0, scale=1, s=np.log(1.06))
		ion_err.set_value(sps.lognorm.rvs(loc=0, scale=1, s=np.log(1.06)))
		ion_err.set_inistd(np.log(1.06))

		self.cs_error_parameters.add_Parameter(ion_err)
		self.parameters.add_Parameter(ion_err)

		bre_err = Parameter(chr(0x03C3) + '_bre', 'cs_error')
		bre_err.set_pdf(sps.lognorm, loc=0, scale=1, s=np.log(1.01))
		bre_err.set_value(sps.lognorm.rvs(loc=0, scale=1, s=np.log(1.01)))
		bre_err.set_inistd(np.log(1.01))

		self.cs_error_parameters.add_Parameter(bre_err)
		self.parameters.add_Parameter(bre_err)

		pai_err = Parameter(chr(0x03C3) + '_pai', 'cs_error')
		pai_err.set_pdf(sps.lognorm, loc=0, scale=1, s=np.log(1.05))
		pai_err.set_value(sps.lognorm.rvs(loc=0, scale=1, s=np.log(1.05)))
		pai_err.set_inistd(np.log(1.05))

		self.cs_error_parameters.add_Parameter(pai_err)
		self.parameters.add_Parameter(pai_err)

		pho_err = Parameter(chr(0x03C3) + '_pho', 'cs_error')
		pho_err.set_pdf(sps.lognorm, loc=0, scale=1, s=np.log(1.3))
		pho_err.set_value(sps.lognorm.rvs(loc=0, scale=1, s=np.log(1.3)))
		pho_err.set_inistd(np.log(1.3))

		self.cs_error_parameters.add_Parameter(pho_err)
		self.parameters.add_Parameter(pho_err)


class SingleConeModelCluster:
	def __init__(self, modeldir):
		self.datadir = op.join(modeldir, 'Data')
		self.matdir = op.join(modeldir, 'Materials')

		self.datafile = op.join(self.datadir, 'data.json')

		self.models = {}

		self.state = ModelState.INITIAL

		self.nwarmupsteps=0
		self.nmainsteps=0

		self.namemap = {}

		# Get all cones
		data = pd.read_json(self.datafile)

		ids = [row['ID'] for index, row in data.iterrows()]

		self.parameters = UniversalParameterGroup()

		for id in ids:
			self.models[id] = SingleConeModel(modeldir, id)
			self.parameters.declare_child(id, self.models[id].parameters)

	def initialize(self):
		for id, model in self.models.items():
			model.initialize()

	def get_detectors(self):
		firstmodel = list(self.models.values())[0]
		return firstmodel.detectors

	def get_data_info(self):
		firstmodel = list(self.models.values())[0]
		return firstmodel.data

	def import_simulations(self, prefix, outdir, resauxdf):
		for id, model in self.models.items():
			file = os.path.join(outdir, prefix + '_' + id + '_results.json')
			model.import_simulations(file, resauxdf)

		firstmodel = list(self.models.values())[0]

		self.nwarmupsteps = firstmodel.nwarmupsteps
		self.nmainsteps = firstmodel.nmainsteps
		self.state = firstmodel.state

		self.parameters.sync_from_children()

		self.parameters.show()

	# def calculate_synthetic_data(self, cs_err_df):
	# 	pass
	#
	# def calculate_model_logprob(self):
	# 	pass
	#
	# def get_material(self, matname):
	# 	pass
	#
	# def propose_model(self):
	# 	pass

	def MH_step(self, nSimulations, popup, mode):
		# Simulation
		icone = 0
		for id, cone in self.models.items():
			icone += 1
			popup.nametowidget('ptcone').config(text='Cone %s / %s' % (icone, len(self.models)))
			for i in range (1, nSimulations + 1):
				popup.nametowidget('ptsim').config(text='Simulating %s / %s' % (i, nSimulations))
				popup.nametowidget('pbsim').step()
				popup.update()
				cone.MH_step()
			popup.nametowidget('pbcone').step()
			popup.update()

		# Update simulation numbers
		if mode == 'WU':
			self.nwarmupsteps += nSimulations
		elif mode == 'MR':
			self.nmainsteps += nSimulations

		# Synchronise local cone parameters with global parameters
		self.parameters.sync_from_children()

		self.parameters.show()

	# def show_models(self):
	# 	pass

	def set_verbosity(self, level):
		for model in self.models.values():
			model.set_verbosity(level)

	def adapt_proposal_dist(self, n_sims=None):
		for model in self.models.values():
			model.adapt_proposal_dist(n_sims)

	# def acc_prob(self):
	# 	pass
	#
	# def plot(self):
	# 	pass

	def get_parameternames(self):
		return list(self.parameters.parameter_df.columns)

	def get_cone_length(self, cone):
		return list(self.models.values())[0].get_cone_length(cone)


	def get_data(self, parameter_name=None):
		if parameter_name:
			return self.parameters.parameter_df[parameter_name]
		else:
			return self.parameters.parameter_df

	def get_mr_data(self, parameter_name=None):
		if self.state == ModelState.MAINRUN:
			if parameter_name:
				return self.parameters.parameter_df[parameter_name].iloc[-self.nmainsteps:]
			else:
				return self.parameters.parameter_df.iloc[-self.nmainsteps:]
		else:
			return None

	def save_to_file(self, savepath):
		savedir = os.path.dirname(savepath)

		chainnr = os.path.basename(savepath).split(sep='_')[0]

		for id, model in self.models.items():
			model.save_to_file(os.path.join(savedir, chainnr + '_' + id + '_results.json'))

	def change_state(self, newModelState):
		for id, model in self.models.items():
			model.change_state(newModelState)

	# def initialise_cs_errors(self):
	# 	pass

	def get_nr_of_cones(self):
		return len(self.models)


class SingleConeModel:
	def __init__(self, modeldir, coneID):
		self.datadir = op.join(modeldir, 'Data')
		self.matdir = op.join(modeldir, 'Materials')
		self.accepted = 0
		self.steps = 0
		self.parameters = 0
		self.detectors = {}
		self.materials = {}

		self.nwarmupsteps = 0
		self.nmainsteps = 0

		self.parameter_names = []
		self.jucova = np.empty(())

		self.coneID = coneID

		self.cones = {}
		self.logprob = 0

		self.conlprobs = pd.Series()

		self.verbosity = 0

		self.initialized = False

		self.parameters = GlobalParameterGroup()
		self.cs_error_parameters = LocalParameterGroup() # Cross-section errors
		self.parameters.declare_child(self.cs_error_parameters)

		self.state = ModelState.INITIAL

		self.initialise_cs_errors()

		self.read_materials()
		self.read_data()

	def __del__(self):
		print('Model destroyed')

	def initialize(self):
		if not self.initialized:
			self.calculate_synthetic_data(self.cs_error_parameters.parameter_df.iloc[-1])
			self.logprob = self.calculate_model_logprob()

			self.parameters.set_jucova_scaling()
			self.initialized = True

	def read_materials(self):
		matlist = os.listdir(self.matdir)

		for mat in matlist:
			material_directory = op.join(self.matdir, mat)
			infoloc = op.join(material_directory, 'info.json')
			compoloc = op.join(material_directory, 'composition.json')
			denloc = op.join(material_directory, 'density.json')

			with open(infoloc) as file:
				info = json.load(file)

			with open(compoloc) as file:
				compo = json.load(file)

			with open(denloc) as file:
				density = json.load(file)

			if info['composition'] == 'Groom element':
				thisele = Element(mat, compo, self, matpath=material_directory,
								  denmode=info['density'], den_dict=density)
				self.materials.update({mat: thisele})
			elif info['composition'] == 'Groom mixture':
				thismix = Mixture(mat, compo, self, matpath=material_directory,
				                  denmode=info['density'], den_dict=density)
				self.materials.update({mat: thismix})
			elif info['composition'] == 'Oxides':
				pass

	def read_data(self):
		data_files = os.listdir(self.datadir)
		data_files.remove('data.json')

		dataloc = op.join(self.datadir, 'data.json')

		for det in data_files:
			detpath = op.join(self.datadir, det)
			detdf = pd.read_json(detpath)
			self.detectors[detdf.at['ID', 'value']] = detdf

		self.data = pd.read_json(dataloc)
		self.data['Detector'] = self.data['ID'].str.split('_').apply(lambda row: row[0])

		# Organise data in cones
		row = self.data.loc[self.data['ID'] == self.coneID].squeeze()

		thiscone = Cone(self, row['ID'], row)
		self.cones[row['ID']] = thiscone

	def calculate_synthetic_data(self, cs_err_df):
		for name, con in self.cones.items():
			con.calculate_flux(cs_err_df)

	def calculate_model_logprob(self):
		lprob = 0
		for matkey in self.materials:
			matlprob = self.materials[matkey].calculate_logprob()
			lprob += matlprob

		for name, con in self.cones.items():
			conelprob = con.calculate_logprob()
			lprob += conelprob

		errlprob = self.cs_error_parameters.calculate_logprob()
		lprob += errlprob

		return lprob

	def get_material(self, matname):
		return self.materials[matname]

	def propose_model(self):
		# Change parameter dataframes and parameter objects
		self.parameters.propose_model()

		for key in self.materials:
			self.materials[key].update()

		for name, con in self.cones.items():
			con.update()

	def MH_step(self):
		if not self.initialized:
			self.initialize()
			self.initialized = True

		if self.state == ModelState.WARMUP:
			self.nwarmupsteps += 1
		elif self.state == ModelState.MAINRUN:
			self.nmainsteps += 1

		self.steps += 1
		# Propose a new model
		self.propose_model()
		if self.verbosity > 0:
			print('---------------------------')
			print('Old/Proposed model:')
			print()
			print(self.parameters.parameter_df.tail(2))
			print()

		# Check if parameters are within bounds (otherwise a calculation is not necessary)
		oob = self.parameters.check_oob()

		if not oob:
			# Perform forward calculation
			self.calculate_synthetic_data(self.cs_error_parameters.parameter_df.iloc[-1])

			# Calculate new log probability
			lprob = self.calculate_model_logprob()

			if self.verbosity > 0:
				print('old: ', self.logprob, ' || new: ', lprob)

			# Calculate transition probability
			transprob = np.exp(lprob - self.logprob)

			if transprob > 1: # Accept upward step
				if self.verbosity > 0:
					print('-- accepted (up) --')
				self.logprob = lprob
				self.accepted += 1
			else:
				alpha = sps.uniform.rvs()
				# print('alpha: ', alpha)
				if alpha < transprob: # Accept downward step
					if self.verbosity > 0:
						print('-- accepted (down) --')
						print('alpha: ', alpha)
					self.logprob = lprob
					self.accepted += 1
				else: # Reject downward step
					if self.verbosity > 0:
						print('-- rejected --')
					self.parameters.parameter_df.iloc[-1, :] = self.parameters.parameter_df.iloc[-2, :]

			if self.verbosity > 0:
				print('Log-Probability: ', lprob)
				print('Transition probability: ', transprob)

		else: # Outright rejection
			if self.verbosity > 0:
				print('-- rejected (oob) --')
			self.parameters.parameter_df.iloc[-1, :] = self.parameters.parameter_df.iloc[-2, :]

	def show_models(self):
		print(self.parameters.parameter_df)

	def set_verbosity(self, level):
		self.verbosity = level

	def adapt_proposal_dist(self, n_sims=None):
		self.parameters.adapt_cov_matrix(n_sims)

	def acc_prob(self):
		return float(self.accepted)/float(self.steps)

	def plot(self):
		axes = self.parameters.parameter_df.plot(subplots=True)
		axes[1].set_yscale('log')
		axes[1].set_xlabel('# of steps')
		plt.show()

	def get_parameternames(self):
		return list(self.parameters.parameter_df.columns)

	def get_cone_length(self, cone):
		conename = '_'.join(cone.split('_')[:-1])
		conelength = self.data.loc[self.data['ID'] == conename]['d_topo'].values[0]

		return conelength

	def get_data(self, parameter_name=None):
		if parameter_name:
			return self.parameters.parameter_df[parameter_name]
		else:
			return self.parameters.parameter_df

	def get_mr_data(self, parameter_name=None):
		if self.state == ModelState.MAINRUN:
			if parameter_name:
				return self.parameters.parameter_df[parameter_name].iloc[-self.nmainsteps:]
			else:
				return self.parameters.parameter_df.iloc[-self.nmainsteps:]
		else:
			return None

	def save_to_file(self, savepath):
		self.parameters.save_to_file(savepath)

	def show_parametrs(self):
		self.parameters.show()
		self.parameters.show_children()

	def import_simulations(self, filepath, resauxdf):
		import_df = pd.read_json(filepath)

		import_header = import_df.columns
		existing_header = self.parameters.parameter_df.columns

		difference = import_header.difference(existing_header)

		if difference.empty:
			self.parameters.replace_df(import_df)
			self.parameters.complete_sync_to_children()

			self.nwarmupsteps = resauxdf['n_warmup_sim']
			self.nmainsteps = resauxdf['n_main_sim']

		if self.nwarmupsteps > 0:
			if self.nmainsteps > 0:
				self.change_state(ModelState.MAINRUN)
			else:
				self.change_state(ModelState.WARMUP)
		else:
			self.change_state(ModelState.INITIAL)

		# print(self.nwarmupsteps)
		# print(self.nmainsteps)

	def change_state(self, newModelState):
		self.state = newModelState

	def initialise_cs_errors(self):
		ion_err = Parameter(chr(0x03C3) + '_ion', 'cs_error')
		ion_err.set_pdf(sps.lognorm, loc=0, scale=1, s=np.log(1.06))
		ion_err.set_value(sps.lognorm.rvs(loc=0, scale=1, s=np.log(1.06)))
		ion_err.set_inistd(np.log(1.06))

		self.cs_error_parameters.add_Parameter(ion_err)
		self.parameters.add_Parameter(ion_err)

		bre_err = Parameter(chr(0x03C3) + '_bre', 'cs_error')
		bre_err.set_pdf(sps.lognorm, loc=0, scale=1, s=np.log(1.01))
		bre_err.set_value(sps.lognorm.rvs(loc=0, scale=1, s=np.log(1.01)))
		bre_err.set_inistd(np.log(1.01))

		self.cs_error_parameters.add_Parameter(bre_err)
		self.parameters.add_Parameter(bre_err)

		pai_err = Parameter(chr(0x03C3) + '_pai', 'cs_error')
		pai_err.set_pdf(sps.lognorm, loc=0, scale=1, s=np.log(1.05))
		pai_err.set_value(sps.lognorm.rvs(loc=0, scale=1, s=np.log(1.05)))
		pai_err.set_inistd(np.log(1.05))

		self.cs_error_parameters.add_Parameter(pai_err)
		self.parameters.add_Parameter(pai_err)

		pho_err = Parameter(chr(0x03C3) + '_pho', 'cs_error')
		pho_err.set_pdf(sps.lognorm, loc=0, scale=1, s=np.log(1.3))
		pho_err.set_value(sps.lognorm.rvs(loc=0, scale=1, s=np.log(1.3)))
		pho_err.set_inistd(np.log(1.3))

		self.cs_error_parameters.add_Parameter(pho_err)
		self.parameters.add_Parameter(pho_err)

# class SingleConeModel:
# 	def __init__(self, modeldir, coneID):
# 		self.datadir = op.join(modeldir, 'Data')
# 		self.matdir = op.join(modeldir, 'Materials')
# 		self.accepted = 0
# 		self.steps = 0
# 		self.parameters = 0
# 		self.detectors = {}
# 		self.materials = {}
#
# 		self.coneID = coneID
#
# 		self.cones = {}
#
# 		self.nwarmupsteps = 0
# 		self.nmainsteps = 0
#
# 		self.parameter_names = []
# 		self.jucova = np.empty(())
#
# 		self.logprob = 0
#
# 		self.conlprobs = pd.Series()
#
# 		self.verbosity = 0
#
# 		self.initialized = False
# 		#
# 		self.parameters = GlobalParameterGroup()
# 		self.cs_error_parameters = LocalParameterGroup() # Cross-section errors
# 		self.parameters.declare_child(self.cs_error_parameters)
#
# 		self.state = ModelState.INITIAL
#
# 		self.initialise_cs_errors()
# 		self.read_materials()
# 		self.read_data()
#
# 	def __del__(self):
# 		print('Model destroyed')
#
# 	def initialize(self):
# 		if not self.initialized:
# 			self.calculate_synthetic_data(self.cs_error_parameters.parameter_df.iloc[-1])
# 			self.logprob = self.calculate_model_logprob()
#
# 			self.parameters.set_jucova_scaling()
# 			self.initialized = True
#
# 	def read_materials(self):
# 		matlist = os.listdir(self.matdir)
#
# 		print(matlist)
#
# 		for mat in matlist:
# 			material_directory = op.join(self.matdir, mat)
# 			infoloc = op.join(material_directory, 'info.json')
# 			compoloc = op.join(material_directory, 'composition.json')
# 			denloc = op.join(material_directory, 'density.json')
#
# 			with open(infoloc) as file:
# 				info = json.load(file)
#
# 			with open(compoloc) as file:
# 				compo = json.load(file)
#
# 			with open(denloc) as file:
# 				density = json.load(file)
#
# 			name = self.coneID + '_' + mat
#
# 			if info['composition'] == 'Groom element':
# 				thisele = Element(name, compo, self, matpath=material_directory,
# 								  denmode=info['density'], den_dict=density)
# 				self.materials.update({name: thisele})
# 			elif info['composition'] == 'Groom mixture':
# 				thismix = Mixture(name, compo, self, matpath=material_directory,
# 				                  denmode=info['density'], den_dict=density)
# 				self.materials.update({name: thismix})
# 			elif info['composition'] == 'Oxides':
# 				pass
#
# 	def read_data(self):
# 		data_files = os.listdir(self.datadir)
# 		data_files.remove('data.json')
#
# 		dataloc = op.join(self.datadir, 'data.json')
#
# 		for det in data_files:
# 			detpath = op.join(self.datadir, det)
# 			detdf = pd.read_json(detpath)
# 			self.detectors[detdf.at['ID', 'value']] = detdf
#
# 		self.data = pd.read_json(dataloc)
# 		self.data['Detector'] = self.data['ID'].str.split('_').apply(lambda row: row[0])
#
# 		# Organise data in cones
# 		row = self.data.loc[self.data['ID'] == self.coneID].squeeze()
#
# 		thiscone = Cone(self, row['ID'], row)
# 		self.cones[row['ID']] = thiscone
#
#
# 	def calculate_synthetic_data(self, cs_err_df):
# 		for name, con in self.cones.items():
# 			con.calculate_flux(cs_err_df)
#
# 	def calculate_model_logprob(self):
# 		lprob = 0
# 		for matkey in self.materials:
# 			matlprob = self.materials[matkey].calculate_logprob()
# 			lprob += matlprob
#
# 		for name, con in self.cones.items():
# 			conelprob = con.calculate_logprob()
# 			lprob += conelprob
#
# 		errlprob = self.cs_error_parameters.calculate_logprob()
# 		lprob += errlprob
#
# 		return lprob
#
# 	def get_material(self, matname):
# 		name = self.coneID + '_' + matname
# 		return self.materials[name]
#
# 	def propose_model(self):
# 		# Change parameter dataframes and parameter objects
# 		self.parameters.propose_model()
#
# 		for key in self.materials:
# 			self.materials[key].update()
#
# 		for name, con in self.cones.items():
# 			con.update()
#
# 		# print(self.parameters.parameter_df)
# 		# for name, con in self.cones.items():
# 		# 	for key in con.parameters.parameters:
# 		# 		print(con.parameters.parameters[key].value)
# 		#
# 		# for mat in self.materials:
# 		# 	for key in self.materials[mat].parameters.parameters:
# 		# 		print(self.materials[mat].parameters.parameters[key].value)
# 			# print(self.materials[mat].parameters.parameter_df)
#
# 	def MH_step(self):
# 		if not self.initialized:
# 			self.initialize()
# 			self.initialized = True
#
# 		if self.state == ModelState.WARMUP:
# 			self.nwarmupsteps += 1
# 		elif self.state == ModelState.MAINRUN:
# 			self.nmainsteps += 1
#
# 		self.steps += 1
# 		# Propose a new model
# 		self.propose_model()
# 		if self.verbosity > 0:
# 			print('---------------------------')
# 			print('Old/Proposed model:')
# 			print()
# 			print(self.parameters.parameter_df.tail(2))
# 			print()
#
# 		# Check if parameters are within bounds (otherwise a calculation is not necessary)
# 		oob = self.parameters.check_oob()
#
# 		if not oob:
# 			# Perform forward calculation
# 			self.calculate_synthetic_data(self.cs_error_parameters.parameter_df.iloc[-1])
#
# 			# Calculate new log probability
# 			lprob = self.calculate_model_logprob()
#
# 			if self.verbosity > 0:
# 				print('old: ', self.logprob, ' || new: ', lprob)
#
# 			# Calculate transition probability
# 			transprob = np.exp(lprob - self.logprob)
#
# 			if transprob > 1: # Accept upward step
# 				if self.verbosity > 0:
# 					print('-- accepted (up) --')
# 				self.logprob = lprob
# 				self.accepted += 1
# 			else:
# 				alpha = sps.uniform.rvs()
# 				# print('alpha: ', alpha)
# 				if alpha < transprob: # Accept downward step
# 					if self.verbosity > 0:
# 						print('-- accepted (down) --')
# 						print('alpha: ', alpha)
# 					self.logprob = lprob
# 					self.accepted += 1
# 				else: # Reject downward step
# 					if self.verbosity > 0:
# 						print('-- rejected --')
# 					self.parameters.parameter_df.iloc[-1, :] = self.parameters.parameter_df.iloc[-2, :]
#
# 			if self.verbosity > 0:
# 				print('Log-Probability: ', lprob)
# 				print('Transition probability: ', transprob)
#
# 		else: # Outright rejection
# 			if self.verbosity > 0:
# 				print('-- rejected (oob) --')
# 			self.parameters.parameter_df.iloc[-1, :] = self.parameters.parameter_df.iloc[-2, :]
#
# 	def show_models(self):
# 		print(self.parameters.parameter_df)
#
# 	def set_verbosity(self, level):
# 		self.verbosity = level
#
# 	def adapt_proposal_dist(self, n_sims=None):
# 		self.parameters.adapt_cov_matrix(n_sims)
#
# 	def acc_prob(self):
# 		return float(self.accepted)/float(self.steps)
#
# 	def plot(self):
# 		axes = self.parameters.parameter_df.plot(subplots=True)
# 		axes[1].set_yscale('log')
# 		axes[1].set_xlabel('# of steps')
# 		plt.show()
#
# 	def get_parameternames(self):
# 		return list(self.parameters.parameter_df.columns)
#
# 	def get_cone_length(self, cone):
# 		conename = '_'.join(cone.split('_')[:-1])
# 		conelength = self.data.loc[self.data['ID'] == conename]['d_topo'].values[0]
#
# 		return conelength
#
# 	def get_data(self, parameter_name=None):
# 		if parameter_name:
# 			return self.parameters.parameter_df[parameter_name]
# 		else:
# 			return self.parameters.parameter_df
#
# 	def get_mr_data(self, parameter_name=None):
# 		if self.state == ModelState.MAINRUN:
# 			if parameter_name:
# 				return self.parameters.parameter_df[parameter_name].iloc[-self.nmainsteps:]
# 			else:
# 				return self.parameters.parameter_df.iloc[-self.nmainsteps:]
# 		else:
# 			return None
#
# 	def save_to_file(self, savepath):
# 		self.parameters.save_to_file(savepath)
#
# 	def show_parametrs(self):
# 		self.parameters.show()
# 		self.parameters.show_children()
#
# 	def import_simulations(self, filepath, resauxdf):
# 		import_df = pd.read_json(filepath)
#
# 		import_header = import_df.columns
# 		existing_header = self.parameters.parameter_df.columns
#
# 		difference = import_header.difference(existing_header)
#
# 		if difference.empty:
# 			self.parameters.replace_df(import_df)
# 			self.parameters.complete_sync_to_children()
#
# 			self.nwarmupsteps = resauxdf['n_warmup_sim']
# 			self.nmainsteps = resauxdf['n_main_sim']
#
# 		if self.nwarmupsteps > 0:
# 			if self.nmainsteps > 0:
# 				self.change_state(ModelState.MAINRUN)
# 			else:
# 				self.change_state(ModelState.WARMUP)
# 		else:
# 			self.change_state(ModelState.INITIAL)
#
# 		# print(self.nwarmupsteps)
# 		# print(self.nmainsteps)
#
# 	def change_state(self, newModelState):
# 		self.state = newModelState
#
# 	def initialise_cs_errors(self):
# 		name = self.coneID + '_'
#
# 		ion_err = Parameter(name + chr(0x03C3) + '_ion', 'cs_error')
# 		ion_err.set_pdf(sps.lognorm, loc=0, scale=1, s=np.log(1.06))
# 		ion_err.set_value(sps.lognorm.rvs(loc=0, scale=1, s=np.log(1.06)))
# 		ion_err.set_inistd(np.log(1.06))
#
# 		self.cs_error_parameters.add_Parameter(ion_err)
# 		self.parameters.add_Parameter(ion_err)
#
# 		bre_err = Parameter(name + chr(0x03C3) + '_bre', 'cs_error')
# 		bre_err.set_pdf(sps.lognorm, loc=0, scale=1, s=np.log(1.01))
# 		bre_err.set_value(sps.lognorm.rvs(loc=0, scale=1, s=np.log(1.01)))
# 		bre_err.set_inistd(np.log(1.01))
#
# 		self.cs_error_parameters.add_Parameter(bre_err)
# 		self.parameters.add_Parameter(bre_err)
#
# 		pai_err = Parameter(name + chr(0x03C3) + '_pai', 'cs_error')
# 		pai_err.set_pdf(sps.lognorm, loc=0, scale=1, s=np.log(1.05))
# 		pai_err.set_value(sps.lognorm.rvs(loc=0, scale=1, s=np.log(1.05)))
# 		pai_err.set_inistd(np.log(1.05))
#
# 		self.cs_error_parameters.add_Parameter(pai_err)
# 		self.parameters.add_Parameter(pai_err)
#
# 		pho_err = Parameter(name + chr(0x03C3) + '_pho', 'cs_error')
# 		pho_err.set_pdf(sps.lognorm, loc=0, scale=1, s=np.log(1.3))
# 		pho_err.set_value(sps.lognorm.rvs(loc=0, scale=1, s=np.log(1.3)))
# 		pho_err.set_inistd(np.log(1.3))
#
# 		self.cs_error_parameters.add_Parameter(pho_err)
# 		self.parameters.add_Parameter(pho_err)


# Single Cone Batch Inversion Model (-> Sicobi)
# class SicobiModel():
# 	def __init__(self, modeldir):
# 		self.datadir = op.join(modeldir, 'Data')
# 		self.matdir = op.join(modeldir, 'Materials')
#
# 		self.detectors = {}
# 		self.cones = {}
#
# 		# self.accepted = 0
# 		# self.steps = 0
# 		# self.parameters = 0
# 		# self.materials = {}
# 		#
# 		# self.nwarmupsteps = 0
# 		# self.nmainsteps = 0
# 		#
# 		# self.parameter_names = []
# 		# self.jucova = np.empty(())
# 		#
# 		# self.logprob = 0
# 		#
# 		# self.conlprobs = pd.Series()
# 		#
# 		# self.verbosity = 0
# 		#
# 		# self.initialized = False
# 		#
# 		# self.parameters = GlobalParameterGroup()
# 		# self.cs_error_parameters = LocalParameterGroup()  # Cross-section errors
# 		# self.parameters.declare_child(self.cs_error_parameters)
# 		#
# 		# self.state = ModelState.INITIAL
# 		#
# 		# self.initialise_cs_errors()
#
# 		self.read_data()
#
# 	def __del__(self):
# 		print('Model destroyed')
#
# 	def initialize(self):
# 		if not self.initialized:
# 			# self.calculate_synthetic_data(self.cs_error_parameters.parameter_df.iloc[-1])
# 			# self.logprob = self.calculate_model_logprob()
# 			#
# 			# self.parameters.set_jucova_scaling()
# 			self.initialized = True
#
# 	def read_materials(self):
# 		matlist = os.listdir(self.matdir)
#
# 		for mat in matlist:
# 			material_directory = op.join(self.matdir, mat)
# 			infoloc = op.join(material_directory, 'info.json')
# 			compoloc = op.join(material_directory, 'composition.json')
# 			denloc = op.join(material_directory, 'density.json')
#
# 			with open(infoloc) as file:
# 				info = json.load(file)
#
# 			with open(compoloc) as file:
# 				compo = json.load(file)
#
# 			with open(denloc) as file:
# 				density = json.load(file)
#
# 			if info['composition'] == 'Groom element':
# 				thisele = Element(mat, compo, self, matpath=material_directory,
# 				                  denmode=info['density'], den_dict=density)
# 				self.materials.update({mat: thisele})
# 			elif info['composition'] == 'Groom mixture':
# 				thismix = Mixture(mat, compo, self, matpath=material_directory,
# 				                  denmode=info['density'], den_dict=density)
# 				self.materials.update({mat: thismix})
# 			elif info['composition'] == 'Oxides':
# 				pass
#
# 	def read_data(self):
# 		data_files = os.listdir(self.datadir)
# 		data_files.remove('data.json')
#
# 		dataloc = op.join(self.datadir, 'data.json')
#
# 		for det in data_files:
# 			detpath = op.join(self.datadir, det)
# 			detdf = pd.read_json(detpath)
# 			self.detectors[detdf.at['ID', 'value']] = detdf
#
# 		self.data = pd.read_json(dataloc)
# 		self.data['Detector'] = self.data['ID'].str.split('_').apply(lambda row: row[0])
#
# 		# Organise data in cones
# 		for index, row in self.data.iterrows():
# 			thiscone = Cone(self, row['ID'], row)
# 			# self.cones.append(thiscone)
# 			self.cones[row['ID']] = thiscone




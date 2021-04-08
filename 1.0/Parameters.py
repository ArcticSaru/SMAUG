from typing import List, Any

import numpy as np
import scipy.linalg as spl
import scipy.stats as sps
import pandas as pd
from functools import partial


class Parameter:

	def __init__(self, name, par_type):
		self.name = name
		self.type = par_type
		self.value = 0
		self.bounds = []
		self.l_pdf = None
		self.rvs = None
		self.logprob = 0
		self.inistd = 0

	def set_value(self, value):
		if isinstance(value, list):
			if len(value) == 1:
				self.value = value[0]

		else:
			self.value = value

	def set_bounds(self, bmin, bmax):
		self.bounds = [bmin, bmax]

	def set_pdf(self, pdfobject, **kwargs):
		# Set pdf
		self.l_pdf = partial(pdfobject.logpdf, **kwargs)
		self.rvs = partial(pdfobject.rvs, **kwargs)

	# # Set initial value
	# self.set_value(pdfobject.rvs(**kwargs))

	def set_inistd(self, initial_sigma):
		self.inistd = initial_sigma

	def calculate_logprob(self):
		self.logprob = self.l_pdf(self.value)
		# print(self.name, ': ', self.logprob)

		return self.logprob

	def rvs(self):
		return self.rvs()

	def check_oob(self):
		return False


class MultiDimParameter:
	def __init__(self, name, par_type):
		self.name = name
		self.type = par_type
		self.value = 0
		self.bounds = []
		self.l_pdf = None
		self.logprob = 0
		self.inistd = 0
		self.oob = False

	def set_value(self, value):
		self.value = value

	def set_bounds(self, bmin, bmax):
		self.bounds = [bmin, bmax]

	def set_pdf(self, pdfobject, **kwargs):
		# Set pdf
		# print(kwargs)
		self.l_pdf = partial(pdfobject.logpdf, **kwargs)

	def calculate_logprob(self):
		self.logprob = self.l_pdf(self.value)
		#print(self.name, ': ', self.logprob)

		return self.logprob

	def check_oob(self):
		if self.bounds[0] < self.value[0] < self.bounds[1]:
			return False
		else:
			#print(self.name, ' is oob: ', self.value)
			return True


class LocalParameterGroup:
	def __init__(self):
		self.parameter_df = pd.DataFrame()
		self.parameters = {}
		self.jucova = np.empty((0, 0))

	def add_Parameter(self, par):
		if isinstance(par, Parameter):
			self.parameter_df[par.name] = [par.value]
			self.parameters[par.name] = par
			self.jucova = spl.block_diag(self.jucova, par.inistd ** 2)
		elif isinstance(par, MultiDimParameter):
			self.parameters[par.name] = par
			for i in range(0, len(par.value)):
				self.parameter_df[par.name + str(i + 1)] = par.value[i]
			#TODO: scaling factor for thickness jucova...
			parcova = np.diag(np.ones((len(par.value), 1))) * 0.1
			# print(len(par.value))
			# print(parcova)
			#self.jucova = spl.block_diag(self.jucova, parcova)

	def check_oob(self):
		for par in self.parameters.values():
			oob = par.check_oob()
			if oob:
				return True

		return False

	def check_oob_par(self, parameter):
		return self.parameters[parameter].check_oob()


	def calculate_logprob(self):
		lprob = 0
		for parkey in self.parameters:
			lprob += self.parameters[parkey].calculate_logprob()
		return lprob

	def update_parameter(self):
		for key in self.parameters:
			boolidx = self.parameter_df.columns.str.contains(key, regex=False)

			self.parameters[key].set_value(self.parameter_df.iloc[-1, boolidx].values.tolist())

	def contains(self, parameter_name):
		if parameter_name in self.parameters.keys():
			return True
		else:
			return False

	def generate_sample(self, parameter_name):
		if self.contains(parameter_name):
			return self.parameters[parameter_name].rvs()
		else:
			return False


class GlobalParameterGroup:
	def __init__(self):
		self.parameter_df = pd.DataFrame()
		self.children = []

		# self.parameters = {}
		self.jucova = np.empty((0, 0))

	def add_Parameter(self, par):
		if isinstance(par, Parameter):
			self.parameter_df[par.name] = [par.value]
			# self.parameters[par.name] = par
			self.jucova = spl.block_diag(self.jucova, par.inistd ** 2)
		elif isinstance(par, MultiDimParameter):
			# self.parameters[par.name] = par
			for i in range(0, len(par.value)):
				self.parameter_df[par.name + str(i + 1)] = par.value[i]
			#TODO: Set factor as user input
			parcova = np.diag(np.ones((len(par.value), 1)))*0.05
			# print(len(par.value))
			# print(parcova)
			self.jucova = spl.block_diag(self.jucova, parcova)

	def check_oob(self):
		for child in self.children:
			oob = child.check_oob()
			if oob:
				return True

		return False

	def propose_model(self):
		delta_model = sps.multivariate_normal.rvs(cov=self.jucova)
		new_row = pd.Series(self.parameter_df.values[-1] + delta_model, index=self.parameter_df.columns)

		# append to Data Frame
		self.parameter_df = self.parameter_df.append(new_row, ignore_index=True)

		# print(self.jucova)
		# print('update: ', delta_model)
		# print(new_row)

		self.sync_to_children()

	def propose_sicobi_model(self):
		not_simulate = []

		gibbspars = [col for col in self.parameter_df.columns if (chr(0x03C3) in col or chr(0x03C1) in col)]
		lenpars = [col for col in self.parameter_df.columns if '_r' in col]
		# print(gibbspars)

		# print('pars:\n', self.parameter_df)

		delta_model = sps.multivariate_normal.rvs(cov=self.jucova)
		new_row = pd.Series(self.parameter_df.values[-1] + delta_model, index=self.parameter_df.columns)

		# print('new_row:\n', new_row)

		# Modify entries
		for par in gibbspars:
			for child in self.children:
				val = child.generate_sample(par)
				if val:
					new_row.at[par] = val
		# print('new_row:\n', new_row)

		self.parameter_df = self.parameter_df.append(new_row, ignore_index=True)

		self.sync_to_children()

	def declare_child(self, child_obj):
		self.children.append(child_obj)

	def sync_to_children(self):
		# Write global entries to local ones
		for child in self.children:
			common_cols = self.parameter_df.columns.intersection(child.parameter_df.columns)
			child.parameter_df = pd.concat([child.parameter_df, self.parameter_df.tail(1)[common_cols]])

			child.update_parameter()

	def complete_sync_to_children(self):
		# Write global entries to local ones
		for child in self.children:
			common_cols = self.parameter_df.columns.intersection(child.parameter_df.columns)
			child.parameter_df = self.parameter_df[common_cols]

			child.update_parameter()

	def adapt_cov_matrix(self, n_sims=None):
		datlen = len(self.parameter_df.index)
		# print('# simulations used & tot: (%s|%s)' % (n_sims, len(self.parameter_df.index)))

		temp_df = self.parameter_df

		if n_sims.isdigit():
			n_sims = int(n_sims)
			if n_sims < datlen: # use onl
				temp_df = self.parameter_df.iloc[-n_sims:]

		self.jucova = temp_df.cov().values
		#factor = (2.4 / np.sqrt(float(len(temp_df.columns))))**2
		factor = 2.4 / np.sqrt(float(len(temp_df.columns)))
		self.jucova = factor * self.jucova
		print(self.jucova)

	def set_jucova_scaling(self):
		#factor = (2.4/np.sqrt(float(len(self.parameter_df.columns))))**2
		factor = 2.4 / np.sqrt(float(len(self.parameter_df.columns)))
		#print('Factor: ', factor)
		self.jucova = factor*self.jucova
		#print('Covariance Matrix: \n', self.jucova)

	def save_to_file(self, savepath):
		self.parameter_df.to_json(savepath)

	def show(self):
		print(self.parameter_df)

	def show_children(self):
		for child in self.children:
			print(child.parameter_df)

	def replace_df(self, newdf):
		self.parameter_df = newdf

	def get_columns(self):
		return list(self.parameter_df.columns)

	def show(self):
		print(self.parameter_df)


class UniversalParameterGroup:
	def __init__(self, columns=None):
		self.parameter_df = pd.DataFrame(columns=columns)
		self.children = {}

	def declare_child(self, id, child_obj):
		self.children[id] = child_obj
		localnames = child_obj.get_columns()

		names = list(self.parameter_df.columns)

		for name in localnames:
			if id in name:
				names.append(name)
			else:
				names.append(id + '_' + name)

		self.parameter_df = pd.DataFrame(columns=names)

	def sync_from_children(self):
		self.parameter_df = self.parameter_df[0:0]
		for id, child in self.children.items():
			parameters = list(child.parameter_df.columns)
			for par in parameters:
				newpar = par
				if id not in par:
					newpar = id + '_' + par

				# print(child.parameter_df[par])

				self.parameter_df[newpar] = child.parameter_df[par]

	def show(self):
		print(self.parameter_df)

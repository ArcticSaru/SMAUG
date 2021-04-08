import tkinter as tk
import tkinter.ttk as ttk
import time
import json
import pandas as pd
from GUI_Elements import *
from Model import *
from tkinter import messagebox


class JSONEncoder(json.JSONEncoder):
	def default(self, obj):
		if hasattr(obj, 'to_json'):
			return obj.to_json(orient='records')
		return json.JSONEncoder.default(self, obj)


class MainWindow(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		self.modeldir = None
		self.outputdir = None
		self.nchains = 1
		self.modelstack = []

		self.runsicobi=False

		self.variogram = None
		self.analysisinfo = None
		self.analysisflag = False

		pd.options.mode.chained_assignment = None

		tk.Frame.__init__(self, self.parent)
		self.resultInfo = None
		self.initUI()

	def initUI(self):
		# Top level container
		self.bottomframe = BottomWindow(self)
		self.centreframe = CentreWindow(self)
		self.leftframe = LeftWindow(self)
		self.rightframe = RightWindow(self)
		self.leftbottomframe = LeftBottomWindow(self)

		# Layout top level frames
		self.leftframe.grid(row=0, column=0, sticky='ns')
		self.leftbottomframe.grid(row=1, column=0, sticky='nsew')
		self.centreframe.grid(row=0, column=1, sticky='nsew')
		self.bottomframe.grid(row=1, column=1, sticky='we')
		self.rightframe.grid(row=0, column=2, rowspan=2, sticky='ns')

		self.rowconfigure(0, weight=1)
		self.columnconfigure(1, weight=1)

	def loadnormal(self, modeldir=None):
		oldmodeldir = self.modeldir
		oldoutdir = self.outputdir
		if not modeldir:
			# Query model folder name
			self.modeldir = filedialog.askdirectory(initialdir='./Models')
		else:
			self.modeldir = modeldir

		if self.modeldir:
			self.outputdir = os.path.join(self.modeldir, 'Output')

			# Check model typei
			invInfofile = os.path.join(self.outputdir, 'invInfo.json')

			answer = True
			if os.path.exists(invInfofile):
				with open(invInfofile) as jf:
					invInfo = json.load(jf)
					if invInfo['Mode'] != 'normal':
						answer = tk.messagebox.askokcancel('Warning','If you continue your data will be erased (stored data is from a Sicobi-model)')

			if answer:
				self.resultInfo = {}
				self.runsicobi = False

				# Reset data container
				self.modelstack.clear()
				self.nchains = len(self.modelstack)
				self.variogram = None
				self.analysisinfo = None
				self.resultInfo = None
				self.analysisflag = False
				self.rightframe.empty_chains()
				self.bottomframe.reset_chains()

				gc.collect()

				# Set model title
				self.centreframe.set_title(os.path.basename(self.modeldir))

				data_error_flag = False
				# If possible load preexisting data
				if os.path.exists(self.outputdir):
					resauxfile = os.path.join(self.outputdir, 'result_chain_info.json')
					if os.path.exists(resauxfile):
						self.resultInfo = pd.read_json(resauxfile)
						result_filelist = [f for f in os.listdir(self.outputdir) if '_results.json' in f]
						if len(result_filelist) == len(self.resultInfo.index):
							for idf, f in enumerate(result_filelist):
								file = os.path.join(self.outputdir, f)
								self.add_chain()
								self.modelstack[idf].import_simulations(file, self.resultInfo.iloc[idf])

								# Update chain widgets
								self.nchains = len(self.modelstack)
								for idm, model in enumerate(self.modelstack, start=1):
									self.rightframe.update_chain_info(idm, model.nwarmupsteps, model.nmainsteps)

						else:
							data_error_flag = True
					else:
						data_error_flag = True
				else:
					data_error_flag = True

				# Add initial model, when not loading
				if data_error_flag:
					self.add_chain()

				self.analysisinfo = pd.DataFrame(columns=self.modelstack[0].get_parameternames(), index=['R_hat', 'n_eff'])
				self.bottomframe.nwustepstx.config(text=self.modelstack[0].nwarmupsteps)
				self.bottomframe.nmrstepstx.config(text=self.modelstack[0].nmainsteps)
				self.leftframe.set_parameter_list()
				self.leftframe.paroptionmenu.configure(state=tk.ACTIVE)
				self.bottomframe.change_modelstate(self.modelstack[0].state)
				self.centreframe.reset_figures()
				self.centreframe.draw()
			else:
				self.modeldir = oldmodeldir
				self.outputdir = oldoutdir

	def loadsicobi(self, modeldir = None):
		oldmodeldir = self.modeldir
		oldoutdir = self.outputdir
		if not modeldir:
			# Query model folder name
			self.modeldir = filedialog.askdirectory(initialdir='./Models')
		else:
			self.modeldir = modeldir

		if self.modeldir:
			self.outputdir = os.path.join(self.modeldir, 'Output')

			# Check model typei
			invInfofile = os.path.join(self.outputdir, 'invInfo.json')

			answer = True
			if os.path.exists(invInfofile):
				with open(invInfofile) as jf:
					invInfo = json.load(jf)
					if invInfo['Mode'] != 'sicobi':
						answer = tk.messagebox.askokcancel('Warning','If you continue your data will be erased (stored data is from a normal model)')

			if answer:
				self.resultInfo = {}
				self.runsicobi = True

				# Reset data container
				self.modelstack.clear()
				self.nchains = len(self.modelstack)
				self.variogram = None
				self.analysisinfo = None
				self.resultInfo = None
				self.analysisflag = False
				self.rightframe.empty_chains()
				self.bottomframe.reset_chains()

				gc.collect()

				# Set model title
				self.centreframe.set_title(os.path.basename(self.modeldir))

				data_error_flag = False
				# If possible load preexisting data
				if os.path.exists(self.outputdir):
					resauxfile = os.path.join(self.outputdir, 'result_chain_info.json')
					if os.path.exists(resauxfile):
						self.resultInfo = pd.read_json(resauxfile)
						result_filelist = [f for f in os.listdir(self.outputdir) if '_results.json' in f]
						for i in range(len(self.resultInfo.index)):
							self.add_chain()
							self.modelstack[i].import_simulations(str(i+1), self.outputdir, self.resultInfo.iloc[i])

						# Update chain widgets
						self.nchains = len(self.modelstack)
						for idm, model in enumerate(self.modelstack, start=1):
							self.rightframe.update_chain_info(idm, model.nwarmupsteps, model.nmainsteps)

					else:
						data_error_flag = True
				else:
					data_error_flag = True

				# Add initial model, when not loading
				if data_error_flag:
					self.add_chain()

				self.analysisinfo = pd.DataFrame(columns=self.modelstack[0].get_parameternames(), index=['R_hat', 'n_eff'])
				self.bottomframe.nwustepstx.config(text=self.modelstack[0].nwarmupsteps)
				self.bottomframe.nmrstepstx.config(text=self.modelstack[0].nmainsteps)
				self.leftframe.set_parameter_list()
				self.leftframe.paroptionmenu.configure(state=tk.ACTIVE)
				self.bottomframe.change_modelstate(self.modelstack[0].state)
				self.centreframe.reset_figures()
				self.centreframe.draw()
			else:
				self.modeldir = oldmodeldir
				self.outputdir = oldoutdir

	def warmup(self):
		# Initialise progress pop-up window
		w = self.winfo_reqwidth()
		h = self.winfo_reqheight()
		ws = self.winfo_screenwidth()
		hs = self.winfo_screenheight()
		x = (ws / 2) - (w / 2)
		y = (hs / 2) - (h / 2)

		nsim = self.bottomframe.get_simulation_number()
		ncha = len(self.modelstack)

		top = tk.Toplevel(self)
		top.title('im warming up')
		top.geometry('+%d+%d' % (x, y))

		# Make main window inactive
		top.grab_set()

		progtext = tk.Label(top, text='Chain %s / %s' % (0, ncha))
		progtext.grid(column=0, row=0, padx=10, pady=10)

		progbar = ttk.Progressbar(top, orient='horizontal', length=250, maximum=ncha)
		progbar.grid(column=0, row=1, padx=10, pady=10)

		if self.runsicobi:
			ncones = self.modelstack[0].get_nr_of_cones()

			progtext3 = tk.Label(top, text='Cone %s / %s' % (0, ncones), name='ptcone')
			progtext3.grid(column=0, row=2, padx=10, pady=10)

			progbar3 = ttk.Progressbar(top, orient='horizontal', length=250, maximum=ncones, name='pbcone')
			progbar3.grid(column=0, row=3, padx=10, pady=10)

		progtext2 = tk.Label(top, text='Initializing', name='ptsim')
		progtext2.grid(column=0, row=4, padx=10, pady=10)

		progbar2 = ttk.Progressbar(top, orient='horizontal', length=250, maximum=nsim, name='pbsim')
		progbar2.grid(column=0, row=5, padx=10, pady=10)

		top.update()

		# stopbtn = tk.Button(top, text='Stop', command=top.destroy)
		# stopbtn.grid(column=0, row=4, padx=10, pady=10)

		# Change state of models
		self.bottomframe.change_modelstate(ModelState.WARMUP)
		for idm, model in enumerate(self.modelstack, start=1):
			progtext.config(text='Chain %s / %s' % (idm, ncha))
			top.update()
			model.initialize()
			model.set_verbosity(0)
			model.change_state(ModelState.WARMUP)
			progbar.step()

		progbar['value'] = 0

		# Perform simulations for Sicobi model
		if self.runsicobi:
			for idm, model in enumerate(self.modelstack, start=1):
				progtext.config(text='Chain %s / %s' % (idm, ncha))
				model.MH_step(nsim, top, 'WU')  # Warmup
				progbar.step()
				top.update()
		else:
			for idm, model in enumerate(self.modelstack, start=1):
				progtext.config(text='Chain %s / %s' % (idm, ncha))
				for i in range(1, nsim + 1):
					if top.winfo_exists():
						progtext2.config(text='Simulating %s / %s' % (i, nsim))
						progbar2.step()
						top.update()
					else:
						top.grab_release()
						break
					model.MH_step()
				progbar.step()

		if top.winfo_exists():
			top.grab_release()
			top.destroy()

		print('Simulation finished')

		self.save_model()

		act_par = self.leftframe.get_lb_selection()
		if act_par:
			self.draw_parameter(act_par)

	def mainrun(self):
		# Initialise progress pop-up window
		w = self.winfo_reqwidth()
		h = self.winfo_reqheight()
		ws = self.winfo_screenwidth()
		hs = self.winfo_screenheight()
		x = (ws / 2) - (w / 2)
		y = (hs / 2) - (h / 2)

		nsim = self.bottomframe.get_simulation_number()
		ncha = len(self.modelstack)

		top = tk.Toplevel(self)
		top.title('im warming up')
		top.geometry('+%d+%d' % (x, y))

		# Make main window inactive
		top.grab_set()

		progtext = tk.Label(top, text='Chain %s / %s' % (0, ncha))
		progtext.grid(column=0, row=0, padx=10, pady=10)

		progbar = ttk.Progressbar(top, orient='horizontal', length=250, maximum=ncha)
		progbar.grid(column=0, row=1, padx=10, pady=10)

		if self.runsicobi:
			ncones = self.modelstack[0].get_nr_of_cones()

			progtext3 = tk.Label(top, text='Cone %s / %s' % (0, ncones), name='ptcone')
			progtext3.grid(column=0, row=2, padx=10, pady=10)

			progbar3 = ttk.Progressbar(top, orient='horizontal', length=250, maximum=ncones, name='pbcone')
			progbar3.grid(column=0, row=3, padx=10, pady=10)

		progtext2 = tk.Label(top, text='Initializing', name='ptsim')
		progtext2.grid(column=0, row=4, padx=10, pady=10)

		progbar2 = ttk.Progressbar(top, orient='horizontal', length=250, maximum=nsim, name='pbsim')
		progbar2.grid(column=0, row=5, padx=10, pady=10)

		top.update()

		# stopbtn = tk.Button(top, text='Stop', command=top.destroy)
		# stopbtn.grid(column=0, row=4, padx=10, pady=10)

		# Change state of models
		self.bottomframe.change_modelstate(ModelState.MAINRUN)
		for idm, model in enumerate(self.modelstack, start=1):
			progtext.config(text='Chain %s / %s' % (idm, ncha))
			top.update()
			model.initialize()
			model.set_verbosity(0)
			model.change_state(ModelState.MAINRUN)
			progbar.step()

		progbar['value'] = 0

		# Perform simulations for Sicobi model
		if self.runsicobi:
			for idm, model in enumerate(self.modelstack, start=1):
				progtext.config(text='Chain %s / %s' % (idm, ncha))
				model.MH_step(nsim, top, 'MR')  # Mainrun
				progbar.step()
				top.update()
		else:
			for idm, model in enumerate(self.modelstack, start=1):
				progtext.config(text='Chain %s / %s' % (idm, ncha))
				for i in range(1, nsim + 1):
					if top.winfo_exists():
						progtext2.config(text='Simulating %s / %s' % (i, nsim))
						progbar2.step()
						top.update()
					else:
						top.grab_release()
						break
					model.MH_step()
				progbar.step()

		if top.winfo_exists():
			top.grab_release()
			top.destroy()

		print('Simulation finished')

		self.save_model()

		act_par = self.leftframe.get_lb_selection()
		if act_par:
			self.draw_parameter(act_par)

	def save_model(self):
		# Update chain widgets and resultInfo Data Frame
		res_list = []
		for idm, model in enumerate(self.modelstack, start=1):
			self.rightframe.update_chain_info(idm, model.nwarmupsteps, model.nmainsteps)

			res_list.append({'n_warmup_sim': model.nwarmupsteps, 'n_main_sim': model.nmainsteps})

		self.resultInfo = pd.DataFrame(res_list)

		# Write resultInfo DF to file
		self.resultInfo.to_json(os.path.join(self.outputdir, 'result_chain_info.json'))

		# Write results to file
		for idm, model in enumerate(self.modelstack, start=1):
			model.save_to_file(os.path.join(self.outputdir, str(idm) + '_results.json'))

		mode = 'normal'
		if self.runsicobi:
			mode = 'sicobi'

		inversionInfo = {'Mode': mode}

		with open(os.path.join(self.outputdir, 'invInfo.json'), 'w') as fp:
			json.dump(inversionInfo, fp)

	def reset_model(self):
		if self.modelstack:
			filelist = [f for f in os.listdir(self.outputdir)]
			for f in filelist:
				os.remove(os.path.join(self.outputdir, f))

			if self.runsicobi:
				self.loadsicobi(self.modeldir)
			else:
				self.loadnormal(self.modeldir)

	def draw_parameter(self, parameter_name):
		# Clear canvas
		self.centreframe.reset_figures()

		# Fetch drawing parameter
		drawpar = self.bottomframe.drawchainvar.get()

		# Draw Chains
		if drawpar == 'All':
			for model in self.modelstack:
				pardat = model.get_data(parameter_name)
				# print(pardat)
				self.centreframe.draw_chains(pardat, model.nwarmupsteps, model.nmainsteps)
		else:
			index = int(drawpar) - 1
			pardat = self.modelstack[index].get_data(parameter_name)
			self.centreframe.draw_chains(pardat, self.modelstack[index].nwarmupsteps, self.modelstack[index].nmainsteps)

		# Draw Histogram
		merged_df = pd.Series(name=parameter_name, dtype=float)
		for model in self.modelstack:
			#print(model.get_mr_data(parameter_name))
			merged_df = merged_df.append(model.get_mr_data(parameter_name), ignore_index=True)

		self.centreframe.draw_histogram(merged_df)

		# Draw Variogram
		if self.variogram is not None:
			self.centreframe.draw_variogram(self.variogram[[parameter_name]])

		self.centreframe.draw(parameter_name)

		# Upsate parameter info in bottom-left frame
		self.leftbottomframe.show_parameter_stats(merged_df)

	def adcovma(self):
		for model in self.modelstack:
			model.adapt_proposal_dist(self.bottomframe.adcovmaentry.get())

	def save_stat_parameters(self):
		# Get active parameters
		pars = self.leftframe.parlistbox.lists[0].get(0, tk.END)

		parrhat = self.leftframe.parlistbox.lists[1].get(0, tk.END)
		parneff = self.leftframe.parlistbox.lists[2].get(0, tk.END)

		pardata = pd.DataFrame()

		if pars:
			for parname in pars:
				merged_ser = pd.Series(name=parname)
				for model in self.modelstack:
					# print(model.get_mr_data(parameter_name))
					merged_ser = merged_ser.append(model.get_mr_data(parname), ignore_index=True)

				if bool(re.search('_r[0-9]', parname)):
					# Get total cone length
					totcole = self.get_cone_length(parname)

					ratios = merged_ser.apply(lambda x: 10 ** x)
					factor = ratios.div(ratios + 1)
					factor2 = 1 / (ratios + 1)

					merged_ser = totcole * factor

				pardata[parname] = merged_ser

			mean = pardata.mean()
			quantiles = pardata.quantile([0.16, 0.84])

			filename = filedialog.asksaveasfilename(filetypes=(("CSV Files", "*.csv"),))

			if filename:
				mean.to_csv(filename + '_mean.csv', header=True, encoding='utf-16')
				quantiles.to_csv(filename + '_quantiles.csv', header=True, encoding='utf-16')

	def save_stat_positions(self):
		if self.leftframe.parvar.get() == 'Thickness':
			# Get active parameters
			pars = self.leftframe.parlistbox.lists[0].get(0, tk.END)

			parrhat = self.leftframe.parlistbox.lists[1].get(0, tk.END)
			parneff = self.leftframe.parlistbox.lists[2].get(0, tk.END)

			pardata = pd.DataFrame()

			if pars:
				for parname in pars:
					merged_ser = pd.Series(name=parname)
					for model in self.modelstack:
						# print(model.get_mr_data(parameter_name))
						merged_ser = merged_ser.append(model.get_mr_data(parname), ignore_index=True)

					if bool(re.search('_r[0-9]', parname)):
						# Get total cone length
						totcole = self.get_cone_length(parname)

						ratios = merged_ser.apply(lambda x: 10 ** x)
						factor = ratios.div(ratios + 1)
						factor2 = 1 / (ratios + 1)

						merged_ser = totcole * factor

					pardata[parname] = merged_ser

				quantiles = pardata.quantile([0.16, 0.84])

				# Get detector information
				detectors = self.modelstack[0].get_detectors()

				# Get data information
				data = self.modelstack[0].get_data_info()

				positions = pd.DataFrame(columns=['x_lower', 'y_lower', 'z_lower', 'x_upper', 'y_upper', 'z_upper'])

				for col in quantiles.columns:
					itentifier_split = col.split(sep='_')
					detname = itentifier_split[0]
					cone_id = itentifier_split[1]
					cone_name = detname + '_' + cone_id

					# Get datainfo row
					cone_info = data.where(data['ID'] == cone_name).dropna()
					# print(cone_info, '\n')

					# Get rest of parameters and calculate coordinates of +/- 1 Sigma bounds for position
					det = detectors[detname]
					R0 = [float(det.loc['E (m, CH1903)'].value), float(det.loc['N (m, CH1903)'].value), float(det.loc['Z (m, CH1903)'].value)]

					theta_geo = cone_info[chr(952)].iloc[0]
					phi_geo = cone_info[chr(966)].iloc[0]

					r0 = quantiles.iloc[0][col]
					r1 = quantiles.iloc[1][col]

					theta = np.deg2rad(theta_geo)
					phi = np.deg2rad(90 - phi_geo)

					dx0 = r0*np.sin(theta)*np.cos(phi)
					dy0 = r0*np.sin(theta)*np.sin(phi)
					dz0 = r0*np.cos(theta)

					dx1 = r1*np.sin(theta)*np.cos(phi)
					dy1 = r1*np.sin(theta)*np.sin(phi)
					dz1 = r1*np.cos(theta)

					new_row = [R0[0]+dx0, R0[1]+dy0, R0[2]+dz0, R0[0]+dx1, R0[1]+dy1, R0[2]+dz1]

					posSeries = pd.Series(new_row, index=positions.columns, name=cone_name)
					positions = positions.append(posSeries)

				filename = filedialog.asksaveasfilename(filetypes=(("CSV Files", "*.csv"),))

				if filename:
					positions.to_csv(filename + '_cone_positions.csv', header=True, encoding='utf-16')

		else:
			tk.messagebox.showinfo('Warning','Please switch parameters to \" Thickness \" ')

	def save_full_parameters(self):
		# Get active parameters
		pars = self.leftframe.parlistbox.lists[0].get(0, tk.END)

		parrhat = self.leftframe.parlistbox.lists[1].get(0, tk.END)
		parneff = self.leftframe.parlistbox.lists[2].get(0, tk.END)

		pardata = pd.DataFrame()

		if pars:
			for parname in pars:
				merged_ser = pd.Series(name=parname)
				for model in self.modelstack:
					# print(model.get_mr_data(parameter_name))
					merged_ser = merged_ser.append(model.get_mr_data(parname), ignore_index=True)

				if bool(re.search('_r[0-9]', parname)):
					# Get total cone length
					totcole = self.get_cone_length(parname)

					ratios = merged_ser.apply(lambda x: 10 ** x)
					factor = ratios.div(ratios + 1)
					factor2 = 1 / (ratios + 1)

					merged_ser = totcole * factor

				pardata[parname] = merged_ser

			filename = filedialog.asksaveasfilename(filetypes=(("CSV Files", "*.csv"),))

			if filename:
				pardata.to_csv(filename + 'full_simulation_data.csv', header=True, encoding='utf-16')

	def save_full_positions(self):
		if self.leftframe.parvar.get() == 'Thickness':
			# Get active parameters
			pars = self.leftframe.parlistbox.lists[0].get(0, tk.END)

			parrhat = self.leftframe.parlistbox.lists[1].get(0, tk.END)
			parneff = self.leftframe.parlistbox.lists[2].get(0, tk.END)

			pardata = pd.DataFrame()

			if pars:
				# Load and merge data
				for parname in pars:
					merged_ser = pd.Series(name=parname)
					for model in self.modelstack:
						# print(model.get_mr_data(parameter_name))
						merged_ser = merged_ser.append(model.get_mr_data(parname), ignore_index=True)

					if bool(re.search('_r[0-9]', parname)):
						# Get total cone length
						totcole = self.get_cone_length(parname)

						ratios = merged_ser.apply(lambda x: 10 ** x)
						factor = ratios.div(ratios + 1)
						factor2 = 1 / (ratios + 1)

						merged_ser = totcole * factor

					pardata[parname] = merged_ser

				cone_dfs = {}

				# print(pardata)

				# Get detector information
				detectors = self.modelstack[0].get_detectors()

				# Get data information
				data = self.modelstack[0].get_data_info()

				# Create Cone dataframes
				for col in pardata.columns:
					itentifier_split = col.split(sep='_')
					detname = itentifier_split[0]
					cone_id = itentifier_split[1]
					cone_name = detname + '_' + cone_id

					datalist = []

					# Get datainfo row
					cone_info = data.where(data['ID'] == cone_name).dropna()
					# print(cone_info, '\n')

					# Get rest of parameters and calculate coordinates of +/- 1 Sigma bounds for position
					det = detectors[detname]
					R0 = [float(det.loc['E (m, CH1903)'].value), float(det.loc['N (m, CH1903)'].value),
					      float(det.loc['Z (m, CH1903)'].value)]

					theta_geo = cone_info[chr(952)].iloc[0]
					phi_geo = cone_info[chr(966)].iloc[0]

					theta = np.deg2rad(theta_geo)
					phi = np.deg2rad(90 - phi_geo)

					for i, row_value in pardata[col].iteritems():
						dx = row_value * np.sin(theta) * np.cos(phi)
						dy = row_value * np.sin(theta) * np.sin(phi)
						dz = row_value * np.cos(theta)

						X = R0[0] + dx
						Y = R0[1] + dy
						Z = R0[2] + dz

						datalist.append([X, Y, Z])

					cone_dfs[cone_name] = pd.DataFrame(datalist, columns=['x', 'y', 'z'])

				filename = filedialog.asksaveasfilename()

				if filename:
					with open(filename + '.json', 'w') as fp:
						json.dump(cone_dfs, fp, cls=JSONEncoder)

		else:
			tk.messagebox.showinfo('Warning', 'Please switch parameters to \" Thickness \" ')

	def analysis(self):
		fullvariogram = bool(self.bottomframe.variogramvar.get())
		print('analyse some stuff')
		# Prepare master DFs
		#columns = self.modelstack[0].parameter_df.columns
		meandf = pd.DataFrame()
		vardf = pd.DataFrame()

		chunksize = 0

		# Loop through all chains
		for model in self.modelstack:
			# Split chain in multiple segments of equal length
			# Get all main run simulations
			mrdf = model.get_mr_data()

			nsplit = self.bottomframe.nsubcha.get()
			chunksize = len(mrdf)//nsplit
			maxindex = chunksize*nsplit

			mrdf.drop(mrdf.tail(len(mrdf)-maxindex).index, inplace=True)

			# print('int divisor: ', chunksize)
			# print('max index: ', maxindex)
			# print('drop last ', len(mrdf)-maxindex, ' rows')

			# for (columnName, columnItem) in anadf.iteritems():
			# length_mr = mrdf[columnName].iloc[-self.model.nmainsteps:].shape[0]
			# print(len(mrdf))

			# Loop thorugh all subchains
			for g, df in mrdf.groupby(np.arange(len(mrdf)) // chunksize):
				# print(df.mean(axis=0))
				# print(df.var(axis=0))
				meandf = meandf.append(df.mean(axis=0), ignore_index=True)
				vardf = vardf.append(df.var(axis=0), ignore_index=True)
				# print(g, ': ', len(df))

		W = vardf.mean(axis=0)
		B = meandf.var(axis=0)*chunksize
		postvar = (chunksize-1)/chunksize * W + 1/chunksize * B
		R_hat2 = postvar/W
		R_hat = R_hat2.transform(np.sqrt)
			# print(mrdf[columnName].iloc[-self.model.nmainsteps:].shape[0])

		print('Mean: \n', meandf)
		print('Variance: \n', vardf)
		print('W: \n', W)
		print('B: \n', B)
		print('PVar: \n', postvar)
		print('R_hat: \n', R_hat)
		print('\n')

		boolser = pd.Series(1, index=mrdf.columns)

		# Compute Variogram
		self.variogram = pd.DataFrame()
		self.analysisinfo = pd.DataFrame(columns=mrdf.columns)

		print(chunksize)

		for lag in range(1, chunksize):
			seriodf = pd.DataFrame()

			for model in self.modelstack:
				mrdf = model.get_mr_data()

				nsplit = self.bottomframe.nsubcha.get()
				chunksize = len(mrdf) // nsplit
				maxindex = chunksize * nsplit

				mrdf.drop(mrdf.tail(len(mrdf) - maxindex).index, inplace=True)

				for g, df in mrdf.groupby(np.arange(len(mrdf)) // chunksize):
					seriodf = seriodf.append(df.diff(periods=lag).dropna().transform(lambda x: x*x).mean(axis=0), ignore_index=True)

			self.variogram = self.variogram.append(seriodf.mean(axis=0).div(-2*postvar).add(1), ignore_index=True)

			tailsum = self.variogram.tail(2).sum() > 0

			if lag % 2 == 0:
				boolser *= tailsum.astype('int')

			print('lag: ', lag)
			if not fullvariogram:
				if not boolser.any():
					break

		# print(self.variogram)

		neffdict = {}
		for value in mrdf.columns:
			autocorrsum = 0
			for idr, valr in self.variogram[value].iteritems():
				autocorrsum += valr
				if idr % 2 != 0:
					checksum = self.variogram[value][idr + 1] + self.variogram[value][idr + 2]
					if checksum < 0 or idr >= len(self.variogram) - 4:
						# print(value, ': ', idr, ' || rho_hat: ', autocorrsum)
						break

			neffdict[value] = len(self.modelstack) * maxindex / (1+2*autocorrsum)

		neffser = pd.Series(data=neffdict)

		# Draw information in left frame
		self.analysisinfo = self.analysisinfo.append(R_hat, ignore_index=True)
		self.analysisinfo = self.analysisinfo.append(neffser, ignore_index=True)
		self.analysisinfo = self.analysisinfo.rename({0: 'R_hat', 1: 'n_eff'})

		self.draw_analysisinfo()

	def draw_analysisinfo(self):
		self.leftframe.draw_analysisinfo(self.analysisinfo)

	def get_model_parameternames(self):
		return self.model.get_parameternames()

	def get_cone_length(self, cone):
		return self.modelstack[0].get_cone_length(cone)

	def add_chain(self):
		if self.modeldir:
			if self.runsicobi:
				self.modelstack.append(SingleConeModelCluster(self.modeldir))
			else:
				self.modelstack.append(Model(self.modeldir))

			self.rightframe.add_chain()
			self.bottomframe.add_chain()

		# print(self.modelstack)

	def remove_chain(self):
		if len(self.modelstack) > 1:
			self.modelstack.pop()
			self.rightframe.remove_chain()
			self.bottomframe.remove_chain()
			gc.collect()

		# print(self.modelstack)


def main():
	root = tk.Tk()
	root.title("Inversion")
	screen_width = root.winfo_screenwidth()
	screen_height = root.winfo_screenheight()
	sizestring = str(int(0.8 * screen_width)) + 'x' + str(int(0.8 * screen_height)) + '+50+50'
	root.geometry(sizestring)

	window = MainWindow(master=root)
	window.pack(side="top", fill="both", expand=True)
	root.mainloop()


if __name__ == '__main__':
	main()

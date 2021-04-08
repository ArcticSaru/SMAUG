import tkinter as tk
import tkinter.ttk as ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.ticker import MaxNLocator
from matplotlib.backend_bases import key_press_handler
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from tkinter import filedialog
from tkinter import ttk
from tkinter import messagebox
import tkinter.font as tkFont
import numpy as np
import pandas as pd
import logging
import re
from Model import ModelState


class BottomWindow(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=7)

		self.drawoptions = ['All']

		self.nchainsvar = tk.IntVar()
		self.nchainsvar.set(1)

		self.cbvar = tk.IntVar()

		self.nsubcha = tk.IntVar()
		self.nsubcha.set(2)

		self.simsteps = tk.StringVar()
		self.simsteps.set(500)

		self.adcovmastr = tk.StringVar()
		self.adcovmastr.set(100)

		self.drawchainvar = tk.StringVar()
		self.drawchainvar.set('All')

		self.variogramvar = tk.IntVar()
		self.variogramvar.set(0)

		self.bold_font = tkFont.Font(self, family='Helvetica', size=12, weight=tkFont.BOLD)

		self.__createlayout()

	def __createlayout(self):
		self.label = tk.Label(self, text='Actions', font=self.bold_font)
		self.loadmodelbutton = tk.Button(self, text='Load model', command=self.parent.loadnormal)
		self.loadsicobibutton = tk.Button(self, text='Load Sicobi', command= self.parent.loadsicobi)
		self.warmupbutton = tk.Button(self, text='Start warmup', state=tk.DISABLED, command=self.parent.warmup)
		self.mainrunbutton = tk.Button(self, text='Start main run', state=tk.DISABLED, command=self.parent.mainrun)
		self.nsimlb = tk.Label(self, text='# simulation steps: ')
		self.simentry = tk.Entry(self, textvariable=self.simsteps)
		self.intseparator = ttk.Separator(self, orient=tk.HORIZONTAL)
		self.adcovmatbtutton = tk.Button(self, text='Adapt Cov-matrix', state=tk.DISABLED, command=self.parent.adcovma)
		self.adcovmatext1 = tk.Label(self, text='acc. to last')
		self.adcovmaentry = tk.Entry(self, textvariable=self.adcovmastr)
		self.adcovmatext2 = tk.Label(self, text='simulations')

		self.separator = ttk.Separator(self, orient=tk.VERTICAL)

		self.modelstatelb = tk.Label(self, text='Model state:')
		self.modelstatetx = tk.Label(self, text='No model')
		self.nwustepslb = tk.Label(self, text='# warmup steps:')
		self.nwustepstx = tk.Label(self, text='0')
		self.nmrstepslb = tk.Label(self, text='# main run steps:')
		self.nmrstepstx = tk.Label(self, text='0')
		self.resetbutton = tk.Button(self, text='Reset simulations', state=tk.DISABLED, command=self.areyousuredialog)

		self.separator2 = ttk.Separator(self, orient=tk.VERTICAL)

		self.analysislb = tk.Label(self, text='Main run analysis')
		self.subchainlb = tk.Label(self, text='Number of chain divisions: ')
		self.subchainentry = tk.Entry(self, textvariable=self.nsubcha)
		self.analysisbutton = tk.Button(self, text='Analyse', state=tk.DISABLED, command=self.parent.analysis)
		self.variogramcheckbox = tk.Checkbutton(self, text='Draw full variogram', variable=self.variogramvar)
		self.variogramtext = tk.Label(self, text='(Can take a long time!) ')
		self.exportlabel = tk.Label(self, text='Export active parameters:')
		self.exportparstat = tk.Button(self, text='Statistical', command=self.parent.save_stat_parameters)
		self.exportparfull = tk.Button(self, text='Full', command=self.parent.save_full_parameters)
		self.exportposlabel = tk.Label(self, text='Export interface positions:')
		self.exportposstat = tk.Button(self, text='Statistical', command=self.parent.save_stat_positions)
		self.exportposfull = tk.Button(self, text='DEM', command=self.parent.save_full_positions)

		self.separator3 = ttk.Separator(self, orient=tk.VERTICAL)

		self.chainlabel = tk.Label(self, text='MCMC Chain management')
		self.addchainbtn = tk.Button(self, text='Add Chain', state=tk.DISABLED, command=self.parent.add_chain)
		self.removechainbtn = tk.Button(self, text='Remove Chain', state=tk.DISABLED, command=self.parent.remove_chain)
		self.drawchaintext = tk.Label(self, text='Draw chains: ')
		self.drawchainoptm = tk.OptionMenu(self, self.drawchainvar, 'All')

		self.label.grid(column=0, row=0, sticky='nw')
		self.loadmodelbutton.grid(column=0, row=1, sticky='w', padx=2, pady=2)
		self.loadsicobibutton.grid(column=0, row=2, sticky='w', padx=2, pady=2)
		self.warmupbutton.grid(column=1, row=1, padx=2, pady=2)
		self.mainrunbutton.grid(column=2, row=1, padx=2, pady=2)
		self.nsimlb.grid(column=0, row=3, padx=2, pady=2)
		self.simentry.grid(column=1, row=3, padx=2, pady=5)
		self.intseparator.grid(column=0, row=4, columnspan=4, sticky='ew', pady=5)
		self.adcovmatbtutton.grid(column=0, row=5, padx=2, pady=2)
		self.adcovmatext1.grid(column=1, row=5, padx=2, pady=2)
		self.adcovmaentry.grid(column=2, row=5, padx=2, pady=2)
		self.adcovmatext2.grid(column=3, row=5, padx=2, pady=2)

		self.separator.grid(column=4, row=0, rowspan=6, sticky='ns', padx=20)

		self.modelstatelb.grid(column=5, row=0, padx=2, pady=2)
		self.modelstatetx.grid(column=6, row=0, padx=2, pady=2)
		# self.nwustepslb.grid(column=5, row=2, padx=2, pady=2)
		# self.nwustepstx.grid(column=6, row=2, padx=2, pady=2)
		# self.nmrstepslb.grid(column=5, row=3, padx=2, pady=2)
		# self.nmrstepstx.grid(column=6, row=3, padx=2, pady=2)
		self.resetbutton.grid(column=6, row=5, padx=2, pady=20)

		self.separator2.grid(column=7, row=0, rowspan=6, sticky='ns', padx=20)

		self.analysislb.grid(column=8, row=0, padx=2, pady=10, sticky='w')
		self.subchainlb.grid(column=8, row=1, padx=2, pady=2, sticky='w')
		self.subchainentry.grid(column=9, row=1, padx=2, pady=2, sticky='w')
		self.analysisbutton.grid(column=8, row=4, padx=2, pady=2, columnspan=2)
		self.variogramcheckbox.grid(column=8, row=2, padx=2, pady=2, sticky='w')
		self.variogramtext.grid(column=9, row=2, sticky='w')
		self.exportlabel.grid(column=8, row=5, padx=2, pady=2)
		self.exportparstat.grid(column=9, row=5, padx=2, pady=2)
		self.exportparfull.grid(column=10, row=5, padx=2, pady=2)
		self.exportposlabel.grid(column=8, row=6, padx=2, pady=2)
		self.exportposstat.grid(column=9, row=6, padx=2, pady=2)
		self.exportposfull.grid(column=10, row=6, padx=2, pady=2)

		self.separator3.grid(column=11, row=0, rowspan=6, sticky='ns', padx=20)

		self.chainlabel.grid(column=12, row=0, padx=2, pady=2)
		self.addchainbtn.grid(column=12, row=2, padx=2, pady=2)
		self.removechainbtn.grid(column=13, row=2, padx=2, pady=2)
		self.drawchaintext.grid(column=12, row=4, padx=2, pady=2)
		self.drawchainoptm.grid(column=12, row=5, padx=10, pady=10)

	def change_modelstate(self, modelstate):
		if modelstate == ModelState.INITIAL:
			self.modelstatetx.config(text='Initial')
			self.warmupbutton.config(state=tk.ACTIVE)
			self.mainrunbutton.config(state=tk.DISABLED)
			self.resetbutton.config(state=tk.DISABLED)
			self.adcovmatbtutton.config(state=tk.DISABLED)
			self.analysisbutton.config(state=tk.DISABLED)
			self.addchainbtn.config(state=tk.ACTIVE)
			self.removechainbtn.config(state=tk.ACTIVE)
		elif modelstate == ModelState.WARMUP:
			self.modelstatetx.config(text='Warm up')
			self.warmupbutton.config(state=tk.ACTIVE)
			self.mainrunbutton.config(state=tk.ACTIVE)
			self.resetbutton.config(state=tk.ACTIVE)
			self.adcovmatbtutton.config(state=tk.ACTIVE)
			self.analysisbutton.config(state=tk.DISABLED)
			self.addchainbtn.config(state=tk.DISABLED)
			self.removechainbtn.config(state=tk.DISABLED)
		elif modelstate == ModelState.MAINRUN:
			self.modelstatetx.config(text='Main run')
			self.warmupbutton.config(state=tk.DISABLED)
			self.mainrunbutton.config(state=tk.ACTIVE)
			self.resetbutton.config(state=tk.ACTIVE)
			self.adcovmatbtutton.config(state=tk.ACTIVE)
			self.analysisbutton.config(state=tk.ACTIVE)
			self.addchainbtn.config(state=tk.DISABLED)
			self.removechainbtn.config(state=tk.DISABLED)

	def areyousuredialog(self):
		result = tk.messagebox.askquestion('Delete', 'Do you really want to delete all simulations from this model?',
										   icon='warning')
		if result == 'yes':
			self.parent.reset_model()

	def get_simulation_number(self):
		return int(self.simsteps.get())

	def add_chain(self):
		self.drawoptions.append(str(len(self.drawoptions)))

		self.drawchainvar.set('All')
		self.drawchainoptm['menu'].delete(0, tk.END)

		for choice in self.drawoptions:
			self.drawchainoptm['menu'].add_command(label=choice, command=tk._setit(self.drawchainvar, choice))

	def remove_chain(self):
		if len(self.drawoptions) > 2:
			self.drawoptions.pop()
			self.drawchainoptm['menu'].delete(tk.END)

	def reset_chains(self):
		self.drawoptions = ['All']


class RightWindow(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)

		self.ncha = 0

		self.__createlayout()

	def __createlayout(self):
		self.label = tk.Label(self, text='MCMC_Chains')
		self.chainlistbox = ChainListbox(self, [('Name', 20), ('#WU', 10), ('#MR', 10)])

		self.label.grid(column=0, row=0, sticky='nw')
		self.chainlistbox.grid(column=0, row=1, sticky='nsew')

		self.rowconfigure(1, weight=1)
		self.columnconfigure(0, weight=1)

	def add_chain(self):
		self.ncha += 1
		name = 'chain ' + str(self.ncha)
		self.chainlistbox.append([name, 0, 0])

	def remove_chain(self):
		self.ncha -= 1
		self.chainlistbox.pop()

	def update_chain_info(self, chain_number, nwus, nmrs):
		name = 'chain ' + str(chain_number)
		self.chainlistbox.delete(chain_number-1)
		self.chainlistbox.insert(chain_number-1, [name, nwus, nmrs])

	def empty_chains(self):
		self.ncha = 0
		self.chainlistbox.delete_all()


class LeftBottomWindow(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1,
		                  highlightcolor='black', highlightbackground='black', highlightthickness=1)

		self.__createlayout()

	def __createlayout(self):
		self.analysislabel = tk.Label(self, text='Parameter information')
		self.meantxt = tk.Label(self, text='Mean: ')
		self.meanval = tk.Label(self, text='')
		self.mediantxt = tk.Label(self, text='Median: ')
		self.medianval = tk.Label(self, text='')
		self.onesigCItxt = tk.Label(self, text='68% central interval: ')
		self.onesigCIval = tk.Label(self, text='[]')
		self.twosigCItxt = tk.Label(self, text='95% central interval: ')
		self.twosigCIval = tk.Label(self, text='[]')
		self.sigtxt = tk.Label(self, text='Std-dev: ')
		self.sigval = tk.Label(self, text='')

		self.analysislabel.grid(column=0, row=0, padx=2, pady=5, columnspan=2)
		self.meantxt.grid(column=0, row=1, padx=2, pady=2, sticky='w')
		self.meanval.grid(column=1, row=1, padx=10, pady=2, sticky='e')
		self.mediantxt.grid(column=0, row=2, padx=2, pady=2, sticky='w')
		self.medianval.grid(column=1, row=2, padx=10, pady=2, sticky='e')
		self.onesigCItxt.grid(column=0, row=3, padx=2, pady=2, sticky='w')
		self.onesigCIval.grid(column=1, row=3, padx=10, pady=2, sticky='e')
		self.twosigCItxt.grid(column=0, row=4, padx=2, pady=2, sticky='w')
		self.twosigCIval.grid(column=1, row=4, padx=10, pady=2, sticky='e')
		self.sigtxt.grid(column=0, row=5, padx=2, pady=2, sticky='w')
		self.sigval.grid(column=1, row=5, padx=10, pady=2, sticky='e')

		self.columnconfigure(0, weight=1)

	def show_parameter_stats(self, parameter_df):
		if bool(re.search('_r[0-9]', parameter_df.name)):
			# Get total cone length
			totcole = self.parent.get_cone_length(parameter_df.name)

			ratios = parameter_df.apply(lambda x: 10 ** x)
			factor = ratios.div(ratios + 1)
			factor2 = 1 / (ratios + 1)

			lengths = totcole * factor

			mean = lengths.mean()
			std = lengths.std()
			quantiles = lengths.quantile([0.025, 0.16, 0.5, 0.84, 0.975])

		else:
			mean = parameter_df.mean()
			std = parameter_df.std()
			quantiles = parameter_df.quantile([0.025, 0.16, 0.5, 0.84, 0.975])

		# Use scientific notation for flux parameters
		if bool(re.search(chr(0x03BB), parameter_df.name)):
			self.meanval.config(text='{:.3E}'.format(mean))
			self.medianval.config(text='{:.3E}'.format(quantiles.loc[0.5]))
			self.onesigCIval.config(text='[' + '{:.3E}'.format(quantiles.loc[0.16]) + ' , '
			                             + '{:.3E}'.format(quantiles.loc[0.84]) + ']')
			self.twosigCIval.config(text='[' + '{:.3E}'.format(quantiles.loc[0.025]) + ' , '
			                             + '{:.3E}'.format(quantiles.loc[0.975]) + ']')
			self.sigval.config(text='{:.3E}'.format(std))
		else:
			self.meanval.config(text='{:.3f}'.format(mean))
			self.medianval.config(text='{:.3f}'.format(quantiles.loc[0.5]))
			self.onesigCIval.config(text='[' + '{:.3f}'.format(quantiles.loc[0.16]) + ' , '
			                             + '{:.3f}'.format(quantiles.loc[0.84]) + ']')
			self.twosigCIval.config(text='[' + '{:.3f}'.format(quantiles.loc[0.025]) + ' , '
			                             + '{:.3f}'.format(quantiles.loc[0.975]) + ']')
			self.sigval.config(text='{:.3f}'.format(std))


class LeftWindow(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)

		self.parvar = tk.StringVar(self)
		self.parvar.set('All')
		self.parvar.trace('w', self.set_parameter_list)

		self.__createlayout()

	def __createlayout(self):
		self.label = tk.Label(self, text='Parameters')
		self.paroptionmenu = tk.OptionMenu(self, self.parvar, 'All', 'Composition', 'Density', 'EL_Error', 'Flux', 'Thickness')
		self.paroptionmenu.configure(state=tk.DISABLED)
		self.parlistbox = AnalysisListbox(self, [('Name', 20), ('R_hat', 10), ('n_eff', 10)])

		# self.separator = ttk.Separator(self, orient=tk.HORIZONTAL)

		self.label.grid(column=0, row=0, sticky='nw')
		self.paroptionmenu.grid(column=1, row=0, sticky='ne')
		self.parlistbox.grid(column=0, row=1, sticky='nsew', columnspan=2)

		# self.separator.grid(column=0, row=2, sticky='ew', pady=20, columnspan=2)

		self.rowconfigure(1, weight=1)

	def set_parameter_list(self, *args):
		parlist = self.parent.analysisinfo.columns
		self.parlistbox.delete_all()

		pava = self.parvar.get()
		sel_df = self.parent.analysisinfo

		if pava == 'Composition':
			sel_df = sel_df.filter(regex='c_')
		elif pava == 'Density':
			sel_df = sel_df.filter(regex=chr(0x03C1))
		elif pava == 'EL_Error':
			sel_df = sel_df.filter(regex=chr(0x03C3))
		elif pava == 'Flux':
			sel_df = sel_df.filter(regex=chr(0x03BB))
		elif pava == 'Thickness':
			sel_df = sel_df.filter(regex='_r[0-9]')

		self.parlistbox.insert_df(sel_df)
		self.parent.update()

	def draw_parameter(self, parameter_name):
		self.parent.draw_parameter(parameter_name)

	def get_lb_selection(self):
		index = self.parlistbox.curselection()
		if index:
			return self.parlistbox.get(self.parlistbox.curselection())[0]
		else:
			return None

	def draw_analysisinfo(self, infodf):
		saveidx = self.parlistbox.curselection()
		self.parlistbox.delete_all()
		self.parlistbox.insert_df(infodf)
		self.parlistbox.selection_reset(saveidx)


class CentreWindow(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)
		self.__createlayout()

	def __createlayout(self):
		self.modelname = tk.Label(self, text='Model: ')

		self.fig = plt.figure()
		self.spec = gridspec.GridSpec(ncols=2, nrows=2, wspace=0.5, hspace=0.5, figure=self.fig)
		self.ax1 = self.fig.add_subplot(self.spec[0, 0])
		self.ax2 = self.fig.add_subplot(self.spec[0, 1])
		self.ax3 = self.fig.add_subplot(self.spec[1, 0])
		self.ax4 = self.fig.add_subplot(self.spec[1, 1])

		self.canvas = FigureCanvasTkAgg(self.fig, master=self)

		self.modelname.pack(side=tk.TOP, anchor=tk.NW)
		self.canvas.get_tk_widget().pack(fill='both', expand=True)

		self.draw()

	def set_title(self, name):
		self.modelname.config(text='Model: {}'.format(name))

	def draw_chains(self, parameter_data, nwustep, nmrstep):

		if nwustep > 0:
			parameter_data.iloc[:nwustep + 1].plot(ax=self.ax1)

		if nmrstep > 0:
			parameter_data.iloc[nwustep + 1:].plot(ax=self.ax2)

	def draw_histogram(self, parameter_data):
		if bool(re.search('_r[0-9]', parameter_data.name)):

			# Get total cone length
			totcole = self.parent.get_cone_length(parameter_data.name)

			ratios = parameter_data.apply(lambda x: 10**x)
			factor = ratios.div(ratios + 1)
			factor2 = 1 / (ratios + 1)

			lengths = totcole * factor

			lengths.hist(bins='fd', edgecolor='black', linewidth=1.2, ax=self.ax3)
		else:
			parameter_data.hist(bins='fd', edgecolor='black', linewidth=1.2, ax=self.ax3)

	def draw_variogram(self, variogram):
		variogram.plot(ax=self.ax4, legend=False)

	def draw(self, name=None):
		self.decorate_figures(name)
		self.canvas.draw()

	def reset_figures(self):
		self.ax1.clear()
		self.ax2.clear()
		self.ax3.clear()
		self.ax4.clear()

	def decorate_figures(self, par_name=None):
		self.fig.suptitle(par_name)

		self.ax1.set_title('Warm up')
		self.ax1.ticklabel_format(axis='both', style='sci', scilimits=(-2, 4))
		self.ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
		self.ax1.grid()

		self.ax2.set_title('Main run')
		self.ax2.ticklabel_format(axis='both', style='sci', scilimits=(-2, 4))
		self.ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
		self.ax2.grid()

		if self.parent.analysisinfo is not None and par_name is not None:
			self.ax3.set_title('Histogram | $n_{eff}$: %.2f | $\hat{R}$: %.3f'
		                   % (self.parent.analysisinfo.at['n_eff', par_name],
		                      self.parent.analysisinfo.at['R_hat', par_name]))
		else:
			self.ax3.set_title('Histogram')

		self.ax3.ticklabel_format(axis='both', style='sci', scilimits=(-2, 4))
		self.ax3.set_axisbelow(True)

		self.ax4.set_title('Variogram')
		self.ax4.ticklabel_format(axis='both', style='sci', scilimits=(-2, 4))
		self.ax4.grid()


class ChainListbox(tk.Frame):
	def __init__(self, master, lists):
		self.parent = master
		tk.Frame.__init__(self, self.parent)
		self.lists = []

		frame = tk.Frame(self)
		frame.pack(side=tk.RIGHT, fill=tk.Y)
		tk.Label(frame, borderwidth=1, relief=tk.RAISED).pack(fill=tk.X)
		sb = tk.Scrollbar(frame, orient=tk.VERTICAL, command=self._onVsb)
		sb.pack(expand=tk.YES, fill=tk.Y)

		for l, w in lists:
			frame = tk.Frame(self)
			frame.pack(side=tk.LEFT, expand=tk.YES, fill=tk.BOTH)
			tk.Label(frame, text=l, borderwidth=1, relief=tk.RAISED).pack(fill=tk.X)
			lb = tk.Listbox(frame, width=w, borderwidth=0, selectborderwidth=0,
			                relief=tk.FLAT, exportselection=tk.FALSE, yscrollcommand=sb.set)
			lb.pack(expand=tk.YES, fill=tk.BOTH)
			self.lists.append(lb)
			# lb.bind('<MouseWheel>', self._onMouseWheel)
			# lb.bind('<Button-1>', lambda e, s=self: s._select(e.y))
			# lb.bind('<B1-Motion>', lambda e, s=self: s._select(e.y))
			lb.bind('<MouseWheel>', self._onMouseWheel)
			lb.bind('<Button-1>', lambda e, s=self: s._highlight(e.y))
			lb.bind('<B1-Motion>', lambda e, s=self: s._highlight(e.y))
			lb.bind('<ButtonRelease-1>', lambda e, s=self: s._select(e.y))
			lb.bind('<Control-1>', lambda e, s=self: s._addselect(e.y))

	def _highlight(self, y):
		row = self.lists[0].nearest(y)
		if row >= 0:
			self.selection_clear(0, tk.END)
			self.selection_set(row)
		return 'break'

	def _select(self, y):
		return 'break'

	def _addselect(self, y):
		row = self.lists[0].nearest(y)
		self.selection_set(row)
		return 'break'

	def _button2(self, x, y):
		for l in self.lists: l.scan_mark(x, y)
		return 'break'

	def _b2motion(self, x, y):
		for l in self.lists: l.scan_dragto(x, y)
		return 'break'

	def _onVsb(self, *args):
		for l in self.lists:
			l.yview(*args)

	def _onMouseWheel(self, event):
		delta = -event.delta
		for l in self.lists:
			l.yview("scroll", delta, "units")
		return 'break'

	def curselection(self):
		return self.lists[0].curselection()

	def get(self, first, last=None):
		result = []
		for l in self.lists:
			result.append(l.get(first, last))
		if last:
			return map([None] + result)

		return result

	def index(self, index):
		return self.lists[0].index(index)

	def append(self, elements):
		if len(elements) == len(self.lists):
			for il, l in enumerate(self.lists):
				l.insert(tk.END, elements[il])

	def pop(self, last=None):
		for l in self.lists:
			l.delete(tk.END, last)

	def insert(self, first, elements):
		if len(elements) == len(self.lists):
			for il, l in enumerate(self.lists):
				l.insert(first, elements[il])

	def delete(self, first, last=None):
		for l in self.lists:
			l.delete(first, last)

	def update_entries(self, row, update):
		print(self.lists)

	def size(self):
		return self.lists[0].size()

	def see(self, index):
		for l in self.lists:
			l.see(index)

	def selection_anchor(self, index):
		for l in self.lists:
			l.selection_anchor(index)

	def selection_clear(self, first, last=None):
		for l in self.lists:
			l.selection_clear(first, last)

	def selection_includes(self, index):
		return self.lists[0].selection_includes(index)

	def selection_set(self, first, last=None):
		for l in self.lists:
			l.selection_set(first, last)

	def delete_all(self):
		for l in self.lists:
			l.delete(0, tk.END)

	def select_all(self):
		for l in self.lists:
			l.selection_set(0, tk.END)


class AnalysisListbox(tk.Frame):
	def __init__(self, master, lists):
		self.parent = master
		tk.Frame.__init__(self, self.parent)
		self.lists = []

		frame = tk.Frame(self)
		frame.pack(side=tk.RIGHT, fill=tk.Y)
		tk.Label(frame, borderwidth=1, relief=tk.RAISED).pack(fill=tk.X)
		sb = tk.Scrollbar(frame, orient=tk.VERTICAL, command=self._onVsb)
		sb.pack(expand=tk.YES, fill=tk.Y)

		for l, w in lists:
			frame = tk.Frame(self)
			frame.pack(side=tk.LEFT, expand=tk.YES, fill=tk.BOTH)
			tk.Label(frame, text=l, borderwidth=1, relief=tk.RAISED).pack(fill=tk.X)
			lb = tk.Listbox(frame, width=w, borderwidth=0, selectborderwidth=0,
			                relief=tk.FLAT, exportselection=tk.FALSE, yscrollcommand=sb.set)
			lb.pack(expand=tk.YES, fill=tk.BOTH)
			self.lists.append(lb)
			# lb.bind('<MouseWheel>', self._onMouseWheel)
			# lb.bind('<Button-1>', lambda e, s=self: s._select(e.y))
			# lb.bind('<B1-Motion>', lambda e, s=self: s._select(e.y))
			lb.bind('<MouseWheel>', self._onMouseWheel)
			lb.bind('<Button-1>', lambda e, s=self: s._highlight(e.y))
			lb.bind('<B1-Motion>', lambda e, s=self: s._highlight(e.y))
			lb.bind('<ButtonRelease-1>', lambda e, s=self: s._select(e.y))
			#lb.bind('<Control-1>', lambda e, s=self: s._addselect(e.y))

	def _highlight(self, y):
		row = self.lists[0].nearest(y)
		if row >= 0:
			self.selection_clear(0, tk.END)
			self.selection_set(row)
		return 'break'

	def _select(self, y):
		row = self.lists[0].nearest(y)
		if row >= 0:
			self.parent.draw_parameter(self.get(row)[0])
			pass
		return 'break'

	def _addselect(self, y):
		row = self.lists[0].nearest(y)
		self.selection_set(row)
		return 'break'

	def _button2(self, x, y):
		for l in self.lists: l.scan_mark(x, y)
		return 'break'

	def _b2motion(self, x, y):
		for l in self.lists: l.scan_dragto(x, y)
		return 'break'

	def _onVsb(self, *args):
		for l in self.lists:
			l.yview(*args)

	def _onMouseWheel(self, event):
		delta = -event.delta
		for l in self.lists:
			l.yview("scroll", delta, "units")
		return 'break'

	def curselection(self):
		return self.lists[0].curselection()

	def delete(self, first, last=None):
		for l in self.lists:
			l.delete(first, last)

	def get(self, first, last=None):
		result = []
		for l in self.lists:
			result.append(l.get(first, last))
		if last: return map([None] + result)
		return result

	def index(self, index):
		return self.lists[0].index(index)

	def insert_df(self, df):
		mylist = [[name] + element.to_list() for name, element in df.iteritems()]
		fmtlst = ['{:s}', '{:.3f}', '{:.2f}']

		for e in mylist:
			i = 0
			for l in self.lists:
				l.insert(tk.END, fmtlst[i].format(e[i]))
				i = i + 1

	def size(self):
		return self.lists[0].size()

	def see(self, index):
		for l in self.lists:
			l.see(index)

	def selection_anchor(self, index):
		for l in self.lists:
			l.selection_anchor(index)

	def selection_clear(self, first, last=None):
		for l in self.lists:
			l.selection_clear(first, last)

	def selection_includes(self, index):
		return self.lists[0].selection_includes(index)

	def selection_set(self, first, last=None):
		for l in self.lists:
			l.selection_set(first, last)

	def selection_reset(self, first, last=None):
		if first:
			for l in self.lists:
				l.selection_set(first, last)

			self.parent.draw_parameter(self.get(first)[0])

	def delete_all(self):
		for l in self.lists:
			l.delete(0, tk.END)

	def select_all(self):
		for l in self.lists:
			l.selection_set(0, tk.END)


class SingleListbox(tk.Frame):
	def __init__(self, master, lists):
		self.parent = master
		tk.Frame.__init__(self, self.parent)
		self.lists = []

		frame = tk.Frame(self)
		frame.pack(side=tk.RIGHT, fill=tk.Y)
		tk.Label(frame, borderwidth=1, relief=tk.RAISED).pack(fill=tk.X)
		sb = tk.Scrollbar(frame, orient=tk.VERTICAL, command=self._onVsb)
		sb.pack(expand=tk.YES, fill=tk.Y)

		frame = tk.Frame(self)
		frame.pack(side=tk.LEFT, expand=tk.YES, fill=tk.BOTH)
		tk.Label(frame, text=lists[0], borderwidth=1, relief=tk.RAISED).pack(fill=tk.X)
		lb = tk.Listbox(frame, width=lists[1], borderwidth=0, selectborderwidth=0, selectmode=tk.EXTENDED,
		                relief=tk.FLAT, exportselection=tk.FALSE, yscrollcommand=sb.set)
		lb.pack(expand=tk.YES, fill=tk.BOTH)
		self.lists.append(lb)
		lb.bind('<MouseWheel>', self._onMouseWheel)
		lb.bind('<Button-1>', lambda e, s=self: s._highlight(e.y))
		lb.bind('<B1-Motion>', lambda e, s=self: s._highlight(e.y))
		lb.bind('<ButtonRelease-1>', lambda e, s=self: s._select(e.y))
		lb.bind('<Control-1>', lambda e, s=self: s._addselect(e.y))

	def _highlight(self, y):
		row = self.lists[0].nearest(y)
		if row >= 0:
			self.selection_clear(0, tk.END)
			self.selection_set(row)
		return 'break'

	def _select(self, y):
		row = self.lists[0].nearest(y)
		if row >= 0:
			self.parent.draw_parameter(self.get(row)[0])
		return 'break'

	def _addselect(self, y):
		row = self.lists[0].nearest(y)
		self.selection_set(row)
		return 'break'

	def _button2(self, x, y):
		for l in self.lists: l.scan_mark(x, y)
		return 'break'

	def _b2motion(self, x, y):
		for l in self.lists: l.scan_dragto(x, y)
		return 'break'

	def _onVsb(self, *args):
		for l in self.lists:
			l.yview(*args)

	def _onMouseWheel(self, event):
		delta = -event.delta
		for l in self.lists:
			l.yview("scroll", delta, "units")
		return 'break'

	def curselection(self):
		return self.lists[0].curselection()

	def delete(self, first, last=None):
		for l in self.lists:
			l.delete(first, last)

	def get(self, first, last=None):
		result = []
		for l in self.lists:
			result.append(l.get(first, last))
		if last: return map([None] + result)
		return result

	def index(self, index):
		return self.lists[0].index(index)

	def insert(self, elements):
		for e in elements:
			for l in self.lists:
				l.insert(tk.END, e)

	def size(self):
		return self.lists[0].size()

	def see(self, index):
		for l in self.lists:
			l.see(index)

	def selection_anchor(self, index):
		for l in self.lists:
			l.selection_anchor(index)

	def selection_clear(self, first, last=None):
		for l in self.lists:
			l.selection_clear(first, last)

	def selection_includes(self, index):
		return self.lists[0].selection_includes(index)

	def selection_set(self, first, last=None):
		for l in self.lists:
			l.selection_set(first, last)

	def delete_all(self):
		for l in self.lists:
			l.delete(0, tk.END)

	def select_all(self):
		self.selection_set(0, tk.END)
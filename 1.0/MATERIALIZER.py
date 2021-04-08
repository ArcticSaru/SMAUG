import tkinter as tk
import tkinter.ttk as ttk

from tkinter import filedialog
from tkinter import messagebox
from tkinter import font
import math
import numpy as np
import scipy.stats as sps
import scipy.optimize as spo
import time
import pandas as pd
import pandas.errors as pde
import json
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.gridspec
from io import StringIO
import os
from myListboxes import MaterialListbox
from Material import *
import threading


class DensityPanelInitial(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)
		self.__createlayout()

	def __createlayout(self):
		self.denlabel = tk.Label(self, text='Density')

		self.denlabel.grid(row=0, sticky='nw')
		self.rowconfigure(0, weight=1)
		self.columnconfigure(0, weight=1)

	def saveDensity(self, spath):
		pass


class DensityPanelFixed(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)

		self.meanvar = tk.DoubleVar()

		self.__createlayout()

	def __createlayout(self):
		self.denlabel = tk.Label(self, text='Density')

		self.meanlbl = tk.Label(self, text='Value:')
		self.meanety = tk.Entry(self, textvariable=self.meanvar)

		self.meanlbl.grid(row=1, column=0)
		self.meanety.grid(row=1, column=1)

		self.columnconfigure(0, weight=1)
		self.rowconfigure(1, weight=1)

	def saveDensity(self, spath):
		dendict = {'value': self.meanvar.get()}

		with open(os.path.normpath(os.path.join(spath, 'density.json')), 'w') as fp:
			json.dump(dendict, fp)



class DensityPanelNormal(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)

		self.meanvar = tk.DoubleVar()
		self.stdvar = tk.DoubleVar()

		self.__createlayout()

	def __createlayout(self):
		self.denlabel = tk.Label(self, text='Density')

		self.meanlbl = tk.Label(self, text='Mean:')
		self.meanety = tk.Entry(self, textvariable=self.meanvar)

		self.stdlbl = tk.Label(self, text='Std-dev:')
		self.stdety = tk.Entry(self, textvariable=self.stdvar)

		self.showbtn = tk.Button(self, text='Show', command=self.showplot)

		self.denlabel.grid(row=0, sticky='nw')
		self.meanlbl.grid(row=1, column=0)
		self.meanety.grid(row=1, column=1)
		self.stdlbl.grid(row=1, column=2)
		self.stdety.grid(row=1, column=3)
		self.showbtn.grid(row=2, column=2, sticky='ew')

		self.columnconfigure(0, weight=1)
		self.rowconfigure(1, weight=1)

	def saveDensity(self, spath):
		dendict = {'mean': self.meanvar.get()}

		if self.stdvar.get() < 0.0:
			messagebox.showerror('Std must be greater than 0, please change!')
			return

		dendict.update({'std': self.stdvar.get()})

		with open(os.path.normpath(os.path.join(spath, 'density.json')), 'w') as fp:
			json.dump(dendict, fp)

	def showplot(self):
		# print(self.meanvar.get())
		# print(self.stdvar.get())
		self.parent.showdensityplot(self.meanvar.get(), self.stdvar.get())

class DensityPanelLogNormal(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)

		self.meanvar = tk.DoubleVar()
		self.stdvar = tk.DoubleVar()

		self.__createlayout()

	def __createlayout(self):
		self.denlabel = tk.Label(self, text='Density')

		self.loglabel = tk.Label(self, text='Mean and standard dev of nonlog-data:')

		self.meanlbl = tk.Label(self, text='mean:')
		self.meanety = tk.Entry(self, textvariable=self.meanvar)

		self.stdlbl = tk.Label(self, text='std:')
		self.stdety = tk.Entry(self, textvariable=self.stdvar)

		self.showbtn = tk.Button(self, text='Show', command=self.showplot)

		self.denlabel.grid(row=0, columnspan=3, sticky='nw')
		self.loglabel.grid(row=1, columnspan=3)
		self.meanlbl.grid(row=2, column=0)
		self.meanety.grid(row=2, column=1)
		self.stdlbl.grid(row=2, column=2)
		self.stdety.grid(row=2, column=3)
		self.showbtn.grid(row=3, column=2, sticky='ew')

		self.columnconfigure(0, weight=1)
		self.rowconfigure(2, weight=1)

	def saveDensity(self, spath):
		dendict = {}

		if self.stdvar.get() > 0.0:
			truevar = np.log(1. + (self.stdvar.get() / self.meanvar.get()) ** 2)
			truesig = truevar ** 0.5
			truemu = np.log(self.meanvar.get()) - 0.5 * truevar
			dendict.update({'scale': np.exp(truemu)})
			dendict.update({'shape': truesig})
		else:
			messagebox.showerror('Std must be greater than 0, please change!')
			return

		with open(os.path.normpath(os.path.join(spath, 'density.json')), 'w') as fp:
			json.dump(dendict, fp)

	def showplot(self):
		self.parent.showdensityplot(self.meanvar.get(), self.stdvar.get())


class DensityPanelData(tk.Frame):
	modes = ['Bulk', 'Skeletal', 'Grain']

	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)

		self.meancheckvar = tk.IntVar()
		self.meancheckvar.set(1)
		self.covcheckvar = tk.IntVar()
		self.covcheckvar.set(0)

		self.meanvar = 0.0
		self.stdvar = 0.0

		self.selectvar = tk.IntVar()
		self.modevar = tk.StringVar()

		self.modevar.set('Choose')

		self.__createlayout()

	def __createlayout(self):
		self.denlabel = tk.Label(self, text='Density')

		self.loaddatabutton = tk.Button(self, text='Load data', command=self.parent.loadDensityData)
		self.meancheckbox = tk.Checkbutton(self, text='Include mean', var=self.meancheckvar, state=tk.DISABLED)
		self.covcheckbox = tk.Checkbutton(self, text='Include Cov-Matrix', var=self.covcheckvar)

		self.rawcheckbox = tk.Radiobutton(self, text='Load raw data', variable=self.selectvar, value=0)
		self.sampcheckbox = tk.Radiobutton(self, text='Load sample data', variable=self.selectvar, value=1)

		self.modeselect = tk.OptionMenu(self, self.modevar, *self.modes)
		self.modeselect.config(bg ='gray97', relief=tk.GROOVE)

		self.denlabel.grid(row=0, column=0, sticky='nw')
		self.loaddatabutton.grid(row=0, column=1)
		self.meancheckbox.grid(row=1, column=2, sticky='e')
		self.covcheckbox.grid(row=2, column=2, sticky='e')
		self.rawcheckbox.grid(row=1, column=0, sticky='w')
		self.sampcheckbox.grid(row=2, column=0, sticky='w')
		self.modeselect.grid(row=2, column=1)

		self.rowconfigure(0, weight=1)
		self.columnconfigure(1, weight=1)

	def saveDensity(self, spath):
		dendict = {'scale': self.meanvar}

		if self.covcheckvar.get() == 1:
			dendict.update({'shape': self.stdvar})

		with open(os.path.normpath(os.path.join(spath, 'density.json')), 'w') as fp:
			json.dump(dendict, fp)


class MaterialPanel(tk.Frame):
	DenOptionList = ['Fixed', 'From data', 'Normal', 'Log-normal']
	#CompOptionList = ['Oxides', 'Formula', 'Groom mixture', 'Groom element']
	CompOptionList = ['Oxides', 'Groom mixture', 'Groom element']

	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)
		self.__createlayout()

	def __createlayout(self):
		self.denvar = tk.StringVar()
		self.compvar = tk.StringVar()

		self.denvar.set('Choose')
		self.compvar.set('Choose')

		self.compvar.trace('w', self.compositionModeChange)
		self.denvar.trace('w', self.densityModeChange)

		self.matlabel = tk.Label(self, text='Material')
		self.denlabel = tk.Label(self, text='Density mode')
		self.complabel = tk.Label(self, text='Composition mode')
		self.densitybox = tk.OptionMenu(self, self.denvar, *self.DenOptionList)
		self.compobox = tk.OptionMenu(self, self.compvar, *self.CompOptionList)
		self.saveButton = tk.Button(self, text='Save Material', command=self.parent.savematerial)
		self.groomButton = tk.Button(self, text='Create Groom Table', command=self.parent.make_groom_table)

		self.matlabel.grid(column=1, row=0, sticky='n', padx=5, pady=5)
		self.densitybox.grid(column=0, row=1, sticky='w', padx=5, pady=5)
		self.compobox.grid(column=2, row=1, sticky='e', padx=5, pady=5)
		self.denlabel.grid(column=0, row=0, sticky='w', padx=5, pady=5)
		self.complabel.grid(column=2, row=0, sticky='e', padx=5, pady=5)
		self.saveButton.grid(column=0, row=2, padx=5, pady=5)
		self.groomButton.grid(column=2, row=2, padx=5, pady=5)

		self.rowconfigure(0, weight=1)
		self.columnconfigure(0, weight=1)
		self.columnconfigure(1, weight=1)
		self.columnconfigure(2, weight=1)

	def compositionModeChange(self, *args):
		mode = self.compvar.get()

		self.parent.compositionModeChange(mode)

	def densityModeChange(self, *args):
		mode = self.denvar.get()

		self.parent.densityModeChange(mode)

	def getElementInfo(self, name):
		return self.parent.getElementInfo(name)


class CompositionPanelInitial(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)
		self.__createlayout()

	def __createlayout(self):
		self.compolabel = tk.Label(self, text='Composition')

		self.compolabel.grid(row=0, sticky='ne')

		self.rowconfigure(0, weight=1)
		self.columnconfigure(0, weight=1)

	def saveComposition(self, spath):
		pass

	def make_groom_table(self, spath):
		pass

class CompositionPanelGroEle(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)
		self.__createlayout()

	def __createlayout(self):
		self.compolabel = tk.Label(self, text='Composition')
		self.componame = tk.Label(self, text='')
		self.formulalbl = tk.Label(self, text='Formula:')
		self.formulaety = tk.Entry(self, text='', state='disabled')
		self.zlbl = tk.Label(self, text='Z')
		self.zety = tk.Entry(self, text='', state='disabled')
		self.albl = tk.Label(self, text='A')
		self.aety = tk.Entry(self, text='', state='disabled')
		self.rholbl = tk.Label(self, text='Density')
		self.rhoety = tk.Entry(self, text='', state='disabled')
		self.Ilbl = tk.Label(self, text='I_eff')
		self.Iety = tk.Entry(self, text='', state='disabled')
		self.Cbarlbl = tk.Label(self, text='Cbar')
		self.Cbarety = tk.Entry(self, text='', state='disabled')
		self.x0lbl = tk.Label(self, text='x0')
		self.x0ety = tk.Entry(self, text='', state='disabled')
		self.x1lbl = tk.Label(self, text='x1')
		self.x1ety = tk.Entry(self, text='', state='disabled')
		self.aalbl = tk.Label(self, text='a')
		self.aaety = tk.Entry(self, text='', state='disabled')
		self.kslbl = tk.Label(self, text='k')
		self.ksety = tk.Entry(self, text='', state='disabled')
		self.d0lbl = tk.Label(self, text='d0')
		self.d0ety = tk.Entry(self, text='', state='disabled')

		self.compolabel.grid(row=0, column=7, sticky='ne')
		self.componame.grid(row=0, column=0, sticky='nw')
		self.formulalbl.grid(row=1, column=0)
		self.formulaety.grid(row=1, column=1)
		self.zlbl.grid(row=1, column=2)
		self.zety.grid(row=1, column=3)
		self.albl.grid(row=1, column=4)
		self.aety.grid(row=1, column=5)
		self.rholbl.grid(row=1, column=6)
		self.rhoety.grid(row=1, column=7)
		self.Ilbl.grid(row=2, column=0)
		self.Iety.grid(row=2, column=1)
		self.Cbarlbl.grid(row=2, column=2)
		self.Cbarety.grid(row=2, column=3)
		self.x0lbl.grid(row=2, column=4)
		self.x0ety.grid(row=2, column=5)
		self.x1lbl.grid(row=2, column=6)
		self.x1ety.grid(row=2, column=7)
		self.aalbl.grid(row=3, column=0)
		self.aaety.grid(row=3, column=1)
		self.kslbl.grid(row=3, column=2)
		self.ksety.grid(row=3, column=3)
		self.d0lbl.grid(row=3, column=4)
		self.d0ety.grid(row=3, column=5)

		self.rowconfigure(0, weight=1)
		for ci in range(0, 6):
			self.columnconfigure(ci, weight=1)

	def setActive(self):
		self.formulaety.config(state='normal')
		self.zety.config(state='normal')
		self.aety.config(state='normal')
		self.rhoety.config(state='normal')
		self.Iety.config(state='normal')
		self.Cbarety.config(state='normal')
		self.x0ety.config(state='normal')
		self.x1ety.config(state='normal')
		self.aaety.config(state='normal')
		self.ksety.config(state='normal')
		self.d0ety.config(state='normal')

	def setDisabled(self):
		self.formulaety.config(state='disabled')
		self.zety.config(state='disabled')
		self.aety.config(state='disabled')
		self.rhoety.config(state='disabled')
		self.Iety.config(state='disabled')
		self.Cbarety.config(state='disabled')
		self.x0ety.config(state='disabled')
		self.x1ety.config(state='disabled')
		self.aaety.config(state='disabled')
		self.ksety.config(state='disabled')
		self.d0ety.config(state='disabled')

	def saveComposition(self, spath):
		name = self.componame.cget('text')

		element_row = self.parent.getElementInfo(name)
		element_info = element_row.squeeze().to_dict()

		print(element_info)

		with open(os.path.normpath(os.path.join(spath, 'composition.json')), 'w') as fp:
			json.dump(element_info, fp)

		ele_obj = ElementGUI(element_info['Formula'], element_info)
		ele_obj.create_radloss_table_gui(spath)

	def make_groom_table(self, spath):
		name = self.componame.cget('text')

		element_row = self.parent.getElementInfo(name)
		element_info = element_row.squeeze().to_dict()

		ele_obj = ElementGUI(element_info['Formula'], element_info)
		ele_obj.make_groom_table(spath)


class CompositionPanelGroMix(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)
		self.__createlayout()

	def __createlayout(self):
		self.compolabel = tk.Label(self, text='Composition')
		self.componame = tk.Label(self, text='')
		self.formulalbl = tk.Label(self, text='Formula:')
		self.formulaety = tk.Entry(self, text='', state='disabled')
		self.zalbl = tk.Label(self, text='Z/A')
		self.zaety = tk.Entry(self, text='', state='disabled')
		self.albl = tk.Label(self, text='A')
		self.aety = tk.Entry(self, text='', state='disabled')
		self.rholbl = tk.Label(self, text='Density')
		self.rhoety = tk.Entry(self, text='', state='disabled')
		self.Ilbl = tk.Label(self, text='I_eff')
		self.Iety = tk.Entry(self, text='', state='disabled')
		self.Cbarlbl = tk.Label(self, text='Cbar')
		self.Cbarety = tk.Entry(self, text='', state='disabled')
		self.x0lbl = tk.Label(self, text='x0')
		self.x0ety = tk.Entry(self, text='', state='disabled')
		self.x1lbl = tk.Label(self, text='x1')
		self.x1ety = tk.Entry(self, text='', state='disabled')
		self.aalbl = tk.Label(self, text='a')
		self.aaety = tk.Entry(self, text='', state='disabled')
		self.kslbl = tk.Label(self, text='k')
		self.ksety = tk.Entry(self, text='', state='disabled')
		self.d0lbl = tk.Label(self, text='d0')
		self.d0ety = tk.Entry(self, text='', state='disabled')

		self.compolabel.grid(row=0, column=7, sticky='ne')
		self.componame.grid(row=0, column=0, sticky='nw')
		self.formulalbl.grid(row=1, column=0)
		self.formulaety.grid(row=1, column=1)
		self.zalbl.grid(row=1, column=2)
		self.zaety.grid(row=1, column=3)
		self.albl.grid(row=1, column=4)
		self.aety.grid(row=1, column=5)
		self.rholbl.grid(row=1, column=6)
		self.rhoety.grid(row=1, column=7)
		self.Ilbl.grid(row=2, column=0)
		self.Iety.grid(row=2, column=1)
		self.Cbarlbl.grid(row=2, column=2)
		self.Cbarety.grid(row=2, column=3)
		self.x0lbl.grid(row=2, column=4)
		self.x0ety.grid(row=2, column=5)
		self.x1lbl.grid(row=2, column=6)
		self.x1ety.grid(row=2, column=7)
		self.aalbl.grid(row=3, column=0)
		self.aaety.grid(row=3, column=1)
		self.kslbl.grid(row=3, column=2)
		self.ksety.grid(row=3, column=3)
		self.d0lbl.grid(row=3, column=4)
		self.d0ety.grid(row=3, column=5)

		self.rowconfigure(0, weight=1)
		for ci in range(0, 6):
			self.columnconfigure(ci, weight=1)

	def setActive(self):
		self.formulaety.config(state='normal')
		self.zaety.config(state='normal')
		self.aety.config(state='normal')
		self.rhoety.config(state='normal')
		self.Iety.config(state='normal')
		self.Cbarety.config(state='normal')
		self.x0ety.config(state='normal')
		self.x1ety.config(state='normal')
		self.aaety.config(state='normal')
		self.ksety.config(state='normal')
		self.d0ety.config(state='normal')

	def setDisabled(self):
		self.formulaety.config(state='disabled')
		self.zaety.config(state='disabled')
		self.aety.config(state='disabled')
		self.rhoety.config(state='disabled')
		self.Iety.config(state='disabled')
		self.Cbarety.config(state='disabled')
		self.x0ety.config(state='disabled')
		self.x1ety.config(state='disabled')
		self.aaety.config(state='disabled')
		self.ksety.config(state='disabled')
		self.d0ety.config(state='disabled')

	def saveComposition(self, spath):
		name = self.componame.cget('text')

		mixture_row = self.parent.getElementInfo(name)
		mixture_info = mixture_row.squeeze().to_dict()

		groom_df = self.parent.getGroomDF()

		with open(os.path.normpath(os.path.join(spath, 'composition.json')), 'w') as fp:
			json.dump(mixture_info, fp)

		mix_obj = MixtureGUI(mixture_info['Name'], groom_df)
		mix_obj.create_radloss_table_gui(spath)

	def make_groom_table(self, spath):
		name = self.componame.cget('text')

		mix_row = self.parent.getElementInfo(name)
		# element_info = element_row.squeeze().to_dict()

		mix_obj = MixtureGUI(name, self.parent.getGroomDF())
		mix_obj.make_groom_table(spath)

class CompositionPanelOxides(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)

		self.meanvar = tk.IntVar()
		self.meanvar.set(1)
		self.covvar = tk.IntVar()
		self.covvar.set(0)

		self.__createlayout()

	def __createlayout(self):
		self.compolabel = tk.Label(self, text='Composition')
		self.table = ttk.Treeview(self)
		self.loaddatabutton = tk.Button(self, text='Load data', command=self.parent.loadOxideData)
		self.meancheckbox = tk.Checkbutton(self, text='Include mean', var=self.meanvar, state=tk.DISABLED)
		self.covcheckbox = tk.Checkbutton(self, text='Include Cov-Matrix', var=self.covvar)

		self.compolabel.grid(row=0, sticky='ne')
		self.loaddatabutton.grid(row=0, sticky='nw')
		self.meancheckbox.grid(row=1, column=0, sticky='w')
		self.covcheckbox.grid(row=2, column=0, sticky='w')

		self.rowconfigure(0, weight=1)
		self.columnconfigure(0, weight=1)

	def saveComposition(self, spath):
		# TODO: implement
		pass


class CompositionPanelFormula(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)
		self.__createlayout()

	def __createlayout(self):
		self.compolabel = tk.Label(self, text='Composition')

		self.compolabel.grid(row=0, sticky='ne')

		self.rowconfigure(0, weight=1)
		self.columnconfigure(0, weight=1)

	def saveComposition(self, spath):
		# TODO: implement
		pass

	def make_groom_table(self, spath):
		pass


class ControlWindow(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)
		self.frames = {}
		self.denframes = {}
		self.__createlayout()

	def __createlayout(self):
		for DF in (DensityPanelInitial, DensityPanelFixed, DensityPanelNormal, DensityPanelLogNormal, DensityPanelData):
			panel_name = DF
			frame = DF(self)
			self.denframes[panel_name] = frame

		self.controlframe = MaterialPanel(self)

		for F in (CompositionPanelInitial, CompositionPanelGroEle, CompositionPanelGroMix, CompositionPanelOxides,
				  CompositionPanelFormula):
			panel_name = F
			frame = F(self)
			self.frames[panel_name] = frame

		self.sep1 = ttk.Separator(self, orient=tk.VERTICAL)
		self.sep2 = ttk.Separator(self, orient=tk.VERTICAL)

		for frame in self.denframes:
			self.denframes[frame].grid(column=0, row=0, padx=2, pady=2, sticky='nsew')

		self.controlframe.grid(column=2, row=0, padx=2, pady=2, sticky='nsew')

		for frame in self.frames:
			self.frames[frame].grid(column=4, row=0, padx=2, pady=2, sticky='nsew')

		self.sep1.grid(column=1, row=0, padx=2, pady=2, sticky='ns')
		self.sep2.grid(column=3, row=0, padx=2, pady=2, sticky='ns')

		self.columnconfigure(0, weight=1)
		self.columnconfigure(2, weight=1)
		self.columnconfigure(4, weight=1)

	def show_compo_frame(self, panel_name):
		frame = self.frames[panel_name]
		frame.tkraise()

	def show_den_frame(self, panel_name):
		frame = self.denframes[panel_name]
		frame.tkraise()

	def compositionModeChange(self, mode):
		self.parent.compositionModeChange(mode)

		if mode == 'Groom element':
			self.show_compo_frame(CompositionPanelGroEle)
		elif mode == 'Groom mixture':
			self.show_compo_frame(CompositionPanelGroMix)
		elif mode == 'Oxides':
			self.show_compo_frame(CompositionPanelOxides)
		elif mode == 'Formula':
			self.show_compo_frame(CompositionPanelFormula)
		else:
			self.show_compo_frame(CompositionPanelInitial)

	def densityModeChange(self, mode):
		self.parent.densityModeChange(mode)

		if mode == 'Normal':
			self.show_den_frame(DensityPanelNormal)
		elif mode == 'Fixed':
			self.show_den_frame(DensityPanelFixed)
		elif mode == 'Log-normal':
			self.show_den_frame(DensityPanelLogNormal)
		elif mode == 'From data':
			self.show_den_frame(DensityPanelData)

	def loadOxideData(self):
		self.parent.loadOxideData()

	def showdensityplot(self, mu, std):
		self.parent.showdensityplot(mu, std)

	def loadDensityData(self):
		self.parent.loadDensityData()

	def savematerial(self):
		if self.controlframe.denvar.get() != 'Choose' and self.controlframe.compvar.get() != 'Choose':
			denmode = self.parent.getActDenMode()
			compmode = self.parent.getActCompMode()

			savepath = tk.filedialog.askdirectory(initialdir='./Materials')

			if savepath:
				infodict = {}

				dclass = DensityPanelInitial
				if denmode == 'Normal':
					dclass = DensityPanelNormal
				elif denmode == 'Log-normal':
					dclass = DensityPanelLogNormal
				elif denmode == 'From data':
					dclass = DensityPanelData


				self.denframes[dclass].saveDensity(savepath)


				cclass = CompositionPanelInitial
				if compmode == 'Oxides':
					cclass = CompositionPanelOxides
				elif compmode == 'Formula':
					cclass = CompositionPanelFormula
				elif compmode == 'Groom mixture':
					cclass = CompositionPanelGroMix
				elif compmode == 'Groom element':
					cclass = CompositionPanelGroEle

				self.frames[cclass].saveComposition(savepath)

				infodict.update({'density': denmode, 'composition': compmode})

				with open(os.path.normpath(os.path.join(savepath,'info.json')), 'w') as fp:
					json.dump(infodict, fp)

		else:
			messagebox.showerror('Error', 'Choose both, Density AND Composition')

	def getElementInfo(self, name):
		return self.parent.getElementInfo(name)

	def getGroomDF(self):
		return self.parent.getGroomDF()

	def make_groom_table(self):
		if self.controlframe.compvar.get() != 'Choose':
			compmode = self.parent.getActCompMode()
			savepath = tk.filedialog.askdirectory(initialdir='./Materials')

			if savepath:
				cclass = CompositionPanelInitial
				if compmode == 'Oxides':
					cclass = CompositionPanelOxides
				elif compmode == 'Formula':
					cclass = CompositionPanelFormula
				elif compmode == 'Groom mixture':
					cclass = CompositionPanelGroMix
				elif compmode == 'Groom element':
					cclass = CompositionPanelGroEle

				self.frames[cclass].make_groom_table(savepath)



class ListboxWindow(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)
		self.__createlayout()

	def __createlayout(self):
		self.matlistbox = MaterialListbox(self, (('Pd_idx', 10), ('Name', 50)))

		self.matlistbox.grid(column=0, row=0, padx=2, pady=2, sticky='ns')

		self.rowconfigure(0, weight=1)

	def changeCompoInfo(self, df_index):
		self.parent.changeCompoInfo(df_index)


class VisualisationWindow(tk.Frame):

	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)
		self.__createlayout()

	def __createlayout(self):
		self.fig = plt.figure(constrained_layout=True)

		self.drawlabel = tk.Label(self, text='Visualisation')
		self.r_progbar = ttk.Progressbar()
		self.drawArea = FigureCanvasTkAgg(self.fig, self)

		self.drawlabel.grid(column=0, row=0, sticky='w', padx=2, pady=2)
		self.drawArea.get_tk_widget().grid(column=0, row=1, sticky='nsew', columnspan=2)

		self.rowconfigure(1, weight=1)
		self.columnconfigure(0, weight=1)


class MainWindow(tk.Frame):

	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent)
		self.initUI()
		self.readGroomData()

		self.actCompoMode = ''
		self.actDenMode = ''

	def initUI(self):
		self.bottomframe = ControlWindow(self)
		self.visframeleft = VisualisationWindow(self)
		self.visframeright = VisualisationWindow(self)
		self.rightframe = ListboxWindow(self)

		self.visframeleft.grid(column=0, row=0, padx=2, pady=2, sticky='nsew')
		self.visframeright.grid(column=1, row=0, padx=2, pady=2, sticky='nsew')
		self.rightframe.grid(column=2, row=0, padx=2, pady=2, sticky='ns')
		self.bottomframe.grid(column=0, row=1, columnspan=3, padx=2, pady=2, sticky='ew')

		self.rowconfigure(0, weight=1)
		self.columnconfigure(0, weight=1)
		self.columnconfigure(1, weight=1)

		self.bottomframe.show_compo_frame(CompositionPanelInitial)
		self.bottomframe.show_den_frame(DensityPanelInitial)

	def getActDenMode(self):
		return self.actDenMode

	def getActCompMode(self):
		return self.actCompoMode

	def readGroomData(self):
		exepath = os.path.split((os.path.abspath(__file__)))[0]
		filepath = os.path.normpath(os.path.join(exepath, 'Input_Material', 'Groom', 'properties.txt'))

		row_list = []
		with open(filepath) as file:
			eleidx = 0
			elelist = []
			matdict = {}
			for idl, line in enumerate(file):
				eleidx += 1

				# print(idl, eleidx, ' :', line)
				spli = line.split()

				if eleidx == 1:
					A = float(spli[3])
					rho = float(spli[-6])
					nele = int(spli[-4])
					nlines = int(spli[-2]) + nele + 4
					thing = spli[-1]
					zovera = float(spli[-8])
				elif eleidx == 2:
					if thing == 'E' or thing == 'R':
						formula = spli[0]
						name = spli[1]
					else:
						formula = ''
						name = spli[0]
				elif eleidx == 4:
					I = float(spli[0])
					Cbar = float(spli[1])
					x0 = float(spli[2])
					x1 = float(spli[3])
					aa = float(spli[4])
					sk = float(spli[5])
					dlt0 = float(spli[6])
				elif 5 <= eleidx < nele + 5:
					elelist.append((int(spli[0]), float(spli[2])))
				elif eleidx > nlines:
					if thing == 'E' or thing == 'R':
						matdict.update(
							{
								'Name': name,
								'Formula': formula,
								'Z': elelist[0][0],
								'A': A,
								'Density': rho,
								'I': I,
								'Cbar': Cbar,
								'x0': x0,
								'x1': x1,
								'a': aa,
								'k': sk,
								'dlt0': dlt0,
								'Type': thing
							}
						)
					else:
						matdict.update(
							{
								'Name': name,
								'Z/A': zovera,
								'A': A,
								'Density': rho,
								'I': I,
								'Cbar': Cbar,
								'x0': x0,
								'x1': x1,
								'a': aa,
								'k': sk,
								'dlt0': dlt0,
								'Constituents': elelist,
								'Type': thing
							}
						)
					row_list.append(matdict)
					matdict = {}
					eleidx = 0
					elelist = []

		self.groom_df = pd.DataFrame(row_list)

	def compositionModeChange(self, mode):

		self.actCompoMode = mode

		if mode == 'Groom element':
			self.rightframe.matlistbox.delete_all()
			self.visframeright.fig.clear()
			self.visframeright.drawArea.draw()

			element = self.groom_df['Type'] == 'E'
			radioactive = self.groom_df['Type'] == 'R'

			for index, row in self.groom_df[element | radioactive].iterrows():
				self.rightframe.matlistbox.insert(tk.END, [index, row['Name']])

			self.rightframe.matlistbox.selection_set(0)
		elif mode == 'Groom mixture':
			self.rightframe.matlistbox.delete_all()
			self.visframeright.fig.clear()
			self.visframeright.drawArea.draw()

			notelement = self.groom_df['Type'] != 'E'
			notradioactive = self.groom_df['Type'] != 'R'

			for index, row in self.groom_df[notelement & notradioactive].iterrows():
				self.rightframe.matlistbox.insert(tk.END, [index, row['Name']])

			self.rightframe.matlistbox.selection_set(0)
		elif mode == 'Oxides':
			self.rightframe.matlistbox.delete_all()
			self.visframeright.fig.clear()
			self.visframeright.drawArea.draw()
		elif mode == 'Formula':
			self.rightframe.matlistbox.delete_all()
			self.visframeright.fig.clear()
			self.visframeright.drawArea.draw()

	def changeCompoInfo(self, df_index):
		#print(self.actCompoMode, df_index)
		if self.actCompoMode == 'Groom element':
			self.bottomframe.frames[CompositionPanelGroEle].setActive()

			self.bottomframe.frames[CompositionPanelGroEle].componame.config(text=self.groom_df.at[df_index, 'Name'])
			self.bottomframe.frames[CompositionPanelGroEle].formulaety.delete(0,tk.END)
			self.bottomframe.frames[CompositionPanelGroEle].formulaety.insert(0, self.groom_df.at[df_index, 'Formula'])
			self.bottomframe.frames[CompositionPanelGroEle].zety.delete(0,tk.END)
			self.bottomframe.frames[CompositionPanelGroEle].zety.insert(0, self.groom_df.at[df_index, 'Z'])
			self.bottomframe.frames[CompositionPanelGroEle].aety.delete(0,tk.END)
			self.bottomframe.frames[CompositionPanelGroEle].aety.insert(0, self.groom_df.at[df_index, 'A'])
			self.bottomframe.frames[CompositionPanelGroEle].rhoety.delete(0,tk.END)
			self.bottomframe.frames[CompositionPanelGroEle].rhoety.insert(0, self.groom_df.at[df_index, 'Density'])
			self.bottomframe.frames[CompositionPanelGroEle].Iety.delete(0,tk.END)
			self.bottomframe.frames[CompositionPanelGroEle].Iety.insert(0, self.groom_df.at[df_index, 'I'])
			self.bottomframe.frames[CompositionPanelGroEle].Cbarety.delete(0,tk.END)
			self.bottomframe.frames[CompositionPanelGroEle].Cbarety.insert(0, self.groom_df.at[df_index, 'Cbar'])
			self.bottomframe.frames[CompositionPanelGroEle].x0ety.delete(0,tk.END)
			self.bottomframe.frames[CompositionPanelGroEle].x0ety.insert(0, self.groom_df.at[df_index, 'x0'])
			self.bottomframe.frames[CompositionPanelGroEle].x1ety.delete(0, tk.END)
			self.bottomframe.frames[CompositionPanelGroEle].x1ety.insert(0, self.groom_df.at[df_index, 'x1'])
			self.bottomframe.frames[CompositionPanelGroEle].aaety.delete(0,tk.END)
			self.bottomframe.frames[CompositionPanelGroEle].aaety.insert(0, self.groom_df.at[df_index, 'a'])
			self.bottomframe.frames[CompositionPanelGroEle].ksety.delete(0,tk.END)
			self.bottomframe.frames[CompositionPanelGroEle].ksety.insert(0, self.groom_df.at[df_index, 'k'])
			self.bottomframe.frames[CompositionPanelGroEle].d0ety.delete(0,tk.END)
			self.bottomframe.frames[CompositionPanelGroEle].d0ety.insert(0, self.groom_df.at[df_index, 'dlt0'])

			self.bottomframe.frames[CompositionPanelGroEle].setDisabled()
		elif self.actCompoMode == 'Groom mixture':
			self.bottomframe.frames[CompositionPanelGroMix].setActive()

			self.bottomframe.frames[CompositionPanelGroMix].componame.config(text=self.groom_df.at[df_index, 'Name'])
			self.bottomframe.frames[CompositionPanelGroMix].formulaety.delete(0, tk.END)
			self.bottomframe.frames[CompositionPanelGroMix].formulaety.insert(0, '')
			self.bottomframe.frames[CompositionPanelGroMix].zaety.delete(0, tk.END)
			self.bottomframe.frames[CompositionPanelGroMix].zaety.insert(0, self.groom_df.at[df_index, 'Z/A'])
			self.bottomframe.frames[CompositionPanelGroMix].aety.delete(0, tk.END)
			self.bottomframe.frames[CompositionPanelGroMix].aety.insert(0, self.groom_df.at[df_index, 'A'])
			self.bottomframe.frames[CompositionPanelGroMix].rhoety.delete(0, tk.END)
			self.bottomframe.frames[CompositionPanelGroMix].rhoety.insert(0, self.groom_df.at[df_index, 'Density'])
			self.bottomframe.frames[CompositionPanelGroMix].Iety.delete(0, tk.END)
			self.bottomframe.frames[CompositionPanelGroMix].Iety.insert(0, self.groom_df.at[df_index, 'I'])
			self.bottomframe.frames[CompositionPanelGroMix].Cbarety.delete(0, tk.END)
			self.bottomframe.frames[CompositionPanelGroMix].Cbarety.insert(0, self.groom_df.at[df_index, 'Cbar'])
			self.bottomframe.frames[CompositionPanelGroMix].x0ety.delete(0, tk.END)
			self.bottomframe.frames[CompositionPanelGroMix].x0ety.insert(0, self.groom_df.at[df_index, 'x0'])
			self.bottomframe.frames[CompositionPanelGroMix].x1ety.delete(0, tk.END)
			self.bottomframe.frames[CompositionPanelGroMix].x1ety.insert(0, self.groom_df.at[df_index, 'x1'])
			self.bottomframe.frames[CompositionPanelGroMix].aaety.delete(0, tk.END)
			self.bottomframe.frames[CompositionPanelGroMix].aaety.insert(0, self.groom_df.at[df_index, 'a'])
			self.bottomframe.frames[CompositionPanelGroMix].ksety.delete(0, tk.END)
			self.bottomframe.frames[CompositionPanelGroMix].ksety.insert(0, self.groom_df.at[df_index, 'k'])
			self.bottomframe.frames[CompositionPanelGroMix].d0ety.delete(0, tk.END)
			self.bottomframe.frames[CompositionPanelGroMix].d0ety.insert(0, self.groom_df.at[df_index, 'dlt0'])

			self.bottomframe.frames[CompositionPanelGroMix].setDisabled()

	def loadOxideData(self):
		filedir = os.path.split(os.path.abspath(__file__))[0]
		data_dir = os.path.normpath(os.path.join(filedir, 'Input_Material', 'OxideData'))

		fileloc = tk.filedialog.askopenfilename(initialdir=data_dir, filetypes=[('CSV files', '*.csv')])

		if fileloc:
			# Check for possible data file errors
			try:
				df = pd.read_csv(fileloc, sep=None, engine='python', index_col=0).transpose()
			except pde.EmptyDataError:
				messagebox.showerror('Error', 'Chosen csv file is empty')
				return

			df = df.div(df.sum(axis=1), axis=0)

			oxides = df.columns.to_list()
			ele_keys = set()
			for oxi in oxides:
				ele, sto = self.parse_formula(oxi)
				ele_keys.update(ele)
			ele_keys = list(ele_keys)

			# Put O at end
			ele_keys.remove('O')
			ele_keys.insert(len(ele_keys), 'O')

			# Compose transformation matrix (oxides -> elements)
			Trans_mat = np.zeros((len(df.columns.to_list()), len(ele_keys)))

			for oxi in oxides:
				element, stochio = self.parse_formula(oxi)
				row_index = oxides.index(oxi)
				oxide_mass = 0
				for ide, ele in enumerate(element):
					elerow = self.groom_df.loc[(self.groom_df['Formula'] == ele) & (self.groom_df['Type'] == 'E')]
					ele_mass = float(stochio[ide]) * elerow['A'].values[0]
					oxide_mass += ele_mass
					Trans_mat[row_index, ele_keys.index(ele)] = ele_mass

				Trans_mat[row_index, :] /= oxide_mass

			element_values = np.matmul(df.values, Trans_mat)
			element_values /= element_values[:, np.size(element_values[1]) - 1, None]
			element_values = np.log(element_values)

			# Log-ratio dataframe
			eledf = pd.DataFrame(data=element_values, index=df.index, columns=ele_keys)

			mean_lr = eledf.mean(axis=0)
			cov_lr_df = eledf.cov()

			print(mean_lr)
			print(cov_lr_df)

			# Plot
			n_subplots = len(ele_keys)

			n_plot_col = int(4)
			n_plot_row = int(np.ceil(n_subplots/n_plot_col))

			self.visframeright.fig.clear()

			self.visframeright.fig.suptitle(r'Plot of $z_r$ vs order statistics')
			gs = self.visframeright.fig.add_gridspec(n_plot_row+1, n_plot_col)

			i = 0
			for row in range(n_plot_row):
				for col in range(n_plot_col):

					if i < n_subplots:
						axi = self.visframeright.fig.add_subplot(gs[row, col])

						dat_size = element_values.shape

						test_statistic = []
						order_statistic = []

						if i < n_subplots-1:
							for j in range(0,dat_size[0]):
								test_statistic.append(sps.norm.cdf((element_values[j,i]-mean_lr[i])/np.sqrt(cov_lr_df.values[i,i])))
								order_statistic.append(float(2*j+1)/float(2*dat_size[0]))

						test_statistic.sort()
						axi.scatter(order_statistic,test_statistic)
						axi.add_line(mlines.Line2D([0,1],[0,1],color='r'))

						axi.set_xlim([0,1])
						axi.set_ylim([0,1])
						axi.set_xlabel('Order statistics')
						axi.set_ylabel(r'$z_r$')
						axi.set_aspect('equal')

						axi.set_title(ele_keys[i])
						axi.grid()
						i += 1


			axmt = self.visframeright.fig.add_subplot(gs[n_plot_row, :])
			axmt.set_axis_off()
			axmt.set_title('Mean log-ratio with last component')

			tabdata = mean_lr.to_numpy()
			td_formt = ['%.2f' % j for j in list(tabdata)]

			columns = mean_lr.index

			ytable = axmt.table([list(td_formt)], colLabels=columns)

			cellDict = ytable.get_celld()
			for i in range(0, len(columns)):
				cellDict[(0, i)].set_height(.3)
				cellDict[(1, i)].set_height(.3)

				for j in range(1, 1):
					cellDict[(j, i)].set_height(.2)

			#self.visframeright.fig.subplots_adjust(left=0.1, bottom=0.16, right=0.90, top=0.92, wspace=0.67, hspace=0.67)
			self.visframeright.drawArea.draw()

	def parse_formula(self, formula):
		ele_list = []
		stochi_idx = []

		end_of_formula = False

		idx = 0

		act_word = ''
		act_numb = ''
		while not end_of_formula:
			act_char = formula[idx]
			if idx + 1 == len(formula):
				end_of_formula = True
			nidx = (idx + 1) % len(formula)
			next_char = formula[nidx]

			act_word += act_char
			idx = idx + 1

			if next_char.isupper():
				ele_list.append(act_word)
				act_word = ''
				stochi_idx.append('1')

			elif next_char.isdigit():
				ele_list.append(act_word)
				act_word = ''
				while True:
					act_char = formula[idx]
					if idx + 1 == len(formula):
						end_of_formula = True
					nidx = (idx + 1) % len(formula)
					next_char = formula[nidx]

					act_numb += act_char
					idx = idx + 1
					if next_char.isupper():
						stochi_idx.append(act_numb)
						act_numb = ''
						break

		return [ele_list, stochi_idx]

	def densityModeChange(self, mode):
		self.actDenMode = mode

		if mode == 'Normal':
			pass
		elif mode == 'Log-normal':
			pass
		elif mode == 'From data':
			pass

	def showdensityplot(self, mu, std):
		if self.actDenMode == 'Normal':
			self.visframeleft.fig.clear()
			gs = self.visframeleft.fig.add_gridspec(2, 1)
			ax = self.visframeleft.fig.add_subplot(gs[0])

			if std == 0:
				ax.axvline(mu)
			else:
				x = np.linspace(mu-5*std, mu+5*std, 200)
				y = sps.norm.pdf(x, mu, std)
				ax.plot(x, y)
				ax.set_xlim([mu - 5 * std, mu + 5 * std])

			ax.set_xlabel('Density')
			ax.set_ylabel('Probability density')
			ax.set_title('Normal distribution')

			ax.grid()

			self.visframeleft.drawArea.draw()

		elif self.actDenMode == 'Log-normal':
			self.visframeleft.fig.clear()
			gs = self.visframeleft.fig.add_gridspec(2, 1)
			ax = self.visframeleft.fig.add_subplot(gs[0])

			if std == 0:
				ax.axvline(mu)
			else:
				varlog = np.log(1 + (std/mu)**2)
				siglog = varlog**0.5
				meanlog = np.log(mu)-0.5*varlog

				x = np.linspace(np.exp(meanlog - 5 * siglog), np.exp(meanlog + 5 * siglog), 200)
				y = sps.lognorm.pdf(x, scale=np.exp(meanlog), s=siglog)
				ax.plot(x, y)
				ax.set_xlim([np.exp(meanlog - 5 * siglog), np.exp(meanlog + 5 * siglog)])

			ax.set_xlabel('Density')
			ax.set_ylabel('Probability density')
			ax.set_title('Lognormal distribution')

			ax.grid()

			self.visframeleft.drawArea.draw()

		elif self.actDenMode == 'From data':
			pass

	def loadDensityData(self):
		def lognorm_func(x, mu, sigma):
			return 1. / ((2. * np.pi * sigma ** 2) ** 0.5 * x) * np.exp(-0.5 * ((np.log(x) - mu) / (sigma)) ** 2)


		filedir = os.path.split(os.path.abspath(__file__))[0]
		data_dir = os.path.normpath(os.path.join(filedir, 'Input_Material', 'DensityData'))

		projectloc = tk.filedialog.askdirectory(initialdir=data_dir, title='Choose project folder')

		if projectloc:
			if self.bottomframe.denframes[DensityPanelData].selectvar.get() == 0:
				mode_str = self.bottomframe.denframes[DensityPanelData].modevar.get()
				modefile = mode_str.lower() + '_density.csv'
				save_plot = os.path.normpath(os.path.join(projectloc, 'plots'))
				save_data = os.path.normpath(os.path.join(projectloc, 'samples', mode_str))
				fileloc = os.path.normpath(os.path.join(projectloc, 'raw_data', modefile))

				data = pd.read_csv(fileloc, sep=None, engine='python')

				# Extract sample names
				sample_names = data['sample'].unique()


				if mode_str == 'Skeletal' or mode_str == 'Grain':
					mode_str = mode_str.lower()

					# Loop through all subsamples
					#TODO: Progressbar
					#progbar = ttk.Progressbar()

					for name in sample_names:
						# Extract partial dataframe
						partial_df = data.loc[data['sample'] == name]
						partial_df.reset_index(drop=True, inplace=True)

						# Create temporary dataframe
						temp_data = []

						# Fill with data and do calculations
						for idx, row in partial_df.iterrows():
							subsample_name = name + '-' + str(idx + 1)

							# m
							var_logm = np.log(1 + row['mass std'] ** 2 / row['mass'] ** 2)
							mu_logm = np.log(row['mass']) - 0.5 * var_logm

							acc_v = 0.0003 * row['sample volume'] + 0.0003 * row[
								'chamber volume']  # accuracy of measurement
							err_v = np.sqrt(acc_v ** 2 + row['volume std'] ** 2)
							var_logv = np.log(1 + err_v ** 2 / row['sample volume'] ** 2)
							mu_logv = np.log(row['sample volume']) - 0.5 * var_logv

							mu_logrho = mu_logm - mu_logv
							var_logrho = var_logm + var_logv
							sig_logrho = np.sqrt(var_logrho)

							temp_data.append((subsample_name, sig_logrho, 0, np.exp(mu_logrho)))

							fig1 = plt.figure()
							ax1 = fig1.add_subplot(1, 1, 1)

							x_dat = np.linspace(np.exp(mu_logrho - 5 * sig_logrho),
												np.exp(mu_logrho + 5 * sig_logrho), 1000)
							y_dat = sps.lognorm.pdf(x_dat, scale=np.exp(mu_logrho), s=sig_logrho)

							ax1.plot(x_dat, y_dat, color='orange')
							ax1.set_xlim([np.exp(mu_logrho - 5 * sig_logrho), np.exp(mu_logrho + 5 * sig_logrho)])
							ax1.set_axisbelow(True)
							ax1.grid(linestyle='dashed')
							ax1.set_title('Lognorm for %s density of %s' % (mode_str, subsample_name))

							fitline = mlines.Line2D([0, 1], [0, 0], color='orange')

							leg_handles = [fitline]
							leg_lables = ['Lognormal pdf']

							ax1.legend(leg_handles, leg_lables)
							ax1.set_xlabel(r'Material density $[g * cm^{-3}]$')
							ax1.set_ylabel(r'Probability density $[g * cm^{-3}]^{-1}$')

							fig1.savefig(os.path.normpath(os.path.join(save_plot, mode_str, subsample_name + '.png')))
							plt.close(fig1)

						# Export subsamples to csv files
						temp_df = pd.DataFrame(temp_data, columns=['sample_nr', 'shape_log_rho', 'loc_log_rho',
																   'scale_log_rho'])
						temp_df.to_csv(os.path.normpath(os.path.join(save_data, name + '.csv')), index=False)
				elif mode_str == 'Bulk':
					mode_str = mode_str.lower()
					n_mc_steps = 10000

					progbar = ttk.Progressbar(self.visframeleft, orient='horizontal', length=50)
					progbar.grid(row=0, column=0, sticky='e')

					# Loop through all subsamples
					for name in sample_names:
						print('here')
						# Extract partial dataframe
						partial_df = data.loc[data['sample'] == name]
						partial_df.reset_index(drop=True, inplace=True)

						# Create temporary dataframe
						temp_data = []

						# Fill with data and do calculations
						for idx, row in partial_df.iterrows():
							subsample_name = name + '-' + str(idx + 1)

							# Monte Carlo simulation
							rho_list = []
							for i in range(0, n_mc_steps):
								# m_sus
								var_log_msus = np.log(1 + row['e_m_sus'] ** 2 / row['m_sus'] ** 2)
								mu_log_msus = np.log(row['m_sus']) - 0.5 * var_log_msus
								msus = sps.lognorm.rvs(scale=np.exp(mu_log_msus), s=np.sqrt(var_log_msus))

								# rho_h2o
								var_log_rhoh2o = np.log(1 + row['e_rho_h2o'] ** 2 / row['rho_h2o'] ** 2)
								mu_log_rhoh2o = np.log(row['rho_h2o']) - 0.5 * var_log_rhoh2o
								rhoh2o = sps.lognorm.rvs(scale=np.exp(mu_log_rhoh2o), s=np.sqrt(var_log_rhoh2o))

								# rho_t
								var_log_rhot = np.log(1 + row['e_rho_t'] ** 2 / row['rho_t'] ** 2)
								mu_log_rhot = np.log(row['rho_t']) - 0.5 * var_log_rhot
								rhot = sps.lognorm.rvs(scale=np.exp(mu_log_rhot), s=np.sqrt(var_log_rhot))

								# rho_p
								var_log_rhop = np.log(1 + row['e_rho_p'] ** 2 / row['rho_p'] ** 2)
								mu_log_rhop = np.log(row['rho_p']) - 0.5 * var_log_rhop
								rhop = sps.lognorm.rvs(scale=np.exp(mu_log_rhop), s=np.sqrt(var_log_rhop))

								mt = -1
								mp = -1

								while mt <= 0 and mp <= 0:
									# m_stp
									var_log_mstp = np.log(1 + row['e_m_stp'] ** 2 / row['m_stp'] ** 2)
									mu_log_mstp = np.log(row['m_stp']) - 0.5 * var_log_mstp
									mstp = sps.lognorm.rvs(scale=np.exp(mu_log_mstp), s=np.sqrt(var_log_mstp))

									# m_st
									var_log_mst = np.log(1 + row['e_m_st'] ** 2 / row['m_st'] ** 2)
									mu_log_mst = np.log(row['m_st']) - 0.5 * var_log_mst
									mst = sps.lognorm.rvs(scale=np.exp(mu_log_mst), s=np.sqrt(var_log_mst))

									# m_s
									var_log_ms = np.log(1 + row['e_m_s'] ** 2 / row['m_s'] ** 2)
									mu_log_ms = np.log(row['m_s']) - 0.5 * var_log_ms
									ms = sps.lognorm.rvs(scale=np.exp(mu_log_ms), s=np.sqrt(var_log_ms))

									# m_rt
									var_log_mrt = np.log(1 + row['e_m_rt'] ** 2 / row['m_rt'] ** 2)
									mu_log_mrt = np.log(row['m_rt']) - 0.5 * var_log_mrt
									mrt = sps.lognorm.rvs(scale=np.exp(mu_log_mrt), s=np.sqrt(var_log_mrt))

									mp = mstp - mst
									mt = mst - ms - mrt

								rho = (rhoh2o * ms) / (mstp - msus - mt * rhoh2o / rhot - mp * rhoh2o / rhop)
								rho_list.append(rho)


							sha, lo, sca = sps.lognorm.fit(rho_list, floc=0)

							histo_edges = np.histogram_bin_edges(rho_list, bins='auto')
							lnorm_func = sps.lognorm.pdf(histo_edges, s=sha, loc=lo, scale=sca)

							fig1 = plt.figure()
							ax1 = fig1.add_subplot(1, 1, 1)
							ax1.set_title('Histogram & Lognorm-fit for %s' % subsample_name)
							ax1.set_axisbelow(True)
							ax1.grid(linestyle='dashed')
							ax1.hist(rho_list, bins=histo_edges, density=True, edgecolor='black', color='green')
							ax1.plot(histo_edges, lnorm_func, color='orange')
							# ax1.plot(histo_edges,norm_func)

							histopatch = mpatches.Rectangle((0, 0), 1, 1, fc='green', ec='black')
							fitline = mlines.Line2D([0, 1], [0, 0], color='orange')

							leg_handles = [histopatch, fitline]
							leg_lables = ['MC simulation', 'Lognormal pdf fit']

							ax1.legend(leg_handles, leg_lables)
							ax1.set_xlabel(r'Material density $[g * cm^{-3}]$')
							ax1.set_ylabel(r'Probability density $[g * cm^{-3}]^{-1}$')

							fig1.savefig(os.path.normpath(os.path.join(save_plot, mode_str, subsample_name + '.png')))
							plt.close(fig1)

							temp_data.append((subsample_name, sha, lo, sca))

						# Export subsamples to csv files
						temp_df = pd.DataFrame(temp_data, columns=['sample_nr', 'shape_log_rho', 'loc_log_rho',
																   'scale_log_rho'])
						temp_df.to_csv(os.path.normpath(os.path.join(save_data, name + '.csv')), index=False)
						progbar.step(1/len(sample_names)*100)
						self.update_idletasks()

					progbar.destroy()


			mode_str = self.bottomframe.denframes[DensityPanelData].modevar.get()
			mode = mode_str.lower()
			fileloc = os.path.normpath(os.path.join(projectloc, 'samples', mode))

			print(fileloc)

			if os.path.exists(fileloc):
				self.visframeleft.fig.clear()

				xlim_sp = (1, 4)
				ylim_sp = (0, 70)

				# Set working directory to path of python script
				absolute_path = os.path.abspath(__file__)
				directory_name = os.path.dirname(absolute_path)
				os.chdir(directory_name)


				############### Samples #################
				# List all sample files in data folder...
				sample_files = os.listdir(fileloc)

				# ... and calculate the number of plots & subplots
				n_plots = len(sample_files)
				n_subplots = math.ceil(np.sqrt(n_plots))

				# Prepare dataframes
				samples = pd.DataFrame(columns=['sample_name', 'shape_log_rho', 'loc_log_rho', 'scale_log_rho'])

				# Prepare x & y axes
				x_dat = np.linspace(xlim_sp[0], xlim_sp[1], 1000)
				y_litho = np.zeros_like(x_dat)

				# Prepare subplots
				supergs = self.visframeleft.fig.add_gridspec(2,1)
				gs = supergs[0].subgridspec(n_subplots+1, n_subplots)#, wspace=0.75, hspace=1.5)
					#matplotlib.gridspec.GridSpec(n_subplots+1, n_subplots)

				merg = []
				fitt = []

				# Initialise progress bar
				n_progress = len(sample_files) + 3

#					with click.progressbar(length=n_progress) as bar:
				progbar = ttk.Progressbar(self.visframeleft, orient='horizontal', length=100)
				progbar.grid(row=0, column=0, sticky='e')

				# Loop through all samples
				y_samp_max = 0

				for isamp, sample in enumerate(sample_files):

					# Retrieve sample name
					sample_name = sample[:-4]
					path_to_file = os.path.normpath(os.path.join(fileloc, sample))

					# Read sample data into dataframe
					data = pd.read_csv(path_to_file)

					# Prepare y-axis for different solutions
					y_dat = np.zeros_like(x_dat)
					y_fit = np.zeros_like(x_dat)

					# Get subplot indices
					jj = isamp % n_subplots
					ii = int(isamp / n_subplots)

					ax = self.visframeleft.fig.add_subplot(gs[ii, jj])

					# Get mean and standard deviation for each subsample
					for iss in range(0, data.shape[0]):
						sca = data.iloc[iss]['scale_log_rho']
						lo = data.iloc[iss]['loc_log_rho']
						sha = data.iloc[iss]['shape_log_rho']

						# Perform a disjunction of pdfs (see Tarantola, 2005)
						for i in range(0, x_dat.shape[0]):
							y_dat[i] += sps.lognorm.pdf(x_dat[i], s=sha, loc=lo, scale=sca)

					# Normalise disjunction of pdf's
					y_dat /= data.shape[0]

					if np.amax(y_dat) > y_samp_max:
						y_samp_max = np.amax(y_dat)

					# Fit a gaussian to each group of subsamples (fit a sample gaussian)
					# popt, pcov = spo.curve_fit(gauss_func, x_dat, y_dat, p0=[data['mean_rho'].mean(), data['std_rho'].mean()])
					popt, pcov = spo.curve_fit(lognorm_func, x_dat, y_dat)

					sca_fit = np.exp(popt[0])
					lo_fit = 0
					sha_fit = popt[1]

					# sca_fit, lo_fit, sha_fit = sps.lognorm.fit(y_dat, floc=0)

					for i in range(0, x_dat.shape[0]):
						# Save gaussian values in y_fit
						y_fit[i] = sps.lognorm.pdf(x_dat[i], s=sha_fit, loc=lo_fit, scale=sca_fit)
						# Perform disjunction of sample pdfs to lithology pdf
						y_litho[i] += y_dat[i] / n_plots

					# Store gaussian parameters for each sample in dataframe
					samples.loc[isamp] = [sample_name, sca_fit, lo_fit, sha_fit]

					# Add plots to list (for legend later) & plot subplots
					merg.append(ax.plot(x_dat, y_dat, color='#1f77b4', linestyle='-'))
					fitt.append(ax.plot(x_dat, y_fit, color='#2ca02c', linestyle='--'))

					# Format subplots
					ax.set_title(sample_name)
					if ii == n_subplots - 1:
						ax.set_xlabel(r'$\rho$')

					if jj == 0:
						ax.set_ylabel('prob. density')

					ax.set_xlim([np.exp(np.log(sca_fit) - 5 * sha_fit), np.exp(np.log(sca_fit) + 5 * sha_fit)])
					ax.set_ylim(ylim_sp[0], ylim_sp[1])
					ax.xaxis.set_major_locator(plt.MaxNLocator(2))
					ax.grid()

					progbar.step(1/float(n_progress)*100)
					self.update_idletasks()

				# # Put subplot legend in first unused subplot
				# # get indices
				# ilegend = n_plots % n_subplots + 1
				# jlegend = int(n_plots / n_subplots)
				#
				# # Width & Height of legend
				# dl = 1. / float(n_subplots)
				#
				# # Get location in figure frame
				# lloc = ilegend * dl
				#
				# # Plot legend
				# pdf_line = mlines.Line2D([], [], color='#1f77b4', linestyle='-', label='merged pdfs')
				# gauss_line = mlines.Line2D([], [], color='#2ca02c', linestyle='--', label='fitted lognormal')

				#self.visframeleft.fig.legend(handles=[pdf_line, gauss_line])
											 # , bbox_to_anchor=(lloc, 0., dl, dl),
											 # framealpha=1.0)
				#self.visframeleft.fig.subplots_adjust(left=0.1, bottom=0.16, right=0.90, top=0.8, wspace=0.75,
				#									  hspace=0.75)

				#progbar.destroy()

				############ Lithology ##################

				# Normalise pdf
				lnorm_pdf = 0.
				for i in range(0, x_dat.shape[0]):
					dx = x_dat[i] - x_dat[i - 1]
					pdf_trp = (y_litho[i] + y_litho[i - 1]) / 2
					lnorm_pdf += pdf_trp * dx

				progbar.step(1/float(n_progress)*100)

				y_litho /= lnorm_pdf

				# Fit gaussian to lithology distribution
				# popt_lith, pcov_lith = spo.curve_fit(gauss_func, x_dat, y_litho, p0=[samples['mean'].mean(), samples['std'].mean()])
				popt, pcov = spo.curve_fit(lognorm_func, x_dat, y_litho)

				sca_fit = np.exp(popt[0])
				lo_fit = 0
				sha_fit = popt[1]

				self.bottomframe.denframes[DensityPanelData].meanvar = sca_fit
				self.bottomframe.denframes[DensityPanelData].stdvar = sha_fit

				# sca_fit, lo_fit, sha_fit = sps.lognorm.fit(y_litho, floc=0)

				# Save gaussian values in y_litho_fit
				y_litho_fit = np.zeros_like(x_dat)
				for i in range(0, x_dat.shape[0]):
					y_litho_fit[i] = sps.lognorm.pdf(x_dat[i], s=sha_fit, loc=lo_fit, scale=sca_fit)

				progbar.step(1/float(n_progress)*100)

				# Calculate empirical cdf
				y_litho_cdf = np.zeros_like(x_dat)
				y_fit_cdf = np.zeros_like(x_dat)

				# perform trapezoidal integration
				# and estimate quantiles and mean
				q25 = 0.
				q50 = 0.
				q75 = 0.
				ypdf25 = 0.
				ypdf50 = 0.
				ypdf75 = 0.
				mean = 0.
				ymean = 0.

				for i in range(1, x_dat.shape[0]):

					# Lognormal cdf
					y_fit_cdf[i] = sps.lognorm.cdf(x_dat[i], s=sha_fit, loc=lo_fit, scale=sca_fit)

					# Trapezoidal integration
					dx = x_dat[i] - x_dat[i - 1]
					pdf_trp = (y_litho[i] + y_litho[i - 1]) / 2
					y_litho_cdf[i] = y_litho_cdf[i - 1] + pdf_trp * dx
					# update mean
					mean_trp = (y_litho[i] * x_dat[i] + y_litho[i - 1] * x_dat[i - 1]) / 2
					mean += mean_trp * dx

					# check for quantile crossings
					# 25% - quantile
					if y_litho_cdf[i] > 0.25 and q25 == 0:
						dy = y_litho_cdf[i] - y_litho_cdf[i - 1]
						q25 = dx / dy * (0.25 - y_litho_cdf[i - 1]) + x_dat[i - 1]
						ypdf25 = y_litho[i - 1] + (y_litho[i] - y_litho[i - 1]) / dx * (q25 - x_dat[i - 1])
					# 50% - quantile (median)
					if y_litho_cdf[i] > 0.5 and q50 == 0:
						dy = y_litho_cdf[i] - y_litho_cdf[i - 1]
						q50 = dx / dy * (0.5 - y_litho_cdf[i - 1]) + x_dat[i - 1]
						ypdf50 = y_litho[i - 1] + (y_litho[i] - y_litho[i - 1]) / dx * (q50 - x_dat[i - 1])
					# 75% - quantile
					if y_litho_cdf[i] > 0.75 and q75 == 0:
						dy = y_litho_cdf[i] - y_litho_cdf[i - 1]
						q75 = dx / dy * (0.75 - y_litho_cdf[i - 1]) + x_dat[i - 1]
						ypdf75 = y_litho[i - 1] + (y_litho[i] - y_litho[i - 1]) / dx * (q75 - x_dat[i - 1])

				progbar.step(1/float(n_progress)*100)

				# # for i in range(1, x_dat.shape[0]):
				# # 	if x_dat[i] > mean:
				# # 		dx = x_dat[i] - x_dat[i - 1]
				# # 		ymean = y_litho[i - 1] + (y_litho[i] - y_litho[i - 1]) / dx * (mean - x_dat[i - 1])
				# # 		break
				#
				# Adapt axes in sample plot
				axes = self.visframeleft.fig.get_axes()
				for ax in axes:
					ax.set_ylim([0, y_samp_max * 1.1])

				gs2 = supergs[1].subgridspec(1, 1)
				ax2 = self.visframeleft.fig.add_subplot(gs2[0])

				# Lognormal lithology quartiles
				gq75 = sps.lognorm.ppf(0.75, s=sha_fit, loc=lo_fit, scale=sca_fit)
				gq50 = sps.lognorm.ppf(0.50, s=sha_fit, loc=lo_fit, scale=sca_fit)
				gq25 = sps.lognorm.ppf(0.25, s=sha_fit, loc=lo_fit, scale=sca_fit)

				# Plot pdfs
				ax2.plot(x_dat, y_litho, label='merged pdfs')
				ax2.plot(x_dat, y_litho_fit, linestyle='--', color='green', label='fitted lognormal')

				# Plot quantiles
				ax2.vlines(x=q50, ymin=0, ymax=ypdf50, color='blue', linestyles='dotted', alpha=0.5,
						   label=r'$q_{25,pdf}$, $q_{50,pdf}$, $q_{75,pdf}$')
				ax2.vlines(x=q75, ymin=0, ymax=ypdf75, color='blue', linestyles='dotted', alpha=0.5)
				ax2.vlines(x=q25, ymin=0, ymax=ypdf25, color='blue', linestyles='dotted', alpha=0.5)

				# ax2.vlines(x=mean, ymin=0, ymax=ymean, color='orange', linestyle='dotted', alpha=0.5, label=r'$\mu_{pdf}$')

				ax2.vlines(x=gq75, ymin=0, ymax=sps.lognorm.pdf(gq75, s=sha_fit, loc=lo_fit, scale=sca_fit),
						   color='green', linestyles='dashdot',
						   alpha=0.5, label=r'$q_{25,LN}$, $q_{50,LN}$, $q_{75,LN}$')
				ax2.vlines(x=gq50, ymin=0, ymax=sps.lognorm.pdf(gq50, s=sha_fit, loc=lo_fit, scale=sca_fit),
						   color='green', linestyles='dashdot', alpha=0.5)
				ax2.vlines(x=gq25, ymin=0, ymax=sps.lognorm.pdf(gq25, s=sha_fit, loc=lo_fit, scale=sca_fit),
						   color='green', linestyles='dashdot',
						   alpha=0.5)

				# Plot format
				ax2.legend(loc='upper right', framealpha=1.0)
				ax2.grid()
				ax2.set_xlabel(r'Material density $[g * cm^{-3}]$')
				ax2.set_ylabel(r'Probability density $[g * cm^{-3}]^{-1}$')
				ax2.set_title('Lithology ¦ PDF')
				ax2.set_xlim([np.exp(np.log(sca_fit) - 5 * sha_fit), np.exp(np.log(sca_fit) + 5 * sha_fit)])

				# # Plot cumulative distribution function
				# fig3 = plt.figure(3)
				# ax3 = fig3.add_subplot(1, 1, 1)
				# ax3.plot(x_dat, y_litho_cdf, label='merged pdfs')
				# ax3.plot(x_dat, y_fit_cdf, linestyle='--', color='green', label='fitted lognormal')
				#
				# ax3.vlines(x=q25, ymin=0, ymax=0.25, color='blue', linestyles='dotted', alpha=0.5,
				# 		   label=r'$q_{25,pdf}$, $q_{50,pdf}$, $q_{75,pdf}$')
				# ax3.vlines(x=q50, ymin=0, ymax=0.5, color='blue', linestyles='dotted', alpha=0.5)
				# ax3.vlines(x=q75, ymin=0, ymax=0.75, color='blue', linestyles='dotted', alpha=0.5)
				# ax3.vlines(x=gq25, ymin=0, ymax=0.25, color='green', linestyles='dashdot',
				# 		   alpha=0.5, label=r'$q_{25,lognormal}$, $q_{50,lognormal}$, $q_{75,lognormal}$')
				# ax3.vlines(x=gq50, ymin=0, ymax=0.5, color='green', linestyle='dashdot', alpha=0.5)
				# ax3.vlines(x=gq75, ymin=0, ymax=0.75, color='green', linestyles='dashdot',
				# 		   alpha=0.5)
				#
				# ax3.legend(loc='upper left', framealpha=1.0)
				# ax3.set_title('Lithology ¦ CDF')
				#
				# ax3.set_xlabel(r'Material density $[g * cm^{-3}]$')
				# ax3.set_ylabel(r'Probability')
				#
				# if shrink_support:
				# 	ax3.set_xlim([np.exp(np.log(sca_fit) - 5 * sha_fit), np.exp(np.log(sca_fit) + 5 * sha_fit)])
				#
				# ax3.grid()

				# # Output to json files
				# litho_pdf_dict = {'pdf_class': 'interpolate', 'parameter': 'density', 'x': x_dat.tolist(),
				# 				  'y': y_litho.tolist()}
				# litho_gauss_dict = {'pdf_class': 'normal', 'parameter': 'density', 'scale': sca_fit,
				# 					'location': lo_fit, 'shape': sha_fit}
				#
				# with open(save_dir + filename + '_pdf.json', 'w') as fp:
				# 	json.dump(litho_pdf_dict, fp)
				#
				# with open(save_dir + filename + '_gauss.json', 'w') as fp:
				# 	json.dump(litho_gauss_dict, fp)

				self.visframeleft.drawArea.draw()
			else:
				messagebox.showerror('Error', mode + ' folder not found')

	def getElementInfo(self, name):
		return self.groom_df.loc[self.groom_df['Name'] == name]

	def getGroomDF(self):
		return self.groom_df



def main():
	root = tk.Tk()
	root.title('Material creation app')

	screen_width = root.winfo_screenwidth()
	screen_height = root.winfo_screenheight()
	sizestring = str(int(0.8 * screen_width)) + 'x' + str(int(0.8 * screen_height)) + '+50+50'
	root.geometry(sizestring)

	window = MainWindow(master=root)
	window.pack(side="top", fill="both", expand=True)
	root.mainloop()


if __name__ == '__main__':
	main()

import tkinter as tk
import tkinter.ttk as ttk

import os
from distutils.dir_util import copy_tree
from shutil import copy2
from collections import OrderedDict
from tkinter import filedialog
from tkinter import ttk
from tkinter import messagebox
from tkinter import font
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection, LineCollection
from matplotlib import collections
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from io import StringIO
import json
from myListboxes import *


class JSONEncoder(json.JSONEncoder):
	def default(self, obj):
		if hasattr(obj, 'to_json'):
			return obj.to_json(orient='records')
		return json.JSONEncoder.default(self, obj)


class LeftWindow(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)
		self.__createlayout()

	def __createlayout(self):
		self.unbinlabel = tk.Label(self, text='Unused bins')
		self.unbinlistbox = BinListbox(self, (('ID', 10), ('# Tracks', 10), (chr(952), 10), (chr(966), 10)))

		self.unbinlabel.grid(column=0, row=0, padx=2, pady=2, sticky='w')
		self.unbinlistbox.grid(column=0, row=1, sticky='nsew')

		self.rowconfigure(1, weight=1)

	def highlightSelection(self, row_idx):
		self.parent.highlightSelection(row_idx)


class RightWindow(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)
		self.__createlayout()

	def __createlayout(self):
		self.grouplabel = tk.Label(self, text='Groups')
		self.grouplistbox = GroupListbox(self, (('Name', 30), ('# Bins', 10)))
		self.materialabel = tk.Label(self, text='Materials')
		self.materialistbox = SingleListbox(self, (('Material (Order: first = near Detector)', 40)))
		self.binlabel = tk.Label(self, text='Attached bins')
		self.usedbinlistbox = BinListbox(self, (('ID', 7), ('# Tracks', 8), (chr(952), 9), (chr(966), 9), ('d_topo', 7)))

		self.grouplabel.grid(column=0, row=0, padx=2, pady=2, sticky='nw')
		self.grouplistbox.grid(column=0, row=1, sticky='ns')
		self.materialabel.grid(column=0, row=2, padx=2, pady=2, sticky='nw')
		self.materialistbox.grid(column=0, row=3, sticky='ns')
		self.binlabel.grid(column=0, row=4, padx=2, pady=2, sticky='nw')
		self.usedbinlistbox.grid(column=0, row=5, sticky='ns')

		self.rowconfigure(1, weight=1)
		self.rowconfigure(3, weight=1)
		self.rowconfigure(5, weight=1)

	def changeBinView(self):
		self.parent.changeBinView()

	def changeMatView(self):
		self.parent.changeMatView()

	def highlightSelection(self, row_idx):
		pass

class CentreWindow(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)
		self.__createlayout()

	def __createlayout(self):
		self.fig = plt.Figure(figsize=(10, 8))

		self.visLabel = tk.Label(self, text='Visualization')
		self.drawArea = FigureCanvasTkAgg(self.fig, self)

		self.visLabel.grid(column=0, row=0, sticky='w', padx=2, pady=2)
		self.drawArea.get_tk_widget().grid(column=0, row=1, sticky='nsew')

		self.rowconfigure(1, weight=1)
		self.columnconfigure(0, weight=1)


class BottomWindow(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)
		self.__createlayout()

	def __createlayout(self):
		self.buttonLabel = tk.Label(self, text='Actions')
		self.loadDataButton = tk.Button(self, text='Load data', command=self.parent.loadData)
		self.loadDemButton = tk.Button(self, text='Load DEM', command=self.parent.loadDEM)
		self.saveIntersButton = tk.Button(self, text='Save selected intersections', command=self.parent.saveIntersections)
		self.selectIntersectionButton = tk.Button(self, text='Select intersection', command=self.parent.intersection)
		self.addGroupButton = tk.Button(self, text='Add Group', command=self.parent.createGroup)
		self.delGroupButton = tk.Button(self, text='Delete Group(s)', command=self.parent.deleteGroup)
		self.addToGroupBtn = tk.Button(self, text='Add to group', command=self.parent.addToGroup)
		self.remFromGroupBtn = tk.Button(self, text='Remove from group', command=self.parent.remFromGroup)
		self.addMaterialBtn = tk.Button(self, text='Add Material', command=self.parent.addMaterial)
		self.remMaterialBtn = tk.Button(self, text='Remove Material', command=self.parent.remMaterial)
		self.saveModelBtn = tk.Button(self, text='Save Model', command=self.parent.saveModel)

		self.buttonLabel.grid(column=0, row=0, sticky='nw')
		self.loadDataButton.grid(column=0, row=1, sticky='w', padx=2, pady=2)
		self.loadDemButton.grid(column=1, row=1, sticky='w', padx=2, pady=2)
		self.saveIntersButton.grid(column=0, row=2, sticky='w', padx=2, pady=2)
		self.selectIntersectionButton.grid(column=3, row=1, padx=2, pady=2)
		self.addGroupButton.grid(column=4, row=1, sticky='e', padx=2, pady=2)
		self.delGroupButton.grid(column=4, row=2, sticky='e', padx=2, pady=2)
		self.addToGroupBtn.grid(column=6, row=1, sticky='e', padx=2, pady=2)
		self.remFromGroupBtn.grid(column=6, row=2, sticky='e', padx=2, pady=2)
		self.addMaterialBtn.grid(column=5, row=1, sticky='e', padx=2, pady=2)
		self.remMaterialBtn.grid(column=5, row=2, sticky='e', padx=2, pady=2)
		self.saveModelBtn.grid(column=3, row=3, padx=2, pady=2)

		self.columnconfigure(3, weight=1)


class MainWindow(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent)
		self.initUI()

		# Parameters
		self.groupList = []
		self.materialList = []

		self.headerLog = set()
		self.actHeader = pd.DataFrame()
		self.actData = pd.DataFrame()

		self.dem_dict = {}
		self.actDEM = np.zeros((1, 1))
		self.alreadyDrawn = False

		self.intersections = OrderedDict()
		self.highlight = 0

		self.nDataPoints = 0

		self.highlight = False

	def initUI(self):
		# Top level container
		self.leftframe = LeftWindow(self)
		self.centreframe = CentreWindow(self)
		self.rightframe = RightWindow(self)
		self.bottomframe = BottomWindow(self)

		# Layout top level frames
		self.leftframe.grid(row=0, column=0, rowspan=2, sticky='ns')
		self.centreframe.grid(row=0, column=1, sticky='nsew')
		self.bottomframe.grid(row=1, column=1, sticky='we')
		self.rightframe.grid(row=0, column=2, rowspan=2, sticky='ns')

		self.rowconfigure(0, weight=1)
		self.columnconfigure(1, weight=1)
		self.dataLoaded = False

	def redrawDEM(self):
		self.centreframe.fig.clear()

		X = self.X
		Y = self.Y

		ax = self.centreframe.fig.add_subplot(1, 1, 1, projection='3d')

		# ax.plot_surface(X, Y, self.actDEM, rstride=2, cstride=2)
		self.dem_lc = ax.plot_wireframe(X, Y, self.actDEM, rstride=10, cstride=10)
		# ax.plot_wireframe(X, Y, self.actDEM, rstride=1, cstride=1)
		# Create cubic bounding box to simulate equal aspect ratio
		max_range = np.nanmax(np.array([np.nanmax(X) - np.nanmin(X), np.nanmax(Y) - np.nanmin(Y),
		                                np.nanmax(self.actDEM) - np.nanmin(self.actDEM)]))
		Xb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].flatten() + 0.5 * (np.nanmax(X) + np.nanmin(X))
		Yb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].flatten() + 0.5 * (np.nanmax(Y) + np.nanmin(Y))
		Zb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][2].flatten() \
		     + 0.5 * (np.nanmax(self.actDEM) + np.nanmin(self.actDEM))

		# Comment or uncomment following both lines to test the fake bounding box:
		for xb, yb, zb in zip(Xb, Yb, Zb):
			ax.plot([xb], [yb], [zb], 'b')

		ax.grid()

		if self.dataLoaded:
			# Draw detector
			det_vec = np.array([float(self.actHeader.loc['E (m, CH1903)', 'value']),
			                    float(self.actHeader.loc['N (m, CH1903)', 'value']),
			                    float(self.actHeader.loc['Z (m, CH1903)', 'value'])])
			self.det_pc = ax.scatter(det_vec[0], det_vec[1], det_vec[2], s=10, c='magenta', marker='D')

		self.centreframe.drawArea.draw()

	def loadData(self):
		fileloc = tk.filedialog.askopenfilename(initialdir='./binnedData/')

		if fileloc:
			if fileloc.endswith('_aux.json'):
				mainfile = fileloc.strip('_aux.json') + '.json'

				try:
					auxfile = fileloc
					auxdf = pd.read_json(auxfile)
				except ValueError:
					messagebox.showerror('Error', fileloc + ' seems to be corrupted')
					return

				try:
					maindf = pd.read_json(mainfile,
										  dtype={'# Tracks': 'int64', chr(952): 'float64', chr(966): 'float64'})
				except ValueError:
					messagebox.showerror('Error', fileloc + ' does not have an accompanying main-file')
					return

			elif fileloc.endswith('.json'):
				auxfile = fileloc.strip('.json') + '_aux.json'
				try:
					maindf = pd.read_json(fileloc,
										  dtype={'# Tracks': 'int64', chr(952): 'float64', chr(966): 'float64'})
				except ValueError:
					messagebox.showerror('Error', fileloc + ' seems to be corrupted')
					return

				try:
					auxdf = pd.read_json(auxfile)
				except ValueError:
					messagebox.showerror('Error', fileloc + ' does not have an accompanying aux-file')
					return

			else:
				messagebox.showerror('Error', 'Check your json files in the binnedData folder!')
				return

			if self.alreadyDrawn:
				self.redrawDEM()

			self.headerLog.add(auxfile)

			self.leftframe.unbinlistbox.delete_all()

			self.actHeader = auxdf
			self.actData = maindf
			id = self.actHeader.loc['ID', 'value']
			self.actData.sort_index(inplace=True)
			self.highlight = 0
			self.intersections = OrderedDict()

			datalist = self.actData.values.tolist()

			indices = [self.nDataPoints + i for i in range(len(datalist))]

			self.nDataPoints += len(datalist)

			idlist = []
			for idx, row in enumerate(datalist):

				idlist.append(id + '_' + str(idx))
				inrow = [id + '_' + str(idx), int(row[0]), '{:,.2f}'.format(row[1]), '{:,.2f}'.format(row[2])]
				self.leftframe.unbinlistbox.insert(tk.END, inrow)
				self.leftframe.unbinlistbox.itemconfig(tk.END, foreground='grey')

			self.actData['d_topo'] = 0.
			self.actData['ID'] = idlist
			self.actData.reset_index(inplace=True)
			self.actData['index'] = indices
			self.actData['Materials'] = ''#[list() for x in range(len(self.actData.index))]

			self.dataLoaded = True

			# print(self.actData.dtypes)
			# print(self.actData)

	def loadDEM(self):

		self.highlight = 0
		self.intersections = OrderedDict()
		# Load DEM txt file
		fileloc = filedialog.askopenfilename(initialdir='./DEMs/')

		if fileloc:
			self.dem_dict = {}
			self.centreframe.fig.clear()
			for i, row in self.actData.iterrows():
				self.leftframe.unbinlistbox.itemconfig(i, foreground='grey')

			with open(fileloc, 'r') as file:
				# Readout header
				for i in range(6):
					a = (file.readline()).split()

					if a[0] == 'ncols' or a[0] == 'nrows' or a[0] == 'cellsize' or a[0] == 'NODATA_value':
						a[1] = int(a[1])
					else:
						a[1] = float(a[1])

					self.dem_dict[a[0]] = a[1]

				# Readout DEM in txt form
				lines = file.readlines()
				self.actDEM = np.zeros((self.dem_dict['nrows'], self.dem_dict['ncols']))

				for idl, line in enumerate(lines):
					self.actDEM[idl, :] = [float(i) for i in line.split()]


			# Mask NODATA values
			mask = self.actDEM == -9999
			self.actDEM[mask] = np.nan

			# Construct grid
			nx = self.dem_dict['cellsize'] * self.dem_dict['ncols']
			ny = self.dem_dict['cellsize'] * self.dem_dict['nrows']

			xend = self.dem_dict['xllcorner'] + nx
			yend = self.dem_dict['yllcorner'] + ny

			x = np.linspace(self.dem_dict['xllcorner'], xend, self.dem_dict['ncols'])
			y = np.linspace(yend, self.dem_dict['yllcorner'], self.dem_dict['nrows'])

			self.X, self.Y = np.meshgrid(x, y)
			X = self.X
			Y = self.Y

			ax = self.centreframe.fig.add_subplot(1, 1, 1, projection='3d')

			#ax.plot_surface(X, Y, self.actDEM, rstride=2, cstride=2)
			self.dem_lc = ax.plot_wireframe(X, Y, self.actDEM, rstride=10, cstride=10)
			#ax.plot_wireframe(X, Y, self.actDEM, rstride=1, cstride=1)
			# Create cubic bounding box to simulate equal aspect ratio
			max_range = np.nanmax(np.array([np.nanmax(X) - np.nanmin(X), np.nanmax(Y) - np.nanmin(Y),
								  np.nanmax(self.actDEM) - np.nanmin(self.actDEM)]))
			Xb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].flatten() + 0.5 * (np.nanmax(X) + np.nanmin(X))
			Yb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].flatten() + 0.5 * (np.nanmax(Y) + np.nanmin(Y))
			Zb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][2].flatten() \
				 + 0.5 * (np.nanmax(self.actDEM) + np.nanmin(self.actDEM))

			# Comment or uncomment following both lines to test the fake bounding box:
			for xb, yb, zb in zip(Xb, Yb, Zb):
				ax.plot([xb], [yb], [zb], 'b')

			ax.grid()

			if self.dataLoaded:
				# Draw detector
				det_vec = np.array([float(self.actHeader.loc['E (m, CH1903)', 'value']),
									float(self.actHeader.loc['N (m, CH1903)', 'value']),
									float(self.actHeader.loc['Z (m, CH1903)', 'value'])])
				self.det_pc = ax.scatter(det_vec[0], det_vec[1], det_vec[2], s=10, c='magenta', marker='D')

			self.centreframe.drawArea.draw()

			self.alreadyDrawn = True

	def intersection(self):

		self.leftframe.unbinlistbox.selection_clear(0, tk.END)

		def doIntersect(Pdem1, Pdem2, Pdem3, Pdem4, PDet, Dirc):
			# New algorithm
			ex = Pdem4 - Pdem3
			ey = Pdem1 - Pdem3

			normal = np.cross(ex, ey)
			a = np.dot(normal, Pdem3)

			lamb = ((a - normal[0]*PDet[0] - normal[1]*PDet[1] - normal[2]*PDet[2])
					/(normal[0]*Dirc[0] + normal[1]*Dirc[1] + normal[2]*Dirc[2]))

			Insect = PDet + lamb * Dirc

			if Pdem3[0] <= Insect[0] <= Pdem4[0] and Pdem3[1] <= Insect[1] <= Pdem1[1]:
				intersection = Insect
			else:
				intersection = np.array([])

			return intersection

		if not self.intersections:

			A = []
			B = []

			for ir in range(0, self.dem_dict['nrows'] - 1):
				for ic in range(0, self.dem_dict['ncols']-1):
					dx = self.X[ir, ic] - self.dem_dict['xllcorner']
					dy = self.Y[ir, ic] - self.dem_dict['yllcorner']
					# Paraboloid
					# a = [(self.X[ir, ic] - self.dem_dict['xllcorner'])**2,
					#      (self.Y[ir, ic] - self.dem_dict['yllcorner'])**2, 1]
					a = [dx**2, dy**2, dx*dy, dx, dy, 1]
					b = self.actDEM[ir, ic]
					if not np.isnan(b):
						A.append(a)
						B.append(b)

			A = np.array(A)
			B = np.array(B)

			ATA = np.matmul(A.transpose(), A)
			ATB = np.matmul(A.transpose(), B)

			sol = np.matmul(np.linalg.inv(ATA), ATB)

			DX = self.X - self.dem_dict['xllcorner']
			DY = self.Y - self.dem_dict['yllcorner']

			ZPlane = (  sol[0]*DX**2
					  + sol[1]*DY**2
					  + sol[2]*DX*DY
					  + sol[3]*DX
					  + sol[4]*DY
					  + sol[5])

			# Paraboloid
			# ZPlane = solution[0]*(self.X - self.dem_dict['xllcorner'])**2 \
			#          + solution[1]*(self.Y - self.dem_dict['yllcorner'])**2 + solution[2]

			ax = self.centreframe.fig.get_axes()[0]
			#ax.plot_wireframe(self.X, self.Y, ZPlane, rstride=10, cstride=10, color='orange')
			self.centreframe.drawArea.draw()

			det_vec = np.array([float(self.actHeader.loc['E (m, CH1903)', 'value']),
					   float(self.actHeader.loc['N (m, CH1903)', 'value']),
					   float(self.actHeader.loc['Z (m, CH1903)', 'value'])])

			angles = zip(self.actData['index'], self.actData[chr(952)], self.actData[chr(966)])

			intersections = []

			MAX = self.actData.shape[0]
			progVar = tk.DoubleVar()
			progBar = ttk.Progressbar(self.centreframe, length=100, variable=progVar)
			progBar.grid(column=0, row=0, sticky='e', padx=2, pady=2)


			for pbidx, theta, phi in angles:
				dfidx = self.actData.loc[self.actData['index'] == pbidx].index[0]
				found = False

				# Prepare variables
				t = np.deg2rad(theta)
				p = np.deg2rad(phi)

				xdir = np.sin(t) * np.sin(p)
				ydir = np.sin(t) * np.cos(p)
				zdir = np.cos(t)

				xp = det_vec[0] - self.dem_dict['xllcorner']
				yp = det_vec[1] - self.dem_dict['yllcorner']
				zp = det_vec[2]

				dir_vec = np.array([xdir, ydir, zdir])

				# Quadratic surface (c2 - quadratic coefficient, c1 - linear coefficient, c0 - constant coefficient)
				c2 = sol[0]*xdir**2+sol[1]*ydir**2+sol[2]*xdir*ydir
				c1 = 2*sol[0]*xp*xdir + 2*sol[1]*yp*ydir + sol[2]*(xp*ydir + yp*xdir) + sol[3]*xdir + sol[4]*ydir - zdir
				c0 = sol[0]*xp**2 + sol[1]*yp**2 + sol[2]*xp*yp + sol[3]*xp + sol[4]*yp + sol[5] - zp

				lamb1 = (-c1 + np.sqrt(c1**2 - 4*c2*c0))/(2*c2)
				lamb2 = (-c1 - np.sqrt(c1**2 - 4*c2*c0))/(2*c2)

				if lamb1 > 0 and lamb2 > 0:
					lambd = min(lamb1, lamb2)
				elif lamb1 > lamb2:
					lambd = lamb1
				elif lamb2 > lamb1:
					lambd = lamb2
				else:
					print('problem Oo')

				xinter = det_vec[0] + lambd*dir_vec[0]
				yinter = det_vec[1] + lambd*dir_vec[1]
				zinter = det_vec[2] + lambd*dir_vec[2]

				xmin = self.dem_dict['xllcorner']
				ymin = self.dem_dict['yllcorner']
				xmax = xmin + self.dem_dict['ncols'] * self.dem_dict['cellsize']
				ymax = ymin + self.dem_dict['nrows'] * self.dem_dict['cellsize']

				# Height correction only if point within DEM
				if xmin <= xinter <= xmax and ymin <= yinter <= ymax:

					ni = round((xinter-xmin)/self.dem_dict['cellsize'])
					nj = round((ymax-yinter)/self.dem_dict['cellsize'])

					buffer = 40
					# Determine search indices
					nimax = ni + buffer
					if nimax > self.dem_dict['ncols']-1 :
						nimax = self.dem_dict['ncols']-1
					nimin = ni - buffer
					if nimin < 0:
						nimin = 0

					njmax = nj + buffer
					if njmax > self.dem_dict['nrows']-1:
						njmax = self.dem_dict['nrows']-1
					njmin = nj - buffer
					if njmin < 0:
						njmin = 0

					z_mean = 0
					cntr = 0
					# for ir in range(int(njmin), int(njmax)):
					#     for ic in range(int(nimin), int(nimax)):
					#         if self.actDEM is not np.nan:
					#             z_mean += self.actDEM[ir, ic]
					#             cntr += 1
					njmin = int(njmin)
					nimin = int(nimin)
					njmax = int(njmax)
					nimax = int(nimax)

					if not np.isnan(self.actDEM[njmin, nimin]):
						z_mean += self.actDEM[njmin, nimin]
						cntr += 1
					if not np.isnan(self.actDEM[njmin, nimax]):
						z_mean += self.actDEM[njmin, nimax]
						cntr += 1
					if not np.isnan(self.actDEM[njmax, nimin]):
						z_mean += self.actDEM[njmax, nimin]
						cntr += 1
					if not np.isnan(self.actDEM[njmax, nimax]):
						z_mean += self.actDEM[njmax, nimax]
						cntr += 1


					if z_mean == 0:
						xcorr = xinter
						ycorr = yinter
						zcorr = zinter
					else:
						zcorr = z_mean / cntr
						alpha = (zcorr - det_vec[2]) / zdir
						xcorr = det_vec[0] + alpha * xdir
						ycorr = det_vec[1] + alpha * ydir
				else:
					xcorr = xinter
					ycorr = yinter
					zcorr = zinter

				# Search buffer (Snakelike)
				if xmin <= xcorr <= xmax and ymin <= ycorr <= ymax:
					stringy = 'Buffer'
					# Get indices of intersection point (height corrected)
					ni = int(round((xcorr - xmin) / self.dem_dict['cellsize']))
					nj = int(round((ymax - ycorr) / self.dem_dict['cellsize']))

					found = False

					# Initial check
					if (np.isnan(self.actDEM[nj, ni]) or np.isnan(self.actDEM[nj + 1, ni])
						or np.isnan(self.actDEM[nj, ni + 1]) or np.isnan(self.actDEM[nj + 1, ni + 1])):
						pass
					else:
						# upper left corner
						ulc = np.array([self.X[nj, ni], self.Y[nj, ni], self.actDEM[nj, ni]])
						# upper right corner
						urc = np.array([self.X[nj, ni + 1], self.Y[nj, ni + 1], self.actDEM[nj, ni + 1]])
						# lower left corner
						llc = np.array([self.X[nj + 1, ni], self.Y[nj + 1, ni], self.actDEM[nj + 1, ni]])
						# lower right corner
						lrc = np.array([self.X[nj + 1, ni + 1], self.Y[nj + 1, ni + 1], self.actDEM[nj + 1, ni + 1]])

						intersection = doIntersect(ulc, urc, llc, lrc, det_vec, dir_vec)

						if intersection.size != 0:
							# print(ir, ic)
							intersections.append(intersection)
							dist = ((intersection[0] - det_vec[0]) ** 2
									+ (intersection[1] - det_vec[1]) ** 2
									+ (intersection[2] - det_vec[2]) ** 2) ** 0.5
							self.actData.at[dfidx, 'd_topo'] = dist
							found = True

					nctr = 0
					rn = 1
					while not found:
						tn = np.pi/2*nctr
						for ri in range(rn):
							njpass = False
							nipass = False
							nj += int((-np.sin(tn)))
							ni += int(np.cos(tn))
							if 0 <= nj < self.dem_dict['nrows']-1:
								njpass = True
							if 0 <= ni < self.dem_dict['ncols']-1:
								nipass = True

							if njpass and nipass:
								if (np.isnan(self.actDEM[nj, ni]) or np.isnan(self.actDEM[nj + 1, ni])
										or np.isnan(self.actDEM[nj, ni + 1]) or np.isnan(self.actDEM[nj + 1, ni + 1])):
									pass
								else:
									# upper left corner
									ulc = np.array([self.X[nj, ni], self.Y[nj, ni], self.actDEM[nj, ni]])
									# upper right corner
									urc = np.array([self.X[nj, ni + 1], self.Y[nj, ni + 1], self.actDEM[nj, ni + 1]])
									# lower left corner
									llc = np.array([self.X[nj + 1, ni], self.Y[nj + 1, ni], self.actDEM[nj + 1, ni]])
									# lower right corner
									lrc = np.array(
										[self.X[nj + 1, ni + 1], self.Y[nj + 1, ni + 1], self.actDEM[nj + 1, ni + 1]])

									intersection = doIntersect(ulc, urc, llc, lrc, det_vec, dir_vec)

									if intersection.size != 0:
										intersections.append(intersection)
										self.intersections[pbidx] = intersection
										dist = ((intersection[0] - det_vec[0]) ** 2
												+ (intersection[1] - det_vec[1]) ** 2
												+ (intersection[2] - det_vec[2]) ** 2) ** 0.5
										self.actData.at[dfidx, 'd_topo'] = dist
										found = True
										break

						if found:
							self.leftframe.unbinlistbox.selection_set(dfidx)

						# Update
						nctr += 1
						# Security check
						if nctr > 1000:
							#print('no°:', pbidx, ' had to break')
							break

						if nctr % 2 == 0: # if even
							rn += 1

				else: # Search whole DEM
					stringy = 'Whole'
					nimin = 0
					nimax = self.dem_dict['ncols']-1
					njmin = 0
					njmax = self.dem_dict['nrows']-1

					# Intersection loop
					for nj in range(int(njmin), int(njmax)):
						found = False
						for ni in range(int(nimin), int(nimax)):
							# upper left corner
							ulc = np.array([self.X[nj, ni], self.Y[nj, ni], self.actDEM[nj, ni]])
							if np.isnan(ulc).any():
								continue
							# upper right corner
							urc = np.array([self.X[nj, ni+1], self.Y[nj, ni+1], self.actDEM[nj, ni+1]])
							if np.isnan(urc).any():
								continue
							# lower left corner
							llc = np.array([self.X[nj+1, ni], self.Y[nj+1, ni], self.actDEM[nj+1, ni]])
							if np.isnan(llc).any():
								continue
							# lower right corner
							lrc = np.array([self.X[nj+1, ni+1], self.Y[nj+1, ni+1], self.actDEM[nj+1, ni+1]])
							if np.isnan(lrc).any():
								continue
							intersection = doIntersect(ulc, urc, llc, lrc, det_vec, dir_vec)

							if intersection.size != 0:
								intersections.append(intersection)
								self.intersections[pbidx] = intersection
								dist = ((intersection[0] - det_vec[0])**2
										+ (intersection[1] - det_vec[1])**2
										+ (intersection[2] - det_vec[2])**2)**0.5
								self.actData.at[dfidx, 'd_topo'] = dist
								found = True
								break

						if found:
							self.leftframe.unbinlistbox.selection_set(dfidx)
							break

				if found:
					#print('no° ', pbidx, ': ', nj, ni, ' | ', stringy, ' | ', intersection, ' | nctr: ', nctr)
					print(self.actData)
					print(pbidx)
					self.leftframe.unbinlistbox.itemconfig(dfidx, foreground='black')
				else:
					#print('no° ', pbidx, ': ', nj, ni, ' | ', stringy, ' | ', 'nothing', ' | nctr: ', nctr)
					pass

				progVar.set(float(pbidx) / float(MAX) * 100)
				self.update_idletasks()

				# pla_inter.append([xinter, yinter, zinter])
				# if xmin <= xcorr <= xmax and ymin <= ycorr <= ymax:
				#     print(xcorr, ycorr, zcorr)
				#     pla_inter.append([xcorr, ycorr, zcorr])

			progBar.destroy()
			# print(self.dem_dict)
			# pla_inter = np.array(pla_inter)

			print(self.intersections)

			itsc = np.array(intersections)
			ax = self.centreframe.fig.get_axes()[0]
			self.inters_pc = ax.scatter(itsc[:, 0], itsc[:, 1], itsc[:, 2], c='r', depthshade=False)
			#ax.scatter(pla_inter[:, 0], pla_inter[:, 1], pla_inter[:, 2], c='magenta', depthshade=False)

			self.centreframe.drawArea.draw()

		else:
			# TODO: Implement security if after selection, one item has moved to attached bins
			for key, value in self.intersections.items():
				index = self.actData[self.actData['index'] == key].index[0]
				self.leftframe.unbinlistbox.selection_set(index)

		self.highlightSelection(self.leftframe.unbinlistbox.curselection())

		# # Slow but works
		#
		# self.intersections = {}
		#
		# det_vec = np.array([float(self.actHeader.loc['E (m, CH1903)', 'value']),
		#            float(self.actHeader.loc['N (m, CH1903)', 'value']),
		#            float(self.actHeader.loc['Z (m, CH1903)', 'value'])])
		#
		# angles = zip(self.actData[chr(952)], self.actData[chr(966)])
		#
		# print(self.actHeader)
		#
		# print('Detector vector: ', det_vec)
		#
		# print(self.dem_dict)
		#
		#
		# MAX = self.actData.shape[0]
		# progVar = tk.DoubleVar()
		# progBar = ttk.Progressbar(self.centreframe, length=100, variable=progVar)
		# progBar.grid(column=0, row=0, sticky='e', padx=2, pady=2)
		# idx = 0
		#
		#
		# intersections = []
		# for theta, phi in angles:
		#     t = np.deg2rad(theta)
		#     p = np.deg2rad(phi)
		#
		#     xdir = np.sin(t)*np.sin(p)
		#     ydir = np.sin(t) * np.cos(p)
		#     zdir = np.cos(t)
		#
		#     dir_vec = np.array([xdir, ydir, zdir])
		#
		#     for ir in range(0, self.dem_dict['nrows']-1):
		#         find = False
		#         for ic in range(0, self.dem_dict['ncols']-1):
		#             # upper left corner
		#             ulc = np.array([self.X[ir, ic], self.Y[ir, ic], self.actDEM[ir, ic]])
		#             # upper right corner
		#             urc = np.array([self.X[ir, ic+1], self.Y[ir, ic+1], self.actDEM[ir, ic+1]])
		#             # lower left corner
		#             llc = np.array([self.X[ir+1, ic], self.Y[ir+1, ic], self.actDEM[ir+1, ic]])
		#             # lower right corner
		#             lrc = np.array([self.X[ir+1, ic+1], self.Y[ir+1, ic+1], self.actDEM[ir+1, ic+1]])
		#
		#             if not(np.nan in ulc or np.nan in urc or np.nan in llc or np.nan in lrc):
		#                 intersection = doIntersect(ulc, urc, llc, lrc, det_vec, dir_vec)
		#
		#                 if intersection.size != 0:
		#                     print(ir, ic)
		#                     intersections.append(intersection)
		#                     self.intersections[idx] = len(intersection)
		#                     find = True
		#                     break
		#         if find:
		#             self.leftframe.unbinlistbox.selection_set(idx)
		#             break
		#     idx += 1
		#
		#     progVar.set(float(idx) / float(MAX) * 100)
		#     self.update_idletasks()
		#
		# itsc = np.array(intersections)
		# ax = self.centreframe.fig.get_axes()[0]
		# ax.scatter(itsc[:, 0], itsc[:, 1], itsc[:, 2], c='r')
		# self.highlight = ax.scatter([], [], [], color='yellow')
		#
		# self.centreframe.drawArea.draw()
		#
		# progBar.destroy()

	def createGroup(self):

		def acceptEntries(event=None):
			entry = gnentry.get()
			if entry:
				self.rightframe.grouplistbox.insert(tk.END, [entry, 0])
				self.groupList.append(pd.DataFrame())
				self.materialList.append([])
			groupWindow.destroy()

		x = self.winfo_rootx()#screenwidth()
		y = self.winfo_rooty()#screenheight()

		ws = self.winfo_width()
		hs = self.winfo_height()

		groupWindow = tk.Toplevel(self)
		groupWindow.title('Create Group')
		groupWindow.transient(self)
		groupWindow.geometry('+%d+%d'%(x + ws/2, y + hs/2))
		groupWindow.bind('<Return>', acceptEntries)

		gnlbl = tk.Label(groupWindow, text='Name: ')
		gnentry = tk.Entry(groupWindow)
		accButton = tk.Button(groupWindow, text='Accept', command=acceptEntries)

		gnlbl.grid(column=0, row=0, sticky='w', padx=5, pady=5)
		gnentry.grid(column=1, row=0, sticky='e', padx=5, pady=5)
		accButton.grid(column=0, row=1, sticky='w', padx=5, pady=5)

		groupWindow.grab_set()
		groupWindow.focus()
		gnentry.focus_set()

		self.wait_window(groupWindow)
		self.rightframe.grouplistbox.selection_clear(0, tk.END)
		self.rightframe.grouplistbox.selection_set(tk.END)

	def deleteGroup(self):
		sel = self.rightframe.grouplistbox.curselection()

		oldsize = self.rightframe.grouplistbox.size()

		for item in reversed(sel):
			self.rightframe.usedbinlistbox.select_all()
			self.remFromGroup()
			self.rightframe.materialistbox.select_all()
			self.remMaterial()
			self.rightframe.grouplistbox.delete(item)
			del self.groupList[item]
			del self.materialList[item]


		if self.rightframe.grouplistbox.size() != 0:
			if sel[0] == oldsize-1:
				self.rightframe.grouplistbox.selection_set(oldsize-2)
			else:
				self.rightframe.grouplistbox.selection_set(sel[0])


	def addToGroup(self):

		grsel = self.rightframe.grouplistbox.curselection()

		if len(grsel) != 1:
			messagebox.showerror('Error', 'Please select 1 group')
		else:
			grselidx = grsel[0]
			sel = self.leftframe.unbinlistbox.curselection()

			rows = self.actData.loc[sel, :]

			self.groupList[grselidx] = self.groupList[grselidx].append(rows, ignore_index=True)
			self.actData.drop(rows.index, inplace=True)
			self.actData.reset_index(drop=True, inplace=True)

			for item in reversed(sel):
				self.leftframe.unbinlistbox.delete(item)

			self.rightframe.usedbinlistbox.delete_all()
			datalist = self.groupList[grselidx].values.tolist()

			for idx, row in enumerate(datalist):
				inrow = ([row[5], int(row[1]), '{:,.2f}'.format(row[2]),
						  '{:,.2f}'.format(row[3]), '{:,.2f}'.format(row[4])])
				self.rightframe.usedbinlistbox.insert(tk.END, inrow)

			self.rightframe.grouplistbox.setSize(self.rightframe.usedbinlistbox.size())

			self.rightframe.usedbinlistbox.sort(self.groupList[grselidx]['index'].to_list())

			self.groupList[grselidx].sort_values(by=['index'], inplace=True)
			self.groupList[grselidx].reset_index(drop=True, inplace=True)

	def remFromGroup(self):

		grsel = self.rightframe.grouplistbox.curselection()

		if len(grsel) != 1:
			messagebox.showerror('Error', 'Please select 1 group')
		else:
			grselidx = grsel[0]
			if not self.groupList[grselidx].empty:
				sel = self.rightframe.usedbinlistbox.curselection()

				rows = self.groupList[grselidx].loc[sel, :]

				self.actData = self.actData.append(rows, ignore_index=True)
				self.groupList[grselidx].drop(rows.index, inplace=True)
				self.groupList[grselidx].reset_index(drop=True, inplace=True)

				for item in reversed(sel):
					self.rightframe.usedbinlistbox.delete(item)

				self.leftframe.unbinlistbox.delete_all()
				datalist = self.actData.values.tolist()
				for idx, row in enumerate(datalist):
					inrow = ([row[5], int(row[1]), '{:,.2f}'.format(row[2]),
							  '{:,.2f}'.format(row[3]), '{:,.2f}'.format(row[4])])
					self.leftframe.unbinlistbox.insert(tk.END, inrow)

				self.rightframe.grouplistbox.setSize(self.rightframe.usedbinlistbox.size())

				self.leftframe.unbinlistbox.sort(self.actData['index'].to_list())

				self.actData.sort_values(by=['index'], inplace=True)
				self.actData.reset_index(drop=True, inplace=True)

	def changeBinView(self):
		grselidx = self.rightframe.grouplistbox.curselection()[0]

		datalist = self.groupList[grselidx].values.tolist()
		self.rightframe.usedbinlistbox.delete_all()
		for idx, row in enumerate(datalist):
			inrow = ([row[5], int(row[1]), '{:,.2f}'.format(row[2]),
					  '{:,.2f}'.format(row[3]), '{:,.2f}'.format(row[4])])
			self.rightframe.usedbinlistbox.insert(tk.END, inrow)

	def changeMatView(self):
		grselidx = self.rightframe.grouplistbox.curselection()[0]

		datalist = self.materialList[grselidx]
		self.rightframe.materialistbox.delete_all()
		for row in datalist:
			self.rightframe.materialistbox.insert(tk.END, [row])


	def highlightSelection(self, row_idx):
		if self.intersections:
			globidx = []
			for ri in row_idx:
				globidx.append(self.actData.at[ri, 'index'])
			# globidx = self.actData.at[row_idx, 'index']

			check = False
			for item in globidx:
				if item in self.intersections:
					check = True

			#if globidx in self.intersections:
			if check:
				self.inters_pc.set_visible(False)

				col = []
				rcolor = (1,0,0,1)
				gcolor = (0,1,0,1)
				colormap = np.array([rcolor, gcolor])
				x = []
				y = []
				z = []

				for i in self.intersections:
					if i in globidx:
						col.append(1)
					else:
						col.append(0)
					x.append(self.intersections[i][0])
					y.append(self.intersections[i][1])
					z.append(self.intersections[i][2])

				ax = self.centreframe.fig.get_axes()[0]
				if self.highlight:
					self.highlight.remove()
				self.highlight = ax.scatter(x, y, z, c=colormap[col], depthshade=False)
				self.highlight.set_visible(True)
			else:
				if self.highlight:
					self.highlight.set_visible(False)
				self.inters_pc.set_visible(True)

		self.centreframe.drawArea.draw()

	def saveIntersections(self):
		if self.intersections:
			selected = self.leftframe.unbinlistbox.curselection()

			globidx = []

			tosave = []
			for ri in selected:
				globidx.append(self.actData.at[ri, 'index'])

			for idx in self.intersections:
				if idx in globidx:
					tosave.append([self.intersections[idx][0], self.intersections[idx][1], self.intersections[idx][2]])

			savefile = filedialog.asksaveasfilename()

			savedf = pd.DataFrame(tosave, columns=['x', 'y', 'z'])

			if savefile:
				with open(savefile + '.json', 'w') as fp:
					json.dump(savedf, fp, cls=JSONEncoder)

	def addMaterial(self):
		matdir = filedialog.askdirectory(initialdir='./Materials')

		if matdir:
			matname = os.path.split(matdir)[1]
			self.rightframe.materialistbox.insert(tk.END, [matname])

			grselidx = self.rightframe.grouplistbox.curselection()[0]
			self.materialList[grselidx].append(matname)



	def remMaterial(self):
		sel = self.rightframe.materialistbox.curselection()

		for item in reversed(sel):
			self.rightframe.materialistbox.delete(item)

			grselidx = self.rightframe.grouplistbox.curselection()[0]
			del self.materialList[grselidx][item]


	def saveModel(self):
		finaldf = pd.DataFrame()

		totmatlist = set()

		for i, idf in enumerate(self.groupList):
			stringy = ''
			for mati in self.materialList[i]:
				totmatlist.add(mati)
				stringy += mati
				stringy += ' '
			idf['Materials'] = stringy

			finaldf = finaldf.append(idf, ignore_index=True)

		totmatlist = list(totmatlist)


		# Query Model name
		newfolder = filedialog.asksaveasfilename(initialdir='./Models')
		rootdir = os.path.split(os.path.split(newfolder)[0])[0]

		# Create Model directory
		os.mkdir(newfolder)

		newmatpath = os.path.normpath(os.path.join(newfolder, 'Materials'))
		os.mkdir(newmatpath)

		oldmatpath = os.path.normpath(os.path.join(rootdir, 'Materials'))

		# Copy Material files
		for mati in totmatlist:
			matipath = os.path.normpath(os.path.join(newmatpath, mati))
			oldmatipath = os.path.normpath(os.path.join(oldmatpath, mati))

			copy_tree(oldmatipath, matipath)

		# Data folder
		datapath = os.path.normpath(os.path.join(newfolder, 'Data'))
		os.mkdir(datapath)

		# Copy header data
		for headi in list(self.headerLog):
			copy2(headi, datapath)

		# Output folder
		outpath = os.path.normpath(os.path.join(newfolder, 'Output'))
		os.mkdir(outpath)

		# Dump json data
		datafile = os.path.normpath(os.path.join(datapath, 'data.json'))
		finaldf.to_json(datafile)


def main():
	root = tk.Tk()
	root.title("Model Builder")
	screen_width = root.winfo_screenwidth()
	screen_height = root.winfo_screenheight()
	sizestring = str(int(0.8*screen_width)) + 'x' + str(int(0.8*screen_height)) + '+50+50'
	root.geometry(sizestring)

	window = MainWindow(master=root)
	window.pack(side="top", fill="both", expand=True)
	root.mainloop()


if __name__ == '__main__':
	main()

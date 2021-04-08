import tkinter as tk
from mpl_toolkits import mplot3d
from collections import OrderedDict
from tkinter import filedialog
import numpy as np
import scipy as sp
import scipy.stats as sps
import scipy.special as spsp
import scipy.ndimage as spn
import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import json
import os


class JSONEncoder(json.JSONEncoder):
	def default(self, obj):
		if hasattr(obj, 'to_json'):
			return obj.to_json(orient='records')
		return json.JSONEncoder.default(self, obj)


class RightWindow(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)

		self.cellsize = tk.DoubleVar()
		self.smoothingval = tk.IntVar()
		self.sigma = tk.DoubleVar()
		self.dampingval = tk.DoubleVar()
		self.smoothafterdamp = tk.BooleanVar()

		self.colors = ['orange', 'g', 'r', 'm', 'y', 'k']
		self.coloriter = iter(self.colors)

		self.sigma.set(0)
		self.cellsize.set(20)
		self.smoothingval.set(0)
		self.dampingval.set(4)
		self.smoothafterdamp.set(False)

		self.interpoldict = {}
		self.fixdict = {}

		self.fixdf = pd.DataFrame()
		self.interpdf = pd.DataFrame()

		self.result = 0
		self.ncols = 0
		self.nrows = 0
		self.xllcorner = 0
		self.yllcorner = 0
		self.cellsize_forsave = 0
		self.nodatavalue = -9999

		self.__createlayout()

	def __createlayout(self):
		self.toplabel = tk.Label(self, text='Interpolated length group')
		self.interpollist = SingleListbox(self, ('Name', 14))
		self.loadinterpolbtn = tk.Button(self, text='Load', command=self.loadinterp)
		self.cellsizetext = tk.Label(self, text='Cellsize')
		self.cellsizeety = tk.Entry(self, textvariable=self.cellsize)
		self.interpoltext = tk.Label(self, text='Cell smoothing | # cells')
		self.interpolentry = tk.Entry(self, textvariable=self.smoothingval)
		self.sigmatext = tk.Label(self, text='+/- sigma (0: mean)')
		self.sigmaentry = tk.Entry(self, textvariable=self.sigma)

		self.botlabel = tk.Label(self, text='Fixed length group')
		self.fixedlist = SingleListbox(self, ('Name', 14))
		self.loadfixedbtn = tk.Button(self, text='Load', command=self.loadfixed)

		self.clearbtn = tk.Button(self, text='Clear all', command=self.clearall)

		self.dampingtxt = tk.Label(self, text='Damping weight fixed:')
		self.dampingval = tk.Entry(self, textvariable=self.dampingval)
		self.dampsmoothcb = tk.Checkbutton(self, text='Damping before Smoothing?', variable=self.smoothafterdamp)

		self.showbtn = tk.Button(self, text='Show surface', command=self.calculate)
		self.savebtn = tk.Button(self, text='Save surface as DEM', command=self.savedem)

		self.toplabel.grid(column=0, row=0, padx=2, pady=2, sticky='w')
		self.interpollist.grid(column=0, row=1, columnspan=2, padx=2, pady=2, sticky='nsew')
		self.loadinterpolbtn.grid(column=0, row=2, padx=2, pady=2)
		self.cellsizetext.grid(column=0, row=3, padx=2, pady=2)
		self.cellsizeety.grid(column=1, row=3, padx=2, pady=2)
		self.interpoltext.grid(column=0, row=4, padx=2, pady=2)
		self.interpolentry.grid(column=1, row=4, padx=2, pady=2)
		self.sigmatext.grid(column=0, row=5, padx=2, pady=2)
		self.sigmaentry.grid(column=1, row=5, padx=2, pady=2)

		self.botlabel.grid(column=0, row=6, padx=2, pady=2, sticky='w')
		self.fixedlist.grid(column=0, row=7, columnspan=2, padx=2, pady=2, sticky='nsew')
		self.loadfixedbtn.grid(column=0, row=8, padx=2, pady=2)

		self.clearbtn.grid(column=1, row=8, padx=2, pady=2)

		self.dampingtxt.grid(column=0, row=9, padx=2, pady=2)
		self.dampingval.grid(column=1, row=9, padx=2, pady=2)
		self.dampsmoothcb.grid(column=0, row=10, columnspan=2, padx=2, pady=2, sticky='w')

		self.showbtn.grid(column=0, row=11, padx=2, pady=20)
		self.savebtn.grid(column=1, row=11, padx=2, pady=20)

		self.rowconfigure(1, weight=1)
		self.rowconfigure(7, weight=1)

	def loadinterp(self):
		files = filedialog.askopenfilenames()

		if files:
			for file in files:
				listname = os.path.basename(file).strip('.json')

				jsonfile = json.load(open(file))

				self.interpollist.insert([listname])
				self.interpoldict[listname] = jsonfile

			# 	print(jsonfile)
			#
			# print(self.interpoldict)

	def loadfixed(self):
		files = filedialog.askopenfilenames()

		if files:
			for file in files:
				listname = os.path.basename(file).strip('.json')

				jsonfile = json.load(open(file))

				self.fixedlist.insert([listname])
				self.fixdict[listname] = jsonfile

			# 	print(jsonfile)
			#
			# print(self.fixdict)

	def __smoothing_calc(self, npx_size):
		dim = npx_size * 2 + 1
		generator = [spsp.binom(dim - 1, i) for i in range(dim)]
		smooth = np.outer(generator, generator)

		smooth /= np.sum(np.concatenate(smooth))

		return smooth

	def smoothing(self, smoothval, smoothing_data):
		smoothing_kernel = self.__smoothing_calc(smoothval)

		result_smoothed_wl = spn.convolve(smoothing_data, smoothing_kernel)
		V = smoothing_data.copy()
		V[np.isnan(smoothing_data)] = 0
		VV = spn.convolve(V, smoothing_kernel)

		W = 0 * smoothing_data.copy() + 1
		W[np.isnan(smoothing_data)] = 0
		WW = spn.convolve(W, smoothing_kernel)

		Z = VV/WW
		result_smoothed = Z.copy()
		result_smoothed[np.isnan(smoothing_data)]=np.nan

		return result_smoothed

	def damping(self, weight, var_dem, fix_dem):
		weight = float(weight)
		Q = var_dem.copy()
		Q[np.isnan(var_dem)] = 0
		WQ = 0 * var_dem.copy() + 1
		WQ[np.isnan(var_dem)] = 0

		P = weight * fix_dem.copy()
		P[np.isnan(fix_dem)] = 0
		WP = 0 * fix_dem.copy() + weight
		WP[np.isnan(fix_dem)] = 0

		result = (Q + P) / (WQ + WP)

		return result

	def calculate(self):
		# Set up grid
		xmin = 0
		xmax = 0
		ymin = 0
		ymax = 0

		# Load Dataframes
		fixdf = pd.DataFrame(columns=['x', 'y', 'z'], dtype=float)
		interpdf = pd.DataFrame(columns=['x', 'y', 'z'], dtype=float)

		for jf in self.fixdict.values():
			fixdf = fixdf.append(pd.read_json(jf), ignore_index=True)

		print(fixdf['z'])

		# Check global boundary values
		for index, row in fixdf.iterrows():
			if xmin == 0 or row['x'].min() < xmin:
				xmin = row['x'].min()
			if ymin == 0 or row['y'].min() < ymin:
				ymin = row['y'].min()
			if row['x'].max() > xmax:
				xmax = row['x'].max()
			if row['y'].max() > ymax:
				ymax = row['y'].max()

		# Load Dataframes and calculate statistic
		sigma = self.sigma.get()
		qu = sps.norm.cdf(sigma)
		for jf in self.interpoldict.values():
			for df in jf.values():
				dfs = pd.read_json(df)

				if sigma == 0:
					interpdf = interpdf.append(dfs.mean(), ignore_index=True)
				else:
					interpdf = interpdf.append(dfs.quantile(qu), ignore_index=True)

		# Check global boundary values
		for index, row in interpdf.iterrows():
			if xmin == 0 or row['x'].min() < xmin:
				xmin = row['x'].min()
			if ymin == 0 or row['y'].min() < ymin:
				ymin = row['y'].min()
			if row['x'].max() > xmax:
				xmax = row['x'].max()
			if row['y'].max() > ymax:
				ymax = row['y'].max()

		# Create Grid
		DX = xmax - xmin
		DY = ymax - ymin

		cellsize = self.cellsize.get()

		nrows = int(np.ceil(DY / cellsize))
		ncols = int(np.ceil(DX / cellsize))

		xllcorner = xmin - (ncols * cellsize - DX) / 2
		yllcorner = ymin - (nrows * cellsize - DY) / 2

		# Create Interpolation Grid
		height_sum = np.zeros((nrows + 1, ncols + 1))
		wt_sum = np.zeros((nrows + 1, ncols + 1))
		result_fix = np.zeros((nrows + 1, ncols + 1))

		cellareafix = cellsize**2

		# Interpolation of fixed points (usually bedrock)
		for index, row in fixdf.iterrows():
			# Get x-index
			xindex = int((row['x'] - xllcorner) // cellsize)

			# Get y-index
			yindex_dem = int((row['y'] - yllcorner) // cellsize)
			yindex_arr = nrows - int((row['y'] - yllcorner) // cellsize) - 1

			dxl = row['x'] - (xindex * cellsize + xllcorner)
			dxu = cellsize - dxl
			dyl = row['y'] - (yindex_dem * cellsize + yllcorner)
			dyu = cellsize - dyl

			wt_00 = dxu * dyu / cellareafix
			wt_p0 = dxl * dyu / cellareafix
			wt_0p = dxu * dyl / cellareafix
			wt_pp = dxl * dyl / cellareafix

			height_sum[yindex_arr, xindex] += row['z'] * wt_00
			height_sum[yindex_arr, xindex + 1] += row['z'] * wt_p0
			height_sum[yindex_arr + 1, xindex] += row['z'] * wt_0p
			height_sum[yindex_arr + 1, xindex + 1] += row['z'] * wt_pp

			wt_sum[yindex_arr, xindex] += wt_00
			wt_sum[yindex_arr, xindex + 1] += wt_p0
			wt_sum[yindex_arr + 1, xindex] += wt_0p
			wt_sum[yindex_arr + 1, xindex + 1] += wt_pp

		for i in range(nrows + 1):
			for j in range(ncols + 1):
				if height_sum[i, j] == 0 or wt_sum[i, j] == 0:
					result_fix[i, j] = np.nan
				else:
					result_fix[i, j] = height_sum[i, j] / wt_sum[i, j]

		# Interpolation for smoothed points (usually glacier)
		smoothing_kernel = self.__smoothing_calc(self.smoothingval.get())

		height_sum = np.zeros((nrows + 1, ncols + 1))
		wt_sum = np.zeros((nrows + 1, ncols + 1))
		result_smooth = np.zeros((nrows + 1, ncols + 1))

		for index, row in interpdf.iterrows():
			# Get x-index
			xindex = int((row['x'] - xllcorner) // cellsize)

			# Get y-index
			yindex_dem = int((row['y'] - yllcorner) // cellsize)
			yindex_arr = nrows - int((row['y'] - yllcorner) // cellsize) - 1

			dxl = row['x'] - (xindex * cellsize + xllcorner)
			dxu = cellsize - dxl
			dyl = row['y'] - (yindex_dem * cellsize + yllcorner)
			dyu = cellsize - dyl

			wt_00 = dxu * dyu / cellareafix
			wt_p0 = dxl * dyu / cellareafix
			wt_0p = dxu * dyl / cellareafix
			wt_pp = dxl * dyl / cellareafix

			height_sum[yindex_arr, xindex] += row['z'] * wt_00
			height_sum[yindex_arr, xindex + 1] += row['z'] * wt_p0
			height_sum[yindex_arr + 1, xindex] += row['z'] * wt_0p
			height_sum[yindex_arr + 1, xindex + 1] += row['z'] * wt_pp

			wt_sum[yindex_arr, xindex] += wt_00
			wt_sum[yindex_arr, xindex + 1] += wt_p0
			wt_sum[yindex_arr + 1, xindex] += wt_0p
			wt_sum[yindex_arr + 1, xindex + 1] += wt_pp

		for i in range(nrows + 1):
			for j in range(ncols + 1):
				if height_sum[i, j] == 0 or wt_sum[i, j] == 0:
					result_smooth[i, j] = np.nan
				else:
					result_smooth[i, j] = height_sum[i, j] / wt_sum[i, j]

		result_smoothed_wl = spn.convolve(result_smooth, smoothing_kernel)

		if self.smoothafterdamp.get():
			res1 = self.damping(self.dampingval.get(), result_smooth, result_fix)
			result = self.smoothing(self.smoothingval.get(), res1)
		else:
			res1 = self.smoothing(self.smoothingval.get(), result_smooth)
			result = self.damping(self.dampingval.get(), res1, result_fix)

		# # Smoothing
		# result_smoothed = self.smoothing(self.smoothingval.get(), result_smooth)
		#
		# # Damping
		# result = self.damping(self.dampingval.get(), result_smoothed, result_fix)

		self.result = result.copy()
		self.result[np.isnan(self.result)] = self.nodatavalue
		self.ncols = ncols
		self.nrows = nrows
		self.xllcorner = xllcorner
		self.yllcorner = yllcorner
		self.cellsize_forsave = cellsize

		fig1 = plt.figure()
		ax1 = fig1.add_subplot(2,3,1)
		ax1.imshow(result_fix)
		ax2 = fig1.add_subplot(2,3,2)
		im = ax2.imshow(result)
		ax3 = fig1.add_subplot(2,3,3)
		ax3.imshow(result_smooth)
		ax4 = fig1.add_subplot(2,3,4)
		ax4.imshow(result_smoothed_wl)

		fig1.subplots_adjust(right=0.8)
		# put colorbar at desire position
		cbar_ax = fig1.add_axes([0.85, 0.15, 0.05, 0.7])
		fig1.colorbar(im, cax=cbar_ax)

		plt.show()

		# Draw
		# Construct grid
		nx = cellsize * ncols
		ny = cellsize * nrows

		xend = xllcorner + nx
		yend = yllcorner + ny

		x = np.linspace(xllcorner, xend, ncols + 1)
		y = np.linspace(yend, yllcorner, nrows + 1)

		X, Y = np.meshgrid(x, y)

		# ax.plot_surface(X, Y, self.actDEM, rstride=2, cstride=2)
		self.parent.centreframe.fig.axes[0].plot_wireframe(X, Y, result, color=next(self.coloriter))

		self.parent.centreframe.drawArea.draw()

	def savedem(self):
		savefile = filedialog.asksaveasfilename()

		if savefile:
			header = 'ncols {}\n' \
					 'nrows {}\n' \
					 'xllcorner {}\n' \
					 'yllcorner {}\n' \
					 'cellsize {}\n' \
					 'NODATA_value {}'.format(self.ncols + 1, self.nrows + 1, self.xllcorner, self.yllcorner,
											  self.cellsize_forsave, self.nodatavalue)

			# header = 'NCOLS {}\n' \
			# 		 'NROWS {}\n' \
			# 		 'XLLCORNER {}\n' \
			# 		 'YLLCORNER {}\n' \
			# 		 'CELLSIZE {}\n' \
			# 		 'NODATA_VALUE {}'.format(ncols + 1, nrows + 1, xllcorner, yllcorner, cellsize, -9999)

			np.savetxt(savefile + '.txt', self.result, delimiter='\t', header=header, comments='')

	def clearall(self):
		self.interpoldict = {}
		self.fixdict = {}
		self.interpollist.delete_all()
		self.fixedlist.delete_all()

	def reset_colors(self):
		self.coloriter = iter(self.colors)


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

		ax = self.fig.add_subplot(1, 1, 1, projection='3d')
		self.drawArea.draw()

class BottomWindow(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)
		self.__createlayout()

	def __createlayout(self):
		self.buttonLabel = tk.Label(self, text='Actions')
		self.loadDataButton = tk.Button(self, text='Load cones', command=self.parent.loadCones)
		self.loadDemButton = tk.Button(self, text='Load DEM', command=self.parent.loadDEM)
		self.clearButton = tk.Button(self, text='Clear all', command= self.parent.clearall)

		self.buttonLabel.grid(column=0, row=0, sticky='nw')
		self.loadDataButton.grid(column=0, row=1, sticky='w', padx=2, pady=2)
		self.loadDemButton.grid(column=1, row=1, sticky='w', padx=2, pady=2)
		self.clearButton.grid(column=3, row=1, sticky='w', padx=2, pady=2)


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

		self.intersections = OrderedDict()
		self.highlight = 0

		self.nDataPoints = 0

		self.highlight = False

	def initUI(self):
		# Top level container
		self.centreframe = CentreWindow(self)
		self.bottomframe = BottomWindow(self)
		self.rightframe = RightWindow(self)

		# Layout top level frames
		self.centreframe.grid(row=0, column=1, sticky='nsew')
		self.bottomframe.grid(row=1, column=1, sticky='we')
		self.rightframe.grid(row=0, column=2, rowspan=2, sticky='ns')

		self.rowconfigure(0, weight=1)
		self.columnconfigure(1, weight=1)

	def loadDEM(self):

		self.highlight = 0
		self.intersections = OrderedDict()
		# Load DEM txt file
		fileloc = filedialog.askopenfilename(initialdir='./DEMs/')

		if fileloc:
			self.dem_dict = {}
			self.centreframe.fig.clear()

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

			# # Draw detector
			# det_vec = np.array([float(self.actHeader.loc['E (m, CH1903)', 'value']),
			#                     float(self.actHeader.loc['N (m, CH1903)', 'value']),
			#                     float(self.actHeader.loc['Z (m, CH1903)', 'value'])])
			# self.det_pc = ax.scatter(det_vec[0], det_vec[1], det_vec[2], s=10, c='magenta', marker='D')

		self.centreframe.drawArea.draw()

	def loadCones(self):
		fileloc = filedialog.askopenfilename()

		if fileloc:
			cones = pd.read_csv(fileloc, encoding='utf-16', index_col=0)

			for i in range(len(cones)):
				xvec = [cones.iloc[i, 0], cones.iloc[i, 3]]
				yvec = [cones.iloc[i, 1], cones.iloc[i, 4]]
				zvec = [cones.iloc[i, 2], cones.iloc[i, 5]]
				self.centreframe.fig.axes[0].plot(xvec, yvec, zvec, color='black')

		self.centreframe.drawArea.draw()

	def clearall(self):
		self.centreframe.fig.axes[0].clear()
		self.centreframe.drawArea.draw()

		self.rightframe.reset_colors()


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


def main():
	root = tk.Tk()
	root.title("Model Viewer")
	screen_width = root.winfo_screenwidth()
	screen_height = root.winfo_screenheight()
	sizestring = str(int(0.5*screen_width)) + 'x' + str(int(0.5*screen_height)) + '+50+50'
	root.geometry(sizestring)

	window = MainWindow(master=root)
	window.pack(side="top", fill="both", expand=True)
	root.mainloop()


if __name__ == '__main__':
	main()

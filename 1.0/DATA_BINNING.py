import tkinter as tk
import tkinter.ttk as ttk
from tkinter import filedialog
from tkinter import messagebox
from tkinter import font
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from mpl_toolkits.axes_grid1 import make_axes_locatable
from io import StringIO
from myListboxes import MultiListbox


class LeftWindow(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)
		self.__createlayout()

	def __createlayout(self):
		self.listlabel = tk.Label(self, text='Bin list')
		self.binListBox = MultiListbox(self, (('Nr°', 5), ('# Tracks', 10), (chr(952), 10), (chr(966), 10)))

		self.listlabel.grid(column=0, row=0, sticky='w')
		self.binListBox.grid(column=0, row=1, sticky='nsew')

		self.rowconfigure(1, weight=1)


class VisualisationWindow(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)
		self.__createlayout()

	def __createlayout(self):
		self.fig = plt.Figure()  # (figsize=(10,8))

		self.drawlabel = tk.Label(self, text='Visualisation')
		self.drawArea = FigureCanvasTkAgg(self.fig, self)

		self.drawlabel.grid(column=0, row=0, sticky='w', padx=2, pady=2)
		self.drawArea.get_tk_widget().grid(column=0, row=1, sticky='nsew', columnspan=2)

		self.rowconfigure(1, weight=1)
		self.columnconfigure(0, weight=1)

	def draw_load(self, locdf, globdf, thetamax, detectorangles):
		self.fig.clear()

		self.fig.suptitle('Detector orientation: Inclination {}° | Azimuth {}°'.format(detectorangles[0], detectorangles[1]))

		ax1 = self.fig.add_subplot(1, 2, 1, projection='polar')
		ax1.set_theta_zero_location('N')
		ax1.set_theta_direction(-1)
		ax1.set_rlim([0, thetamax])
		ax1.set_title('Local track distribution (Detector POV)')
		ax1.scatter(locdf.values[:, 1], locdf.values[:, 0]*180/np.pi, edgecolors='black', s=0.5)

		ax2 = self.fig.add_subplot(1, 2, 2, projection='polar')
		ax2.set_theta_zero_location('N')
		ax2.set_theta_direction(-1)
		ax2.set_rlim([0, 100])
		ax2.set_title('Global track distribution (Top view)')
		ax2.scatter(globdf.values[:, 1], globdf.values[:, 0]*180/np.pi, edgecolors='black', s=0.5)

		self.drawArea.draw()

	def draw_binning(self,bindf, locdf, globdf, polhisto, fluxhisto, thetamax, thetagrid, phigrid):
		self.fig.clear()

		ax1 = self.fig.add_subplot(2, 3, 1, projection='polar')
		ax1.set_title('Track distribution (Detector POV)')
		ax1.set_theta_zero_location('N')
		ax1.set_theta_direction(-1)
		ax1.set_rlim([0, thetamax])
		ax1.scatter(locdf.values[:, 1], locdf.values[:, 0] * 180 / np.pi, edgecolors='black', s=0.5)

		ax4 = self.fig.add_subplot(2, 3, 4, projection='polar')
		ax4.set_title('Global track distribution')
		ax4.set_theta_zero_location('N')
		ax4.set_theta_direction(-1)
		ax4.set_rlim([0, 100])
		ax4.scatter(globdf.values[:, 1], globdf.values[:, 0] * 180 / np.pi, edgecolors='black', s=0.5)

		colmap = plt.cm.jet
		norm = mpl.colors.Normalize(vmin=np.min(polhisto), vmax=np.max(polhisto))

		ax2 = self.fig.add_subplot(2, 3, 2, polar=True)
		ax2.set_title('Binned tracks')
		ax2.set_theta_zero_location('N')
		ax2.set_theta_direction(-1)
		# im2 = ax2.hist2d(self.actData.values[:, 1], self.actData.values[:, 0]*180/np.pi, bins=(phi, theta), cmin=self.minVal,
		#            cmap=colmap)
		im2 = ax2.pcolormesh(phigrid, thetagrid, polhisto, cmap=colmap)
		plt.colorbar(im2, ax=ax2)

		ax3 = self.fig.add_subplot(2, 3, 3, polar=True)
		ax3.set_title('Binned fluxes')
		ax3.set_theta_zero_location('N')
		ax3.set_theta_direction(-1)

		minval = np.power(10, np.floor(np.log10(np.nanmin(fluxhisto))))
		maxval = np.power(10, np.ceil(np.log10(np.nanmax(fluxhisto))))

		im3 = ax3.pcolormesh(phigrid, thetagrid, fluxhisto, cmap=colmap,
		                     norm=mpl.colors.LogNorm(vmin=minval, vmax=maxval))
		plt.colorbar(im3, ax=ax3)

		ax5 = self.fig.add_subplot(2, 3, 5, projection='polar')
		ax5.set_title('Bin distribution')
		ax5.set_theta_zero_location('N')
		ax5.set_theta_direction(-1)
		ax5.set_rlim([0, 100])
		ax5.scatter(bindf.values[:, 2]/180*np.pi, bindf.values[:, 1], edgecolors='red', s=1, marker='s')

		# #ax2 = self.fig.add_subplot(2, 2, 2, projection='polar')
		# ax2 = self.visualisationframe.fig.add_subplot(2, 3, 2, polar=True)
		# ax2.set_title('Binned tracks')
		# ax2.set_theta_zero_location('N')
		# ax2.set_theta_direction(-1)
		# im2 = ax2.hist2d(self.actData.values[:, 1], self.actData.values[:, 0], bins=(phi, theta), cmin=self.minVal,
		#            cmap=colmap)
		# plt.colorbar(im2[3], ax=ax2)

		# ax3 = self.visualisationframe.fig.add_subplot(2, 3, 3, polar=True)
		# ax3.set_title('Binned fluxes')
		# ax3.set_theta_zero_location('N')
		# ax3.set_theta_direction(-1)
		#
		# minval = np.power(10, np.floor(np.log10(np.nanmin(fluxhisto))))
		# maxval = np.power(10, np.ceil(np.log10(np.nanmax(fluxhisto))))
		#
		# im3 = ax3.pcolormesh(phigrid, thetagrid, fluxhisto, cmap=colmap,
		#                norm=mpl.colors.LogNorm(vmin=minval, vmax=maxval))
		# plt.colorbar(im3, ax=ax3)
		#
		#
		# ax4 = self.visualisationframe.fig.add_subplot(2, 3, 4)
		# ax4.hist(binlist, bins='auto', facecolor='green', edgecolor='black')
		# ax4.set_ylabel('# of bins')
		# ax4.set_xlabel('# tracks/bin')
		#
		#
		# ax5 = self.visualisationframe.fig.add_subplot(2, 3, 5)
		# ax5.hist2d(self.actData.values[:, 1]*180/np.pi, self.actData.values[:, 0], bins=(phi*180/np.pi, theta),
		# 		   cmin=self.minVal)
		# ax5.invert_yaxis()
		# ax5.set_ylabel('Zenith angle')
		# ax5.set_xlabel('Azimuth')
		#
		#
		# print(self.actData)
		# datval = self.actData.values.astype(float)
		# datval[:, 0] = datval[:, 0]/180*np.pi
		# datval[:, 1] = datval[:, 1]
		# datval = datval.T
		#
		# print(datval)
		#
		#
		# datshape = np.shape(datval)
		# qubit = np.zeros(datshape, dtype=complex)
		# for i in range(datshape[1]):
		# 	qubit[0, i] = np.cos(datval[0, i]/2)
		# 	qubit[1, i] = np.exp(datval[1, i]*1j)*np.sin(datval[0, i]/2)
		#
		# print('qubit:\n', qubit)
		#
		# thetaDdeg = (float(self.actHeader.loc['inclination angle (deg)'].value))
		# phiDdeg = float(self.actHeader.loc['azimuth (deg)'].value)
		#
		# thetaD = np.deg2rad(thetaDdeg)
		# phiD = np.deg2rad(phiDdeg)
		#
		# print('thetaD:\n', thetaD)
		#
		# Ry180 = np.array([[0, -1], [1, 0]])
		# qubity = np.matmul(Ry180, qubit)
		#
		# flipped = np.zeros_like(datval)
		# for i in range(datshape[1]):
		# 	flipped[0, i] = np.pi - 2*np.arctan(np.abs(qubity[1, i]/np.abs(qubity[0, i])))
		# 	flipped[1, i] = np.angle(qubity[1, i]) - np.angle(qubity[0, i])
		# 	if flipped[1, i] < 0:
		# 		flipped[1, i] = flipped[1, i] + np.pi *2
		#
		# print('flipped:\n', flipped)
		#
		# qubitnew = np.zeros(datshape, dtype=complex)
		# for i in range(datshape[1]):
		# 	qubitnew[0, i] = np.cos(flipped[0, i]/2)
		# 	qubitnew[1, i] = np.exp(flipped[1, i]*1j)*np.sin(flipped[0, i]/2)
		#
		#
		# Ry = np.array([[np.cos(thetaD/2), -np.sin(thetaD/2)], [np.sin(thetaD/2), np.cos(thetaD/2)]])
		# Rz = np.array([[np.exp(-1j*phiD/2), 0], [0, np.exp(1j*phiD/2)]])
		#
		# qubityy = np.matmul(Ry, qubitnew)
		# qubitxyz = np.matmul(Rz, qubityy)
		#
		# print('theta: ', thetaDdeg)
		# print('phi: ', phiDdeg)
		# print('X-rotation:\n', Ry)
		# print('Z-rotation:\n', Rz)
		#
		# rotatedData = np.zeros_like(datval)
		# for i in range(datshape[1]):
		# 	rotatedData[0, i] = 2*np.arctan(np.abs(qubitxyz[1, i]/np.abs(qubitxyz[0, i])))/np.pi*180
		# 	rotatedData[1, i] = np.angle(qubitxyz[1, i]) - np.angle(qubitxyz[0, i])
		# 	if rotatedData[1, i] < 0:
		# 		rotatedData[1, i] = rotatedData[1, i] + np.pi *2
		#
		# rotatedData = rotatedData.T
		# print('final:\n', rotatedData)
		#
		#

		# ax6.scatter(flipped[1, :], flipped[0, :]*180/np.pi, edgecolors='black', s=0.5)

		self.drawArea.draw()

	def setTheta(self, maxtheta):
		axes = self.fig.get_axes()
		axes[0].set_rlim([0, maxtheta])

		self.drawArea.draw()

class ButtonWindow(tk.Frame):
	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent, borderwidth=1)

		self.binoptions = tk.StringVar()
		self.binoptions.set('Detector centric')

		self.themax = tk.DoubleVar()
		self.themax.set(90)

		self.__createlayout()

	def __createlayout(self):
		self.buttonlabel = tk.Label(self, text='Loading')
		self.loadDatabutton = tk.Button(self, text='Load data', command=self.parent.loadData)

		self.separator = ttk.Separator(self, orient=tk.VERTICAL)

		self.bintext = tk.Label(self, text='Binning')
		self.binoptionmenu = tk.OptionMenu(self, self.binoptions, 'Detector centric', 'Global')
		self.thetamaxtext = tk.Label(self, text='Set maximum Theta (degree):')
		self.thetamaxentry = tk.Entry(self, textvariable=self.themax)
		self.setthetamaxbtn = tk.Button(self, text='Set', command=self.parent.setTheta)
		self.binButton = tk.Button(self, text='Bin data', command=self.parent.binData)

		self.separator2 = ttk.Separator(self, orient=tk.VERTICAL)

		self.savelabel = tk.Label(self, text='Saving')
		self.saveFigButton = tk.Button(self, text='Save figure', command=self.parent.saveFigure)
		self.saveDataButton = tk.Button(self, text='Save data', command=self.parent.saveData)


		self.buttonlabel.grid(column=0, row=0, sticky='w')
		self.loadDatabutton.grid(column=0, row=1, sticky='w', padx=2, pady=2)

		self.separator.grid(column=1, row=0, sticky='ns', padx=20, pady=2, rowspan=3)

		self.bintext.grid(column=3, row=0, padx=2, pady=2)
		self.binoptionmenu.grid(column=3, row=1, padx=2, pady=2)
		self.thetamaxtext.grid(column=3, row=2, padx=2, pady=2)
		self.thetamaxentry.grid(column=4, row=2, padx=2, pady=2)
		self.setthetamaxbtn.grid(column=5, row=2, padx=10, pady=2)
		self.binButton.grid(column=6, row=1, padx=2, pady=2, rowspan=3)

		self.separator2.grid(column=7, row=0, sticky='ns', padx=20, pady=2, rowspan=3)

		self.savelabel.grid(column=8, row=0, sticky='w', padx=2, pady=2)
		self.saveFigButton.grid(column=8, row=1, sticky='w', padx=2, pady=2)
		self.saveDataButton.grid(column=8, row=2, sticky='w', padx=2, pady=2)


class MainWindow(tk.Frame):

	def __init__(self, master=None):
		self.parent = master
		tk.Frame.__init__(self, self.parent)
		self.initUI()

		# Initial Parameters
		self.actHeader = pd.DataFrame()
		self.actData = pd.DataFrame()
		self.saveHeader = pd.DataFrame()
		self.saveData = pd.DataFrame()
		self.globData = pd.DataFrame()

		self.thetabins = 0
		self.phibins = 0
		self.minVal = 1
		self.ID = ''

		self.actfile = None

		self.maxtheta = 90

	def initUI(self):
		self.leftframe = LeftWindow(self)
		self.visualisationframe = VisualisationWindow(self)
		self.buttonframe = ButtonWindow(self)

		self.leftframe.grid(row=0, column=0, rowspan=2, sticky='ns')
		self.visualisationframe.grid(row=0, column=1, sticky='nsew')
		self.buttonframe.grid(row=1, column=1, sticky='we')

		self.rowconfigure(0, weight=1)
		self.columnconfigure(1, weight=1)

		self.buttonframe.thetamaxentry.config(state=tk.DISABLED)
		self.buttonframe.setthetamaxbtn.config(state=tk.DISABLED)

	def setTheta(self, *args):
		self.maxtheta = self.buttonframe.themax.get()

		if self.maxtheta > 90:
			self.maxtheta = 90
			self.buttonframe.themax.set(90)
		elif self.maxtheta < 0:
			self.maxtheta = 0
			self.buttonframe.themax.set(0)

		self.visualisationframe.setTheta(self.maxtheta)

	def reload(self, *args):
		if self.actfile is not None:
			self.loadData(self.actfile)

	def loadData(self, file=None):
		if file is None:
			fileloc = filedialog.askopenfilename(initialdir='./trackData', title='Select detector file')
			self.actfile = fileloc
		else:
			fileloc = file

		if fileloc:
			self.buttonframe.thetamaxentry.config(state=tk.NORMAL)
			self.buttonframe.setthetamaxbtn.config(state=tk.ACTIVE)

			self.phibins = 0
			self.thetabins = 0

			self.leftframe.binListBox.delete_all()


			with open(fileloc, 'r') as file:
				header_parse = pd.read_csv(StringIO(file.read()), sep='#', skipinitialspace=True, header=None)

			hd = header_parse.iloc[:, 1].dropna()
			self.actHeader = hd.str.extract('\s*(?P<key>[^;]+)\s*;\s*(?P<value>.+)', expand=True).dropna()
			self.actHeader.set_index('key', inplace=True)

			with open(fileloc, 'r') as file:
				self.actData = pd.read_csv(StringIO(file.read()), sep=';', skiprows=7)

			self.actData['z'] = 1

			self.visualisationframe.drawlabel.configure(text='Visualisation: (#Tracks: %s)' % self.actData.shape[0])

			# Calculation of detector rest frame data
			new_data = pd.DataFrame(columns=['theta', 'phi'])

			MAX = self.actData.shape[0]

			progVar = tk.DoubleVar()
			progBar = ttk.Progressbar(self.visualisationframe, length=100, variable=progVar)
			progBar.grid(column=0, row=0, sticky='e')
			self.update()

			for i in range(0, MAX):
				x = self.actData.iloc[i, 0]
				y = self.actData.iloc[i, 1]
				z = self.actData.iloc[i, 2]

				new_data.loc[i, 'theta'] = np.arccos(z / np.sqrt(x ** 2 + y ** 2 + z ** 2))
				# new_data.loc[i, 'theta'] = np.arccos(z / np.sqrt(x ** 2 + y ** 2 + z ** 2))
				new_data.loc[i, 'phi'] = np.pi / 2 - np.arctan2(y, x)
				if new_data.loc[i, 'phi'] < 0:
					new_data.loc[i, 'phi'] = new_data.loc[i, 'phi'] + 2 * np.pi

				progVar.set(float(i) / float(MAX) * 100)
				self.update_idletasks()

			self.actData = new_data

			# Calculation of globalt rest frame data (rotation on two spherical coordinate systems <-> rotations on bloch sphere)
			# Get data as numpy array
			datval = self.actData.values.astype(float)
			datval = datval.T

			datshape = np.shape(datval)

			# Transform data into a qubit
			qubit = np.zeros(datshape, dtype=complex)
			for i in range(datshape[1]):
				qubit[0, i] = np.cos(datval[0, i] / 2)
				qubit[1, i] = np.exp(datval[1, i] * 1j) * np.sin(datval[0, i] / 2)

			# Get detector rotation angles
			thetaDdeg = (float(self.actHeader.loc['inclination angle (deg)'].value))
			phiDdeg = float(self.actHeader.loc['azimuth (deg)'].value)

			thetaD = np.deg2rad(thetaDdeg)
			phiD = np.deg2rad(phiDdeg)

			# Perform first a flip of the data (behind the detector to front of the detector
			# This can be achieved by rotation the spherical data 180° in x direction (which is now on the bottom half
			# of the sphere, and then projecting it again to the upper half

			# Rotate
			Rx180 = np.array([[0, -1], [1, 0]])
			qubity = np.matmul(Rx180, qubit)

			# Backconversion from qubit to angles & projection to upper half
			flipped = np.zeros_like(datval)
			for i in range(datshape[1]):
				flipped[0, i] = np.pi - 2 * np.arctan(np.abs(qubity[1, i] / np.abs(qubity[0, i])))
				flipped[1, i] = np.angle(qubity[1, i]) - np.angle(qubity[0, i])
				if flipped[1, i] < 0:
					flipped[1, i] = flipped[1, i] + np.pi * 2

			# Convert back to qubit for detector rotations
			qubitnew = np.zeros(datshape, dtype=complex)
			for i in range(datshape[1]):
				qubitnew[0, i] = np.cos(flipped[0, i] / 2)
				qubitnew[1, i] = np.exp(flipped[1, i] * 1j) * np.sin(flipped[0, i] / 2)

			# Rotation
			Ry = np.array([[np.cos(thetaD / 2), -np.sin(thetaD / 2)], [np.sin(thetaD / 2), np.cos(thetaD / 2)]])
			Rz = np.array([[np.exp(-1j * phiD / 2), 0], [0, np.exp(1j * phiD / 2)]])

			qubityy = np.matmul(Ry, qubitnew)
			qubitxyz = np.matmul(Rz, qubityy)

			# Backconversion to angles
			rotatedData = np.zeros_like(datval)
			for i in range(datshape[1]):
				rotatedData[0, i] = 2 * np.arctan(np.abs(qubitxyz[1, i] / np.abs(qubitxyz[0, i])))
				rotatedData[1, i] = np.angle(qubitxyz[1, i]) - np.angle(qubitxyz[0, i])
				if rotatedData[1, i] < 0:
					rotatedData[1, i] = rotatedData[1, i] + np.pi * 2

			rotatedData = rotatedData.T
			self.globData = pd.DataFrame(rotatedData, columns=('theta', 'phi'))

			thetavar = float(self.buttonframe.thetamaxentry.get())
			if thetavar > 90:
				thetavar = 90
			elif thetavar < 0:
				thetavar = 0

			# Draw on screen
			self.visualisationframe.draw_load(self.actData, self.globData, thetavar, (thetaDdeg, phiDdeg))

			progBar.destroy()

	def binData(self):
		if self.actData.empty:
			messagebox.showerror("Error", "No data loaded. Load data first!")
		else:
			def acceptEntries(self):
				tbins = tentry.get()
				pbins = pentry.get()
				miva = minentry.get()

				if tbins.isdigit() and pbins.isdigit():

					tbins = int(tbins)
					pbins = int(pbins)

					if tbins > 0 and pbins > 0:
						self.thetabins = tbins
						self.phibins = pbins
						self.ID = identry.get()

						if miva.isdigit():
							self.minVal = int(miva)

						self.buttonframe.thetamaxentry.config(state=tk.DISABLED)
						self.buttonframe.setthetamaxbtn.config(state=tk.DISABLED)
						binWindow.destroy()
					else:
						messagebox.showerror("Error", "# Bins have to be greater than 0!")
				else:
					messagebox.showerror("Error", "# Bins must be a positive integer!")

			def clearEntries():
				tentry.delete(0, 'end')
				pentry.delete(0, 'end')

			binWindow = tk.Toplevel(self)
			binWindow.title('Bin selection')

			tlbl = tk.Label(binWindow, text='# bins Zenith: ')
			plbl = tk.Label(binWindow, text='# bins Azimuth: ')
			minlbl = tk.Label(binWindow, text='Min # of tracks: ')
			idlbl = tk.Label(binWindow, text='Detector ID: ')

			tentry = tk.Entry(binWindow)
			pentry = tk.Entry(binWindow)
			minentry = tk.Entry(binWindow)
			identry = tk.Entry(binWindow)

			accButton = tk.Button(binWindow, text='Accept', command= lambda: acceptEntries(self))
			clearButton = tk.Button(binWindow, text='Clear', command=clearEntries)

			tlbl.grid(column=0, row=0, sticky='w', padx=5, pady=5)
			plbl.grid(column=0, row=1, sticky='w', padx=5, pady=5)
			minlbl.grid(column=0, row=2, sticky='w', padx=5, pady=5)
			idlbl.grid(column=0, row=3, sticky='w', padx=5, pady=5)

			tentry.grid(column=1, row=0, sticky='e', padx=5, pady=5)
			pentry.grid(column=1, row=1, sticky='e', padx=5, pady=5)
			minentry.grid(column=1, row=2, sticky='e', padx=5, pady=5)
			identry.grid(column=1, row= 3, sticky='e', padx=5, pady=5)

			accButton.grid(column=0, row=4, sticky='w', padx=5, pady=5)
			clearButton.grid(column=1, row=4, sticky='e', padx=5, pady=5)

			x = int((self.winfo_screenwidth() - self.winfo_reqwidth()) / 2)
			y = int((self.winfo_screenheight() - self.winfo_reqheight()) / 2)
			binWindow.geometry('+{}+{}'.format(x, y))

			binWindow.transient(self)
			binWindow.grab_set()
			self.wait_window(binWindow)

			if self.thetabins > 0 and self.phibins > 0:
				self.buttonframe.thetamaxentry.config(state=tk.DISABLED)
				self.buttonframe.setthetamaxbtn.config(state=tk.DISABLED)

				self.maxtheta = float(self.buttonframe.thetamaxentry.get())

				self.saveData = pd.DataFrame()
				self.saveHeader = pd.DataFrame()

				binopt = self.buttonframe.binoptions.get()

				costhemax = np.cos(float(self.buttonframe.thetamaxentry.get())/180*np.pi)

				phi = np.linspace(0, 2 * np.pi, self.phibins + 1)
				phideg = phi * 180/np.pi
				if binopt == 'Global':
					theta = np.arccos(np.linspace(1, 0, self.thetabins + 1))/np.pi*180
				else:
					theta = np.arccos(np.linspace(1, costhemax, self.thetabins + 1)) / np.pi * 180

				solAng = (phi[1]-phi[0])*(np.cos(np.deg2rad(theta[0]))-np.cos(np.deg2rad(theta[1])))

				if binopt == 'Global':
					polhisto, _, _ = np.histogram2d(self.globData.values[:, 1], self.globData.values[:, 0]*180/np.pi, bins=(phi, theta))
				else:
					polhisto, _, _ = np.histogram2d(self.actData.values[:, 1],
					                                self.actData.values[:, 0] * 180 / np.pi, bins=(phi, theta))

				polhisto = polhisto.T
				fluxhisto = np.zeros_like(polhisto)

				# print(solAng)
				# print(polhisto)
				# print(fluxhisto)

				phigrid, thetagrid = np.meshgrid(phi, theta)

				time = float(self.actHeader.loc['exposure time (sec)'].value)
				area = float(self.actHeader.loc['effective area (cm2)'].value)

				binlist = []
				self.leftframe.binListBox.delete_all()
				self.saveData = pd.DataFrame(columns=['# Tracks', chr(952), chr(966)])
				for i in range(0, np.shape(polhisto)[1]):
					for j in range(0, np.shape(polhisto)[0]):
						if self.minVal <= polhisto[j, i]:
							binlist.append(polhisto[j, i])
							phimean = (phideg[i+1]+phideg[i])/2
							thetamean = (theta[j+1]+theta[j])/2

							if binopt == 'Detector centric':
								thetamean = np.deg2rad(thetamean)
								phimean = np.deg2rad(phimean)

								# Transform data into a qubit
								qubit = np.zeros((2,1), dtype=complex)

								qubit[0] = np.cos(thetamean / 2)
								qubit[1] = np.exp(phimean * 1j) * np.sin(thetamean / 2)

								# Get detector rotation angles
								thetaDdeg = (float(self.actHeader.loc['inclination angle (deg)'].value))
								phiDdeg = float(self.actHeader.loc['azimuth (deg)'].value)

								thetaD = np.deg2rad(thetaDdeg)
								phiD = np.deg2rad(phiDdeg)

								# Perform first a flip of the data (behind the detector to front of the detector
								# This can be achieved by rotation the spherical data 180° in x direction (which is now on the bottom half
								# of the sphere, and then projecting it again to the upper half

								# Rotate
								Rx180 = np.array([[0, -1], [1, 0]])
								qubity = np.matmul(Rx180, qubit)

								# Backconversion from qubit to angles & projection to upper half
								flipped = np.zeros_like(qubit)
								flipped[0] = np.pi - 2 * np.arctan(
									np.abs(qubity[1] / np.abs(qubity[0])))
								flipped[1] = np.angle(qubity[1]) - np.angle(qubity[0])
								if flipped[1] < 0:
									flipped[1] += np.pi * 2

								# Convert back to qubit for detector rotations
								qubitnew = np.zeros_like(qubit)
								qubitnew[0] = np.cos(flipped[0] / 2)
								qubitnew[1] = np.exp(flipped[1] * 1j) * np.sin(flipped[0] / 2)

								# Rotation
								Ry = np.array([[np.cos(thetaD / 2), -np.sin(thetaD / 2)],
								               [np.sin(thetaD / 2), np.cos(thetaD / 2)]])
								Rz = np.array([[np.exp(-1j * phiD / 2), 0], [0, np.exp(1j * phiD / 2)]])

								qubityy = np.matmul(Ry, qubitnew)
								qubitxyz = np.matmul(Rz, qubityy)

								# Backconversion to angles
								rotatedData = np.zeros_like(qubit, dtype=float)
								rotatedData[0] = 2 * np.arctan(
									np.abs(qubitxyz[1] / np.abs(qubitxyz[0])))
								rotatedData[1] = np.angle(qubitxyz[1]) - np.angle(qubitxyz[0])
								if rotatedData[1] < 0:
									rotatedData[1] = rotatedData[1] + np.pi * 2

								thetameanrot = rotatedData[0,0] *180/np.pi
								phimeanrot = rotatedData[1,0] *180/np.pi

								self.saveData.loc[len(self.saveData)] = [polhisto[j, i], thetameanrot, phimeanrot]
								self.leftframe.binListBox.insert(tk.END, (len(binlist), int(polhisto[j, i]),
								                                          f'{thetameanrot:9.2f}', f'{phimeanrot:9.2f}'))
							else:
								self.saveData.loc[len(self.saveData)] = [polhisto[j, i], thetamean, phimean]
								self.leftframe.binListBox.insert(tk.END, (len(binlist), int(polhisto[j, i]),
								                                          f'{thetamean:9.2f}', f'{phimean:9.2f}'))


							fluxhisto[j, i] = polhisto[j, i] / (solAng * time * area * np.cos(thetamean/180*np.pi))
						else:
							polhisto[j, i] = np.nan
							fluxhisto[j, i] = np.nan

				print(self.saveData)

				self.visualisationframe.draw_binning(self.saveData, self.actData, self.globData, polhisto, fluxhisto, self.maxtheta, thetagrid, phigrid)

				self.visualisationframe.drawlabel.configure(text='Visualisation: (#Tracks: %s | #Bins: %s)'
			                              % (self.actData.shape[0], len(binlist)))

				self.saveHeader = self.actHeader.copy()
				self.saveHeader.reset_index(inplace=True)

				self.saveHeader = self.saveHeader.append({'key': 'ID', 'value': self.ID}, ignore_index=True)
				self.saveHeader = self.saveHeader.append({'key': 'Solid angle (sr)', 'value': solAng}, ignore_index=True)

				self.saveHeader.set_index('key', inplace=True)


	def saveFigure(self):
		if not self.visualisationframe.fig.get_axes():
			messagebox.showerror('Error', 'Figure is empty!')
		else:
			fileloc = filedialog.asksaveasfilename(initialdir='./images/', title='Save figure',
												   filetypes=[('Portable network graphic','png'),('SVG-file','svg'),
															   ('Portable file document','pdf')], defaultextension='.png')
			if fileloc:
				self.visualisationframe.fig.savefig(fileloc)

	def saveData(self):
		if self.leftframe.binListBox.index(tk.END) == 0:
			messagebox.showerror('Error', 'There are currently no bins!')
		else:
			fileloc = filedialog.asksaveasfilename(initialdir='./binnedData/', title='Save data')

			if fileloc:

				if fileloc.endswith('_aux.json'):
					fileloc = fileloc.strip('_aux.json')
				elif fileloc.endswith('.json'):
					fileloc = fileloc.strip('.json')

				self.saveHeader.to_json(fileloc + '_aux.json')
				self.saveData.to_json(fileloc + '.json', force_ascii=True)


def main():
	root = tk.Tk()
	root.title('Data binning app')

	screen_width = root.winfo_screenwidth()
	screen_height = root.winfo_screenheight()
	sizestring = str(int(0.8 * screen_width)) + 'x' + str(int(0.8 * screen_height)) + '+50+50'
	root.geometry(sizestring)

	window = MainWindow(master=root)
	window.pack(side="top", fill="both", expand=True)
	root.mainloop()


if __name__ == '__main__':
	main()

import os
import os.path as op
import pandas as pd
import json
import numpy as np
import scipy as sp
import scipy.stats as sps
import functools
import math
import Material

from Parameters import *
from cosmicflux import getIntegratedFlux


class MultiUniformDistribution:
	def __init__(self):
		pass

	def logpdf(self, x, dimension, min, max):
		if any(xi < min or xi > max for xi in x):
			return -np.inf
		else:
			return -dimension * np.log(max - min)


class Cone:
	def __init__(self, parent, name, data):
		# Initialise parameters
		self.name = name
		self.parent = parent
		self.data = data
		self.parameter_names = []
		self.hstep = 20  # [cm]
		self.flux = 0
		self.logprob = 0
		self.fluerr = 0.15  # 15% error on flux model

		self.thickvar_flag = False

		self.ntracks = self.data['# Tracks']

		self.segments = []

		self.parameters = LocalParameterGroup()
		self.parent.parameters.declare_child(self.parameters)

		mats = self.data['Materials'].split()

		time = float(self.parent.detectors[self.data['Detector']].at['exposure time (sec)', 'value'])
		solidangle = float(self.parent.detectors[self.data['Detector']].at['Solid angle (sr)', 'value'])
		# TODO: calculate total area exposed to directional flux
		area = float(self.parent.detectors[self.data['Detector']].at['effective area (cm2)', 'value'])

		self.exposure = time * solidangle * area

		for mat in mats:
			self.segments.append([0, self.parent.get_material(mat)])

		for i in range(2, 1 + len(self.segments)):
			self.parameter_names.append(name + '_r' + str(i - 1))

		if len(self.segments) == 1:
			self.segments[0][0] = self.data['d_topo']

			self.lp_thicknesses = self.__ldelta_function
			self.lenpar = []
		else:
			self.thickvar_flag = True

			# Initialise thickness parameter
			thickratiopar = MultiDimParameter(name + '_r', 'thickness_ratio')
			thickratiopar.set_bounds(-2, 2)
			thickratiopar.set_value(sps.uniform.rvs(loc=thickratiopar.bounds[0],
													scale=thickratiopar.bounds[1] - thickratiopar.bounds[0],
													size=len(self.segments) - 1))
			thickratiopar.set_pdf(MultiUniformDistribution(), dimension=len(self.segments) - 1,
								  min=thickratiopar.bounds[0], max=thickratiopar.bounds[1])

			# thickratiopar = Parameter(name + '_dm', 'detector_side_material')
			# thickratiopar.set_bounds(0, self.data['d_topo'])
			# thickratiopar.set_value(sps.uniform.rvs(loc=thickratiopar.bounds[0],
			#                                         scale=thickratiopar.bounds[1] - thickratiopar.bounds[0],
			#                                         size=len(self.segments) - 1))
			# thickratiopar.set_pdf(MultiUniformDistribution(), dimension=len(self.segments) - 1,
			#                       min=thickratiopar.bounds[0], max=thickratiopar.bounds[1])

			self.parameters.add_Parameter(thickratiopar)
			self.parent.parameters.add_Parameter(thickratiopar)

			self.update_lengths(thickratiopar.value)

	def __ldelta_function(self, x):
		return 0

	def __lmuniform_function(self, x, mi, ma):
		lprob = 0
		for i in x:
			lprob += sps.uniform.logpdf(i, mi, ma)
		return lprob

	def update_lengths(self, ratios):
		explen = np.power(10, ratios)
		factor = self.data['d_topo'] / (math.fsum(explen) + 1)

		lengths = [xi * factor for xi in explen]
		lengths.append(factor)

		for il, litem in enumerate(lengths):
			self.segments[il][0] = litem

	# def update_lengths(self, distance):
	#     dist = distance[0]
	#     self.segments[0][0] = dist
	#     self.segments[1][0] = self.data['d_topo'] - dist

	def get_parameters(self):
		return self.parameters.iloc[[-1]]

	def get_logprob(self):
		return self.logprob

	def set_logprob(self, lprob):
		self.logprob = lprob

	def update(self):
		if self.thickvar_flag:
			self.update_lengths(self.parameters.parameters[self.name + '_r'].value)

	def show_segments(self):
		for seg in self.segments:
			print(seg)

	def get_name(self):
		return self.name

	def calculate_logprob(self):
		lprob = 0

		# Parameters
		parlprob = self.parameters.calculate_logprob()
		lprob += parlprob

		# Likelihood ("Convolution" of Lognormal with Poisson, Integration)
		llike = np.log(self.calculate_loglikelihood())
		lprob += llike

		# print('Cone: ', self.name, ' | ', parlprob, ' | ', llike)

		return lprob

	def calculate_loglikelihood(self):
		varflux = np.log(1 + self.fluerr ** 2)
		stdflux = varflux ** 0.5
		muflux = np.log(self.flux) - 0.5 * varflux

		(u, h) = np.linspace(-3, -1, 1000, retstep=True)

		integral = 0

		# Marginalise flux parameter
		for ui in u:
			intpoint = np.exp(ui - np.exp(-ui))
			# intpoint = ui

			poi = sps.poisson.pmf(self.ntracks, intpoint * self.exposure)
			logn = sps.lognorm.pdf(intpoint, s=stdflux, scale=self.flux)

			weight = (1 + np.exp(-ui)) * (np.exp(ui - np.exp(-ui)))
			# weight = 1

			integral += poi * logn * weight * h
			#print(intpoint, logn, stdflux, self.flux)

		return integral


	def calculate_cutoff_energy(self, cs_err_df):
		# Employ runge kutta sceme
		# Starting values
		ei = 1000  # [MeV]
		xi = 0
		# print(self.data['d_topo'])
		# print(self.segments)
		for seg in self.segments:
			# print(seg[1].name)
			# print(seg[0])
			xmax = seg[0] * 100
			n_rk = int(round(xmax / self.hstep))  # Approximate number of segments

			for i in range(0, n_rk):
				k1 = seg[1].get_energy_loss(ei, cs_err_df)
				k2 = seg[1].get_energy_loss(ei + self.hstep / 2. * k1, cs_err_df)
				k3 = seg[1].get_energy_loss(ei + self.hstep / 2. * k2, cs_err_df)
				k4 = seg[1].get_energy_loss(ei + self.hstep * k3, cs_err_df)

				ei = ei + self.hstep / 6. * (k1 + 2. * k2 + 2. * k3 + k4)
				xi = xi + self.hstep

		return ei

	def calculate_flux(self, cs_err_df):
		h0 = float(self.parent.detectors[self.data['Detector']].at['Z (m, CH1903)', 'value'])
		dist = float(self.data['d_topo'])
		angle = np.deg2rad(self.data[chr(0x03B8)])

		height = h0 + dist * np.sin(angle)

		self.flux = getIntegratedFlux(self.calculate_cutoff_energy(cs_err_df), height, angle)
		# print('Flux: ', self.flux)


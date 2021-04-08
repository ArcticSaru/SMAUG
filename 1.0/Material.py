import os
import numpy as np
import scipy.stats as sps
import pandas as pd
import csv
import math
import sys
import tkinter as tk
from tkinter import ttk
from functools import partial

from numpy.core._multiarray_umath import ndarray

from Parameters import *

# TODO: Probability mode for lognormal and from data


class Material:
	"""description of class"""

	# Class variables
	# Physical Constants
	K = 0.307075  # [MeV*cm^2*mol/g^2]
	massMuon = 105.6583568  # [MeV/c^2]
	massProton = 938.2720813  # [MeV/c^2]
	massElectron = 0.510998902  # [MeV/c^2]
	radiusElectron = 2.817940385E-13  # [cm]
	alpha = 1. / 137.03599976  # Fine structure constant []
	euler = 2.7182818284
	pi = 3.1415926535
	avogadro = 6.0221412927E23  # [mol^-1]
	hbarc = 0.1973269788  # h-bar * speed of light [GeV*fm]

	# Mathematical Constants
	numberIntegrationSegments = 63
	rangeNumberIntegrationSegments = 40
	numberRadlossTableEnergies = 100
	rungeKuttaStep = 20  # [cm]

	ionError = 1
	bremsError = 1
	pairError = 1
	photoError = 1

	def __init__(self, name):
		self.name = name

	def get_cutoff_energy(self, distance):

		h = self.rungeKuttaStep

		nsteps = int(np.around(distance / h))

		ei = 10.
		xi = 0.

		for i in range(1, nsteps + 1):
			k1 = self.get_energy_loss(ei)
			k2 = self.get_energy_loss(ei + h / 2. * k1)
			k3 = self.get_energy_loss(ei + h / 2. * k2)
			k4 = self.get_energy_loss(ei + h * k3)

			ei = ei + h / 6. * (k1 + 2. * k2 + 2. * k3 + k4)
			xi = xi + h

		return ei

	def get_energy_loss(self, energy, cs_err_df):

		engyloss = self.density * (cs_err_df[chr(0x03C3) + '_ion'] * self.ionisation_loss(energy) +
								   self.radiation_loss(energy, cs_err_df))

		return engyloss

	def create_radloss_table(self, spath):
		# Set up energies from 10 MeV - 100 TeV
		energies = np.logspace(1.0, 8.0, num=100, base=10.0)

		rlTable = np.zeros_like(energies)

		for ide, engy in enumerate(energies):
			print(ide)
			rlTable[ide] = (self.radiation_loss_calc(engy))

		self.rlTable = np.vstack((energies, rlTable)).T
		savefile = os.path.normpath(os.path.join(spath, self.name + '_RLT.txt'))
		np.savetxt(savefile, self.rlTable)

	def ionisation_loss(self, energy):
		# Define kinetic variables
		gamma = (energy + self.massMuon) / self.massMuon
		gammaSquared = gamma * gamma
		betaSquared = ((energy * (2. * self.massMuon + energy))
					   / ((energy + self.massMuon) * (energy + self.massMuon)))
		beta = betaSquared ** 0.5
		betagamma = beta * gamma
		betagammaSquared = betagamma * betagamma
		xCompare = np.log10(betagamma)

		# Dfeine Ionisation Loss variables
		qMax = ((2. * self.massElectron * betagammaSquared)
				/ (1. + 2. * gamma * self.massElectron / self.massMuon
				   + (self.massElectron / self.massMuon) * (self.massElectron / self.massMuon)))

		correction = (self.K / (4. * self.pi) * self.Z / self.A * self.alpha
					  * (np.log(2. * (energy + self.massMuon) / self.massMuon)
						 - np.log(2 * qMax / self.massElectron) / 3) * np.power(np.log(2 * qMax / self.massElectron),
																				2))

		shellCorrection = ((0.422377 / betagammaSquared + 0.0304043 / (betagammaSquared * betagammaSquared)
							- 0.00038106 / (
										betagammaSquared * betagammaSquared * betagammaSquared)) * 1E-6 * self.I * self.I
						   + (3.858016 / betagammaSquared - 0.1667989 / (betagammaSquared * betagammaSquared)
							  + 0.00157955 / (
										  betagammaSquared * betagammaSquared * betagammaSquared)) * 1E-9 * self.I * self.I * self.I)

		if (xCompare < self.x0):
			delta = self.delta * np.power(10., 2. * (xCompare - self.x0))
		elif (xCompare > self.x1):
			delta = 2 * np.log(10) * xCompare - self.C
		else:
			delta = 2 * np.log(10) * xCompare - self.C + self.a * np.power(self.x1 - xCompare, self.k)

		# Return Ionisation dE/dX
		first = 0.5 * np.log((2. * self.massElectron * betagammaSquared * qMax * 1E12) / (self.I * self.I))
		second = - betaSquared
		third = - delta / 2.
		fourth = (qMax * qMax) / (8. * (gamma * gamma * self.massMuon * self.massMuon))
		fifth = - shellCorrection / self.Z
		sixth = correction

		result = self.K / betaSquared * self.Z / self.A * (first + second + third + fourth + fifth) + sixth

		return result

	def radiation_loss(self, energy, cs_err_df):
		# Interpolate on a logarithmic grid
		# Get correct energy index
		logx = np.log10(energy)

		if logx < 1 or logx > 8:
			print('Energy exceeds the range from 10 MeV to 100 TeV')
			exit(1)

		k = math.floor((logx - 1) * 99. / 7.)

		bre1 = self.rlTable[k, 1]
		bre2 = self.rlTable[k + 1, 1]
		pai1 = self.rlTable[k, 2]
		pai2 = self.rlTable[k + 1, 2]
		pho1 = self.rlTable[k, 3]
		pho2 = self.rlTable[k + 1, 3]

		lx1 = np.log10(self.rlTable[k, 0])
		lx2 = np.log10(self.rlTable[k + 1, 0])

		exponent = (logx - lx1) / (lx2 - lx1)

		bre = bre1 * np.power(bre2 / bre1, exponent) * cs_err_df[chr(0x03C3) + '_bre']
		pai = pai1 * np.power(pai2 / pai1, exponent) * cs_err_df[chr(0x03C3) + '_pai']
		pho = pho1 * np.power(pho2 / pho1, exponent) * cs_err_df[chr(0x03C3) + '_pho']

		return (bre + pai + pho) * energy

	def radiation_loss_calc(self, energy):
		result = (self.bremsstrahlung_loss_calc(energy) + self.pairproduction_loss_calc(
			energy) + self.photonuclear_loss_calc(energy)) * energy
		return result

	def bremsstrahlung_loss_calc(self, energy):  # After Groom

		# Prepare variables
		bremsstrahlungLoss = 0.
		bremsIntegrationLowerBound = 0.
		bremsIntegrationUpperBound = 1.
		bremsIntervalBegin = -3.
		bremsIntervalEnd = 3.
		bremsIntegrationSegmentLength = (bremsIntervalEnd - bremsIntervalBegin) / self.numberIntegrationSegments
		bremsDn = 1.54 * np.power(self.A, 0.27)

		# Check for Hydrogen
		if (self.Z == 1):
			bremsB = 202.4
			bremsB2 = 446.
		else:
			bremsB = 182.7
			bremsB2 = 1429.

		# Bremsstrahlung Integration
		for segmentIterator in range(0, self.numberIntegrationSegments + 1):
			if (bremsIntegrationUpperBound > 0):
				# exponential integration point and weights
				bremsDoubleExponentialIntegrationPoint = np.tanh(
					self.pi / 2. * np.sinh(bremsIntervalBegin + segmentIterator * bremsIntegrationSegmentLength))

				bremsFractionalEnergyTransfer = ((bremsIntegrationLowerBound + bremsIntegrationUpperBound) / 2.
												 + (
															 bremsIntegrationUpperBound - bremsIntegrationLowerBound) / 2. * bremsDoubleExponentialIntegrationPoint)

				bremsIntegrationWeight = ((bremsIntegrationUpperBound - bremsIntegrationLowerBound) / 2.
										  * (self.pi / 2. * np.cosh(
							bremsIntervalBegin + segmentIterator * bremsIntegrationSegmentLength))
										  / np.power(np.cosh(self.pi / 2. * np.sinh(
							bremsIntervalBegin + segmentIterator * bremsIntegrationSegmentLength)), 2.))

				# Parameters
				delta = self.massMuon * self.massMuon * bremsFractionalEnergyTransfer / (
							2. * energy * (1. - bremsFractionalEnergyTransfer))

				# Bremsstrahlung Electronic Crosssection
				phiin = (np.log(
					(self.massMuon / delta) / (self.massMuon * delta / (self.massElectron * self.massElectron)
											   + np.sqrt(self.euler))) - np.log(
					1. + self.massElectron / (delta * bremsB2 * np.sqrt(self.euler) * np.power(self.Z, -2. / 3.))))

				bremsCrosssectionElec = (
							4 * self.alpha * self.Z * np.power(self.massElectron * self.radiusElectron / self.massMuon,
															   2.)
							* (
										4. / 3. - 4. * bremsFractionalEnergyTransfer / 3. + bremsFractionalEnergyTransfer * bremsFractionalEnergyTransfer)
							* phiin * self.avogadro / self.A)

				# Bremsstrahlung Nuclear Crosssection
				DelN = np.log(bremsDn / (1 + delta * (bremsDn * np.sqrt(self.euler) - 2.) / self.massMuon))

				Phi = (np.log((bremsB * self.massMuon * np.power(self.Z, -1. / 3.) / self.massElectron)
							  / (1. + delta * np.sqrt(self.euler) * bremsB * np.power(self.Z,
																					  -1. / 3.) / self.massElectron)) - DelN)

				bremsCrosssectionNucl = (
							4. * self.alpha * np.power(self.Z * self.massElectron * self.radiusElectron / self.massMuon,
													   2.)
							* (
										4. / 3. - 4. * bremsFractionalEnergyTransfer / 3. + bremsFractionalEnergyTransfer * bremsFractionalEnergyTransfer)
							* Phi * self.avogadro / self.A)

				# std::cout << bremsCrosssectionElec << " | " << bremsCrosssectionNucl << std::endl

				# Calculate Bremsstrahlung crosssection
				bremsTotalCrosssection = bremsCrosssectionNucl + bremsCrosssectionElec

				# Integrate
				bremsstrahlungLoss += bremsIntegrationWeight * bremsIntegrationSegmentLength * bremsTotalCrosssection

		# std::cout << bremsstrahlungEnergyLossG4(energy) << " | "
		if (bremsstrahlungLoss < 0):
			bremsstrahlungLoss = 0.

		return bremsstrahlungLoss

	def pairproduction_loss_calc(self, energy):  # After GEANT 4

		# Prepare variables
		pairproductionLoss = 0.
		pairIntegrationLowerBound = 4. * self.massElectron / energy
		pairIntegrationUpperBound = 1. - (3. * np.sqrt(self.euler) * self.massMuon * np.power(self.Z, 1. / 3.)) / (
					4. * energy)
		pairIntervalBegin = -3.
		pairIntervalEnd = 3.
		pairIntegrationSegmentLength = (pairIntervalEnd - pairIntervalBegin) / self.numberIntegrationSegments
		pairTheta = np.power(self.massMuon / self.massElectron, 2.)

		# Check for Hydrogen
		if (self.Z == 1.):
			pairAstar = 202.4
			pairGamma1 = 4.4E-5
			pairGamma2 = 4.8E-5
		else:
			pairAstar = 183.
			pairGamma1 = 1.95E-5
			pairGamma2 = 5.3E-5

		if (pairIntegrationUpperBound > pairIntegrationLowerBound):

			# Pair production Integration
			for segmentIterator in range(0, self.numberIntegrationSegments + 1):
				if (pairIntegrationUpperBound > 0):

					# exponential integration point and weights
					pairDoubleExponentialIntegrationPoint = np.tanh(
						self.pi / 2. * np.sinh(pairIntervalBegin + segmentIterator * pairIntegrationSegmentLength))

					pairEnergyTransfer = ((pairIntegrationLowerBound + pairIntegrationUpperBound) / 2.
										  + (
													  pairIntegrationUpperBound - pairIntegrationLowerBound) / 2. * pairDoubleExponentialIntegrationPoint)

					pairIntegrationWeight = ((pairIntegrationUpperBound - pairIntegrationLowerBound) / 2.
											 * (self.pi / 2. * np.cosh(
								pairIntervalBegin + segmentIterator * pairIntegrationSegmentLength))
											 / np.power(np.cosh(self.pi / 2. * np.sinh(
								pairIntervalBegin + segmentIterator * pairIntegrationSegmentLength)), 2.))

					# Changing variables
					# For fractional crosssection
					pairBeta = pairEnergyTransfer * pairEnergyTransfer / (2. * (1. - pairEnergyTransfer))

					# For rho integration
					pairG = 0.
					pairRhoIntervalBegin = -3.
					pairRhoIntervalEnd = 3.
					pairRhoIntegrationSegmentLength = (
																  pairRhoIntervalEnd - pairRhoIntervalBegin) / self.numberIntegrationSegments
					pairRhoIntegrationLowerBound = 0.

					if (abs(pairIntegrationLowerBound - pairEnergyTransfer) <= 1E-16):
						pairEnergyTransfer = pairIntegrationLowerBound

					pairRhoIntegrationUpperBound = (
								(1. - 6. * np.power(self.massMuon, 2.) / (energy * energy * (1. - pairEnergyTransfer)))
								* np.sqrt(1. - 4. * self.massElectron / (energy * pairEnergyTransfer)))

					# std::cout << segmentIterator << " | " << pairIntegrationLowerBound << " | " << pairIntegrationUpperBound << " | " << pairEnergyTransfer << std::endl

					# Do rho integration
					for rhoIterator in range(0, self.numberIntegrationSegments + 1):

						pairRhoDEIntegrationPoint = np.tanh(self.pi / 2. * np.sinh(
							pairRhoIntervalBegin + rhoIterator * pairRhoIntegrationSegmentLength))

						pairRhoIntegrationWeight = ((pairRhoIntegrationUpperBound - pairRhoIntegrationLowerBound) / 2.
													* (self.pi / 2. * np.cosh(
									pairRhoIntervalBegin + rhoIterator * pairRhoIntegrationSegmentLength))
													/ np.power(np.cosh(self.pi / 2. * np.sinh(
									pairRhoIntervalBegin + rhoIterator * pairRhoIntegrationSegmentLength)), 2.))

						pairRhoRho = ((pairRhoIntegrationLowerBound + pairRhoIntegrationUpperBound) / 2.
									  + (
												  pairRhoIntegrationUpperBound - pairRhoIntegrationLowerBound) / 2. * pairRhoDEIntegrationPoint)

						pairRhoRhoSquared = pairRhoRho * pairRhoRho

						pairRhoChi = np.power(self.massMuon * pairEnergyTransfer / (2. * self.massElectron), 2.) * (
									1. - (pairRhoRhoSquared)) / (1. - pairEnergyTransfer)

						pairRhoYe = ((5. - pairRhoRhoSquared + 4. * pairBeta * (1. + pairRhoRhoSquared))
									 / (2. * (1. + 3. * pairBeta) * np.log(
									3. + 1. / pairRhoChi) - pairRhoRhoSquared - 2. * pairBeta * (
													2. - pairRhoRhoSquared)))

						pairRhoYu = ((4. + pairRhoRhoSquared + 3. * pairBeta * (1. + pairRhoRhoSquared))
									 / ((1. + pairRhoRhoSquared) * (1.5 + 2. * pairBeta) * np.log(
									3. + pairRhoChi) + 1. - 1.5 * pairRhoRhoSquared))

						pairRhoLe = (np.log(
							(pairAstar * np.power(self.Z, -1. / 3.) * np.sqrt((1. + pairRhoChi) * (1. + pairRhoYe)))
							/ (1. + (2. * self.massElectron * np.sqrt(self.euler) * pairAstar * np.power(self.Z,
																										 -1. / 3.)
									 * (1. + pairRhoChi) * (1 + pairRhoYe)) / (
										   energy * pairEnergyTransfer * (1. - pairRhoRhoSquared)))) - 1. / 2.
									 * np.log(
									1. + (9. * self.massElectron * self.massElectron * np.power(self.Z, 2. / 3.)
										  / (4. * self.massMuon * self.massMuon)) * (1. + pairRhoChi) * (
												1. + pairRhoYe)))

						pairRhoLu = (np.log((self.massMuon / self.massElectron * pairAstar * np.power(self.Z, -1. / 3.)
											 * np.sqrt((1. + 1. / pairRhoChi) * (1. + pairRhoYu)))
											/ (1. + (
									2. * self.massElectron * np.sqrt(self.euler) * pairAstar * np.power(self.Z,
																										-1. / 3.) * (
												1. + pairRhoChi) * (1. + pairRhoYu))
											   / (energy * pairEnergyTransfer * (1. - pairRhoRhoSquared))))
									 - np.log(
									1.5 * np.power(self.Z, 1. / 3.) * np.sqrt((1 + 1. / pairRhoChi) * (1 + pairRhoYu))))

						if (pairRhoChi >= 1E3):

							pairRhoBe = 1. / (2. * pairRhoChi) * (
										(3. - pairRhoRhoSquared) + 2. * pairBeta * (1. + pairRhoRhoSquared))

							pairRhoBu = (((1. + pairRhoRhoSquared) * (1. + 3. * pairBeta / 2.)
										  - 1. / pairRhoChi * (1. + 2. * pairBeta) * (1. - pairRhoRhoSquared)) * np.log(
								1. + pairRhoChi)
										 + (pairRhoChi * (1. - pairRhoRhoSquared - pairBeta)) / (1. + pairRhoChi) + (
													 1. + 2. * pairBeta) * (1. - pairRhoRhoSquared))

						elif (pairRhoChi <= 1E-3):

							pairRhoBe = (((2. + pairRhoRhoSquared) * (1. + pairBeta) + pairRhoChi * (
										3. + pairRhoRhoSquared))
										 * np.log(1. + 1. / pairRhoChi) + (1. - pairRhoRhoSquared - pairBeta) / (
													 1. + pairRhoChi) - (3. + pairRhoRhoSquared))

							pairRhoBu = pairRhoChi / 2. * (
										(5. - pairRhoRhoSquared) + pairBeta * (3. + pairRhoRhoSquared))

						else:

							pairRhoBe = (((2. + pairRhoRhoSquared) * (1. + pairBeta) + pairRhoChi * (
										3. + pairRhoRhoSquared))
										 * np.log(1. + 1. / pairRhoChi) + (1. - pairRhoRhoSquared - pairBeta) / (
													 1. + pairRhoChi) - (3. + pairRhoRhoSquared))

							pairRhoBu = (((1. + pairRhoRhoSquared) * (1. + 3. * pairBeta / 2.) - 1. / pairRhoChi * (
										1. + 2. * pairBeta) * (1. - pairRhoRhoSquared))
										 * np.log(1. + pairRhoChi) + (pairRhoChi * (1. - pairRhoRhoSquared - pairBeta))
										 / (1. + pairRhoChi) + (1. + 2. * pairBeta) * (1. - pairRhoRhoSquared))

						pairRhoPhie = pairRhoBe * pairRhoLe

						if (pairRhoPhie < 0):
							pairRhoPhie = 0.

						pairRhoPhiu = pairRhoBu * pairRhoLu

						if (pairRhoPhiu < 0):
							pairRhoPhiu = 0.

						pairGincremental = pairRhoPhie + np.power(self.massElectron / self.massMuon, 2.) * pairRhoPhiu

						pairG += pairRhoIntegrationWeight * pairRhoIntegrationSegmentLength * pairGincremental

					pairZetanum = (0.073 * np.log((energy / self.massMuon) / (
								1. + pairGamma1 * np.power(self.Z, 2. / 3.) * energy / self.massMuon)) - 0.26)

					if (pairZetanum < 0 or energy <= 35 * self.massMuon):
						pairZeta = 0.
					else:
						pairZeta = pairZetanum / (0.058 * np.log((energy / self.massMuon) / (
									1. + pairGamma2 * np.power(self.Z, 1. / 3.) * energy / self.massMuon)) - 0.14)

					pairCrosssection = (4. / (3. * self.pi) * self.Z * (self.Z + pairZeta) / self.A * self.avogadro
										* np.power(self.alpha * self.radiusElectron, 2.) * (
													1. - pairEnergyTransfer) * pairG)

					# For safety (Check for NaN)
					if (pairCrosssection != pairCrosssection):
						pairCrosssection = 0

					# Do pair production crosssection integration
					pairproductionLoss += pairIntegrationSegmentLength * pairIntegrationWeight * pairCrosssection

		return pairproductionLoss

	def photonuclear_loss_calc(self, energy):  # After Groom

		# Prepare variables
		photonuclearLoss = 0.
		photoIntegrationLowerBound = 0.
		photoIntegrationUpperBound = 1.
		photoIntervalBegin = -3.
		photoIntervalEnd = 3.
		photoIntegrationSegmentLength = (photoIntervalEnd - photoIntervalBegin) / self.numberIntegrationSegments
		photoM1 = 0.54
		photoM2 = 1.8
		photoMuonMassGeV = self.massMuon / 1000.
		photoEnergyGeV = energy / 1000.

		# Photonuclear Integration
		for segmentIterator in range(0, self.numberIntegrationSegments + 1):
			# exponential integration point and weights
			photoDoubleExponentialIntegrationPoint = np.tanh(
				self.pi / 2. * np.sinh(photoIntervalBegin + segmentIterator * photoIntegrationSegmentLength))

			photoFractionalEnergyTransfer = ((photoIntegrationLowerBound + photoIntegrationUpperBound) / 2
											 + (
														 photoIntegrationUpperBound - photoIntegrationLowerBound) / 2 * photoDoubleExponentialIntegrationPoint)

			photoIntegrationWeight = ((photoIntegrationUpperBound - photoIntegrationLowerBound) / 2
									  * (self.pi / 2 * np.cosh(
						photoIntervalBegin + segmentIterator * photoIntegrationSegmentLength))
									  / np.power(np.cosh(
						self.pi / 2 * np.sinh(photoIntervalBegin + segmentIterator * photoIntegrationSegmentLength)),
												 2))

			# Changing variables
			photoEnergyTransfer = photoFractionalEnergyTransfer * photoEnergyGeV

			photoSigmaGNu = 114.3 + 1.647 * np.power(np.log(0.0213 * photoEnergyTransfer), 2.)

			photoXg = 0.00282 * np.power(self.A, 1. / 3.) * photoSigmaGNu

			photoG = 3. / np.power(photoXg, 3.) * (photoXg * photoXg / 2. - 1. + np.exp(-photoXg) * (1. + photoXg))

			photoKappa = 1. - 2. / photoFractionalEnergyTransfer + 2. / np.power(photoFractionalEnergyTransfer, 2.)

			photoT = photoMuonMassGeV * photoMuonMassGeV * photoFractionalEnergyTransfer * photoFractionalEnergyTransfer / (
						1. - photoFractionalEnergyTransfer)

			photoTerm1 = (0.75 * photoG * (
						photoKappa * np.log(1. + photoM1 / photoT) - (photoKappa * photoM1) / (photoM1 + photoT)
						- (2. * photoMuonMassGeV * photoMuonMassGeV) / photoT))
			#                  *(photoM1 + 2. * photoT)
			# / (photoM1 + photoT) + 4. * photoMuonMassGeV*photoMuonMassGeV / photoM1*np.log(1. + photoM1 / photoT)))

			photoTerm2 = 0.25 * (
						photoKappa * np.log(1. + photoM2 / photoT) - 2. * photoMuonMassGeV * photoMuonMassGeV / photoT)
			#                 (((photoKappa + 2. * photoMuonMassGeV*photoMuonMassGeV / photoM2)*np.log(1. + photoM2 / photoT)
			# - (2. * photoMuonMassGeV*photoMuonMassGeV) / photoT) / 4.)

			photoTerm3 = (photoMuonMassGeV * photoMuonMassGeV / (2. * photoT)
						  * (0.75 * photoG * photoM1 / (photoM1 + photoT) + 0.25 * photoM2 / photoT * np.log(
						1. + photoT / photoM2)))

			# Calculate Crosssection
			photoCrosssection = (self.alpha / (2. * self.pi) * self.A * photoSigmaGNu * photoFractionalEnergyTransfer
								 * (photoTerm1 + photoTerm2 + photoTerm3))

			photoIntegrand = (self.avogadro * 1e-30 / self.A) * photoFractionalEnergyTransfer * photoCrosssection

			# Integrate
			photonuclearLoss += photoIntegrationWeight * photoIntegrationSegmentLength * photoIntegrand

		if (photonuclearLoss < 0):
			photonuclearLoss = 0

		return photonuclearLoss


class Oxides(Material):
	def __init__(self, name):
		super().__init__(name)
		pass


class MixtureGUI(Material):
	def __init__(self, name, groom_df, parent=None, matpath=None, denmode=None, den_dict=None):
		super().__init__(name)

		self.groom_df = groom_df

		ele_dict = groom_df.loc[groom_df['Name'] == name].squeeze().to_dict()
		self.ele_dict = ele_dict

		self.formula = ele_dict['Formula']
		self.ZoverA = float(ele_dict['Z/A'])
		self.rho0 = ele_dict['Density']
		self.I = ele_dict['I']
		self.C = ele_dict['Cbar']
		self.x0 = ele_dict['x0']
		self.x1 = ele_dict['x1']
		self.a = ele_dict['a']
		self.k = ele_dict['k']
		self.delta = ele_dict['dlt0']
		self.density = self.rho0
		self.neles = len(self.ele_dict['Constituents'])

	def make_groom_table(self, spath):
		# Get number of elements
		ele_list = self.ele_dict['Constituents']

		# Set up energies from 10 MeV - 100 TeV
		energies = [10e0,14e0,20e0,30e0,40e0,80e0,
					10e1,14e1,20e1,30e1,40e1,80e1,
					10e2,14e2,20e2,30e2,40e2,80e2,
					10e3,14e3,20e3,30e3,40e3,80e3,
					10e4,14e4,20e4,30e4,40e4,80e4,
					10e5,14e5,20e5,30e5,40e5,80e5,
					10e6,14e6,20e6,30e6,40e6,80e6,
					10e7]


		popup = tk.Toplevel(takefocus=True)
		proglbl = tk.Label(popup, text='Progress')
		eleprogbar = ttk.Progressbar(popup, length=200, maximum=self.neles)
		progbar = ttk.Progressbar(popup, length=200, maximum=np.alen(energies))

		eletxt = tk.Label(popup, text=str(0) + '/' + str(self.neles))
		progtxt = tk.Label(popup, text=str(0) + ' %')
		popup.update_idletasks()

		proglbl.grid(column=0, row=0, columnspan=2)
		eleprogbar.grid(column=0, row=1, padx=10, pady=10)
		eletxt.grid(column=1, row=1, padx=10, pady=10)
		progbar.grid(column=0, row=2, padx=10, pady=10)
		progtxt.grid(column=1, row=2, padx=10, pady=10)
		popup.update()


		ioniTable = np.zeros_like(energies)
		bremTable = np.zeros_like(energies)
		pairTable = np.zeros_like(energies)
		photTable = np.zeros_like(energies)

		for ide, engy in enumerate(energies):
			ioniTable[ide] = (self.ionisation_loss(engy))


		for iele, ele in enumerate(ele_list):

			eletxt.config(text=str(iele+1) + '/' + str(self.neles))
			this_ele_dict = self.groom_df.iloc[ele[0]-1]
			ele_obj = ElementGUI(this_ele_dict['Formula'], this_ele_dict)

			for ide, engy in enumerate(energies):
				progbar.step(1)
				frac = int(ide / len(energies) * 100)

				progtxt.config(text=str('{:2.2f}'.format(ide/np.alen(energies)*100.)) + '%')
				popup.update()

				bremTable[ide] += (ele_obj.bremsstrahlung_loss_calc(engy))*ele[1]
				pairTable[ide] += (ele_obj.pairproduction_loss_calc(engy))*ele[1]
				photTable[ide] += (ele_obj.photonuclear_loss_calc(engy))*ele[1]

			eleprogbar.step(1)
			popup.update()


		table = pd.DataFrame({'Energy': energies, 'Ionization': ioniTable, 'Bremsstrahlung': bremTable,
							  'Pair-Production': pairTable, 'Photonuclear': photTable})

		table.to_csv(os.path.join(spath, self.name + '_GT.csv'))

		popup.destroy()


	def create_radloss_table_gui(self, spath):
		# Get number of elements
		ele_list = self.ele_dict['Constituents']

		# Set up energies from 10 MeV - 100 TeV
		energies = np.logspace(1.0, 8.0, num=100, base=10.0)


		popup = tk.Toplevel(takefocus=True)
		proglbl = tk.Label(popup, text='Progress')
		eleprogbar = ttk.Progressbar(popup, length=200, maximum=self.neles)
		progbar = ttk.Progressbar(popup, length=200, maximum=np.alen(energies))

		eletxt = tk.Label(popup, text=str(0) + '/' + str(self.neles))
		progtxt = tk.Label(popup, text=str(0) + ' %')
		popup.update_idletasks()

		proglbl.grid(column=0, row=0, columnspan=2)
		eleprogbar.grid(column=0, row=1, padx=10, pady=10)
		eletxt.grid(column=1, row=1, padx=10, pady=10)
		progbar.grid(column=0, row=2, padx=10, pady=10)
		progtxt.grid(column=1, row=2, padx=10, pady=10)
		popup.update()

		rlTable = np.zeros_like(energies)

		for iele, ele in enumerate(ele_list):

			eletxt.config(text=str(iele+1) + '/' + str(self.neles))

			this_ele_dict = self.groom_df.iloc[ele[0]-1]

			ele_obj = ElementGUI(this_ele_dict['Formula'], this_ele_dict)

			for ide, engy in enumerate(energies):
				progbar.step(1)
				progtxt.config(text=str('{:2.2f}'.format(ide/np.alen(energies)*100.)) + '%')
				popup.update()
				rlTable[ide] += (ele_obj.radiation_loss_calc(engy))*ele[1]

			eleprogbar.step(1)
			popup.update()

		popup.destroy()

		self.rlTable = np.vstack((energies, rlTable)).T
		savefile = os.path.normpath(os.path.join(spath, self.name + '_RLT.txt'))
		np.savetxt(savefile, self.rlTable)


	def ionisation_loss(self, energy):
		# Define kinetic variables
		gamma = (energy + self.massMuon) / self.massMuon
		gammaSquared = gamma * gamma
		betaSquared = ((energy * (2. * self.massMuon + energy))
					   / ((energy + self.massMuon) * (energy + self.massMuon)))
		beta = betaSquared ** 0.5
		betagamma = beta * gamma
		betagammaSquared = betagamma * betagamma
		xCompare = np.log10(betagamma)

		# Dfeine Ionisation Loss variables
		qMax = ((2. * self.massElectron * betagammaSquared)
				/ (1. + 2. * gamma * self.massElectron / self.massMuon
				   + (self.massElectron / self.massMuon) * (self.massElectron / self.massMuon)))

		correction = (self.K / (4. * self.pi) * self.ZoverA * self.alpha
					  * (np.log(2. * (energy + self.massMuon) / self.massMuon)
						 - np.log(2 * qMax / self.massElectron) / 3) * np.power(np.log(2 * qMax / self.massElectron),
																				2))

		shellCorrection = ((0.422377 / betagammaSquared + 0.0304043 / (betagammaSquared * betagammaSquared)
							- 0.00038106 / (
									betagammaSquared * betagammaSquared * betagammaSquared)) * 1E-6 * self.I * self.I
						   + (3.858016 / betagammaSquared - 0.1667989 / (betagammaSquared * betagammaSquared)
							  + 0.00157955 / (
									  betagammaSquared * betagammaSquared * betagammaSquared)) * 1E-9 * self.I * self.I * self.I)

		if (xCompare < self.x0):
			delta = self.delta * np.power(10., 2. * (xCompare - self.x0))
		elif (xCompare > self.x1):
			delta = 2 * np.log(10) * xCompare - self.C
		else:
			delta = 2 * np.log(10) * xCompare - self.C + self.a * np.power(self.x1 - xCompare, self.k)

		# Return Ionisation dE/dX
		first = 0.5 * np.log((2. * self.massElectron * betagammaSquared * qMax * 1E12) / (self.I * self.I))
		second = - betaSquared
		third = - delta / 2.
		fourth = (qMax * qMax) / (8. * (gamma * gamma * self.massMuon * self.massMuon))
		fifth = 0
		sixth = correction

		result = self.K / betaSquared * self.ZoverA * (first + second + third + fourth + fifth) + sixth

		return result


class ElementGUI(Material):
	def __init__(self, name, ele_dict, parent=None, matpath=None, denmode=None, den_dict=None):
		super().__init__(name)

		self.formula = ele_dict['Formula']
		self.Z = float(ele_dict['Z'])
		self.A = float(ele_dict['A'])
		self.rho0 = ele_dict['Density']
		self.I = ele_dict['I']
		self.C = ele_dict['Cbar']
		self.x0 = ele_dict['x0']
		self.x1 = ele_dict['x1']
		self.a = ele_dict['a']
		self.k = ele_dict['k']
		self.delta = ele_dict['dlt0']
		self.density = self.rho0

	def create_radloss_table_gui(self, spath):
		popup = tk.Toplevel(takefocus=True)
		proglbl = tk.Label(popup, text='Progress')
		progbar = ttk.Progressbar(popup, length=200, maximum=100)
		progtxt = tk.Label(popup, text=str(0) + ' %')
		popup.update_idletasks()

		proglbl.grid(column=0, row=0, columnspan=2)
		progbar.grid(column=0, row=1, padx=10, pady=10)
		progtxt.grid(column=1, row=1, padx=10, pady=10)
		popup.update()

		# Set up energies from 10 MeV - 100 TeV
		energies = np.logspace(1.0, 8.0, num=100, base=10.0)

		rlTable = np.zeros_like(energies)

		for ide, engy in enumerate(energies):
			progbar.step(1)
			progtxt.config(text=str(ide) + '%')
			popup.update()
			rlTable[ide] = (self.radiation_loss_calc(engy))

		self.rlTable = np.vstack((energies, rlTable)).T
		savefile = os.path.normpath(os.path.join(spath, self.name + '_RLT.txt'))
		np.savetxt(savefile, self.rlTable)

		popup.destroy()

	def make_groom_table(self, spath):
		# Set up energies from 10 MeV - 100 TeV
		energies = [10e0,14e0,20e0,30e0,40e0,80e0,
					10e1,14e1,20e1,30e1,40e1,80e1,
					10e2,14e2,20e2,30e2,40e2,80e2,
					10e3,14e3,20e3,30e3,40e3,80e3,
					10e4,14e4,20e4,30e4,40e4,80e4,
					10e5,14e5,20e5,30e5,40e5,80e5,
					10e6,14e6,20e6,30e6,40e6,80e6,
					10e7]

		popup = tk.Toplevel(takefocus=True)
		proglbl = tk.Label(popup, text='Progress')
		progbar = ttk.Progressbar(popup, length=200, maximum=len(energies))
		progtxt = tk.Label(popup, text=str(0) + ' %')
		popup.update_idletasks()

		proglbl.grid(column=0, row=0, columnspan=2)
		progbar.grid(column=0, row=1, padx=10, pady=10)
		progtxt.grid(column=1, row=1, padx=10, pady=10)
		popup.update()


		ioniTable = np.zeros_like(energies)
		bremTable = np.zeros_like(energies)
		pairTable = np.zeros_like(energies)
		photTable = np.zeros_like(energies)

		for ide, engy in enumerate(energies):
			progbar.step(1)
			frac = int(ide/len(energies)*100)

			progtxt.config(text=str(frac) + '%')
			popup.update()
			ioniTable[ide] = (self.ionisation_loss(engy))
			bremTable[ide] = (self.bremsstrahlung_loss_calc(engy))
			pairTable[ide] = (self.pairproduction_loss_calc(engy))
			photTable[ide] = (self.photonuclear_loss_calc(engy))

		table = pd.DataFrame({'Energy': energies, 'Ionization': ioniTable, 'Bremsstrahlung': bremTable,
							  'Pair-Production': pairTable, 'Photonuclear': photTable})

		table.to_csv(os.path.join(spath, self.name + '_GT.csv'))

		popup.destroy()


class Mixture(Material):
	def __init__(self, name, mix_dict, parent, matpath, denmode, den_dict):
		super().__init__(name)

		self.parent = parent

		#print(name, mix_dict, parent)

		# Initialise parameters
		self.name = name
		self.filename = mix_dict['Name']
		self.ZoverA = float(mix_dict['Z/A'])
		self.rho0 = mix_dict['Density']
		self.I = mix_dict['I']
		self.C0 = mix_dict['Cbar']
		self.x00 = mix_dict['x0']
		self.x10 = mix_dict['x1']
		self.a = mix_dict['a']
		self.k = mix_dict['k']
		self.delta = mix_dict['dlt0']
		self.density = self.rho0
		self.x0 = self.x00
		self.x1 = self.x10
		self.C = self.C0

		self.rlTable = np.empty(())
		self.set_radloss_table(matpath)

		self.parameters = LocalParameterGroup()
		self.parent.parameters.declare_child(self.parameters)

		# Define lp for densities and set initial density model
		denpar = Parameter(self.name + '_' + chr(0x03C1), 'density')

		# print(den_dict['mean'], den_dict['std'])

		if denmode == 'Normal':
			denpar.set_pdf(sps.norm, loc=den_dict['mean'], scale=den_dict['std'])
			denpar.set_value(sps.norm.rvs(loc=den_dict['mean'], scale=den_dict['std']))
			denpar.set_inistd(den_dict['std'])
			self.parameters.add_Parameter(denpar)
			self.parent.parameters.add_Parameter(denpar)
		elif denmode == 'Log-normal':
			pass
		elif denmode == 'From data':
			pass

		#print(self.parameters.parameter_df)

	def set_radloss_table(self, spath):
		tabloc = os.path.join(spath, self.filename + '_RLT.txt')
		self.rlTable = np.loadtxt(tabloc)

	def ionisation_loss(self, energy):
		# Define kinetic variables
		gamma = (energy + self.massMuon) / self.massMuon
		gammaSquared = gamma * gamma
		betaSquared = ((energy * (2. * self.massMuon + energy))
					   / ((energy + self.massMuon) * (energy + self.massMuon)))
		beta = betaSquared ** 0.5
		betagamma = beta * gamma
		betagammaSquared = betagamma * betagamma
		xCompare = np.log10(betagamma)

		# Dfeine Ionisation Loss variables
		qMax = ((2. * self.massElectron * betagammaSquared)
				/ (1. + 2. * gamma * self.massElectron / self.massMuon
				   + (self.massElectron / self.massMuon) * (self.massElectron / self.massMuon)))

		correction = (self.K / (4. * self.pi) * self.ZoverA * self.alpha
					  * (np.log(2. * (energy + self.massMuon) / self.massMuon)
						 - np.log(2 * qMax / self.massElectron) / 3) * np.power(np.log(2 * qMax / self.massElectron),
																				2))

		shellCorrection = ((0.422377 / betagammaSquared + 0.0304043 / (betagammaSquared * betagammaSquared)
							- 0.00038106 / (
									betagammaSquared * betagammaSquared * betagammaSquared)) * 1E-6 * self.I * self.I
						   + (3.858016 / betagammaSquared - 0.1667989 / (betagammaSquared * betagammaSquared)
							  + 0.00157955 / (
									  betagammaSquared * betagammaSquared * betagammaSquared)) * 1E-9 * self.I * self.I * self.I)

		if (xCompare < self.x0):
			delta = self.delta * np.power(10., 2. * (xCompare - self.x0))
		elif (xCompare > self.x1):
			delta = 2 * np.log(10) * xCompare - self.C
		else:
			delta = 2 * np.log(10) * xCompare - self.C + self.a * np.power(self.x1 - xCompare, self.k)

		# Return Ionisation dE/dX
		first = 0.5 * np.log((2. * self.massElectron * betagammaSquared * qMax * 1E12) / (self.I * self.I))
		second = - betaSquared
		third = - delta / 2.
		fourth = (qMax * qMax) / (8. * (gamma * gamma * self.massMuon * self.massMuon))
		fifth = 0
		sixth = correction

		result = self.K / betaSquared * self.ZoverA * (first + second + third + fourth + fifth) + sixth

		return result

	def calculate_logprob(self):
		return self.parameters.calculate_logprob()

	def update(self):
		self.update_density(self.parameters.parameters[self.name + '_' + chr(0x03C1)].value)

	def update_density(self, density):
		self.density = density
		r = self.density / self.rho0
		self.C = self.C0 - np.log(r)
		self.x0 = self.x00 - 0.5 * np.log10(r)
		self.x1 = self.x10 - 0.5 * np.log10(r)


class Element(Material):
	rlTable: ndarray

	def __init__(self, name, ele_dict, parent, matpath, denmode, den_dict):
		super().__init__(name)

		self.parent = parent
		#print(ele_dict)

		# Initialise parameters
		self.name = name
		self.formula = ele_dict['Formula']
		self.Z = float(ele_dict['Z'])
		self.A = float(ele_dict['A'])
		self.rho0 = ele_dict['Density']
		self.I = ele_dict['I']
		self.C0 = ele_dict['Cbar']
		self.x00 = ele_dict['x0']
		self.x10 = ele_dict['x1']
		self.a = ele_dict['a']
		self.k = ele_dict['k']
		self.delta = ele_dict['dlt0']
		self.density = self.rho0
		self.x0 = self.x00
		self.x1 = self.x10
		self.C = self.C0
		# self.logprob = 0

		self.rlTable = np.empty(())
		self.set_radloss_table(matpath)

		self.parameters = LocalParameterGroup()
		self.parent.parameters.declare_child(self.parameters)

		# Define lp for densities and set initial density model
		denpar = Parameter(name + '_' + chr(0x03C1), 'density')
		if denmode == 'Normal':
			denpar.set_pdf(sps.norm, loc=den_dict['mean'], scale=den_dict['std'])
			denpar.set_value(sps.norm.rvs(loc=den_dict['mean'], scale=den_dict['std']))
			denpar.set_inistd(den_dict['std'])
			self.parameters.add_Parameter(denpar)
			self.parent.parameters.add_Parameter(denpar)
		elif denmode == 'Log-normal':
			pass
		elif denmode == 'From data':
			pass

		self.update()

	def update(self):
		self.update_density(self.parameters.parameters[self.name + '_' + chr(0x03C1)].value)

	def update_density(self, density):
		self.density = density
		r = self.density / self.rho0
		self.C = self.C0 - np.log(r)
		self.x0 = self.x00 - 0.5 * np.log10(r)
		self.x1 = self.x10 - 0.5 * np.log10(r)

	def set_radloss_table(self, matdirectory):
		tabloc = os.path.join(matdirectory, self.formula + '_RLT.txt')
		self.rlTable = np.loadtxt(tabloc)

	def update_logprob(self):
		self.logprob = 0
		self.logprob += self.lp_density(self.density)

	def get_parameters(self):
		return self.parameters.iloc[[-1]]

	def get_logprob(self):
		return self.logprob

	def calculate_logprob(self):
		return self.parameters.calculate_logprob()

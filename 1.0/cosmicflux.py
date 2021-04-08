import numpy as np


def getVerticalFlux (momentum:float) -> float:
	y = np.log10(momentum)
	y2 = y*y
	y3 = y*y2
	exponent = -(0.0209*y3-0.2555*y2+1.288*y+0.2455)
	
	result = 0.00253*np.power(momentum,exponent)

	return result

def getDifferentialFlux (momentum:float, height:float, zenithAngle:float) -> float:
	cosTheta = np.cos(zenithAngle)
	redmom = momentum*cosTheta
	inclinedFlux = getVerticalFlux(momentum*cosTheta)*cosTheta*cosTheta*cosTheta
	result = inclinedFlux*np.exp(height/(3400.+1100.*momentum*cosTheta))

	return result # returns in [cm-2 s-1 sr-1]/[GeV]
	#return result * 10	# returns in [m-2 s-1 sr-1]/[MeV]

def getIntegratedFlux(lowerEnergy:float, height:float, zenithAngle:float) -> float:
	flux = 0.0
	lowerIntegrationBoundary = lowerEnergy / 1000. # Convert MeV to GeV
	fluxIntervalBegin = -3.
	fluxIntervalEnd = 4
	fluxIntegrationSegmentLength = (fluxIntervalEnd - fluxIntervalBegin) / 200.

	for i in range(0,201):
		t = fluxIntervalBegin+i*fluxIntegrationSegmentLength
		fluxDEIntegrationPoint = np.exp(np.pi/2*np.sinh(t))
		fluxEnergyPoint = fluxDEIntegrationPoint + lowerIntegrationBoundary
		fluxIntegrationWeight = np.pi/2*np.cosh(t)*fluxDEIntegrationPoint   

		if (fluxEnergyPoint < 0.10567):
			momentum = 0
			flincr = 0
		else:
			momentum = np.sqrt(fluxEnergyPoint*fluxEnergyPoint-0.10567*0.10567)
			flincr = getDifferentialFlux(momentum, height, zenithAngle)*fluxIntegrationWeight*fluxIntegrationSegmentLength

		#print(flincr)
		flux = flux + flincr

	#return flux * 10000   # returns [m-2 s-1 sr-1]
	#return flux * 60*60*24  # returns [cm-2 day-1 sr-1]
	return flux # returns [cm-2 s-1 sr-1]

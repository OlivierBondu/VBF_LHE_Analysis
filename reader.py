from math import log, tan, acos, pi
import gzip

class Event(object):
	def __init__(self, particleList, status):
		self.plist = particleList
		self.status = status

	def getByPid(self, pid, isAbs=False, isNot=False):
		plistByPid = []
		for particle in self.plist:
			if not isAbs and not isNot and particle.pid == pid:
				plistByPid.append(particle)
			if isAbs and not isNot and abs(particle.pid) == pid:
				plistByPid.append(particle)
			if isAbs and isNot and abs(particle.pid) != pid:
				plistByPid.append(particle)
			if not isAbs and isNot and particle.pid != pid:
				plistByPid.append(particle)
		return plistByPid



class LorentzVector( object ):
	def __init__(self, pid, listPxPyPzE ):
		self.pid = pid
		self.px = listPxPyPzE[0]
		self.py = listPxPyPzE[1]
		self.pz = listPxPyPzE[2]
		self.e = listPxPyPzE[3]
		self.p2 = self.px**2 + self.py**2 + self.pz**2
		self.m2 = self.e**2 - self.p2
		self.p = self.p2 ** 0.5
		self.pt = (self.px**2 + self.py**2)**0.5 
		try:
			self.m = self.m2 ** 0.5
		except ValueError:
			self.m = 0
		self.theta = acos(self.pz / self.p)
		self.phi = acos(self.py / self.pt)
		self.eta = -log(tan(self.theta / 2.))

	def listPxPyPzE(self):
		return [self.px, self.py, self.pz, self.e]

def mass(a, b):
	return ((a.e + b.e)**2 - ((a.px + b.px)**2 + (a.py + b.py)**2 + (a.pz + b.pz)**2))**.5

def dr(a, b):
	deta = a.eta - b.eta
	dphi = a.phi - b.phi
	if dphi > pi:
		dphi -= pi
	if dphi < -pi:
		dphi += pi
	return (deta**2 + dphi**2)**.5

lheDir="files/"

def getFinalParticles( lheFile ):
	data = {}
	processingEvent = False
	particleList = []
	ievent = 0
	for line in lheFile:
#		if "<event>" in line: # standard lhe format
#			processingEvent = True
#		if "</event>" in line: # standard lhe format
#			processingEvent = False
#			data[ievent] = Event(particleList, 0)
#			ievent += 1
#			particleList = []			
#		if processingEvent:
#			fullParticle = line.split()
#			if len(fullParticle) >= 13 and int(fullParticle[1]) == 1:
#				particle = LorentzVector( int(fullParticle[0]), map(float, fullParticle[6:10]) )
#				particleList.append(particle)
		if "#" in line: # non-standard Xanda's lhe format
			processingEvent = True
		if "8" == line.rstrip() and len(particleList) > 0: # non-standard Xanda's lhe format
			processingEvent = False
			data[ievent] = Event(particleList, 0)
			ievent += 1
			particleList = []			
		if processingEvent:
			fullParticle = line.split()
#			if len(fullParticle) >= 13 and int(fullParticle[1]) == 1:
#				particle = LorentzVector( int(fullParticle[0]), map(float, fullParticle[6:10]) )
#				particleList.append(particle)
			if len(fullParticle) >= 5 and int(fullParticle[0]) != 25:
				particle = LorentzVector( int(fullParticle[0]), map(float, fullParticle[1:5]) )
#				print map(float, fullParticle[1:5])
				particleList.append(particle)
	return data

# ##### ##### #####	#####
# START THE MAIN PROCESS
# ##### ##### #####	#####

process = "zbbbbjj"
for irun in xrange(1,21):
#		fileName = process + "_10k_run" + str(irun) + "_unweighted_events.lhe"
#		fileName = "Events_4b2j/13tev/parton/pp_hh_vbf_SM_13tev_nocuts.lhe"
	fileName = "MGraviton_850.lhe.decayed" 
	print "reading file ", fileName
	with open(lheDir + fileName, 'r') as lheFile:
		data = getFinalParticles( lheFile )
	print "nProcessed: ", len(data)
#		for particle in data[0].plist:
#			print particle.pid, particle.m, particle.pt, particle.eta

	for ievent in data:
		# Acceptance cuts
#		for particle in data[ievent].plist:
#			if particle.pt < 30 or abs(particle.eta) > 4.5:
#				data[ievent].status = 1
#		if data[ievent].status: continue
		for particle in data[ievent].getByPid(5, isAbs=True):
			if particle.pt < 10 or abs(particle.eta) > 2.5:
				data[ievent].status = 2

		listOfNonB = data[ievent].getByPid(5, isAbs=True, isNot=True)
		listOfB = data[ievent].getByPid(5, isAbs=True)
		# M(jj) > 400 GeV
		if data[ievent].status: continue
		p1 = data[ievent].plist[0]; p2 = data[ievent].plist[0]
		p1 = listOfNonB[0]
		p2 = listOfNonB[1]
		if mass(p1, p2) < 400:
			data[ievent].status = 3
			continue
		# |Deta(jj)| < 4
		if abs(p1.eta - p2.eta) > 4:
			data[ievent].status = 4
			continue

		# prepare the two higgs candidates
		Higgs1 = [0, -1]
		Higgs2 = [-1, -1]
		jet1 = listOfB[0]
		minmass = 14000
		for ijet, jet2 in enumerate(listOfB, 1):
			for jjet, jet3 in enumerate(listOfB, ijet+1):
				for kjet, jet4 in enumerate(listOfB, jjet+1):
					if abs(mass(jet1, jet2) - mass(jet3, jet4)) < minmass:
						minmass = abs(mass(jet1, jet2) - mass(jet3, jet4))
						Higgs1[1] = ijet
						Higgs2[0] = jjet
						Higgs2[1] = kjet

#FIXME			
#			higgs1 = Event(25, [])	

		# Each higgs candidate must lie within 10% of the Higgs mass
		if not 0.9 * 125 < mass(listOfB[Higgs1[0]], listOfB[Higgs1[1]]) < 1.1 * 125:
			data[ievent].status = 5
			continue
		if not 0.9 * 125 < mass(listOfB[Higgs2[0]], listOfB[Higgs2[1]]) < 1.1 * 125:
			data[ievent].status = 6
			continue

	ncuts = 6
	nevents = [len(data) for x in range(ncuts+1)]
	for cut in range(1,ncuts+1):
		nevents[cut] = nevents[cut-1]
		for ievent in data:
			if data[ievent].status == cut:
				nevents[cut] -= 1

	print "cutflow: icut, nevent, efficiency, expected yield @ 3ab-1"
	for icut, nevent in enumerate(nevents):
		print icut, nevent, float(nevent) / nevents[0], nevent * 0.8222 * 3000 / nevents[0]
	break
#	break

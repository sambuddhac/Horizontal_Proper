#Definition for Marketover class public Member Methods
import numpy as np
from Python_src.nettran import Nettran
from Python_src.log import log
from Python_src.profiler import Profiler
#define BIGM 1000000000000000000 // Defines the value of the Big M for transforming the bilinear terms to linear constraints

class Marketover(object):
	def __init__(self, zoneCount, totNodes, SECount, candCount, vectorOfSubnet, listOfSharedNode, listOfGlobalSharedNode, listOfSharedZone, ListOfSESerial, ListOfSEFromRank, ListOfSEToRank, ListOfSEReact, ListOfSECap, ListOfCandSerial, ListOfCandFromRank, ListOfCandToRank, ListOfCandReact, ListOfCandCap, milpAlgoChoice, scenarios, compareBase): #constructor
		self.zoneNumber = zoneCount #Initialize the total number of load zones
		self.subNetVector = vectorOfSubnet #Create handle to the vector of subnetworks
		self.zoneList = listOfSharedZone #List of zones between which exsiting and candidate shared lines exist
		self.nodeList = listOfSharedNode #List of nodes at the ends of the existing and candidate shared lines
		self.sharedGlobalList = listOfGlobalSharedNode 
		self.nodeNumber = totNodes #total number of nodes, which form the ends of the shared lines
		self.candLineNumber = candCount #Number of shared candidate lines
		self.sharedELineNumber = SECount #Number of shared existing Transmission lines
		self.SESerial = ListOfSESerial #Serial list of shared existing lines
		self.SEFromRank = ListOfSEFromRank #List of rank of from nodes of shared existing lines
		self.SEToRank = ListOfSEToRank #List of rank of to nodes of shared existing lines
		self.SEReactance = ListOfSEReact #List of reactances of shared existing lines
		self.SECapacity = ListOfSECap #List of line flow limits of shared existing lines
		self.candSerial = ListOfCandSerial #Serial list of shared candidate lines
		self.candFromRank = ListOfCandFromRank #List of rank of from nodes of shared candidate lines
		self.candToRank = ListOfCandToRank #List of rank of to nodes of shared candidate lines
		self.candReactance = ListOfCandReact #List of reactances of shared candidate lines
		self.candCapacity = ListOfCandCap #List of line flow limits of shared candidate lines
		self.lpSolveAlgo = milpAlgoChoice #Simplex for 1 and IPM for 2
		self.countOfScenarios = scenarios #Number of total random scenarios
		self.compareBasis = compareBase
		self.lineInterDecVec = np.zeros(int, float)
		self.phAngDecVec = np.zeros(int, float)
		"""
		int skip = 0;
		for (candserialiterator = candSerial.begin(); candserialiterator != candSerial.end(); ++candserialiterator) {
			if (skip > 0) {
				cout << "Candidate Serial is " << *candserialiterator << endl;
			}
			++skip;
		}
		for (candfromiterator = candFromRank.begin(); candfromiterator != candFromRank.end(); ++candfromiterator) {
			cout << "From rank is " << *candfromiterator << endl;
		}
		for (candtoiterator = candToRank.begin(); candtoiterator != candToRank.end(); ++candtoiterator) {
			cout << "To rank is " << *candtoiterator << endl;
		}
		"""
		#end constructor

	#def __del__(): #destructor
		#log.info("\nSimulation ended")
		#destructor ends

	def MILPMarketover(self, LagMultXi, LagMultPi, totalCandLineNum, totalSharedNodeNum): #Function MILPMarketover() implements the Mixed Integer Linear Programming Solver routine by calling GLPK routines for average heat rate objective
		#CREATION OF THE MIP SOLVER INSTANCE */
		dimRow = self.countOfScenarios*(4 * self.candLineNumber + 2 * self.sharedELineNumber) #Total number of rows of the A matrix (number of structural constraints of the LP) first term to account for lower and upper line limits & lower and upper definition limits of candidate shared lines, and second term for lower and upper line limits for shared existing lines
		dimCol = self.countOfScenarios*self.nodeNumber+(self.countOfScenarios+1)*self.candLineNumber #Total number of columns of the LP (number of Decision Variables) first term to account for voltage phase angles for inter-zonal lines' nodes, and second term for power flow values and binary integer decision variable values for shared candidate lines

	#Function MILP() ends

	def LBMarketover(self, LagMultXi, LagMultPi, totalCandLineNum, totalSharedNodeNum): #Function LBMarketover() calculates the lower bound of the Mixed Integer Linear Programming Solver routine by calling GLPK routines for average heat rate objective
		dimRow = self.countOfScenarios*(4 * self.candLineNumber + 2 * self.sharedELineNumber) #Total number of rows of the A matrix (number of structural constraints of the LP) first term to account for lower and upper line limits & lower and upper definition limits of candidate shared lines, and second term for lower and upper line limits for shared existing lines
		dimCol = self.countOfScenarios*self.nodeNumber+(self.countOfScenarios+1)*self.candLineNumber #Total number of columns of the LP (number of Decision Variables) first term to account for voltage phase angles for inter-zonal lines' nodes, and second term for power flow values and binary integer decision variable values for shared candidate lines
	#Function MILP() ends

	def bufferintermediateDecision(self, iterCountOut):
		if iterCountOut==0: 
			for scenarioTrack in range(self.countOfScenarios): 
				for nodeTrack in range(self.nodeNumber):
					self.interimContDecVarPrev.append(0)
	
			for lineTrack in range(self.candLineNumber):
				self.interimIntDecVarPrev.append(0)
		else: 
			interimContDecVarPrev = self.interimContDecVar
			interimIntDecVarPrev = self.interimIntDecVar
	

	def getGlobalUpper(self, LagMultXi, LagMultPi, regionalUpper, numberOfZones): #Function getGlobalUpper() returns the global upper bound for the investment coordination problem
		#Sum up the regional upper bounds, which are the tentative regional minima, at the end of every iteration
		revisedUpperBound = 0 #total upper bound initialized
		for i in range(numberOfZones):
			revisedUpperBound += regionalUpper[i]
		#Calculate the MO upper bound, which is the minimum value of the MO objective, calculated with a given mix of the intermediate, shared continuous variable values from different regions/zones
		for scenPos in range(self.countOfScenarios):
			rankPresChecker = [] #Vector to check if a particular rank has been accounted for 
			length = len(self.angleDecIndex[scenPos]) #length of the rank vector
			for i in range(length): #Put as many zeroes on the rankPresChecker vector as is the length of the rank vector 
				rankPresChecker.append(0)
			interAngleTermVec = [] #Intermediate vector for storing the costs of the angle terms from each sub-network
			tracker = 0 #tracker to track the position of the angleDecIndexIterator
			for angleDecIndexIterator in self.angleDecIndex[scenPos]: #Iterate through rank vector
				interNodeRank = angleDecIndexIterator #Store the value of the rank of the present node iterate in the list 
				if rankPresChecker[tracker] == 0: #If this node rank hasn't been already accounted for
					indexList = [pos1 for pos1 in range(len(self.angleDecIndex[scenPos])) if self.angleDecIndex[scenPos][pos1] == interNodeRank] # Get all the indices of interNodeRank in the angleDecIndex vector
					for indexListIterator in indexList: # while all the different positions of this rank hasn't been accounted for
						rankPresChecker[indexListIterator] = 1 #Mark the element in the corresponding position of checker vector 1, to indicate this rank has been accounted for 
						interAngleTermVec.append(-(self.phaseAngleDecision[scenPos][indexListIterator])*(LagMultXi[scenPos*self.nodeNumber+(angleDecIndexIterator-1)])) #Calculate cost term
					smallest_element = min(interAngleTermVec)
					revisedUpperBound += smallest_element
					interAngleTermVec = []
					indexList = []
				tracker += 1 #Increment the tracker
		#Calculate the MO upper bound, which is the minimum value of the MO objective, calculated with a given mix of the intermediate, shared discrete variable values from different regions/zones
		rankPresCheckerInt = [] #Vector to check if a particular rank has been accounted for 
		lengthInt = len(self.lineDecIndex) #length of the rank vector
		for i in range(lengthInt):#Put as many zeroes on the rankPresChecker vector as is the length of the rank vector 
			self.rankPresCheckerInt.append(0)
		interLineTermVec = [] #Intermediate vector for storing the costs of the line building decisions from each sub-network
		trackerInt = 0 #tracker to track the position of the lineDecIndexIterator
		for lineDecIndexIterator in self.lineDecIndex: #Iterate through rank vector
			interLineRank = lineDecIndexIterator #Store the value of the rank of the present line iterate in the list 
			if rankPresCheckerInt[trackerInt] == 0: #If this line rank hasn't been already accounted for
				indexList = [pos1 for pos1 in range(len(self.lineDecIndex)) if self.lineDecIndex[pos1] == interLineRank] # Get all the indices of interNodeRank in the angleDecIndex vector
				for indexListIterator in indexList: # while all the different positions of this rank hasn't been accounted for
					rankPresCheckerInt[indexListIterator] = 1 #Mark the element in the corresponding position of checker vector 1, to indicate this rank has been accounted for 
					interLineTermVec.append(-(self.lineInterDecision[indexListIterator])*(LagMultPi[(*lineDecIndexIterator-1)])) #Calculate cost term
				smallest_element = min(interLineTermVec)
				revisedUpperBound += smallest_element
				interLineTermVec = []
				indexList = []
			trackerInt += 1 #Increment the tracker
		return revisedUpperBound
	#Function getGlobalUpper ends

	def getGlobalLower(regionalLower, numberOfZones): #Function getGlobalLower() returns the global lower bound for the investment coordination problem
		revisedLowerBound = 0 #total lower bound initialized
		for i in range(numberOfZones):
			revisedLowerBound += regionalLower[i]
		return revisedLowerBound
	#Function getGlobalLower ends

	def getGlobalConsensus(self): #Function getGlobalConsensus() returns the global consensus for the investment coordination problem
		globConsensus = 0 #total consensus
		#Calculate the MO upper bound, which is the minimum value of the MO objective, calculated with a given mix of the intermediate, shared continuous variable values from different regions/zones
		for scenPos in range(self.countOfScenarios):
			rankPresChecker = [] #Vector to check if a particular rank has been accounted for 
			length = len(self.angleDecIndex[scenPos]) #length of the rank vector
			for i in range(length): #Put as many zeroes on the rankPresChecker vector as is the length of the rank vector 
				rankPresChecker.append(0)
			tracker = 0 #tracker to track the position of the angleDecIndexIterator
			for angleDecIndexIterator in self.angleDecIndex[scenPos]: #Iterate through rank vector
				interNodeRank = angleDecIndexIterator #Store the value of the rank of the present node iterate in the list 
				if rankPresChecker[tracker] == 0: #If this node rank hasn't been already accounted for
					indexList = [pos1 for pos1 in range(len(self.angleDecIndex[scenPos])) if self.angleDecIndex[scenPos][pos1] == interNodeRank] # Get all the indices of interNodeRank in the angleDecIndex vector
					for indexListIterator in indexList: # while all the different positions of this rank hasn't been accounted for
						rankPresChecker[indexListIterator] = 1 #Mark the element in the corresponding position of checker vector 1, to indicate this rank has been accounted for 
						globConsensus += (self.phaseAngleDecision[scenPos][indexListIterator]-self.interimContDecVar[angleDecIndexIterator-1]) ** 2 #Calculate cost term
					indexList = []
				tracker+=1 #Increment the tracker
		#Calculate the MO upper bound, which is the minimum value of the MO objective, calculated with a given mix of the intermediate, shared discrete variable values from different regions/zones
		rankPresCheckerInt = [] #Vector to check if a particular rank has been accounted for 
		lengthInt = len(self.lineDecIndex) #length of the rank vector
		for i in range(lengthInt): #Put as many zeroes on the rankPresChecker vector as is the length of the rank vector 
			rankPresCheckerInt.append(0)
		trackerInt = 0 #tracker to track the position of the lineDecIndexIterator
		for lineDecIndexIterator in self.lineDecIndex: #Iterate through rank vector
			interLineRank = lineDecIndexIterator #Store the value of the rank of the present line iterate in the list 
			if rankPresCheckerInt[trackerInt] == 0: #If this line rank hasn't been already accounted for
				indexList = [pos1 for pos1 in range(len(self.lineDecIndex)) if self.lineDecIndex[pos1] == interLineRank] # Get all the indices of interNodeRank in the angleDecIndex vector
				for indexListIterator in indexList: # while all the different positions of this rank hasn't been accounted for
					rankPresCheckerInt[indexListIterator] = 1 #Mark the element in the corresponding position of checker vector 1, to indicate this rank has been accounted for 
					globConsensus += (self.lineInterDecision[indexListIterator]-self.interimIntDecVar[lineDecIndexIterator-1]) ** 2 #Calculate cost term
				indexList = []
			trackerInt += 1 #Increment the tracker
		return globConsensus**0.5
	
	#Function getGlobalConsensus ends

	def finDecLineConstr(self): #Final decisions on the construction stauses of candidate lines taking into account the decisions of different zones
		rankPresCheckerInt = [] #Vector to check if a particular rank has been accounted for 
		lengthInt = len(self.lineDecIndex) #length of the rank vector
		for i in range(lengthInt): #Put as many zeroes on the rankPresChecker vector as is the length of the rank vector 
			rankPresCheckerInt.append(0)
		trackerInt = 0 #tracker to track the position of the lineDecIndexIterator
		for lineDecIndexIterator in self.lineDecIndex: #Iterate through rank vector
			interLineRank = lineDecIndexIterator #Store the value of the rank of the present line iterate in the list 
			if rankPresCheckerInt[trackerInt] == 0: #If this line rank hasn't been already accounted for
				first = 0 #Initialize a flag to indicate the first occurence of a rank
				firstInSeries = 0 #flag to indicate if the first verdict for this rank is 1 or 0
				indexList = [pos1 for pos1 in range(len(self.lineDecIndex)) if self.lineDecIndex[pos1] == interLineRank] # Get all the indices of interNodeRank in the angleDecIndex vector
				for indexListIterator in indexList: # while all the different positions of this rank hasn't been accounted for
					rankPresCheckerInt[indexListIterator] = 1 #Mark the element in the corresponding position of checker vector 1, to indicate this rank has been accounted for 
					if self.lineInterDecision[indexListIterator]==1 and first == 0: #Check if the decision of the different subnets is 1
						self.constructedRanks.append(interLineRank)
						firstInSeries = 1 #Set the flag to indicate that the first verdict for this rank is 1
					elif self.lineInterDecision[indexListIterator]==0 and first != 0 and firstInSeries == 1:
						self.constructedRanks.pop()
					first += 1 #Increment the flag to indicate the recurrence of the rank
				indexList = []
			trackerInt += 1 #Increment the tracker
		"""
		for construcIterator in self.constructedRanks:	
			log.info(" First message ")
			log.info("Constructed shared candidate lines are {}".format(construcIterator))
			log.info(" Second message ")
		log.info("end of finDecLineConstr ")
		"""

	def scanBuiltLinesList(self, rankGlob): #Scans the list of built candidate lines for the input argument of global rank to check if this line is built
		pos = self.constructedRanks.index(rankGlob) #find the first position of this rank in the vector
		if pos != self.constructedRanks[len(self.constructedRanks)-1]: #if the rank is found in the list
    			return 1 #return 1
		else:
			return 0 #Otherwise, return 0

	def populateLineDec(self, interimDecBin, zonalIndex, globalRank):#Method to pass the intermediate message for the zonal line building decision to the MO
		self.lineInterDecision.append(interimDecBin)
		self.lineDecIndex.append(globalRank)
		self.zonalIndVectorInt.append(zonalIndex)
		self.lineInterDecVec[zonalIndex*self.candLineNumber+globalRank]=interimDecBin

	def populateAngleDec(self, interimDecAngle, zonalIndex, scenario, globalRank): #Method to pass the intermediate message for the zonal node angle decision to the MO
		(self.phaseAngleDecision[scenario]).append(interimDecAngle)
		(self.angleDecIndex[scenario]).append(globalRank)
		(self.zonalIndVectorCont[scenario]).append(zonalIndex)
		self.phAngDecVec[scenario*self.zoneNumber*self.nodeNumber+zonalIndex*self.nodeNumber+globalRank]=interimDecAngle

	def rewardPenaltyCont(self, lagrangeMult, matrixIndex, iteration): #Function getGlobalUpper() returns the global upper bound for the investment coordination problem
		#cout << "\nITERATION COUNT : " << iteration << endl;
		for scenPos in range(self.countOfScenarios):
			tracker = 0
			for angleDecIndexIterator in self.angleDecIndex[scenPos]:
				if matrixIndex==(scenPos*self.zoneNumber*self.nodeNumber+(self.zonalIndVectorCont[scenPos][tracker]-1)*self.nodeNumber+(angleDecIndexIterator)):		
					if self.compareBasis == 1:
						lagrangeMult += (self.phaseAngleDecision[scenPos][tracker]-self.interimContDecVar[scenPos*self.nodeNumber+angleDecIndexIterator-1])
						#cout << "From MO Continuous Lagrange Multiplier update: " << "Angle from Zones. Tracker: " << tracker << " Scenario: " << scenPos << " Zonal Angle value " << ((phaseAngleDecision[scenPos]).at(tracker)) << " MO Angle index: " << scenPos*nodeNumber+(*angleDecIndexIterator-1) << " MO Angle value: " << (interimContDecVar[scenPos*nodeNumber+(*angleDecIndexIterator-1)]) << " Updated Lagrange Multiplier: " << lagrangeMult << endl;
					else:
						lagrangeMult += (self.phaseAngleDecision[scenPos][tracker]-self.interimContDecVarPrev[scenPos*self.nodeNumber+angleDecIndexIterator-1])
						#log.info("From MO Continuous Lagrange Multiplier update: Angle from Zones. Tracker: {} Scenario: {} Zonal Angle value {} MO Angle index: {} MO Angle value: {} Updated Lagrange Multiplier: {}".format(tracker, scenPos, self.phaseAngleDecision[scenPos][tracker], scenPos*self.nodeNumber+(angleDecIndexIterator-1), interimContDecVar[scenPos*nodeNumber+(angleDecIndexIterator-1)], lagrangeMult))
				tracker += 1
			#log.info("Tracker {}".format(tracker))
		return lagrangeMult
	#Function getGlobalUpper ends

	def rewardPenaltyInteger(self, lagrangeMult, matrixIndex, iteration): #Function getGlobalUpper() returns the global upper bound for the investment coordination problem
		#log.info("\nITERATION COUNT : {}".format(iteration))
		tracker = 0
		for lineDecIndexIterator in self.lineDecIndex:
			if matrixIndex==(self.zonalIndVectorInt[tracker]-1)*self.candLineNumber+lineDecIndexIterator:
				if self.compareBasis == 1:	
					lagrangeMult += (self.lineInterDecision[tracker]-self.interimIntDecVar[lineDecIndexIterator-1])
					#cout << "From MO Integer Lagrange Multiplier update: " << "Line Decision from Zones. Tracker: " << tracker << " Zonal Decision value " << lineInterDecision[tracker] << " MO decision index: " << (*lineDecIndexIterator-1) << " MO Decision value: " << interimIntDecVar[(*lineDecIndexIterator-1)] << " Updated Lagrange Multiplier: " << lagrangeMult << endl
				else:
					lagrangeMult += (self.lineInterDecision[tracker]-self.interimIntDecVarPrev[lineDecIndexIterator-1])
					#log.info("From MO Integer Lagrange Multiplier update: " << "Line Decision from Zones. Tracker: " << tracker << " Zonal Decision value " << lineInterDecision[tracker] << " MO decision index: " << (*lineDecIndexIterator-1) << " MO Decision value: " << interimIntDecVar[(*lineDecIndexIterator-1)] << " Updated Lagrange Multiplier: " << lagrangeMult << endl
			tracker+=1
		#log.info("Tracker {}".format(tracker))
		return lagrangeMult	
	#Function getGlobalUpper ends

	def clearVectors(self): #Clears the different interim vectors for making them ready for the next iteration
		for scenPos in range(self.countOfScenarios):
			#Clear vectors for interim continuous decision variables
			self.phaseAngleDecision[scenPos] = []
			self.angleDecIndex[scenPos] = []
			self.zonalIndVectorCont[scenPos] = []
		self.interimContDecVar = []
		self.interimContDecVarPrev = []
		#Clear vectors for interim integer decision variables
		self.lineInterDecision = []
		self.lineDecIndex = []
		self.zonalIndVectorInt = []
		self.interimIntDecVar = []
		self.interimIntDecVarPrev = []

	def clearDelayedVectors(self): #Clears the different interim vectors only buffer vectors
		self.interimContDecVarPrev = []
		self.interimIntDecVarPrev = []

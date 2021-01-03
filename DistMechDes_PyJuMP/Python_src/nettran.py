#Definition for Nettran class
import julia
import os
import subprocess
import pandas as pd
import numpy as np
import json
import sys
import traceback
from Python_src.log import log
from Python_src.profiler import Profiler
import gurobipy as gp
from gurobipy import GRB
from Python_src.powergenerator import Powergenerator
from Python_src.transl import transmissionLine
from Python_src.load import Load
from Python_src.node import Node
from Python_src.sharedLine import SELine
from Python_src.candidateLine import candLine
from Python_src.intcandidateLine import intCandLine

profiler = Profiler()

if sys.platform in ["darwin", "linux"]:
	log.info("Using Julia executable in {}".format(str(subprocess.check_output(["which", "julia"]), 'utf-8').strip('\n')))
elif sys.platform in ["win32", "win64", "cygwin"]:
	log.info("Using Julia executable in {}".format(str(subprocess.check_output(["which", "julia"]), 'utf-8').strip('\n')))

log.info("Loading Julia...")
profiler.start()
julSol = julia.Julia()
julSol.using("Pkg")
julSol.eval('Pkg.activate(".")')
julSol.include(os.path.join("JuMP_src", "HorMILPDistMech.jl")) # definition of Gensolver class for base case scenario first interval
log.info(("Julia took {:..2f} seconds to start and include Horizontal Investment Coordination mechanism design models.".format(profiler.get_interval())))

class Nettran(object):
	def __init__(self, jSONIndex, zoneIndex, zoneCount, objChoice, milpAlgoChoice): #constructor
		self.otherNodeCount = 0
		self.realizedCLines = 0
		self.realizedIntCLines = 0
		self.zonalCount = zoneCount
		self.zonalIndex = zoneIndex #Assigns the ID number of this zone
		self.netFile = jSONIndex['Network File'] #String for storing the name of the network file
		self.genFile = jSONIndex['Generator File'] #String for storing the name of the generator file
		self.sharedLineFile = jSONIndex['Shared Lines File'] # String for storing the name of the shared existing lines file
		self.tranFile = jSONIndex['Transmission Lines File'] #String for storing the name of the transmission line file
		self.loadFile = jSONIndex['Load File'] #String for storing the name of the load file
		self.candLineFile = jSONIndex['Candidate Lines File'] #String for storing the name of the candidate lines file
		self.intCandLineFile = jSONIndex['Intra Candidate Lines File'] #String for storing the name of the candidate lines file
		self.containsFlag = 0 #Indicates whether a particular node in a particular connected zone has been accounted for or not; 0 for no, 1 for yes
		#Specify the type of the curve
		self.nodeObject = []
		self.genObject = []
		self.simMode = objChoice
		self.diffZoneNodeID = []; self.diffZoneNodeID.append(0) #Initialize the list of external-zone connected node ID's so that it's not empty
		self.diffZoneID = []; self.diffZoneID.append(0) #Initialize the list of external connected zone ID's so that it's not empty
		self.globalRankDiffNode = []; self.globalRankDiffNode.append(0) #Initialize the list of external-zone connected node global rank so that it's not empty
		self.diffZoneNodeExistingID = []; self.diffZoneNodeExistingID.append(0) #initialize the list of external-zone existing node ID's so that it's not empty
		self.diffZoneExistingID = []; self.diffZoneExistingID.append(0) #Initialize the list of external-zone existing zone ID's so that it's not empty.
		self.globalExistingRank = []; self.globalExistingRank.append(0) #Initialize the list os external-zone connected node global rank so that it's not empty
		self.lpSolveAlgo = milpAlgoChoice #Simplex for 1 and IPM for 2
		#/* Nodes */################################################################################################################################################################################
		matrixNetFile = json.load(open(os.path.join("data", self.netFile)))
		self.nodeNumber = matrixNetFile['nodeNumber'] 
		self.sharedELines = matrixNetFile['sharedELines'] 
		self.sharedCLines = matrixNetFile['sharedCLines'] 
		self.genNumber = matrixNetFile['genNumber'] 
		self.loadNumber = matrixNetFile['loadNumber'] 
		self.tranNumber = matrixNetFile['tranNumber']
		self.internalCLines = matrixNetFile['internalCLines'] #get the dimensions of the Network
		self.countOfScenarios = matrixNetFile['loadStochScenarios'] #get the number of stochastic scenarios of load
		for l in range(self.nodeNumber):
			#log.info("\nCreating the {} -th Node:".format(l + 1))
		
			nodeInstance = Node(l + 1, self.zonalIndex) #creates nodeInstance object with ID l + 1

			self.nodeObject.append( nodeInstance ) #pushes the nodeInstance object into the vector
		#end initialization for Nodes
		matrixNetFile.close() #close the network file
		#/* Generators */##########################################################################################################################################################################
		#/* Instantiate Generators */
		# Open the .json file to read the Powergenerator parameter values
		matrixGen = json.load(open(os.path.join("data", self.genFile))) #ifstream constructor opens the file of Generators
		j = 0 #counter for generators
		for matrixGenFile in matrixGen:
			gNodeID = matrixGenFile['genNodeID'] #node object ID to which the particular generator object is connected
			#log.info("\nConnection Node defined.\n")
			#Parameters for Generator
			#Quadratic Coefficient: 
			c2 = matrixGenFile['quadCostCoeff'] * (100**2)
			#Linear coefficient: 
			c1 = matrixGenFile['linCostCoeff'] * 100
			#Constant term: 
			c0 = matrixGenFile['noLoadCost']
			#Maximum Limit: 
			PgMax = matrixGenFile['PgMax'] / 100
			#Minimum Limit: 
			PgMin = matrixGenFile['PgMin'] / 100
			#/* Secant Approximation of the Quadratic Cost curve */
			#Tangent Ratio of the secant approximation of the intercepted cost curve
			tanTheta = (c2*(PgMax**2)+c1*PgMax-c2*(PgMin**2)-c1*PgMin)/(PgMax-PgMin)
			minCost = c2*(PgMin**2)+c1*PgMin+c0-tanTheta*PgMin
			#Intercept value or cost at minimum power level
			#check the bounds and validity of the parameter values
			genInstance = Powergenerator(j+1, self.nodeObject[ gNodeID - 1 ],  tanTheta, minCost, PgMax, PgMin)
			self.genObject.append(genInstance) #push the generator object into the array
			j +=1 #increment counter
		matrixGen.close() #Close the generator file
		self.genDF = pd.DataFrame([g.__dict__ for g in self.genObject])
		#/* Transmission Lines */###################################################################################################################################################################
		matrixTranFile = json.load(open(os.path.join("data", self.tranFile))) #ifstream constructor opens the file of Transmission lines
		#/* Instantiate Transmission Lines */
		k=0 #counter for transmission lines
		for matrixTran in matrixTranFile:
			#node IDs of the node objects to which this transmission line is connected.
			tNodeID1 = matrixTran['fromNode'] #From end
			tNodeID2 = matrixTran['toNode'] #To end
			#Reactance
			reacT = matrixTran['Reactance']
			#values of maximum allowable power flow on line in the forward and reverse direction:
			ptMax = matrixTran['lineLimit']/100
			#creates transLineInstance object with ID k + 1
			transLineInstance = transmissionLine(k + 1, self.nodeObject[ tNodeID1 - 1 ], self.nodeObject[ tNodeID2 - 1 ], ptMax, reacT) 
			self.translObject.append( transLineInstance ) #pushes the transLineInstance object into the vector
			k +=1 #increment the counter
		#end initialization for Transmission Lines 
		matrixTranFile.close() #Close the transmission line file 
		self.tranDF = pd.DataFrame([t.__dict__ for t in self.translObject])
		#/* Shared Existing Transmission Lines */####################################################################################################################################################
		matrixSETranFile = json.load(open(os.path.join("data", self.sharedLineFile))) #ifstream constructor opens the file of Transmission lines
		
		#/* Instantiate Shared Existing Transmission Lines */
		for matrixSETran in matrixSETranFile:
			#log.info("Tran File Test Message 1 from line {} before creation1".format(k))
			#node IDs of the node objects to which this transmission line is connected.
			serNum = matrixSETran['globalSerial'] #global serial number of the shared existing transmission line
			tNodeID1 = matrixSETran['fromNode'] #From end node 
			nodeZone1 = matrixSETran['fromZone'] #From end zone number
			tNodeID2 = matrixSETran['toNode'] #To end node
			nodeZone2 = matrixSETran['toZone'] #To end zone number
			#Reactance:
			reacT = matrixSETran['Reactance']
			#values of maximum allowable power flow on line in the forward and reverse direction:
			ptMax = matrixSETran['lineLimit']/100
			#creates SELine object with ID k + 1
			if nodeZone1 == self.zonalIndex: #If the node 1 belongs to this zone
				SELineInstance = SELine(k + 1, serNum, self.nodeObject[ tNodeID1 - 1 ], tNodeID1, nodeZone1, tNodeID2, nodeZone2, self.zonalIndex, ptMax, reacT) #Create the shared existing transmission line object with node 1 
				indCount = 0 #Initialize a counter for tracking the position in the vector of the iterator
				toFromFlag = 1 #Indicates that the from node is the intra-zonal node 
				for diffZNIt in self.diffZoneNodeID:
					if diffZNIt == tNodeID2 and self.diffZoneID[indCount] == nodeZone2: #Check whether the other-zone node is already present in the list
						self.containsFlag = 1 #Loop through the diffZoneNodeID and diffZoneID vectors, to see if the other zone-node is already present or not, if yes, set the containsFlag to 1
						SELineInstance.outerNodeIndex(indCount, toFromFlag) #Calling the SELine method, which in turn calls the intra-zonal node method, to create a list of rankings of all the outer-zone nodes that are connected to this node, if the intra-zonal node is already present
					indCount += 1 #Increment the counter
				if self.containsFlag == 0: #If the diffZoneNodeID and diffZoneID vectors do not contain the other zone-node, then push the node and the other zone index in the respective vectors
					self.diffZoneNodeID.append(tNodeID2)
					self.diffZoneID.append(nodeZone2)
					self.diffZoneNodeExistingID.append(tNodeID2) #initialize the list of external-zone existing node ID's so that it's not empty
					self.diffZoneExistingID.append(nodeZone2) #Initialize the list of external-zone existing zone ID's so that it's not empty.
					self.otherNodeCount += 1 #Increment the counter to account for the total number of other-zone nodes
					SELineInstance.outerNodeIndex(self.otherNodeCount, toFromFlag) #Calling the SELine method, which in turn calls the intra-zonal node method, to create a list of rankings of all the outer-zone nodes that are connected to this node, if the intra-zonal node is newly added
				self.SELineObject.append(SELineInstance) #pushes the transLineInstance object into the vector
			else: #Otherwise, if the node 2 belongs to this zone
				SELineInstance = SELine(k + 1, serNum, self.nodeObject[ tNodeID2 - 1 ], tNodeID1, nodeZone1, tNodeID2, nodeZone2, self.zonalIndex, ptMax, reacT) #Create the shared existing transmission line object with node 2
				indCount = 0 #Initialize a counter for tracking the position in the vector of the iterator
				toFromFlag = -1 #Indicates that the To node is the intra-zonal node
				for diffZNIt in self.diffZoneNodeID:
					if diffZNIt == tNodeID1 and self.diffZoneID[indCount] == nodeZone1: #Check whether the other-zone node is already present in the list
						self.containsFlag = 1 #Loop through the diffZoneNodeID and diffZoneID vectors, to see if the other zone-node is already present or not, if yes, set the containsFlag to 1
						SELineInstance.outerNodeIndex(indCount, toFromFlag) #Calling the SELine method, which in turn calls the intra-zonal node method, to create a list of rankings of all the outer-zone nodes that are connected to this node, if the intra-zonal node is already present
					indCount += 1 #Increment the counter
				if self.containsFlag == 0: #If the diffZoneNodeID and diffZoneID vectors do not contain the other zone-node, then push the node and the other zone index in the respective vectors
					self.diffZoneNodeID.append(tNodeID1)
					self.diffZoneID.append(nodeZone1)
					self.diffZoneNodeExistingID.append(tNodeID1) #initialize the list of external-zone existing node ID's so that it's not empty
					self.diffZoneExistingID.append(nodeZone1) #Initialize the list of external-zone existing zone ID's so that it's not empty.
					self.otherNodeCount += 1  #Increment the counter to account for the total number of other-zone nodes
					SELineInstance.outerNodeIndex(self.otherNodeCount, toFromFlag) #Calling the SELine method, which in turn calls the intra-zonal node method, to create a list of rankings of all the outer-zone nodes that are connected to this node, if the intra-zonal node is newly added
				self.SELineObject.append( SELineInstance ) # pushes the transLineInstance object into the vector	
			self.containsFlag = 0 #Reset the containsFlag for matching the next item
		#end initialization for Shared Existing Transmission Lines
		matrixSETranFile.close() #Close the shared existing file
		self.sharedExistingDF = pd.DataFrame([seDF.__dict__ for seDF in self.SELineObject])
		#/* Shared Candidate Transmission Lines */######################################################################################################################################################
		matrixCETranFile = json.load(open(os.path.join("data", self.candLineFile))) #ifstream constructor opens the file of candidate Transmission lines

		#/* Instantiate Shared Candidate Transmission Lines */
		for matrixCETran in matrixCETranFile:
			# node object IDs to which the particular transmission line object is connected
			# node IDs of the node objects to which this transmission line is connected.
			serNum = matrixCETran['globalSerial'] #global serial number of the shared existing transmission line
			tNodeID1 = matrixCETran['fromNode'] #From end node 
			nodeZone1 = matrixCETran['fromZone'] #From end zone number
			tNodeID2 = matrixCETran['toNode'] #To end node
			nodeZone2 = matrixCETran['toZone'] #To end zone number
			#Parameters for Transmission Line
			#Reactance
			reacT = matrixCETran['Reactance']
			#values of maximum allowable power flow on line in the forward and reverse direction:
			#Forward direction:
			ptMax = matrixCETran['lineLimit']/100
			lifeTime = matrixCETran['lifeTime'] #life time of the candidate line
			interestRate = matrixCETran['interestRate'] #interest rate of the investment 
			costPerCap = matrixCETran['costPerCap']*ptMax #capital cost for the construction 
			presAbsence = matrixCETran['presAbsence'] #status of the construction 
			ownership = matrixCETran['ownership'] #ownership of the candidate line for this zone
			
			#creates candLineInstance object with ID k + 1
			if nodeZone1 == self.zonalIndex:
				candLineInstance = candLine( k + 1, serNum, self.nodeObject[ tNodeID1 - 1 ], tNodeID1, nodeZone1, tNodeID2, nodeZone2, self.zonalIndex, ptMax, reacT, interestRate, lifeTime, costPerCap, presAbsence, ownership ) #Create the shared candidate transmission line object with node 1 
				indCount = 0 #Initialize a counter for tracking the position in the vector of the iterator
				toFromFlag = 1 #Indicates that the from node is the intra-zonal node 
				for diffZNIt in self.diffZoneNodeID:
					if diffZNIt == tNodeID2 and self.diffZoneID[indCount] == nodeZone2:
						containsFlag = 1 #Loop through the diffZoneNodeID and diffZoneID vectors, to see if the other zone-node is already present or not, if yes, set the containsFlag to 1
						candLineInstance.outerNodeIndex(indCount, toFromFlag) #Calling the SELine method, which in turn calls the intra-zonal node method, to create a list of rankings of all the outer-zone nodes that are connected to this node, if the intra-zonal node is already present
					indCount +=1 #Increment the counter
				if self.containsFlag == 0: #If the diffZoneNodeID and diffZoneID vectors do not contain the other zone-node, then push the node and the other zone index in the respective vectors
					self.diffZoneNodeID.append(tNodeID2)
					self.diffZoneID.append(nodeZone2)
					self.otherNodeCount += 1 #Increment the counter to account for the total number of other-zone nodes
					candLineInstance.outerNodeIndex(self.otherNodeCount, toFromFlag) #Calling the SELine method, which in turn calls the intra-zonal node method, to create a list of rankings of all the outer-zone nodes that are connected to this node, if the intra-zonal node is newly added
				self.candLineObject.append( candLineInstance ) #pushes the transLineInstance object into the vector
			else:
				candLineInstance = candLine( k + 1, serNum, self.nodeObject[ tNodeID2 - 1 ], tNodeID1, nodeZone1, tNodeID2, nodeZone2, self.zonalIndex, ptMax, reacT, interestRate, lifeTime, costPerCap, presAbsence, ownership ) #Create the shared candidate transmission line object with node 2
				indCount = 0 #Initialize a counter for tracking the position in the vector of the iterator
				toFromFlag = -1 #Indicates that the To node is the intra-zonal node 
				for diffZNIt in self.diffZoneNodeID:
					if diffZNIt == tNodeID1 and self.diffZoneID[indCount] == nodeZone1:
						containsFlag = 1 #Loop through the diffZoneNodeID and diffZoneID vectors, to see if the other zone-node is already present or not, if yes, set the containsFlag to 1
						candLineInstance.outerNodeIndex(indCount, toFromFlag) #Calling the SELine method, which in turn calls the intra-zonal node method, to create a list of rankings of all the outer-zone nodes that are connected to this node, if the intra-zonal node is already present
					indCount += 1 #Increment the counter
				if self.containsFlag == 0: #If the diffZoneNodeID and diffZoneID vectors do not contain the other zone-node, then push the node and the other zone index in the respective vectors
					self.diffZoneNodeID.append(tNodeID1)
					self.diffZoneID.append(nodeZone1)
					self.otherNodeCount+=1 #Increment the counter to account for the total number of other-zone nodes
					candLineInstance.outerNodeIndex(self.otherNodeCount, toFromFlag) #Calling the SELine method, which in turn calls the intra-zonal node method, to create a list of rankings of all the outer-zone nodes that are connected to this node, if the intra-zonal node is newly added
				self.candLineObject.append( candLineInstance ) #pushes the transLineInstance object into the vector
			containsFlag = 0 #Reset the containsFlag for matching the next item
		#end initialization for candidate Transmission Lines
		matrixCETranFile.close() #Close the candidate lines file
		if self.internalCLines>0:
		#/* Internal Candidate Transmission Lines */################################################################################################################################################
			matrixIntCETranFile = json.load(open(os.path.join("data", self.intCandLineFile))) #ifstream constructor opens the file of internal candidate Transmission lines
			#/* Instantiate Internal Candidate Transmission Lines */
			for matrixIntCETran in matrixIntCETranFile:
				#node object IDs to which the particular transmission line object is connected
				#node IDs of the node objects to which this transmission line is connected.
				tNodeID1 = matrixIntCETran['fromNode'] #From end node 
				tNodeID2 = matrixIntCETran['toNode'] #To end node
				#Parameters for Transmission Line
				#Reactance:
				reacT = matrixIntCETran['Reactance']
				#values of maximum allowable power flow on line in the forward and reverse direction:
				#Forward direction:
				ptMax = matrixIntCETran['lineLimit']/100
				lifeTime = matrixIntCETran['lifeTime'] #life time of the candidate line
				interestRate = matrixIntCETran['interestRate'] #interest rate of the investment 
				costPerCap = matrixIntCETran['costPerCap']*ptMax #capital cost for the construction 
				presAbsence = matrixIntCETran['presAbsence'] #status of the construction
			
				#creates intCandLineInstance object with ID k + 1
				intCandLineInstance = intCandLine( k + 1, self.nodeObject[ tNodeID1 - 1 ], self.nodeObject[ tNodeID2 - 1 ], ptMax, reacT, interestRate, lifeTime, costPerCap, presAbsence ) #Create the internal candidate transmission line object with node 1 
				self.intCandLineObject.append( intCandLineInstance ) #pushes the transLineInstance object into the vector
			#end initialization for candidate Transmission Lines
			matrixIntCETranFile.close() #Close the candidate lines file
		#/* Loads */################################################################################################################################################################################
		matrixLoadFile = json.load(open(os.path.join("data", self.loadFile))) #ifstream constructor opens the file of Loads

		int loadFields; // Number of columns in the load file
		matrixLoadFile >> loadFields; // get the dimensions of the Load matrix
		countOfScenarios = loadFields-1
		#Initialize the default loads on all nodes to zero
		for l in range(self.nodeNumber):
			(self.nodeObject[l]).initLoad(self.countOfScenarios) #Initialize the default loads on all nodes to zero
		#end initialization for Nodes		
		#/* Instantiate Loads */
		for matrixLoad in matrixLoadFile:
			#node object ID to which the particular load object is connected
			#node ID of the node object to which this load object is connected.
			lNodeID = matrixLoad['loadNodeID']
			for f in range(self.countOfScenarios):
				#value of allowable power consumption capability of load in pu with a negative sign to indicate consumption:
				#Power Consumption:
				P_Load[f] = matrixLoad['loadMW']/100

			loadInstance = Load( j + 1, self.nodeObject[ lNodeID - 1 ], loadFields-1, P_Load ) #creates loadInstance object object with ID number j + 1

			self.loadObject.append( loadInstance ) #pushes the loadInstance object into the vector
		#end initialization for Loads
		matrixLoadFile.close() #Closes the load file
		for f in range(self.countOfScenarios):
			self.probability.append(1/(self.countOfScenarios))

		#check the bounds and validity of the parameter values
		otherZoneNodeIter = self.diffZoneNodeID.begin() #Initialize the otherZoneNodeIter iterator to point to the beginning of the diffZoneNodeID vector
		otherZoneNodeIter+=1 #Increment the pointer to point to the actual non-zero entry
		otherZoneIter = self.diffZoneID.begin() #Initialize the otherZoneIter iterator to point to the beginning of the diffZoneID vector
		otherZoneIter+=1 #Increment the pointer to point to the actual non-zero entry	
		sharedELineIt = self.SELineObject.begin() #Initialize the sharedELineIt iterator to point to the beginning of the SELineObject vector 
		sharedELineFromIt = self.SELineObject.begin() #Initialize the sharedELineFromIt iterator to point to the beginning of the SELineObject vector
		sharedELineFZoneIt = self.SELineObject.begin() #Initialize the sharedELineFZoneIt iterator to point to the beginning of the SELineObject vector
		sharedELineToIt = self.SELineObject.begin() #Initialize the sharedELineToIt iterator to point to the beginning of the SELineObject vector
		sharedELineTZoneIt = self.SELineObject.begin() #Initialize the sharedELineTZoneIt iterator to point to the beginning of the SELineObject vector
		sharedELineReactIt = self.SELineObject.begin() #Initialize the sharedELineReactIt iterator to point to the beginning of the SELineObject vector
		sharedELineCapIt = self.SELineObject.begin() #Initialize the sharedELineCapIt iterator to point to the beginning of the SELineObject vector
		sharedCandLineIt = self.candLineObject.begin() #Initialize the sharedCandLineIt iterator to point to the beginning of the candLineObject vector 
		sharedCandLineFromIt = self.candLineObject.begin() #Initialize the sharedCandLineFromIt iterator to point to the beginning of the candLineObject vector
		sharedCandLineFZoneIt = self.candLineObject.begin() #Initialize the sharedCandLineFZoneIt iterator to point to the beginning of the candLineObject vector
		sharedCandLineToIt = self.candLineObject.begin() #Initialize the sharedCandLineToIt iterator to point to the beginning of the candLineObject vector
		sharedCandLineTZoneIt = self.candLineObject.begin() #Initialize the sharedCandLineTZoneIt iterator to point to the beginning of the candLineObject vector
		sharedCandLineReactIt = self.candLineObject.begin() #Initialize the sharedCandLineReactIt iterator to point to the beginning of the candLineObject vector
		sharedCandLineCapIt = self.candLineObject.begin() #Initialize the sharedCandLineCapIt iterator to point to the beginning of the candLineObject vector
	#end constructor

	def getNumberOfScenarios(self): #Returns the number of scenarios
		return self.countOfScenarios

	def MILPPiecewiseLin(): #Function MILPPiecewiseLin() implements the Mixed Integer Linear Programming Unit Commitment Solver routine by calling GLPK routines for piecewise linear objective
		log.info("\nUnder Construction")

	def MILPPolynomial(): #Function MILPPolynomial() implements the Mixed Integer Linear Programming Unit Commitment Solver routine by calling GLPK routines for polynomial convex objective
		log.info("\nUnder Construction")

	def MILPAvgHR(self, coordInstanceRef, LagMultXi, LagMultPi, totalCandLineNum, totalSharedNodeNum): #Function MILPAvgHR() implements the Mixed Integer Linear Programming Unit Commitment Solver routine by calling GLPK routines for average heat rate objective for Horizontal Coordination Investment decision making
		zonal_dec_stageI = julSol.HorMILPDistMech(coordInstanceRef, LagMultXi, LagMultPi, totalCandLineNum, totalSharedNodeNum)
	#/* CREATION OF THE MIP SOLVER INSTANCE */
	
	#Function MILP() ends

	def calcMILPBounds(self, LagMultXi, LagMultPi, totalCandLineNum, totalSharedNodeNum): #Function MILPAvgHR() implements the Mixed Integer Linear Programming Unit Commitment Solver routine by calling GLPK routines for average heat rate objective for Horizontal Coordination Investment decision making
		zonal_dec_stageI = julSol.HorMILPDistMech(LagMultXi, LagMultPi, totalCandLineNum, totalSharedNodeNum)


	def getConnZone(self, i): #returns the pointer to the base of the vector, diffZoneID
		if self.otherZoneIter != self.diffZoneID.end(): #check to see if the end of the diffZoneID vector has been reached
			self.otherZoneIter += 1 #if not, increment otherZoneIter iterator for the next call
			return self.diffZoneID[i] #return the zone ID present at the position i
		else:
			return -1 #if end is reached, return -1
	#Function getConnZone() ends

	def getConnNode(self, i): #returns the pointer to the base of the vector, diffZoneNodeID
		if self.otherZoneNodeIter != self.diffZoneNodeID.end(): #check to see if the end of the diffZoneNodeID vector has been reached
			self.otherZoneNodeIter+=1 #if not, increment otherZoneNodeIter iterator for the next call
			return self.diffZoneNodeID[i] #return the node ID present at the position i
		else:
			return -1 #if end is reached, return -1
	#Function getConnNode() ends

	def getSESerial(self, i): #returns the pointer to the base of the vector, SELineObject
		if self.sharedELineIt != self.SELineObject.end(): #check to see if the end of the SELineObject vector has been reached
			self.sharedELineIt+=1 #if not, increment sharedELineIt iterator for the next call
			return (self.SELineObject[i]).getTranslID() #return the global serial number of the present SE line at the position i
		else:
			return -1 #if end is reached, return -1
	#Function getSESerial() ends

	def getSEFromNode(self, i): #returns the ID number of the from node of the shared existing line
		if self.sharedELineFromIt != self.SELineObject.end(): #check to see if the end of the SELineObject vector has been reached
			self.sharedELineFromIt+=1 #if not, increment iterator for the next call
			if (self.SELineObject[i]).getFlowDir()==1: #If the intrazonal node is the from node
				return (self.SELineObject[i]).getIntlNodeID() #return the ID of the internal node of the SE line at the position i
			else:
				return (self.SELineObject[i]).getExtNodeID() #else, return the ID of the external node of the SE line at the position i
		else:
			return -1 #if end is reached, return -1 
	#Function getSEFromNode() ends

	def getSEFromZone(self, i): #returns the ID number of the from zone of the shared existing line
		if self.sharedELineFZoneIt != self.SELineObject.end(): #check to see if the end of the SELineObject vector has been reached
			self.sharedELineFZoneIt+=1 #if not, increment iterator for the next call
			if (self.SELineObject[i]).getFlowDir()==1: #If the intrazonal node is the from node
				return (self.SELineObject[i]).getIntlZoneID() #return the ID of current zone of the SE line at the position i
			else:
				return (self.SELineObject[i]).getExtZoneID() #else, return the ID of the external zone of the SE line at the position i
		else:
			return -1 #if end is reached, return -1 
	#Function getSEFromZone() ends

	def getSEToNode(self, i): #returns the ID number of the to node of the shared existing line
		if self.sharedELineToIt != self.SELineObject.end(): #check to see if the end of the SELineObject vector has been reached
			self.sharedELineToIt+=1 #if not, increment iterator for the next call
			if (self.SELineObject[i]).getFlowDir()==1: #If the intrazonal node is the from node
				return (self.SELineObject[i]).getExtNodeID() #return the ID of the external node of the SE line at the position i
			else:
				return (self.SELineObject[i]).getIntlNodeID() #else, return the ID of the internal node of the SE line at the position i
		else:
			return -1 #if end is reached, return -1 
	#Function getSEToNode() ends

	def getSEToZone(self, i): #returns the ID number of the to zone of the shared existing line
		if self.sharedELineTZoneIt != self.SELineObject.end(): #check to see if the end of the SELineObject vector has been reached
			self.sharedELineTZoneIt+=1 #if not, increment iterator for the next call
			if (self.SELineObject[i]).getFlowDir()==1: #If the intrazonal node is the from node
				return (self.SELineObject[i]).getExtZoneID() #return the ID of external zone of the SE line at the position i
			else:
				return (self.SELineObject[i]).getIntlZoneID() #else, return the ID of the current zone of the SE line at the position i
		else:
			return -1 #if end is reached, return -1 
	#Function getSEToZone() ends

	def getSEReactance(self, i): #returns the reactance of the shared existing line
		if self.sharedELineReactIt != self.SELineObject.end(): #check to see if the end of the SELineObject vector has been reached
			self.sharedELineReactIt+=1 #if not, increment iterator for the next call
			return (self.SELineObject[i]).getReactance() #return the reactance of the SE line at the position i
		else:
			return -1 #if end is reached, return -1 
	#Function getConnNode() ends

	def getSECapacity(self, i): #returns the capacity of the shared existing line
		if self.sharedELineCapIt != self.SELineObject.end(): #check to see if the end of the SELineObject vector has been reached
			self.sharedELineCapIt+=1 #if not, increment iterator for the next call
			return (self.SELineObject[i]).getFlowLimit() #return the capacity of the SE line at the position i
		else:
			return -1 #if end is reached, return -1 
	 #Function getConnNode() ends

	def getCandSerial(self, i): #returns the pointer to the base of the vector, SELineObject
		if self.sharedCandLineIt != self.candLineObject.end(): #check to see if the end of the SELineObject vector has been reached
			self.sharedCandLineIt+=1 #if not, increment sharedELineIt iterator for the next call
			return (self.candLineObject[i]).getTranslID() #return the global serial number of the present SE line at the position i
		else:
			return -1 #if end is reached, return -1
	#Function getSESerial() ends

	def getCandFromNode(self, i): #returns the ID number of the from node of the shared existing line
		if (sharedCandLineFromIt != candLineObject.end()) { // check to see if the end of the SELineObject vector has been reached
			sharedCandLineFromIt++; // if not, increment iterator for the next call
			if ((candLineObject.at(i))->getFlowDir()==1) // If the intrazonal node is the from node
				return (candLineObject.at(i))->getIntlNodeID(); // return the ID of the internal node of the SE line at the position i
			else
				return (candLineObject.at(i))->getExtNodeID(); // else, return the ID of the external node of the SE line at the position i
		}
		else
			return -1; // if end is reached, return -1 
	#Function getSEFromNode() ends

	def getCandFromZone(self, i) // returns the ID number of the from zone of the shared existing line
{
	if (sharedCandLineFZoneIt != candLineObject.end()) { // check to see if the end of the SELineObject vector has been reached
		sharedCandLineFZoneIt++; // if not, increment iterator for the next call
		if ((candLineObject.at(i))->getFlowDir()==1) // If the intrazonal node is the from node
			return (candLineObject.at(i))->getIntlZoneID(); // return the ID of current zone of the SE line at the position i
		else
			return (candLineObject.at(i))->getExtZoneID(); // else, return the ID of the external zone of the SE line at the position i
	}
	else
		return -1; // if end is reached, return -1 
} // Function getSEFromZone() ends

	def getCandToNode(self, i) // returns the ID number of the to node of the shared existing line
{
	if (sharedCandLineToIt != candLineObject.end()) { // check to see if the end of the SELineObject vector has been reached
		sharedCandLineToIt++; // if not, increment iterator for the next call
		if ((candLineObject.at(i))->getFlowDir()==1) // If the intrazonal node is the from node
			return (candLineObject.at(i))->getExtNodeID(); // return the ID of the external node of the SE line at the position i
		else
			return (candLineObject.at(i))->getIntlNodeID(); // else, return the ID of the internal node of the SE line at the position i
	}
	else
		return -1; // if end is reached, return -1 
} // Function getSEToNode() ends

	def getCandToZone(self, i) // returns the ID number of the to zone of the shared existing line
{
	if (sharedCandLineTZoneIt != candLineObject.end()) { // check to see if the end of the SELineObject vector has been reached
		sharedCandLineTZoneIt++; // if not, increment iterator for the next call
		if ((candLineObject.at(i))->getFlowDir()==1) // If the intrazonal node is the from node
			return (candLineObject.at(i))->getExtZoneID(); // return the ID of external zone of the SE line at the position i
		else
			return (candLineObject.at(i))->getIntlZoneID(); // else, return the ID of the current zone of the SE line at the position i
	}
	else
		return -1; // if end is reached, return -1 
} // Function getSEToZone() ends

	def getCandReactance(self, i) // returns the reactance of the shared existing line
{
	if (sharedCandLineReactIt != candLineObject.end()) { // check to see if the end of the SELineObject vector has been reached
		sharedCandLineReactIt++; // if not, increment iterator for the next call
		return (candLineObject.at(i))->getReactance(); // return the reactance of the SE line at the position i
	}
	else
		return -1; // if end is reached, return -1 
} // Function getConnNode() ends

	def getCandCapacity(self, i) // returns the capacity of the shared existing line
{
	if (sharedCandLineCapIt != candLineObject.end()) { // check to see if the end of the SELineObject vector has been reached
		sharedCandLineCapIt++; // if not, increment iterator for the next call
		return (candLineObject.at(i))->getFlowLimit(); // return the capacity of the SE line at the position i
	}
	else
		return -1; // if end is reached, return -1 
} // Function getConnNode() ends

	def setSEFromRank(self, existingCounter, rankFrom) // sets the from rank for the internal zone node end of the shared existing line
{
	 (*(SELineObject.at(existingCounter))).assignRank(rankFrom);
}

	def setSEToRank(self, existingCounter, rankTo) // sets the to rank for the internal zone node end of the shared existing line
{
	(*(SELineObject.at(existingCounter))).assignRank(rankTo);
}

	def setCandFromRank(self, candidateCounter, rankFrom) // sets the from rank for the internal zone node end of the shared candidate line
{
	(*(candLineObject.at(candidateCounter))).assignRank(rankFrom);
}

	def setCandToRank(self, candidateCounter, rankTo) // sets the to rank for the internal zone node end of the shared candidate line
{
	(*(candLineObject.at(candidateCounter))).assignRank(rankTo);
}

	def setSEFromRankConn(self, existingCounter, rankFrom) // populates the vector of the global ranks of all the external from nodes of SE lines
{
	for (globalIterator = globalRankDiffNode.begin(); globalIterator != globalRankDiffNode.end(); ++globalIterator) {
		if (*globalIterator == rankFrom) { // Check whether the other-zone node is already present in the list
			containsFlag = 1;  // Loop through the diffZoneNodeID and diffZoneID vectors, to see if the other zone-node is already present or not, if yes, set the containsFlag to 1
		}
	}
	if (containsFlag == 0) { // If globalRankDiffNode vector does not contain the other zone-node rank, then push the node rank in the vector
		globalRankDiffNode.push_back(rankFrom);
		globalExistingRank.push_back(rankFrom); // list of global ranking of external zone nodes for only SE lines and shared candidate lines
		++existingOtherZoneNodeNum; // Increment the number of other zone nodes that are ends of SE lines, only when a new one is added
	}
	containsFlag = 0; // reset
	(*(SELineObject.at(existingCounter))).connectRank(rankFrom);			
}

	def setSEToRankConn(self, existingCounter, rankTo) // populates the vector of the global ranks of all the external to nodes of SE lines
{
	for (globalIterator = globalRankDiffNode.begin(); globalIterator != globalRankDiffNode.end(); ++globalIterator) {
		if (*globalIterator == rankTo) { // Check whether the other-zone node is already present in the list
			containsFlag = 1;  // Loop through the diffZoneNodeID and diffZoneID vectors, to see if the other zone-node is already present or not, if yes, set the containsFlag to 1
		}
	}
	if (containsFlag == 0) { // If globalRankDiffNode vector does not contain the other zone-node rank, then push the node rank in the vector
		globalRankDiffNode.push_back(rankTo);
		globalExistingRank.push_back(rankTo); // list of global ranking of external zone nodes for only SE lines and shared candidate lines
		++existingOtherZoneNodeNum; // Increment the number of other zone nodes that are ends of SE lines, only when a new one is added
	}
	containsFlag = 0; // reset
	(*(SELineObject.at(existingCounter))).connectRank(rankTo);
}

	def setCandFromRankConn(self, candidateCounter, rankFrom): #populates the vector of the global ranks of all the external from nodes of cand lines
		for globalIterator in self.globalRankDiffNode:
			if globalIterator == rankFrom: #Check whether the other-zone node is already present in the list
				self.containsFlag = 1 #Loop through the diffZoneNodeID and diffZoneID vectors, to see if the other zone-node is already present or not, if yes, set the containsFlag to 1
		if self.containsFlag == 0: #If globalRankDiffNode vector does not contain the other zone-node rank, then push the node rank in the vector
			self.globalRankDiffNode.append(rankFrom)
		self.containsFlag = 0 #reset
		(self.candLineObject[candidateCounter]).connectRank(rankFrom)

	def setCandToRankConn(self, candidateCounter, rankTo): #populates the vector of the global ranks of all the external to nodes of cand lines
		for globalIterator in self.globalRankDiffNode:
			if globalIterator == rankTo: #Check whether the other-zone node is already present in the list
				self.containsFlag = 1 #Loop through the diffZoneNodeID and diffZoneID vectors, to see if the other zone-node is already present or not, if yes, set the containsFlag to 1
		if self.containsFlag == 0: #If globalRankDiffNode vector does not contain the other zone-node rank, then push the node rank in the vector
			self.globalRankDiffNode.append(rankTo)
		self.containsFlag = 0 #reset
		(self.candLineObject[candidateCounter]).connectRank(rankTo)

	def assignCandGlobalRank(self, candidateCounter, candGlobalRank): #assigns the global rank of the shared candidate line
		(self.candLineObject[candidateCounter]).assignLineRank(candGlobalRank)

	def setRealizedCLines(self, coordInstanceRef): #Assigns the realizedCLine variable, the number of candidate lines that are actually built
		#Built shared candidate lines
		for candIterator in self.candLineObject:
			candGlobRank=candIterator.getGlobalRank() #Get the global rank of the candidate line 		
			statusPresAbs = coordInstanceRef.scanBuiltLinesList(candGlobRank) #Passing on the shared node angle decision message to the MO
			if statusPresAbs==1:
				candIterator.setPresAbsStatus()
			if candIterator.returnPresAbsStatus()==1:
				self.realCandLine.append(candIterator) #Populate the vector of realized candidate lines
				#log.info("Constructed shared line ID: {}".format(candIterator.getTranslID()))
				candIterator.modifyNodeReact() #Adjust the connected nodal reactances after the lines have been decided to be built
				nodeZone2 = candIterator.getExtNodeID()
				tNodeID2 = candIterator.getExtZoneID()
				indCount = 0 #Initialize a counter for tracking the position in the vector of the iterator
				rankBuilt = candIterator.getExtNodeGlobalRank()
				for diffZNIt in self.diffZoneNodeExistingID:
					if (diffZNIt == tNodeID2) and (self.diffZoneExistingID[indCount] == nodeZone2): #Check whether the other-zone node is already present in the list
						self.containsFlag = 1 #Loop through the diffZoneNodeID and diffZoneID vectors, to see if the other zone-node is already present or not, if yes, set the containsFlag to 1
					indCount += 1 #Increment the counter
				if self.containsFlag == 0: #If the diffZoneNodeID and diffZoneID vectors do not contain the other zone-node, then push the node and the other zone index in the respective vectors
					self.diffZoneNodeExistingID.append(tNodeID2) #initialize the list of external-zone existing node ID's so that it's not empty
					#log.info(" The node ID of outer zonal node is : {}".format(tNodeID2))
					self.diffZoneExistingID.push_back(nodeZone2) #Initialize the list of external-zone existing zone ID's so that it's not empty.
					#log.info(" The zone ID of outer zonal node is : {}".format(nodeZone2))
				for globalIterator in self.globalRankDiffNode:
					if globalIterator == rankBuilt: #Check whether the other-zone node is already present in the list
						self.containsFlagGlob = 1 #Loop through the diffZoneNodeID and diffZoneID vectors, to see if the other zone-node is already present or not, if yes, set the containsFlag to 1
				if self.containsFlagGlob == 0: #If globalRankDiffNode vector does not contain the other zone-node rank, then push the node rank in the vector
					self.globalExistingRank.append(rankBuilt) #list of global ranking of external zone nodes for only SE lines and shared candidate lines
					#cout << " The rank of outer zonal node is : " << rankBuilt << endl;
					self.existingOtherZoneNodeNum += 1 #Increment the number of other zone nodes that are ends of SE lines, only when a new one is added
					#cout << " The number of existing other zone node is : " << existingOtherZoneNodeNum << endl
				self.containsFlag = 0 #reset
				self.containsFlagGlob = 0 #reset
				#(SELineObject[existingCounter]).connectRank(rankFrom)
				self.realizedCLines += 1
		#cout << "Number of constructed shared lines in network " << zonalIndex << " is " << realizedCLines << endl; 

		#Built internal candidate lines
		#Iterator for candidate lines
		for intCandIterator in self.realIntCandLineObject:
			intCandIterator.setPresAbsStatus()
			#log.info("Constructed internal line ID: {}".format(intCandIterator.getTranslID()))
			realizedIntCLines += 1
		for intCandIterator in self.intCandLineObject:
			if (intCandIterator.returnPresAbsStatus()==1) and (find(realIntCandLineObject.begin(), realIntCandLineObject.end(), (*intCandIterator))==realIntCandLineObject.end()):
				realIntCandLineObject.push_back((*intCandIterator)); // Populate the vector of realized candidate lines
				//cout << "Constructed internal line ID: " << (*intCandIterator)->getTranslID() << endl;
				(*intCandIterator)->modifyNodeReact(); // Adjust the connected nodal reactances after the lines have been decided to be built
				realizedIntCLines += 1
		#cout << "Number of constructed internal lines in network " << zonalIndex << " is " << realizedIntCLines << endl

	def TestBuiltExternalNodes(self):
		#cout << " Subnetwork # : " << zonalIndex << endl;
		indCount = 0
		cancelFirstEntry = 0
		for diffZNIt in self.diffZoneNodeExistingID:
			if cancelFirstEntry > 0:
				log.info("\nOuter node number : {} Outer node zone : {} Outer node global rank : {}".format(diffZNIt, self.diffZoneExistingID[indCount], self.globalExistingRank[indCount]))
			indCount += 1 #Increment the counter
			cancelFirstEntry += 1
		#log.info(" Total number of outer zone nodes for this subnetwork to which existing or built lines are connected: {}".format(existingOtherZoneNodeNum))

	def returnMultiplicity(self): #Returns the total multiplicity of the shared nodes
		connMult = 0
		for nodeIterator in self.nodeObject:
			connMult += nodeIterator.getNodeMultiplicity()
		return connMult

	def APPQPAvgHR(Marketover &coordInstanceRef, double LagMultXi[], int totalCandLineNum, int totalSharedNodeNum, GRBEnv* environmentGUROBI, int iterCount): #Calls the GUROBI solver object and solver method to solve the problem of determining the values of the continuous variables


	def getZonalDecision(self): #Returns the intermediate decision variable values from APP
		return self.thetaBuffer

	def getZonalRanks(self): #Returns the global ranks of the shared decision variables from APP
		return self.globRankBuffer
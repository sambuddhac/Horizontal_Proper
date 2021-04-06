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
#define AVERAGE_HEAT 1 // Defines the Average Heat generator cost function mode
#define PIECEWISE_LINEAR 2 // Defines the Piecewise Linear generator cost function mode
#define POLYNOMIAL 3 // Defines the Convex Polynomial generator cost function mode
#define BIGM 1000000000000000000 // Defines the value of the Big M for transforming the bilinear terms to linear constraints

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
	def __init__(self, jSONList, zoneCount, objChoice): # constructor
		self.simMode = objChoice #Specify the type of the curve
		self.zonalCount = zoneCount
		self.nodeNumVector = []
		self.nodeNumVector.append(0) #Initialize the node number vector to indicate there no nodes in the fictitous 0-th zone
		for zonalIndex in jSONList: #Iterate through the zones 
			self.netFile = zonalIndex['Network File'] #String for storing the name of the network file
			self.genFile = zonalIndex['Generator File'] #String for storing the name of the generator file
			self.tranFile = zonalIndex['Transmission Lines File'] #String for storing the name of the transmission line file
			self.loadFile = zonalIndex['Load File'] #String for storing the name of the load file
			self.intCandLineFile = zonalIndex['Intra Candidate Lines File'] #String for storing the name of the candidate lines file
			#/* Nodes */
			matrixNetFile = json.load(open(os.path.join("data", self.netFile)))
			self.nodeNumber = matrixNetFile['nodeNumber']
			self.genNumber = matrixNetFile['genNumber'] 
			self.loadNumber = matrixNetFile['loadNumber'] 
			self.tranNumber = matrixNetFile['tranNumber']
			self.countOfScenarios = matrixNetFile['loadStochScenarios'] #get the number of stochastic scenarios of load#get the dimensions of the Network
			for l in range(self.nodeNumber):
				#log.info("Creating the {} -th Node".format(l + 1)) 
		
				nodeInstance = Node( self.univNodeNum + l + 1, l + 1, zonalIndex ) #creates nodeInstance object with ID l + 1

				self.nodeObject.append( nodeInstance ) #pushes the nodeInstance object into the vector
			#end initialization for Nodes
##/* Generators */
			#/* Instantiate Generators */
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
				genInstance = self.Powergenerator( univGenNum + j+1, nodeObject[ univNodeNum + gNodeID - 1 ],  tanTheta, minCost, PgMax, PgMin )
				self.genObject.append(genInstance) #push the generator object into the array
				j +=1 #increment counter
			matrixGen.close() #Close the generator file
			self.genDF = pd.DataFrame([g.__dict__ for g in self.genObject])
##/* Transmission Lines */
        		if self.tranNumber > 0:
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
					transLineInstance = self.transmissionLine( univTranNum + k + 1, self.nodeObject[ univNodeNum + tNodeID1 - 1 ], self.nodeObject[ univNodeNum + tNodeID2 - 1 ], ptMax, reacT ) 
					self.translObject.append( transLineInstance ) #pushes the transLineInstance object into the vector
					k +=1 #increment the counter
		#end initialization for Transmission Lines 
				matrixTranFile.close() #Close the transmission line file 
				self.tranDF = pd.DataFrame([t.__dict__ for t in self.translObject])
##/* Internal Candidate Transmission Lines */
			if self.internalCLines > 0:
        			matrixIntCETranFile = json.load(open(os.path.join("data", self.intCandLineFile))) #ifstream constructor opens the file of internal candidate Transmission lines
        			k = 0 #Counter for internal candidate lines
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
					intCandLineInstance = self.intCandLine( univIntCandNum + k + 1, self.nodeObject[ univNodeNum + tNodeID1 - 1 ], self.nodeObject[ univNodeNum + tNodeID2 - 1 ], ptMax, reacT, interestRate, lifeTime, costPerCap, presAbsence ) #Create the internal candidate transmission line object with node 1 
					self.intCandLineObject.append( intCandLineInstance ) #pushes the transLineInstance object into the vector
            				k +=1 #increment the counter
			#end initialization for candidate Transmission Lines
				matrixIntCETranFile.close() #Close the candidate lines file
				self.intCandtranDF = pd.DataFrame([t.__dict__ for t in elf.intCandLineObject])
##/* Loads */
        		matrixLoadFile = json.load(open(os.path.join("data", self.loadFile))) #ifstream constructor opens the file of Loads
			countOfScenarios = loadFields-1
			for l in range(self.nodeNumber+self.nodeNumber):
		    		(self.nodeObject[l]).initLoad(self.countOfScenarios) #Initialize the default loads on all nodes to zero
		#end initialization for Nodes		
		#/* Instantiate Loads */
			for j in loadNumber:
				for matrixLoad in matrixLoadFile:
			#node object ID to which the particular load object is connected
			#node ID of the node object to which this load object is connected.
					lNodeID = matrixLoad['loadNodeID']
					for f in range(self.countOfScenarios):
				#value of allowable power consumption capability of load in pu with a negative sign to indicate consumption:
				#Power Consumption:
						P_Load[f] = matrixLoad['loadMW']/100

				loadInstance = self.Load( univLoadNum + j + 1, self.nodeObject[ univNodeNum + lNodeID - 1 ], loadFields-1, P_Load ) #creates loadInstance object object with ID number j + 1

				self.loadObject.append( loadInstance ) #pushes the loadInstance object into the vector
		#end initialization for Loads
			matrixLoadFile.close() #Closes the load file
			for f in range(self.countOfScenarios):
				self.probability.append(1/(self.countOfScenarios))
########################################################################################################################	
		univNodeNum = univNodeNum + nodeNumber #Increment the universal node number
		nodeNumVector.append(nodeNumber)
		univGenNum = univGenNum + genNumber #Increment the universal generator number
		univTranNum = univTranNum + tranNumber #Increment the universal transmission line number
		univLoadNum = univLoadNum + loadNumber #Increment the universal load number
		univIntCandNum = univIntCandNum + internalCLines #Increment the universal intra zonal candidate line number
######################################################################################################################## 
		for zonalIndex in jSONList: #Iterate through the zones 
			self.netFile = zonalIndex['Network File'] 
			self.sharedLineFile = jSONIndex['Shared Lines File'] 
			self.tranFile = zonalIndex['Transmission Lines File'] 
			self.candLineFile = jSONIndex['Candidate Lines File'] 
## Shared Existing Transmission lines
			matrixSETranFile = json.load(open(os.path.join("data", self.sharedLineFile))) #ifstream constructor opens the file of Transmission lines
			k = 0
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
				if serNum in SELineSerList:
					newNodeBase = 0
					for aggrCount in nodeZone1:
						newNodeBase = newNodeBase + nodeNumVector[aggrCount]
					tNodeID1 = tNodeID1 + newNodeBase
					newNodeBase = 0
					for aggrCount in nodeZone2:
						newNodeBase = newNodeBase + nodeNumVector[aggrCount]
					tNodeID2 = tNodeID2 + newNodeBase
					SELineInstance = self.SELine( k + 1, serNum, self.nodeObject[ tNodeID1 - 1 ], self.nodeObject[ tNodeID2 - 1 ], ptMax, reacT )
					SELineObject.append( SELineInstance )
					SELineSerList.append(serNum)
					univSELineNum = univSELineNum + 1
					k +=1
			matrixSETranFile.close() #Close the shared existing file
			self.sharedExistingDF = pd.DataFrame([seDF.__dict__ for seDF in self.SELineObject])
##/* Shared Candidate Transmission Lines */
        		matrixCETranFile = json.load(open(os.path.join("data", self.candLineFile))) #ifstream constructor opens the file of candidate Transmission lines
		#/* Instantiate Shared Candidate Transmission Lines */
			k=0
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
				if serNum in SELineSerList:
					newNodeBase = 0
					for aggrCount in nodeZone1:
						newNodeBase = newNodeBase + nodeNumVector[aggrCount]
					tNodeID1 = tNodeID1 + newNodeBase
					newNodeBase = 0
					for aggrCount in nodeZone2:
						newNodeBase = newNodeBase + nodeNumVector[aggrCount]
					tNodeID2 = tNodeID2 + newNodeBase
					candLineInstance = self.candLine( k + 1, serNum, self.nodeObject[ tNodeID1 - 1 ], self.nodeObject[ tNodeID2 - 1 ], ptMax, reacT, interestRate, lifeTime, costPerCap, presAbsence )
					candLineObject.append( candLineInstance )		# pushes the transLineInstance object into the vector
					candLineSerList.append(serNum)  
					univCandLineNum = univCandLineNum + 1
					k +=1
			matrixCETranFile.close() #Close 
##End of my translation
	assignProb()
	#end constructor
	def assignProb(self):
		for f in range(self.countOfScenarios):
			self.probability.append(1/self.countOfScenarios)

	def MILPAvgHRGUROBI(self): #Function MILPAvgHRGUROBI() implements the Mixed Integer Linear Programming Unit Commitment Solver routine by calling GUROBI routines for average heat rate objective for Horizontal Coordination Investment decision making
		milp_avg_hr_central(matrixNetFile, matrixGenFile, matrixTran, )
 		#Function MILP() ends

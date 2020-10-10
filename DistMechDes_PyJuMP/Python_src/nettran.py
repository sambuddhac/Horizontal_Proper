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
		#/* Nodes */
		matrixNetFile = json.load(open(os.path.join("data", netFile)))
		self.nodeNumber = matrixNetFile['nodeNumber'] 
		self.sharedELines = matrixNetFile['sharedELines'] 
		self.sharedCLines = matrixNetFile['sharedCLines'] 
		self.genNumber = matrixNetFile['genNumber'] 
		self.loadNumber = matrixNetFile['loadNumber'] 
		self.tranNumber = matrixNetFile['tranNumber'] 
		self.internalCLines = matrixNetFile['internalCLines'] #get the dimensions of the Network
		for l in range(self.nodeNumber):
			#log.info("\nCreating the {} -th Node:".format(l + 1))
		
			nodeInstance = Node(l + 1, self.zonalIndex) #creates nodeInstance object with ID l + 1

			self.nodeObject.append( nodeInstance ) #pushes the nodeInstance object into the vector
		#end initialization for Nodes
		matrixNetFile.close() #close the network file
		#/* Generators */
		#/* Instantiate Generators */
		# Open the .json file to read the Powergenerator parameter values
		matrixGen = json.load(open(os.path.join("data", genFile))) #ifstream constructor opens the file of Generators
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
			genInstance = Powergenerator(j+1, nodeObject[ gNodeID - 1 ],  tanTheta, minCost, PgMax, PgMin)
			self.genObject.append(genInstance) #push the generator object into the array
		matrixGenFile.close() #Close the generator file
		#/* Transmission Lines */
		matrixTranFile = json.load(open(os.path.join("data", tranFile))) #ifstream constructor opens the file of Transmission lines
		#/* Instantiate Transmission Lines */
		for matrixTran in matrixTranFile:
			#node IDs of the node objects to which this transmission line is connected.
			tNodeID1 = matrixTran['fromNode'] #From end
			tNodeID2 = matrixTran['toNode'] #To end
			#Reactance
			reacT = matrixTran['Reactance']
			#values of maximum allowable power flow on line in the forward and reverse direction:
			ptMax = matrixTran['lineLimit']/100
			#creates transLineInstance object with ID k + 1
			transLineInstance = transmissionLine(k + 1, nodeObject[ tNodeID1 - 1 ], nodeObject[ tNodeID2 - 1 ], ptMax, reacT) 
			self.translObject.append( transLineInstance ) #pushes the transLineInstance object into the vector
		#end initialization for Transmission Lines 
		matrixTranFile.close() #Close the transmission line file 

		#/* Shared Existing Transmission Lines */
		matrixSETranFile = json.load(open(os.path.join("data", sharedLineFile))) #ifstream constructor opens the file of Transmission lines
		
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
				SELineInstance = SELine(k + 1, serNum, nodeObject[ tNodeID1 - 1 ], tNodeID1, nodeZone1, tNodeID2, nodeZone2, zonalIndex, ptMax, reacT) #Create the shared existing transmission line object with node 1 
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
					otherNodeCount += 1 #Increment the counter to account for the total number of other-zone nodes
					SELineInstance.outerNodeIndex(otherNodeCount, toFromFlag) #Calling the SELine method, which in turn calls the intra-zonal node method, to create a list of rankings of all the outer-zone nodes that are connected to this node, if the intra-zonal node is newly added
				self.SELineObject.append(SELineInstance) #pushes the transLineInstance object into the vector
			else: #Otherwise, if the node 2 belongs to this zone
				SELineInstance = SELine(k + 1, serNum, nodeObject[ tNodeID2 - 1 ], tNodeID1, nodeZone1, tNodeID2, nodeZone2, zonalIndex, ptMax, reacT) #Create the shared existing transmission line object with node 2
				indCount = 0 #Initialize a counter for tracking the position in the vector of the iterator
				toFromFlag = -1 #Indicates that the To node is the intra-zonal node
				for diffZNIt in diffZoneNodeID:
					if diffZNIt == tNodeID1 and self.diffZoneID[indCount] == nodeZone1: #Check whether the other-zone node is already present in the list
						self.containsFlag = 1 #Loop through the diffZoneNodeID and diffZoneID vectors, to see if the other zone-node is already present or not, if yes, set the containsFlag to 1
						SELineInstance.outerNodeIndex(indCount, toFromFlag) #Calling the SELine method, which in turn calls the intra-zonal node method, to create a list of rankings of all the outer-zone nodes that are connected to this node, if the intra-zonal node is already present
					indCount += 1 #Increment the counter
				if self.containsFlag == 0: #If the diffZoneNodeID and diffZoneID vectors do not contain the other zone-node, then push the node and the other zone index in the respective vectors
					self.diffZoneNodeID.append(tNodeID1)
					self.diffZoneID.append(nodeZone1)
					self.diffZoneNodeExistingID.append(tNodeID1) #initialize the list of external-zone existing node ID's so that it's not empty
					self.diffZoneExistingID.append(nodeZone1) #Initialize the list of external-zone existing zone ID's so that it's not empty.
					otherNodeCount += 1  #Increment the counter to account for the total number of other-zone nodes
					SELineInstance.outerNodeIndex(otherNodeCount, toFromFlag) #Calling the SELine method, which in turn calls the intra-zonal node method, to create a list of rankings of all the outer-zone nodes that are connected to this node, if the intra-zonal node is newly added
				self.SELineObject.append( SELineInstance ) # pushes the transLineInstance object into the vector	
			self.containsFlag = 0 #Reset the containsFlag for matching the next item
		#end initialization for Shared Existing Transmission Lines
		matrixSETranFile.close() #Close the shared existing file
		#/* Shared Candidate Transmission Lines */
		matrixCETranFile = json.load(open(os.path.join("data", candLineFile))) #ifstream constructor opens the file of candidate Transmission lines

		#/* Instantiate Shared Candidate Transmission Lines */
		for ( int k = 0; k < sharedCLines; ++k ) {
			int serNum, tNodeID1, tNodeID2, nodeZone1, nodeZone2, presAbsence, lifeTime, ownership; // node object IDs to which the particular transmission line object is connected
			do {
				//node IDs of the node objects to which this transmission line is connected.
				serNum = matrixCETran[ k ][ 0 ]; // global serial number of the shared existing transmission line
				tNodeID1 = matrixCETran[ k ][ 1 ]; // From end node 
				nodeZone1 = matrixCETran[ k ][ 2 ]; // From end zone number
				tNodeID2 = matrixCETran[ k ][ 3 ]; // To end node
				nodeZone2 = matrixCETran[ k ][ 4 ]; // To end zone number 
			} while ( (( nodeZone1 != zonalIndex ) && ( nodeZone2 != zonalIndex )) || (( nodeZone1 == zonalIndex ) && ( nodeZone2 == zonalIndex )) ); // validity check
			double reacT, ptMax, interestRate, costPerCap; // Parameters for Transmission Line
			do {
				//Reactance:
				reacT = matrixCETran[ k ][ 5 ];
				//values of maximum allowable power flow on line in the forward and reverse direction:
				//Forward direction:
				ptMax = matrixCETran[ k ][ 6 ]/100;
				lifeTime = matrixCETran[ k ][ 7 ]; // life time of the candidate line
				interestRate = matrixCETran[ k ][ 8 ]; // interest rate of the investment 
				costPerCap = matrixCETran[ k ][ 9 ]*ptMax; // capital cost for the construction 
				presAbsence = matrixCETran[ k ][ 10 ]; // status of the construction 
				ownership = matrixCETran[ k ][ 11 ]; // ownership of the candidate line for this zone
			} while ( ( reacT <= 0 ) ); // check the bounds and validity of the parameter values
			
			// creates candLineInstance object with ID k + 1
			if ( nodeZone1 == zonalIndex ) {
				candLine *candLineInstance = new candLine( k + 1, serNum, nodeObject[ tNodeID1 - 1 ], tNodeID1, nodeZone1, tNodeID2, nodeZone2, zonalIndex, ptMax, reacT, interestRate, lifeTime, costPerCap, presAbsence, ownership );// Create the shared candidate transmission line object with node 1 
				int indCount = 0; // Initialize a counter for tracking the position in the vector of the iterator
				int toFromFlag = 1; // Indicates that the from node is the intra-zonal node 
				for (diffZNIt = diffZoneNodeID.begin(); diffZNIt != diffZoneNodeID.end(); ++diffZNIt) {
					if ((*diffZNIt == tNodeID2) && (diffZoneID[indCount] == nodeZone2)) {
						containsFlag = 1;  // Loop through the diffZoneNodeID and diffZoneID vectors, to see if the other zone-node is already present or not, if yes, set the containsFlag to 1
						candLineInstance->outerNodeIndex(indCount, toFromFlag); // Calling the SELine method, which in turn calls the intra-zonal node method, to create a list of rankings of all the outer-zone nodes that are connected to this node, if the intra-zonal node is already present
					}
					++indCount; // Increment the counter
				}
				if (containsFlag == 0) { // If the diffZoneNodeID and diffZoneID vectors do not contain the other zone-node, then push the node and the other zone index in the respective vectors
					diffZoneNodeID.push_back(tNodeID2);
					diffZoneID.push_back(nodeZone2);
					++otherNodeCount; // Increment the counter to account for the total number of other-zone nodes
					candLineInstance->outerNodeIndex(otherNodeCount, toFromFlag); // Calling the SELine method, which in turn calls the intra-zonal node method, to create a list of rankings of all the outer-zone nodes that are connected to this node, if the intra-zonal node is newly added
				}
				candLineObject.push_back( candLineInstance ); // pushes the transLineInstance object into the vector
			}
			else {
				candLine *candLineInstance = new candLine( k + 1, serNum, nodeObject[ tNodeID2 - 1 ], tNodeID1, nodeZone1, tNodeID2, nodeZone2, zonalIndex, ptMax, reacT, interestRate, lifeTime, costPerCap, presAbsence, ownership );// Create the shared candidate transmission line object with node 2
				int indCount = 0; // Initialize a counter for tracking the position in the vector of the iterator
				int toFromFlag = -1; // Indicates that the To node is the intra-zonal node 
				for (diffZNIt = diffZoneNodeID.begin(); diffZNIt != diffZoneNodeID.end(); ++diffZNIt) {
					if ((*diffZNIt == tNodeID1) && (diffZoneID[indCount] == nodeZone1)) {
						containsFlag = 1;   // Loop through the diffZoneNodeID and diffZoneID vectors, to see if the other zone-node is already present or not, if yes, set the containsFlag to 1
						candLineInstance->outerNodeIndex(indCount, toFromFlag); // Calling the SELine method, which in turn calls the intra-zonal node method, to create a list of rankings of all the outer-zone nodes that are connected to this node, if the intra-zonal node is already present
					}
					++indCount; // Increment the counter
				}
				if (containsFlag == 0) { // If the diffZoneNodeID and diffZoneID vectors do not contain the other zone-node, then push the node and the other zone index in the respective vectors
					diffZoneNodeID.push_back(tNodeID1);
					diffZoneID.push_back(nodeZone1);
					++otherNodeCount; // Increment the counter to account for the total number of other-zone nodes
					candLineInstance->outerNodeIndex(otherNodeCount, toFromFlag); // Calling the SELine method, which in turn calls the intra-zonal node method, to create a list of rankings of all the outer-zone nodes that are connected to this node, if the intra-zonal node is newly added
				}
				candLineObject.push_back( candLineInstance ); // pushes the transLineInstance object into the vector
			}
			containsFlag = 0; // Reset the containsFlag for matching the next item
		} // end initialization for candidate Transmission Lines
		matrixCETranFile.close(); // Close the candidate lines file
		if(internalCLines>0) {
		/* Internal Candidate Transmission Lines */
		ifstream matrixIntCETranFile( intCandLineFile, ios::in ); // ifstream constructor opens the file of internal candidate Transmission lines
		// exit program if ifstream could not open file
		if ( !matrixIntCETranFile ) {
			cerr << "\nFile for Internal Candidate Transmission lines could not be opened\n" << endl;
			exit( 1 );
		} // end if
		int intCFields; // Number of columns in the internal candidate transmission lines file
		matrixIntCETranFile >> intCFields; // get the dimensions of the internal candidate Transmission line matrix
		double matrixIntCETran[ internalCLines ][ intCFields ]; // internal candidate Transmission line matrix
		for ( int i = 0; i < internalCLines; ++i ) {
			for ( int j = 0; j < intCFields; ++j ) {
				matrixIntCETranFile >> matrixIntCETran[ i ][ j ]; // read the internal candidate Transmission line matrix
			}
		}
		//cout << "\nFile for Internal Candidate Transmission lines read\n" << endl;
		/* Instantiate Internal Candidate Transmission Lines */
		for ( int k = 0; k < internalCLines; ++k ) {
			int serNum, tNodeID1, tNodeID2, presAbsence, lifeTime; // node object IDs to which the particular transmission line object is connected
			do {
				//node IDs of the node objects to which this transmission line is connected.
				tNodeID1 = matrixIntCETran[ k ][ 0 ]; // From end node 
				tNodeID2 = matrixIntCETran[ k ][ 1 ]; // To end node
			} while ( ( tNodeID1 <= 0 ) || ( tNodeID1 > nodeNumber ) || ( tNodeID2 <= 0 ) || ( tNodeID2 > nodeNumber ) || ( tNodeID1 == tNodeID2) ); // validity check // validity check
			double reacT, ptMax, interestRate, costPerCap; // Parameters for Transmission Line
			do {
				//Reactance:
				reacT = matrixIntCETran[ k ][ 2 ];
				//values of maximum allowable power flow on line in the forward and reverse direction:
				//Forward direction:
				ptMax = matrixIntCETran[ k ][ 3 ]/100;
				lifeTime = matrixIntCETran[ k ][ 4 ]; // life time of the candidate line
				interestRate = matrixIntCETran[ k ][ 5 ]; // interest rate of the investment 
				costPerCap = matrixIntCETran[ k ][ 6 ]*ptMax; // capital cost for the construction 
				presAbsence = matrixIntCETran[ k ][ 7 ]; // status of the construction 
			} while ( ( reacT <= 0 ) ); // check the bounds and validity of the parameter values
			
			// creates intCandLineInstance object with ID k + 1
			intCandLine *intCandLineInstance = new intCandLine( k + 1, nodeObject[ tNodeID1 - 1 ], nodeObject[ tNodeID2 - 1 ], ptMax, reacT, interestRate, lifeTime, costPerCap, presAbsence );// Create the internal candidate transmission line object with node 1 
			intCandLineObject.push_back( intCandLineInstance ); // pushes the transLineInstance object into the vector
		} // end initialization for candidate Transmission Lines
		matrixIntCETranFile.close(); // Close the candidate lines file
		}
		/* Loads */
		ifstream matrixLoadFile( loadFile, ios::in ); // ifstream constructor opens the file of Loads

		// exit program if ifstream could not open file
		if ( !matrixLoadFile ) {
			cerr << "\nFile for Loads could not be opened\n" << endl;
			exit( 1 );
		} // end if
		int loadFields; // Number of columns in the load file
		matrixLoadFile >> loadFields; // get the dimensions of the Load matrix
		countOfScenarios = loadFields-1;
		double matrixLoad[ loadNumber ][ loadFields ]; // Load matrix
		for ( int i = 0; i < loadNumber; ++i ) {
			for ( int j = 0; j < loadFields; ++j ) {
				matrixLoadFile >> matrixLoad[ i ][ j ]; // read the Load matrix
			}
		}
		// Initialize the default loads on all nodes to zero
		for ( int l = 0; l < nodeNumber; ++l ) {

			(nodeObject[l])->initLoad( countOfScenarios ); // Initialize the default loads on all nodes to zero

		} // end initialization for Nodes		
		/* Instantiate Loads */
		for ( int j = 0; j < loadNumber; ++j ) {
			int lNodeID; // node object ID to which the particular load object is connected
			do {
				//node ID of the node object to which this load object is connected.
				lNodeID = matrixLoad[ j ][ 0 ]; 
			} while ( ( lNodeID <= 0 ) || ( lNodeID > nodeNumber ) ); // validity check

			double P_Load[loadFields-1]; // Parameters for Load
			for (int f=0;f<(loadFields-1); ++f) {
				do {
					//value of allowable power consumption capability of load in pu with a negative sign to indicate consumption:
					//Power Consumption:
					P_Load[f] = matrixLoad[ j ][ 1+f ]/100;
				} while ( -P_Load[f] <= 0 ); // check the bounds and validity of the parameter values
			}		

			Load *loadInstance = new Load( j + 1, nodeObject[ lNodeID - 1 ], loadFields-1, P_Load ); // creates loadInstance object object with ID number j + 1

			loadObject.push_back( loadInstance ); // pushes the loadInstance object into the vector

		} // end initialization for Loads
		matrixLoadFile.close(); // Closes the load file
	} while ( (genNumber <= 0 ) || ( nodeNumber <= 0 ) || ( loadNumber <= 0 )); //|| ( tranNumber <= 0 ) 
	for (int f=0; f < countOfScenarios; ++f) {
		probability.push_back((static_cast<double>(1)/(countOfScenarios)));
	}
	// check the bounds and validity of the parameter values
	otherZoneNodeIter = diffZoneNodeID.begin(); // Initialize the otherZoneNodeIter iterator to point to the beginning of the diffZoneNodeID vector
	otherZoneNodeIter++; // Increment the pointer to point to the actual non-zero entry
	otherZoneIter = diffZoneID.begin(); // Initialize the otherZoneIter iterator to point to the beginning of the diffZoneID vector
	otherZoneIter++; // Increment the pointer to point to the actual non-zero entry	
	sharedELineIt = SELineObject.begin(); // Initialize the sharedELineIt iterator to point to the beginning of the SELineObject vector 
	sharedELineFromIt = SELineObject.begin(); // Initialize the sharedELineFromIt iterator to point to the beginning of the SELineObject vector
	sharedELineFZoneIt = SELineObject.begin(); // Initialize the sharedELineFZoneIt iterator to point to the beginning of the SELineObject vector
	sharedELineToIt = SELineObject.begin(); // Initialize the sharedELineToIt iterator to point to the beginning of the SELineObject vector
	sharedELineTZoneIt = SELineObject.begin(); // Initialize the sharedELineTZoneIt iterator to point to the beginning of the SELineObject vector
	sharedELineReactIt = SELineObject.begin(); // Initialize the sharedELineReactIt iterator to point to the beginning of the SELineObject vector
	sharedELineCapIt = SELineObject.begin(); // Initialize the sharedELineCapIt iterator to point to the beginning of the SELineObject vector
	sharedCandLineIt = candLineObject.begin(); // Initialize the sharedCandLineIt iterator to point to the beginning of the candLineObject vector 
	sharedCandLineFromIt = candLineObject.begin(); // Initialize the sharedCandLineFromIt iterator to point to the beginning of the candLineObject vector
	sharedCandLineFZoneIt = candLineObject.begin(); // Initialize the sharedCandLineFZoneIt iterator to point to the beginning of the candLineObject vector
	sharedCandLineToIt = candLineObject.begin(); // Initialize the sharedCandLineToIt iterator to point to the beginning of the candLineObject vector
	sharedCandLineTZoneIt = candLineObject.begin(); // Initialize the sharedCandLineTZoneIt iterator to point to the beginning of the candLineObject vector
	sharedCandLineReactIt = candLineObject.begin(); // Initialize the sharedCandLineReactIt iterator to point to the beginning of the candLineObject vector
	sharedCandLineCapIt = candLineObject.begin(); // Initialize the sharedCandLineCapIt iterator to point to the beginning of the candLineObject vector
} // end constructor

	def getNumberOfScenarios(self): #Returns the number of scenarios
		return self.countOfScenarios

	def MILPPiecewiseLin(): #Function MILPPiecewiseLin() implements the Mixed Integer Linear Programming Unit Commitment Solver routine by calling GLPK routines for piecewise linear objective
		log.info("\nUnder Construction")

	def MILPPolynomial(): #Function MILPPolynomial() implements the Mixed Integer Linear Programming Unit Commitment Solver routine by calling GLPK routines for polynomial convex objective
		log.info("\nUnder Construction")

	def MILPAvgHR(self, coordInstanceRef, LagMultXi, LagMultPi, totalCandLineNum, totalSharedNodeNum): #Function MILPAvgHR() implements the Mixed Integer Linear Programming Unit Commitment Solver routine by calling GLPK routines for average heat rate objective for Horizontal Coordination Investment decision making
		zonal_dec_stageI = julSol.HorMILPDistMech(coordInstanceRef, LagMultXi, LagMultPi, totalCandLineNum, totalSharedNodeNum)
	#/* CREATION OF THE MIP SOLVER INSTANCE */
	
} // Function MILP() ends

double Nettran::calcMILPBounds(double LagMultXi[], double LagMultPi[], int totalCandLineNum, int totalSharedNodeNum) // Function MILPAvgHR() implements the Mixed Integer Linear Programming Unit Commitment Solver routine by calling GLPK routines for average heat rate objective for Horizontal Coordination Investment decision making
{
	/* CREATION OF THE MIP SOLVER INSTANCE */
	clock_t begin = clock(); // start the timer
	vector<int>::iterator diffZNIt; // Iterator for diffZoneNodeID
	vector<Powergenerator*>::iterator genIterator; // Iterator for Powergenerator objects
	vector<transmissionLine*>::iterator tranIterator; // Iterator for Transmission line objects
	vector<Load*>::iterator loadIterator; // Iterator for load objects
	vector<Node*>::iterator nodeIterator; // Iterator for node objects
	vector<candLine*>::iterator candIterator; // Iterator for candidate lines
	vector<SELine*>::iterator exsharedIterator; // Iterator for shared existing lines
        vector<intCandLine*>::iterator intCandIterator; // Iterator for candidate lines

	string outSummaryFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputSummaryGLPKBounds/OutBoundSummaryGLPK" + to_string(zonalIndex) + ".txt";
	ofstream outPutFile(outSummaryFileName, ios::out); // Create Output File to output the Summary of Results
	if (!outPutFile){
		cerr << "\nCouldn't open the file" << endl;
		exit(1);
	}

        int dimRow =  countOfScenarios*(2 * genNumber + 4 * sharedCLines + 2 * sharedELines + 2 * tranNumber + nodeNumber + 4*internalCLines); // Total number of rows of the A matrix (number of structural constraints of the LP) first term to account for lower and upper generating limits, second term for lower and upper line limits & lower and upper definition limits of candidate shared lines, third term for lower and upper line limits for shared existing lines, fourth term for lower and upper line limits for internal zonal lines, the fifth term to account for nodal power balance constraints, and sixth term to account for the internal candidate lines
        int dimCol = countOfScenarios*(genNumber+nodeNumber+sharedCLines+otherNodeCount+internalCLines)+sharedCLines+internalCLines; // Total number of columns of the LP (number of Decision Variables) first term to account for power generation MW outputs, second term for voltage phase angles for internal zonal nodes, third term for power flow values and binary integer decision variable values for shared candidate lines, fourth term for the voltage phase angles of other-zone nodes connected through shared existing and candidate lines, and fifth term for the decision variables for internal candidate lines
	outPutFile << "\nTotal Number of Structural Constraints (Rows) is: " << dimRow << endl;
	outPutFile << "\nTotal Number of Decision Variables (Columns) is: " << dimCol << endl;
	glp_prob *milp; // Instantiate GLPK Problem Object pointer
	glp_iocp *ipControlParam = new glp_iocp; // Instantiate the Control parameters for the Integer Programming Problem
	glp_init_iocp(ipControlParam); // Initialize the Control Parameters for the Integer Programming Problem with default values
	ipControlParam->mip_gap = 1e-1; // Set the tolerance for the Integer Programming Problem

	// arrays to store the row and column index combinations for the coefficient matrix
	vector<int> ia; // array to store the non zero element row indices of A matrix
	vector<int> ja; // array to store the non-zero element column indices of A matrix
	vector<double> ar; // array to store the coefficients of the A matrix
	double z; // variable to store the objective value
	vector<double> x; // Coefficient matrix entires, objective function, and decision variables

	milp = glp_create_prob(); // Creates the GLPK MILP Problem
	glp_set_prob_name(milp, "zonalTransDec"); // Names the particular problem instance 
	glp_set_obj_dir(milp, GLP_MIN); // Set direction (Declares the MILP Problem as a Minimization Problem)

	/* SPECIFICATION OF PROBLEM PARAMETERS */
	/*Row Definitions: Specification of RHS or b vector of b<=Ax<=b*/
	glp_add_rows(milp, dimRow);
	//Row Definitions and Bounds Corresponding to Constraints/

	/*******************************************************************************************/

	// Constraints corresponding to supply-demand balance
	string outPGenFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputPowerGLPKBounds/OutBoundPowerGenGLPK" + to_string(zonalIndex) + ".txt"; 
	ofstream powerGenOut(outPGenFileName, ios::out);
	if (!powerGenOut){
		cerr << "\nCouldn't open the file" << endl;
		exit(1);
	}
	// Vectors for storing the output data
	vector<int> busCount; // vector for storing the node/bus serial
	outPutFile << "Constraints corresponding to Supply-Demand Balance right hand side" << endl;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (int rCount = 1; rCount <= nodeNumber; ++rCount){
			busCount.push_back(rCount);
			glp_set_row_name(milp, scenCounter*nodeNumber+rCount, NULL); // Specify the particular supply-demand balance constraint for the particular node in GLPK
			if (((nodeObject[rCount-1])->devpinitMessage(scenCounter))==0)
				glp_set_row_bnds(milp, scenCounter*nodeNumber+rCount, GLP_FX, ((nodeObject[rCount-1])->devpinitMessage(scenCounter)), 0.0); // Specify the right hand side which is the net demand for the particular node
			else 
				glp_set_row_bnds(milp, scenCounter*nodeNumber+rCount, GLP_FX, -((nodeObject[rCount-1])->devpinitMessage(scenCounter)), 0.0); // Specify the right hand side which is the net demand for the particular node
			outPutFile << "Connected load to node " << rCount << " in scenario " << scenCounter+1 << " is " << (nodeObject[rCount-1])->devpinitMessage(scenCounter)*100 << " MW" << endl;
			outPutFile << rCount << "\t";
			if (((nodeObject[rCount-1])->devpinitMessage(scenCounter))==0)
				outPutFile << ((nodeObject[rCount-1])->devpinitMessage(scenCounter))*100 << " MW" << endl;
			else
				outPutFile << -((nodeObject[rCount-1])->devpinitMessage(scenCounter))*100 << " MW" << endl;
		}
	}
	/*******************************************************************************************/

	// Constraints corresponding to Powergenerator Lower Bound
	outPutFile << "Constraints corresponding to Powergenerator Lower Bound" << endl;
	int genRun = 0;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator) {
			++genRun;
			int rCount = countOfScenarios*nodeNumber+genRun;
			glp_set_row_name(milp, rCount, NULL);
			glp_set_row_bnds(milp, rCount, GLP_LO, ((*genIterator)->getPMin()), 0.0);
			outPutFile << rCount << "\t";
			outPutFile << ((*genIterator)->getPMin())*100 << " MW" << endl;
		}
	}
	// Constraints corresponding to Powergenerator Upper Bound
	outPutFile << "Constraints corresponding to Powergenerator Upper Bound" << endl;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator) {
			++genRun;
			int rCount = countOfScenarios*nodeNumber+genRun;
			glp_set_row_name(milp, rCount, NULL);
			glp_set_row_bnds(milp, rCount, GLP_UP, 0.0, ((*genIterator)->getPMax()));
			outPutFile << rCount << "\t";
			outPutFile << ((*genIterator)->getPMax())*100 << " MW" << endl;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for Constraints corresponding to Powergenerator Bounds: " << genRun << endl;
	/*******************************************************************************************/

	// Constraints corresponding to Line Forward Flow limits for internal zonal transmission lines
	outPutFile << "\nTesting of intrazonal transmission line forward flow limits" << endl; 
	int rowCount = countOfScenarios*(nodeNumber+2*genNumber)+1;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (tranIterator = translObject.begin(); tranIterator != translObject.end(); ++tranIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_UP, 0.0, (*tranIterator)->getFlowLimit());
			outPutFile << rowCount << "\t";
			outPutFile << ((*tranIterator)->getFlowLimit())*100 << " MW" << endl;
			++rowCount;
		}
	}
	outPutFile << endl;
	// Constraints corresponding to Line Reverse Flow limits for internal zonal transmission lines
	outPutFile << "\nTesting of intrazonal transmission line reverse flow limits" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (tranIterator = translObject.begin(); tranIterator != translObject.end(); ++tranIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_LO, -((*tranIterator)->getFlowLimit()), 0.0);
			outPutFile << rowCount << "\t";
			outPutFile << -((*tranIterator)->getFlowLimit())*100 << " MW" << endl;
			++rowCount;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for Line Flow limits for internal zonal transmission lines: " << rowCount << endl;

	/*******************************************************************************************/

	// Constraints corresponding to Line Forward Flow limits for shared existing lines
	outPutFile << "\nTesting of Shared Existing transmission line forward flow limits" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (exsharedIterator = SELineObject.begin(); exsharedIterator != SELineObject.end(); ++exsharedIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_UP, 0.0, (*exsharedIterator)->getFlowLimit());
			outPutFile << rowCount << "\t";
			outPutFile << ((*exsharedIterator)->getFlowLimit())*100 << " MW" << endl;
			++rowCount;
		}
	}
	outPutFile << endl;
	// Constraints corresponding to Line Reverse Flow limits for shared existing lines
	outPutFile << "\nTesting of Shared Existing transmission line reverse flow limits" << endl;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (exsharedIterator = SELineObject.begin(); exsharedIterator != SELineObject.end(); ++exsharedIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_LO, -((*exsharedIterator)->getFlowLimit()), 0.0);
			outPutFile << rowCount << "\t";
			outPutFile << -((*exsharedIterator)->getFlowLimit())*100 << " MW" << endl;
			++rowCount;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for Line Flow limits for shared existing lines: " << rowCount << endl;

	/*******************************************************************************************/

	// Constraints corresponding to Line Forward Flow limits for shared candidate lines
	outPutFile << "\nTesting of Shared Candidate transmission line forward flow limits" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_UP, 0.0, 0.0);
			outPutFile << rowCount << "\t";
			outPutFile << 0.0 << endl;
			++rowCount;
		}
	}
	outPutFile << endl;
	// Constraints corresponding to Line Reverse Flow limits for shared candidate lines
	outPutFile << "\nTesting of Shared Candidate transmission line reverse flow limits" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_LO, 0.0, 0.0);
			outPutFile << rowCount << "\t";
			outPutFile << 0.0 << endl;
			++rowCount;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for Line Flow limits for shared candidate lines: " << rowCount << endl;

	/*******************************************************************************************/

	// Constraints corresponding to shared candidate lines flow definition upper bound
	outPutFile << "\nTesting of Definition of Flows on Shared Candidate transmission lines upper bound" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_UP, 0.0, 2.5*((*candIterator)->getFlowLimit()));
			outPutFile << rowCount << "\t";
			outPutFile << BIGM << endl;
			++rowCount;
		}
	}
	outPutFile << endl;
	// Constraints corresponding to shared candidate lines flow definition lower bound
	outPutFile << "\nTesting of Definition of Flows on Shared Candidate transmission lines lower bound" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_LO, -2.5*((*candIterator)->getFlowLimit()), 0.0);
			outPutFile << rowCount << "\t";
			outPutFile << -BIGM << endl;
			++rowCount;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for shared candidate lines flow definition: " << rowCount << endl;
	/*******************************************************************************************/
	// Constraints corresponding to Line Forward Flow limits for internal candidate lines
	outPutFile << "\nTesting of Internal Candidate transmission line forward flow limits" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_UP, 0.0, 0.0);
			outPutFile << rowCount << "\t";
			outPutFile << 0.0 << endl;
			++rowCount;
		}
	}
	outPutFile << endl;
	// Constraints corresponding to Line Reverse Flow limits for internal candidate lines
	outPutFile << "\nTesting of Internal Candidate transmission line reverse flow limits" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_LO, 0.0, 0.0);
			outPutFile << rowCount << "\t";
			outPutFile << 0.0 << endl;
			++rowCount;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for Line Flow limits for internal candidate lines: " << rowCount << endl;

	/*******************************************************************************************/

	// Constraints corresponding to internal candidate lines flow definition upper bound
	outPutFile << "\nTesting of Definition of Flows on Internal Candidate transmission lines upper bound" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_UP, 0.0, 2.5*((*intCandIterator)->getFlowLimit()));
			outPutFile << rowCount << "\t";
			outPutFile << BIGM << endl;
			++rowCount;
		}
	}
	outPutFile << endl;
	// Constraints corresponding to internal candidate lines flow definition lower bound
	outPutFile << "\nTesting of Definition of Flows on Internal Candidate transmission lines lower bound" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_LO, -2.5*((*intCandIterator)->getFlowLimit()), 0.0);
			outPutFile << rowCount << "\t";
			outPutFile << -BIGM << endl;
			++rowCount;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for internal candidate lines flow definition: " << rowCount << endl;
	outPutFile << "\nConstraint bounds (rows) Specified" << endl;
	outPutFile << "\nTotal number of rows: " << rowCount - 1 << endl;
	/*******************************************************************************************/

	/*******************************************************************************************/

	/*Column Definitions, Bounds, and Objective Function Co-efficients*/
	glp_add_cols(milp, dimCol);
	int colCount = 1;
	//Columns corresponding to Power Generation continuous variables for different generators//
	outPutFile << "\nCoefficients of Power generator variables" << endl;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
			glp_set_col_kind(milp, colCount, GLP_CV);
			glp_set_col_name(milp, colCount, NULL);
			glp_set_col_bnds(milp, colCount, GLP_LO, 0.0, 0.0);
			glp_set_obj_coef(milp, colCount, probability.at(scenCounter)*((*genIterator)->getLinCoeff()));
			outPutFile << colCount << "\t";
			outPutFile << ((*genIterator)->getLinCoeff())/100 << " $/MW" << endl;
			++colCount;
		}
	}
	outPutFile << "\nTotal number of columns after accounting for Power Generation continuous variables for different generators: " << colCount << endl;
	/*******************************************************************************************/

	//Columns corresponding to Voltage Phase Angles continuous variables for different intrazonal nodes//
	outPutFile << "\nCoefficients of Voltage Phase Angles continuous variables for different intrazonal nodes" << endl;
	outPutFile << "\nVariable Count\tShared Node\tGlobal Rank\tLagMultXiIndex\tLagMultXiValue" << endl;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {	
		for (nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); ++nodeIterator){
			glp_set_col_kind(milp, colCount, GLP_CV);
			glp_set_col_name(milp, colCount, NULL);
			glp_set_col_bnds(milp, colCount, GLP_DB, 0, (44/7));
			if ( (((*nodeIterator)->getSharedFlag()) == 1) || (((*nodeIterator)->getCandFlag()) == 1) ) {
				glp_set_obj_coef(milp, colCount, LagMultXi[scenCounter*zonalCount*totalSharedNodeNum+(zonalIndex-1)*totalSharedNodeNum+((*nodeIterator)->getGlobalRank())]);
				outPutFile << colCount << "\tYes\t" << ((*nodeIterator)->getGlobalRank()) << "\t" << scenCounter*zonalCount*totalSharedNodeNum+(zonalIndex-1)*totalSharedNodeNum+((*nodeIterator)->getGlobalRank()) << "\t" << LagMultXi[scenCounter*zonalCount*totalSharedNodeNum+(zonalIndex-1)*totalSharedNodeNum+((*nodeIterator)->getGlobalRank())] << endl;	
			}
			else {
				glp_set_obj_coef(milp, colCount, 0.0);
				outPutFile << colCount << "\tNo\t" << "-" << "\t" << "-" << "\t" << "-" << endl;	
			}
			++colCount;
		}
	}
	outPutFile << "\nTotal number of columns after accounting for Voltage Phase Angles continuous variables for different intrazonal nodes: " << colCount << endl;
	/*******************************************************************************************/

	//Columns corresponding to Shared Candidate Line Flows continuous variables//
	outPutFile << "\nCoefficients corresponding to Shared Candidate Line Flows continuous variables" << endl;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			glp_set_col_kind(milp, colCount, GLP_CV);
			glp_set_col_name(milp, colCount, NULL);
			glp_set_col_bnds(milp, colCount, GLP_FR, 0.0, 0.0);
			glp_set_obj_coef(milp, colCount, 0.0);
			outPutFile << colCount << "\t";
			outPutFile << 0 << endl;
			++colCount;
		}
	}
	outPutFile << "\nTotal number of columns after accounting for Shared Candidate Line Flows continuous variables: " << colCount << endl;
	/*******************************************************************************************/

	//Columns corresponding to Shared Candidate Line Construction Decision Binary Integer variables//
	outPutFile << "\nCoefficients corresponding to Shared Candidate Line Construction Decision Binary Integer variables" << endl;
	outPutFile << "\nVariable Count\tGlobal Rank\tLagMultPiIndex\tLagMultPiValue\tInvestment Cost" << endl;
	for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
		glp_set_col_kind(milp, colCount, GLP_CV);
		glp_set_col_name(milp, colCount, NULL);
		glp_set_col_bnds(milp, colCount, GLP_DB, 0.0, 1.0);
		glp_set_obj_coef(milp, colCount, LagMultPi[(zonalIndex-1)*totalCandLineNum+((*candIterator)->getGlobalRank())]+((*candIterator)->returnOwnership())*0.5*((*candIterator)->getInvestCost()));
		outPutFile << colCount << "\t" << ((*candIterator)->getGlobalRank()) << "\t" << (zonalIndex-1)*totalCandLineNum+((*candIterator)->getGlobalRank()) << "\t" << LagMultPi[(zonalIndex-1)*totalCandLineNum+((*candIterator)->getGlobalRank())] << "\t" << ((*candIterator)->returnOwnership())*0.5*((*candIterator)->getInvestCost()) << endl;	
		++colCount;
	}
	outPutFile << "\nTotal number of columns after accounting for Shared Candidate Line Construction Decision Binary Integer variables: " << colCount << endl;
	/*******************************************************************************************/
	/*int diffNodeTestCounter = 0;
	for (globalIterator = globalRankDiffNode.begin(); globalIterator != globalRankDiffNode.end(); ++globalIterator){
		if (diffNodeTestCounter > 0) { // Skip the first element, since it's a dummy "0"	
			cout << " Column count after the shared node test is " << diffNodeTestCounter << endl;
		}
		++diffNodeTestCounter;
	}*/
	//Columns corresponding to Voltage Phase Angles continuous variables for other zone nodes for shared lines//
	outPutFile << "\nCoefficients corresponding to Voltage Phase Angles continuous variables for other zone nodes for shared lines" << endl;
	outPutFile << "\nVariable Count\tGlobal Rank\tLagMultXiIndex\tLagMultXiValue" << endl;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int diffNodeCounter = 0; // counter flag to indicate the first element of the diffZoneNodeID list
		for (globalIterator = globalRankDiffNode.begin(); globalIterator != globalRankDiffNode.end(); ++globalIterator){
			if (diffNodeCounter > 0) { // Skip the first element, since it's a dummy "0"
				glp_set_col_kind(milp, colCount, GLP_CV);
				glp_set_col_name(milp, colCount, NULL);
				glp_set_col_bnds(milp, colCount, GLP_DB, 0, (44/7));
				glp_set_obj_coef(milp, colCount, LagMultXi[scenCounter*zonalCount*totalSharedNodeNum+(zonalIndex-1)*totalSharedNodeNum+(*globalIterator)]);
				outPutFile << colCount << (*globalIterator) << "\t" << scenCounter*zonalCount*totalSharedNodeNum+(zonalIndex-1)*totalSharedNodeNum+(*globalIterator) << "\t" << LagMultXi[scenCounter*zonalCount*totalSharedNodeNum+(zonalIndex-1)*totalSharedNodeNum+(*globalIterator)] << endl;		
				++colCount;
				//cout << " Column count after the shared node " << diffNodeCounter << " is " << colCount << endl;
			}
			++diffNodeCounter;
		}
	}
	outPutFile << "\nTotal number of columns after accounting for Voltage Phase Angles continuous variables for other zone nodes for shared lines: " << colCount << endl;
	/*******************************************************************************************/

	//Columns corresponding to Internal Candidate Line Flows continuous variables//
	outPutFile << "\nCoefficients corresponding to Internal Candidate Line Flows continuous variables" << endl;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			glp_set_col_kind(milp, colCount, GLP_CV);
			glp_set_col_name(milp, colCount, NULL);
			glp_set_col_bnds(milp, colCount, GLP_FR, 0.0, 0.0);
			glp_set_obj_coef(milp, colCount, 0.0);
			outPutFile << colCount << "\t";
			outPutFile << 0 << endl;
			++colCount;
		}
	}
	outPutFile << "\nTotal number of columns after accounting for Internal Candidate Line Flows continuous variables: " << colCount << endl;
	/*******************************************************************************************/

	//Columns corresponding to Internal Candidate Line Construction Decision Binary Integer variables//
	outPutFile << "\nCoefficients corresponding to Internal Candidate Line Construction Decision Binary Integer variables" << endl;
	outPutFile << "\nVariable Count\tInvestment Cost" << endl;
	for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
		glp_set_col_kind(milp, colCount, GLP_CV);
		glp_set_col_name(milp, colCount, NULL);
		glp_set_col_bnds(milp, colCount, GLP_DB, 0.0, 1.0);
		glp_set_obj_coef(milp, colCount, ((*intCandIterator)->getInvestCost()));
		outPutFile << colCount << "\t" << ((*intCandIterator)->getInvestCost()) << endl;	
		++colCount;
	}
	outPutFile << "\nTotal number of columns after accounting for Internal Candidate Line Construction Decision Binary Integer variables: " << colCount << endl;
	outPutFile << "\nTotal Number of columns for generation, angles, integer variables, and flows: " << colCount - 1 << endl;
	outPutFile << "\nDecision Variables and Objective Function defined" << endl;
	outPutFile << "\nTotal Number of columns: " << colCount - 1 << endl;
	/*******************************************************************************************/

	/*******************************************************************************************/

	//Non-Zero entries of A matrix (Constraint/Coefficient matrix entries)//
	int index = 1;
	ia.push_back(0), ja.push_back(0), ar.push_back(0.0);
	// Coefficients for the supply-demand balance constraints
	outPutFile << "\nNon-zero elements of A matrix" << endl;
	outPutFile << "\nRow Number\tColumn Number\tNon-zero Entry\tFrom Reactance\tToReactance" << endl;
	outPutFile << "\nCoefficients for the supply-demand balance constraints" << endl;
	int rCount = 1; // Initialize the row count
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); ++nodeIterator){
			outPutFile << "\nGeneration\t" << rCount << "\n";
			int genListLength = (*nodeIterator)->getGenLength(); // get the number
			cout << "\nScenario Count: " << scenCounter << " node count: " << (*nodeIterator)->getNodeID() << " Conn Gen: " << genListLength << endl;
			for (int cCount = 1; cCount <= genListLength; ++cCount){
				cout << "\nSerial: " << cCount << " Generator Serial: " << (*nodeIterator)->getGenSer(cCount) << endl;
				ia.push_back(rCount), ja.push_back(scenCounter*genNumber+(*nodeIterator)->getGenSer(cCount)), ar.push_back(1.0);
				outPutFile << "\n" << rCount << "\t" << scenCounter*genNumber+(*nodeIterator)->getGenSer(cCount) << "\t" << 1.0 << endl;
				++index;
			}
			outPutFile << "\nIntrazonal Node Angles\t" << rCount << "\n";
			ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber+rCount), ar.push_back(((*nodeIterator)->getToReact())-((*nodeIterator)->getFromReact()));
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber+rCount << "\t" << ((*nodeIterator)->getToReact())-((*nodeIterator)->getFromReact()) << "\t" << ((*nodeIterator)->getFromReact()) << ((*nodeIterator)->getToReact()) << endl;
			++index;
			outPutFile << "\nConnected Intrazonal Node Angles\t" << rCount << "\n";
			int connNodeListLength = (*nodeIterator)->getConNodeLength(); // get the number of intra-zonal nodes connected to this node
			for (int cCount = 1; cCount <= connNodeListLength; ++cCount){
				if (((*nodeIterator)->getConnReact(cCount))<=0) {
					ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber+scenCounter*nodeNumber+((*nodeIterator)->getConnSer(cCount))), ar.push_back(-((*nodeIterator)->getConnReact(cCount)));
					outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber+scenCounter*nodeNumber+((*nodeIterator)->getConnSer(cCount)) << "\t" <<  (((*nodeIterator)->getConnReact(cCount))) << "\n";
                                }
				else {
					ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber+scenCounter*nodeNumber+((*nodeIterator)->getConnSer(cCount))), ar.push_back(((*nodeIterator)->getConnReact(cCount)));
					outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber+scenCounter*nodeNumber+((*nodeIterator)->getConnSer(cCount)) << "\t" <<  (((*nodeIterator)->getConnReact(cCount))) << "\n";
                                }
				++index;

			}
			outPutFile << "\nConnected Outer zonal Node Angles\t" << rCount << "\n";
			int connOutNodeLength = (*nodeIterator)->getExtraNodeLength(); // get the number of outer-zonal nodes connected to this node
			for (int cCount = 1; cCount <= connOutNodeLength; ++cCount){
                                if (((*nodeIterator)->getExtConnReact(cCount))<=0) {
					ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount+(*nodeIterator)->getExtConnSer(cCount)), ar.push_back(-((*nodeIterator)->getExtConnReact(cCount)));
					outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount+(*nodeIterator)->getExtConnSer(cCount) << "\t" <<  (-((*nodeIterator)->getExtConnReact(cCount))) << "\n";
				}
                                else {
					ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount+(*nodeIterator)->getExtConnSer(cCount)), ar.push_back(((*nodeIterator)->getExtConnReact(cCount)));
					outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount+(*nodeIterator)->getExtConnSer(cCount) << "\t" <<  (((*nodeIterator)->getExtConnReact(cCount))) << "\n";
				}
				++index;
			}
			outPutFile << "\nConnected Candidate Lines for which this is the From node\t" << rCount << "\n";
			int connCandListLengthF = (*nodeIterator)->getCandLineLengthF(); // get the number of candidate lines connected to this from node 
			for (int cCount = 1; cCount <= connCandListLengthF; ++cCount){
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber)+scenCounter*sharedCLines+(*nodeIterator)->getCandSerF(cCount)), ar.push_back(-1);
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber)+scenCounter*sharedCLines+(*nodeIterator)->getCandSerF(cCount) << "\t" << -1.0 << "\n";
				++index;
			}
			outPutFile << "\nConnected Candidate Lines for which this is the To node\t" << rCount << "\n";
			int connCandListLengthT = (*nodeIterator)->getCandLineLengthT(); // get the number of candidate lines connected to this to node 
			for (int cCount = 1; cCount <= connCandListLengthT; ++cCount){
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber)+scenCounter*sharedCLines+(*nodeIterator)->getCandSerT(cCount)), ar.push_back(1);
				outPutFile << "\n" << rCount << "\t" << genNumber+nodeNumber+(*nodeIterator)->getCandSerT(cCount) << "\t" << 1.0 << "\n";
				++index;
			}
                	outPutFile << "\nConnected Internal Candidate Lines for which this is the From node\t" << rCount << "\n";
                	int connintCandListLengthF = (*nodeIterator)->getIntCandLineLengthF(); // get the number of internal candidate lines connected to this from node 
                	for (int cCount = 1; cCount <= connintCandListLengthF; ++cCount){
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines+otherNodeCount)+sharedCLines+scenCounter*internalCLines+(*nodeIterator)->getIntCandSerF(cCount)), ar.push_back(-1);			
                        	outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines+otherNodeCount)+sharedCLines+scenCounter*internalCLines+(*nodeIterator)->getIntCandSerF(cCount) << "\t" << -1.0 << "\n";
				++index;
               		}
               		outPutFile << "\nConnected Internal Candidate Lines for which this is the To node\t" << rCount << "\n";
                	int connintCandListLengthT = (*nodeIterator)->getIntCandLineLengthT(); // get the number of internal candidate lines connected to this to node 
                	for (int cCount = 1; cCount <= connintCandListLengthT; ++cCount){
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines+otherNodeCount)+sharedCLines+scenCounter*internalCLines+(*nodeIterator)->getIntCandSerT(cCount)), ar.push_back(1);
                        	outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines+otherNodeCount)+sharedCLines+scenCounter*internalCLines+(*nodeIterator)->getIntCandSerT(cCount) << "\t" << 1.0 << "\n";
				++index;
                	}
			++rCount; // Increment the row count to point to the next node object
		}
	}
	/*******************************************************************************************/

	// Coefficients corresponding to lower generation limits
	outPutFile << "\nCoefficients corresponding to lower generation limits\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
			ia.push_back(rCount), ja.push_back(rCount - countOfScenarios*nodeNumber), ar.push_back(1.0);
			++index;
			outPutFile << rCount << "\t" << (rCount - countOfScenarios*nodeNumber) << "\t" << 1.0 << "\t" << (*genIterator)->getPMin() << endl;
			++rCount; // Increment the row count to point to the next generator object
		}
	}
	// Coefficients corresponding to upper generation limits
	outPutFile << "\nCoefficients corresponding to lower generation limits\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
			ia.push_back(rCount), ja.push_back(rCount - countOfScenarios*(genNumber + nodeNumber)), ar.push_back(1.0);
			++index;
			outPutFile << rCount << "\t" << (rCount - countOfScenarios*(genNumber + nodeNumber)) << "\t" << 1.0 << "\t" << ((*genIterator)->getPMax()) << endl;
			++rCount; // Increment the row count to point to the next generator object
		}
	}
	/*******************************************************************************************/
	
	// Coefficients corresponding to intra-zone Line Forward Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to intra-zone Line Forward Flow Limit Constraints\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (tranIterator = translObject.begin(); tranIterator != translObject.end(); ++tranIterator){
			ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*tranIterator)->getTranslNodeID1()), ar.push_back(1/((*tranIterator)->getReactance()));
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*tranIterator)->getTranslNodeID1() << "\t" << 1/((*tranIterator)->getReactance()) << "\t" << 1/((*tranIterator)->getReactance()) << "\n";
			++index;
			ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*tranIterator)->getTranslNodeID2()), ar.push_back(-1/((*tranIterator)->getReactance()));
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*tranIterator)->getTranslNodeID2() << "\t" << -1/((*tranIterator)->getReactance()) << "\t" << "-" << "\t" << -1/((*tranIterator)->getReactance()) << "\n";
			++index;
			++rCount; // Increment the row count to point to the next transmission line object
		
		}
	}
	// Coefficients corresponding to intra-zone Line Reverse Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to intra-zone Line Reverse Flow Limit Constraints\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (tranIterator = translObject.begin(); tranIterator != translObject.end(); ++tranIterator){
			ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*tranIterator)->getTranslNodeID1()), ar.push_back(1/((*tranIterator)->getReactance()));
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*tranIterator)->getTranslNodeID1() << "\t" << 1/((*tranIterator)->getReactance()) << "\t" << 1/((*tranIterator)->getReactance()) << "\n";
			++index;
			ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*tranIterator)->getTranslNodeID2()), ar.push_back(-1/((*tranIterator)->getReactance()));
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*tranIterator)->getTranslNodeID2() << "\t" << -1/((*tranIterator)->getReactance()) << "\t" << "-" << "\t" << -1/((*tranIterator)->getReactance()) << "\n";
			++index;
			++rCount; // Increment the row count to point to the next transmission line object
		}
	}
	/*******************************************************************************************/
	
	// Coefficients corresponding to shared existing Line Forward Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to shared existing Line Forward Flow Limit Constraints\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (exsharedIterator = SELineObject.begin(); exsharedIterator != SELineObject.end(); ++exsharedIterator){
			if ((*exsharedIterator)->getFlowDir() == 1) {
				ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*exsharedIterator)->getIntlNodeID()), ar.push_back(1/((*exsharedIterator)->getReactance()));
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*exsharedIterator)->getIntlNodeID() << "\t" << 1/((*exsharedIterator)->getReactance()) << "\n";
				++index;
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount+(*exsharedIterator)->getExtNodeRank()), ar.push_back(-1/((*exsharedIterator)->getReactance()));
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*exsharedIterator)->getExtNodeRank() << "\t" << -1/((*exsharedIterator)->getReactance()) << "\n";
				++index;
				++rCount; // Increment the row count to point to the next transmission line object
			}
			else {
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*exsharedIterator)->getExtNodeRank()), ar.push_back(1/((*exsharedIterator)->getReactance()));
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*exsharedIterator)->getExtNodeRank() << "\t" << 1/((*exsharedIterator)->getReactance()) << "\n";
				++index;
				ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber + (*exsharedIterator)->getIntlNodeID()), ar.push_back(-1/((*exsharedIterator)->getReactance()));
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber + (*exsharedIterator)->getIntlNodeID() << "\t" << -1/((*exsharedIterator)->getReactance()) << "\n";
				++index;
				++rCount; // Increment the row count to point to the next transmission line object
			}
		}
	}
	// Coefficients corresponding to shared existing Line Reverse Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to shared existing Line Reverse Flow Limit Constraints\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (exsharedIterator = SELineObject.begin(); exsharedIterator != SELineObject.end(); ++exsharedIterator){
			if ((*exsharedIterator)->getFlowDir() == 1) {
				ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber + (*exsharedIterator)->getIntlNodeID()), ar.push_back(1/((*exsharedIterator)->getReactance()));
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber + (*exsharedIterator)->getIntlNodeID() << "\t" << 1/((*exsharedIterator)->getReactance()) << "\n";
				++index;
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*exsharedIterator)->getExtNodeRank()), ar.push_back(-1/((*exsharedIterator)->getReactance()));
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*exsharedIterator)->getExtNodeRank() << "\t" << -1/((*exsharedIterator)->getReactance()) << "\n";
				++index;
				++rCount; // Increment the row count to point to the next transmission line object
			}
			else {
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*exsharedIterator)->getExtNodeRank()), ar.push_back(1/((*exsharedIterator)->getReactance()));
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*exsharedIterator)->getExtNodeRank() << "\t" << 1/((*exsharedIterator)->getReactance()) << "\n";
				++index;
				ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber + (*exsharedIterator)->getIntlNodeID()), ar.push_back(-1/((*exsharedIterator)->getReactance()));
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber + (*exsharedIterator)->getIntlNodeID() << "\t" << -1/((*exsharedIterator)->getReactance()) << "\n";
				++index;
				++rCount; // Increment the row count to point to the next transmission line object
			}
		}
	}
	/*******************************************************************************************/
	
	// Coefficients corresponding to shared candidate Line Forward Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to shared candidate Line Forward Flow Limit Constraints\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines)), ar.push_back(1);
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines) << "\t" << 1 << "\n";
			++index;
			ia.push_back(rCount), ja.push_back(countOfScenarios*sharedCLines+rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines)-scenCounter*sharedCLines), ar.push_back(-((*candIterator)->getFlowLimit()));
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*sharedCLines+rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines)-scenCounter*sharedCLines << "\t" << -((*candIterator)->getFlowLimit()) << "\n";
			++index;
			++rCount;
		}
	}
	// Coefficients corresponding to shared candidate Line Reverse Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to shared candidate Line Reverse Flow Limit Constraints\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+sharedCLines)), ar.push_back(1);
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+sharedCLines) << "\t" << 1 << "\n";
			++index;
			ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines)-scenCounter*sharedCLines), ar.push_back(((*candIterator)->getFlowLimit()));
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines)-scenCounter*sharedCLines << "\t" << ((*candIterator)->getFlowLimit()) << "\n";
			++index;
			++rCount;
		}
	}
	/*******************************************************************************************/
	
	// Coefficients corresponding to shared candidate Line Definition upper bound
	outPutFile << "\nCoefficients corresponding to shared candidate Line Definition upper bound\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			if ((*candIterator)->getFlowDir() == 1) {
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+2*sharedCLines)), ar.push_back(1);
				++index;
				outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+2*sharedCLines) << "\t" << 1 << "\n";
				ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*candIterator)->getIntlNodeID()), ar.push_back(-1/((*candIterator)->getReactance()));
				++index;
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*candIterator)->getIntlNodeID() << "\t" << -1/((*candIterator)->getReactance()) << "\n";
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*candIterator)->getExtNodeRank()), ar.push_back(1/((*candIterator)->getReactance()));
				++index;
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*candIterator)->getExtNodeRank() << "\t" << 1/((*candIterator)->getReactance()) << "\n";
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+sharedCLines)-scenCounter*sharedCLines), ar.push_back(2.5*((*candIterator)->getFlowLimit()));
				++index;
				outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+sharedCLines)-scenCounter*sharedCLines << "\t" << BIGM << "\n";
				++rCount; // Increment the row count to point to the next transmission line object
			}
			else {
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+2*sharedCLines)), ar.push_back(1);
				++index;
				outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+2*sharedCLines) << "\t" << 1 << "\n";
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*candIterator)->getExtNodeRank()), ar.push_back(-1/((*candIterator)->getReactance()));
				++index;
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*candIterator)->getExtNodeRank() << "\t" << -1/((*candIterator)->getReactance()) << "\n";
				ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*candIterator)->getIntlNodeID()), ar.push_back(1/((*candIterator)->getReactance()));
				++index;
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*candIterator)->getIntlNodeID() << "\t" << 1/((*candIterator)->getReactance()) << "\n";
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+sharedCLines)-scenCounter*sharedCLines), ar.push_back(2.5*((*candIterator)->getFlowLimit()));
				++index;
				outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+sharedCLines)-scenCounter*sharedCLines << "\t" << BIGM << "\n";
				++rCount; // Increment the row count to point to the next transmission line object
			}
		}
	}
	// Coefficients corresponding to shared candidate Line Definition lower bound
	outPutFile << "\nCoefficients corresponding to shared candidate Line Definition lower bound\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			if ((*candIterator)->getFlowDir() == 1) {
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines)), ar.push_back(1);
				++index;
				outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines) << "\t" << 1 << "\n";
				ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*candIterator)->getIntlNodeID()), ar.push_back(-1/((*candIterator)->getReactance()));
				++index;
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*candIterator)->getIntlNodeID() << "\t" << -1/((*candIterator)->getReactance()) << "\n";
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*candIterator)->getExtNodeRank()), ar.push_back(1/((*candIterator)->getReactance()));
				++index;
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*candIterator)->getExtNodeRank() << "\t" << 1/((*candIterator)->getReactance()) << "\n";
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+2*sharedCLines)-scenCounter*sharedCLines), ar.push_back(-2.5*((*candIterator)->getFlowLimit()));
				++index;
				outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+2*sharedCLines)-scenCounter*sharedCLines << "\t" << -BIGM << "\n";
				++rCount; // Increment the row count to point to the next transmission line object
			}
			else {
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines)), ar.push_back(1);
				++index;
				outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines) << "\t" << 1 << "\n";
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*candIterator)->getExtNodeRank()), ar.push_back(-1/((*candIterator)->getReactance()));
				++index;
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*candIterator)->getExtNodeRank() << "\t" << -1/((*candIterator)->getReactance()) << "\n";
				ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*candIterator)->getIntlNodeID()), ar.push_back(1/((*candIterator)->getReactance()));
				++index;
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*candIterator)->getIntlNodeID() << "\t" << 1/((*candIterator)->getReactance()) << "\n";
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+2*sharedCLines)-scenCounter*sharedCLines), ar.push_back(-2.5*((*candIterator)->getFlowLimit()));
				++index;
				outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+2*sharedCLines)-scenCounter*sharedCLines << "\t" << -BIGM << "\n";
				++rCount; // Increment the row count to point to the next transmission line object
			}
		}
	}
	/*******************************************************************************************/

        // Coefficients corresponding to Internal candidate Line Forward Flow Limit Constraints
        outPutFile << "\nCoefficients corresponding to Internal candidate Line Forward Flow Limit Constraints\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
        	for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines)+sharedCLines+countOfScenarios*otherNodeCount), ar.push_back(1);
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines)+sharedCLines+countOfScenarios*otherNodeCount << "\t" << 1 << "\n";
			++index;
			ia.push_back(rCount), ja.push_back(countOfScenarios*internalCLines+rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines)+sharedCLines+countOfScenarios*otherNodeCount-scenCounter*internalCLines), ar.push_back(-((*intCandIterator)->getFlowLimit()));
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*internalCLines+rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines)+sharedCLines+countOfScenarios*otherNodeCount-scenCounter*internalCLines << "\t" << -((*intCandIterator)->getFlowLimit()) << "\n";
			++index;
			++rCount;
        	}
	}
        // Coefficients corresponding to Internal candidate Line Reverse Flow Limit Constraints\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
        	for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines+internalCLines)+sharedCLines+countOfScenarios*otherNodeCount), ar.push_back(1);
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines+internalCLines)+sharedCLines+countOfScenarios*otherNodeCount << "\t" << 1 << "\n";
			++index;
			ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines)+sharedCLines+countOfScenarios*otherNodeCount-scenCounter*internalCLines), ar.push_back(((*intCandIterator)->getFlowLimit()));
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines)+sharedCLines+countOfScenarios*otherNodeCount-scenCounter*internalCLines << "\t" << ((*intCandIterator)->getFlowLimit()) << "\n";
			++index;
			++rCount;
       		}
	}
	/*******************************************************************************************/

        // Coefficients corresponding to Internal candidate Line Definition upper bound
        outPutFile << "\nCoefficients corresponding to Internal candidate Line Definition upper bound\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
        	for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines+2*internalCLines)+sharedCLines+countOfScenarios*otherNodeCount), ar.push_back(1);
			++index;
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines+2*internalCLines)+sharedCLines+countOfScenarios*otherNodeCount << "\t" << 1 << "\n";
			ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*intCandIterator)->getTranslNodeID1()), ar.push_back(-1/((*intCandIterator)->getReactance()));
			++index;
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*intCandIterator)->getTranslNodeID1() << "\t" << -1/((*intCandIterator)->getReactance()) << "\n";
			ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*intCandIterator)->getTranslNodeID2()), ar.push_back(1/((*intCandIterator)->getReactance()));
			++index;
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*intCandIterator)->getTranslNodeID2() << "\t" << 1/((*intCandIterator)->getReactance()) << "\n";
			ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines+internalCLines)+sharedCLines+countOfScenarios*otherNodeCount-scenCounter*internalCLines), ar.push_back(2.5*((*intCandIterator)->getFlowLimit()));
			++index;
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines+internalCLines)+sharedCLines+countOfScenarios*otherNodeCount-scenCounter*internalCLines << "\t" << BIGM << "\n";
			++rCount; // Increment the row count to point to the next transmission line object
        	}
	}
        // Coefficients corresponding to Internal candidate Line Definition lower bound
        outPutFile << "\nCoefficients corresponding to Internal candidate Line Definition lower bound\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
        	for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines+3*internalCLines)+sharedCLines+countOfScenarios*otherNodeCount), ar.push_back(1);
			++index;
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines+3*internalCLines)+sharedCLines+countOfScenarios*otherNodeCount << "\t" << 1 << "\n";
			ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*intCandIterator)->getTranslNodeID1()), ar.push_back(-1/((*intCandIterator)->getReactance()));
			++index;
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*intCandIterator)->getTranslNodeID1() << "\t" << -1/((*intCandIterator)->getReactance()) << "\n";
			ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*intCandIterator)->getTranslNodeID2()), ar.push_back(1/((*intCandIterator)->getReactance()));
			++index;
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*intCandIterator)->getTranslNodeID2() << "\t" << 1/((*intCandIterator)->getReactance()) << "\n";
			ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines+2*internalCLines)+sharedCLines+countOfScenarios*otherNodeCount-scenCounter*internalCLines), ar.push_back(-2.5*((*intCandIterator)->getFlowLimit()));
			++index;
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines+2*internalCLines)+sharedCLines+countOfScenarios*otherNodeCount-scenCounter*internalCLines << "\t" << -BIGM << "\n";
			++rCount; // Increment the row count to point to the next transmission line object
        	}
	}
	/*******************************************************************************************/

	outPutFile << "\nCoefficient Matrix specified" << endl;
	clock_t end1 = clock(); // stop the timer
	double elapsed_secs1 = double(end1 - begin) / CLOCKS_PER_SEC; // Calculate the time required to populate the constraint matrix and objective coefficients
	outPutFile << "\nTotal time taken to define the rows, columns, objective and populate the coefficient matrix = " << elapsed_secs1 << " s " << endl;
	/* RUN THE OPTIMIZATION SIMULATION ALGORITHM */
	cout << "\nSimulation in Progress. Wait !!! ....." << endl;
	glp_load_matrix(milp, --index, &ia[0], &ja[0], &ar[0]); // Loads the Problem
	int lpMethodChoice; // Choice between Simplex or Interior Point Methods to solve the LP relaxation problem
	lpMethodChoice = lpSolveAlgo;
	lpMethodChoice = 1;
	glp_init_iocp(ipControlParam); // Initialize the Mixed Integer Programming Control Parameters with the default setting values
	ipControlParam->mip_gap = 1e-1; // Adjust the tolerance of the Integer Oprtimization part for faster convergence 
	cout << "\nTolerence for MIP is " << ipControlParam->mip_gap;
	switch (lpMethodChoice){
	case 1:
		glp_simplex(milp, NULL);
		break;
	case 2:
		glp_interior(milp, NULL);
		break;
	default:
		cout << "\nInvalid Choice" << endl;
		break;
	}
	string outLPLogFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputSummaryGLPKBounds/OutLPLogBoundGLPK" + to_string(zonalIndex) + ".txt";
	string outMILPLogFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputSummaryGLPKBounds/OutMIPLogBoundGLPK" + to_string(zonalIndex) + ".txt";
	glp_print_sol(milp, outLPLogFileName.c_str()); // write the solution log to the relaxed LP problem
	glp_intopt(milp, ipControlParam); // Calls the GLPK Branch & Bound-Cutting Plane and Simplex Method to solve the integer optimization problem
	glp_print_sol(milp, outMILPLogFileName.c_str()); // write the solution log to the integer optimization problem
	int stat = glp_mip_status(milp); // Outputs the solution status of the problem 

	/* DISPLAY THE SOLUTION DETAILS */
	switch (stat){
	case 1:
		outPutFile << "\nThe solution to the problem is UNDEFINED." << endl;
		cout << "\nThe solution to the problem is UNDEFINED." << endl;
		break;
	case 2:
		outPutFile << "\nThe solution to the problem is FEASIBLE." << endl;
		cout << "\nThe solution to the problem is FEASIBLE." << endl;
		break;
	case 3:
		outPutFile << "\nThe solution to the problem is INFEASIBLE." << endl;
		cout << "\nThe solution to the problem is INFEASIBLE." << endl;
		break;
	case 4:
		outPutFile << "\nNO FEASIBLE solution to the problem exists." << endl;
		cout << "\nNO FEASIBLE solution to the problem exists." << endl;
		break;
	case 5:
		outPutFile << "\nThe solution to the problem is OPTIMAL." << endl;
		cout << "\nThe solution to the problem is OPTIMAL." << endl;
		break;
	case 6:
		outPutFile << "\nThe solution to the problem is UNBOUNDED." << endl;
		cout << "\nThe solution to the problem is UNBOUNDED." << endl;
		break;
	default:
		outPutFile << "\nERROR in Solution Status." << endl;
		cout << "\nERROR in Solution Status." << endl;
		break;
	}

	//Get the Optimal Objective Value results//
	z = glp_mip_obj_val(milp);
	for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
		z += (*genIterator)->getNLCost();
	}

	// Open separate output files for writing results of different variables
	string outIntAngFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputAnglesGLPKBounds/internalAngleBoundGLPK" + to_string(zonalIndex) + ".txt";
	string outCandFlowFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/candidateLinesGLPKBounds/candFlowMWBoundGLPK" + to_string(zonalIndex) + ".txt";
	string outCandDecFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/candidateLinesGLPKBounds/candLineDecisionBoundGLPK" + to_string(zonalIndex) + ".txt";
	string outExtAngFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputAnglesGLPKBounds/externalAngleBoundGLPK" + to_string(zonalIndex) + ".txt";
        string outintCandFlowFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/intcandLinesGLPKBounds/intcandFlowMWBoundGLPK" + to_string(zonalIndex) + ".txt";
        string outintCandLineDecFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/intcandLinesGLPKBounds/intcandLineDecisionBoundGLPK" + to_string(zonalIndex) + ".txt";
        string outTranFlowFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/tranLinesGLPKBounds/tranLineFlowBoundGLPK" + to_string(zonalIndex) + ".txt";
        string outSEFlowFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/SELinesGLPKBounds/SELineFlowGLPK" + to_string(zonalIndex) + ".txt";
	ofstream internalAngleOut(outIntAngFileName, ios::out); //switchStateOut
	ofstream candFlowMWOut(outCandFlowFileName, ios::out); //switchOnOut
	ofstream candLineDecisionOut(outCandDecFileName, ios::out); //switchOffOut
	ofstream externalAngleOut(outExtAngFileName, ios::out);
        ofstream intCandFlowMWOut(outintCandFlowFileName, ios::out);
        ofstream intCandLineDecisionOut(outintCandLineDecFileName, ios::out);
        ofstream tranFlowOut(outTranFlowFileName, ios::out);
        ofstream SEFlowOut(outSEFlowFileName, ios::out);
	outPutFile << "\nThe Optimal Objective value (Generation Dispatch  and Line Building Decision cost) is: " << z << endl;
	powerGenOut << "\nThe Optimal Objective value (Generation Dispatch  and Line Building Decision cost) is: " << z << endl;
	cout << "\nThe Optimal Objective value (Generation Dispatch  and Line Building Decision cost) is: " << z << endl;
	x.push_back(0); // Initialize the decision Variable vector

	//Display Power Generation
	powerGenOut << "\n****************** GENERATORS' POWER GENERATION LEVELS (MW) *********************" << endl;
	powerGenOut << "GENERATOR ID" << "\t" << "GENERATOR MW" << "\n";
	int arrayInd = 1;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
			x.push_back(glp_mip_col_val(milp, arrayInd));
			powerGenOut << (*genIterator)->getGenID() << "\t" << (glp_mip_col_val(milp, arrayInd))*100 << " MW" << endl;
			++arrayInd;
		}
	}
	powerGenOut << "Finished writing Power Generation" << endl;

	// Display Internal node voltage phase angle variables
	internalAngleOut << "\n****************** INTERNAL NODE VOLTAGE PHASE ANGLE VALUES *********************" << endl;
	internalAngleOut << "NODE ID" << "\t" << "VOLTAGE PHASE ANGLE" << "\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); ++nodeIterator){
			if ( (((*nodeIterator)->getSharedFlag()) == 1) || (((*nodeIterator)->getCandFlag()) == 1) ) {
				x.push_back(glp_mip_col_val(milp, arrayInd));
				internalAngleOut << (*nodeIterator)->getNodeID() << "\t" << glp_mip_col_val(milp, arrayInd) << endl;	
				++arrayInd;
			}
			else {
				x.push_back(glp_mip_col_val(milp, arrayInd));
				internalAngleOut << (*nodeIterator)->getNodeID() << "\t" << glp_mip_col_val(milp, arrayInd) << endl;		
				++arrayInd;	
			}		
		}
	}
	internalAngleOut << "Finished writing Internal Node Voltage Phase Angles" << endl;

	// Display Shared Candidate lines' Power Flow variables
	candFlowMWOut << "\n****************** SHARED CANDIDATE LINES POWER FLOW VALUES *********************" << endl;
	candFlowMWOut << "CANDIDATE LINE ID" << "\t" << "MW FLOW VALUE" << "\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			x.push_back(glp_mip_col_val(milp, arrayInd));
			candFlowMWOut << (*candIterator)->getTranslID() << "\t" << (glp_mip_col_val(milp, arrayInd))*100 << " MW" << endl;
			++arrayInd;
		}
	}
	candFlowMWOut << "Finished writing Shared Candidate lines' Power Flow variables" << endl;

	// Display Shared Candidate lines' Construction Decisions
	candLineDecisionOut << "\n****************** SHARED CANDIDATE LINES CONSTRUCTION DECISIONS *********************" << endl;
	candLineDecisionOut << "CANDIDATE LINE ID" << "\t" << "CONSTRUCTION DECISIONS" << "\n";
	for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
		x.push_back(glp_mip_col_val(milp, arrayInd));
		candLineDecisionOut << (*candIterator)->getTranslID() << "\t" << (glp_mip_col_val(milp, arrayInd)) << endl;
		++arrayInd;
	}
	candLineDecisionOut << "Finished writing Shared Candidate lines' Construction decisions" << endl;

	// Display Outer Zonal node angles
	externalAngleOut << "\n****************** OUTER ZONAL NODE VOLTAGE PHASE ANGLE VALUES *********************" << endl;
	externalAngleOut << "EXTERNAL NODE ID" << "\t" << "VOLTAGE PHASE ANGLE" << "\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int diffNodeC = 0;
		for (globalIterator = globalRankDiffNode.begin(); globalIterator != globalRankDiffNode.end(); ++globalIterator){
			if (diffNodeC > 0) { // Skip the dummy element "0" at the beginning
				x.push_back(glp_mip_col_val(milp, arrayInd));
				externalAngleOut << (*globalIterator) << "\t" << (glp_mip_col_val(milp, arrayInd)) << endl;
				++arrayInd;
			}
			++diffNodeC;
		}
	}
	externalAngleOut << "Finished writing outer zonal node voltage phase angle values" << endl;

        // Display Internal Candidate lines' Power Flow variables
        intCandFlowMWOut << "\n****************** INTERNAL CANDIDATE LINES POWER FLOW VALUES *********************" << endl;
        intCandFlowMWOut << "CANDIDATE LINE ID" << "\t" << "MW FLOW VALUE" << "\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
        	for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			x.push_back(glp_mip_col_val(milp, arrayInd));
			intCandFlowMWOut << (*intCandIterator)->getTranslID() << "\t" << (glp_mip_col_val(milp, arrayInd))*100 << " MW" << endl;
			++arrayInd;
        	}
	}
        intCandFlowMWOut << "Finished writing Internal Candidate lines' Power Flow variables" << endl;

        // Display Internal Candidate lines' Construction Decisions
        intCandLineDecisionOut << "\n****************** INTERNAL CANDIDATE LINES CONSTRUCTION DECISIONS *********************" << endl;
        intCandLineDecisionOut << "CANDIDATE LINE ID" << "\t" << "CONSTRUCTION DECISIONS" << "\n";
        for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
		x.push_back(glp_mip_col_val(milp, arrayInd));
		intCandLineDecisionOut << (*intCandIterator)->getTranslID() << "\t" << (glp_mip_col_val(milp, arrayInd)) << endl;
                ++arrayInd;
        }
        intCandLineDecisionOut << "Finished writing Internal Candidate lines' Construction decisions" << endl;
        // Display Internal Transmission lines' Flows
        tranFlowOut << "\n****************** INTERNAL TRANSMISSION LINES FLOWS *********************" << endl;
        tranFlowOut << "SCENARIO ID" << "\t" << "TRANSMISSION LINE ID" << "\t" << "MW FLOW" << "\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
                for (tranIterator = translObject.begin(); tranIterator != translObject.end(); ++tranIterator){
                        tranFlowOut << scenCounter << "\t" << (*tranIterator)->getTranslID() << "\t" << (1/((*tranIterator)->getReactance()))*((glp_mip_col_val(milp, (countOfScenarios*genNumber + scenCounter*nodeNumber+(*tranIterator)->getTranslNodeID1())))-(glp_mip_col_val(milp, (countOfScenarios*genNumber + scenCounter*nodeNumber+(*tranIterator)->getTranslNodeID2()))))*100 << " MW" << endl;
                }
        }
        tranFlowOut << "Finished writing Internal Transmission lines' MW Flows" << endl;

        // Display SE lines' Flows
        SEFlowOut << "\n****************** SHARED EXISTING TRANSMISSION LINES FLOWS *********************" << endl;
        SEFlowOut << "SCENARIO ID" << "\t" << "SHARED EXISTING TRANSMISSION LINE ID" << "\t" << "MW FLOW" << "\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
                for (exsharedIterator = SELineObject.begin(); exsharedIterator != SELineObject.end(); ++exsharedIterator){
                        if ((*exsharedIterator)->getFlowDir() == 1) {
                                SEFlowOut << scenCounter << "\t" << (*exsharedIterator)->getTranslID() << "\t" << (1/((*exsharedIterator)->getReactance()))*((glp_mip_col_val(milp, (countOfScenarios*genNumber + scenCounter*nodeNumber+(*exsharedIterator)->getIntlNodeID())))-(glp_mip_col_val(milp, (countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount+(*exsharedIterator)->getExtNodeRank()))))*100 << " MW" << endl;
                        }
                        else {
                                SEFlowOut << scenCounter << "\t" << (*exsharedIterator)->getTranslID() << "\t" << (1/((*exsharedIterator)->getReactance()))*((glp_mip_col_val(milp, (countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount+(*exsharedIterator)->getExtNodeRank())))-(glp_mip_col_val(milp, (countOfScenarios*genNumber + scenCounter*nodeNumber+(*exsharedIterator)->getIntlNodeID()))))*100 << " MW" << endl;
                        }
                }
        }
        tranFlowOut << "Finished writing Shared Existing Transmission lines' MW Flows" << endl;
	delete ipControlParam; // free the memory of the Integer Programming Control Parameter struct
	glp_delete_prob(milp); // Free the memory of the GLPK Problem Object
	clock_t end2 = clock(); // stop the timer
	double elapsed_secs2 = double(end2 - begin) / CLOCKS_PER_SEC; // Calculate the Total Time
	outPutFile << "\nTotal time taken to solve the MILP Line Construction Decision Making Problem instance and retrieve the results = " << elapsed_secs2 << " s " << endl;
	cout << "\nTotal time taken to solve the MILP Line Construction Decision Making Problem instance and retrieve the results = " << elapsed_secs2 << " s " << endl;

	// Close the different output files
	outPutFile.close();
	powerGenOut.close();
	internalAngleOut.close();
	candFlowMWOut.close();
	candLineDecisionOut.close();
	externalAngleOut.close();
        intCandFlowMWOut.close();
        intCandLineDecisionOut.close();
        tranFlowOut.close();
        SEFlowOut.close();
	cout << "\nSimulation Completed.\nResults written on the different output files" << endl;
	return z;
} // Function MILP() ends

int Nettran::getConnZone(int i) // returns the pointer to the base of the vector, diffZoneID
{
	if (otherZoneIter != diffZoneID.end()) { // check to see if the end of the diffZoneID vector has been reached
		otherZoneIter++; // if not, increment otherZoneIter iterator for the next call
		return diffZoneID.at(i); // return the zone ID present at the position i
	}
	else
		return -1; // if end is reached, return -1
} // Function getConnZone() ends

int Nettran::getConnNode(int i) // returns the pointer to the base of the vector, diffZoneNodeID
{
	if (otherZoneNodeIter != diffZoneNodeID.end()) { // check to see if the end of the diffZoneNodeID vector has been reached
		otherZoneNodeIter++; // if not, increment otherZoneNodeIter iterator for the next call
		return diffZoneNodeID.at(i); // return the node ID present at the position i
	}
	else
		return -1; // if end is reached, return -1
} // Function getConnNode() ends

int Nettran::getSESerial(int i) // returns the pointer to the base of the vector, SELineObject
{
	if (sharedELineIt != SELineObject.end()) { // check to see if the end of the SELineObject vector has been reached
		sharedELineIt++; // if not, increment sharedELineIt iterator for the next call
		return (SELineObject.at(i))->getTranslID(); // return the global serial number of the present SE line at the position i
	}
	else
		return -1; // if end is reached, return -1
} // Function getSESerial() ends

int Nettran::getSEFromNode(int i) // returns the ID number of the from node of the shared existing line
{
	if (sharedELineFromIt != SELineObject.end()) { // check to see if the end of the SELineObject vector has been reached
		sharedELineFromIt++; // if not, increment iterator for the next call
		if ((SELineObject.at(i))->getFlowDir()==1) // If the intrazonal node is the from node
			return (SELineObject.at(i))->getIntlNodeID(); // return the ID of the internal node of the SE line at the position i
		else
			return (SELineObject.at(i))->getExtNodeID(); // else, return the ID of the external node of the SE line at the position i
	}
	else
		return -1; // if end is reached, return -1 
} // Function getSEFromNode() ends

int Nettran::getSEFromZone(int i) // returns the ID number of the from zone of the shared existing line
{
	if (sharedELineFZoneIt != SELineObject.end()) { // check to see if the end of the SELineObject vector has been reached
		sharedELineFZoneIt++; // if not, increment iterator for the next call
		if ((SELineObject.at(i))->getFlowDir()==1) // If the intrazonal node is the from node
			return (SELineObject.at(i))->getIntlZoneID(); // return the ID of current zone of the SE line at the position i
		else
			return (SELineObject.at(i))->getExtZoneID(); // else, return the ID of the external zone of the SE line at the position i
	}
	else
		return -1; // if end is reached, return -1 
} // Function getSEFromZone() ends

int Nettran::getSEToNode(int i) // returns the ID number of the to node of the shared existing line
{
	if (sharedELineToIt != SELineObject.end()) { // check to see if the end of the SELineObject vector has been reached
		sharedELineToIt++; // if not, increment iterator for the next call
		if ((SELineObject.at(i))->getFlowDir()==1) // If the intrazonal node is the from node
			return (SELineObject.at(i))->getExtNodeID(); // return the ID of the external node of the SE line at the position i
		else
			return (SELineObject.at(i))->getIntlNodeID(); // else, return the ID of the internal node of the SE line at the position i
	}
	else
		return -1; // if end is reached, return -1 
} // Function getSEToNode() ends

int Nettran::getSEToZone(int i) // returns the ID number of the to zone of the shared existing line
{
	if (sharedELineTZoneIt != SELineObject.end()) { // check to see if the end of the SELineObject vector has been reached
		sharedELineTZoneIt++; // if not, increment iterator for the next call
		if ((SELineObject.at(i))->getFlowDir()==1) // If the intrazonal node is the from node
			return (SELineObject.at(i))->getExtZoneID(); // return the ID of external zone of the SE line at the position i
		else
			return (SELineObject.at(i))->getIntlZoneID(); // else, return the ID of the current zone of the SE line at the position i
	}
	else
		return -1; // if end is reached, return -1 
} // Function getSEToZone() ends

double Nettran::getSEReactance(int i) // returns the reactance of the shared existing line
{
	if (sharedELineReactIt != SELineObject.end()) { // check to see if the end of the SELineObject vector has been reached
		sharedELineReactIt++; // if not, increment iterator for the next call
		return (SELineObject.at(i))->getReactance(); // return the reactance of the SE line at the position i
	}
	else
		return -1; // if end is reached, return -1 
} // Function getConnNode() ends

double Nettran::getSECapacity(int i) // returns the capacity of the shared existing line
{
	if (sharedELineCapIt != SELineObject.end()) { // check to see if the end of the SELineObject vector has been reached
		sharedELineCapIt++; // if not, increment iterator for the next call
		return (SELineObject.at(i))->getFlowLimit(); // return the capacity of the SE line at the position i
	}
	else
		return -1; // if end is reached, return -1 
} // Function getConnNode() ends

int Nettran::getCandSerial(int i) // returns the pointer to the base of the vector, SELineObject
{
	if (sharedCandLineIt != candLineObject.end()) { // check to see if the end of the SELineObject vector has been reached
		sharedCandLineIt++; // if not, increment sharedELineIt iterator for the next call
		return (candLineObject.at(i))->getTranslID(); // return the global serial number of the present SE line at the position i
	}
	else
		return -1; // if end is reached, return -1
} // Function getSESerial() ends

int Nettran::getCandFromNode(int i) // returns the ID number of the from node of the shared existing line
{
	if (sharedCandLineFromIt != candLineObject.end()) { // check to see if the end of the SELineObject vector has been reached
		sharedCandLineFromIt++; // if not, increment iterator for the next call
		if ((candLineObject.at(i))->getFlowDir()==1) // If the intrazonal node is the from node
			return (candLineObject.at(i))->getIntlNodeID(); // return the ID of the internal node of the SE line at the position i
		else
			return (candLineObject.at(i))->getExtNodeID(); // else, return the ID of the external node of the SE line at the position i
	}
	else
		return -1; // if end is reached, return -1 
} // Function getSEFromNode() ends

int Nettran::getCandFromZone(int i) // returns the ID number of the from zone of the shared existing line
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

int Nettran::getCandToNode(int i) // returns the ID number of the to node of the shared existing line
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

int Nettran::getCandToZone(int i) // returns the ID number of the to zone of the shared existing line
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

double Nettran::getCandReactance(int i) // returns the reactance of the shared existing line
{
	if (sharedCandLineReactIt != candLineObject.end()) { // check to see if the end of the SELineObject vector has been reached
		sharedCandLineReactIt++; // if not, increment iterator for the next call
		return (candLineObject.at(i))->getReactance(); // return the reactance of the SE line at the position i
	}
	else
		return -1; // if end is reached, return -1 
} // Function getConnNode() ends

double Nettran::getCandCapacity(int i) // returns the capacity of the shared existing line
{
	if (sharedCandLineCapIt != candLineObject.end()) { // check to see if the end of the SELineObject vector has been reached
		sharedCandLineCapIt++; // if not, increment iterator for the next call
		return (candLineObject.at(i))->getFlowLimit(); // return the capacity of the SE line at the position i
	}
	else
		return -1; // if end is reached, return -1 
} // Function getConnNode() ends

void Nettran::setSEFromRank(int existingCounter, int rankFrom) // sets the from rank for the internal zone node end of the shared existing line
{
	 (*(SELineObject.at(existingCounter))).assignRank(rankFrom);
}

void Nettran::setSEToRank(int existingCounter, int rankTo) // sets the to rank for the internal zone node end of the shared existing line
{
	(*(SELineObject.at(existingCounter))).assignRank(rankTo);
}

void Nettran::setCandFromRank(int candidateCounter, int rankFrom) // sets the from rank for the internal zone node end of the shared candidate line
{
	(*(candLineObject.at(candidateCounter))).assignRank(rankFrom);
}

void Nettran::setCandToRank(int candidateCounter, int rankTo) // sets the to rank for the internal zone node end of the shared candidate line
{
	(*(candLineObject.at(candidateCounter))).assignRank(rankTo);
}

void Nettran::setSEFromRankConn(int existingCounter, int rankFrom) // populates the vector of the global ranks of all the external from nodes of SE lines
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

void Nettran::setSEToRankConn(int existingCounter, int rankTo) // populates the vector of the global ranks of all the external to nodes of SE lines
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

void Nettran::setCandFromRankConn(int candidateCounter, int rankFrom) // populates the vector of the global ranks of all the external from nodes of cand lines
{
	for (globalIterator = globalRankDiffNode.begin(); globalIterator != globalRankDiffNode.end(); ++globalIterator) {
		if (*globalIterator == rankFrom) { // Check whether the other-zone node is already present in the list
			containsFlag = 1;  // Loop through the diffZoneNodeID and diffZoneID vectors, to see if the other zone-node is already present or not, if yes, set the containsFlag to 1
		}
	}
	if (containsFlag == 0) { // If globalRankDiffNode vector does not contain the other zone-node rank, then push the node rank in the vector
		globalRankDiffNode.push_back(rankFrom);
	}
	containsFlag = 0; // reset
	(*(candLineObject.at(candidateCounter))).connectRank(rankFrom);
}

void Nettran::setCandToRankConn(int candidateCounter, int rankTo) // populates the vector of the global ranks of all the external to nodes of cand lines
{
	for (globalIterator = globalRankDiffNode.begin(); globalIterator != globalRankDiffNode.end(); ++globalIterator) {
		if (*globalIterator == rankTo) { // Check whether the other-zone node is already present in the list
			containsFlag = 1;  // Loop through the diffZoneNodeID and diffZoneID vectors, to see if the other zone-node is already present or not, if yes, set the containsFlag to 1
		}
	}
	if (containsFlag == 0) { // If globalRankDiffNode vector does not contain the other zone-node rank, then push the node rank in the vector
		globalRankDiffNode.push_back(rankTo);
	}
	containsFlag = 0; // reset
	(*(candLineObject.at(candidateCounter))).connectRank(rankTo);
}

void Nettran::assignCandGlobalRank(int candidateCounter, int candGlobalRank) // assigns the global rank of the shared candidate line
{
	(*(candLineObject.at(candidateCounter))).assignLineRank(candGlobalRank);
}

void Nettran::setRealizedCLines(Marketover &coordInstanceRef) // Assigns the realizedCLine variable, the number of candidate lines that are actually built
{
	// *** Built shared candidate lines *** //
	vector<candLine*>::iterator candIterator; // Iterator for candidate lines
	for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator) {
		int candGlobRank=((*candIterator)->getGlobalRank()); // Get the global rank of the candidate line 		
		int statusPresAbs = coordInstanceRef.scanBuiltLinesList(candGlobRank); // Passing on the shared node angle decision message to the MO
		if (statusPresAbs==1)
			(*candIterator)->setPresAbsStatus();
		if ((*candIterator)->returnPresAbsStatus()==1) {
			realCandLine.push_back((*candIterator)); // Populate the vector of realized candidate lines
			//cout << "Constructed shared line ID: " << (*candIterator)->getTranslID() << endl;
			(*candIterator)->modifyNodeReact(); // Adjust the connected nodal reactances after the lines have been decided to be built
			int nodeZone2 = (*candIterator)->getExtNodeID();
			int tNodeID2 = (*candIterator)->getExtZoneID();
			int indCount = 0; // Initialize a counter for tracking the position in the vector of the iterator
			int rankBuilt = (*candIterator)->getExtNodeGlobalRank();
			vector<int>::iterator diffZNIt;
			for (diffZNIt = diffZoneNodeExistingID.begin(); diffZNIt != diffZoneNodeExistingID.end(); ++diffZNIt) {
				if ((*diffZNIt == tNodeID2) && (diffZoneExistingID[indCount] == nodeZone2)) { // Check whether the other-zone node is already present in the list
					containsFlag = 1;  // Loop through the diffZoneNodeID and diffZoneID vectors, to see if the other zone-node is already present or not, if yes, set the containsFlag to 1
				}
				++indCount; // Increment the counter
			}
			if (containsFlag == 0) { // If the diffZoneNodeID and diffZoneID vectors do not contain the other zone-node, then push the node and the other zone index in the respective vectors
				diffZoneNodeExistingID.push_back(tNodeID2); // initialize the list of external-zone existing node ID's so that it's not empty
				//cout << " The node ID of outer zonal node is : " << tNodeID2 << endl;
				diffZoneExistingID.push_back(nodeZone2); // Initialize the list of external-zone existing zone ID's so that it's not empty.
				//cout << " The zone ID of outer zonal node is : " << nodeZone2 << endl;
			}
			for (globalIterator = globalRankDiffNode.begin(); globalIterator != globalRankDiffNode.end(); ++globalIterator) {
				if (*globalIterator == rankBuilt) { // Check whether the other-zone node is already present in the list
					containsFlagGlob = 1;  // Loop through the diffZoneNodeID and diffZoneID vectors, to see if the other zone-node is already present or not, if yes, set the containsFlag to 1
				}
			}
			if (containsFlagGlob == 0) { // If globalRankDiffNode vector does not contain the other zone-node rank, then push the node rank in the vector
				globalExistingRank.push_back(rankBuilt); // list of global ranking of external zone nodes for only SE lines and shared candidate lines
				//cout << " The rank of outer zonal node is : " << rankBuilt << endl;
				++existingOtherZoneNodeNum; // Increment the number of other zone nodes that are ends of SE lines, only when a new one is added
				//cout << " The number of existing other zone node is : " << existingOtherZoneNodeNum << endl;
			}
			containsFlag = 0; // reset
			containsFlagGlob = 0; // reset
			//(*(SELineObject.at(existingCounter))).connectRank(rankFrom);	
			++realizedCLines;
		}
	}
	//cout << "Number of constructed shared lines in network " << zonalIndex << " is " << realizedCLines << endl; 

	// *** Built internal candidate lines *** //
	vector<intCandLine*>::iterator intCandIterator; // Iterator for candidate lines
	for (intCandIterator = realIntCandLineObject.begin(); intCandIterator != realIntCandLineObject.end(); ++intCandIterator) {
		(*intCandIterator)->setPresAbsStatus();
		//cout << "Constructed internal line ID: " << (*intCandIterator)->getTranslID() << endl;
		++realizedIntCLines;
	}
	for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator) {
		if (((*intCandIterator)->returnPresAbsStatus()==1) && (find(realIntCandLineObject.begin(), realIntCandLineObject.end(), (*intCandIterator))==realIntCandLineObject.end())) {
			realIntCandLineObject.push_back((*intCandIterator)); // Populate the vector of realized candidate lines
			//cout << "Constructed internal line ID: " << (*intCandIterator)->getTranslID() << endl;
			(*intCandIterator)->modifyNodeReact(); // Adjust the connected nodal reactances after the lines have been decided to be built
			++realizedIntCLines;
		}
	}
	//cout << "Number of constructed internal lines in network " << zonalIndex << " is " << realizedIntCLines << endl; 
}

void Nettran::TestBuiltExternalNodes()
{
	//cout << " Subnetwork # : " << zonalIndex << endl;
	int indCount = 0;
	vector<int>::iterator diffZNIt;
	int cancelFirstEntry = 0; 
	for (diffZNIt = diffZoneNodeExistingID.begin(); diffZNIt != diffZoneNodeExistingID.end(); ++diffZNIt) {
		if (cancelFirstEntry > 0) {
			//cout << "\nOuter node number : " << *diffZNIt << " Outer node zone : " << diffZoneExistingID[indCount] << " Outer node global rank : " << globalExistingRank[indCount] << endl;
		}
		++indCount; // Increment the counter
		++cancelFirstEntry;
	}
	//cout << " Total number of outer zone nodes for this subnetwork to which existing or built lines are connected : " << existingOtherZoneNodeNum << endl;
}

int Nettran::returnMultiplicity() // Returns the total multiplicity of the shared nodes
{
	vector<Node*>::iterator nodeIterator; // Iterator for node objects
	int connMult = 0;
	for (nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); ++nodeIterator){
		connMult += (*nodeIterator)->getNodeMultiplicity();
	}
	return connMult;
}

double Nettran::APPQPAvgHR(Marketover &coordInstanceRef, double LagMultXi[], int totalCandLineNum, int totalSharedNodeNum, GRBEnv* environmentGUROBI, int iterCount) // Calls the GUROBI solver object and solver method to solve the problem of determining the values of the continuous variables
{
	// CREATION OF THE QP SOLVER INSTANCE //
	double gammaPathLength, betaPathLength;
	clock_t begin = clock(); // start the timer
	vector<int>::iterator diffZNIt; // Iterator for diffZoneNodeID
	vector<Powergenerator*>::iterator genIterator; // Iterator for Powergenerator objects
	vector<transmissionLine*>::iterator tranIterator; // Iterator for Transmission line objects
	vector<Load*>::iterator loadIterator; // Iterator for load objects
	vector<Node*>::iterator nodeIterator; // Iterator for node objects
	vector<candLine*>::iterator candIterator; // Iterator for candidate lines
	vector<SELine*>::iterator exsharedIterator; // Iterator for shared existing lines
	vector<intCandLine*>::iterator intCandIterator; // Iterator for internal candidate lines
	string outSummaryFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outSummaryAPP/OutSumAPPGUROBI" + to_string(zonalIndex) + ".txt";
	ofstream outPutFile(outSummaryFileName, ios::out); // Create Output File to output the Summary of Results
	if (!outPutFile){
		cerr << "\nCouldn't open the file" << endl;
		exit(1);
	}

	int dimRow = (2 * genNumber + 2 * realizedCLines + 2 * sharedELines + 2*realizedIntCLines + 2 * tranNumber + nodeNumber); // Total number of rows of the A matrix (number of structural constraints of the QP) first term to account for lower and upper generating limits, second term for lower and upper line limits of candidate shared lines, third term for lower and upper line limits for shared existing lines, fourth term for lower and upper line limits for internal candidate lines, fifth term for lower and upper line limits for internal lines, and the sixth term to account for nodal power balance constraints
	int dimCol = genNumber+nodeNumber+otherNodeCount; // Total number of columns of the LP (number of Decision Variables) first term to account for power generation MW outputs, second term for voltage phase angles for internal zonal nodes, and third term for the voltage phase angles of other-zone nodes connected through shared existing and candidate lines
	vector<double> thetaDiff;
	vector<double> thetaPrev;
	if (iterCount==0) {
		for (nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); ++nodeIterator){
			if ( (((*nodeIterator)->getSharedFlag()) == 1) || (((*nodeIterator)->getBuiltCandFlag()) == 1) ) {
				thetaPrev.push_back(0);
				thetaDiff.push_back(0);
				globRankBuffer.push_back((*nodeIterator)->getGlobalRank());
			}	
		}
		int diffNodeCount = 0;
		for (globalIterator = globalRankDiffNode.begin(); globalIterator != globalRankDiffNode.end(); ++globalIterator){
			if (diffNodeCount > 0) { // Skip the dummy element "0" at the beginning
				thetaPrev.push_back(0);
				thetaDiff.push_back(0);
				globRankBuffer.push_back((*globalIterator));
			}
			++diffNodeCount;
		}
	}
	else {
		thetaPrev=thetaBuffer;
		thetaDiff=diffBuffer;
		thetaBuffer.clear();
		diffBuffer.clear();
		globRankBuffer.clear();		
	}	
	outPutFile << "\nTotal Number of Structural Constraints (Rows) is: " << dimRow << endl;
	outPutFile << "\nTotal Number of Decision Variables (Columns) is: " << dimCol << endl;
	// Instantiate GUROBI Problem model
	GRBModel *modelSubnetQP = new GRBModel(*environmentGUROBI);
	cout << "\nGurobi model created" << endl;
    	modelSubnetQP->set(GRB_StringAttr_ModelName, "subQP" + to_string(zonalIndex));
	cout << "\nGurobi model created and name set" << endl;
	GRBVar decvar[dimCol+1];
	cout << "\nGurobi decision variables created" << endl;
	double z; // variable to store the objective value

	// SPECIFICATION OF PROBLEM PARAMETERS //
	// Dummy Decision Variable //
	cout << "\nGurobi decision variables to be assigned" << endl;
	decvar[0] = modelSubnetQP->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
	//Decision Variable Definitions, Bounds, and Objective Function Co-efficients//
	cout << "\nGurobi dummy decision variable created" << endl;
	int colCount = 1;
	//Columns corresponding to Power Generation continuous variables for different generators//
	outPutFile << "\nCoefficients of Power generator variables" << endl;
	for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
		decvar[colCount] = modelSubnetQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		outPutFile << colCount << "\t";
		outPutFile << ((*genIterator)->getLinCoeff())/100 << " $/MW" << endl;
		++colCount;
	}
	outPutFile << "\nTotal number of columns after accounting for Power Generation continuous variables for different generators: " << colCount << endl;

	//Columns corresponding to Voltage Phase Angles continuous variables for different intrazonal nodes//
	outPutFile << "\nCoefficients of Voltage Phase Angles continuous variables for different intrazonal nodes" << endl;	
	for (nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); ++nodeIterator){
		decvar[colCount] = modelSubnetQP->addVar(0, (44/7), 0.0, GRB_CONTINUOUS);
		++colCount;
	}
	outPutFile << "\nTotal number of columns after accounting for Voltage Phase Angles continuous variables for different intrazonal nodes: " << colCount << endl;

	//Columns corresponding to Voltage Phase Angles continuous variables for other zone nodes for shared lines//
	outPutFile << "\nCoefficients corresponding to Voltage Phase Angles continuous variables for other zone nodes for shared lines" << endl;
	int diffNodeCounter = 0; // counter flag to indicate the first element of the diffZoneNodeID list
	for (globalIterator = globalRankDiffNode.begin(); globalIterator != globalRankDiffNode.end(); ++globalIterator){
		if (diffNodeCounter > 0) { // Skip the first element, since it's a dummy "0"
			decvar[colCount] = modelSubnetQP->addVar(0, (44/7), 0.0, GRB_CONTINUOUS);	
			++colCount;
			//cout << " Column count after the shared node " << diffNodeCounter << " is " << colCount << endl;
		}
		++diffNodeCounter;
	}

	outPutFile << "\nTotal number of columns after accounting for Voltage Phase Angles continuous variables for other zone nodes for shared lines: " << colCount << endl;
	outPutFile << "\nTotal Number of columns for generation, angles, integer variables, and flows: " << colCount-1 << endl;
	outPutFile << "\nDecision Variables and Objective Function defined" << endl;
	outPutFile << "\nTotal Number of columns: " << colCount-1 << endl;
	//Setting Objective//
	GRBQuadExpr obj = 0.0;
	// Objective Contribution from Dummy Decision Variable //
	obj += 0*(decvar[0]);
	colCount = 1;
	//Columns corresponding to Power Generation continuous variables for different generators//
	for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
		obj += ((*genIterator)->getLinCoeff())*(decvar[colCount]);
		++colCount;
	}

	//Columns corresponding to Voltage Phase Angles continuous variables for different intrazonal nodes//	
	for (nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); ++nodeIterator){
		if ( (((*nodeIterator)->getSharedFlag()) == 1) || (((*nodeIterator)->getCandFlag()) == 1) ) {
			obj += (LagMultXi[(zonalIndex-1)*totalSharedNodeNum+((*nodeIterator)->getGlobalRank())]+gammaPathLength*thetaDiff[(zonalIndex-1)*totalSharedNodeNum+((*nodeIterator)->getGlobalRank())])*(decvar[colCount]);
			obj += (betaPathLength/2)*(decvar[colCount]-thetaPrev[colCount])*(decvar[colCount]-thetaPrev[colCount]);				
		}
		else {
			obj += 0*(decvar[colCount]);	
		}
		++colCount;
	}

	//Columns corresponding to Voltage Phase Angles continuous variables for other zone nodes for shared lines//
	diffNodeCounter = 0; // counter flag to indicate the first element of the diffZoneNodeID list
	for (globalIterator = globalRankDiffNode.begin(); globalIterator != globalRankDiffNode.end(); ++globalIterator){
		if (diffNodeCounter > 0) { // Skip the first element, since it's a dummy "0"
			obj += (LagMultXi[(zonalIndex-1)*totalSharedNodeNum+(*globalIterator)]+gammaPathLength*thetaDiff[(zonalIndex-1)*totalSharedNodeNum+(*globalIterator)])*(decvar[colCount]);	
			obj += (betaPathLength/2)*(decvar[colCount]-thetaPrev[colCount])*(decvar[colCount]-thetaPrev[colCount]);	
			++colCount;
		}
		++diffNodeCounter;
	}

	modelSubnetQP->setObjective(obj, GRB_MINIMIZE);
	cout << " Objective Function and Decision Variables have been defined and the colCount is " << colCount-1 << endl;
	//Row Definitions: Specification of b<=Ax<=b//
	GRBLinExpr lhs[dimRow+1];
	//Row Definitions and Bounds Corresponding to Constraints/


	// Constraints corresponding to supply-demand balance
	string outPGenFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outPowerAPP/OutPGenAPPGUROBI" + to_string(zonalIndex) + ".txt"; 
	ofstream powerGenOut(outPGenFileName, ios::out);
	if (!powerGenOut){
		cerr << "\nCouldn't open the file" << endl;

		exit(1);
	}
	//Non-Zero entries of A matrix (Constraint/Coefficient matrix entries)//
	// Coefficients for the supply-demand balance constraints
	outPutFile << "\nNon-zero elements of A matrix" << endl;
	outPutFile << "\nRow Number\tColumn Number\tNon-zero Entry\tFrom Reactance\tToReactance" << endl;
	outPutFile << "\nCoefficients for the supply-demand balance constraints" << endl;
	// Dummy Constraint //
	lhs[0] = 0*(decvar[0]);
	modelSubnetQP->addConstr(lhs[0], GRB_EQUAL, 0);
	int rCount = 1; // Initialize the row count
	vector<int> busCount; // vector for storing the node/bus serial
	outPutFile << "Constraints corresponding to Supply-Demand Balance right hand side" << endl;
	for (nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); ++nodeIterator){
		outPutFile << "\nGeneration\t" << rCount << "\n";
		int genListLength = (*nodeIterator)->getGenLength(); // get the number
		lhs[rCount]=0;
		for (int cCount = 1; cCount <= genListLength; ++cCount){
			lhs[rCount] += 1*(decvar[((*nodeIterator)->getGenSer(cCount))]);
			outPutFile << "\n" << rCount << "\t" << (*nodeIterator)->getGenSer(cCount) << "\t" << 1.0 << endl;
		}
		outPutFile << "\nIntrazonal Node Angles\t" << rCount << "\n";
		lhs[rCount] += (((*nodeIterator)->getToReact())-((*nodeIterator)->getFromReact()))*(decvar[genNumber+rCount]);
		outPutFile << "\n" << rCount << "\t" << genNumber+rCount << "\t" << -((*nodeIterator)->getToReact())-((*nodeIterator)->getFromReact()) << "\t" << -((*nodeIterator)->getFromReact()) << "\t" << -((*nodeIterator)->getToReact()) << endl;
		outPutFile << "\nConnected Intrazonal Node Angles\t" << rCount << "\n";
		int connNodeListLength = (*nodeIterator)->getConNodeLength(); // get the number of intra-zonal nodes connected to this node
		for (int cCount = 1; cCount <= connNodeListLength; ++cCount){
			if (((*nodeIterator)->getConnReact(cCount))<=0)
				lhs[rCount] -= (((*nodeIterator)->getConnReact(cCount)))*(decvar[genNumber+((*nodeIterator)->getConnSer(cCount))]);
			else
				lhs[rCount] += (((*nodeIterator)->getConnReact(cCount)))*(decvar[genNumber+((*nodeIterator)->getConnSer(cCount))]);
			outPutFile << "\n" << rCount << "\t" << genNumber+((*nodeIterator)->getConnSer(cCount)) << "\t" <<  (-((*nodeIterator)->getConnReact(cCount))) << "\n";

		}
		outPutFile << "\nConnected Outer zonal Node Angles\t" << rCount << "\n";
		int connOutNodeLength = (*nodeIterator)->getExtraNodeLength(); // get the number of outer-zonal nodes connected to this node
		for (int cCount = 1; cCount <= connOutNodeLength; ++cCount){
			lhs[rCount] += (-((*nodeIterator)->getExtConnReact(cCount)))*(decvar[genNumber+nodeNumber+(*nodeIterator)->getExtConnSer(cCount)]);
			outPutFile << "\n" << rCount << "\t" << genNumber+nodeNumber+(*nodeIterator)->getExtConnSer(cCount) << "\t" <<  (-((*nodeIterator)->getExtConnReact(cCount))) << "\n";
		}
		busCount.push_back(rCount);
		if (((nodeObject[rCount-1])->devpinitMessage(0))==0)
			modelSubnetQP->addConstr(lhs[rCount], GRB_EQUAL, ((nodeObject[rCount-1])->devpinitMessage(0)));
		else 
			modelSubnetQP->addConstr(lhs[rCount], GRB_EQUAL, -((nodeObject[rCount-1])->devpinitMessage(0)));
		outPutFile << "Connected load to node " << rCount << " is " << (nodeObject[rCount-1])->devpinitMessage(0)*100 << " MW" << endl;
		outPutFile << rCount << "\t";
		if (((nodeObject[rCount-1])->devpinitMessage(0))==0)

			outPutFile << ((nodeObject[rCount-1])->devpinitMessage(0))*100 << " MW" << endl;
		else
			outPutFile << -((nodeObject[rCount-1])->devpinitMessage(0))*100 << " MW" << endl;
		++rCount; // Increment the row count to point to the next node object
	}


	// Coefficients corresponding to lower generation limits
	outPutFile << "\nCoefficients corresponding to lower generation limits\n";
	for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
		lhs[rCount] = 0;
		lhs[rCount] += decvar[rCount - nodeNumber];
		modelSubnetQP->addConstr(lhs[rCount] >= ((*genIterator)->getPMin()));
		outPutFile << rCount << "\t" << (rCount - nodeNumber) << "\t" << 1.0 << "\t" << (*genIterator)->getPMin() << endl;
		outPutFile << rCount << "\t";
		outPutFile << ((*genIterator)->getPMin())*100 << " MW" << endl;
		++rCount; // Increment the row count to point to the next generator object
	}

	// Coefficients corresponding to upper generation limits
	outPutFile << "\nCoefficients corresponding to lower generation limits\n";
	for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
		lhs[rCount] = 0;
		lhs[rCount] += decvar[rCount - (genNumber + nodeNumber)];
		modelSubnetQP->addConstr(lhs[rCount] <= ((*genIterator)->getPMax()));
		outPutFile << rCount << "\t" << (rCount - (genNumber + nodeNumber)) << "\t" << 1.0 << "\t" << ((*genIterator)->getPMax()) << endl;
		outPutFile << rCount << "\t";
		outPutFile << ((*genIterator)->getPMax())*100 << " MW" << endl;
		++rCount; // Increment the row count to point to the next generator object
	}
	
	// Coefficients corresponding to intra-zone Line Forward Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to intra-zone Line Forward Flow Limit Constraints\n";
	for (tranIterator = translObject.begin(); tranIterator != translObject.end(); ++tranIterator){
		lhs[rCount] = 0;
		lhs[rCount] += (1/((*tranIterator)->getReactance()))*(decvar[genNumber + (*tranIterator)->getTranslNodeID1()]);
		outPutFile << "\n" << rCount << "\t" << genNumber + (*tranIterator)->getTranslNodeID1() << "\t" << 1/((*tranIterator)->getReactance()) << "\t" << 1/((*tranIterator)->getReactance()) << "\n";
		lhs[rCount] += (-1/((*tranIterator)->getReactance()))*(decvar[genNumber + (*tranIterator)->getTranslNodeID2()]);
		outPutFile << "\n" << rCount << "\t" << genNumber + (*tranIterator)->getTranslNodeID2() << "\t" << -1/((*tranIterator)->getReactance()) << "\t" << "-" << "\t" << -1/((*tranIterator)->getReactance()) << "\n";
		modelSubnetQP->addConstr(lhs[rCount] <= ((*tranIterator)->getFlowLimit()));
		outPutFile << rCount << "\t";
		outPutFile << ((*tranIterator)->getFlowLimit())*100 << " MW" << endl;
		++rCount; // Increment the row count to point to the next transmission line object
		
	}
	// Coefficients corresponding to intra-zone Line Reverse Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to intra-zone Line Reverse Flow Limit Constraints\n";
	for (tranIterator = translObject.begin(); tranIterator != translObject.end(); ++tranIterator){
		lhs[rCount] = 0;
		lhs[rCount] += (1/((*tranIterator)->getReactance()))*(decvar[genNumber + (*tranIterator)->getTranslNodeID1()]);
		outPutFile << "\n" << rCount << "\t" << genNumber + (*tranIterator)->getTranslNodeID1() << "\t" << 1/((*tranIterator)->getReactance()) << "\t" << 1/((*tranIterator)->getReactance()) << "\n";
		lhs[rCount] += (-1/((*tranIterator)->getReactance()))*(decvar[genNumber + (*tranIterator)->getTranslNodeID2()]);
		outPutFile << "\n" << rCount << "\t" << genNumber + (*tranIterator)->getTranslNodeID2() << "\t" << -1/((*tranIterator)->getReactance()) << "\t" << "-" << "\t" << -1/((*tranIterator)->getReactance()) << "\n";
		modelSubnetQP->addConstr(lhs[rCount] >= -((*tranIterator)->getFlowLimit()));
		outPutFile << rCount << "\t";
		outPutFile << -((*tranIterator)->getFlowLimit())*100 << " MW" << endl;
		++rCount; // Increment the row count to point to the next transmission line object
	}

	
	// Coefficients corresponding to shared existing Line Forward Flow Limit Constraints

	outPutFile << "\nCoefficients corresponding to shared existing Line Forward Flow Limit Constraints\n";
	for (exsharedIterator = SELineObject.begin(); exsharedIterator != SELineObject.end(); ++exsharedIterator){
		lhs[rCount] = 0;
		if ((*exsharedIterator)->getFlowDir() == 1) {
			lhs[rCount] += (1/((*exsharedIterator)->getReactance()))*(decvar[(genNumber + (*exsharedIterator)->getIntlNodeID())]);
			outPutFile << "\n" << rCount << "\t" << genNumber + (*exsharedIterator)->getIntlNodeID() << "\t" << 1/((*exsharedIterator)->getReactance()) << "\n";
			lhs[rCount] += (-1/((*exsharedIterator)->getReactance()))*(decvar[(genNumber + nodeNumber + (*exsharedIterator)->getExtNodeRank())]);
			outPutFile << "\n" << rCount << "\t" << genNumber + nodeNumber + (*exsharedIterator)->getExtNodeRank() << "\t" << -1/((*exsharedIterator)->getReactance()) << "\n";
			modelSubnetQP->addConstr(lhs[rCount] <= (*exsharedIterator)->getFlowLimit());
			outPutFile << rCount << "\t";
			outPutFile << ((*exsharedIterator)->getFlowLimit())*100 << " MW" << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		}
		else {
			lhs[rCount] += (1/((*exsharedIterator)->getReactance()))*(decvar[(genNumber + nodeNumber + (*exsharedIterator)->getExtNodeRank())]);
			outPutFile << "\n" << rCount << "\t" << genNumber + nodeNumber + (*exsharedIterator)->getExtNodeRank() << "\t" << 1/((*exsharedIterator)->getReactance()) << "\n";
			lhs[rCount] += (-1/((*exsharedIterator)->getReactance()))*(decvar[(genNumber + (*exsharedIterator)->getIntlNodeID())]);
			outPutFile << "\n" << rCount << "\t" << genNumber + (*exsharedIterator)->getIntlNodeID() << "\t" << -1/((*exsharedIterator)->getReactance()) << "\n";
			modelSubnetQP->addConstr(lhs[rCount] <= (*exsharedIterator)->getFlowLimit());
			outPutFile << rCount << "\t";
			outPutFile << ((*exsharedIterator)->getFlowLimit())*100 << " MW" << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		}
	}
	// Coefficients corresponding to shared existing Line Reverse Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to shared existing Line Reverse Flow Limit Constraints\n";
	for (exsharedIterator = SELineObject.begin(); exsharedIterator != SELineObject.end(); ++exsharedIterator){
		lhs[rCount] = 0;
		if ((*exsharedIterator)->getFlowDir() == 1) {
			lhs[rCount] += (1/((*exsharedIterator)->getReactance()))*(decvar[(genNumber + (*exsharedIterator)->getIntlNodeID())]);
			outPutFile << "\n" << rCount << "\t" << genNumber + (*exsharedIterator)->getIntlNodeID() << "\t" << 1/((*exsharedIterator)->getReactance()) << "\n";
			lhs[rCount] += (-1/((*exsharedIterator)->getReactance()))*(decvar[(genNumber + nodeNumber + (*exsharedIterator)->getExtNodeRank())]);
			outPutFile << "\n" << rCount << "\t" << genNumber + nodeNumber + (*exsharedIterator)->getExtNodeRank() << "\t" << -1/((*exsharedIterator)->getReactance()) << "\n";
			modelSubnetQP->addConstr(lhs[rCount] >= -((*exsharedIterator)->getFlowLimit()));
			outPutFile << rCount << "\t";

			outPutFile << -((*exsharedIterator)->getFlowLimit())*100 << " MW" << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		}
		else {
			lhs[rCount] += (1/((*exsharedIterator)->getReactance()))*(decvar[(genNumber + nodeNumber + (*exsharedIterator)->getExtNodeRank())]);
			outPutFile << "\n" << rCount << "\t" << genNumber + nodeNumber + (*exsharedIterator)->getExtNodeRank() << "\t" << 1/((*exsharedIterator)->getReactance()) << "\n";
			lhs[rCount] += (-1/((*exsharedIterator)->getReactance()))*(decvar[(genNumber + (*exsharedIterator)->getIntlNodeID())]);
			outPutFile << "\n" << rCount << "\t" << genNumber + (*exsharedIterator)->getIntlNodeID() << "\t" << -1/((*exsharedIterator)->getReactance()) << "\n";
			modelSubnetQP->addConstr(lhs[rCount] >= -((*exsharedIterator)->getFlowLimit()));
			outPutFile << rCount << "\t";
			outPutFile << -((*exsharedIterator)->getFlowLimit())*100 << " MW" << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		}
	}
	// Coefficients corresponding to realized shared candidate Line Forward Flow Limit Constraints

	outPutFile << "\nCoefficients corresponding to realized shared candidate Line Forward Flow Limit Constraints\n";
	for (candIterator = realCandLine.begin(); candIterator != realCandLine.end(); ++candIterator){
		lhs[rCount] = 0;
		if ((*candIterator)->getFlowDir() == 1) {
			lhs[rCount] += (1/((*candIterator)->getReactance()))*(decvar[(genNumber + (*candIterator)->getIntlNodeID())]);
			outPutFile << "\n" << rCount << "\t" << genNumber + (*exsharedIterator)->getIntlNodeID() << "\t" << 1/((*candIterator)->getReactance()) << "\n";
			lhs[rCount] += (-1/((*candIterator)->getReactance()))*(decvar[(genNumber + nodeNumber + (*candIterator)->getExtNodeRank())]);
			outPutFile << "\n" << rCount << "\t" << genNumber + nodeNumber + (*candIterator)->getExtNodeRank() << "\t" << -1/((*candIterator)->getReactance()) << "\n";
			modelSubnetQP->addConstr(lhs[rCount] <= (*candIterator)->getFlowLimit());
			outPutFile << rCount << "\t";
			outPutFile << ((*candIterator)->getFlowLimit())*100 << " MW" << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		}
		else {
			lhs[rCount] += (1/((*candIterator)->getReactance()))*(decvar[(genNumber + nodeNumber + (*candIterator)->getExtNodeRank())]);
			outPutFile << "\n" << rCount << "\t" << genNumber + nodeNumber + (*candIterator)->getExtNodeRank() << "\t" << 1/((*candIterator)->getReactance()) << "\n";
			lhs[rCount] += (-1/((*candIterator)->getReactance()))*(decvar[(genNumber + (*candIterator)->getIntlNodeID())]);
			outPutFile << "\n" << rCount << "\t" << genNumber + (*candIterator)->getIntlNodeID() << "\t" << -1/((*candIterator)->getReactance()) << "\n";
			modelSubnetQP->addConstr(lhs[rCount] <= (*candIterator)->getFlowLimit());
			outPutFile << rCount << "\t";
			outPutFile << ((*candIterator)->getFlowLimit())*100 << " MW" << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		}
	}
	// Coefficients corresponding to realized shared candidate Line Reverse Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to realized shared candidate Line Reverse Flow Limit Constraints\n";
	for (candIterator = realCandLine.begin(); candIterator != realCandLine.end(); ++candIterator){
		lhs[rCount] = 0;
		if ((*candIterator)->getFlowDir() == 1) {
			lhs[rCount] += (1/((*candIterator)->getReactance()))*(decvar[(genNumber + (*candIterator)->getIntlNodeID())]);
			outPutFile << "\n" << rCount << "\t" << genNumber + (*candIterator)->getIntlNodeID() << "\t" << 1/((*candIterator)->getReactance()) << "\n";
			lhs[rCount] += (-1/((*candIterator)->getReactance()))*(decvar[(genNumber + nodeNumber + (*candIterator)->getExtNodeRank())]);
			outPutFile << "\n" << rCount << "\t" << genNumber + nodeNumber + (*candIterator)->getExtNodeRank() << "\t" << -1/((*candIterator)->getReactance()) << "\n";
			modelSubnetQP->addConstr(lhs[rCount] >= -((*candIterator)->getFlowLimit()));
			outPutFile << rCount << "\t";

			outPutFile << -((*candIterator)->getFlowLimit())*100 << " MW" << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		}
		else {
			lhs[rCount] += (1/((*candIterator)->getReactance()))*(decvar[(genNumber + nodeNumber + (*candIterator)->getExtNodeRank())]);
			outPutFile << "\n" << rCount << "\t" << genNumber + nodeNumber + (*candIterator)->getExtNodeRank() << "\t" << 1/((*candIterator)->getReactance()) << "\n";
			lhs[rCount] += (-1/((*candIterator)->getReactance()))*(decvar[(genNumber + (*candIterator)->getIntlNodeID())]);
			outPutFile << "\n" << rCount << "\t" << genNumber + (*candIterator)->getIntlNodeID() << "\t" << -1/((*candIterator)->getReactance()) << "\n";
			modelSubnetQP->addConstr(lhs[rCount] >= -((*candIterator)->getFlowLimit()));
			outPutFile << rCount << "\t";
			outPutFile << -((*candIterator)->getFlowLimit())*100 << " MW" << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		}
	}
	// Coefficients corresponding to realized intra-zone candidate Line Forward Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to realized intra-zone candidate Line Forward Flow Limit Constraints\n";
	for (intCandIterator = realIntCandLineObject.begin(); intCandIterator != realIntCandLineObject.end(); ++intCandIterator){
		lhs[rCount] = 0;
		lhs[rCount] += (1/((*intCandIterator)->getReactance()))*(decvar[genNumber + (*intCandIterator)->getTranslNodeID1()]);
		outPutFile << "\n" << rCount << "\t" << genNumber + (*intCandIterator)->getTranslNodeID1() << "\t" << 1/((*intCandIterator)->getReactance()) << "\t" << 1/((*intCandIterator)->getReactance()) << "\n";
		lhs[rCount] += (-1/((*intCandIterator)->getReactance()))*(decvar[genNumber + (*intCandIterator)->getTranslNodeID2()]);
		outPutFile << "\n" << rCount << "\t" << genNumber + (*intCandIterator)->getTranslNodeID2() << "\t" << -1/((*intCandIterator)->getReactance()) << "\t" << "-" << "\t" << -1/((*intCandIterator)->getReactance()) << "\n";
		modelSubnetQP->addConstr(lhs[rCount] <= ((*intCandIterator)->getFlowLimit()));
		outPutFile << rCount << "\t";
		outPutFile << ((*intCandIterator)->getFlowLimit())*100 << " MW" << endl;
		++rCount; // Increment the row count to point to the next transmission line object
	}
	// Coefficients corresponding to realized intra-zone candidate Line Reverse Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to realized intra-zone candidate Line Reverse Flow Limit Constraints\n";
	for (intCandIterator = realIntCandLineObject.begin(); intCandIterator != realIntCandLineObject.end(); ++intCandIterator){
		lhs[rCount] = 0;
		lhs[rCount] += (1/((*intCandIterator)->getReactance()))*(decvar[genNumber + (*intCandIterator)->getTranslNodeID1()]);
		outPutFile << "\n" << rCount << "\t" << genNumber + (*intCandIterator)->getTranslNodeID1() << "\t" << 1/((*intCandIterator)->getReactance()) << "\t" << 1/((*intCandIterator)->getReactance()) << "\n";
		lhs[rCount] += (-1/((*intCandIterator)->getReactance()))*(decvar[genNumber + (*intCandIterator)->getTranslNodeID2()]);
		outPutFile << "\n" << rCount << "\t" << genNumber + (*intCandIterator)->getTranslNodeID2() << "\t" << -1/((*intCandIterator)->getReactance()) << "\t" << "-" << "\t" << -1/((*intCandIterator)->getReactance()) << "\n";
		modelSubnetQP->addConstr(lhs[rCount] >= -((*intCandIterator)->getFlowLimit()));
		outPutFile << rCount << "\t";
		outPutFile << -((*intCandIterator)->getFlowLimit())*100 << " MW" << endl;
		++rCount; // Increment the row count to point to the next transmission line object
	}

	outPutFile << "\nConstraint bounds (rows) Specified" << endl;
	outPutFile << "\nTotal number of rows: " << rCount - 1 << endl;

	outPutFile << "\nCoefficient Matrix specified" << endl;
	clock_t end1 = clock(); // stop the timer
	double elapsed_secs1 = double(end1 - begin) / CLOCKS_PER_SEC; // Calculate the time required to populate the constraint matrix and objective coefficients
	outPutFile << "\nTotal time taken to define the rows, columns, objective and populate the coefficient matrix = " << elapsed_secs1 << " s " << endl;
	// RUN THE OPTIMIZATION SIMULATION ALGORITHM //
	cout << "\nSimulation in Progress. Wait !!! ....." << endl;
	modelSubnetQP->optimize(); // Solves the optimization problem
	int stat = modelSubnetQP->get(GRB_IntAttr_Status); // Outputs the solution status of the problem 

	// DISPLAY THE SOLUTION DETAILS //
	if (stat == GRB_INFEASIBLE){
		outPutFile << "\nThe solution to the problem is INFEASIBLE." << endl;
		cout << "\nThe solution to the problem is INFEASIBLE." << endl;
		delete modelSubnetQP; // Free the memory of the GUROBI Problem Model
		return 0;
	} else if (stat == GRB_INF_OR_UNBD) {

		outPutFile << "\nNO FEASIBLE or BOUNDED solution to the problem exists." << endl;
		cout << "\nNO FEASIBLE or BOUNDED solution to the problem exists." << endl;
		delete modelSubnetQP; // Free the memory of the GUROBI Problem Model
		return 0;
	} else if (stat == GRB_UNBOUNDED) {
		outPutFile << "\nThe solution to the problem is UNBOUNDED." << endl;
		cout << "\nThe solution to the problem is UNBOUNDED." << endl;
		delete modelSubnetQP; // Free the memory of the GUROBI Problem Model
		return 0;
	} else if (stat == GRB_OPTIMAL) {
		outPutFile << "\nThe solution to the problem is OPTIMAL." << endl;

		cout << "\nThe solution to the problem is OPTIMAL." << endl;
	}

	//Get the Optimal Objective Value results//
	z = modelSubnetQP->get(GRB_DoubleAttr_ObjVal);
	for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
		z += (*genIterator)->getNLCost();
	}

	// Open separate output files for writing results of different variables
	string outIntAngFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outIntAnglesAPP/intAngleAPPGUROBI" + to_string(zonalIndex) + ".txt";
	string outCandFlowFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/sharedLineFlowAPP/candFlowMWAPPGUROBI" + to_string(zonalIndex) + ".txt";
	string outExtAngFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outExtAnglesAPP/extAngleAPPGUROBI" + to_string(zonalIndex) + ".txt";
	ofstream internalAngleOut(outIntAngFileName, ios::out);
	ofstream candFlowMWOut(outCandFlowFileName, ios::out);
	ofstream externalAngleOut(outExtAngFileName, ios::out);
	outPutFile << "\nThe Optimal Objective value (Generation Dispatch  and Line Building Decision cost) is: " << z << endl;
	powerGenOut << "\nThe Optimal Objective value (Generation Dispatch  and Line Building Decision cost) is: " << z << endl;
	cout << "\nThe Optimal Objective value (Generation Dispatch  and Line Building Decision cost) is: " << z << endl;
	vector<double> x; // Vector for storing decision variable output 
	x.push_back(0); // Initialize the decision Variable vector
	//thetaBuffer.push_back(0);
	//diffBuffer.push_back(0);

	//Display Power Generation
	powerGenOut << "\n****************** GENERATORS' POWER GENERATION LEVELS (MW) *********************" << endl;

	powerGenOut << "GENERATOR ID" << "\t" << "GENERATOR MW" << "\n";
	int arrayInd = 1;
	for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
		x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
		powerGenOut << (*genIterator)->getGenID() << "\t" << ((decvar[arrayInd]).get(GRB_DoubleAttr_X))*100 << " MW" << endl;
		++arrayInd;
	}
	powerGenOut << "Finished writing Power Generation" << endl;

	// Display Internal node voltage phase angle variables
	internalAngleOut << "\n****************** INTERNAL NODE VOLTAGE PHASE ANGLE VALUES *********************" << endl;
	internalAngleOut << "NODE ID" << "\t" << "VOLTAGE PHASE ANGLE" << "\n";
	for (nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); ++nodeIterator){
		if ( (((*nodeIterator)->getSharedFlag()) == 1) || (((*nodeIterator)->getCandFlag()) == 1) ) {
			x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			thetaBuffer.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			//diffBuffer.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			globRankBuffer.push_back((*nodeIterator)->getGlobalRank());
			internalAngleOut << (*nodeIterator)->getNodeID() << "\t" << ((decvar[arrayInd]).get(GRB_DoubleAttr_X)) << endl;
			//coordInstanceRef.populateAngleDec(((decvar[arrayInd]).get(GRB_DoubleAttr_X)), (zonalIndex-1), ((*nodeIterator)->getGlobalRank())); // Passing on the shared node angle decision message to the MO		
			++arrayInd;
		}
		else {
			x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			internalAngleOut << (*nodeIterator)->getNodeID() << "\t" << ((decvar[arrayInd]).get(GRB_DoubleAttr_X)) << endl;		
			++arrayInd;	
		}		
	}

	internalAngleOut << "Finished writing Internal Node Voltage Phase Angles" << endl;

	// Display Outer Zonal node angles
	externalAngleOut << "\n****************** OUTER ZONAL NODE VOLTAGE PHASE ANGLE VALUES *********************" << endl;
	externalAngleOut << "EXTERNAL NODE ID" << "\t" << "VOLTAGE PHASE ANGLE" << "\n";
	int diffNodeC = 0;
	for (globalIterator = globalRankDiffNode.begin(); globalIterator != globalRankDiffNode.end(); ++globalIterator){
		if (diffNodeC > 0) { // Skip the dummy element "0" at the beginning
			x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			thetaBuffer.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			//diffBuffer.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			globRankBuffer.push_back((*globalIterator));
			externalAngleOut << (*globalIterator) << "\t" << ((decvar[arrayInd]).get(GRB_DoubleAttr_X)) << endl;
			//coordInstanceRef.populateAngleDec(((decvar[arrayInd]).get(GRB_DoubleAttr_X)), (zonalIndex-1), (*globalIterator)); // Passing on the shared node angle decision message to the MO
			++arrayInd;
		}
		++diffNodeC;
	}
	externalAngleOut << "Finished writing outer zonal node voltage phase angle values" << endl;
	// Display Shared Candidate lines' Power Flow values
	candFlowMWOut << "\n****************** INTERNAL TRANSMISSION LINES POWER FLOW VALUES *********************" << endl;
	candFlowMWOut << "INTERNAL TRANSMISSION LINE ID" << "\t" << "FROM" << "\t" << "TO" << "\t" << "MW FLOW VALUE" << "\n";
	for (tranIterator = translObject.begin(); tranIterator != translObject.end(); ++tranIterator){
		candFlowMWOut << (*tranIterator)->getTranslID() << "\t" << (*tranIterator)->getTranslNodeID1() << "\t" << (*tranIterator)->getTranslNodeID2() << "\t" << (1/((*tranIterator)->getReactance()))*(((decvar[genNumber + (*tranIterator)->getTranslNodeID1()]).get(GRB_DoubleAttr_X))-((decvar[genNumber + (*tranIterator)->getTranslNodeID2()]).get(GRB_DoubleAttr_X)))*100 << " MW" << endl;;
		++arrayInd; // Increment the row count to point to the next transmission line object
	}
	candFlowMWOut << "\n****************** SHARED EXISTING TRANSMISSION LINES POWER FLOW VALUES *********************" << endl;
	candFlowMWOut << "SHARED EXISTING TRANSMISSION LINE ID" << "FROM" << "\t" << "TO" << "\t" << "MW FLOW VALUE" << "\n";
	for (exsharedIterator = SELineObject.begin(); exsharedIterator != SELineObject.end(); ++exsharedIterator){
		if ((*exsharedIterator)->getFlowDir() == 1) {
			candFlowMWOut << (*exsharedIterator)->getTranslID() << "\t" << (*exsharedIterator)->getIntlNodeID() << "\t" << (*exsharedIterator)->getExtNodeRank() << " (ext. node)" << "\t" << (1/((*exsharedIterator)->getReactance()))*(((decvar[(genNumber + (*exsharedIterator)->getIntlNodeID())]).get(GRB_DoubleAttr_X))-((decvar[(genNumber + nodeNumber + (*exsharedIterator)->getExtNodeRank())]).get(GRB_DoubleAttr_X)))*100 << " MW" << endl;;
		}
		else {
			candFlowMWOut << (*exsharedIterator)->getTranslID() << "\t" << (*exsharedIterator)->getExtNodeRank() << " (ext. node)" << "\t" << (*exsharedIterator)->getIntlNodeID() << "\t" << (1/((*exsharedIterator)->getReactance()))*(((decvar[(genNumber + nodeNumber + (*exsharedIterator)->getExtNodeRank())]).get(GRB_DoubleAttr_X))-((decvar[(genNumber + (*exsharedIterator)->getIntlNodeID())]).get(GRB_DoubleAttr_X)))*100 << " MW" << endl;;
		}
	}
	candFlowMWOut << "\n****************** REALIZED SHARED CANDIDATE TRANSMISSION LINES POWER FLOW VALUES *********************" << endl;
	candFlowMWOut << "REALIZED SHARED CANDIDATE LINE ID" << "FROM" << "\t" << "TO" << "\t" << "MW FLOW VALUE" << "\n";
	for (candIterator = realCandLine.begin(); candIterator != realCandLine.end(); ++candIterator){
		if ((*candIterator)->getFlowDir() == 1) {
			candFlowMWOut << (*candIterator)->getTranslID() << "\t" << (*candIterator)->getIntlNodeID() << "\t" << (*candIterator)->getExtNodeRank() << " (ext. node)" << "\t" << (1/((*candIterator)->getReactance()))*(((decvar[(genNumber + (*candIterator)->getIntlNodeID())]).get(GRB_DoubleAttr_X))-((decvar[(genNumber + nodeNumber + (*candIterator)->getExtNodeRank())]).get(GRB_DoubleAttr_X)))*100 << " MW" << endl;;
		}
		else {
			candFlowMWOut << (*candIterator)->getTranslID() << "\t" << (*candIterator)->getExtNodeRank() << " (ext. node)" << "\t" << (*candIterator)->getIntlNodeID() << "\t" << (1/((*candIterator)->getReactance()))*(((decvar[(genNumber + nodeNumber + (*candIterator)->getExtNodeRank())]).get(GRB_DoubleAttr_X))-((decvar[(genNumber + (*candIterator)->getIntlNodeID())]).get(GRB_DoubleAttr_X)))*100 << " MW" << endl;;
		}
	}
	candFlowMWOut << "\n****************** REALIZED INTERNAL CANDIDATE TRANSMISSION LINES POWER FLOW VALUES *********************" << endl;
	candFlowMWOut << "REALIZED INTERNAL CANDIDATE LINE ID" << "FROM" << "\t" << "TO" << "\t" << "MW FLOW VALUE" << "\n";
	for (intCandIterator = realIntCandLineObject.begin(); intCandIterator != realIntCandLineObject.end(); ++intCandIterator){
		candFlowMWOut << (*intCandIterator)->getTranslID() << "\t" << (*intCandIterator)->getTranslNodeID1() << "\t" << (*intCandIterator)->getTranslNodeID2() << "\t" << (1/((*intCandIterator)->getReactance()))*(((decvar[genNumber + (*intCandIterator)->getTranslNodeID1()]).get(GRB_DoubleAttr_X))-((decvar[genNumber + (*intCandIterator)->getTranslNodeID2()]).get(GRB_DoubleAttr_X)))*100 << " MW" << endl;;
	}
	candFlowMWOut << "Finished writing lines' Power Flows" << endl;

	delete modelSubnetQP; // Free the memory of the GUROBI Problem Model
	clock_t end2 = clock(); // stop the timer
	double elapsed_secs2 = double(end2 - begin) / CLOCKS_PER_SEC; // Calculate the Total Time
	outPutFile << "\nTotal time taken to solve the QP Line Construction Decision Making Problem instance and retrieve the results = " << elapsed_secs2 << " s " << endl;
	cout << "\nTotal time taken to solve the QP Line Construction Decision Making Problem instance and retrieve the results = " << elapsed_secs2 << " s " << endl;

	// Close the different output files
	outPutFile.close();
	powerGenOut.close();
	internalAngleOut.close();
	candFlowMWOut.close();
	externalAngleOut.close();
	cout << "\nSimulation Completed.\nResults written on the different output files" << endl;
	return z;
}

vector<double> Nettran::getZonalDecision() // Returns the intermediate decision variable values from APP
{
	return thetaBuffer;
}

vector<int> Nettran::getZonalRanks() // Returns the global ranks of the shared decision variables from APP
{
	return globRankBuffer;
}
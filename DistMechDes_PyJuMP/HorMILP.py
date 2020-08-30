# HorMILP.cpp : Defines the entry point for Stage-I and Stage-II of the Horizontal Investment Coordination MILP Market Mechanism Design Simulation application.
# Main Method for running the Horizontal Investment Coordination Stage-I and Stage-II MILP Market Mechanism Design Simulation based on Distributed Stochastic UC and APP; Parts of the code follows the code design philosophy of Nick Laws of NREL
# Added Hosna Khajeh and Dr. Mohammad Reza Hesamzadeh as collaborators on 30th August, 2020
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
from Python_src.nettran import Nettran
from Python_src.marketoverseer import marketOverseer

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
def HorMILPCentral(): # Main method begins program execution
    	'''Future Work
	#Choose the type of objective function
    	'''
	systemChoice = int(input("Choose the type of System to be simulated: 1 for Simple two bus/two region, 2 for system combined of IEEE 14, 30, and 5 node systems"))
	basisForComparison = int(input("Choose, for updating the Lagrange multipliers, 0 when the current zonal updates are compared to previous iteration MO update;else, 1"))
	curveChoice = 1 # Number to indicate the type of Objective function among average heat rate, piecewise linear, or polynomial; Assume Average Heat Rate for now
	# Read the master zones file, for deciding upon which other files to read for building the model
	#Read the master zones file, for deciding upon which other files to read for building the model
	int numberOfFields; // Number of rows or individual file types for each of the zones
	string inputMasterFile;
	if (systemChoice==1)
		inputMasterFile = "masterZonesSummary.json";
	else
		inputMasterFile = "masterZonesSummaryRevised.json";
	UBIterate = 0.0 #Initial value for the global upper bound iterates at the end of every iteration
	LBIterate = 0.0 #Initial value for the global lower bound iterates at the end of every iteration
	UBItVec = [] #Vector for storing the values of the global upper bound iterates for every iteration
	LBItVec = [] #Vector for storing the values of the global lower bound iterates for every iteration
	ratioIterate = [] #Vector for storing the values of UB/LB iterates for every iteration
	globalCons = [] #Vector for storing the values of global consensus for every iteration
	ifstream zoneSummaryFile( inputMasterFile, ios::in ); // ifstream constructor opens the master zones summary file
	stringstream buffer; // stringstream object to store the read information from the summary file
	// exit program if ifstream could not open file
	if ( !zoneSummaryFile ) {
		cerr << "\nMaster file for Zones Summary could not be opened\n" << endl;
		exit( 1 );
	} // end if

	//zoneSummaryFile >> numberOfZones >> numberOfFields; // get the number of zones and the number of fields: Future expansion
	#Number of zones between which horizontal investment coordination for transmission lines to be built is considered
	numberOfZones = int(input("\nEnter the number of zones")) #User input the number of zones/regions
	solverChoice = int(input("\nChoose either the GLPK (1) or GUROBI (2) as the Solver. ")) #Choice of the solver
	if solverChoice==1:
		lpMethodChoice = int(input("\nChoose either the Simplex LP Rlaxation (1) or Interior Point Method LP Relaxation (2) as the method to provide the initial basis to the Mixed Integer Unit Commitment Problem. ")) #Simplex or interior method algorithm for MILP
	}
	#GRBEnv* environmentGUROBI = new GRBEnv("GUROBILogFile.log"); // GUROBI Environment object for storing the different optimization models
	numberOfFields = 7; // Number of fields
   	buffer << zoneSummaryFile.rdbuf(); // reads the data in the summary file 
   	string test = buffer.str(); // Extract the strings from the buffer to "test"

   	//create variables that will act as "cursors". we'll take everything between them.
   	size_t pos1 = 0;
   	size_t pos2;
   	//create the array to store the strings.
   	string str[numberOfFields*numberOfZones];
	#Read the summary input file
   	for i in range(numberOfFields):
		for j in range(numberOfZones):
			if j==numberOfZones-1:
				pos2 = test.find("\n", pos1); //search for the bar "\n". pos2 will be where the bar was found.
        			str[i*numberOfZones+j] = test.substr(pos1, (pos2-pos1)); //make a substring, wich is nothing more 
                                              //than a copy of a fragment of the big string.
        			pos1 = pos2+1; // sets pos1 to the next character after pos2.
			else:
        			pos2 = test.find(" ", pos1); //search for the bar " ". pos2 will be where the bar was found.
        			str[i*numberOfZones+j] = test.substr(pos1, (pos2-pos1)); //make a substring, wich is nothing more 
                                              //than a copy of a fragment of the big string.
        			pos1 = pos2+1; // sets pos1 to the next character after pos2. 
    			}
    		}
 	}	

	zonalNetVector = [] #Vector of zonal network objects
	log.info("\n*** NETWORK INITIALIZATION STAGE BEGINS ***\n")
	upperBoundVector = np.zeros(float, numberOfZones) #Vector of upper bounds by iteration
	lowerBoundVector = np.zeros(float, numberOfZones) #Vector of lower bounds by iteration
	for i in range(numberOfZones):
		zonalNetVector.append(Nettran( str, (i+1), numberOfZones, curveChoice, lpMethodChoice )) #push to the vector of zonal networks
	log.info("\n*** NETWORK INITIALIZATION STAGE ENDS: ZONAL SUB-NETWORKS CREATED ***\n")
	"""
	#Initialize the Lagrange Multipliers/Dual Variables/Rewards/Penalties to zero for the first iteration
	LagMultXi = np.zeros(float, 1000) #Lagrange Multipliers/Dual Variables/Rewards/Penalties corresponding to bus voltage phase angle consensus
	LagMultPi = np.zeros(float, 1000) #Lagrange Multipliers/Dual Variables/Rewards/Penalties corresponding to line building decision integer variable consensus
	LagCombXi = np.zeros(float, 1000) #Lagrange Multipliers/Dual Variables/Rewards/Penalties corresponding to bus voltage phase angle consensus for MO
	LagCombPi = np.zeros(float, 1000) #Lagrange Multipliers/Dual Variables/Rewards/Penalties corresponding to line building decision integer variable consensus MO
	"""
	#Stitching the different nodes of shared existing & shared candidate lines and zones in one list
	sharedNodeList = [] #List of all the node IDs of the shared existing and/or shared candidate transmission lines
	sharedZoneList = [] #List of all the corresponding zone IDs of the respective shared line end nodes
	sharedGlobalNodeList = [] #Global ranking of the shared nodes		
	sharedNodeList.append(0) #Initialize the list of shared line end node ID's so that it's not empty
	sharedZoneList.append(0) #Initialize the list of shared line end zone ID's so that it's not empty
	sharedGlobalNodeList.append(0) #Initialize the list of shared global line end zone ID's so that it's not empty
	containsFlag = 0 #Initializes the flag to test if the list so far contains a particular node, to avoid duplication
	otherNodeCount = 0 #Initializes the count of total number of shared nodes
	sharedZone = 0 #Temporary integer for storing the shared zone index
	sharedNode = 0 #Temporary integer for storing the shared node index
	sharedSESerList = [] #List of all the unique global serial numbers of the shared existing lines
	sharedSESerList.append(0) #Initialize the list of shared existing line global serial numbers so that it's not empty
	sharedSEFromList = [] #List of the From node ranks of the shared existing lines
	sharedSEToList = [] #List of the To node ranks of the shared existing lines
	sharedSEReactList = [] #List of all the reactances of the shared existing lines
	sharedSECapList = [] #List of all the flow capacities of the shared existing lines
	containsSEFlag = 0 #Initializes the flag to test if the list so far contains the shared existing line, to avoid duplication
	SELineCount = 0 #Initializes the count of total number of shared existing lines
	sharedELine = 0 #Temporary integer for storing the shared existing line global serial number
	sharedCandSerList  = [] #List of all the unique global serial numbers of the shared candidate lines
	sharedCandSerList.append(0) #Initialize the list of shared candidate line global serial numbers so that it's not empty
	sharedCandFromList  = [] #List of the From node ranks of the shared candidate lines
	sharedCandToList  = [] #List of the To node ranks of the shared candidate lines
	sharedCandReactList  = [] #List of all the reactances of the shared candidate lines
	sharedCandCapList  = [] #List of all the flow capacities of the shared candidate lines
	containsCandFlag = 0 #Initializes the flag to test if the list so far contains the shared candidate line, to avoid duplication
	candLineCount = 0 #Initializes the count of total number of shared candidate lines
	sharedCandLine = 0 #Temporary integer for storing the shared candidate line global serial number
	#Stitching together the shared line end nodes from different zones into one master list
	#log.info("Stitched set of nodes at the ends of shared existing and shared candidate lines with global rankings")
	for i in range(numberOfZones): #Run the loop on all the zones
		#log.info("Zone Considered is {}".format(i))
		sharedZone = 0 #Temporary integer for storing the shared zone index
		sharedNode = 0 #Temporary integer for storing the shared node index
		sharedNodeCounter = 1
		while (sharedZone != -1) and (sharedNode != -1): #Run the loop on shared lists
			sharedZone = zonalNetVector[i].getConnZone(sharedNodeCounter) #Gets the zone ID of the outer zonal node
			sharedNode = zonalNetVector[i].getConnNode(sharedNodeCounter) #Gets the node ID of the outer zonal node
			indCount = 0 #Initialize a counter for tracking the position in the vector of the iterator
			if (sharedZone != -1) and (sharedNode != -1): #as long as end of the zonal vectors are not reached
				for diffZNIt in sharedNodeList: #Iterate over the vector of connected external zone nodes
					if diffZNIt == sharedNode and sharedZoneList[indCount] == sharedZone: #Check whether the other-zone node is already present in the list
						containsFlag = 1 #Loop through the diffZoneNodeID and diffZoneID vectors, to see if the other zone-node is already present or not, if yes, set the containsFlag to 1
					indCount += 1 #Increment the counter
				if containsFlag == 0: #If the diffZoneNodeID and diffZoneID vectors do not contain the other zone-node, then push the node and the other zone index in the respective vectors
					sharedNodeList.append(sharedNode)
					sharedZoneList.append(sharedZone)
					#log.info("Local Node Index: {}".format(sharedNode))
					#log.info("Zone Index: {}".format(sharedZone))
					otherNodeCount += 1 #Increment the counter to account for the total number of other-zone nodes
					#log.info("Rank: {}".format(otherNodeCount))
					sharedGlobalNodeList.append(otherNodeCount)		
				containsFlag = 0 #Reset the containsFlag for matching the next item
			sharedNodeCounter += 1
	#Stitching together the shared existing line list from different zones into one master list
	for i in range(numberOfZones): #Run the loop on all the zones
		sharedELine = 0 #Temporary integer for storing the shared existing line global serial number
		sharedELineCounter = 0
		while sharedELine != -1: #Run the loop on shared existing lines
			sharedELine = zonalNetVector[i].getSESerial(sharedELineCounter) #Gets the global serial number of the shared existing line
			SEFromNode = zonalNetVector[i].getSEFromNode(sharedELineCounter) #Gets the from node of the shared existing line
			SEFromZone = zonalNetVector[i].getSEFromZone(sharedELineCounter) #Gets the from zone of the shared existing line
			SEToNode = zonalNetVector[i].getSEToNode(sharedELineCounter) #Gets the to node of the shared existing line
			SEToZone = zonalNetVector[i].getSEToZone(sharedELineCounter) #Gets the to zone of the shared existing line
			SEImp = zonalNetVector[i].getSEReactance(sharedELineCounter) #Gets the reactance of the shared existing line
			SECap = zonalNetVector[i].getSECapacity(sharedELineCounter) #Gets the MW flow limit of the shared existing line
			if sharedELine != -1: #as long as end of the zonal shared existing line vector is not reached
				for SESerList in sharedSESerList: #Iterate over the list of SELine
					if SESerList == sharedELine: #Check whether the SE line is already present in the list
						containsSEFlag = 1 #Loop through the diffZoneNodeID and diffZoneID vectors, to see if the other zone-node is already present or not, if yes, set the containsFlag to 1
						indCount = 0 #Initialize a counter for tracking the position in the vector of the iterator
						for diffZNIt in sharedNodeList:
							if diffZNIt == SEFromNode and sharedZoneList[indCount] == SEFromZone: #Check whether the from node is already present in the list
								SEFromRank = indCount #if yes, get the rank of the from node from the list (incrementing indCount will point to the rank, since sharedNodeList and sharedZoneList both have a dummy 0 as their fist element
								if i+1 == SEFromZone: #if the present zone has the from node
									zonalNetVector[i].setSEFromRank(sharedELineCounter, SEFromRank) #assign the rank to the internal from node
								else: #else
									zonalNetVector[i].setSEFromRankConn(sharedELineCounter, SEFromRank) #assign the rank of the external from node and also store the rank to the list of rank of connected nodes to the internal to node
							elif diffZNIt == SEToNode and sharedZoneList[indCount] == SEToZone: #Check whether the to node is already present in the list
								SEToRank = indCount #if yes, get the rank of the to node from the list (incrementing indCount will point to the rank, since sharedNodeList and sharedZoneList both have a dummy 0 as their fist element
								if i+1 == SEToZone: #if the present zone has the to node
									zonalNetVector[i].setSEToRank(sharedELineCounter, SEToRank) #assign the rank to the internal to node
								else: #else
									zonalNetVector[i].setSEToRankConn(sharedELineCounter, SEToRank) #assign the rank of the external to node and also store the rank to the list of rank of connected nodes to the internal to node
							indCount += 1 #Increment the counter
				if containsSEFlag == 0: #If the sharedSESerList vector does not contain the SE line, then push the SE line in the vector
					#log.info("From node of SE line {} is {} From zone is {} To node is {} To zone is {}".format(sharedELine, SEFromNode, SEFromZone, SEToNode, SEToZone))
					sharedSESerList.append(sharedELine)
					indCount = 0 #Initialize a counter for tracking the position in the vector of the iterator
					for diffZNIt in sharedNodeList:
						if diffZNIt == SEFromNode and sharedZoneList[indCount] == SEFromZone: #Check whether the from node is already present in the list
							SEFromRank = indCount #if yes, get the rank of the from node from the list (incrementing indCount will point to the rank, since sharedNodeList and sharedZoneList both have a dummy 0 as their fist element
							if i+1 == SEFromZone: #if the present zone has the from node
								zonalNetVector[i].setSEFromRank(sharedELineCounter, SEFromRank) #assign the rank to the internal from node
							else: #else
								zonalNetVector[i].setSEFromRankConn(sharedELineCounter, SEFromRank) #assign the rank of the external from node and also store the rank to the list of rank of connected nodes to the internal to node
							sharedSEFromList.append(SEFromRank)
							#log.info("From node of SE line {} is {} From zone is {} and from rank is {}".format(sharedELine,SEFromNode,SEFromZone,SEFromRank))
						elif diffZNIt == SEToNode and sharedZoneList[indCount] == SEToZone: #Check whether the to node is already present in the list
							SEToRank = indCount #if yes, get the rank of the to node from the list (incrementing indCount will point to the rank, since sharedNodeList and sharedZoneList both have a dummy 0 as their fist element
							if i+1 == SEToZone : #if the present zone has the to node
								zonalNetVector[i].setSEToRank(sharedELineCounter, SEToRank) #assign the rank to the internal to node
							else: #else
								zonalNetVector[i].setSEToRankConn(sharedELineCounter, SEToRank) #assign the rank of the external to node and also store the rank to the list of rank of connected nodes to the internal to node
							sharedSEToList.append(SEToRank)
							#log.info("To node of SE line {} is {} To zone is {} and to rank is {}".format(sharedELine, SEToNode, SEToZone, SEToRank))
						indCount+=1 #Increment the counter
					sharedSEReactList.append(SEImp)
					sharedSECapList.append(SECap)
					SELineCount += 1 #Increment the counter to account for the total number of SE Lines
				containsSEFlag = 0 #Reset the containsFlag for matching the next item
		sharedELineCounter += 1
	#Stitching together the shared candidate line list from different zones into one master list
	for i in range(numberOfZones): #Run the loop on all the zones
		sharedCandLine = 0 #Temporary integer for storing the shared candidate line global serial number
		candLineCounter = 0
		while sharedCandLine != -1: #Run the loop on shared candidate lines
			sharedCandLine = zonalNetVector[i].getCandSerial(candLineCounter) #Gets the global serial number of the shared candidate line
			CandFromNode = zonalNetVector[i].getCandFromNode(candLineCounter) #Gets the from node of the shared candidate line
			CandFromZone = zonalNetVector[i].getCandFromZone(candLineCounter) #Gets the from zone of the shared candidate line
			CandToNode = zonalNetVector[i].getCandToNode(candLineCounter) #Gets the to node of the shared candidate line
			CandToZone = zonalNetVector[i].getCandToZone(candLineCounter) #Gets the to zone of the shared candidate line
			CandImp = zonalNetVector[i].getCandReactance(candLineCounter) #Gets the reactance of the shared candidate line
			CandCap = zonalNetVector[i].getCandCapacity(candLineCounter) #Gets the MW flow limit of the shared candidate line
			if sharedCandLine != -1: #as long as end of the zonal candidate line vector is not reached
				candGlobalRank = 0 #Global ranking of the candidate line in the list of shared candidate lines
				for CandSerList in sharedCandSerList: #Iterate over candline
					if CandSerList == sharedCandLine: #Check whether the Cand line is already present in the list
						containsCandFlag = 1 #Loop through the diffZoneNodeID and diffZoneID vectors, to see if the other zone-node is already present or not, if yes, set the containsFlag to 1 
						zonalNetVector[i].assignCandGlobalRank(candLineCounter, candGlobalRank) #assigns the global rank of the shared candidate line
						indCount = 0 #Initialize a counter for tracking the position in the vector of the iterator
						for diffZNIt in sharedNodeList:
							if diffZNIt == CandFromNode and sharedZoneList[indCount] == CandFromZone: #Check whether the from node is already present in the list
								CandFromRank = indCount #if yes, get the rank of the from node from the list (incrementing indCount will point to the rank, since sharedNodeList and sharedZoneList both have a dummy 0 as their fist element
								if i+1 == CandFromZone: #if the present zone has the from node
									zonalNetVector[i].setCandFromRank(candLineCounter, CandFromRank) #assign the rank to the internal from node
								else: #else
									zonalNetVector[i].setCandFromRankConn(candLineCounter, CandFromRank) #assign the rank of the external from node and also store the rank to the list of rank of connected nodes to the internal to node
							elif diffZNIt == CandToNode and sharedZoneList[indCount] == CandToZone: #Check whether the to node is already present in the list
								CandToRank = indCount #if yes, get the rank of the to node from the list (incrementing indCount will point to the rank, since sharedNodeList and sharedZoneList both have a dummy 0 as their fist element
								if i+1 == CandToZone: #if the present zone has the to node
									zonalNetVector[i].setCandToRank(candLineCounter, CandToRank) #assign the rank to the internal to node
								else: #else
									zonalNetVector[i].setCandToRankConn(candLineCounter, CandToRank) #assign the rank of the external to node and also store the rank to the list of rank of connected nodes to the internal to node
							indCount+=1 #Increment the counter
					candGlobalRank+=1 #Increment the shared candidate line rank
				if containsCandFlag == 0: #If the diffZoneNodeID and diffZoneID vectors do not contain the other zone-node, then push the node and the other zone index in the respective vectors
					#log.info("From node of candidate line {" << sharedCandLine << "} is {" << CandFromNode << "} From zone is {" << CandFromZone << "} To node is {" << CandToNode << "} To zone is {}".format())
					sharedCandSerList.append(sharedCandLine)
					zonalNetVector[i].assignCandGlobalRank(candLineCounter, candLineCount+1) #assigns the global rank of the shared candidate line
					indCount = 0 #Initialize a counter for tracking the position in the vector of the iterator
					for diffZNIt in sharedNodeList:
						if diffZNIt == CandFromNode and sharedZoneList[indCount] == CandFromZone: #Check whether the from node is already present in the list
							CandFromRank = indCount #if yes, get the rank of the from node from the list (incrementing indCount will point to the rank, since sharedNodeList and sharedZoneList both have a dummy 0 as their fist element
							if i+1 == CandFromZone: #if the present zone has the from node
								zonalNetVector[i].setCandFromRank(candLineCounter, CandFromRank) #assign the rank to the internal from node
							else: #else
								zonalNetVector[i].setCandFromRankConn(candLineCounter, CandFromRank) #assign the rank of the external from node and also store the rank to the list of rank of connected nodes to the internal to node
							sharedCandFromList.append(CandFromRank)
							#log.info("From node of candidate line {} is {} From zone is {} and from rank is {}".format(sharedCandLine,CandFromNode,CandFromZone,CandFromRank))
						elif diffZNIt == CandToNode and sharedZoneList[indCount] == CandToZone: #Check whether the to node is already present in the list
							CandToRank = indCount #if yes, get the rank of the to node from the list (incrementing indCount will point to the rank, since sharedNodeList and sharedZoneList both have a dummy 0 as their fist element
							if i+1 == CandToZone: #if the present zone has the to node
								zonalNetVector[i].setCandToRank(candLineCounter, CandToRank) #assign the rank to the internal to node
							else: #else
								zonalNetVector[i].setCandToRankConn(candLineCounter, CandToRank) #assign the rank of the external to node and also store the rank to the list of rank of connected nodes to the internal to node
							sharedCandToList.pappend(CandToRank)
							#log.info("To node of candidate line {} is {} To zone is {} and to rank is {}".format(sharedCandLine,CandToNode,CandToZone,CandToRank))
						indCount += 1 #Increment the counter
					sharedCandReactList.append(CandImp)
					sharedCandCapList.append(CandCap)
					candLineCount += 1 #Increment the counter to account for the total number of SE Lines	
				containsCandFlag = 0 #Reset the containsFlag for matching the next item
		candLineCounter += 1
	log.info("\n*** MARKET OVERSEER INITIALIZATION STAGE BEGINS ***\n")
	scenarios = zonalNetVector[0].getNumberOfScenarios() #Returns the number of scenarios
	marketoverInstance = Marketover( numberOfZones, otherNodeCount, SELineCount, candLineCount, zonalNetVector, sharedNodeList, sharedGlobalNodeList, sharedZoneList, sharedSESerList, sharedSEFromList, sharedSEToList, sharedSEReactList, sharedSECapList, sharedCandSerList, sharedCandFromList, sharedCandToList, sharedCandReactList, sharedCandCapList, lpMethodChoice, scenarios, basisForComparison ) #create the market overseer instance//%%
	log.info("\n*** MARKET OVERSEER INITIALIZATION STAGE ENDS ***\n")
	log.info("\nShared Node List")
	for sharedNodeIterator in sharedNodeList:
		log.info("\nRank of the nodes in serial order is: {}".format(sharedNodeIterator))
	for sharedCandLineIterator in sharedCandSerList:
		log.info("\nRank of the candidate lines in serial order is: {}".format(sharedCandLineIterator))
	#Initialize the Lagrange Multipliers/Dual Variables/Rewards/Penalties to zero for the first iteration
	LagMultXi = np.zeros(float, (scenarios+1)*numberOfZones*otherNodeCount) #Lagrange Multipliers/Dual Variables/Rewards/Penalties corresponding to bus voltage phase angle consensus
	#Initialize the Lagrange Multipliers/Dual Variables/Rewards/Penalties to zero for the first iteration
	LagMultPi = np.zeros(float, (numberOfZones*candLineCount)) #Lagrange Multipliers/Dual Variables/Rewards/Penalties corresponding to line building decision integer variable consensus
	LagCombXi = np.zeros(float, (scenarios+1)*otherNodeCount) #Lagrange Multipliers/Dual Variables/Rewards/Penalties corresponding to bus voltage phase angle consensus for MO
	LagCombPi= np.zeros(float, candLineCount) #Lagrange Multipliers/Dual Variables/Rewards/Penalties corresponding to line building decision integer variable consensus MO
	iterationCounter = 1 #Initialize the iteration counter

	log.info("\n*** DISTRIBUTED STOCHASTIC OPTIMIZATION ALGORITHMIC MARKET MECHANISM DESIGN FIRST STAGE BEGINS ***\n")
	#do {//%%
	for iterCountOut in range(100):
		if iterCountOut != 0:
			marketoverInstance.clearDelayedVectors() #clear the different interim delayed vectors for making them ready for next iteration 
			marketoverInstance.bufferintermediateDecision(iterCountOut) #Buffering the previous iterations' MO decision values, (for comparison basis equal to 0)			
			marketoverInstance.clearVectors() #clear the different interim vectors for making them ready for next iteration
		else:
			marketoverInstance.bufferintermediateDecision(iterCountOut) #Buffering the previous iterations' MO decision values, (for comparison basis equal to 0)			 
		UBIterate = 0.0
		LBIterate = 0.0 
		log.info("\n*** ITERATION {} BEGINS ***\n".format(iterCountOut+1))
		for i  in range(numberOfZones): #Each region solves its own MILP optimization problem 
			log.info("\n*** MIXED INTEGER LINEAR PROGRAMMING FOR ZONE " << i+1 << " BEGINS ***\n")
			log.info("\nZonal Calculations of Beliefs about the Investment decision MILP begins")
			log.info("\nSOLVING MILP")
			upperBoundVector[i]=zonalNetVector[i].MILPAvgHR(solverChoice, marketoverInstance, LagMultXi, LagMultPi, candLineCount, otherNodeCount) #Perform unit commitment for average heat rate objective
			UBIterate += upperBoundVector[i]
			log.info("\nMILP SOLVED")
			log.info("\nESTIMATING LOWER BOUND")
			lowerBoundVector[i]=zonalNetVector[i].calcMILPBounds(solverChoice, LagMultXi, LagMultPi, candLineCount, otherNodeCount) #Calculate the bounds//%%
			LBIterate += lowerBoundVector[i]
			log.info("\nLOWER BOUND ESTIMATED")
		for k in range(scenarios+1):
			for i in range(otherNodeCount): #Initialize the Lagrange Multipliers/Dual Variables/Rewards/Penalties to zero for the MO
				tempXi = 0
				for j in range(numberOfZones):
					tempXi += LagMultXi[k*numberOfZones*otherNodeCount+j*otherNodeCount+i]
				LagCombXi[k*otherNodeCount+i] = tempXi
		for i in range(candLineCount): #Initialize the Lagrange Multipliers/Dual Variables/Rewards/Penalties to zero for the MO
			tempPi = 0
			for j in range(numberOfZones):
				tempPi += LagMultPi[j*candLineCount+i]
			LagCombPi[i] = tempPi
		log.info("\nSOLVING MILP FOR MARKET OVERSEER/TRANSMISSION PLANNING COORDINATOR")
		marketoverInstance.MILPMarketover(solverChoice, LagCombXi, LagCombPi, candLineCount, otherNodeCount) #The MO solves its own optimization problem and also updates the rewards/penalties//%%
		log.info("\nMILP SOLVED")
		log.info("\nESTIMATING UPPER BOUND FOR MARKET OVERSEER/TRANSMISSION PLANNING COORDINATOR")
		upperBound = marketoverInstance.getGlobalUpper(LagCombXi, LagCombPi, upperBoundVector, numberOfZones) #MO calculates the global upper bound after every iteration//%%
		log.info("\nESTIMATED UPPER BOUND FOR MARKET OVERSEER/TRANSMISSION PLANNING COORDINATOR")
		log.info("\nESTIMATING LOWER BOUND FOR MARKET OVERSEER/TRANSMISSION PLANNING COORDINATOR")
		lowerBound = marketoverInstance.LBMarketover(solverChoice, LagCombXi, LagCombPi, candLineCount, otherNodeCount) #MO calculates the lower bound for itself after every iteration//%%
		LBIterate += lowerBound
		#lowerBound += marketoverInstance.getGlobalLower(lowerBoundVector, numberOfZones) #MO calculates the global lower bound after every iteration//%%
		log.info("\nESTIMATED LOWER BOUND FOR MARKET OVERSEER/TRANSMISSION PLANNING COORDINATOR")

		consensus = marketoverInstance.getGlobalConsensus() #MO calculates the global consensus after every iteration//%%
		UBItVec.append(upperBound)
		LBItVec.append(LBIterate)
		ratioIterate.append(upperBound/LBIterate)
		globalCons.append(consensus)
		log.info("\nUPDATING THE VALUES OF THE REWARDS/PENALTIES/LAGRANGE MULTIPLIERS/DUAL VARIABLES")
		for i in range(10): #Initialize the Lagrange Multipliers/Dual Variables/Rewards/Penalties to zero for the first iteration
			tempLagMultXi = LagMultXi[i]
			tempLagMultPi = LagMultPi[i]
			LagMultXi[i] = marketoverInstance.rewardPenaltyCont(tempLagMultXi, i, iterCountOut)
			LagMultPi[i] = marketoverInstance.rewardPenaltyInteger(tempLagMultPi, i, iterCountOut)
		#marketoverInstance.rewardPenaltyUpdate(LagMultXi, LagMultPi, i, iterCountOut)
		log.info("\nUPDATED THE VALUES OF THE REWARDS/PENALTIES/LAGRANGE MULTIPLIERS/DUAL VARIABLES")
		#iterationCounter += 1 #Increment the iteration counter before next iteration
	#} while ((((1-abs(lowerBound/upperBound))<=0.05) or ((1-abs(lowerBound/upperBound))>=-0.05)))
	#and (consensus<=0.05)) # while the tolerance is reached//%%
	log.info("\nZonal Calculations of Beliefs about the Investment decision MILP ends")
	ofstream UBIteratesOut("/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputMetrics/UpperBoundIterates.txt", ios::out);
	// exit program if ifstream could not open file
	if ( !UBIteratesOut ) {
		cerr << "\nUpperBoundIterates file could not be opened\n" << endl;
		exit( 1 );
	} // end if
	ofstream LBIteratesOut("/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputMetrics/LowerBoundIterates.txt", ios::out);
	// exit program if ifstream could not open file
	if ( !LBIteratesOut ) {
		cerr << "\nLowerBoundIterates file could not be opened\n" << endl;
		exit( 1 );
	} // end if
	ofstream ratioIteratesOut("/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputMetrics/RatioIterates.txt", ios::out);
	// exit program if ifstream could not open file
	if ( !ratioIteratesOut ) {
		cerr << "\nRatioIterates file could not be opened\n" << endl;
		exit( 1 );
	} // end if
	ofstream globalConsensusOut("/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputMetrics/globalConsensus.txt", ios::out);
	// exit program if ifstream could not open file
	if ( !globalConsensusOut ) {
		cerr << "\nglobalConsensus file could not be opened\n" << endl;
		exit( 1 );
	} // end if
	int countOfIterate = 1; 
	for (UBItIterator = UBItVec.begin(); UBItIterator != UBItVec.end(); ++UBItIterator){
		UBIteratesOut << countOfIterate << "\t" << *UBItIterator << endl;
		++countOfIterate;
	}
	countOfIterate = 1;
	for (LBItIterator = LBItVec.begin(); LBItIterator != LBItVec.end(); ++LBItIterator){
		LBIteratesOut << countOfIterate << "\t" << *LBItIterator << endl;
		++countOfIterate;
	}
	countOfIterate = 1;
	for (ratioIterator = ratioIterate.begin(); ratioIterator != ratioIterate.end(); ++ratioIterator){
		ratioIteratesOut << countOfIterate << "\t" << *ratioIterator << endl;
		++countOfIterate;
	}
	countOfIterate = 1;
	for (globalConsIterator = globalCons.begin(); globalConsIterator != globalCons.end(); ++globalConsIterator){
		globalConsensusOut << countOfIterate << "\t" << *globalConsIterator << endl;
		++countOfIterate;
	}
	cout << endl << "\n*** DISTRIBUTED STOCHASTIC OPTIMIZATION ALGORITHMIC MARKET MECHANISM DESIGN FIRST STAGE ENDS ***\n" << endl << endl;
	marketoverInstance->finDecLineConstr();	// Populate the list of constructed shared candidate lines, according to Stage-I
	int dimAPPLagArray = 0; // Initialize the Dimension of the Lagrange multipliers for APP
	for ( int i = 0; i < numberOfZones; ++i ) { // Each region solves its own MILP optimization problem
		zonalNetVector[i]->setRealizedCLines(*marketoverInstance); // Perform unit commitment for average heat rate objective
		dimAPPLagArray += zonalNetVector[i]->returnMultiplicity(); // 
	} 
	for ( int i = 0; i < numberOfZones; ++i ) {
		zonalNetVector[i]->TestBuiltExternalNodes();
	} 
	double APPLagMultipliers[dimAPPLagArray];
	for (int i = 0; i < dimAPPLagArray; ++i) 
		APPLagMultipliers[i] = 0;
/*
	log.info("\n*** AUXILIARY PROBLEM PRINCIPLE (APP) OPTIMIZATION ALGORITHMIC MARKET MECHANISM DESIGN SECOND STAGE BEGINS ***\n")
	interAngleMessage = [np.zeros(float, numberOfZones)] #Array of the zonal vectors of the intermediate values of Lagrange Multipliers
	zonalGlobRank = [np.zeros(float, numberOfZones)] #Array of the global ranks of the shared nodes for each zone
	ObjIterate = np.zeros(float, numberOfZones) #Array of zonal optimum objectives
	//do {//%%
	for iterCountOut in range(1):
		log.info("\n*** ITERATION {} BEGINS ***\n".format(iterCountOut+1))
		for i in range(numberOfZones): #Each region solves its own MILP optimization problem 
			log.info("\n*** APP QUADRATIC PROGRAMMING FOR ZONE {} BEGINS ***\n".format(i+1))
			log.info("\nZonal Calculations of Beliefs about the generation and flow decision APP-QP begins")
			log.info("\nSOLVING THE OPTIMIZATION SUB-PROBLEM")
			ObjIterate[i]=zonalNetVector[i].APPQPAvgHR(marketoverInstance, APPLagMultipliers, candLineCount, otherNodeCount, environmentGUROBI, iterCountOut) #Perform unit commitment for average heat rate objective
			interAngleMessage[i] = zonalNetVector[i].getZonalDecision();
			zonalGlobRank[i] = zonalNetVector[i].getZonalRanks();
			log.info("\nOPTIMIZATION SUB-PROBLEM SOLVED")
		for i in range(otherNodeCount+1): #Initialize the Lagrange Multipliers/Dual Variables/Rewards/Penalties to zero for the MO
			tempXi = 0.0
			for j in range(numberOfZones):
				tempXi += LagMultXi[j*otherNodeCount+i]
			LagCombXi[i] = tempXi
		for i in range(candLineCount+1): #Initialize the Lagrange Multipliers/Dual Variables/Rewards/Penalties to zero for the MO
			tempPi = 0
			for j in range(numberOfZones):
				tempPi += LagMultPi[j*candLineCount+i]
			LagCombPi[i] = tempPi

		consensus = marketoverInstance.getGlobalConsensus() #MO calculates the global consensus after every iteration//%%
		UBItVec.append(upperBound)
		LBItVec.append(LBIterate)
		ratioIterate.append(upperBound/LBIterate)
		globalCons.append(consensus)
		log.info("\nUPDATING THE VALUES OF THE REWARDS/PENALTIES/LAGRANGE MULTIPLIERS/DUAL VARIABLES")
		for i in range(1000): #Initialize the Lagrange Multipliers/Dual Variables/Rewards/Penalties to zero for the first iteration
			tempLagMultXi = LagMultXi[i]
			tempLagMultPi = LagMultPi[i]
			LagMultXi[i] = marketoverInstance.rewardPenaltyCont(tempLagMultXi, i)
			LagMultPi[i] = marketoverInstance.rewardPenaltyInteger(tempLagMultPi, i)
		log.info("\nUPDATED THE VALUES OF THE REWARDS/PENALTIES/LAGRANGE MULTIPLIERS/DUAL VARIABLES")
		#iterationCounter+=1 #Increment the iteration counter before next iteration
		marketoverInstance.clearVectors() #clear the different interim vectors for making them ready for next iteration
	#while ((((1-abs(lowerBound/upperBound))<=0.05) or ((1-abs(lowerBound/upperBound))>=-0.05)))
	#and (consensus<=0.05)) #while the tolerance is reached//%%
	log.info("\nZonal Calculations of Beliefs about the Investment decision MILP ends")
	ofstream UBIteratesOut("UpperBoundIterates.txt", ios::out);
	// exit program if ifstream could not open file
	if ( !UBIteratesOut ) {
		cerr << "\nUpperBoundIterates file could not be opened\n" << endl;
		exit( 1 );
	} // end if
	ofstream LBIteratesOut("LowerBoundIterates.txt", ios::out);
	// exit program if ifstream could not open file
	if ( !LBIteratesOut ) {
		cerr << "\nLowerBoundIterates file could not be opened\n" << endl;
		exit( 1 );
	} // end if
	ofstream ratioIteratesOut("RatioIterates.txt", ios::out);
	// exit program if ifstream could not open file
	if ( !ratioIteratesOut ) {
		cerr << "\nRatioIterates file could not be opened\n" << endl;
		exit( 1 );
	} // end if
	ofstream globalConsensusOut("globalConsensus.txt", ios::out);
	#exit program if ifstream could not open file
	if ( !globalConsensusOut ):
		log.info("\nglobalConsensus file could not be opened\n")
	countOfIterate = 1 
	for UBItIterator in UBItVec:
		UBIteratesOut << countOfIterate << "\t" << *UBItIterator << endl;
		countOfIterate += 1
	countOfIterate = 1
	for LBItIterator in LBItVec:
		LBIteratesOut << countOfIterate << "\t" << *LBItIterator << endl;
		countOfIterate += 1
	countOfIterate = 1
	for ratioIterator in ratioIterate:
		ratioIteratesOut << countOfIterate << "\t" << *ratioIterator << endl;
		countOfIterate += 1
	countOfIterate = 1
	for globalConsIterator in globalCons:
		globalConsensusOut << countOfIterate << "\t" << globalConsIterator << endl;
		countOfIterate += 1
	log.info("\n*** DISTRIBUTED STOCHASTIC OPTIMIZATION ALGORITHMIC MARKET MECHANISM DESIGN FIRST STAGE ENDS ***\n")
*/


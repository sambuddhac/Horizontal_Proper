#Definition for Marketover class public Member Methods
import numpy as np
from Python_src.nettran import Network
#include <glpk.h> // Includes the GLPK (GNU Linear Programming Kit) header file
#include "gurobi_c++.h"
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
		"""
		clock_t begin = clock() #start the timer
		vector<int>::iterator diffZNIt; // Iterator for diffZoneNodeID
		vector<Node*>::iterator nodeIterator; // Iterator for node objects
		vector<double>::iterator candCapIterator; // Iterator for candidate lines for capacity
		vector<int>::iterator candIterator; // Iterator for candidate lines	
		vector<double>::iterator exsharedCapIterator; // Iterator for shared existing lines for capacity
		vector<int>::iterator exsharedIterator;

		string outSummaryFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputSummaryGLPK/OutSummaryMOGLPK.txt";
		ofstream outPutFile(outSummaryFileName, ios::out); // Create Output File to output the Summary of Results
		if (!outPutFile){
		cerr << "\nCouldn't open the file" << endl;
		exit(1);
		}
		"""
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
					while interNodeRank in self.angleDecIndex[scenPos]: # while all the different positions of this rank hasn't been accounted for
      						pos1 = (self.angleDecIndex[scenPos]).index(interNodeRank) #get the location in the vector, of this rank
						rankPresChecker[pos1] = 1 #Mark the element in the corresponding position of checker vector 1, to indicate this rank has been accounted for 
						interAngleTermVec.append(-(self.phaseAngleDecision[scenPos][pos1])*(LagMultXi[scenPos*self.nodeNumber+(angleDecIndexIterator-1)])) #Calculate cost term
      						pos = std::find(pos + 1, (self.angleDecIndex[scenPos]).index(interNodeRank), ) #Find position of the next occurence of this rank
					smallest_element = *min_element(interAngleTermVec.begin(), interAngleTermVec.end())
					revisedUpperBound += smallest_element
					interAngleTermVec = []
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
			if self.rankPresCheckerInt[trackerInt] == 0: #If this line rank hasn't been already accounted for
				pos = self.lineDecIndex.index(interLineRank) #find the first position of this rank in the vector
				while(pos != self.lineDecIndex[lengthInt-1]): #while all the different positions of this rank hasn't been accounted for
      					pos1 = std::distance(lineDecIndex.begin(), pos) #get the location in the vector, of this rank
					rankPresCheckerInt.at(pos1) = 1 #Mark the element in the corresponding position of checker vector 1, to indicate this rank has been accounted for 
					interLineTermVec.push_back(-(lineInterDecision[pos1])*(LagMultPi[(*lineDecIndexIterator-1)])) #Calculate cost term
      					pos = std::find(pos + 1, lineDecIndex.end(), interLineRank) #Find position of the next occurence of this rank
				double smallest_element = *min_element(interLineTermVec.begin(), interLineTermVec.end())
				revisedUpperBound += smallest_element
				interLineTermVec = []
			trackerInt += 1 #Increment the tracker
		return revisedUpperBound
	#Function getGlobalUpper ends

	def getGlobalLower(regionalLower, numberOfZones): #Function getGlobalLower() returns the global lower bound for the investment coordination problem
		revisedLowerBound = 0 #total lower bound initialized
		for i in range(numberOfZones):
			revisedLowerBound += regionalLower[i]
		return revisedLowerBound
	#Function getGlobalLower ends

	def getGlobalConsensus(): #Function getGlobalConsensus() returns the global consensus for the investment coordination problem
	globConsensus = 0 #total consensus
	#Calculate the MO upper bound, which is the minimum value of the MO objective, calculated with a given mix of the intermediate, shared continuous variable values from different regions/zones
	for scenPos in range(countOfScenarios):
		int interNodeRank; // Intermediate variable for storing the rank of the node 
		vector<int> rankPresChecker; // Vector to check if a particular rank has been accounted for 
		int length = (angleDecIndex[scenPos]).size(); // length of the rank vector
		vector<int>::iterator angleDecIndexIterator; // Iterator for the rank vector
		for (int i = 0; i < length; ++i) // Put as many zeroes on the rankPresChecker vector as is the length of the rank vector 
			rankPresChecker.push_back(0);
		int tracker = 0; // tracker to track the position of the angleDecIndexIterator
		for (angleDecIndexIterator=(angleDecIndex[scenPos]).begin(); angleDecIndexIterator!=(angleDecIndex[scenPos]).end(); ++angleDecIndexIterator) { // Iterate through rank vector
			interNodeRank = (*angleDecIndexIterator); // Store the value of the rank of the present node iterate in the list 
			if ( rankPresChecker.at(tracker) == 0 ) { // If this node rank hasn't been already accounted for
				auto pos = std::find((angleDecIndex[scenPos]).begin(), (angleDecIndex[scenPos]).end(), interNodeRank); // find the first position of this rank in the vector
				while(pos != (angleDecIndex[scenPos]).end()) // while all the different positions of this rank hasn't been accounted for 
    				{
      					auto pos1 = std::distance((angleDecIndex[scenPos]).begin(), pos); // get the location in the vector, of this rank
					rankPresChecker.at(pos1) = 1; // Mark the element in the corresponding position of checker vector 1, to indicate this rank has been accounted for 
					globConsensus += pow((((phaseAngleDecision[scenPos]).at(pos1))-(interimContDecVar[(*angleDecIndexIterator-1)])), 2); // Calculate cost term
      					pos = std::find(pos + 1, (angleDecIndex[scenPos]).end(), interNodeRank); // Find position of the next occurence of this rank
   				}					
			}
			tracker+=1 #Increment the tracker
		}
	}
	// Calculate the MO upper bound, which is the minimum value of the MO objective, calculated with a given mix of the intermediate, shared discrete variable values from different regions/zones
	int interLineRank; // Intermediate variable for storing the rank of the line 
	vector<int> rankPresCheckerInt; // Vector to check if a particular rank has been accounted for 
	int lengthInt = lineDecIndex.size(); // length of the rank vector
	vector<int>::iterator lineDecIndexIterator; // Iterator for the rank vector
	for (int i = 0; i < lengthInt; ++i) // Put as many zeroes on the rankPresChecker vector as is the length of the rank vector 
		rankPresCheckerInt.push_back(0);
	int trackerInt = 0; // tracker to track the position of the lineDecIndexIterator
	for (lineDecIndexIterator=lineDecIndex.begin(); lineDecIndexIterator!=lineDecIndex.end(); ++lineDecIndexIterator) { // Iterate through rank vector
		interLineRank = (*lineDecIndexIterator); // Store the value of the rank of the present line iterate in the list 
		if ( rankPresCheckerInt.at(trackerInt) == 0 ) { // If this line rank hasn't been already accounted for
			auto pos = std::find(lineDecIndex.begin(), lineDecIndex.end(), interLineRank); // find the first position of this rank in the vector
			while(pos != lineDecIndex.end()) // while all the different positions of this rank hasn't been accounted for 
    			{
      				auto pos1 = std::distance(lineDecIndex.begin(), pos); // get the location in the vector, of this rank
				rankPresCheckerInt.at(pos1) = 1; // Mark the element in the corresponding position of checker vector 1, to indicate this rank has been accounted for 
				globConsensus += pow(((lineInterDecision[pos1])-(interimIntDecVar[(*lineDecIndexIterator-1)])), 2); // Calculate cost term
      				pos = std::find(pos + 1, lineDecIndex.end(), interLineRank); // Find position of the next occurence of this rank
   			}					
		}
		++trackerInt; // Increment the tracker
	} 
	return sqrt(globConsensus);
	
} // Function getGlobalConsensus ends

	def finDecLineConstr(): #Final decisions on the construction stauses of candidate lines taking into account the decisions of different zones
{
	int interLineRank; // Intermediate variable for storing the rank of the line 
	vector<int> rankPresCheckerInt; // Vector to check if a particular rank has been accounted for 
	int lengthInt = lineDecIndex.size(); // length of the rank vector
	vector<int>::iterator lineDecIndexIterator; // Iterator for the rank vector
	for (int i = 0; i < lengthInt; ++i) // Put as many zeroes on the rankPresChecker vector as is the length of the rank vector 
		rankPresCheckerInt.push_back(0);
	int trackerInt = 0; // tracker to track the position of the lineDecIndexIterator
	for (lineDecIndexIterator=lineDecIndex.begin(); lineDecIndexIterator!=lineDecIndex.end(); ++lineDecIndexIterator) { // Iterate through rank vector
		interLineRank = (*lineDecIndexIterator); // Store the value of the rank of the present line iterate in the list 
		if ( rankPresCheckerInt.at(trackerInt) == 0 ) { // If this line rank hasn't been already accounted for
			int first = 0; // Initialize a flag to indicate the first occurence of a rank
			int firstInSeries = 0; // flag to indicate if the first verdict for this rank is 1 or 0
			auto pos = std::find(lineDecIndex.begin(), lineDecIndex.end(), interLineRank); // find the first position of this rank in the vector
			while(pos != lineDecIndex.end()) // while all the different positions of this rank hasn't been accounted for 
    			{
      				auto pos1 = std::distance(lineDecIndex.begin(), pos); // get the location in the vector, of this rank
				rankPresCheckerInt.at(pos1) = 1; // Mark the element in the corresponding position of checker vector 1, to indicate this rank has been accounted for 
				if((lineInterDecision[pos1]==1) && (first == 0)) {// Check if the decision of the different subnets is 1
					constructedRanks.push_back(interLineRank);
					firstInSeries = 1; // Set the flag to indicate that the first verdict for this rank is 1
				}
				else if((lineInterDecision[pos1]==0) && (first != 0) && (firstInSeries == 1))
					constructedRanks.pop_back(); 
      				pos = std::find(pos + 1, lineDecIndex.end(), interLineRank); // Find position of the next occurence of this rank
				++first; // Increment the flag to indicate the recurrence of the rank
   			}					
		}
		++trackerInt; // Increment the tracker
	} 
	vector<int>::iterator construcIterator;
	for (construcIterator=constructedRanks.begin(); construcIterator!=constructedRanks.end(); ++construcIterator) {	
		//cout << " First message " << endl;
		//cout << "Constructed shared candidate lines are " << *construcIterator << endl;
		//cout << " Second message " << endl;
	}
	//cout << "end of finDecLineConstr " << endl; 	
}

	def scanBuiltLinesList(int rankGlob): #Scans the list of built candidate lines for the input argument of global rank to check if this line is built
{
	auto pos = std::find(constructedRanks.begin(), constructedRanks.end(), rankGlob); // find the first position of this rank in the vector
	if(pos != constructedRanks.end()) // if the rank is found in the list
    		return 1; // return 1
	else 
		return 0; // Otherwise, return 0
}

	def populateLineDec(int interimDecBin, int zonalIndex, int globalRank):#Method to pass the intermediate message for the zonal line building decision to the MO
{
	lineInterDecision.push_back(interimDecBin);
	lineDecIndex.push_back(globalRank);
	zonalIndVectorInt.push_back(zonalIndex);
	lineInterDecVec.at(zonalIndex*candLineNumber+globalRank)=interimDecBin;
}

	def populateAngleDec(double interimDecAngle, int zonalIndex, int scenario, int globalRank): #Method to pass the intermediate message for the zonal node angle decision to the MO
{
	(phaseAngleDecision[scenario]).push_back(interimDecAngle);
	(angleDecIndex[scenario]).push_back(globalRank);
	(zonalIndVectorCont[scenario]).push_back(zonalIndex);
	phAngDecVec.at(scenario*zoneNumber*nodeNumber+zonalIndex*nodeNumber+globalRank)=interimDecAngle;
}

	def rewardPenaltyUpdate(double lagrangeMultXi[], double lagrangeMultPi[], int matrixIndex, int iteration): #Function getGlobalUpper() returns the global upper bound for the investment coordination problem
{
"""	for (int scenPos = 0; scenPos < countOfScenarios; ++scenPos) {
		for (int zonalIndex=0;zonalIndex<zoneNumber;++zonalIndex) {
			for (int globalRank=1;globalRank<=nodeNumber;++globalRank) {			
				if(compareBasis == 1) {
					lagrangeMultXi[scenPos*zoneNumber*nodeNumber+zonalIndex*nodeNumber+globalRank] += (((phAngDecVec).at(scenPos*zoneNumber*nodeNumber+zonalIndex*nodeNumber+globalRank))-(interimContDecVar[scenPos*nodeNumber+globalRank]));
					//cout << "From MO Continuous Lagrange Multiplier update: " << "Angle from Zones. Tracker: " << tracker << " Scenario: " << scenPos << " Zonal Angle value " << ((phaseAngleDecision[scenPos]).at(tracker)) << " MO Angle index: " << scenPos*nodeNumber+(*angleDecIndexIterator-1) << " MO Angle value: " << (interimContDecVar[scenPos*nodeNumber+(*angleDecIndexIterator-1)]) << " Updated Lagrange Multiplier: " << lagrangeMult << endl; 
				}
				else {
					lagrangeMultXi[scenPos*zoneNumber*nodeNumber+zonalIndex*nodeNumber+globalRank] += (((phAngDecVec).at(scenPos*zoneNumber*nodeNumber+zonalIndex*nodeNumber+globalRank))-(interimContDecVarPrev[scenPos*nodeNumber+globalRank]));
					//cout << "From MO Continuous Lagrange Multiplier update: " << "Angle from Zones. Tracker: " << tracker << " Scenario: " << scenPos << " Zonal Angle value " << ((phaseAngleDecision[scenPos]).at(tracker)) << " MO Angle index: " << scenPos*nodeNumber+(*angleDecIndexIterator-1) << " MO Angle value: " << (interimContDecVar[scenPos*nodeNumber+(*angleDecIndexIterator-1)]) << " Updated Lagrange Multiplier: " << lagrangeMult << endl; 
				}
			}
		}
	}
	for (int zonalIndex=0;zonalIndex<zoneNumber;++zonalIndex) {
			for (int globalRank=1;globalRank<=candLineNumber;++globalRank) {
				if (compareBasis == 1) {		
					lagrangeMultPi[zonalIndex*candLineNumber+globalRank] += (lineInterDecVec.at(zonalIndex*candLineNumber+globalRank)-interimIntDecVar[globalRank]);
					//cout << "From MO Integer Lagrange Multiplier update: " << "Line Decision from Zones. Tracker: " << tracker << " Zonal Decision value " << lineInterDecision[tracker] << " MO decision index: " << (*lineDecIndexIterator-1) << " MO Decision value: " << interimIntDecVar[(*lineDecIndexIterator-1)] << " Updated Lagrange Multiplier: " << lagrangeMult << endl; 
				}
				else {
					lagrangeMultPi[zonalIndex*candLineNumber+globalRank] += (lineInterDecVecat(zonalIndex*candLineNumber+globalRank)-interimIntDecVarPrev[(*globalRank]);
					//cout << "From MO Integer Lagrange Multiplier update: " << "Line Decision from Zones. Tracker: " << tracker << " Zonal Decision value " << lineInterDecision[tracker] << " MO decision index: " << (*lineDecIndexIterator-1) << " MO Decision value: " << interimIntDecVar[(*lineDecIndexIterator-1)] << " Updated Lagrange Multiplier: " << lagrangeMult << endl; 
				}
			}
	}	"""
}

	def rewardPenaltyCont(double lagrangeMult, int matrixIndex, int iteration): #Function getGlobalUpper() returns the global upper bound for the investment coordination problem
{
	//cout << "\nITERATION COUNT : " << iteration << endl;
	for (int scenPos = 0; scenPos < countOfScenarios; ++scenPos) {
		int tracker = 0;
		vector<int>::iterator angleDecIndexIterator;
		for (angleDecIndexIterator=(angleDecIndex[scenPos]).begin();angleDecIndexIterator!=(angleDecIndex[scenPos]).end();++angleDecIndexIterator) {
			if (matrixIndex==(scenPos*zoneNumber*nodeNumber+(zonalIndVectorCont[scenPos]).at(tracker)-1)*nodeNumber+(*angleDecIndexIterator)) {			
				if(compareBasis == 1) {
					lagrangeMult += (((phaseAngleDecision[scenPos]).at(tracker))-(interimContDecVar[scenPos*nodeNumber+(*angleDecIndexIterator-1)]));
					//cout << "From MO Continuous Lagrange Multiplier update: " << "Angle from Zones. Tracker: " << tracker << " Scenario: " << scenPos << " Zonal Angle value " << ((phaseAngleDecision[scenPos]).at(tracker)) << " MO Angle index: " << scenPos*nodeNumber+(*angleDecIndexIterator-1) << " MO Angle value: " << (interimContDecVar[scenPos*nodeNumber+(*angleDecIndexIterator-1)]) << " Updated Lagrange Multiplier: " << lagrangeMult << endl; 
				}
				else {
					lagrangeMult += (((phaseAngleDecision[scenPos]).at(tracker))-(interimContDecVarPrev[scenPos*nodeNumber+(*angleDecIndexIterator-1)]));
					//cout << "From MO Continuous Lagrange Multiplier update: " << "Angle from Zones. Tracker: " << tracker << " Scenario: " << scenPos << " Zonal Angle value " << ((phaseAngleDecision[scenPos]).at(tracker)) << " MO Angle index: " << scenPos*nodeNumber+(*angleDecIndexIterator-1) << " MO Angle value: " << (interimContDecVar[scenPos*nodeNumber+(*angleDecIndexIterator-1)]) << " Updated Lagrange Multiplier: " << lagrangeMult << endl; 
				}
			}
			++tracker;
		}
		//cout << "Tracker " << tracker << endl;
	}
	return lagrangeMult;	
} // Function getGlobalUpper ends

	def rewardPenaltyInteger(lagrangeMult, matrixIndex, iteration): #Function getGlobalUpper() returns the global upper bound for the investment coordination problem
		#log.info("\nITERATION COUNT : {}".format(iteration))
		tracker = 0
		for lineDecIndexIterator in self.lineDecIndex:
			if matrixIndex==(zonalIndVectorInt[tracker]-1)*candLineNumber+(*lineDecIndexIterator):
				if self.compareBasis == 1:	
					lagrangeMult += (lineInterDecision[tracker]-interimIntDecVar[(*lineDecIndexIterator-1)])
					#cout << "From MO Integer Lagrange Multiplier update: " << "Line Decision from Zones. Tracker: " << tracker << " Zonal Decision value " << lineInterDecision[tracker] << " MO decision index: " << (*lineDecIndexIterator-1) << " MO Decision value: " << interimIntDecVar[(*lineDecIndexIterator-1)] << " Updated Lagrange Multiplier: " << lagrangeMult << endl
				else:
					lagrangeMult += (lineInterDecision[tracker]-interimIntDecVarPrev[(*lineDecIndexIterator-1)])
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

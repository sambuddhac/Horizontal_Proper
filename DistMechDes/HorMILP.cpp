// HorMILP.cpp : Defines the entry point for Stage-I and Stage-II of the Horizontal Investment Coordination MILP Market Mechanism Design Simulation application.
// Main Method for running the Horizontal Investment Coordination Stage-I and Stage-II MILP Market Mechanism Design Simulation based on Distributed Stochastic UC and APP
#include <iostream>
#include <cstring>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "gurobi_c++.h" // includes definition of the GUROBI solver header file
#include "nettran.h" // includes definition of Nettran class for sub-networks for MILP based investment decisions for building transmission lines
#include "marketOverseer.h" // includes the definition of marketOverseer class for attaining consensus//%%
using namespace std;

int main() // Main method begins program execution
{
	int curveChoice; // Number to indicate the type of Objective function among average heat rate, piecewise linear, or polynomial
	/* Future Work
	// Choose the type of objective function
	//cout << "\nChoose the type of Objective Function: 1 for Average Heat Rate, 2 for Piecewise Linear, and 3 for Polynomial (Cubic) Convex function" << endl;
	//cin >> curveChoice; */
	int basisForComparison; // While updating the Lagrange multipliers, 0 when the current zonal updates are compared to previous iteration MO update;else, 1
	int solverChoice; // Choice of the solver
	curveChoice = 1; //Assume Average Heat Rate for now
	int systemChoice; // Choice of the system to be simulated
	cout << "\nChoose the type of System to be simulated: 1 for Simple two bus/two region, 2 for system combined of IEEE 14, 30, and 5 node systems" << endl;
	cin >> systemChoice;
	cout << "\nChoose, for updating the Lagrange multipliers, 0 when the current zonal updates are compared to previous iteration MO update;else, 1" << endl;
	cin >> basisForComparison;
	// Read the master zones file, for deciding upon which other files to read for building the model
	int numberOfZones; // Number of zones between which horizontal investment coordination for transmission lines to be built is considered
	int numberOfFields; // Number of rows or individual file types for each of the zones
	int lpMethodChoice; // Simplex or interior method algorithm for MILP
	string inputMasterFile;
	if (systemChoice==1)
		inputMasterFile = "masterZonesSummary.txt";
	else
		inputMasterFile = "masterZonesSummaryRevised.txt";
	double UBIterate = 0.0; // Initial value for the global upper bound iterates at the end of every iteration
	double LBIterate = 0.0; // Initial value for the global lower bound iterates at the end of every iteration
	vector <double> UBItVec; // Vector for storing the values of the global upper bound iterates for every iteration
	vector <double> LBItVec; // Vector for storing the values of the global lower bound iterates for every iteration
	vector <double> ratioIterate; // Vector for storing the values of UB/LB iterates for every iteration
	vector <double> globalCons; // Vector for storing the values of global consensus for every iteration
	vector <double>::iterator UBItIterator; // Vector iterator for iterating the UBItVec
	vector <double>::iterator LBItIterator; // Vector iterator for iterating the LBItVec
	vector <double>::iterator ratioIterator; // Vector iterator for iterating the ratioIterate vector
	vector <double>::iterator globalConsIterator; // Vector iterator for iterating the globalCons vector
	ifstream zoneSummaryFile( inputMasterFile, ios::in ); // ifstream constructor opens the master zones summary file
	stringstream buffer; // stringstream object to store the read information from the summary file
	// exit program if ifstream could not open file
	if ( !zoneSummaryFile ) {
		cerr << "\nMaster file for Zones Summary could not be opened\n" << endl;
		exit( 1 );
	} // end if

	//zoneSummaryFile >> numberOfZones >> numberOfFields; // get the number of zones and the number of fields: Future expansion
	cout << "\nEnter the number of zones" << endl;
	cin >> numberOfZones; // User input the number of zones/regions
	cout << "\nChoose either the GLPK (1) or GUROBI (2) as the Solver. " << endl;
	cin >> solverChoice;
	if (solverChoice==1) {
		cout << "\nChoose either the Simplex LP Rlaxation (1) or Interior Point Method LP Relaxation (2) as the method to provide the initial basis to the Mixed Integer Unit Commitment Problem. " << endl;
	cin >> lpMethodChoice;
	}
	GRBEnv* environmentGUROBI = new GRBEnv("GUROBILogFile.log"); // GUROBI Environment object for storing the different optimization models
	numberOfFields = 7; // Number of fields
   	buffer << zoneSummaryFile.rdbuf(); // reads the data in the summary file 
   	string test = buffer.str(); // Extract the strings from the buffer to "test"

   	//create variables that will act as "cursors". we'll take everything between them.
   	size_t pos1 = 0;
   	size_t pos2;
   	//create the array to store the strings.
   	string str[numberOfFields*numberOfZones];
	//Read the summary input file
   	for ( int i = 0; i < numberOfFields; ++i ) {
		for ( int j = 0; j < numberOfZones; ++j ) {
			if (j==numberOfZones-1){
				pos2 = test.find("\n", pos1); //search for the bar "\n". pos2 will be where the bar was found.
        			str[i*numberOfZones+j] = test.substr(pos1, (pos2-pos1)); //make a substring, wich is nothing more 
                                              //than a copy of a fragment of the big string.
        			pos1 = pos2+1; // sets pos1 to the next character after pos2. 
    			}
			else {
        			pos2 = test.find(" ", pos1); //search for the bar " ". pos2 will be where the bar was found.
        			str[i*numberOfZones+j] = test.substr(pos1, (pos2-pos1)); //make a substring, wich is nothing more 
                                              //than a copy of a fragment of the big string.
        			pos1 = pos2+1; // sets pos1 to the next character after pos2. 
    			}
    		}
 	}	

	vector< Nettran* > zonalNetVector; // Vector of zonal network objects
	cout << endl << "\n*** NETWORK INITIALIZATION STAGE BEGINS ***\n" << endl << endl;
	double upperBoundVector[numberOfZones]; // Vector of upper bounds by iteration
	double lowerBoundVector[numberOfZones]; // Vector of lower bounds by iteration
	for ( int i = 0; i < numberOfZones; ++i ) {
		Nettran *nettranInstance = new Nettran( str, (i+1), numberOfZones, curveChoice, lpMethodChoice ); // create the network instances for the different zones
		zonalNetVector.push_back( nettranInstance ); // push to the vector of zonal networks
	}
	cout << "\n*** NETWORK INITIALIZATION STAGE ENDS: ZONAL SUB-NETWORKS CREATED ***\n" << endl;
/*
	double LagMultXi[1000]; // Lagrange Multipliers/Dual Variables/Rewards/Penalties corresponding to bus voltage phase angle consensus
	double LagMultPi[1000]; // Lagrange Multipliers/Dual Variables/Rewards/Penalties corresponding to line building decision integer variable consensus
	double LagCombXi[1000]; // Lagrange Multipliers/Dual Variables/Rewards/Penalties corresponding to bus voltage phase angle consensus for MO
	double LagCombPi[1000]; // Lagrange Multipliers/Dual Variables/Rewards/Penalties corresponding to line building decision integer variable consensus MO
	for ( int i = 0; i < 1000; ++i ) { // Initialize the Lagrange Multipliers/Dual Variables/Rewards/Penalties to zero for the first iteration
		LagMultXi[i] = 0.0;
		LagMultPi[i] = 0.0;
	}
*/
	// Stitching the different nodes of shared existing & shared candidate lines and zones in one list
	vector<int> sharedNodeList; // List of all the node IDs of the shared existing and/or shared candidate transmission lines
	vector<int> sharedZoneList; // List of all the corresponding zone IDs of the respective shared line end nodes
	vector<int> sharedGlobalNodeList; // Global ranking of the shared nodes		
	sharedNodeList.push_back(0); // Initialize the list of shared line end node ID's so that it's not empty
	sharedZoneList.push_back(0); // Initialize the list of shared line end zone ID's so that it's not empty
	sharedGlobalNodeList.push_back(0); // Initialize the list of shared global line end zone ID's so that it's not empty
	vector<int>::iterator diffZNIt; // Iterator for sharedNodeList
	int containsFlag = 0; // Initializes the flag to test if the list so far contains a particular node, to avoid duplication
	int otherNodeCount = 0; // Initializes the count of total number of shared nodes
	int sharedZone = 0; // Temporary integer for storing the shared zone index
	int sharedNode = 0; // Temporary integer for storing the shared node index
	vector<int> sharedSESerList; // List of all the unique global serial numbers of the shared existing lines
	sharedSESerList.push_back(0); // Initialize the list of shared existing line global serial numbers so that it's not empty
	vector<int> sharedSEFromList; // List of the From node ranks of the shared existing lines
	vector<int> sharedSEToList; // List of the To node ranks of the shared existing lines
	vector<double> sharedSEReactList; // List of all the reactances of the shared existing lines
	vector<double> sharedSECapList; // List of all the flow capacities of the shared existing lines
	vector<int>::iterator SESerList; // Iterator for sharedSESerList	
	int containsSEFlag = 0; // Initializes the flag to test if the list so far contains the shared existing line, to avoid duplication
	int SELineCount = 0; // Initializes the count of total number of shared existing lines
	int sharedELine = 0; // Temporary integer for storing the shared existing line global serial number
	int SEFromNode; // Node ID of the from node of the shared existing line
	int SEFromZone; // Zone ID of the from zone of the shared existing line
	int SEToNode; // Node ID of the to node of the shared existing line
	int SEToZone; // Zone ID of the to zone of the shared existing line
	int SEFromRank; // Node Rank of the from node of the shared existing line
	int SEToRank; // Node Rank of the to node of the shared existing line
	double SEImp; // Temporary variable for storing the reactance of the shared existing line
	double SECap; // Temporary variable for storing the MW flow limit of the shared existing line
	vector<int> sharedCandSerList; // List of all the unique global serial numbers of the shared candidate lines
	sharedCandSerList.push_back(0); // Initialize the list of shared candidate line global serial numbers so that it's not empty
	vector<int> sharedCandFromList; // List of the From node ranks of the shared candidate lines
	vector<int> sharedCandToList; // List of the To node ranks of the shared candidate lines
	vector<double> sharedCandReactList; // List of all the reactances of the shared candidate lines
	vector<double> sharedCandCapList; // List of all the flow capacities of the shared candidate lines
	vector<int>::iterator CandSerList; // Iterator for sharedCandSerList
	int containsCandFlag = 0; // Initializes the flag to test if the list so far contains the shared candidate line, to avoid duplication
	int candLineCount = 0; // Initializes the count of total number of shared candidate lines
	int sharedCandLine = 0; // Temporary integer for storing the shared candidate line global serial number
	int CandFromNode; // Node ID of the from node of the shared candidate line
	int CandFromZone; // Zone ID of the from zone of the shared candidate line
	int CandToNode; // Node ID of the to node of the shared candidate line
	int CandToZone; // Zone ID of the to zone of the shared candidate line
	int CandFromRank; // Node Rank of the from node of the shared candidate line
	int CandToRank; // Node Rank of the to node of the shared candidate line
	double CandImp; // Temporary variable for storing the reactance of the shared candidate line
	double CandCap; // Temporary variable for storing the MW flow limit of the shared candidate line
	// Stitching together the shared line end nodes from different zones into one master list
	//cout << "Stitched set of nodes at the ends of shared existing and shared candidate lines with global rankings" << endl;
	for ( int i = 0; i < numberOfZones; ++i ) { // Run the loop on all the zones
		//cout << "Zone Considered is " << i << endl;
		sharedZone = 0; // Temporary integer for storing the shared zone index
		sharedNode = 0; // Temporary integer for storing the shared node index
		for ( int sharedNodeCounter = 1; ((sharedZone != -1) && (sharedNode != -1)); ++sharedNodeCounter) { // Run the loop on shared lists
			sharedZone = zonalNetVector[i]->getConnZone(sharedNodeCounter); // Gets the zone ID of the outer zonal node
			sharedNode = zonalNetVector[i]->getConnNode(sharedNodeCounter); // Gets the node ID of the outer zonal node
			int indCount = 0; // Initialize a counter for tracking the position in the vector of the iterator
			if ((sharedZone != -1) && (sharedNode != -1)) { // as long as end of the zonal vectors are not reached
				for (diffZNIt = sharedNodeList.begin(); diffZNIt != sharedNodeList.end(); ++diffZNIt) { // Iterate over the vector of connected external zone nodes
					if ((*diffZNIt == sharedNode) && (sharedZoneList[indCount] == sharedZone)) { // Check whether the other-zone node is already present in the list
						containsFlag = 1;  // Loop through the diffZoneNodeID and diffZoneID vectors, to see if the other zone-node is already present or not, if yes, set the containsFlag to 1
					}
					++indCount; // Increment the counter
				}
				if (containsFlag == 0) { // If the diffZoneNodeID and diffZoneID vectors do not contain the other zone-node, then push the node and the other zone index in the respective vectors
					sharedNodeList.push_back(sharedNode);
					sharedZoneList.push_back(sharedZone);
					//cout << "Local Node Index: " << sharedNode << "\t";
					//cout << "Zone Index: " << sharedZone << "\t";
					++otherNodeCount; // Increment the counter to account for the total number of other-zone nodes
					//cout << "Rank: " << otherNodeCount << endl;
					sharedGlobalNodeList.push_back(otherNodeCount);
					
				} 			
				containsFlag = 0; // Reset the containsFlag for matching the next item
			}
		}
	}
	// Stitching together the shared existing line list from different zones into one master list
	for ( int i = 0; i < numberOfZones; ++i ) { // Run the loop on all the zones
		sharedELine = 0; // Temporary integer for storing the shared existing line global serial number
		for ( int sharedELineCounter = 0; (sharedELine != -1); ++sharedELineCounter) { // Run the loop on shared existing lines
			sharedELine = zonalNetVector[i]->getSESerial(sharedELineCounter); // Gets the global serial number of the shared existing line
			SEFromNode = zonalNetVector[i]->getSEFromNode(sharedELineCounter); // Gets the from node of the shared existing line
			SEFromZone = zonalNetVector[i]->getSEFromZone(sharedELineCounter); // Gets the from zone of the shared existing line
			SEToNode = zonalNetVector[i]->getSEToNode(sharedELineCounter); // Gets the to node of the shared existing line
			SEToZone = zonalNetVector[i]->getSEToZone(sharedELineCounter); // Gets the to zone of the shared existing line
			SEImp = zonalNetVector[i]->getSEReactance(sharedELineCounter); // Gets the reactance of the shared existing line
			SECap = zonalNetVector[i]->getSECapacity(sharedELineCounter); // Gets the MW flow limit of the shared existing line
			if (sharedELine != -1) { // as long as end of the zonal shared existing line vector is not reached
				for (SESerList = sharedSESerList.begin(); SESerList != sharedSESerList.end(); ++SESerList) { // Iterate over the list of SELine
					if (*SESerList == sharedELine) { // Check whether the SE line is already present in the list
						containsSEFlag = 1;  // Loop through the diffZoneNodeID and diffZoneID vectors, to see if the other zone-node is already present or not, if yes, set the containsFlag to 1
						int indCount = 0; // Initialize a counter for tracking the position in the vector of the iterator
						for (diffZNIt = sharedNodeList.begin(); diffZNIt != sharedNodeList.end(); ++diffZNIt) {
							if ((*diffZNIt == SEFromNode) && (sharedZoneList[indCount] == SEFromZone)) { // Check whether the from node is already present in the list
								SEFromRank = indCount; // if yes, get the rank of the from node from the list (incrementing indCount will point to the rank, since sharedNodeList and sharedZoneList both have a dummy 0 as their fist element
								if (i+1 == SEFromZone) // if the present zone has the from node
									zonalNetVector[i]->setSEFromRank(sharedELineCounter, SEFromRank); // assign the rank to the internal from node
								else // else
									zonalNetVector[i]->setSEFromRankConn(sharedELineCounter, SEFromRank); // assign the rank of the external from node and also store the rank to the list of rank of connected nodes to the internal to node
							}
							if ((*diffZNIt == SEToNode) && (sharedZoneList[indCount] == SEToZone)) { // Check whether the to node is already present in the list
								SEToRank = indCount; // if yes, get the rank of the to node from the list (incrementing indCount will point to the rank, since sharedNodeList and sharedZoneList both have a dummy 0 as their fist element
								if (i+1 == SEToZone) // if the present zone has the to node
									zonalNetVector[i]->setSEToRank(sharedELineCounter, SEToRank); // assign the rank to the internal to node
								else // else
									zonalNetVector[i]->setSEToRankConn(sharedELineCounter, SEToRank); // assign the rank of the external to node and also store the rank to the list of rank of connected nodes to the internal to node
							}
							++indCount; // Increment the counter
						}
					}
				}
				if (containsSEFlag == 0) { // If the sharedSESerList vector does not contain the SE line, then push the SE line in the vector
					//cout << "From node of SE line " << sharedELine << " is " << SEFromNode << " From zone is " << SEFromZone << " To node is " << SEToNode << " To zone is " << SEToZone << endl;
					sharedSESerList.push_back(sharedELine);
					int indCount = 0; // Initialize a counter for tracking the position in the vector of the iterator
					for (diffZNIt = sharedNodeList.begin(); diffZNIt != sharedNodeList.end(); ++diffZNIt) {
						if ((*diffZNIt == SEFromNode) && (sharedZoneList[indCount] == SEFromZone)) { // Check whether the from node is already present in the list
							SEFromRank = indCount; // if yes, get the rank of the from node from the list (incrementing indCount will point to the rank, since sharedNodeList and sharedZoneList both have a dummy 0 as their fist element
							if (i+1 == SEFromZone) // if the present zone has the from node
								zonalNetVector[i]->setSEFromRank(sharedELineCounter, SEFromRank); // assign the rank to the internal from node
							else // else
								zonalNetVector[i]->setSEFromRankConn(sharedELineCounter, SEFromRank); // assign the rank of the external from node and also store the rank to the list of rank of connected nodes to the internal to node
							sharedSEFromList.push_back(SEFromRank);
							//cout << "From node of SE line " << sharedELine << " is " << SEFromNode << " From zone is " << SEFromZone << " and from rank is " << SEFromRank << endl;
						}
						if ((*diffZNIt == SEToNode) && (sharedZoneList[indCount] == SEToZone)) { // Check whether the to node is already present in the list
							SEToRank = indCount; // if yes, get the rank of the to node from the list (incrementing indCount will point to the rank, since sharedNodeList and sharedZoneList both have a dummy 0 as their fist element
							if (i+1 == SEToZone) // if the present zone has the to node
								zonalNetVector[i]->setSEToRank(sharedELineCounter, SEToRank); // assign the rank to the internal to node
							else // else
								zonalNetVector[i]->setSEToRankConn(sharedELineCounter, SEToRank); // assign the rank of the external to node and also store the rank to the list of rank of connected nodes to the internal to node
							sharedSEToList.push_back(SEToRank);
							//cout << "To node of SE line " << sharedELine << " is " << SEToNode << " To zone is " << SEToZone << " and to rank is " << SEToRank << endl;
						}
						++indCount; // Increment the counter
					}				
					sharedSEReactList.push_back(SEImp);
					sharedSECapList.push_back(SECap);
					++SELineCount; // Increment the counter to account for the total number of SE Lines
				} 			
				containsSEFlag = 0; // Reset the containsFlag for matching the next item
			}
		}
	}
	// Stitching together the shared candidate line list from different zones into one master list
	for ( int i = 0; i < numberOfZones; ++i ) { // Run the loop on all the zones
		sharedCandLine = 0; // Temporary integer for storing the shared candidate line global serial number
		for ( int candLineCounter = 0; (sharedCandLine != -1); ++candLineCounter) { // Run the loop on shared candidate lines
			sharedCandLine = zonalNetVector[i]->getCandSerial(candLineCounter); // Gets the global serial number of the shared candidate line
			CandFromNode = zonalNetVector[i]->getCandFromNode(candLineCounter); // Gets the from node of the shared candidate line
			CandFromZone = zonalNetVector[i]->getCandFromZone(candLineCounter); // Gets the from zone of the shared candidate line
			CandToNode = zonalNetVector[i]->getCandToNode(candLineCounter); // Gets the to node of the shared candidate line
			CandToZone = zonalNetVector[i]->getCandToZone(candLineCounter); // Gets the to zone of the shared candidate line
			CandImp = zonalNetVector[i]->getCandReactance(candLineCounter); // Gets the reactance of the shared candidate line
			CandCap = zonalNetVector[i]->getCandCapacity(candLineCounter); // Gets the MW flow limit of the shared candidate line
			if (sharedCandLine != -1) { // as long as end of the zonal candidate line vector is not reached
				int candGlobalRank = 0; // Global ranking of the candidate line in the list of shared candidate lines
				for (CandSerList = sharedCandSerList.begin(); CandSerList != sharedCandSerList.end(); ++CandSerList) { // Iterate over candline
					if (*CandSerList == sharedCandLine) { // Check whether the Cand line is already present in the list
						containsCandFlag = 1;  // Loop through the diffZoneNodeID and diffZoneID vectors, to see if the other zone-node is already present or not, if yes, set the containsFlag to 1 
						zonalNetVector[i]->assignCandGlobalRank(candLineCounter, candGlobalRank); // assigns the global rank of the shared candidate line
						int indCount = 0; // Initialize a counter for tracking the position in the vector of the iterator
						for (diffZNIt = sharedNodeList.begin(); diffZNIt != sharedNodeList.end(); ++diffZNIt) {
							if ((*diffZNIt == CandFromNode) && (sharedZoneList[indCount] == CandFromZone)) { // Check whether the from node is already present in the list
								CandFromRank = indCount; // if yes, get the rank of the from node from the list (incrementing indCount will point to the rank, since sharedNodeList and sharedZoneList both have a dummy 0 as their fist element
								if (i+1 == CandFromZone) // if the present zone has the from node
									zonalNetVector[i]->setCandFromRank(candLineCounter, CandFromRank); // assign the rank to the internal from node
								else // else
									zonalNetVector[i]->setCandFromRankConn(candLineCounter, CandFromRank); // assign the rank of the external from node and also store the rank to the list of rank of connected nodes to the internal to node
							}
							if ((*diffZNIt == CandToNode) && (sharedZoneList[indCount] == CandToZone)) { // Check whether the to node is already present in the list
								CandToRank = indCount; // if yes, get the rank of the to node from the list (incrementing indCount will point to the rank, since sharedNodeList and sharedZoneList both have a dummy 0 as their fist element
								if (i+1 == CandToZone) // if the present zone has the to node
									zonalNetVector[i]->setCandToRank(candLineCounter, CandToRank); // assign the rank to the internal to node
								else // else
									zonalNetVector[i]->setCandToRankConn(candLineCounter, CandToRank); // assign the rank of the external to node and also store the rank to the list of rank of connected nodes to the internal to node
							}
							++indCount; // Increment the counter
						}
					}
					++candGlobalRank; // Increment the shared candidate line rank
				}
				if (containsCandFlag == 0) { // If the diffZoneNodeID and diffZoneID vectors do not contain the other zone-node, then push the node and the other zone index in the respective vectors
					//cout << "From node of candidate line " << sharedCandLine << " is " << CandFromNode << " From zone is " << CandFromZone << " To node is " << CandToNode << " To zone is " << CandToZone << endl;
					sharedCandSerList.push_back(sharedCandLine);
					zonalNetVector[i]->assignCandGlobalRank(candLineCounter, candLineCount+1); // assigns the global rank of the shared candidate line
					int indCount = 0; // Initialize a counter for tracking the position in the vector of the iterator
					for (diffZNIt = sharedNodeList.begin(); diffZNIt != sharedNodeList.end(); ++diffZNIt) {
						if ((*diffZNIt == CandFromNode) && (sharedZoneList[indCount] == CandFromZone)) { // Check whether the from node is already present in the list
							CandFromRank = indCount; // if yes, get the rank of the from node from the list (incrementing indCount will point to the rank, since sharedNodeList and sharedZoneList both have a dummy 0 as their fist element
							if (i+1 == CandFromZone) // if the present zone has the from node
								zonalNetVector[i]->setCandFromRank(candLineCounter, CandFromRank); // assign the rank to the internal from node
							else // else
								zonalNetVector[i]->setCandFromRankConn(candLineCounter, CandFromRank); // assign the rank of the external from node and also store the rank to the list of rank of connected nodes to the internal to node
							sharedCandFromList.push_back(CandFromRank);
							//cout << "From node of candidate line " << sharedCandLine << " is " << CandFromNode << " From zone is " << CandFromZone << " and from rank is " << CandFromRank << endl;
						}
						if ((*diffZNIt == CandToNode) && (sharedZoneList[indCount] == CandToZone)) { // Check whether the to node is already present in the list
							CandToRank = indCount; // if yes, get the rank of the to node from the list (incrementing indCount will point to the rank, since sharedNodeList and sharedZoneList both have a dummy 0 as their fist element
							if (i+1 == CandToZone) // if the present zone has the to node
								zonalNetVector[i]->setCandToRank(candLineCounter, CandToRank); // assign the rank to the internal to node
							else // else
								zonalNetVector[i]->setCandToRankConn(candLineCounter, CandToRank); // assign the rank of the external to node and also store the rank to the list of rank of connected nodes to the internal to node
							sharedCandToList.push_back(CandToRank);
							//cout << "To node of candidate line " << sharedCandLine << " is " << CandToNode << " To zone is " << CandToZone << " and to rank is " << CandToRank << endl;
						}
						++indCount; // Increment the counter
					}				
					sharedCandReactList.push_back(CandImp);
					sharedCandCapList.push_back(CandCap);
					++candLineCount; // Increment the counter to account for the total number of SE Lines
				} 			
				containsCandFlag = 0; // Reset the containsFlag for matching the next item
			}
		}
	}
	
	cout << endl << "\n*** MARKET OVERSEER INITIALIZATION STAGE BEGINS ***\n" << endl << endl;
	int scenarios = zonalNetVector[0]->getNumberOfScenarios(); // Returns the number of scenarios
	Marketover *marketoverInstance = new Marketover( numberOfZones, otherNodeCount, SELineCount, candLineCount, &zonalNetVector, &sharedNodeList, &sharedGlobalNodeList, &sharedZoneList, &sharedSESerList, &sharedSEFromList, &sharedSEToList, &sharedSEReactList, &sharedSECapList, &sharedCandSerList, &sharedCandFromList, &sharedCandToList, &sharedCandReactList, &sharedCandCapList, lpMethodChoice, scenarios, basisForComparison ); // create the market overseer instance//%%
	cout << "\n*** MARKET OVERSEER INITIALIZATION STAGE ENDS ***\n" << endl;
	cout << "\nShared Node List" << endl;
	vector<int>::iterator sharedNodeIterator;
	vector<int>::iterator sharedCandLineIterator;
	for (sharedNodeIterator=sharedNodeList.begin(); sharedNodeIterator!=sharedNodeList.end(); ++sharedNodeIterator)
		cout << "\nRank of the nodes in serial order is: " << *sharedNodeIterator << endl;
	for (sharedCandLineIterator=sharedCandSerList.begin(); sharedCandLineIterator!=sharedCandSerList.end(); ++sharedCandLineIterator)
		cout << "\nRank of the candidate lines in serial order is: " << *sharedCandLineIterator << endl;
	double LagMultXi[(scenarios+1)*numberOfZones*otherNodeCount]; // Lagrange Multipliers/Dual Variables/Rewards/Penalties corresponding to bus voltage phase angle consensus
	double LagMultPi[numberOfZones*candLineCount]; // Lagrange Multipliers/Dual Variables/Rewards/Penalties corresponding to line building decision integer variable consensus
	double LagCombXi[(scenarios+1)*otherNodeCount]; // Lagrange Multipliers/Dual Variables/Rewards/Penalties corresponding to bus voltage phase angle consensus for MO
	double LagCombPi[candLineCount]; // Lagrange Multipliers/Dual Variables/Rewards/Penalties corresponding to line building decision integer variable consensus MO
	for ( int i = 0; i < (scenarios+1)*numberOfZones*otherNodeCount; ++i ) { // Initialize the Lagrange Multipliers/Dual Variables/Rewards/Penalties to zero for the first iteration
		LagMultXi[i] = 0.0;
	}
	for ( int i = 0; i < numberOfZones*candLineCount; ++i ) { // Initialize the Lagrange Multipliers/Dual Variables/Rewards/Penalties to zero for the first iteration
		LagMultPi[i] = 0.0;
	}
	double lowerBound; // Lower bound of the investment coordination optimization problem
	double upperBound; // Upper bound of the investment coordination optimization problem
	double consensus; // Magnitude of lack of consensus among the different regions, with each other as well as with the MO
	int iterationCounter = 1; // Initialize the iteration counter

	cout << endl << "\n*** DISTRIBUTED STOCHASTIC OPTIMIZATION ALGORITHMIC MARKET MECHANISM DESIGN FIRST STAGE BEGINS ***\n" << endl << endl;
	//do {//%%
	for ( int iterCountOut = 0; iterCountOut < 100; ++iterCountOut ) {
		if (iterCountOut != 0) {
			marketoverInstance->clearDelayedVectors(); // clear the different interim delayed vectors for making them ready for next iteration 
			marketoverInstance->bufferintermediateDecision(iterCountOut); // Buffering the previous iterations' MO decision values, (for comparison basis equal to 0)			
			marketoverInstance->clearVectors(); // clear the different interim vectors for making them ready for next iteration
		}
		else
			marketoverInstance->bufferintermediateDecision(iterCountOut); // Buffering the previous iterations' MO decision values, (for comparison basis equal to 0)			 
		UBIterate = 0;
		LBIterate = 0; 
		cout << "\n*** ITERATION " << iterCountOut+1 << " BEGINS ***\n" << endl << endl;
		for ( int i = 0; i < numberOfZones; ++i ) { // Each region solves its own MILP optimization problem 
			cout << endl << "\n*** MIXED INTEGER LINEAR PROGRAMMING FOR ZONE " << i+1 << " BEGINS ***\n" << endl << endl;
			cout << "\nZonal Calculations of Beliefs about the Investment decision MILP begins" << endl;
			switch (curveChoice) {
				case 1:
					cout << "\nSOLVING MILP" << endl;
					if (solverChoice == 1)
						upperBoundVector[i]=zonalNetVector[i]->MILPAvgHR(*marketoverInstance, LagMultXi, LagMultPi, candLineCount, otherNodeCount); // Perform unit commitment for average heat rate objective
					else
						upperBoundVector[i]=zonalNetVector[i]->MILPAvgHRGUROBI(*marketoverInstance, LagMultXi, LagMultPi, candLineCount, otherNodeCount, environmentGUROBI); // Perform unit commitment for average heat rate objective
					UBIterate += upperBoundVector[i];
					cout << "\nMILP SOLVED" << endl;
					cout << "\nESTIMATING LOWER BOUND" << endl;
					if (solverChoice == 1)
						lowerBoundVector[i]=zonalNetVector[i]->calcMILPBounds(LagMultXi, LagMultPi, candLineCount, otherNodeCount); // Calculate the bounds//%%
					else
						lowerBoundVector[i]=zonalNetVector[i]->calcMILPBoundsGUROBI(LagMultXi, LagMultPi, candLineCount, otherNodeCount, environmentGUROBI); // Calculate the bounds//%%
					LBIterate += lowerBoundVector[i];
					cout << "\nLOWER BOUND ESTIMATED" << endl;
					break;
				case 2:
					zonalNetVector[i]->MILPPiecewiseLin(); // Perform unit commitment for piecewise linear objective
					break;
				case 3:
					zonalNetVector[i]->MILPPolynomial(); // Perform unit commitment for polynomial objective
					break;
				default:
					cout << "\nInvalid choice of Objective function" << endl;
					break;
			}
		}
		for (int k = 0; k <= scenarios; ++k) {
			for ( int i = 0; i < otherNodeCount; ++i ) { // Initialize the Lagrange Multipliers/Dual Variables/Rewards/Penalties to zero for the MO
				double tempXi = 0;
				for ( int j = 0; j < numberOfZones; ++j ) {
					tempXi += LagMultXi[k*numberOfZones*otherNodeCount+j*otherNodeCount+i];
				}
				LagCombXi[k*otherNodeCount+i] = tempXi;
			}
		}
		for ( int i = 0; i < candLineCount; ++i ) { // Initialize the Lagrange Multipliers/Dual Variables/Rewards/Penalties to zero for the MO
			double tempPi = 0;
			for ( int j = 0; j < numberOfZones; ++j ) {
				tempPi += LagMultPi[j*candLineCount+i];
			}
			LagCombPi[i] = tempPi;
		}
		cout << "\nSOLVING MILP FOR MARKET OVERSEER/TRANSMISSION PLANNING COORDINATOR" << endl;
		if (solverChoice == 1)
			marketoverInstance->MILPMarketover(LagCombXi, LagCombPi, candLineCount, otherNodeCount); // The MO solves its own optimization problem and also updates the rewards/penalties//%%
		else
			marketoverInstance->MILPMarketoverGUROBI(LagCombXi, LagCombPi, candLineCount, otherNodeCount, environmentGUROBI); // The MO solves its own optimization problem and also updates the rewards/penalties//%%
		cout << "\nMILP SOLVED" << endl;
		cout << "\nESTIMATING UPPER BOUND FOR MARKET OVERSEER/TRANSMISSION PLANNING COORDINATOR" << endl;
		upperBound = marketoverInstance->getGlobalUpper(LagCombXi, LagCombPi, upperBoundVector, numberOfZones); // MO calculates the global upper bound after every iteration//%%
		cout << "\nESTIMATED UPPER BOUND FOR MARKET OVERSEER/TRANSMISSION PLANNING COORDINATOR" << endl;
		cout << "\nESTIMATING LOWER BOUND FOR MARKET OVERSEER/TRANSMISSION PLANNING COORDINATOR" << endl;
		if (solverChoice == 1)
			lowerBound = marketoverInstance->LBMarketover(LagCombXi, LagCombPi, candLineCount, otherNodeCount); // MO calculates the lower bound for itself after every iteration//%%
		else
			lowerBound = marketoverInstance->LBMarketoverGUROBI(LagCombXi, LagCombPi, candLineCount, otherNodeCount, environmentGUROBI); // MO calculates the lower bound for itself after every iteration//%%
		LBIterate += lowerBound;
		//lowerBound += marketoverInstance->getGlobalLower(lowerBoundVector, numberOfZones); // MO calculates the global lower bound after every iteration//%%
		cout << "\nESTIMATED LOWER BOUND FOR MARKET OVERSEER/TRANSMISSION PLANNING COORDINATOR" << endl;

		consensus = marketoverInstance->getGlobalConsensus(); // MO calculates the global consensus after every iteration//%%
		UBItVec.push_back(upperBound);
		LBItVec.push_back(LBIterate);
		ratioIterate.push_back(upperBound/LBIterate);
		globalCons.push_back(consensus);
		cout << "\nUPDATING THE VALUES OF THE REWARDS/PENALTIES/LAGRANGE MULTIPLIERS/DUAL VARIABLES" << endl;
		for ( int i = 0; i < 10; ++i ) { // Initialize the Lagrange Multipliers/Dual Variables/Rewards/Penalties to zero for the first iteration
			double tempLagMultXi = LagMultXi[i];
			double tempLagMultPi = LagMultPi[i];
			LagMultXi[i] = marketoverInstance->rewardPenaltyCont(tempLagMultXi, i, iterCountOut);
			LagMultPi[i] = marketoverInstance->rewardPenaltyInteger(tempLagMultPi, i, iterCountOut);
		}
		//marketoverInstance->rewardPenaltyUpdate(LagMultXi, LagMultPi, i, iterCountOut);
		cout << "\nUPDATED THE VALUES OF THE REWARDS/PENALTIES/LAGRANGE MULTIPLIERS/DUAL VARIABLES" << endl;
		//++iterationCounter; // Increment the iteration counter before next iteration
	}
	//} while ((((1-abs(lowerBound/upperBound))<=0.05) || ((1-abs(lowerBound/upperBound))>=-0.05)));
	//&& (consensus<=0.05)); // while the tolerance is reached//%%
	cout << "\nZonal Calculations of Beliefs about the Investment decision MILP ends" << endl;
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
	cout << endl << "\n*** AUXILIARY PROBLEM PRINCIPLE (APP) OPTIMIZATION ALGORITHMIC MARKET MECHANISM DESIGN SECOND STAGE BEGINS ***\n" << endl << endl;
	vector<double> interAngleMessage[numberOfZones]; // Array of the zonal vectors of the intermediate values of Lagrange Multipliers
	vector<int> zonalGlobRank[numberOfZones]; // Array of the global ranks of the shared nodes for each zone
	double ObjIterate[numberOfZones]; // Array of zonal optimum objectives
	//do {//%%
	for ( int iterCountOut = 0; iterCountOut < 1; ++iterCountOut ) {
		cout << "\n*** ITERATION " << iterCountOut+1 << " BEGINS ***\n" << endl << endl;
		for ( int i = 0; i < numberOfZones; ++i ) { // Each region solves its own MILP optimization problem 
			cout << endl << "\n*** APP QUADRATIC PROGRAMMING FOR ZONE " << i+1 << " BEGINS ***\n" << endl << endl;
			cout << "\nZonal Calculations of Beliefs about the generation and flow decision APP-QP begins" << endl;
			switch (curveChoice) {
				case 1:
					cout << "\nSOLVING THE OPTIMIZATION SUB-PROBLEM" << endl;
					ObjIterate[i]=zonalNetVector[i]->APPQPAvgHR(*marketoverInstance, APPLagMultipliers, candLineCount, otherNodeCount, environmentGUROBI, iterCountOut); // Perform unit commitment for average heat rate objective
					interAngleMessage[i] = zonalNetVector[i]->getZonalDecision();
					zonalGlobRank[i] = zonalNetVector[i]->getZonalRanks();
					cout << "\nOPTIMIZATION SUB-PROBLEM SOLVED" << endl;
					break;
				case 2:
					zonalNetVector[i]->MILPPiecewiseLin(); // Perform unit commitment for piecewise linear objective
					break;
				case 3:
					zonalNetVector[i]->MILPPolynomial(); // Perform unit commitment for polynomial objective
					break;
				default:
					cout << "\nInvalid choice of Objective function" << endl;
					break;
			}
		}
		for ( int i = 0; i <= otherNodeCount; ++i ) { // Initialize the Lagrange Multipliers/Dual Variables/Rewards/Penalties to zero for the MO
			double tempXi = 0;
			for ( int j = 0; j < numberOfZones; ++j ) {
				tempXi += LagMultXi[j*otherNodeCount+i];
			}
			LagCombXi[i] = tempXi;
		}
		for ( int i = 0; i <= candLineCount; ++i ) { // Initialize the Lagrange Multipliers/Dual Variables/Rewards/Penalties to zero for the MO
			double tempPi = 0;
			for ( int j = 0; j < numberOfZones; ++j ) {
				tempPi += LagMultPi[j*candLineCount+i];
			}
			LagCombPi[i] = tempPi;
		}

		consensus = marketoverInstance->getGlobalConsensus(); // MO calculates the global consensus after every iteration//%%
		UBItVec.push_back(upperBound);
		LBItVec.push_back(LBIterate);
		ratioIterate.push_back(upperBound/LBIterate);
		globalCons.push_back(consensus);
		cout << "\nUPDATING THE VALUES OF THE REWARDS/PENALTIES/LAGRANGE MULTIPLIERS/DUAL VARIABLES" << endl;
		for ( int i = 0; i < 1000; ++i ) { // Initialize the Lagrange Multipliers/Dual Variables/Rewards/Penalties to zero for the first iteration
			double tempLagMultXi = LagMultXi[i];
			double tempLagMultPi = LagMultPi[i];
			LagMultXi[i] = marketoverInstance->rewardPenaltyCont(tempLagMultXi, i);
			LagMultPi[i] = marketoverInstance->rewardPenaltyInteger(tempLagMultPi, i);
		}
		cout << "\nUPDATED THE VALUES OF THE REWARDS/PENALTIES/LAGRANGE MULTIPLIERS/DUAL VARIABLES" << endl;
		//++iterationCounter; // Increment the iteration counter before next iteration
		marketoverInstance->clearVectors(); // clear the different interim vectors for making them ready for next iteration 
	}
	//} while ((((1-abs(lowerBound/upperBound))<=0.05) || ((1-abs(lowerBound/upperBound))>=-0.05)));
	//&& (consensus<=0.05)); // while the tolerance is reached//%%
	cout << "\nZonal Calculations of Beliefs about the Investment decision MILP ends" << endl;
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
	delete marketoverInstance; // Free the memory of the Marketoverseer class object//%%
	for ( int i = 0; i < numberOfZones; ++i ) {
		delete zonalNetVector[i]; // Free the memory of the Nettran class objects
	}
*/
	delete environmentGUROBI; // Free the memory of the GUROBI environment object
	return 0; // Indicates successful Program Completion
} // End of main method


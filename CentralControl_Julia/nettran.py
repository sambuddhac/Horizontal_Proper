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
			matrixNetFile >> nodeNumber >> sharedELines >> sharedCLines >> genNumber >> loadNumber >> tranNumber >> internalCLines; // get the dimensions of the Network
			for ( int l = 0; l < nodeNumber; ++l ) {
				//cout << "\nCreating the " << l + 1 << " -th Node:\n";
		
				Node *nodeInstance = new Node( univNodeNum + l + 1, l + 1, zonalIndex ); // creates nodeInstance object with ID l + 1

				nodeObject.push_back( nodeInstance ); // pushes the nodeInstance object into the vector

			} // end initialization for Nodes
			matrixNetFile.close(); // close the network file

			#/* Generators */

			#/* Instantiate Generators */
			// Open the .txt file to read the Powergenerator parameter values
			ifstream matrixGenFile( genFile, ios::in ); // ifstream constructor opens the file of Generators

			// exit program if ifstream could not open file
			if ( !matrixGenFile ) {
				cerr << "\nFile for Generators could not be opened\n" << endl;
				exit( 1 );
			} // end if
			int genFields; // Number of columns in the generator database
			matrixGenFile >> genFields; // get the dimensions of the Generator matrix
			double matrixGen[ genNumber ][ genFields ]; // Generator matrix
			for ( int i = 0; i < genNumber; ++i ) {
				for ( int j = 0; j < genFields; ++j ) {
					matrixGenFile >> matrixGen[ i ][ j ]; // read the Generator matrix
				}
			}
			for (int j = 0; j < genNumber; ++j) {
				int gNodeID; // node object ID to which the particular generator object is connected
				do {
					gNodeID = matrixGen[ j ][ 0 ];
				} while ( ( gNodeID <= 0 ) || ( gNodeID > nodeNumber ) ); // validity check

				double c2, c1, c0, PgMax, PgMin, tanTheta, minCost; // Parameters for Generator
				do {
					//Quadratic Coefficient: 
					c2 = matrixGen[ j ][ 1 ] * (pow(100, 2.0));
					//Linear coefficient: 
					c1 = matrixGen[ j ][ 2 ] * 100;
					//Constant term: 
					c0 = matrixGen[ j ][ 3 ];
					//Maximum Limit: 
					PgMax = matrixGen[ j ][ 4 ] / 100;
					//Minimum Limit: 
					PgMin = matrixGen[ j ][ 5 ] / 100;
					/* Secant Approximation of the Quadratic Cost curve */
					//Tangent Ratio of the secant approximation of the intercepted cost curve
					tanTheta = (c2*(pow(PgMax, 2.0))+c1*PgMax-c2*(pow(PgMin, 2.0))-c1*PgMin)/(PgMax-PgMin);
					minCost = c2*(pow(PgMin, 2.0))+c1*PgMin+c0-tanTheta*PgMin;
					//Intercept value or cost at minimum power level
				
				} while ( (c2 < 0 ) || ( c1 < 0 ) || ( PgMax <= 0 ) || ( PgMin < 0 ) || ( PgMax <= PgMin ) ); 
				// check the bounds and validity of the parameter values
				Powergenerator *genInstance = new Powergenerator( univGenNum + j+1, nodeObject[ univNodeNum + gNodeID - 1 ],  tanTheta, minCost, PgMax, PgMin );
				genObject.push_back(genInstance); // push the generator object into the array
			}
			matrixGenFile.close(); // Close the generator file
 			if (tranNumber > 0) {
				#/* Transmission Lines */
				ifstream matrixTranFile( tranFile, ios::in ); // ifstream constructor opens the file of Transmission lines

				// exit program if ifstream could not open file
				if ( !matrixTranFile ) {
					cerr << "\nFile for Transmission lines could not be opened\n" << endl;
					exit( 1 );
				} // end if
				int translFields; // Number of columns in the transmission lines file
				matrixTranFile >> translFields; // get the dimensions of the Transmission line matrix
				double matrixTran[ tranNumber ][ translFields ]; // Transmission line matrix
				for ( int i = 0; i < tranNumber; ++i ) {
					for ( int j = 0; j < translFields; ++j ) {
						matrixTranFile >> matrixTran[ i ][ j ]; // read the Transmission line matrix
					}
				}

				#/* Instantiate Transmission Lines */
				for ( int k = 0; k < tranNumber; ++k ) {
					int tNodeID1, tNodeID2; // node object IDs to which the particular transmission line object is connected
					do {
						//node IDs of the node objects to which this transmission line is connected.
						tNodeID1 = matrixTran[ k ][ 0 ]; //From end
						tNodeID2 = matrixTran[ k ][ 1 ]; //To end
					} while ( ( tNodeID1 <= 0 ) || ( tNodeID1 > nodeNumber ) || ( tNodeID2 <= 0 ) || ( tNodeID2 > nodeNumber ) || ( tNodeID1 == tNodeID2) ); // validity check
					double reacT, ptMax; // Parameters for Transmission Line
					do {
						//Reactance:
						reacT = matrixTran[ k ][ 2 ];
						//values of maximum allowable power flow on line in the forward and reverse direction:
						ptMax = matrixTran[ k ][ 3 ]/100;
					} while ( reacT <= 0 ); // check the bounds and validity of the parameter values
					// creates transLineInstance object with ID k + 1
					transmissionLine *transLineInstance = new transmissionLine( univTranNum + k + 1, nodeObject[ univNodeNum + tNodeID1 - 1 ], nodeObject[ univNodeNum + tNodeID2 - 1 ], ptMax, reacT ); 
					translObject.push_back( transLineInstance ); // pushes the transLineInstance object into the vector

				} // end initialization for Transmission Lines 
				matrixTranFile.close(); // Close the transmission line file 
			}
			if (internalCLines>0) {
				#/* Internal Candidate Transmission Lines */
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
				#/* Instantiate Internal Candidate Transmission Lines */
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
						costPerCap = matrixIntCETran[ k ][ 6 ]*ptMax*100; // capital cost for the construction 
						presAbsence = matrixIntCETran[ k ][ 7 ]; // status of the construction 
					} while ( ( reacT <= 0 ) ); // check the bounds and validity of the parameter values
			
					// creates intCandLineInstance object with ID k + 1
					intCandLine *intCandLineInstance = new intCandLine( univIntCandNum + k + 1, nodeObject[ univNodeNum + tNodeID1 - 1 ], nodeObject[ univNodeNum + tNodeID2 - 1 ], ptMax, reacT, interestRate, lifeTime, costPerCap, presAbsence );// Create the internal candidate transmission line object with node 1 
					intCandLineObject.push_back( intCandLineInstance ); // pushes the transLineInstance object into the vector
				} // end initialization for candidate Transmission Lines
				matrixIntCETranFile.close(); // Close the candidate lines file
			}
			#/* Loads */
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
			for ( int l = univNodeNum; l < univNodeNum+nodeNumber; ++l ) {
				(nodeObject[l])->initLoad( countOfScenarios ); // Initialize the default loads on all nodes to zero

			} // end initialization for Nodes		
			#/* Instantiate Loads */
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

				Load *loadInstance = new Load( univLoadNum + j + 1, nodeObject[ univNodeNum + lNodeID - 1 ], loadFields-1, P_Load ); // creates loadInstance object object with ID number j + 1

				loadObject.push_back( loadInstance ); // pushes the loadInstance object into the vector

			} // end initialization for Loads
			matrixLoadFile.close(); // Closes the load file
		univNodeNum += nodeNumber #Increment the universal node number
		nodeNumVector.append(nodeNumber)
		univGenNum += genNumber #Increment the universal generator number
		univTranNum += tranNumber #Increment the universal transmission line number
		univLoadNum += loadNumber #Increment the universal load number
		univIntCandNum += internalCLines #Increment the universal intra zonal candidate line number
	}
	for ( int zonalIndex = 1; zonalIndex <= zonalCount; ++zonalIndex ) { // Iterate through the zones 
		strcpy( netFile, initSummary[(zonalIndex-1)].c_str() ); // String for storing the name of the network file
		strcpy( sharedLineFile, initSummary[(zonalIndex-1)+2*zoneCount].c_str() ); // String for storing the name of the shared existing lines file
		strcpy( candLineFile, initSummary[(zonalIndex-1)+5*zoneCount].c_str() ); // String for storing the name of the candidate lines file
		do {
			ifstream matrixNetFile( netFile, ios::in ); // ifstream constructor opens the file of Network	
			// exit program if ifstream could not open file
			if ( !matrixNetFile ) {
				cerr << "\nFile for Network could not be opened\n" << endl;
				exit( 1 );
			} // end if
			matrixNetFile >> nodeNumber >> sharedELines >> sharedCLines >> genNumber >> loadNumber >> tranNumber >> internalCLines; // get the dimensions of the Network
			/* Shared Existing Transmission Lines */
			ifstream matrixSETranFile( sharedLineFile, ios::in ); // ifstream constructor opens the file of Transmission lines

			// exit program if ifstream could not open file
			if ( !matrixSETranFile ) {
				cerr << "\nFile for Shared Existing Transmission lines could not be opened\n" << endl;
				exit( 1 );
			} // end if
			int tranSEFields; // Number of columns in the transmission lines file
			matrixSETranFile >> tranSEFields; // get the dimensions of the Transmission line matrix
			double matrixSETran[ sharedELines ][ tranSEFields ]; // Transmission line matrix
			for ( int i = 0; i < sharedELines; ++i ) {
				for ( int j = 0; j < tranSEFields; ++j ) {
					matrixSETranFile >> matrixSETran[ i ][ j ]; // read the Transmission line matrix
				}
			}

			/* Instantiate Shared Existing Transmission Lines */
			for ( int k = 0; k < sharedELines; ++k ) {
				int serNum, tNodeID1, tNodeID2, nodeZone1, nodeZone2; // node object IDs to which the particular transmission line object is connected
				//cout << "Tran File Test Message 1 from line " << k << " before creation1" << endl;  
				do {
					//node IDs of the node objects to which this transmission line is connected.
					serNum = matrixSETran[ k ][ 0 ]; // global serial number of the shared existing transmission line
					tNodeID1 = matrixSETran[ k ][ 1 ]; // From end node 
					nodeZone1 = matrixSETran[ k ][ 2 ]; // From end zone number
					tNodeID2 = matrixSETran[ k ][ 3 ]; // To end node
					nodeZone2 = matrixSETran[ k ][ 4 ]; // To end zone number 
				} while ( (( nodeZone1 != zonalIndex ) && ( nodeZone2 != zonalIndex )) || (( nodeZone1 == zonalIndex ) && ( nodeZone2 == zonalIndex )) ); // validity check
				double reacT, ptMax; // Parameters for Transmission Line
				do {
					//Reactance:
					reacT = matrixSETran[ k ][ 5 ];
					cout << "Reactance read " << reacT;
					//values of maximum allowable power flow on line in the forward and reverse direction:
					ptMax = matrixSETran[ k ][ 6 ]/100;
				} while ( reacT <= 0 ); // check the bounds and validity of the parameter values
				if(find(SELineSerList.begin(), SELineSerList.end(), serNum)==SELineSerList.end()) { // To avoid duplication, check if the SE Line has already been created, if not
					// creates SELine object with ID k + 1
					// Renumbering the node IDs with the aggregated count
					int newNodeBase = 0; // 
					for (int aggrCount = 0; aggrCount < nodeZone1; ++aggrCount) {
						newNodeBase+=nodeNumVector[aggrCount];
					}
					tNodeID1 += newNodeBase;
					newNodeBase = 0; // 
					for (int aggrCount = 0; aggrCount < nodeZone2; ++aggrCount) {
						newNodeBase+=nodeNumVector[aggrCount];
					}
					tNodeID2 += newNodeBase;
					SELine *SELineInstance = new SELine( k + 1, serNum, nodeObject[ tNodeID1 - 1 ], nodeObject[ tNodeID2 - 1 ], ptMax, reacT ); // Create the shared existing transmission line object with node 1 
					SELineObject.push_back( SELineInstance ); // pushes the transLineInstance object into the vector
					SELineSerList.push_back(serNum);
					++univSELineNum; // Increment to account for the total number of shared existing lines
				}		
			} // end initialization for Shared Existing Transmission Lines
			matrixSETranFile.close(); // Close the shared existing file

			/* Shared Candidate Transmission Lines */
			ifstream matrixCETranFile( candLineFile, ios::in ); // ifstream constructor opens the file of candidate Transmission lines

			// exit program if ifstream could not open file
			if ( !matrixCETranFile ) {
				cerr << "\nFile for Shared Candidate Transmission lines could not be opened\n" << endl;
				exit( 1 );
			} // end if
			int tranCFields; // Number of columns in the candidate transmission lines file
			matrixCETranFile >> tranCFields; // get the dimensions of the Transmission line matrix
			double matrixCETran[ sharedCLines ][ tranCFields ]; // candidate Transmission line matrix
			for ( int i = 0; i < sharedCLines; ++i ) {
				for ( int j = 0; j < tranCFields; ++j ) {
					matrixCETranFile >> matrixCETran[ i ][ j ]; // read the candidate Transmission line matrix
				}
			}

			/* Instantiate Shared Candidate Transmission Lines */
			for ( int k = 0; k < sharedCLines; ++k ) {
				int serNum, tNodeID1, tNodeID2, nodeZone1, nodeZone2, presAbsence, lifeTime; // node object IDs to which the particular transmission line object is connected
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
					costPerCap = matrixCETran[ k ][ 9 ]*ptMax*100; // capital cost for the construction 
					presAbsence = matrixCETran[ k ][ 10 ]; // status of the construction 
				} while ( ( reacT <= 0 ) ); // check the bounds and validity of the parameter values
				if(find(candLineSerList.begin(), candLineSerList.end(), serNum)==candLineSerList.end()) { // To avoid duplication, check if the shared candidate Line has already been created, if not
					int newNodeBase = 0; // 
					for (int aggrCount = 0; aggrCount < nodeZone1; ++aggrCount) {
						newNodeBase+=nodeNumVector[aggrCount];
					}
					tNodeID1 += newNodeBase;
					newNodeBase = 0; // 
					for (int aggrCount = 0; aggrCount < nodeZone2; ++aggrCount) {
						newNodeBase+=nodeNumVector[aggrCount];
					}
					tNodeID2 += newNodeBase;			
					// creates candLineInstance object with ID k + 1
					candLine *candLineInstance = new candLine( k + 1, serNum, nodeObject[ tNodeID1 - 1 ], nodeObject[ tNodeID2 - 1 ], ptMax, reacT, interestRate, lifeTime, costPerCap, presAbsence );// Create the shared candidate transmission line object with node 1 
					candLineObject.push_back( candLineInstance ); // pushes the transLineInstance object into the vector
					candLineSerList.push_back(serNum);
					++univCandLineNum; // Increment to account for the total number of shared candidate lines
				}
			} // end initialization for candidate Transmission Lines
			matrixCETranFile.close(); // Close the candidate lines file
		} while ( (sharedELines <= 0 ) || ( sharedCLines <= 0 ) );
		// check the bounds and validity of the parameter values
	}
	assignProb();
} // end constructor

Nettran::~Nettran() // destructor
{
	cout << "\nSimulation ended" << endl;
} // destructor ends

void Nettran::MILPPiecewiseLin(void) // Function MILPPiecewiseLin() implements the Mixed Integer Linear Programming Unit Commitment Solver routine by calling GLPK routines for piecewise linear objective
{
	cout << "\nUnder Construction" << endl;
}

void Nettran::assignProb() {
	for (int f=0; f < countOfScenarios; ++f) {
		probability.push_back((static_cast<double>(1)/(countOfScenarios)));
	}
}

void Nettran::MILPPolynomial() // Function MILPPolynomial() implements the Mixed Integer Linear Programming Unit Commitment Solver routine by calling GLPK routines for polynomial convex objective
{
	cout << "\nUnder Construction" << endl;
}

double Nettran::MILPAvgHRGUROBI(GRBEnv* environmentGUROBI) // Function MILPAvgHRGUROBI() implements the Mixed Integer Linear Programming Unit Commitment Solver routine by calling GUROBI routines for average heat rate objective for Horizontal Coordination Investment decision making
{
	// CREATION OF THE MIP SOLVER INSTANCE //
	clock_t begin = clock(); // start the timer
	vector<Powergenerator*>::iterator genIterator; // Iterator for Powergenerator objects
	vector<transmissionLine*>::iterator tranIterator; // Iterator for Transmission line objects
	vector<Load*>::iterator loadIterator; // Iterator for load objects
	vector<Node*>::iterator nodeIterator; // Iterator for node objects
	vector<candLine*>::iterator candIterator; // Iterator for candidate lines
	vector<SELine*>::iterator exsharedIterator; // Iterator for shared existing lines
	vector<intCandLine*>::iterator intCandIterator; // Iterator for candidate lines
	vector<double>::iterator probIterator; 
	for (probIterator = probability.begin(); probIterator != probability.end(); ++probIterator) {
		cout << "The value of probability is " << *probIterator << endl; 
	}

	string outSummaryFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/outputCent/outputSummaryResults/OutSummaryGUROBI.txt";
	ofstream outPutFile(outSummaryFileName, ios::out); // Create Output File to output the Summary of Results
	if (!outPutFile){
		cerr << "\nCouldn't open the file" << endl;
		exit(1);
	}

        int dimRow = countOfScenarios*(2 * univGenNum + 4 * univCandLineNum + 2 * univSELineNum + 2 * univTranNum + univNodeNum + 4*univIntCandNum); // Total number of rows of the A matrix (number of structural constraints of the LP) first term to account for lower and upper generating limits, second term for lower and upper line limits & lower and upper definition limits of candidate shared lines, third term for lower and upper line limits for shared existing lines, fourth term for lower and upper line limits for internal zonal lines, the fifth term to account for nodal power balance constraints, and sixth term to account for the internal candidate lines
        int dimCol = countOfScenarios*(univGenNum+univNodeNum+univCandLineNum+univIntCandNum)+univCandLineNum+univIntCandNum; // Total number of columns of the LP (number of Decision Variables) first term to account for power generation MW outputs, second term for voltage phase angles for internal zonal nodes, third term for power flow values and binary integer decision variable values for shared candidate lines, fourth term for the voltage phase angles of other-zone nodes connected through shared existing and candidate lines, and fifth term for the decision variables for internal candidate lines
	outPutFile << "\nTotal Number of Structural Constraints (Rows) is: " << dimRow << endl;
	outPutFile << "\nTotal Number of Decision Variables (Columns) is: " << dimCol << endl;
	// Instantiate GUROBI Problem model
	GRBModel *modelSubnetMILP = new GRBModel(*environmentGUROBI);
	cout << "\nGurobi model created" << endl;
    	modelSubnetMILP->set(GRB_StringAttr_ModelName, "assignment");
	cout << "\nGurobi model created and name set" << endl;
	GRBVar decvar[dimCol+1];
	cout << "\nGurobi decision variables created" << endl;
	double z; // variable to store the objective value

	// SPECIFICATION OF PROBLEM PARAMETERS //
	// Dummy Decision Variable //
	cout << "\nGurobi decision variables to be assigned" << endl;
	decvar[0] = modelSubnetMILP->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
	//Decision Variable Definitions, Bounds, and Objective Function Co-efficients//
	cout << "\nGurobi dummy decision variable created" << endl;
	int colCount = 1;
	//Columns corresponding to Power Generation continuous variables for different generators//
	outPutFile << "\nCoefficients of Power generator variables" << endl;
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
			decvar[colCount] = modelSubnetMILP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
			outPutFile << colCount << "\t";
			outPutFile << ((*genIterator)->getLinCoeff())/100 << " $/MW" << endl;
			++colCount;
		}
	}
	outPutFile << "\nTotal number of columns after accounting for Power Generation continuous variables for different generators: " << colCount << endl;

	//Columns corresponding to Voltage Phase Angles continuous variables for different nodes//
	outPutFile << "\nCoefficients of Voltage Phase Angles continuous variables for different nodes" << endl;
	outPutFile << "\nVariable Count\tShared Node\tGlobal Rank" << endl;
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {	
		for (nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); ++nodeIterator){
			if ( (((*nodeIterator)->getSharedFlag()) == 1) || (((*nodeIterator)->getCandFlag()) == 1) ) {
				decvar[colCount] = modelSubnetMILP->addVar((0), (44/7), 0.0, GRB_CONTINUOUS);
				outPutFile << colCount << "\tYes\t" << ((*nodeIterator)->getNodeID()) << endl;	
			}
			else {
				decvar[colCount] = modelSubnetMILP->addVar((0), (44/7), 0.0, GRB_CONTINUOUS);
				outPutFile << colCount << "\tNo\t" << ((*nodeIterator)->getNodeID()) << endl;	
			}
			++colCount;
		}
	}
	outPutFile << "\nTotal number of columns after accounting for Voltage Phase Angles continuous variables for different nodes: " << colCount << endl;

	//Columns corresponding to Shared Candidate Line Flows continuous variables//
	outPutFile << "\nCoefficients corresponding to Shared Candidate Line Flows continuous variables" << endl;
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			decvar[colCount] = modelSubnetMILP->addVar((-GRB_INFINITY), (GRB_INFINITY), 0.0, GRB_CONTINUOUS);
			outPutFile << colCount << "\t";
			outPutFile << 0 << endl;
			++colCount;
		}
	}
	outPutFile << "\nTotal number of columns after accounting for Shared Candidate Line Flows continuous variables: " << colCount << endl;

	//Columns corresponding to Shared Candidate Line Construction Decision Binary Integer variables//
	outPutFile << "\nCoefficients corresponding to Shared Candidate Line Construction Decision Binary Integer variables" << endl;
	outPutFile << "\nVariable Count\tGlobal Rank\tInvestment Cost" << endl;
	for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
		decvar[colCount] = modelSubnetMILP->addVar(0, 1, 0.0, GRB_BINARY);
		outPutFile << colCount << "\t" << ((*candIterator)->getTranslID()) << "\t" << ((*candIterator)->getInvestCost()) << endl;	
		++colCount;
	}
	outPutFile << "\nTotal number of columns after accounting for Shared Candidate Line Construction Decision Binary Integer variables: " << colCount << endl;

	//Columns corresponding to Internal Candidate Line Flows continuous variables//
	outPutFile << "\nCoefficients corresponding to Internal Candidate Line Flows continuous variables" << endl;
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			decvar[colCount] = modelSubnetMILP->addVar((-GRB_INFINITY), (GRB_INFINITY), 0.0, GRB_CONTINUOUS);
			outPutFile << colCount << "\t";
			outPutFile << 0 << endl;
			++colCount;
		}
	}
	outPutFile << "\nTotal number of columns after accounting for Internal Candidate Line Flows continuous variables: " << colCount << endl;

	//Columns corresponding to Internal Candidate Line Construction Decision Binary Integer variables//
	outPutFile << "\nCoefficients corresponding to Internal Candidate Line Construction Decision Binary Integer variables" << endl;
	outPutFile << "\nVariable Count\tInvestment Cost" << endl;
	for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
		decvar[colCount] = modelSubnetMILP->addVar(0, 1, 0.0, GRB_BINARY);
		outPutFile << colCount << "\t" << ((*intCandIterator)->getInvestCost()) << endl;	
		++colCount;
	}
	outPutFile << "\nTotal number of columns after accounting for Internal Candidate Line Construction Decision Binary Integer variables: " << colCount << endl;
	outPutFile << "\nTotal Number of columns for generation, angles, integer variables, and flows: " << colCount-1 << endl;
	outPutFile << "\nDecision Variables and Objective Function defined" << endl;
	outPutFile << "\nTotal Number of columns: " << colCount-1 << endl;
	//Setting Objective//
	for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
		cout << "Gen Max cost : " << ((*genIterator)->getLinCoeff())*((*genIterator)->getPMax()) << endl;
	}
	for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
		cout << "Shared Transmission investment cost : " << (2*((*candIterator)->getInvestCost())) << endl;	
	}
	for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
		cout << "Internal Transmission investment cost : " << ((*intCandIterator)->getInvestCost()) << endl;	
	}
	GRBLinExpr obj = 0.0;
	// Objective Contribution from Dummy Decision Variable //
	obj += 0*(decvar[0]);
	colCount = 1;
	//Columns corresponding to Power Generation continuous variables for different generators//
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
			obj += (probability.at(scenCounter))*(((*genIterator)->getLinCoeff()))*(decvar[colCount]);
			++colCount;
		}
	}
	//Columns corresponding to Voltage Phase Angles continuous variables for different nodes//
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {	
		for (nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); ++nodeIterator){
			obj += 0*(decvar[colCount]);
			++colCount;
		}
	}
	//Columns corresponding to Shared Candidate Line Flows continuous variables//
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			obj += 0*(decvar[colCount]);
			++colCount;
		}
	}
	//Columns corresponding to Shared Candidate Line Construction Decision Binary Integer variables//
	for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
		obj += (2*((*candIterator)->getInvestCost()))*(decvar[colCount]);	
		++colCount;
	}

	//Columns corresponding to Internal Candidate Line Flows continuous variables//
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			obj += 0*(decvar[colCount]);
			++colCount;
		}
	}
	//Columns corresponding to Internal Candidate Line Construction Decision Binary Integer variables//
	for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
		obj += ((*intCandIterator)->getInvestCost())*(decvar[colCount]);	
		++colCount;
	}

	modelSubnetMILP->setObjective(obj, GRB_MINIMIZE);
	cout << " Objective Function and Decision Variables have been defined and the colCount is " << colCount-1 << endl;
	//Row Definitions: Specification of b<=Ax<=b//
	GRBLinExpr lhs[dimRow+1];
	//Row Definitions and Bounds Corresponding to Constraints/
	// Constraints corresponding to supply-demand balance
	string outPGenFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/outputCent/outputPowerResults/OutPowerGenGUROBI.txt"; 
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
	modelSubnetMILP->addConstr(lhs[0], GRB_EQUAL, 0);
	int rCount = 1; // Initialize the row count
	vector<int> busCount; // vector for storing the node/bus serial
	outPutFile << "Constraints corresponding to Supply-Demand Balance right hand side" << endl;
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); ++nodeIterator){
			outPutFile << "\nGeneration\t" << rCount << "\n";
			int genListLength = (*nodeIterator)->getGenLength(); // get the number
			lhs[rCount]=0;
			//cout << "\nScenario Count: " << scenCounter << " node count: " << (*nodeIterator)->getNodeID() << " Conn Gen: " << genListLength << endl;
			for (int cCount = 1; cCount <= genListLength; ++cCount){
				//cout << "\nSerial: " << cCount << " Generator Serial: " << (*nodeIterator)->getGenSer(cCount) << endl;
				lhs[rCount] += 1*(decvar[scenCounter*univGenNum+(*nodeIterator)->getGenSer(cCount)]);
				outPutFile << "\n" << rCount << "\t" << scenCounter*univGenNum+(*nodeIterator)->getGenSer(cCount) << "\t" << 1.0 << endl;
			}
			outPutFile << "\nIntrazonal Node Angles\t" << rCount << "\n";
			lhs[rCount] += (((*nodeIterator)->getToReact())-((*nodeIterator)->getFromReact()))*(decvar[countOfScenarios*univGenNum+rCount]);
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*univGenNum+rCount << "\t" << (((*nodeIterator)->getToReact())-((*nodeIterator)->getFromReact())) << "\t" << ((*nodeIterator)->getFromReact()) << "\t" << ((*nodeIterator)->getToReact()) << endl;
			outPutFile << "\nConnected Intrazonal Node Angles\t" << rCount << "\n";
			int connNodeListLength = (*nodeIterator)->getConNodeLength(); // get the number of nodes connected to this node
			for (int cCount = 1; cCount <= connNodeListLength; ++cCount){
				if (((*nodeIterator)->getConnReact(cCount))<=0)
					lhs[rCount] -= (((*nodeIterator)->getConnReact(cCount)))*(decvar[countOfScenarios*univGenNum+scenCounter*univNodeNum+((*nodeIterator)->getConnSer(cCount))]);
				else
					lhs[rCount] += (((*nodeIterator)->getConnReact(cCount)))*(decvar[countOfScenarios*univGenNum+scenCounter*univNodeNum+((*nodeIterator)->getConnSer(cCount))]);		
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*univGenNum+scenCounter*univNodeNum+((*nodeIterator)->getConnSer(cCount)) << "\t" <<  (((*nodeIterator)->getConnReact(cCount))) << "\n";

			}
			outPutFile << "\nConnected Candidate Lines for which this is the From node\t" << rCount << "\n";
			int connCandListLengthF = (*nodeIterator)->getCandLineLengthF(); // get the number of candidate lines connected to this from node 
			for (int cCount = 1; cCount <= connCandListLengthF; ++cCount){
				lhs[rCount] += (-1)*(decvar[countOfScenarios*(univGenNum+univNodeNum)+scenCounter*univCandLineNum+(*nodeIterator)->getCandSerF(cCount)]);
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*(univGenNum+univNodeNum)+scenCounter*univCandLineNum+(*nodeIterator)->getCandSerF(cCount) << "\t" << -1.0 << "\n";
			}
			outPutFile << "\nConnected Candidate Lines for which this is the To node\t" << rCount << "\n";
			int connCandListLengthT = (*nodeIterator)->getCandLineLengthT(); // get the number of candidate lines connected to this to node 
			for (int cCount = 1; cCount <= connCandListLengthT; ++cCount){
				lhs[rCount] += decvar[countOfScenarios*(univGenNum+univNodeNum)+scenCounter*univCandLineNum+(*nodeIterator)->getCandSerT(cCount)];
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*(univGenNum+univNodeNum)+scenCounter*univCandLineNum+(*nodeIterator)->getCandSerT(cCount) << "\t" << 1.0 << "\n";
			}
			outPutFile << "\nConnected Internal Candidate Lines for which this is the From node\t" << rCount << "\n";
			int connintCandListLengthF = (*nodeIterator)->getIntCandLineLengthF(); // get the number of internal candidate lines connected to this from node 
			for (int cCount = 1; cCount <= connintCandListLengthF; ++cCount){
				lhs[rCount] += (-1)*(decvar[countOfScenarios*(univGenNum+univNodeNum+univCandLineNum)+univCandLineNum+scenCounter*univIntCandNum+(*nodeIterator)->getIntCandSerF(cCount)]);
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*(univGenNum+univNodeNum+univCandLineNum)+univCandLineNum+scenCounter*univIntCandNum+(*nodeIterator)->getIntCandSerF(cCount) << "\t" << -1.0 << "\n";
			}
			outPutFile << "\nConnected Internal Candidate Lines for which this is the To node\t" << rCount << "\n";
			int connintCandListLengthT = (*nodeIterator)->getIntCandLineLengthT(); // get the number of internal candidate lines connected to this to node 
			for (int cCount = 1; cCount <= connintCandListLengthT; ++cCount){
				lhs[rCount] += decvar[countOfScenarios*(univGenNum+univNodeNum+univCandLineNum)+univCandLineNum+scenCounter*univIntCandNum+(*nodeIterator)->getIntCandSerT(cCount)];
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*(univGenNum+univNodeNum+univCandLineNum)+univCandLineNum+scenCounter*univIntCandNum+(*nodeIterator)->getIntCandSerT(cCount) << "\t" << 1.0 << "\n";
			}
			busCount.push_back(rCount);
			if (((*nodeIterator)->devpinitMessage(scenCounter))==0) {
				modelSubnetMILP->addConstr(lhs[rCount], GRB_EQUAL, ((*nodeIterator)->devpinitMessage(scenCounter)));
			}
			else {
				modelSubnetMILP->addConstr(lhs[rCount], GRB_EQUAL, -((*nodeIterator)->devpinitMessage(scenCounter)));
			}
			outPutFile << "Connected load to node " << rCount << " is " << (*nodeIterator)->devpinitMessage(scenCounter)*100 << " MW" << endl;
			outPutFile << rCount << "\t";
			if (((*nodeIterator)->devpinitMessage(scenCounter))==0)
				outPutFile << ((*nodeIterator)->devpinitMessage(scenCounter))*100 << " MW" << endl;
			else
				outPutFile << -((*nodeIterator)->devpinitMessage(scenCounter))*100 << " MW" << endl;
			++rCount; // Increment the row count to point to the next node object
		}
	}
	// Coefficients corresponding to lower generation limits
	outPutFile << "\nCoefficients corresponding to lower generation limits\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
			lhs[rCount] = 0;
			lhs[rCount] += decvar[rCount - countOfScenarios*univNodeNum];
			modelSubnetMILP->addConstr(lhs[rCount] >= ((*genIterator)->getPMin()));
			outPutFile << rCount << "\t" << (rCount - countOfScenarios*univNodeNum) << "\t" << 1.0 << "\t" << (*genIterator)->getPMin() << endl;
			outPutFile << rCount << "\t";
			outPutFile << ((*genIterator)->getPMin())*100 << " MW" << endl;
			++rCount; // Increment the row count to point to the next generator object
		}
	}
	// Coefficients corresponding to upper generation limits
	outPutFile << "\nCoefficients corresponding to upper generation limits\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
			lhs[rCount] = 0;
			lhs[rCount] += decvar[rCount - countOfScenarios*(univGenNum + univNodeNum)];
			modelSubnetMILP->addConstr(lhs[rCount] <= ((*genIterator)->getPMax()));
			outPutFile << rCount << "\t" << (rCount - countOfScenarios*(univGenNum + univNodeNum)) << "\t" << 1.0 << "\t" << ((*genIterator)->getPMax()) << endl;
			outPutFile << rCount << "\t";
			outPutFile << ((*genIterator)->getPMax())*100 << " MW" << endl;
			++rCount; // Increment the row count to point to the next generator object
		}
	}
	// Coefficients corresponding to intra-zone Line Forward Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to intra-zone Line Forward Flow Limit Constraints\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (tranIterator = translObject.begin(); tranIterator != translObject.end(); ++tranIterator){
			lhs[rCount] = 0;
			lhs[rCount] += (1/((*tranIterator)->getReactance()))*(decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*tranIterator)->getTranslNodeID1()]);
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*tranIterator)->getTranslNodeID1() << "\t" << 1/((*tranIterator)->getReactance()) << "\t" << 1/((*tranIterator)->getReactance()) << "\n";
			lhs[rCount] += (-1/((*tranIterator)->getReactance()))*(decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*tranIterator)->getTranslNodeID2()]);
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*tranIterator)->getTranslNodeID2() << "\t" << -1/((*tranIterator)->getReactance()) << "\t" << "-" << "\t" << -1/((*tranIterator)->getReactance()) << "\n";
			modelSubnetMILP->addConstr(lhs[rCount] <= ((*tranIterator)->getFlowLimit()));
			outPutFile << rCount << "\t";
			outPutFile << ((*tranIterator)->getFlowLimit())*100 << " MW" << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		
		}
	}
	// Coefficients corresponding to intra-zone Line Reverse Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to intra-zone Line Reverse Flow Limit Constraints\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (tranIterator = translObject.begin(); tranIterator != translObject.end(); ++tranIterator){
			lhs[rCount] = 0;
			lhs[rCount] += (1/((*tranIterator)->getReactance()))*(decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*tranIterator)->getTranslNodeID1()]);
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*tranIterator)->getTranslNodeID1() << "\t" << 1/((*tranIterator)->getReactance()) << "\t" << 1/((*tranIterator)->getReactance()) << "\n";
			lhs[rCount] += (-1/((*tranIterator)->getReactance()))*(decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*tranIterator)->getTranslNodeID2()]);
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*tranIterator)->getTranslNodeID2() << "\t" << -1/((*tranIterator)->getReactance()) << "\t" << "-" << "\t" << -1/((*tranIterator)->getReactance()) << "\n";
			modelSubnetMILP->addConstr(lhs[rCount] >= -((*tranIterator)->getFlowLimit()));
			outPutFile << rCount << "\t";
			outPutFile << -((*tranIterator)->getFlowLimit())*100 << " MW" << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		}
	}	
	// Coefficients corresponding to shared existing Line Forward Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to shared existing Line Forward Flow Limit Constraints\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (exsharedIterator = SELineObject.begin(); exsharedIterator != SELineObject.end(); ++exsharedIterator){
			lhs[rCount] = 0;
			lhs[rCount] += (1/((*exsharedIterator)->getReactance()))*(decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*exsharedIterator)->getFromNodeID()]);
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*exsharedIterator)->getFromNodeID() << "\t" << 1/((*exsharedIterator)->getReactance()) << "\n";
			lhs[rCount] += (-1/((*exsharedIterator)->getReactance()))*(decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*exsharedIterator)->getToNodeID()]);
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*exsharedIterator)->getToNodeID() << "\t" << -1/((*exsharedIterator)->getReactance()) << "\n";
			modelSubnetMILP->addConstr(lhs[rCount] <= (*exsharedIterator)->getFlowLimit());
			outPutFile << rCount << "\t";
			outPutFile << ((*exsharedIterator)->getFlowLimit())*100 << " MW" << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		}
	}
	// Coefficients corresponding to shared existing Line Reverse Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to shared existing Line Reverse Flow Limit Constraints\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (exsharedIterator = SELineObject.begin(); exsharedIterator != SELineObject.end(); ++exsharedIterator){
			lhs[rCount] = 0;
			lhs[rCount] += (1/((*exsharedIterator)->getReactance()))*(decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*exsharedIterator)->getFromNodeID()]);
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*exsharedIterator)->getFromNodeID() << "\t" << 1/((*exsharedIterator)->getReactance()) << "\n";
			lhs[rCount] += (-1/((*exsharedIterator)->getReactance()))*(decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*exsharedIterator)->getToNodeID()]);
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*exsharedIterator)->getToNodeID() << "\t" << -1/((*exsharedIterator)->getReactance()) << "\n";
			modelSubnetMILP->addConstr(lhs[rCount] >= -((*exsharedIterator)->getFlowLimit()));
			outPutFile << rCount << "\t";
			outPutFile << -((*exsharedIterator)->getFlowLimit())*100 << " MW" << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		}
	}	
	// Coefficients corresponding to shared candidate Line Forward Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to shared candidate Line Forward Flow Limit Constraints\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			lhs[rCount] = 0;
			lhs[rCount] += decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum)];
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum) << "\t" << 1 << "\n";
			lhs[rCount] += (-((*candIterator)->getFlowLimit()))*(decvar[countOfScenarios*univCandLineNum+rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum)-scenCounter*univCandLineNum]);
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*univCandLineNum+rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum)-scenCounter*univCandLineNum << "\t" << -((*candIterator)->getFlowLimit()) << "\n";
			modelSubnetMILP->addConstr(lhs[rCount] <= 0.0);
			outPutFile << rCount << "\t";
			outPutFile << 0.0 << endl;
			++rCount;
		}
	}
	// Coefficients corresponding to shared candidate Line Reverse Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to shared candidate Line Reverse Flow Limit Constraints\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			lhs[rCount] = 0;
			lhs[rCount] += decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+univCandLineNum)];
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+univCandLineNum) << "\t" << 1 << "\n";
			lhs[rCount] += (((*candIterator)->getFlowLimit()))*(decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum)-scenCounter*univCandLineNum]);
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum)-scenCounter*univCandLineNum << "\t" << ((*candIterator)->getFlowLimit()) << "\n";
			modelSubnetMILP->addConstr(lhs[rCount] >= 0.0);
			outPutFile << rCount << "\t";
			outPutFile << 0.0 << endl;
			++rCount;
		}
	}	
	// Coefficients corresponding to shared candidate Line Definition upper bound
	outPutFile << "\nCoefficients corresponding to shared candidate Line Definition upper bound\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			lhs[rCount] = 0;
			lhs[rCount] += decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+2*univCandLineNum)];
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+2*univCandLineNum) << "\t" << 1 << "\n";
			lhs[rCount] += (-1/((*candIterator)->getReactance()))*(decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*candIterator)->getFromNodeID()]);
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*candIterator)->getFromNodeID() << "\t" << -1/((*candIterator)->getReactance()) << "\n";
			lhs[rCount] += (1/((*candIterator)->getReactance()))*(decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*candIterator)->getToNodeID()]);
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*candIterator)->getToNodeID() << "\t" << 1/((*candIterator)->getReactance()) << "\n";
			lhs[rCount] += (2.5*((*candIterator)->getFlowLimit()))*(decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+univCandLineNum)-scenCounter*univCandLineNum]);
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+univCandLineNum)-scenCounter*univCandLineNum << "\t" << BIGM << "\n";
			modelSubnetMILP->addConstr(lhs[rCount] <= (2.5*((*candIterator)->getFlowLimit())));//
			outPutFile << rCount << "\t";
			outPutFile << BIGM << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		}
	}
	// Coefficients corresponding to shared candidate Line Definition lower bound
	outPutFile << "\nCoefficients corresponding to shared candidate Line Definition lower bound\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			lhs[rCount] = 0;
			lhs[rCount] += decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum)];
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum) << "\t" << 1 << "\n";
			lhs[rCount] += (-1/((*candIterator)->getReactance()))*(decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*candIterator)->getFromNodeID()]);
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*candIterator)->getFromNodeID() << "\t" << -1/((*candIterator)->getReactance()) << "\n";
			lhs[rCount] += (1/((*candIterator)->getReactance()))*(decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*candIterator)->getToNodeID()]);
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*candIterator)->getToNodeID() << "\t" << 1/((*candIterator)->getReactance()) << "\n";
			lhs[rCount] += (-(2.5*((*candIterator)->getFlowLimit())))*(decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+2*univCandLineNum)-scenCounter*univCandLineNum]);
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+2*univCandLineNum)-scenCounter*univCandLineNum << "\t" << -BIGM << "\n";
			modelSubnetMILP->addConstr(lhs[rCount] >=  -(2.5*((*candIterator)->getFlowLimit())));
			outPutFile << rCount << "\t";
			outPutFile << -BIGM << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		}
	}
	// Coefficients corresponding to Internal candidate Line Forward Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to Internal candidate Line Forward Flow Limit Constraints\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			lhs[rCount] = 0;
			lhs[rCount] += decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum)+univCandLineNum];
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum)+univCandLineNum << "\t" << 1 << "\n";
			lhs[rCount] += (-((*intCandIterator)->getFlowLimit()))*(decvar[countOfScenarios*univIntCandNum+rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum)+univCandLineNum-scenCounter*univIntCandNum]);
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*univIntCandNum+rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum)+univCandLineNum-scenCounter*univIntCandNum << "\t" << -((*intCandIterator)->getFlowLimit()) << "\n";
			modelSubnetMILP->addConstr(lhs[rCount] <= 0.0);
			outPutFile << rCount << "\t";
			outPutFile << 0.0 << endl;
			++rCount;
		}
	}
	// Coefficients corresponding to Internal candidate Line Reverse Flow Limit Constraints\n";
	outPutFile << "\nCoefficients corresponding to Internal candidate Line Reverse Flow Limit Constraints\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			lhs[rCount] = 0;
			lhs[rCount] += decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum+univIntCandNum)+univCandLineNum];
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum+univIntCandNum)+univCandLineNum << "\t" << 1 << "\n";
			lhs[rCount] += (((*intCandIterator)->getFlowLimit()))*(decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum)+univCandLineNum-scenCounter*univIntCandNum]);
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum)+univCandLineNum-scenCounter*univIntCandNum << "\t" << ((*intCandIterator)->getFlowLimit()) << "\n";
			modelSubnetMILP->addConstr(lhs[rCount] >= 0.0);
			outPutFile << rCount << "\t";
			outPutFile << 0.0 << endl;
			++rCount;
		}
	}
	// Coefficients corresponding to Internal candidate Line Definition upper bound
	outPutFile << "\nCoefficients corresponding to Internal candidate Line Definition upper bound\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			lhs[rCount] = 0;
			lhs[rCount] += decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum+2*univIntCandNum)+univCandLineNum];
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum+2*univIntCandNum)+univCandLineNum << "\t" << 1 << "\n"; 
			lhs[rCount] += (-1/((*intCandIterator)->getReactance()))*(decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*intCandIterator)->getTranslNodeID1()]);
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*intCandIterator)->getTranslNodeID1() << "\t" << -1/((*intCandIterator)->getReactance()) << "\n";
			lhs[rCount] += (1/((*intCandIterator)->getReactance()))*(decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*intCandIterator)->getTranslNodeID2()]);
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*intCandIterator)->getTranslNodeID2() << "\t" << 1/((*intCandIterator)->getReactance()) << "\n";
			lhs[rCount] += ((2.5*((*intCandIterator)->getFlowLimit())))*(decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum+univIntCandNum)+univCandLineNum-scenCounter*univIntCandNum]);
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum+univIntCandNum)+univCandLineNum-scenCounter*univIntCandNum << "\t" << BIGM << "\n";
			modelSubnetMILP->addConstr(lhs[rCount] <= (2.5*((*intCandIterator)->getFlowLimit())));
			outPutFile << rCount << "\t";
			outPutFile << BIGM << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		}
	}
	// Coefficients corresponding to Internal candidate Line Definition lower bound
	outPutFile << "\nCoefficients corresponding to Internal candidate Line Definition lower bound\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			lhs[rCount] = 0;
			lhs[rCount] += decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum+3*univIntCandNum)+univCandLineNum];
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum+3*univIntCandNum)+univCandLineNum << "\t" << 1 << "\n";
			lhs[rCount] += (-1/((*intCandIterator)->getReactance()))*(decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*intCandIterator)->getTranslNodeID1()]);
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*intCandIterator)->getTranslNodeID1() << "\t" << -1/((*intCandIterator)->getReactance()) << "\n";
			lhs[rCount] += (1/((*intCandIterator)->getReactance()))*(decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*intCandIterator)->getTranslNodeID2()]);
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*intCandIterator)->getTranslNodeID2() << "\t" << 1/((*intCandIterator)->getReactance()) << "\n";
			lhs[rCount] += ((2.5*((*intCandIterator)->getFlowLimit())))*(decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum+2*univIntCandNum)+univCandLineNum-scenCounter*univIntCandNum]);
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum+2*univIntCandNum)+univCandLineNum-scenCounter*univIntCandNum << "\t" << -BIGM << "\n";
			modelSubnetMILP->addConstr(lhs[rCount] >=  -(2.5*((*intCandIterator)->getFlowLimit())));
			outPutFile << rCount << "\t";
			outPutFile << -BIGM << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		}
	}
	outPutFile << "\nConstraint bounds (rows) Specified" << endl;
	outPutFile << "\nTotal number of rows: " << rCount - 1 << endl;
	outPutFile << "\nCoefficient Matrix specified" << endl;
	clock_t end1 = clock(); // stop the timer
	double elapsed_secs1 = double(end1 - begin) / CLOCKS_PER_SEC; // Calculate the time required to populate the constraint matrix and objective coefficients
	outPutFile << "\nTotal time taken to define the rows, columns, objective and populate the coefficient matrix = " << elapsed_secs1 << " s " << endl;
	// RUN THE OPTIMIZATION SIMULATION ALGORITHM //
	cout << "\nSimulation in Progress. Wait !!! ....." << endl;
	modelSubnetMILP->optimize(); // Solves the optimization problem
	int stat = modelSubnetMILP->get(GRB_IntAttr_Status); // Outputs the solution status of the problem 

	// DISPLAY THE SOLUTION DETAILS //
	if (stat == GRB_INFEASIBLE){
		outPutFile << "\nThe solution to the problem is INFEASIBLE." << endl;
		cout << "\nThe solution to the problem is INFEASIBLE." << endl;
		delete modelSubnetMILP; // Free the memory of the GUROBI Problem Model
		return 0;
	} else if (stat == GRB_INF_OR_UNBD) {
		outPutFile << "\nNO FEASIBLE or BOUNDED solution to the problem exists." << endl;
		cout << "\nNO FEASIBLE or BOUNDED solution to the problem exists." << endl;
		delete modelSubnetMILP; // Free the memory of the GUROBI Problem Model
		return 0;
	} else if (stat == GRB_UNBOUNDED) {
		outPutFile << "\nThe solution to the problem is UNBOUNDED." << endl;
		cout << "\nThe solution to the problem is UNBOUNDED." << endl;
		delete modelSubnetMILP; // Free the memory of the GUROBI Problem Model
		return 0;
	} else if (stat == GRB_OPTIMAL) {
		outPutFile << "\nThe solution to the problem is OPTIMAL." << endl;
		cout << "\nThe solution to the problem is OPTIMAL." << endl;
	}

	//Get the Optimal Objective Value results//
	z = modelSubnetMILP->get(GRB_DoubleAttr_ObjVal);
	cout << "Objective value " << z << endl;
	for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
		z += (*genIterator)->getNLCost();
	}
	cout << "Objective value after adding NL cost " << z << endl;
	// Open separate output files for writing results of different variables
	string outIntAngFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/outputCent/outputAnglesResults/internalAngleGUROBI.txt";
	string outCandFlowFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/outputCent/candidateLinesResults/candFlowMWGUROBI.txt";
	string outCandDecFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/outputCent/candidateLinesResults/candLineDecisionGUROBI.txt";
	string outintCandFlowFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/outputCent/intcandLinesResults/intcandFlowMWGUROBI.txt";
	string outintCandLineDecFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/outputCent/intcandLinesResults/intcandLineDecisionGUROBI.txt";
	string outSEFlowFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/outputCent/SELinesResults/SELineDecisionGUROBI.txt";
	string outConstraintSatisfaction = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/outputCent/ConstraintSat/constraintSatGUROBI.txt";
	ofstream internalAngleOut(outIntAngFileName, ios::out); //switchStateOut
	ofstream candFlowMWOut(outCandFlowFileName, ios::out); //switchOnOut
	ofstream candLineDecisionOut(outCandDecFileName, ios::out); //switchOffOut
	ofstream intCandFlowMWOut(outintCandFlowFileName, ios::out);
	ofstream intCandLineDecisionOut(outintCandLineDecFileName, ios::out);
	ofstream SEFlowOut(outSEFlowFileName, ios::out);
	ofstream constraintOut(outConstraintSatisfaction, ios::out);
	outPutFile << "\nThe Optimal Objective value (Generation Dispatch  and Line Building Decision cost) is: " << z << endl;
	powerGenOut << "\nThe Optimal Objective value (Generation Dispatch  and Line Building Decision cost) is: " << z << endl;
	cout << "\nThe Optimal Objective value (Generation Dispatch  and Line Building Decision cost) is: " << z << endl;
	vector<double> x; // Vector for storing decision variable output 
	x.push_back(0); // Initialize the decision Variable vector

	//Display Power Generation
	powerGenOut << "\n****************** GENERATORS' POWER GENERATION LEVELS (MW) *********************" << endl;
	powerGenOut << "GENERATOR ID" << "\t" << "GENERATOR MW" << "\n";
	int arrayInd = 1;
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
			x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			powerGenOut << (*genIterator)->getGenID() << "\t" << ((decvar[arrayInd]).get(GRB_DoubleAttr_X))*100 << " MW" << endl;
			++arrayInd;
		}
	}
	powerGenOut << "Finished writing Power Generation" << endl;

	// Display Internal node voltage phase angle variables
	internalAngleOut << "\n****************** NODE VOLTAGE PHASE ANGLE VALUES *********************" << endl;
	internalAngleOut << "NODE ID" << "\t" << "VOLTAGE PHASE ANGLE" << "\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); ++nodeIterator){
			if ( (((*nodeIterator)->getSharedFlag()) == 1) || (((*nodeIterator)->getCandFlag()) == 1) ) {
				x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
				internalAngleOut << (*nodeIterator)->getNodeID() << "\t" << ((decvar[arrayInd]).get(GRB_DoubleAttr_X)) << endl;		
				++arrayInd;
			}
			else {
				x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
				internalAngleOut << (*nodeIterator)->getNodeID() << "\t" << ((decvar[arrayInd]).get(GRB_DoubleAttr_X)) << endl;		
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
			x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			candFlowMWOut << (*candIterator)->getTranslID() << "\t" << ((decvar[arrayInd]).get(GRB_DoubleAttr_X))*100 << " MW" << endl;
			++arrayInd;
		}
	}
	candFlowMWOut << "Finished writing Shared Candidate lines' Power Flow variables" << endl;

	// Display Shared Candidate lines' Construction Decisions
	candLineDecisionOut << "\n****************** SHARED CANDIDATE LINES CONSTRUCTION DECISIONS *********************" << endl;
	candLineDecisionOut << "CANDIDATE LINE ID" << "\t" << "CONSTRUCTION DECISIONS" << "\n";
	for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
		x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
		candLineDecisionOut << (*candIterator)->getTranslID() << "\t" << ((decvar[arrayInd]).get(GRB_DoubleAttr_X)) << endl;
		++arrayInd;
	}
	candLineDecisionOut << "Finished writing Shared Candidate lines' Construction decisions" << endl;

	// Display Internal Candidate lines' Power Flow variables
	intCandFlowMWOut << "\n****************** INTERNAL CANDIDATE LINES POWER FLOW VALUES *********************" << endl;
	intCandFlowMWOut << "CANDIDATE LINE ID" << "\t" << "MW FLOW VALUE" << "\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			intCandFlowMWOut << (*intCandIterator)->getTranslID() << "\t" << ((decvar[arrayInd]).get(GRB_DoubleAttr_X))*100 << " MW" << endl;
			++arrayInd;
		}
	}
	intCandFlowMWOut << "Finished writing Internal Candidate lines' Power Flow variables" << endl;

	// Display Internal Candidate lines' Construction Decisions
	intCandLineDecisionOut << "\n****************** INTERNAL CANDIDATE LINES CONSTRUCTION DECISIONS *********************" << endl;
	intCandLineDecisionOut << "CANDIDATE LINE ID" << "\t" << "CONSTRUCTION DECISIONS" << "\n";
	for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
		x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
		intCandLineDecisionOut << (*intCandIterator)->getTranslID() << "\t" << ((decvar[arrayInd]).get(GRB_DoubleAttr_X)) << endl;
		++arrayInd;
	}
	intCandLineDecisionOut << "Finished writing Internal Candidate lines' Construction decisions" << endl;

	// Display Shared Existing Transmission lines' Flows
	SEFlowOut << "\n****************** SHARED EXISTING TRANSMISSION LINES FLOWS *********************" << endl;
	SEFlowOut << "SCENARIO ID" << "\t" << "SHARED EXISTING TRANSMISSION LINE ID" << "\t" << "MW FLOW" << "\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (exsharedIterator = SELineObject.begin(); exsharedIterator != SELineObject.end(); ++exsharedIterator){
			SEFlowOut << scenCounter << "\t" << (*exsharedIterator)->getTranslID() << "\t" << (1/((*exsharedIterator)->getReactance()))*((decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*exsharedIterator)->getFromNodeID()]).get(GRB_DoubleAttr_X)-(decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*exsharedIterator)->getToNodeID()]).get(GRB_DoubleAttr_X))*100 << " MW" << endl;
		}
	}
	SEFlowOut << "Finished writing Shared Existing Transmission lines' MW Flows" << endl;

	// Display constraint satisfaction for testing purposes
	constraintOut << "\n****************** CONSTRAINT SATISFACTION FOR TESTING *********************" << endl;
	constraintOut << "POWER BALANCE CONSTRAINTS" << "\n";
	double LHS[dimRow+1];
	//Non-Zero entries of A matrix (Constraint/Coefficient matrix entries)//
	// Coefficients for the supply-demand balance constraints
	constraintOut << "\nNon-zero elements of A matrix" << endl;
	constraintOut << "\nRow Number\tColumn Number\tNon-zero Entry\tFrom Reactance\tToReactance" << endl;
	constraintOut << "\nCoefficients for the supply-demand balance constraints" << endl;
	// Dummy Constraint //
	LHS[0] = 0*((decvar[0]).get(GRB_DoubleAttr_X));
	rCount = 1; // Initialize the row count
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); ++nodeIterator){
			constraintOut << "\nGeneration\t" << rCount << "\n";
			int genListLength = (*nodeIterator)->getGenLength(); // get the number
			LHS[rCount]=0;
			constraintOut << "\nScenario Count: " << scenCounter << " node count: " << (*nodeIterator)->getNodeID() << " Conn Gen: " << genListLength << endl;
			for (int cCount = 1; cCount <= genListLength; ++cCount){
				constraintOut << "\nSerial: " << cCount << " Generator Serial: " << (*nodeIterator)->getGenSer(cCount) << endl;
				LHS[rCount] += 1*((decvar[scenCounter*univGenNum+(*nodeIterator)->getGenSer(cCount)]).get(GRB_DoubleAttr_X));
				constraintOut << "\nRowCount" << "\t" << "Column Count" << "\t" << "LHS Value" << endl;
				constraintOut << "\n" << rCount << "\t" << scenCounter*univGenNum+(*nodeIterator)->getGenSer(cCount) << "\t" << LHS[rCount] << endl;
			}
			constraintOut << "\nIntrazonal Node Angles\t" << rCount << "\n";
			if ((((*nodeIterator)->getToReact())+((*nodeIterator)->getFromReact()))>=0)
				LHS[rCount] -= (((*nodeIterator)->getToReact())+((*nodeIterator)->getFromReact()))*((decvar[countOfScenarios*univGenNum+rCount]).get(GRB_DoubleAttr_X));
			else
				LHS[rCount] += (((*nodeIterator)->getToReact())+((*nodeIterator)->getFromReact()))*((decvar[countOfScenarios*univGenNum+rCount]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Count" << "\t" << "Total Reactance Reciprocal" << "\t" << "From Reactance Reciprocal" << "\t" << "To Reactance Reciprocal" << "\t" << "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << countOfScenarios*univGenNum+rCount << "\t" << ((*nodeIterator)->getToReact())+((*nodeIterator)->getFromReact()) << "\t" << ((*nodeIterator)->getFromReact()) << "\t" << ((*nodeIterator)->getToReact()) << "\t" << LHS[rCount] << endl;
			constraintOut << "\nConnected Intrazonal Node Angles\t" << rCount << "\n";
			int connNodeListLength = (*nodeIterator)->getConNodeLength(); // get the number of nodes connected to this node
			for (int cCount = 1; cCount <= connNodeListLength; ++cCount){
				if ((((*nodeIterator)->getToReact())+((*nodeIterator)->getFromReact()))>=0)
					LHS[rCount] -= (((*nodeIterator)->getConnReact(cCount)))*((decvar[countOfScenarios*univGenNum+scenCounter*univNodeNum+((*nodeIterator)->getConnSer(cCount))]).get(GRB_DoubleAttr_X));
				else
					LHS[rCount] += (((*nodeIterator)->getConnReact(cCount)))*((decvar[countOfScenarios*univGenNum+scenCounter*univNodeNum+((*nodeIterator)->getConnSer(cCount))]).get(GRB_DoubleAttr_X));
				constraintOut << "\nRow Count" << "\t" << "Column Number" << "\t" <<  "Reciprocal Reactance" << "\t" << "LHS Value" << endl;
				constraintOut << "\n" << rCount << "\t" << countOfScenarios*univGenNum+scenCounter*univNodeNum+((*nodeIterator)->getConnSer(cCount)) << "\t" <<  (((*nodeIterator)->getConnReact(cCount))) << "\t" << LHS[rCount] << endl;

			}
			constraintOut << "\nConnected Candidate Lines for which this is the From node\t" << rCount << "\n";
			int connCandListLengthF = (*nodeIterator)->getCandLineLengthF(); // get the number of candidate lines connected to this from node 
			for (int cCount = 1; cCount <= connCandListLengthF; ++cCount){
				LHS[rCount] += (-1)*((decvar[countOfScenarios*(univGenNum+univNodeNum)+scenCounter*univCandLineNum+(*nodeIterator)->getCandSerF(cCount)]).get(GRB_DoubleAttr_X));
				constraintOut << "\nRow Count" << "\t" << "Column Number" << "\t" <<  "LHS Value" << endl;
				constraintOut << "\n" << rCount << "\t" << countOfScenarios*(univGenNum+univNodeNum)+scenCounter*univCandLineNum+(*nodeIterator)->getCandSerF(cCount) << "\t" << LHS[rCount] << endl;
			}
			constraintOut << "\nConnected Candidate Lines for which this is the To node\t" << rCount << "\n";
			int connCandListLengthT = (*nodeIterator)->getCandLineLengthT(); // get the number of candidate lines connected to this to node 
			for (int cCount = 1; cCount <= connCandListLengthT; ++cCount){
				LHS[rCount] += ((decvar[countOfScenarios*(univGenNum+univNodeNum)+scenCounter*univCandLineNum+(*nodeIterator)->getCandSerT(cCount)]).get(GRB_DoubleAttr_X));
				constraintOut << "\nRow Count" << "\t" << "Column Number" << "\t" <<  "LHS Value" << endl;
				constraintOut << "\n" << rCount << "\t" << countOfScenarios*(univGenNum+univNodeNum)+scenCounter*univCandLineNum+(*nodeIterator)->getCandSerT(cCount) << "\t" << LHS[rCount] << endl;
			}
			constraintOut << "\nConnected Internal Candidate Lines for which this is the From node\t" << rCount << "\n";
			int connintCandListLengthF = (*nodeIterator)->getIntCandLineLengthF(); // get the number of internal candidate lines connected to this from node 
			for (int cCount = 1; cCount <= connintCandListLengthF; ++cCount){
				LHS[rCount] += (-1)*((decvar[countOfScenarios*(univGenNum+univNodeNum+univCandLineNum)+univCandLineNum+scenCounter*univIntCandNum+(*nodeIterator)->getIntCandSerF(cCount)]).get(GRB_DoubleAttr_X));
				constraintOut << "\nRow Count" << "\t" << "Column Number" << "\t" <<  "LHS Value" << endl;
				constraintOut << "\n" << rCount << "\t" << countOfScenarios*(univGenNum+univNodeNum+univCandLineNum)+univCandLineNum+scenCounter*univIntCandNum+(*nodeIterator)->getIntCandSerF(cCount) << "\t" << LHS[rCount] << endl;
			}
			constraintOut << "\nConnected Internal Candidate Lines for which this is the To node\t" << rCount << "\n";
			int connintCandListLengthT = (*nodeIterator)->getIntCandLineLengthT(); // get the number of internal candidate lines connected to this to node 
			for (int cCount = 1; cCount <= connintCandListLengthT; ++cCount){
				LHS[rCount] += ((decvar[countOfScenarios*(univGenNum+univNodeNum+univCandLineNum)+univCandLineNum+scenCounter*univIntCandNum+(*nodeIterator)->getIntCandSerT(cCount)]).get(GRB_DoubleAttr_X));
				constraintOut << "\nRow Count" << "\t" << "Column Number" << "\t" <<  "LHS Value" << endl;
				constraintOut << "\n" << rCount << "\t" << countOfScenarios*(univGenNum+univNodeNum+univCandLineNum)+univCandLineNum+scenCounter*univIntCandNum+(*nodeIterator)->getIntCandSerT(cCount) << "\t" << LHS[rCount] << endl;
			}
			constraintOut << "Connected load to node " << rCount << " is " << (*nodeIterator)->devpinitMessage(scenCounter)*100 << " MW" << endl;
			constraintOut << "RHS Value" << endl;
			if (((*nodeIterator)->devpinitMessage(scenCounter))==0)
				constraintOut << ((*nodeIterator)->devpinitMessage(scenCounter)) << endl;
			else
				constraintOut << -((*nodeIterator)->devpinitMessage(scenCounter)) << endl;
			++rCount; // Increment the row count to point to the next node object
		}
	}
	// Coefficients corresponding to lower generation limits
	constraintOut << "LOWER GENERATION LIMIT CONSTRAINTS" << "\n";
	constraintOut << "\nCoefficients corresponding to lower generation limits\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
			LHS[rCount] = 0;
			LHS[rCount] += ((decvar[rCount - countOfScenarios*univNodeNum]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number" << "\t" <<  "LHS Value" << "\t" <<  "RHS Value" << endl;
			constraintOut << rCount << "\t" << (rCount - countOfScenarios*univNodeNum) << "\t" << LHS[rCount] << ">=" << (*genIterator)->getPMin() << endl;
			constraintOut << "PMin" << "\t";
			constraintOut << ((*genIterator)->getPMin())*100 << " MW" << endl;
			++rCount; // Increment the row count to point to the next generator object
		}
	}
	// Coefficients corresponding to upper generation limits
	constraintOut << "UPPER GENERATION LIMIT CONSTRAINTS" << "\n";
	constraintOut << "\nCoefficients corresponding to upper generation limits\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
			LHS[rCount] = 0;
			LHS[rCount] += ((decvar[rCount - countOfScenarios*(univGenNum + univNodeNum)]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number" << "\t" <<  "LHS Value" << "\t" <<  "RHS Value" << endl;
			constraintOut << rCount << "\t" << (rCount - countOfScenarios*(univGenNum + univNodeNum)) << "\t" << LHS[rCount] << "<=" << ((*genIterator)->getPMax()) << endl;
			constraintOut << "PMax" << "\t";
			constraintOut << ((*genIterator)->getPMax())*100 << " MW" << endl;
			++rCount; // Increment the row count to point to the next generator object
		}
	}
	// Coefficients corresponding to intra-zone Line Forward Flow Limit Constraints
	constraintOut << "INTRA-ZONE TRANSMISSION LINES' FORWARD FLOW LIMITS CONSTRAINTS" << "\n";
	constraintOut << "\nCoefficients corresponding to intra-zone Line Forward Flow Limit Constraints\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (tranIterator = translObject.begin(); tranIterator != translObject.end(); ++tranIterator){
			LHS[rCount] = 0;
			LHS[rCount] += (1/((*tranIterator)->getReactance()))*((decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*tranIterator)->getTranslNodeID1()]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(FromEnd)" << "\t" <<  "Reciprocal Reactance" << "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*tranIterator)->getTranslNodeID1() << "\t" << 1/((*tranIterator)->getReactance()) << "\t" << LHS[rCount] << "\n";
			LHS[rCount] += (-1/((*tranIterator)->getReactance()))*((decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*tranIterator)->getTranslNodeID2()]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(ToEnd)" << "\t" <<  "Reciprocal Reactance" << "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*tranIterator)->getTranslNodeID2() << "\t" << -1/((*tranIterator)->getReactance()) << "\t" << LHS[rCount] << "\n";
			constraintOut << "Flow Limit" << "\t";
			constraintOut << ((*tranIterator)->getFlowLimit())*100 << " MW" << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		
		}
	}
	// Coefficients corresponding to intra-zone Line Reverse Flow Limit Constraints
	constraintOut << "INTRA-ZONE TRANSMISSION LINES' REVERSE FLOW LIMITS CONSTRAINTS" << "\n";
	constraintOut << "\nCoefficients corresponding to intra-zone Line Reverse Flow Limit Constraints\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (tranIterator = translObject.begin(); tranIterator != translObject.end(); ++tranIterator){
			LHS[rCount] = 0;
			LHS[rCount] += (1/((*tranIterator)->getReactance()))*((decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*tranIterator)->getTranslNodeID1()]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(FromEnd)" << "\t" <<  "Reciprocal Reactance" << "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*tranIterator)->getTranslNodeID1() << "\t" << 1/((*tranIterator)->getReactance()) << "\t" << LHS[rCount] << "\n";
			LHS[rCount] += (-1/((*tranIterator)->getReactance()))*((decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*tranIterator)->getTranslNodeID2()]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(ToEnd)" << "\t" <<  "Reciprocal Reactance" << "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*tranIterator)->getTranslNodeID2() << "\t" << -1/((*tranIterator)->getReactance()) << "\t" << LHS[rCount] << "\n";
			constraintOut << "Reverse Flow Limit" << "\t";
			constraintOut << -((*tranIterator)->getFlowLimit())*100 << " MW" << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		}
	}	
	// Coefficients corresponding to shared existing Line Forward Flow Limit Constraints
	constraintOut << "SHARED EXISTING TRANSMISSION LINES' FORWARD FLOW LIMITS CONSTRAINTS" << "\n";
	constraintOut << "\nCoefficients corresponding to shared existing Line Forward Flow Limit Constraints\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (exsharedIterator = SELineObject.begin(); exsharedIterator != SELineObject.end(); ++exsharedIterator){
			LHS[rCount] = 0;
			LHS[rCount] += (1/((*exsharedIterator)->getReactance()))*((decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*exsharedIterator)->getFromNodeID()]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(FromEnd)" << "\t" <<  "Reciprocal Reactance" << "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*exsharedIterator)->getFromNodeID() << "\t" << 1/((*exsharedIterator)->getReactance()) << "\t" << LHS[rCount] << "\n";
			LHS[rCount] += (-1/((*exsharedIterator)->getReactance()))*((decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*exsharedIterator)->getToNodeID()]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(ToEnd)" << "\t" <<  "Reciprocal Reactance" << "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*exsharedIterator)->getToNodeID() << "\t" << -1/((*exsharedIterator)->getReactance()) << "\t" << LHS[rCount] << "\n";
			constraintOut << "Flow Limit" << "\t";
			constraintOut << ((*exsharedIterator)->getFlowLimit())*100 << " MW" << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		}
	}
	// Coefficients corresponding to shared existing Line Reverse Flow Limit Constraints
	constraintOut << "SHARED EXISTING TRANSMISSION LINES' REVERSE FLOW LIMITS CONSTRAINTS" << "\n";
	constraintOut << "\nCoefficients corresponding to shared existing Line Reverse Flow Limit Constraints\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (exsharedIterator = SELineObject.begin(); exsharedIterator != SELineObject.end(); ++exsharedIterator){
			LHS[rCount] = 0;
			LHS[rCount] += (1/((*exsharedIterator)->getReactance()))*((decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*exsharedIterator)->getFromNodeID()]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(FromEnd)" << "\t" <<  "Reciprocal Reactance" << "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*exsharedIterator)->getFromNodeID() << "\t" << 1/((*exsharedIterator)->getReactance()) << "\t" << LHS[rCount] << "\n";
			LHS[rCount] += (-1/((*exsharedIterator)->getReactance()))*((decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*exsharedIterator)->getToNodeID()]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(ToEnd)" << "\t" <<  "Reciprocal Reactance" << "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*exsharedIterator)->getToNodeID() << "\t" << -1/((*exsharedIterator)->getReactance()) << "\t" << LHS[rCount] << "\n";
			constraintOut << "Reverse Flow Limit" << "\t";
			constraintOut << -((*exsharedIterator)->getFlowLimit())*100 << " MW" << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		}
	}	
	// Coefficients corresponding to shared candidate Line Forward Flow Limit Constraints
	constraintOut << "SHARED CANDIDATE TRANSMISSION LINES' FORWARD FLOW LIMITS CONSTRAINTS" << "\n";
	constraintOut << "\nCoefficients corresponding to shared candidate Line Forward Flow Limit Constraints\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			LHS[rCount] = 0;
			LHS[rCount] += ((decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum)]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(Flow)" << "\t" <<  "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum) << "\t" << LHS[rCount] << "\n";
			LHS[rCount] += (-((*candIterator)->getFlowLimit()))*((decvar[countOfScenarios*univCandLineNum+rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum)-scenCounter*univCandLineNum]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(Integer)" << "\t" <<  "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << countOfScenarios*univCandLineNum+rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum)-scenCounter*univCandLineNum << "\t" << LHS[rCount] << "\n";
			constraintOut << LHS[rCount] << "<=" << 0.0 << endl;
			++rCount;
		}
	}
	// Coefficients corresponding to shared candidate Line Reverse Flow Limit Constraints
	constraintOut << "SHARED CANDIDATE TRANSMISSION LINES' REVERSE FLOW LIMITS CONSTRAINTS" << "\n";
	constraintOut << "\nCoefficients corresponding to shared candidate Line Reverse Flow Limit Constraints\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			LHS[rCount] = 0;
			LHS[rCount] += ((decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+univCandLineNum)]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(Flow)" << "\t" <<  "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+univCandLineNum) << "\t" << LHS[rCount] << "\n";
			LHS[rCount] += (((*candIterator)->getFlowLimit()))*((decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum)-scenCounter*univCandLineNum]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(Integer)" << "\t" <<  "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum)-scenCounter*univCandLineNum << "\t" << LHS[rCount] << "\n";
			constraintOut << LHS[rCount] << ">=" << 0.0 << endl;
			++rCount;
		}
	}	
	// Coefficients corresponding to shared candidate Line Definition upper bound
	constraintOut << "SHARED CANDIDATE TRANSMISSION LINES' DEFINITION UPPER BOUNDS CONSTRAINTS" << "\n";
	constraintOut << "\nCoefficients corresponding to shared candidate Line Definition upper bound\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			LHS[rCount] = 0;
			LHS[rCount] += ((decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+2*univCandLineNum)]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(Flow)" << "\t" <<  "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+2*univCandLineNum) << "\t" << LHS[rCount] << "\n";
			LHS[rCount] += (-1/((*candIterator)->getReactance()))*((decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*candIterator)->getFromNodeID()]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(From Node)" << "\t" << "From Reciprocal Reactance" << "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*candIterator)->getFromNodeID() << "\t" << -1/((*candIterator)->getReactance()) << "\t" << LHS[rCount] << "\n";
			LHS[rCount] += (1/((*candIterator)->getReactance()))*((decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*candIterator)->getToNodeID()]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(To Node)" << "\t" << "To Reciprocal Reactance" << "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*candIterator)->getToNodeID() << "\t" << 1/((*candIterator)->getReactance()) << "\t" << LHS[rCount] << "\n";
			LHS[rCount] += (2.5*((*candIterator)->getFlowLimit()))*((decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+univCandLineNum)-scenCounter*univCandLineNum]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(Integer)" << "\t" << "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+univCandLineNum)-scenCounter*univCandLineNum << "\t" << LHS[rCount] << "\n";
			constraintOut << LHS[rCount] << "<=" << (2.5*((*candIterator)->getFlowLimit())) << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		}
	}
	// Coefficients corresponding to shared candidate Line Definition lower bound
	constraintOut << "SHARED CANDIDATE TRANSMISSION LINES' DEFINITION LOWER BOUNDS CONSTRAINTS" << "\n";
	constraintOut << "\nCoefficients corresponding to shared candidate Line Definition lower bound\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			LHS[rCount] = 0;
			LHS[rCount] += ((decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum)]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(Flow)" << "\t" <<  "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum) << "\t" << LHS[rCount] << "\n";
			LHS[rCount] += (-1/((*candIterator)->getReactance()))*((decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*candIterator)->getFromNodeID()]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(From Node)" << "\t" << "From Reciprocal Reactance" << "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*candIterator)->getFromNodeID() << "\t" << -1/((*candIterator)->getReactance()) << "\t" << LHS[rCount] << "\n";
			LHS[rCount] += (1/((*candIterator)->getReactance()))*((decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*candIterator)->getToNodeID()]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(Integer)" << "\t" <<  "To Reciprocal Reactance" << "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*candIterator)->getToNodeID() << "\t" << 1/((*candIterator)->getReactance()) << "\t" << LHS[rCount] << "\n";
			LHS[rCount] += (-(2.5*((*candIterator)->getFlowLimit())))*((decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+2*univCandLineNum)-scenCounter*univCandLineNum]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(Integer)" << "\t" << "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+2*univCandLineNum)-scenCounter*univCandLineNum << "\t" << LHS[rCount] << "\n";
			constraintOut << LHS[rCount] << ">=" << -(2.5*((*candIterator)->getFlowLimit())) << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		}
	}
	// Coefficients corresponding to Internal candidate Line Forward Flow Limit Constraints
	constraintOut << "INTERNAL CANDIDATE TRANSMISSION LINES' FORWARD FLOW LIMITS CONSTRAINTS" << "\n";
	constraintOut << "\nCoefficients corresponding to Internal candidate Line Forward Flow Limit Constraints\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			LHS[rCount] = 0;
			LHS[rCount] += ((decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum)+univCandLineNum]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(Flow)" << "\t" <<  "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum)+univCandLineNum << "\t" << LHS[rCount] << "\n";
			LHS[rCount] += (-((*intCandIterator)->getFlowLimit()))*((decvar[countOfScenarios*univIntCandNum+rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum)+univCandLineNum-scenCounter*univIntCandNum]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(Integer)" << "\t" <<  "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << countOfScenarios*univIntCandNum+rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum)+univCandLineNum-scenCounter*univIntCandNum << "\t" << LHS[rCount] << "\n";
			constraintOut << LHS[rCount] << "<=" << 0.0 << endl;
			++rCount;
		}
	}
	// Coefficients corresponding to Internal candidate Line Reverse Flow Limit Constraints\n";
	constraintOut << "INTERNAL CANDIDATE TRANSMISSION LINES' REVERSE FLOW LIMITS CONSTRAINTS" << "\n";
	constraintOut << "\nCoefficients corresponding to Internal candidate Line Reverse Flow Limit Constraints\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			LHS[rCount] = 0;
			LHS[rCount] += ((decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum+univIntCandNum)+univCandLineNum]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(Flow)" << "\t" <<  "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum+univIntCandNum)+univCandLineNum << "\t" << LHS[rCount] << "\n";
			LHS[rCount] += (((*intCandIterator)->getFlowLimit()))*((decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum)+univCandLineNum-scenCounter*univIntCandNum]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(Integer)" << "\t" <<  "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum)+univCandLineNum-scenCounter*univIntCandNum << "\t" << LHS[rCount] << "\n";
			constraintOut << LHS[rCount] << ">=" << 0.0 << endl;
			++rCount;
		}
	}
	// Coefficients corresponding to Internal candidate Line Definition upper bound
	constraintOut << "INTERNAL CANDIDATE TRANSMISSION LINES' DEFINITION UPPER BOUNDS CONSTRAINTS" << "\n";
	constraintOut << "\nCoefficients corresponding to Internal candidate Line Definition upper bound\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			LHS[rCount] = 0;
			LHS[rCount] += ((decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum+2*univIntCandNum)+univCandLineNum]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(Flow)" << "\t" <<  "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum+2*univIntCandNum)+univCandLineNum << "\t" << LHS[rCount] << "\n"; 
			LHS[rCount] += (-1/((*intCandIterator)->getReactance()))*((decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*intCandIterator)->getTranslNodeID1()]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(From Node)" << "\t" << "From Reciprocal Reactance" << "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*intCandIterator)->getTranslNodeID1() << "\t" << -1/((*intCandIterator)->getReactance()) << "\t" << LHS[rCount] << "\n";
			LHS[rCount] += (1/((*intCandIterator)->getReactance()))*((decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*intCandIterator)->getTranslNodeID2()]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(To Node)" << "\t" << "To Reciprocal Reactance" << "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*intCandIterator)->getTranslNodeID2() << "\t" << 1/((*intCandIterator)->getReactance()) << "\t" << LHS[rCount] << "\n";
			LHS[rCount] += ((2.5*((*intCandIterator)->getFlowLimit())))*((decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum+univIntCandNum)+univCandLineNum-scenCounter*univIntCandNum]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(Integer)" << "\t" << "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum+univIntCandNum)+univCandLineNum-scenCounter*univIntCandNum << "\t" << LHS[rCount] << "\n";
			constraintOut << LHS[rCount] << "<=" << (2.5*((*intCandIterator)->getFlowLimit())) << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		}
	}
	// Coefficients corresponding to Internal candidate Line Definition lower bound
	constraintOut << "INTERNAL CANDIDATE TRANSMISSION LINES' DEFINITION LOWER BOUNDS CONSTRAINTS" << "\n";
	constraintOut << "\nCoefficients corresponding to Internal candidate Line Definition lower bound\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			LHS[rCount] = 0;
			LHS[rCount] += ((decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum+3*univIntCandNum)+univCandLineNum]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(Flow)" << "\t" <<  "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum+3*univIntCandNum)+univCandLineNum << "\t" << LHS[rCount] << "\n"; 
			LHS[rCount] += (-1/((*intCandIterator)->getReactance()))*((decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*intCandIterator)->getTranslNodeID1()]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(From Node)" << "\t" << "From Reciprocal Reactance" << "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*intCandIterator)->getTranslNodeID1() << "\t" << -1/((*intCandIterator)->getReactance()) << "\t" << LHS[rCount] << "\n"; 
			LHS[rCount] += (1/((*intCandIterator)->getReactance()))*((decvar[countOfScenarios*univGenNum + scenCounter*univNodeNum+(*intCandIterator)->getTranslNodeID2()]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(To Node)" << "\t" << "To Reciprocal Reactance" << "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << countOfScenarios*univGenNum + scenCounter*univNodeNum+(*intCandIterator)->getTranslNodeID2() << "\t" << 1/((*intCandIterator)->getReactance()) << "\t" << LHS[rCount] << "\n"; 
			LHS[rCount] += ((2.5*((*intCandIterator)->getFlowLimit())))*((decvar[rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum+2*univIntCandNum)+univCandLineNum-scenCounter*univIntCandNum]).get(GRB_DoubleAttr_X));
			constraintOut << "\nRow Count" << "\t" << "Column Number(Integer)" << "\t" << "LHS Value" << endl;
			constraintOut << "\n" << rCount << "\t" << rCount-countOfScenarios*(univGenNum+2*univTranNum+2*univSELineNum+3*univCandLineNum+2*univIntCandNum)+univCandLineNum-scenCounter*univIntCandNum << "\t" << LHS[rCount] << "\n"; 
			constraintOut << LHS[rCount] << ">=" << -(2.5*((*intCandIterator)->getFlowLimit())) << endl;
			++rCount; // Increment the row count to point to the next transmission line object
		}
	}
	constraintOut << "OBJECTIVE TESTING" << "\n";
	double OBJ = 0.0;
	// Objective Contribution from Dummy Decision Variable //
	OBJ += 0*((decvar[0]).get(GRB_DoubleAttr_X));
	colCount = 1;
	constraintOut << "Columns corresponding to Power Generation continuous variables for different generators" << endl;
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
			OBJ += (probability.at(scenCounter))*(((*genIterator)->getLinCoeff()))*((decvar[colCount]).get(GRB_DoubleAttr_X));
			constraintOut << "Generator Objective Value at scenario " << scenCounter << " for generator " << (*genIterator)->getGenID() << " is " << (probability.at(scenCounter))*(((*genIterator)->getLinCoeff()))*((decvar[colCount]).get(GRB_DoubleAttr_X)) << " and cumulative objective: " << OBJ << " Column count: " << colCount << endl;
			constraintOut << "Scenario Probability is " << (probability.at(scenCounter)) << " Generator Cost Coefficient is " << (((*genIterator)->getLinCoeff())) << " Generation Decision is " << ((decvar[colCount]).get(GRB_DoubleAttr_X)) << endl;
			++colCount;
		}
	}
	constraintOut << "Columns corresponding to Voltage Phase Angles continuous variables for different nodes" << endl;
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {	
		for (nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); ++nodeIterator){
			OBJ += 0*((decvar[colCount]).get(GRB_DoubleAttr_X));
			constraintOut << "Angle Value at scenario " << scenCounter << " for node " << (*nodeIterator)->getNodeID() << " is " << ((decvar[colCount]).get(GRB_DoubleAttr_X)) << " and cumulative objective: " << OBJ << " Column count: " << colCount << endl;
			++colCount;
		}
	}
	constraintOut << "Columns corresponding to Shared Candidate Line Flows continuous variables" << endl;
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			OBJ += 0*((decvar[colCount]).get(GRB_DoubleAttr_X));
			constraintOut << "shared candidate Line Flow at scenario " << scenCounter << " for line " << (*candIterator)->getTranslID() << " is " << ((decvar[colCount]).get(GRB_DoubleAttr_X)) << " and cumulative objective: " << OBJ << " Column count: " << colCount << endl;
			++colCount;
		}
	}
	constraintOut << "Columns corresponding to Shared Candidate Line Construction Decision Binary Integer variables" << endl;
	for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
		OBJ += (2*((*candIterator)->getInvestCost()))*((decvar[colCount]).get(GRB_DoubleAttr_X));
		constraintOut << "shared candidate Line construction decision for line " << (*candIterator)->getTranslID() << " is " << ((decvar[colCount]).get(GRB_DoubleAttr_X)) << " Cost is " << (2*((*candIterator)->getInvestCost()))*((decvar[colCount]).get(GRB_DoubleAttr_X)) << " and cumulative objective: " << OBJ << " Column count: " << colCount << endl;	
		++colCount;
	}

	constraintOut << "Columns corresponding to Internal Candidate Line Flows continuous variables" << endl;
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			OBJ += 0*((decvar[colCount]).get(GRB_DoubleAttr_X));
			constraintOut << "internal candidate Line Flow at scenario " << scenCounter << " for line " << (*intCandIterator)->getTranslID() << " is " << ((decvar[colCount]).get(GRB_DoubleAttr_X)) << " and cumulative objective: " << OBJ << " Column count: " << colCount << endl;
			++colCount;
		}
	}
	constraintOut << "Columns corresponding to Internal Candidate Line Construction Decision Binary Integer variables" << endl;
	for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
		OBJ += ((*intCandIterator)->getInvestCost())*((decvar[colCount]).get(GRB_DoubleAttr_X));
		constraintOut << "internal candidate Line construction decision for line " << (*intCandIterator)->getTranslID() << " is " << ((decvar[colCount]).get(GRB_DoubleAttr_X)) << " Cost is " << ((*intCandIterator)->getInvestCost())*((decvar[colCount]).get(GRB_DoubleAttr_X)) << " and cumulative objective: " << OBJ << " Column count: " << colCount << endl;		
		++colCount;
	}

	delete modelSubnetMILP; // Free the memory of the GUROBI Problem Model
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
	intCandFlowMWOut.close();
	intCandLineDecisionOut.close();
	SEFlowOut.close();
	constraintOut.close();
	cout << "\nSimulation Completed.\nResults written on the different output files" << endl;
	return z;
} // Function MILP() ends



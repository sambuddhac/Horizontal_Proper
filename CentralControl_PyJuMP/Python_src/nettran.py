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

			#/* Generators */

			#/* Instantiate Generators */
			#Open the .txt file to read the Powergenerator parameter values
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
	assignProb()
	#end constructor
	def assignProb(self):
		for f in range(self.countOfScenarios):
			self.probability.append(1/self.countOfScenarios)

	def MILPAvgHRGUROBI(self): #Function MILPAvgHRGUROBI() implements the Mixed Integer Linear Programming Unit Commitment Solver routine by calling GUROBI routines for average heat rate objective for Horizontal Coordination Investment decision making
 		#Function MILP() ends



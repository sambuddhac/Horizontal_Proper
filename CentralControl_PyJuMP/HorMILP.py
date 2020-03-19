# HorMILP.cpp : Defines the entry point for the centralized/merged regions/Single control area Horizontal Investment Coordination MILP Simulation application.
# Main Method for running the Horizontal Investment Coordination MILP Simulation; Parts of the code follows the code design philosophy of Nick Laws of NREL
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
julSol.include(os.path.join("JuMP_src", "HorMILPCentralized.jl")) # definition of Gensolver class for base case scenario first interval
log.info(("Julia took {:..2f} seconds to start and include Horizontal Investment Coordination centralized models.".format(profiler.get_interval())))

def HorMILPCentral(): # Main method begins program execution
	'''Future Work
	#Choose the type of objective function
    '''
	systemChoice = int(input("Choose the type of System to be simulated: 1 for Simple two bus/two region, 2 for system combined of IEEE 14, 30, and 5 node systems"))
	curveChoice = 1 # Number to indicate the type of Objective function among average heat rate, piecewise linear, or polynomial; Assume Average Heat Rate for now
	# Read the master zones file, for deciding upon which other files to read for building the model
	int numberOfZones; // Number of zones between which horizontal investment coordination for transmission lines to be built is considered
	int numberOfFields; // Number of rows or individual file types for each of the zones
	if systemChoice==1:
		inputMasterFile = "masterZonesSummary.json"
	else:
		inputMasterFile = "masterZonesSummaryRevised.json"
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
	cout << endl << "\n*** NETWORK INITIALIZATION STAGE BEGINS ***\n" << endl << endl;
	Nettran *nettranInstance = new Nettran( str, numberOfZones, curveChoice ); // create the network instances for the different zones
	cout << "\n*** NETWORK INITIALIZATION STAGE ENDS: ZONAL SUB-NETWORKS CREATED ***\n" << endl;
	cout << endl << "\n*** SOLUTION OF SINGLE AREA MILP HORIZONTAL COORDINATION BEGINS ***\n" << endl << endl;
	switch (curveChoice) {
		case 1:
			{cout << "\nSOLVING MILP" << endl;
			double solution = nettranInstance->MILPAvgHRGUROBI(environmentGUROBI); // Perform unit commitment for average heat rate objective
			cout << "\nMILP SOLVED" << endl;
			break;}
		case 2:
			nettranInstance->MILPPiecewiseLin(); // Perform unit commitment for piecewise linear objective
			break;
		case 3:
			nettranInstance->MILPPolynomial(); // Perform unit commitment for polynomial objective
			break;
		default:
			cout << "\nInvalid choice of Objective function" << endl;
			break;
	}
	cout << endl << "\n*** SOLUTION OF SINGLE AREA MILP HORIZONTAL COORDINATION ENDS ***\n" << endl << endl;
	delete nettranInstance; // Free the memory of the Nettran class object
	delete environmentGUROBI; // Free the memory of the GUROBI environment object
	return 0; // Indicates successful Program Completion
} // End of main method


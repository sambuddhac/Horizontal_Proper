// Definition for Marketover class public Member Methods
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cmath>
#include <cstring>
#include <glpk.h> // Includes the GLPK (GNU Linear Programming Kit) header file
#include "gurobi_c++.h"
#include "marketOverseer.h" // includes definition of Nettran class
#include "nettran.h"
#define BIGM 1000000000000000000 // Defines the value of the Big M for transforming the bilinear terms to linear constraints

using namespace std;

Marketover::Marketover(int zoneCount, int totNodes, int SECount, int candCount, vector <Nettran*>* vectorOfSubnet, vector <int>* listOfSharedNode, vector <int>* listOfGlobalSharedNode, vector <int>* listOfSharedZone, vector <int>* ListOfSESerial, vector <int>* ListOfSEFromRank, vector <int>* ListOfSEToRank, vector <double>* ListOfSEReact, vector <double>* ListOfSECap, vector <int>* ListOfCandSerial, vector <int>* ListOfCandFromRank, vector <int>* ListOfCandToRank, vector <double>* ListOfCandReact, vector <double>* ListOfCandCap, int milpAlgoChoice, int scenarios, int compareBase) // constructor
{
	zoneNumber = zoneCount; // Initialize the total number of load zones
	subNetVector = *vectorOfSubnet; // Create handle to the vector of subnetworks
	zoneList = *listOfSharedZone; // List of zones between which exsiting and candidate shared lines exist
	nodeList = *listOfSharedNode; // List of nodes at the ends of the existing and candidate shared lines
	sharedGlobalList = *listOfGlobalSharedNode; 
	nodeNumber = totNodes; // total number of nodes, which form the ends of the shared lines
	candLineNumber = candCount; // Number of shared candidate lines
	sharedELineNumber = SECount; // Number of shared existing Transmission lines
	SESerial = *ListOfSESerial; // Serial list of shared existing lines
	SEFromRank = *ListOfSEFromRank; // List of rank of from nodes of shared existing lines
	SEToRank = *ListOfSEToRank; // List of rank of to nodes of shared existing lines
	SEReactance = *ListOfSEReact; // List of reactances of shared existing lines
	SECapacity = *ListOfSECap; // List of line flow limits of shared existing lines
	candSerial = *ListOfCandSerial; // Serial list of shared candidate lines
	candFromRank = *ListOfCandFromRank; // List of rank of from nodes of shared candidate lines
	candToRank = *ListOfCandToRank; // List of rank of to nodes of shared candidate lines
	candReactance = *ListOfCandReact; // List of reactances of shared candidate lines
	candCapacity = *ListOfCandCap; // List of line flow limits of shared candidate lines
	lpSolveAlgo = milpAlgoChoice; // Simplex for 1 and IPM for 2
	countOfScenarios = scenarios; // Number of total random scenarios
	//vector<double> upperBoundVector = ; // Vector of upper bounds by iteration
	vector<int>::iterator candserialiterator;
	vector<int>::iterator candfromiterator;
	vector<int>::iterator candtoiterator;
	compareBasis = compareBase;
	for (int i = 0; i < 1000; ++i) {
		lineInterDecVec.push_back(0);
		phAngDecVec.push_back(0);
	}
	/*
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
	*/

} // end constructor

Marketover::~Marketover() // destructor
{
	cout << "\nSimulation ended" << endl;
} // destructor ends

void Marketover::MILPMarketover(double LagMultXi[], double LagMultPi[], int totalCandLineNum, int totalSharedNodeNum) // Function MILPMarketover() implements the Mixed Integer Linear Programming Solver routine by calling GLPK routines for average heat rate objective
{
	/* CREATION OF THE MIP SOLVER INSTANCE */
	clock_t begin = clock(); // start the timer
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

	int dimRow = countOfScenarios*(4 * candLineNumber + 2 * sharedELineNumber); // Total number of rows of the A matrix (number of structural constraints of the LP) first term to account for lower and upper line limits & lower and upper definition limits of candidate shared lines, and second term for lower and upper line limits for shared existing lines
	int dimCol = countOfScenarios*nodeNumber+(countOfScenarios+1)*candLineNumber; // Total number of columns of the LP (number of Decision Variables) first term to account for voltage phase angles for inter-zonal lines' nodes, and second term for power flow values and binary integer decision variable values for shared candidate lines
	outPutFile << "\nTotal Number of Structural Constraints (Rows) is: " << dimRow << endl;
	outPutFile << "\nTotal Number of Decision Variables (Columns) is: " << dimCol << endl;
	//cout << "\nTotal Number of Structural Constraints (Rows) is: " << dimRow << endl;
	//cout << "\nTotal Number of Decision Variables (Columns) is: " << dimCol << endl;
	//cout << "\nTotal Number of Shared Nodes is: " << nodeNumber << endl;	
	//cout << "\nTotal Number of shared candidate lines is: " << candLineNumber << endl;
	//cout << "\nTotal Number of shared existing lines is: " << sharedELineNumber << endl;
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
	glp_set_prob_name(milp, "MOTransDec"); // Names the particular problem instance 
	glp_set_obj_dir(milp, GLP_MAX); // Set direction (Declares the MILP Problem as a Maximization Problem)

	/* SPECIFICATION OF PROBLEM PARAMETERS */
	/*Row Definitions: Specification of RHS or b vector of b<=Ax<=b*/
	glp_add_rows(milp, dimRow);
	//Row Definitions and Bounds Corresponding to Constraints/
	string outPGenFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputPowerGLPK/OutPowerGenMOGLPK.txt"; 
	ofstream powerGenOut(outPGenFileName, ios::out);
	if (!powerGenOut){
		cerr << "\nCouldn't open the file" << endl;
		exit(1);
	}
	outPutFile << "\nTesting of interzonal transmission limits" << endl; 
	int rowCount = 1;

	/*******************************************************************************************/

	// Constraints corresponding to Line Forward Flow limits for shared existing lines
	outPutFile << "\nConstraints corresponding to Line Forward Flow limits for shared existing lines" << endl;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) { 
		for (exsharedCapIterator = SECapacity.begin(); exsharedCapIterator != SECapacity.end(); ++exsharedCapIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_UP, 0.0, (*exsharedCapIterator));
			outPutFile << rowCount << "\t";
			outPutFile << ((*exsharedCapIterator)) << endl;
			++rowCount;
		}
	}
	outPutFile << endl;
	// Constraints corresponding to Line Reverse Flow limits for shared existing lines
	outPutFile << "\nConstraints corresponding to Line Reverse Flow limits for shared existing lines" << endl;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (exsharedCapIterator = SECapacity.begin(); exsharedCapIterator != SECapacity.end(); ++exsharedCapIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_LO, -((*exsharedCapIterator)), 0.0);
			outPutFile << rowCount << "\t";
			outPutFile << -((*exsharedCapIterator)) << endl;
			++rowCount;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for Line Flow limits for shared existing lines: " << rowCount << endl;

	/*******************************************************************************************/

	// Constraints corresponding to Line Forward Flow limits for shared candidate lines
	outPutFile << "\nConstraints corresponding to Line Forward Flow limits for shared candidate lines" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candCapIterator = candCapacity.begin(); candCapIterator != candCapacity.end(); ++candCapIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_UP, 0.0, 0.0);
			outPutFile << rowCount << "\t";
			outPutFile << 0.0 << endl;
			++rowCount;
		}
	}
	outPutFile << endl;
	// Constraints corresponding to Line Reverse Flow limits for shared candidate lines
	outPutFile << "\nConstraints corresponding to Line Reverse Flow limits for shared candidate lines" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candCapIterator = candCapacity.begin(); candCapIterator != candCapacity.end(); ++candCapIterator){
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
	outPutFile << "\nTesting of Definition of Flows on Shared Candidate transmission lines" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candCapIterator = candCapacity.begin(); candCapIterator != candCapacity.end(); ++candCapIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_UP, 0.0, (2.5*(*candCapIterator)));
			outPutFile << rowCount << "\t";
			outPutFile << BIGM << endl;
			++rowCount;
		}
	}
	outPutFile << endl;
	// Constraints corresponding to shared candidate lines flow definition lower bound
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candCapIterator = candCapacity.begin(); candCapIterator != candCapacity.end(); ++candCapIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_LO, -(2.5*(*candCapIterator)), 0.0);
			outPutFile << rowCount << "\t";
			outPutFile << -BIGM << endl;
			++rowCount;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for shared candidate lines flow definition: " << rowCount << endl;
	outPutFile << "\nConstraint bounds (rows) Specified" << endl;
	outPutFile << "\nTotal number of rows: " << rowCount - 1 << endl;
	//cout << "\nConstraint bounds (rows) Specified" << endl;
	//cout << "\nTotal number of rows: " << rowCount - 1 << endl;
	/*******************************************************************************************/

	/*******************************************************************************************/

	/*Column Definitions, Bounds, and Objective Function Co-efficients*/
	glp_add_cols(milp, dimCol);
	int colCount = 1;

	/*******************************************************************************************/

	//Columns corresponding to Shared Candidate Line Flows continuous variables//
	outPutFile << "\nColumns corresponding to Shared Candidate Line Flows continuous variables" << endl;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candCapIterator = candCapacity.begin(); candCapIterator != candCapacity.end(); ++candCapIterator){
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
	outPutFile << "\nColumns corresponding to Shared Candidate Line Construction Decision Binary Integer variables" << endl;
	int candInd = 0; // Initialize the counter for indexing the candidate lines
	int diffCandCounter = 0; // counter flag to indicate the first element of the candSerial list
	for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
		if (diffCandCounter > 0) { // Skip the first element, since it's a dummy "0"
			glp_set_col_kind(milp, colCount, GLP_BV);
			glp_set_col_name(milp, colCount, NULL);
			glp_set_col_bnds(milp, colCount, GLP_DB, 0.0, 1.0);
			glp_set_obj_coef(milp, colCount, LagMultPi[candInd]);
			outPutFile << colCount << "\t";
			outPutFile << LagMultPi[candInd] << "\t" << candInd << endl;
			++candInd;
			++colCount;
		}
		++diffCandCounter;
	}
	outPutFile << "\nTotal number of columns after accounting for Shared Candidate Line Construction Decision Binary Integer variables: " << colCount << endl;
	/*******************************************************************************************/

	//Columns corresponding to Voltage Phase Angles continuous variables for other zone nodes for shared lines//
	outPutFile << "\nColumns corresponding to Voltage Phase Angles continuous variables for other zone nodes for shared lines" << endl;
	int otherInd = 0; // Initialize the counter for indexing the other-zone nodes
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int columnCount = 1;
		int diffNodeCounter = 0; // counter flag to indicate the first element of the diffZoneNodeID list
		/*for (diffZNIt = diffZoneNodeID.begin(); diffZNIt != diffZoneNodeID.end(); ++diffZNIt){
			cout << " Column count after the shared node is " << columnCount << endl;	
			++columnCount;
		}*/
		for (diffZNIt = nodeList.begin(); diffZNIt != nodeList.end(); ++diffZNIt){
			if (diffNodeCounter > 0) { // Skip the first element, since it's a dummy "0"
				glp_set_col_kind(milp, colCount, GLP_CV);
				glp_set_col_name(milp, colCount, NULL);
				glp_set_col_bnds(milp, colCount, GLP_DB, 0, (44/7));
				glp_set_obj_coef(milp, colCount, LagMultXi[otherInd]);
				outPutFile << colCount << "\t";
				outPutFile << LagMultXi[otherInd] << "\t" << otherInd << endl;
				++otherInd;	
				++colCount;
			}
			++diffNodeCounter;
		}
	}
	outPutFile << "\nDecision Variables and Objective Function defined" << endl;
	outPutFile << "\nTotal Number of columns: " << colCount - 1 << endl;
	//cout << "\nDecision Variables and Objective Function defined" << endl;
	//cout << "\nTotal Number of columns: " << colCount - 1 << endl;
	/*******************************************************************************************/

	//Non-Zero entries of A matrix (Constraint/Coefficient matrix entries)//
	int index = 1;
	ia.push_back(0), ja.push_back(0), ar.push_back(0.0);
	int rCount = 1; // Initialize the row count

	/*******************************************************************************************/
	
	// Coefficients corresponding to shared existing Line Forward Flow Limit Constraints
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared existing Line Forward Flow Limit Constraints" << endl;
		for (exsharedIterator = SESerial.begin(); exsharedIterator != SESerial.end(); ++exsharedIterator){
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				//cout << "Test Message 2" << endl;
				//cout << "SE From Rank" << SEFromRank.at(jockey) << endl;
				//cout << "SE Reactance" << 1/(SEReactance.at(jockey)) << endl;
				ia.push_back(rCount), ja.push_back((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEFromRank.at(jockey)), ar.push_back(1/(SEReactance.at(jockey)));
				++index;
				outPutFile << "\n" << rCount << "\t" << (countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEFromRank.at(jockey) << "\t" << 1/(SEReactance.at(jockey)) << endl;
				ia.push_back(rCount), ja.push_back((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEToRank.at(jockey)), ar.push_back(-1/(SEReactance.at(jockey)));
				++index;
				outPutFile << "\n" << rCount << "\t" << (countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEToRank.at(jockey) << "\t" << -1/(SEReactance.at(jockey)) << endl;
				++rCount; // Increment the row count to point to the next transmission line object
				++jockey;
				//cout << "\nTest after shared existing Line Forward Flow Limit Constraint " << diffSerCounter << endl;
			}
			++diffSerCounter;	
		}
	}
	// Coefficients corresponding to shared existing Line Reverse Flow Limit Constraints
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared existing Line Reverse Flow Limit Constraints" << endl;
		for (exsharedIterator = SESerial.begin(); exsharedIterator != SESerial.end(); ++exsharedIterator){
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				ia.push_back(rCount), ja.push_back((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEFromRank.at(jockey)), ar.push_back(1/(SEReactance.at(jockey)));
				++index;
				outPutFile << "\n" << rCount << "\t" << (countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEFromRank.at(jockey) << "\t" << 1/(SEReactance.at(jockey)) << endl;
				ia.push_back(rCount), ja.push_back((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEToRank.at(jockey)), ar.push_back(-1/(SEReactance.at(jockey)));
				++index;
				outPutFile << "\n" << rCount << "\t" << (countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEToRank.at(jockey) << "\t" << -1/(SEReactance.at(jockey)) << endl;
				++rCount; // Increment the row count to point to the next transmission line object
				++jockey;
				//cout << "\nTest after shared existing Line Reverse Flow Limit Constraint " << diffSerCounter << endl;
			}
			++diffSerCounter;
		}
	}
	/*******************************************************************************************/
	
	// Coefficients corresponding to shared candidate Line Forward Flow Limit Constraints
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared candidate Line Forward Flow Limit Constraints" << endl;
		for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				ia.push_back(rCount), ja.push_back(rCount-2*countOfScenarios*sharedELineNumber), ar.push_back(1);
				++index;
				outPutFile << "\n" << rCount << "\t" << rCount-2*countOfScenarios*sharedELineNumber << "\t" << 1 << endl;
				ia.push_back(rCount), ja.push_back(countOfScenarios*candLineNumber+rCount-2*countOfScenarios*sharedELineNumber-scenCounter*candLineNumber), ar.push_back(-candCapacity.at(jockey));
				++index;
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*candLineNumber+rCount-2*countOfScenarios*sharedELineNumber-scenCounter*candLineNumber << "\t" << -candCapacity.at(jockey) << endl;
				++rCount;
				++jockey;
				//cout << "\nTest after shared candidate Line Forward Flow Limit Constraint " << diffSerCounter << endl;
			}
			++diffSerCounter;
		}
	}
	// Coefficients corresponding to shared candidate Line Reverse Flow Limit Constraints
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared candidate Line Reverse Flow Limit Constraints" << endl;
		for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber)), ar.push_back(1);
				++index;
				outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber) << "\t" << 1 << endl;
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(2*sharedELineNumber)-scenCounter*candLineNumber), ar.push_back(candCapacity.at(jockey));
				++index;
				outPutFile << "\n" << rCount << "\t" << (rCount-countOfScenarios*(2*sharedELineNumber)-scenCounter*candLineNumber) << "\t" << candCapacity.at(jockey) << endl;
				++rCount;
				++jockey;
				//cout << "\nTest after shared candidate Line Reverse Flow Limit Constraint " << diffSerCounter << endl;
			}
			++diffSerCounter;
		}
	}
	/*******************************************************************************************/
	
	// Coefficients corresponding to shared candidate Line Definition upper bound
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared candidate Line Definition upper bound" << endl;
		for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(2*sharedELineNumber+2*candLineNumber)), ar.push_back(1);
				++index;
				outPutFile << rCount << "\t" << rCount-countOfScenarios*(2*sharedELineNumber+2*candLineNumber) << " Power term for candidate line " << diffSerCounter << "-th candidate line" << endl;
				ia.push_back(rCount), ja.push_back((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candFromRank.at(jockey)), ar.push_back(-1/(candReactance.at(jockey)));
				++index;
				outPutFile << rCount << "\t" << (countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candFromRank.at(jockey) << "\t" << -1/(candReactance.at(jockey)) << endl;
				ia.push_back(rCount), ja.push_back((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candToRank.at(jockey)), ar.push_back(1/(candReactance.at(jockey)));
				++index;
				outPutFile << rCount << "\t" << (countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candToRank.at(jockey) << "\t" << 1/(candReactance.at(jockey)) << endl;
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber)-scenCounter*candLineNumber), ar.push_back(2.5*(candCapacity.at(jockey)));
				++index;
				outPutFile << rCount << "\t" << rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber)-scenCounter*candLineNumber << "\t" << BIGM << endl;
				++rCount; // Increment the row count to point to the next transmission line object
				++jockey;
				//cout << "\nTest after shared candidate Line Definition upper bound " << diffSerCounter << endl;
			}
			++diffSerCounter;
		}
	}
	// Coefficients corresponding to shared candidate Line Definition lower bound
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared candidate Line Definition lower bound" << endl;
		for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(2*sharedELineNumber+3*candLineNumber)), ar.push_back(1);
				++index;
				outPutFile << rCount << "\t" << rCount-countOfScenarios*(2*sharedELineNumber+3*candLineNumber) << " Power term for candidate line " << diffSerCounter << "-th candidate line" << endl;
				ia.push_back(rCount), ja.push_back((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candFromRank.at(jockey)), ar.push_back(-1/(candReactance.at(jockey)));
				++index;
				outPutFile << rCount << "\t" << (countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candFromRank.at(jockey) << "\t" << -1/(candReactance.at(jockey)) << endl;
				ia.push_back(rCount), ja.push_back((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candToRank.at(jockey)), ar.push_back(1/(candReactance.at(jockey)));
				++index;
				outPutFile << rCount << "\t" << (countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candToRank.at(jockey) << "\t" << 1/(candReactance.at(jockey)) << endl;
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(2*sharedELineNumber+2*candLineNumber)-scenCounter*candLineNumber), ar.push_back(-(2.5*(candCapacity.at(jockey))));
				++index;
				outPutFile << rCount << "\t" << rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber)-scenCounter*candLineNumber << "\t" << -BIGM << endl;
				++rCount; // Increment the row count to point to the next transmission line object
				++jockey;
				//cout << "\nTest after shared candidate Line Definition lower bound " << diffSerCounter << endl;
			}
			++diffSerCounter;
		}
	}
	/*******************************************************************************************/

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
	string outLPLogFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputSummaryGLPK/OutLPLogMOGLPK.txt";
	string outMILPLogFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputSummaryGLPK/OutMIPLogMOGLPK.txt";
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

	// Open separate output files for writing results of different variables
	string outCandFlowFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/candidateLinesGLPK/candFlowMWMOGLPK.txt";
	string outCandDecFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/candidateLinesGLPK/candLineDecisionMOGLPK.txt";
	string outExtAngFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputAnglesGLPK/externalAngleMOGLPK.txt";
	ofstream candFlowMWOut(outCandFlowFileName, ios::out); //switchOnOut
	ofstream candLineDecisionOut(outCandDecFileName, ios::out); //switchOffOut
	ofstream externalAngleOut(outExtAngFileName, ios::out);
	outPutFile << "\nThe Optimal Objective value (Line Building Decision cost) is: " << -z << endl;
	powerGenOut << "\nThe Optimal Objective value (Line Building Decision cost) is: " << -z << endl;
	cout << "\nThe Optimal Objective value (Line Building Decision cost) is: " << -z << endl;
	x.push_back(0); // Initialize the decision Variable vector

	// Display Shared Candidate lines' Power Flow variables
	candFlowMWOut << "\n****************** SHARED CANDIDATE LINES POWER FLOW VALUES *********************" << endl;
	candFlowMWOut << "CANDIDATE LINE ID" << "\t" << "MW FLOW VALUE" << "\n";
	int arrayInd = 1;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int diffCandCounter2 = 0; // counter flag to indicate the first element of the candSerial list
		for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
			if (diffCandCounter2 > 0) { // Skip the first element, since it's a dummy "0"
				x.push_back(glp_mip_col_val(milp, arrayInd));
				candFlowMWOut << (*candIterator) << "\t" << (glp_mip_col_val(milp, arrayInd))*100 << " MW" << endl;
				++arrayInd;
			}
			++diffCandCounter2;
		}
	}
	candFlowMWOut << "Finished writing Shared Candidate lines' Power Flow variables" << endl;

	// Display Shared Candidate lines' Construction Decisions
	candLineDecisionOut << "\n****************** SHARED CANDIDATE LINES CONSTRUCTION DECISIONS *********************" << endl;
	candLineDecisionOut << "CANDIDATE LINE ID" << "\t" << "CONSTRUCTION DECISIONS" << "\n";
	int candInd1 = 0; // Initialize the counter for indexing the candidate lines
	int diffCandCounter1 = 0; // counter flag to indicate the first element of the candSerial list
	for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
		if (diffCandCounter1 > 0) { // Skip the first element, since it's a dummy "0"
			x.push_back(glp_mip_col_val(milp, arrayInd));
			interimIntDecVar.push_back(glp_mip_col_val(milp, arrayInd));
			candLineDecisionOut << (*candIterator) << "\t" << (glp_mip_col_val(milp, arrayInd)) << endl;
			++arrayInd;
			++candInd1;
		}
		++diffCandCounter1;
	}
	candLineDecisionOut << "Finished writing Shared Candidate lines' Construction decisions" << endl;

	// Display shared node angles
	externalAngleOut << "\n****************** OUTER ZONAL NODE VOLTAGE PHASE ANGLE VALUES *********************" << endl;
	externalAngleOut << "EXTERNAL NODE ID" << "\t" << "VOLTAGE PHASE ANGLE" << "\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int otherInd1 = 0; // Initialize the counter for indexing the other-zone nodes
		int dummyStart = 0;
		for (diffZNIt = nodeList.begin(); diffZNIt != nodeList.end(); ++diffZNIt){
			if (dummyStart > 0) { // Skip the dummy element "0" at the beginning
				x.push_back(glp_mip_col_val(milp, arrayInd));
				interimContDecVar.push_back(glp_mip_col_val(milp, arrayInd));
				externalAngleOut << (*diffZNIt) << "\t" << (glp_mip_col_val(milp, arrayInd)) << endl;
				++arrayInd;
				++otherInd1;
			}
			++dummyStart;
		}
	}
	externalAngleOut << "Finished writing shared node voltage phase angle values" << endl;

	delete ipControlParam; // free the memory of the Integer Programming Control Parameter struct
	glp_delete_prob(milp); // Free the memory of the GLPK Problem Object
	clock_t end2 = clock(); // stop the timer
	double elapsed_secs2 = double(end2 - begin) / CLOCKS_PER_SEC; // Calculate the Total Time
	outPutFile << "\nTotal time taken to solve the MILP Line Construction Decision Making Problem instance and retrieve the results = " << elapsed_secs2 << " s " << endl;
	cout << "\nTotal time taken to solve the MILP Line Construction Decision Making Problem instance and retrieve the results = " << elapsed_secs2 << " s " << endl;

	// Close the different output files
	outPutFile.close();
	powerGenOut.close();
	candFlowMWOut.close();
	candLineDecisionOut.close();
	externalAngleOut.close();
	cout << "\nSimulation Completed.\nResults written on the different output files" << endl;
	//%%calcMILPBounds(); // Calculate the bounds
} // Function MILP() ends

double Marketover::LBMarketover(double LagMultXi[], double LagMultPi[], int totalCandLineNum, int totalSharedNodeNum) // Function LBMarketover() calculates the lower bound of the Mixed Integer Linear Programming Solver routine by calling GLPK routines for average heat rate objective
{
	/* CREATION OF THE MIP SOLVER INSTANCE */
	clock_t begin = clock(); // start the timer
	vector<int>::iterator diffZNIt; // Iterator for diffZoneNodeID
	vector<Node*>::iterator nodeIterator; // Iterator for node objects
	vector<double>::iterator candCapIterator; // Iterator for candidate lines for capacity
	vector<int>::iterator candIterator; // Iterator for candidate lines	
	vector<double>::iterator exsharedCapIterator; // Iterator for shared existing lines for capacity
	vector<int>::iterator exsharedIterator;

	string outSummaryFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputSummaryGLPKBounds/OutSummaryMOBoundGLPK.txt";
	ofstream outPutFile(outSummaryFileName, ios::out); // Create Output File to output the Summary of Results
	if (!outPutFile){
		cerr << "\nCouldn't open the file" << endl;
		exit(1);
	}

	int dimRow = countOfScenarios*(4 * candLineNumber + 2 * sharedELineNumber); // Total number of rows of the A matrix (number of structural constraints of the LP) first term to account for lower and upper line limits & lower and upper definition limits of candidate shared lines, and second term for lower and upper line limits for shared existing lines
	int dimCol = countOfScenarios*nodeNumber+(countOfScenarios+1)*candLineNumber; // Total number of columns of the LP (number of Decision Variables) first term to account for voltage phase angles for inter-zonal lines' nodes, and second term for power flow values and binary integer decision variable values for shared candidate lines
	outPutFile << "\nTotal Number of Structural Constraints (Rows) is: " << dimRow << endl;
	outPutFile << "\nTotal Number of Decision Variables (Columns) is: " << dimCol << endl;
	//cout << "\nTotal Number of Structural Constraints (Rows) is: " << dimRow << endl;
	//cout << "\nTotal Number of Decision Variables (Columns) is: " << dimCol << endl;
	//cout << "\nTotal Number of Shared Nodes is: " << nodeNumber << endl;	
	//cout << "\nTotal Number of shared candidate lines is: " << candLineNumber << endl;
	//cout << "\nTotal Number of shared existing lines is: " << sharedELineNumber << endl;
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
	glp_set_prob_name(milp, "MOTransDec"); // Names the particular problem instance 
	glp_set_obj_dir(milp, GLP_MAX); // Set direction (Declares the MILP Problem as a Minimization Problem)

	/* SPECIFICATION OF PROBLEM PARAMETERS */
	/*Row Definitions: Specification of RHS or b vector of b<=Ax<=b*/
	glp_add_rows(milp, dimRow);
	//Row Definitions and Bounds Corresponding to Constraints/
	string outPGenFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputPowerGLPKBounds/OutPowerGenMOBoundGLPK.txt"; 
	ofstream powerGenOut(outPGenFileName, ios::out);
	if (!powerGenOut){
		cerr << "\nCouldn't open the file" << endl;
		exit(1);
	}
	int rowCount = 1;

	/*******************************************************************************************/

	// Constraints corresponding to Line Forward Flow limits for shared existing lines
	outPutFile << "\nConstraints corresponding to Line Forward Flow limits for shared existing lines" << endl;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) { 
		for (exsharedCapIterator = SECapacity.begin(); exsharedCapIterator != SECapacity.end(); ++exsharedCapIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_UP, 0.0, (*exsharedCapIterator));
			outPutFile << rowCount << "\t";
			outPutFile << ((*exsharedCapIterator)) << endl;
			++rowCount;
		}
	}
	outPutFile << endl;
	// Constraints corresponding to Line Reverse Flow limits for shared existing lines
	outPutFile << "\nConstraints corresponding to Line Reverse Flow limits for shared existing lines" << endl;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (exsharedCapIterator = SECapacity.begin(); exsharedCapIterator != SECapacity.end(); ++exsharedCapIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_LO, -((*exsharedCapIterator)), 0.0);
			outPutFile << rowCount << "\t";
			outPutFile << -((*exsharedCapIterator)) << endl;
			++rowCount;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for Line Flow limits for shared existing lines: " << rowCount << endl;

	/*******************************************************************************************/

	// Constraints corresponding to Line Forward Flow limits for shared candidate lines
	outPutFile << "\nConstraints corresponding to Line Forward Flow limits for shared candidate lines" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candCapIterator = candCapacity.begin(); candCapIterator != candCapacity.end(); ++candCapIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_UP, 0.0, 0.0);
			outPutFile << rowCount << "\t";
			outPutFile << 0.0 << endl;
			++rowCount;
		}
	}
	outPutFile << endl;
	// Constraints corresponding to Line Reverse Flow limits for shared candidate lines
	outPutFile << "\nConstraints corresponding to Line Reverse Flow limits for shared candidate lines" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candCapIterator = candCapacity.begin(); candCapIterator != candCapacity.end(); ++candCapIterator){
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
	outPutFile << "\nTesting of Definition of Flows on Shared Candidate transmission lines" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candCapIterator = candCapacity.begin(); candCapIterator != candCapacity.end(); ++candCapIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_UP, 0.0, (2.5*(*candCapIterator)));
			outPutFile << rowCount << "\t";
			outPutFile << BIGM << endl;
			++rowCount;
		}
	}
	outPutFile << endl;
	// Constraints corresponding to shared candidate lines flow definition lower bound
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candCapIterator = candCapacity.begin(); candCapIterator != candCapacity.end(); ++candCapIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_LO, -(2.5*(*candCapIterator)), 0.0);
			outPutFile << rowCount << "\t";
			outPutFile << -BIGM << endl;
			++rowCount;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for shared candidate lines flow definition: " << rowCount << endl;
	outPutFile << "\nConstraint bounds (rows) Specified" << endl;
	outPutFile << "\nTotal number of rows: " << rowCount - 1 << endl;
	//cout << "\nConstraint bounds (rows) Specified" << endl;
	//cout << "\nTotal number of rows: " << rowCount - 1 << endl;
	/*******************************************************************************************/

	/*******************************************************************************************/

	/*Column Definitions, Bounds, and Objective Function Co-efficients*/
	glp_add_cols(milp, dimCol);
	int colCount = 1;

	/*******************************************************************************************/

	//Columns corresponding to Shared Candidate Line Flows continuous variables//
	outPutFile << "\nColumns corresponding to Shared Candidate Line Flows continuous variables" << endl;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candCapIterator = candCapacity.begin(); candCapIterator != candCapacity.end(); ++candCapIterator){
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
	outPutFile << "\nColumns corresponding to Shared Candidate Line Construction Decision Binary Integer variables" << endl;
	int candInd = 0; // Initialize the counter for indexing the candidate lines
	int diffCandCounter = 0; // counter flag to indicate the first element of the candSerial list
	for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
		if (diffCandCounter > 0) { // Skip the first element, since it's a dummy "0"
			glp_set_col_kind(milp, colCount, GLP_CV);
			glp_set_col_name(milp, colCount, NULL);
			glp_set_col_bnds(milp, colCount, GLP_DB, 0.0, 1.0);
			glp_set_obj_coef(milp, colCount, LagMultPi[candInd]);
			outPutFile << colCount << "\t";
			outPutFile << LagMultPi[candInd] << "\t" << candInd << endl;
			++candInd;
			++colCount;
		}
		++diffCandCounter;
	}
	outPutFile << "\nTotal number of columns after accounting for Shared Candidate Line Construction Decision Binary Integer variables: " << colCount << endl;
	/*******************************************************************************************/

	//Columns corresponding to Voltage Phase Angles continuous variables for other zone nodes for shared lines//
	outPutFile << "\nColumns corresponding to Voltage Phase Angles continuous variables for other zone nodes for shared lines" << endl;
	int otherInd = 0; // Initialize the counter for indexing the other-zone nodes
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int columnCount = 1;
		int diffNodeCounter = 0; // counter flag to indicate the first element of the diffZoneNodeID list
		/*for (diffZNIt = diffZoneNodeID.begin(); diffZNIt != diffZoneNodeID.end(); ++diffZNIt){
			cout << " Column count after the shared node is " << columnCount << endl;	
			++columnCount;
		}*/
		for (diffZNIt = nodeList.begin(); diffZNIt != nodeList.end(); ++diffZNIt){
			if (diffNodeCounter > 0) { // Skip the first element, since it's a dummy "0"
				glp_set_col_kind(milp, colCount, GLP_CV);
				glp_set_col_name(milp, colCount, NULL);
				glp_set_col_bnds(milp, colCount, GLP_DB, 0, (44/7));
				glp_set_obj_coef(milp, colCount, LagMultXi[otherInd]);
				outPutFile << colCount << "\t";
				outPutFile << LagMultXi[otherInd] << "\t" << otherInd << endl;
				++otherInd;	
				++colCount;
			}
			++diffNodeCounter;
		}
	}
	outPutFile << "\nDecision Variables and Objective Function defined" << endl;
	outPutFile << "\nTotal Number of columns: " << colCount - 1 << endl;
	//cout << "\nDecision Variables and Objective Function defined" << endl;
	//cout << "\nTotal Number of columns: " << colCount - 1 << endl;
	/*******************************************************************************************/

	//Non-Zero entries of A matrix (Constraint/Coefficient matrix entries)//
	int index = 1;
	ia.push_back(0), ja.push_back(0), ar.push_back(0.0);
	int rCount = 1; // Initialize the row count

	/*******************************************************************************************/
	
	// Coefficients corresponding to shared existing Line Forward Flow Limit Constraints
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared existing Line Forward Flow Limit Constraints" << endl;
		for (exsharedIterator = SESerial.begin(); exsharedIterator != SESerial.end(); ++exsharedIterator){
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				//cout << "Test Message 2" << endl;
				//cout << "SE From Rank" << SEFromRank.at(jockey) << endl;
				//cout << "SE Reactance" << 1/(SEReactance.at(jockey)) << endl;
				ia.push_back(rCount), ja.push_back((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEFromRank.at(jockey)), ar.push_back(1/(SEReactance.at(jockey)));
				++index;
				outPutFile << "\n" << rCount << "\t" << (countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEFromRank.at(jockey) << "\t" << 1/(SEReactance.at(jockey)) << endl;
				ia.push_back(rCount), ja.push_back((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEToRank.at(jockey)), ar.push_back(-1/(SEReactance.at(jockey)));
				++index;
				outPutFile << "\n" << rCount << "\t" << (countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEToRank.at(jockey) << "\t" << -1/(SEReactance.at(jockey)) << endl;
				++rCount; // Increment the row count to point to the next transmission line object
				++jockey;
				//cout << "\nTest after shared existing Line Forward Flow Limit Constraint " << diffSerCounter << endl;
			}
			++diffSerCounter;	
		}
	}
	// Coefficients corresponding to shared existing Line Reverse Flow Limit Constraints
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared existing Line Reverse Flow Limit Constraints" << endl;
		for (exsharedIterator = SESerial.begin(); exsharedIterator != SESerial.end(); ++exsharedIterator){
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				ia.push_back(rCount), ja.push_back((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEFromRank.at(jockey)), ar.push_back(1/(SEReactance.at(jockey)));
				++index;
				outPutFile << "\n" << rCount << "\t" << (countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEFromRank.at(jockey) << "\t" << 1/(SEReactance.at(jockey)) << endl;
				ia.push_back(rCount), ja.push_back((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEToRank.at(jockey)), ar.push_back(-1/(SEReactance.at(jockey)));
				++index;
				outPutFile << "\n" << rCount << "\t" << (countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEToRank.at(jockey) << "\t" << -1/(SEReactance.at(jockey)) << endl;
				++rCount; // Increment the row count to point to the next transmission line object
				++jockey;
				//cout << "\nTest after shared existing Line Reverse Flow Limit Constraint " << diffSerCounter << endl;
			}
			++diffSerCounter;
		}
	}
	/*******************************************************************************************/
	
	// Coefficients corresponding to shared candidate Line Forward Flow Limit Constraints
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared candidate Line Forward Flow Limit Constraints" << endl;
		for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				ia.push_back(rCount), ja.push_back(rCount-2*countOfScenarios*sharedELineNumber), ar.push_back(1);
				++index;
				outPutFile << "\n" << rCount << "\t" << rCount-2*countOfScenarios*sharedELineNumber << "\t" << 1 << endl;
				ia.push_back(rCount), ja.push_back(countOfScenarios*candLineNumber+rCount-2*countOfScenarios*sharedELineNumber-scenCounter*candLineNumber), ar.push_back(-candCapacity.at(jockey));
				++index;
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*candLineNumber+rCount-2*countOfScenarios*sharedELineNumber-scenCounter*candLineNumber << "\t" << -candCapacity.at(jockey) << endl;
				++rCount;
				++jockey;
				//cout << "\nTest after shared candidate Line Forward Flow Limit Constraint " << diffSerCounter << endl;
			}
			++diffSerCounter;
		}
	}
	// Coefficients corresponding to shared candidate Line Reverse Flow Limit Constraints
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared candidate Line Reverse Flow Limit Constraints" << endl;
		for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber)), ar.push_back(1);
				++index;
				outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber) << "\t" << 1 << endl;
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(2*sharedELineNumber)-scenCounter*candLineNumber), ar.push_back(candCapacity.at(jockey));
				++index;
				outPutFile << "\n" << rCount << "\t" << (rCount-countOfScenarios*(2*sharedELineNumber)-scenCounter*candLineNumber) << "\t" << candCapacity.at(jockey) << endl;
				++rCount;
				++jockey;
				//cout << "\nTest after shared candidate Line Reverse Flow Limit Constraint " << diffSerCounter << endl;
			}
			++diffSerCounter;
		}
	}
	/*******************************************************************************************/
	
	// Coefficients corresponding to shared candidate Line Definition upper bound
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared candidate Line Definition upper bound" << endl;
		for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(2*sharedELineNumber+2*candLineNumber)), ar.push_back(1);
				++index;
				outPutFile << rCount << "\t" << rCount-countOfScenarios*(2*sharedELineNumber+2*candLineNumber) << " Power term for candidate line " << diffSerCounter << "-th candidate line" << endl;
				ia.push_back(rCount), ja.push_back((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candFromRank.at(jockey)), ar.push_back(-1/(candReactance.at(jockey)));
				++index;
				outPutFile << rCount << "\t" << (countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candFromRank.at(jockey) << "\t" << -1/(candReactance.at(jockey)) << endl;
				ia.push_back(rCount), ja.push_back((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candToRank.at(jockey)), ar.push_back(1/(candReactance.at(jockey)));
				++index;
				outPutFile << rCount << "\t" << (countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candToRank.at(jockey) << "\t" << 1/(candReactance.at(jockey)) << endl;
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber)-scenCounter*candLineNumber), ar.push_back(2.5*(candCapacity.at(jockey)));
				++index;
				outPutFile << rCount << "\t" << rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber)-scenCounter*candLineNumber << "\t" << BIGM << endl;
				++rCount; // Increment the row count to point to the next transmission line object
				++jockey;
				//cout << "\nTest after shared candidate Line Definition upper bound " << diffSerCounter << endl;
			}
			++diffSerCounter;
		}
	}
	// Coefficients corresponding to shared candidate Line Definition lower bound
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared candidate Line Definition lower bound" << endl;
		for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(2*sharedELineNumber+3*candLineNumber)), ar.push_back(1);
				++index;
				outPutFile << rCount << "\t" << rCount-countOfScenarios*(2*sharedELineNumber+3*candLineNumber) << " Power term for candidate line " << diffSerCounter << "-th candidate line" << endl;
				ia.push_back(rCount), ja.push_back((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candFromRank.at(jockey)), ar.push_back(-1/(candReactance.at(jockey)));
				++index;
				outPutFile << rCount << "\t" << (countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candFromRank.at(jockey) << "\t" << -1/(candReactance.at(jockey)) << endl;
				ia.push_back(rCount), ja.push_back((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candToRank.at(jockey)), ar.push_back(1/(candReactance.at(jockey)));
				++index;
				outPutFile << rCount << "\t" << (countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candToRank.at(jockey) << "\t" << 1/(candReactance.at(jockey)) << endl;
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(2*sharedELineNumber+2*candLineNumber)-scenCounter*candLineNumber), ar.push_back(-(2.5*(candCapacity.at(jockey))));
				++index;
				outPutFile << rCount << "\t" << rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber)-scenCounter*candLineNumber << "\t" << -BIGM << endl;
				++rCount; // Increment the row count to point to the next transmission line object
				++jockey;
				//cout << "\nTest after shared candidate Line Definition lower bound " << diffSerCounter << endl;
			}
			++diffSerCounter;
		}
	}
	/*******************************************************************************************/

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
	string outLPLogFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputSummaryGLPKBounds/OutLPLogMOBoundGLPK.txt";
	string outMILPLogFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputSummaryGLPKBounds/OutMIPLogMOBoundGLPK.txt";
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

	// Open separate output files for writing results of different variables
	string outCandFlowFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/candidateLinesGLPKBounds/candFlowMWMOBoundGLPK.txt";
	string outCandDecFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/candidateLinesGLPKBounds/candLineDecisionMOBoundGLPK.txt";
	string outExtAngFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputAnglesGLPKBounds/externalAngleMOBoundGLPK.txt";
	ofstream candFlowMWOut(outCandFlowFileName, ios::out); //switchOnOut
	ofstream candLineDecisionOut(outCandDecFileName, ios::out); //switchOffOut
	ofstream externalAngleOut(outExtAngFileName, ios::out);
	outPutFile << "\nThe Optimal Objective value (Line Building Decision cost) is: " << -z << endl;
	powerGenOut << "\nThe Optimal Objective value (Line Building Decision cost) is: " << -z << endl;
	cout << "\nThe Optimal Objective value (Line Building Decision cost) is: " << -z << endl;
	x.push_back(0); // Initialize the decision Variable vector

	// Display Shared Candidate lines' Power Flow variables
	candFlowMWOut << "\n****************** SHARED CANDIDATE LINES POWER FLOW VALUES *********************" << endl;
	candFlowMWOut << "CANDIDATE LINE ID" << "\t" << "MW FLOW VALUE" << "\n";
	int arrayInd = 1;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int diffCandCounter2 = 0; // counter flag to indicate the first element of the candSerial list
		for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
			if (diffCandCounter2 > 0) { // Skip the first element, since it's a dummy "0"
				x.push_back(glp_mip_col_val(milp, arrayInd));
				candFlowMWOut << (*candIterator) << "\t" << (glp_mip_col_val(milp, arrayInd))*100 << " MW" << endl;
				++arrayInd;
			}
			++diffCandCounter2;
		}
	}
	candFlowMWOut << "Finished writing Shared Candidate lines' Power Flow variables" << endl;

	// Display Shared Candidate lines' Construction Decisions
	candLineDecisionOut << "\n****************** SHARED CANDIDATE LINES CONSTRUCTION DECISIONS *********************" << endl;
	candLineDecisionOut << "CANDIDATE LINE ID" << "\t" << "CONSTRUCTION DECISIONS" << "\n";
	int candInd1 = 0; // Initialize the counter for indexing the candidate lines
	int diffCandCounter1 = 0; // counter flag to indicate the first element of the candSerial list
	for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
		if (diffCandCounter1 > 0) { // Skip the first element, since it's a dummy "0"
			x.push_back(glp_mip_col_val(milp, arrayInd));
			candLineDecisionOut << (*candIterator) << "\t" << (glp_mip_col_val(milp, arrayInd)) << endl;
			++arrayInd;
			++candInd1;
		}
		++diffCandCounter1;
	}
	candLineDecisionOut << "Finished writing Shared Candidate lines' Construction decisions" << endl;

	// Display shared node angles
	externalAngleOut << "\n****************** OUTER ZONAL NODE VOLTAGE PHASE ANGLE VALUES *********************" << endl;
	externalAngleOut << "EXTERNAL NODE ID" << "\t" << "VOLTAGE PHASE ANGLE" << "\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int otherInd1 = 0; // Initialize the counter for indexing the other-zone nodes
		int dummyStart = 0;
		for (diffZNIt = nodeList.begin(); diffZNIt != nodeList.end(); ++diffZNIt){
			if (dummyStart > 0) { // Skip the dummy element "0" at the beginning
				x.push_back(glp_mip_col_val(milp, arrayInd));
				externalAngleOut << (*diffZNIt) << "\t" << (glp_mip_col_val(milp, arrayInd)) << endl;
				++arrayInd;
				++otherInd1;
			}
			++dummyStart;
		}
	}
	externalAngleOut << "Finished writing shared node voltage phase angle values" << endl;

	delete ipControlParam; // free the memory of the Integer Programming Control Parameter struct
	glp_delete_prob(milp); // Free the memory of the GLPK Problem Object
	clock_t end2 = clock(); // stop the timer
	double elapsed_secs2 = double(end2 - begin) / CLOCKS_PER_SEC; // Calculate the Total Time
	outPutFile << "\nTotal time taken to solve the MILP Line Construction Decision Making Problem instance and retrieve the results = " << elapsed_secs2 << " s " << endl;
	cout << "\nTotal time taken to solve the MILP Line Construction Decision Making Problem instance and retrieve the results = " << elapsed_secs2 << " s " << endl;

	// Close the different output files
	outPutFile.close();
	powerGenOut.close();
	candFlowMWOut.close();
	candLineDecisionOut.close();
	externalAngleOut.close();
	cout << "\nSimulation Completed.\nResults written on the different output files" << endl;
	//%%calcMILPBounds(); // Calculate the bounds
	return -z;
} // Function MILP() ends

void Marketover::bufferintermediateDecision(int iterCountOut)
{
	if (iterCountOut==0) {
		for (int scenarioTrack = 0; scenarioTrack < countOfScenarios; ++scenarioTrack) {
			for (int nodeTrack = 0; nodeTrack < nodeNumber; ++nodeTrack) {
				interimContDecVarPrev.push_back(0);
			}
		}
		for (int lineTrack = 0; lineTrack < candLineNumber; ++lineTrack) {
			interimIntDecVarPrev.push_back(0);
		}		
	}
	else {
		interimContDecVarPrev = interimContDecVar;
		interimIntDecVarPrev = interimIntDecVar;
	}
}

double Marketover::getGlobalUpper(double LagMultXi[], double LagMultPi[], double regionalUpper[], int numberOfZones) // Function getGlobalUpper() returns the global upper bound for the investment coordination problem
{
	// Sum up the regional upper bounds, which are the tentative regional minima, at the end of every iteration
	double revisedUpperBound = 0; // total upper bound initialized
	for ( int i = 0; i < numberOfZones; ++i ) {
		revisedUpperBound += regionalUpper[i];
	}
	// Calculate the MO upper bound, which is the minimum value of the MO objective, calculated with a given mix of the intermediate, shared continuous variable values from different regions/zones
	for (int scenPos = 0; scenPos < countOfScenarios; ++scenPos) {
		int interNodeRank; // Intermediate variable for storing the rank of the node 
		vector<int> rankPresChecker; // Vector to check if a particular rank has been accounted for 
		int length = (angleDecIndex[scenPos]).size(); // length of the rank vector
		vector<int>::iterator angleDecIndexIterator; // Iterator for the rank vector
		for (int i = 0; i < length; ++i) // Put as many zeroes on the rankPresChecker vector as is the length of the rank vector 
			rankPresChecker.push_back(0);
		vector<double> interAngleTermVec; // Intermediate vector for storing the costs of the angle terms from each sub-network
		int tracker = 0; // tracker to track the position of the angleDecIndexIterator
		for (angleDecIndexIterator=(angleDecIndex[scenPos]).begin(); angleDecIndexIterator!=(angleDecIndex[scenPos]).end(); ++angleDecIndexIterator) { // Iterate through rank vector
			interNodeRank = (*angleDecIndexIterator); // Store the value of the rank of the present node iterate in the list 
			if ( rankPresChecker.at(tracker) == 0 ) { // If this node rank hasn't been already accounted for
				auto pos = std::find((angleDecIndex[scenPos]).begin(), (angleDecIndex[scenPos]).end(), interNodeRank); // find the first position of this rank in the vector
				while(pos != (angleDecIndex[scenPos]).end()) // while all the different positions of this rank hasn't been accounted for 
    				{
      					auto pos1 = std::distance((angleDecIndex[scenPos]).begin(), pos); // get the location in the vector, of this rank
					rankPresChecker.at(pos1) = 1; // Mark the element in the corresponding position of checker vector 1, to indicate this rank has been accounted for 
					interAngleTermVec.push_back(-((phaseAngleDecision[scenPos]).at(pos1))*(LagMultXi[scenPos*nodeNumber+(*angleDecIndexIterator-1)])); // Calculate cost term
      					pos = std::find(pos + 1, (angleDecIndex[scenPos]).end(), interNodeRank); // Find position of the next occurence of this rank
   				}
				double smallest_element = *min_element(interAngleTermVec.begin(), interAngleTermVec.end());
				revisedUpperBound += smallest_element;
				interAngleTermVec.clear();					
			}
			++tracker; // Increment the tracker
		}
	}
	// Calculate the MO upper bound, which is the minimum value of the MO objective, calculated with a given mix of the intermediate, shared discrete variable values from different regions/zones
	int interLineRank; // Intermediate variable for storing the rank of the line 
	vector<int> rankPresCheckerInt; // Vector to check if a particular rank has been accounted for 
	int lengthInt = lineDecIndex.size(); // length of the rank vector
	vector<int>::iterator lineDecIndexIterator; // Iterator for the rank vector
	for (int i = 0; i < lengthInt; ++i) // Put as many zeroes on the rankPresChecker vector as is the length of the rank vector 
		rankPresCheckerInt.push_back(0);
	vector<double> interLineTermVec; // Intermediate vector for storing the costs of the line building decisions from each sub-network
	int trackerInt = 0; // tracker to track the position of the lineDecIndexIterator
	for (lineDecIndexIterator=lineDecIndex.begin(); lineDecIndexIterator!=lineDecIndex.end(); ++lineDecIndexIterator) { // Iterate through rank vector
		interLineRank = (*lineDecIndexIterator); // Store the value of the rank of the present line iterate in the list 
		if ( rankPresCheckerInt.at(trackerInt) == 0 ) { // If this line rank hasn't been already accounted for
			auto pos = std::find(lineDecIndex.begin(), lineDecIndex.end(), interLineRank); // find the first position of this rank in the vector
			while(pos != lineDecIndex.end()) // while all the different positions of this rank hasn't been accounted for 
    			{
      				auto pos1 = std::distance(lineDecIndex.begin(), pos); // get the location in the vector, of this rank
				rankPresCheckerInt.at(pos1) = 1; // Mark the element in the corresponding position of checker vector 1, to indicate this rank has been accounted for 
				interLineTermVec.push_back(-(lineInterDecision[pos1])*(LagMultPi[(*lineDecIndexIterator-1)])); // Calculate cost term
      				pos = std::find(pos + 1, lineDecIndex.end(), interLineRank); // Find position of the next occurence of this rank
   			}
			double smallest_element = *min_element(interLineTermVec.begin(), interLineTermVec.end());
			revisedUpperBound += smallest_element;
			interLineTermVec.clear();					
		}
		++trackerInt; // Increment the tracker
	}
	return revisedUpperBound;
	
} // Function getGlobalUpper ends

double Marketover::getGlobalLower(double regionalLower[], int numberOfZones) // Function getGlobalLower() returns the global lower bound for the investment coordination problem
{
	double revisedLowerBound = 0; // total lower bound initialized
	for ( int i = 0; i < numberOfZones; ++i ) {
		revisedLowerBound += regionalLower[i];
	}
	return revisedLowerBound;
	
} // Function getGlobalLower ends

double Marketover::getGlobalConsensus() // Function getGlobalConsensus() returns the global consensus for the investment coordination problem
{
	double globConsensus = 0; // total consensus
	// Calculate the MO upper bound, which is the minimum value of the MO objective, calculated with a given mix of the intermediate, shared continuous variable values from different regions/zones
	for (int scenPos = 0; scenPos < countOfScenarios; ++scenPos) {
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
			++tracker; // Increment the tracker
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

void Marketover::finDecLineConstr() // Final decisions on the construction stauses of candidate lines taking into account the decisions of different zones
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

int Marketover::scanBuiltLinesList(int rankGlob) // Scans the list of built candidate lines for the input argument of global rank to check if this line is built
{
	auto pos = std::find(constructedRanks.begin(), constructedRanks.end(), rankGlob); // find the first position of this rank in the vector
	if(pos != constructedRanks.end()) // if the rank is found in the list
    		return 1; // return 1
	else 
		return 0; // Otherwise, return 0
}

void Marketover::populateLineDec(int interimDecBin, int zonalIndex, int globalRank) // Method to pass the intermediate message for the zonal line building decision to the MO
{
	lineInterDecision.push_back(interimDecBin);
	lineDecIndex.push_back(globalRank);
	zonalIndVectorInt.push_back(zonalIndex);
	lineInterDecVec.at(zonalIndex*candLineNumber+globalRank)=interimDecBin;
}

void Marketover::populateAngleDec(double interimDecAngle, int zonalIndex, int scenario, int globalRank) // Method to pass the intermediate message for the zonal node angle decision to the MO
{
	(phaseAngleDecision[scenario]).push_back(interimDecAngle);
	(angleDecIndex[scenario]).push_back(globalRank);
	(zonalIndVectorCont[scenario]).push_back(zonalIndex);
	phAngDecVec.at(scenario*zoneNumber*nodeNumber+zonalIndex*nodeNumber+globalRank)=interimDecAngle;
}

void Marketover::rewardPenaltyUpdate(double lagrangeMultXi[], double lagrangeMultPi[], int matrixIndex, int iteration) // Function getGlobalUpper() returns the global upper bound for the investment coordination problem
{
/*	for (int scenPos = 0; scenPos < countOfScenarios; ++scenPos) {
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
	}	*/
}

double Marketover::rewardPenaltyCont(double lagrangeMult, int matrixIndex, int iteration) // Function getGlobalUpper() returns the global upper bound for the investment coordination problem
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

double Marketover::rewardPenaltyInteger(double lagrangeMult, int matrixIndex, int iteration) // Function getGlobalUpper() returns the global upper bound for the investment coordination problem
{
	//cout << "\nITERATION COUNT : " << iteration << endl;
	int tracker = 0;
	vector<int>::iterator lineDecIndexIterator;
	for (lineDecIndexIterator=lineDecIndex.begin();lineDecIndexIterator!=lineDecIndex.end();++lineDecIndexIterator) {
		if (matrixIndex==(zonalIndVectorInt[tracker]-1)*candLineNumber+(*lineDecIndexIterator))	{
			if (compareBasis == 1) {		
				lagrangeMult += (lineInterDecision[tracker]-interimIntDecVar[(*lineDecIndexIterator-1)]);
				//cout << "From MO Integer Lagrange Multiplier update: " << "Line Decision from Zones. Tracker: " << tracker << " Zonal Decision value " << lineInterDecision[tracker] << " MO decision index: " << (*lineDecIndexIterator-1) << " MO Decision value: " << interimIntDecVar[(*lineDecIndexIterator-1)] << " Updated Lagrange Multiplier: " << lagrangeMult << endl; 
			}
			else {
				lagrangeMult += (lineInterDecision[tracker]-interimIntDecVarPrev[(*lineDecIndexIterator-1)]);
				//cout << "From MO Integer Lagrange Multiplier update: " << "Line Decision from Zones. Tracker: " << tracker << " Zonal Decision value " << lineInterDecision[tracker] << " MO decision index: " << (*lineDecIndexIterator-1) << " MO Decision value: " << interimIntDecVar[(*lineDecIndexIterator-1)] << " Updated Lagrange Multiplier: " << lagrangeMult << endl; 
			}
		}
		++tracker;
	}
	//cout << "Tracker " << tracker << endl;
	return lagrangeMult;	
} // Function getGlobalUpper ends

void Marketover::clearVectors() // Clears the different interim vectors for making them ready for the next iteration
{
	for (int scenPos = 0; scenPos < countOfScenarios; ++scenPos) {
		// Clear vectors for interim continuous decision variables
		(phaseAngleDecision[scenPos]).clear();
		(angleDecIndex[scenPos]).clear();
		(zonalIndVectorCont[scenPos]).clear();
	}
	interimContDecVar.clear();
	interimContDecVarPrev.clear();
	// Clear vectors for interim integer decision variables
	lineInterDecision.clear();
	lineDecIndex.clear();
	zonalIndVectorInt.clear();
	interimIntDecVar.clear();
	interimIntDecVarPrev.clear();
}

void Marketover::clearDelayedVectors() // Clears the different interim vectors only buffer vectors
{
	interimContDecVarPrev.clear();
	interimIntDecVarPrev.clear();
}

void Marketover::MILPMarketoverGUROBI(double LagMultXi[], double LagMultPi[], int totalCandLineNum, int totalSharedNodeNum, GRBEnv* environmentGUROBI) // Function MILPMarketover() implements the Mixed Integer Linear Programming Solver routine by calling GUROBI routines for average heat rate objective
{
	// CREATION OF THE MIP SOLVER INSTANCE //
	clock_t begin = clock(); // start the timer
	vector<int>::iterator diffZNIt; // Iterator for diffZoneNodeID
	vector<Node*>::iterator nodeIterator; // Iterator for node objects
	vector<double>::iterator candCapIterator; // Iterator for candidate lines for capacity
	vector<int>::iterator candIterator; // Iterator for candidate lines	
	vector<double>::iterator exsharedCapIterator; // Iterator for shared existing lines for capacity
	vector<int>::iterator exsharedIterator;

	string outSummaryFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputSummaryGUROBI/OutSummaryMOGUROBI.txt";
	ofstream outPutFile(outSummaryFileName, ios::out); // Create Output File to output the Summary of Results
	if (!outPutFile){
		cerr << "\nCouldn't open the file" << endl;
		exit(1);
	}

        int dimRow = countOfScenarios*(4 * candLineNumber + 2 * sharedELineNumber); // Total number of rows of the A matrix (number of structural constraints of the LP) first term to account for lower and upper line limits & lower and upper definition limits of candidate shared lines, and second term for lower and upper line limits for shared existing lines
        int dimCol = countOfScenarios*nodeNumber+(countOfScenarios+1)*candLineNumber; // Total number of columns of the LP (number of Decision Variables) first term to account for voltage phase angles for inter-zonal lines' nodes, and second term for power flow values and binary integer decision variable values for shared candidate lines

	outPutFile << "\nTotal Number of Structural Constraints (Rows) is: " << dimRow << endl;
	outPutFile << "\nTotal Number of Decision Variables (Columns) is: " << dimCol << endl;
	//cout << "\nTotal Number of Structural Constraints (Rows) is: " << dimRow << endl;
	//cout << "\nTotal Number of Decision Variables (Columns) is: " << dimCol << endl;
	//cout << "\nTotal Number of Shared Nodes is: " << nodeNumber << endl;	
	//cout << "\nTotal Number of shared candidate lines is: " << candLineNumber << endl;
	//cout << "\nTotal Number of shared existing lines is: " << sharedELineNumber << endl;
	GRBModel *modelMOMILP = new GRBModel(*environmentGUROBI);
    	modelMOMILP->set(GRB_StringAttr_ModelName, "assignmentMO");
	GRBVar decvar[dimCol+1];
	double z; // variable to store the objective value
	vector<double> x; // Coefficient matrix entires, objective function, and decision variables

	// SPECIFICATION OF PROBLEM PARAMETERS //
	// Dummy Decision Variable //
	decvar[0] = modelMOMILP->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);

	//Decision Variable Definitions, Bounds, and Objective Function Co-efficients//
	int colCount = 1;

	//Columns corresponding to Shared Candidate Line Flows continuous variables//
	outPutFile << "\nColumns corresponding to Shared Candidate Line Flows continuous variables" << endl;
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candCapIterator = candCapacity.begin(); candCapIterator != candCapacity.end(); ++candCapIterator){
			decvar[colCount] = modelMOMILP->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
			outPutFile << colCount << "\t";
			outPutFile << 0 << endl;
			++colCount;
		}
	}
	outPutFile << "\nTotal number of columns after accounting for Shared Candidate Line Flows continuous variables: " << colCount << endl;

	//Columns corresponding to Shared Candidate Line Construction Decision Binary Integer variables//
	outPutFile << "\nColumns corresponding to Shared Candidate Line Construction Decision Binary Integer variables" << endl;
	int candInd = 0; // Initialize the counter for indexing the candidate lines
	int diffCandCounter = 0; // counter flag to indicate the first element of the candSerial list
	for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
		if (diffCandCounter > 0) { // Skip the first element, since it's a dummy "0"
			decvar[colCount] = modelMOMILP->addVar(0, 1, 0.0, GRB_BINARY);
			outPutFile << colCount << "\t";
			outPutFile << LagMultPi[candInd] << "\t" << candInd << endl;
			++candInd;
			++colCount;
		}
		++diffCandCounter;
	}
	outPutFile << "\nTotal number of columns after accounting for Shared Candidate Line Construction Decision Binary Integer variables: " << colCount << endl;

	//Columns corresponding to Voltage Phase Angles continuous variables for other zone nodes for shared lines//
	outPutFile << "\nColumns corresponding to Voltage Phase Angles continuous variables for other zone nodes for shared lines" << endl;
	int otherInd = 0; // Initialize the counter for indexing the other-zone nodes
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int columnCount = 1;
		int diffNodeCounter = 0; // counter flag to indicate the first element of the diffZoneNodeID list
		//for (diffZNIt = diffZoneNodeID.begin(); diffZNIt != diffZoneNodeID.end(); ++diffZNIt){
			//cout << " Column count after the shared node is " << columnCount << endl;	
			//++columnCount;
		//}//
		for (diffZNIt = nodeList.begin(); diffZNIt != nodeList.end(); ++diffZNIt){
			if (diffNodeCounter > 0) { // Skip the first element, since it's a dummy "0"
				decvar[colCount] = modelMOMILP->addVar(0, (44/7), 0.0, GRB_CONTINUOUS);
				outPutFile << colCount << "\t";
				outPutFile << LagMultXi[otherInd] << "\t" << otherInd << endl;
				++otherInd;	
				++colCount;
			}
			++diffNodeCounter;
		}
	}
	outPutFile << "\nDecision Variables and Objective Function defined" << endl;
	outPutFile << "\nTotal Number of columns: " << colCount - 1 << endl;
	//cout << "\nDecision Variables and Objective Function defined" << endl;
	//cout << "\nTotal Number of columns: " << colCount - 1 << endl;
	//Setting Objective//
	GRBLinExpr obj = 0.0;
	// Objective Contribution from Dummy Decision Variable //
	obj += 0*(decvar[0]);
	colCount = 1;

	//Columns corresponding to Shared Candidate Line Flows continuous variables//
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candCapIterator = candCapacity.begin(); candCapIterator != candCapacity.end(); ++candCapIterator){
			obj += 0*(decvar[colCount]);
			++colCount;
		}
	}
	//Columns corresponding to Shared Candidate Line Construction Decision Binary Integer variables//
	candInd = 0; // Initialize the counter for indexing the candidate lines
	diffCandCounter = 0; // counter flag to indicate the first element of the candSerial list
	for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
		if (diffCandCounter > 0) { // Skip the first element, since it's a dummy "0"
			obj += (LagMultPi[candInd])*(decvar[colCount]);
			++candInd;
			++colCount;
		}
		++diffCandCounter;
	}

	//Columns corresponding to Voltage Phase Angles continuous variables for other zone nodes for shared lines//
	otherInd = 0; // Initialize the counter for indexing the other-zone nodes
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int columnCount = 1;
		int diffNodeCounter = 0; // counter flag to indicate the first element of the diffZoneNodeID list
		for (diffZNIt = nodeList.begin(); diffZNIt != nodeList.end(); ++diffZNIt){
			if (diffNodeCounter > 0) { // Skip the first element, since it's a dummy "0"
				obj += (LagMultXi[otherInd])*(decvar[colCount]);
				++otherInd;	
				++colCount;
			}
			++diffNodeCounter;
		}
	}
	//cout << "\nDecision Variables and Objective Function defined" << endl;
	//cout << "\nTotal Number of columns: " << colCount - 1 << endl;

	modelMOMILP->setObjective(obj, GRB_MAXIMIZE);
	//cout << " Objective Function and Decision Variables have been defined and the colCount is " << colCount-1 << endl;

	//Row Definitions: Specification of RHS or b vector of b<=Ax<=b//
	GRBLinExpr lhs[dimRow+1];
	//Row Definitions and Bounds Corresponding to Constraints/
	string outPGenFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputPowerGUROBI/OutPowerGenMOGUROBI.txt"; 
	ofstream powerGenOut(outPGenFileName, ios::out);
	if (!powerGenOut){
		cerr << "\nCouldn't open the file" << endl;
		exit(1);
	}
	outPutFile << "\nTesting of interzonal transmission limits" << endl; 
	// Dummy Constraint //
	lhs[0] = 0*(decvar[0]);
	modelMOMILP->addConstr(lhs[0], GRB_EQUAL, 0);


	//Non-Zero entries of A matrix (Constraint/Coefficient matrix entries)//
	int rCount = 1; // Initialize the row count
	
	// Coefficients corresponding to shared existing Line Forward Flow Limit Constraints
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared existing Line Forward Flow Limit Constraints" << endl;
		for (exsharedIterator = SESerial.begin(); exsharedIterator != SESerial.end(); ++exsharedIterator){
			lhs[rCount] = 0;
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				//cout << "Test Message 2" << endl;
				//cout << "SE From Rank" << SEFromRank.at(jockey) << endl;
				//cout << "SE Reactance" << 1/(SEReactance.at(jockey)) << endl;
				lhs[rCount] += (1/(SEReactance.at(jockey)))*(decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEFromRank.at(jockey))]);
				outPutFile << "\n" << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEFromRank.at(jockey)) << "\t" << 1/(SEReactance.at(jockey)) << endl;
				lhs[rCount] += (-1/(SEReactance.at(jockey)))*(decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEToRank.at(jockey))]);
				outPutFile << "\n" << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEToRank.at(jockey)) << "\t" << -1/(SEReactance.at(jockey)) << endl;
				modelMOMILP->addConstr(lhs[rCount] <= (SECapacity.at(jockey)));
				outPutFile << rCount << "\t";
				outPutFile << (SECapacity.at(jockey)) << endl;
				++rCount; // Increment the row count to point to the next transmission line object
				++jockey;
				//cout << "\nTest after shared existing Line Forward Flow Limit Constraint " << diffSerCounter << endl;
			}
			++diffSerCounter;	
		}
	}
	// Coefficients corresponding to shared existing Line Reverse Flow Limit Constraints
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared existing Line Reverse Flow Limit Constraints" << endl;
		for (exsharedIterator = SESerial.begin(); exsharedIterator != SESerial.end(); ++exsharedIterator){
			lhs[rCount] = 0;
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				lhs[rCount] += (1/(SEReactance.at(jockey)))*(decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEFromRank.at(jockey))]);
				outPutFile << "\n" << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEFromRank.at(jockey)) << "\t" << 1/(SEReactance.at(jockey)) << endl;
				lhs[rCount] += (-1/(SEReactance.at(jockey)))*(decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEToRank.at(jockey))]);
				outPutFile << "\n" << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEToRank.at(jockey)) << "\t" << -1/(SEReactance.at(jockey)) << endl;
				modelMOMILP->addConstr(lhs[rCount] >= -(SECapacity.at(jockey)));
				outPutFile << rCount << "\t";
				outPutFile << -(SECapacity.at(jockey)) << endl;
				++rCount; // Increment the row count to point to the next transmission line object
				++jockey;
				//cout << "\nTest after shared existing Line Reverse Flow Limit Constraint " << diffSerCounter << endl;
			}
			++diffSerCounter;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for Line Flow limits for shared existing lines: " << rCount << endl;
	
	// Coefficients corresponding to shared candidate Line Forward Flow Limit Constraints
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared candidate Line Forward Flow Limit Constraints" << endl;
		for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
			lhs[rCount] = 0;
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				lhs[rCount] += (1)*(decvar[(rCount-2*countOfScenarios*sharedELineNumber)]);
				outPutFile << "\n" << rCount << "\t" << (rCount-2*countOfScenarios*sharedELineNumber) << "\t" << 1 << endl;
				lhs[rCount] += (-candCapacity.at(jockey))*(decvar[(countOfScenarios*candLineNumber+rCount-2*countOfScenarios*sharedELineNumber-scenCounter*candLineNumber)]);
				outPutFile << "\n" << rCount << "\t" << (countOfScenarios*candLineNumber+rCount-2*countOfScenarios*sharedELineNumber-scenCounter*candLineNumber) << "\t" << -candCapacity.at(jockey) << endl;
				modelMOMILP->addConstr(lhs[rCount] <= 0.0);
				outPutFile << rCount << "\t";
				outPutFile << 0.0 << endl;
				++rCount;
				++jockey;
				//cout << "\nTest after shared candidate Line Forward Flow Limit Constraint " << diffSerCounter << endl;
			}
			++diffSerCounter;
		}
	}
	// Coefficients corresponding to shared candidate Line Reverse Flow Limit Constraints
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared candidate Line Reverse Flow Limit Constraints" << endl;
		for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
			lhs[rCount] = 0;
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				lhs[rCount] += (1)*(decvar[(rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber))]);
				outPutFile << "\n" << rCount << "\t" << (rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber)) << "\t" << 1 << endl;
				lhs[rCount] += (candCapacity.at(jockey))*(decvar[(rCount-countOfScenarios*(2*sharedELineNumber)-scenCounter*candLineNumber)]);
				outPutFile << "\n" << rCount << "\t" << (rCount-countOfScenarios*(2*sharedELineNumber)-scenCounter*candLineNumber) << "\t" << candCapacity.at(jockey) << endl;
				modelMOMILP->addConstr(lhs[rCount] >= 0.0);
				outPutFile << rCount << "\t";
				outPutFile << 0.0 << endl;
				++rCount;
				++jockey;
				//cout << "\nTest after shared candidate Line Reverse Flow Limit Constraint " << diffSerCounter << endl;
			}
			++diffSerCounter;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for Line Flow limits for shared candidate lines: " << rCount << endl;
	
	// Coefficients corresponding to shared candidate Line Definition upper bound
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared candidate Line Definition upper bound" << endl;
		for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
			lhs[rCount] = 0;
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				lhs[rCount] += (1)*(decvar[(rCount-countOfScenarios*(2*sharedELineNumber+2*candLineNumber))]);
				outPutFile << rCount << "\t" << (rCount-countOfScenarios*(2*sharedELineNumber+2*candLineNumber)) << " Power term for candidate line " << diffSerCounter << "-th candidate line" << endl;
				lhs[rCount] += (-1/(candReactance.at(jockey)))*(decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candFromRank.at(jockey))]);
				outPutFile << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candFromRank.at(jockey)) << "\t" << -1/(candReactance.at(jockey)) << endl;
				lhs[rCount] += (1/(candReactance.at(jockey)))*(decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candToRank.at(jockey))]);
				outPutFile << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candToRank.at(jockey)) << "\t" << 1/(candReactance.at(jockey)) << endl;
				lhs[rCount] += (2.5*(candCapacity.at(jockey)))*(decvar[(rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber)-scenCounter*candLineNumber)]);
				outPutFile << rCount << "\t" << (rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber)-scenCounter*candLineNumber) << "\t" << BIGM << endl;
				modelMOMILP->addConstr(lhs[rCount] <= 2.5*(candCapacity.at(jockey)));
				outPutFile << rCount << "\t";
				outPutFile << BIGM << endl;
				++rCount; // Increment the row count to point to the next transmission line object
				++jockey;
				//cout << "\nTest after shared candidate Line Definition upper bound " << diffSerCounter << endl;
			}
			++diffSerCounter;
		}
	}
	// Coefficients corresponding to shared candidate Line Definition lower bound
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared candidate Line Definition lower bound" << endl;
		for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
			lhs[rCount] = 0;
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				lhs[rCount] += (1)*(decvar[(rCount-countOfScenarios*(2*sharedELineNumber+3*candLineNumber))]);
				outPutFile << rCount << "\t" << (rCount-countOfScenarios*(2*sharedELineNumber+3*candLineNumber)) << " Power term for candidate line " << diffSerCounter << "-th candidate line" << endl;
				lhs[rCount] += (-1/(candReactance.at(jockey)))*(decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candFromRank.at(jockey))]);
				outPutFile << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candFromRank.at(jockey)) << "\t" << -1/(candReactance.at(jockey)) << endl;
				lhs[rCount] += (1/(candReactance.at(jockey)))*(decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candToRank.at(jockey))]);
				outPutFile << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candToRank.at(jockey)) << "\t" << 1/(candReactance.at(jockey)) << endl;
				lhs[rCount] += (-(2.5*(candCapacity.at(jockey))))*(decvar[(rCount-countOfScenarios*(2*sharedELineNumber+2*candLineNumber)-scenCounter*candLineNumber)]);
				outPutFile << rCount << "\t" << (rCount-countOfScenarios*(2*sharedELineNumber+2*candLineNumber)-scenCounter*candLineNumber) << "\t" << -BIGM << endl;
				modelMOMILP->addConstr(lhs[rCount] >= -(2.5*(candCapacity.at(jockey))));
				outPutFile << rCount << "\t";
				outPutFile << -BIGM << endl;
				++rCount; // Increment the row count to point to the next transmission line object
				++jockey;
				//cout << "\nTest after shared candidate Line Definition lower bound " << diffSerCounter << endl;
			}
			++diffSerCounter;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for shared candidate lines flow definition: " << rCount << endl;
	outPutFile << "\nConstraint bounds (rows) Specified" << endl;
	outPutFile << "\nTotal number of rows: " << rCount - 1 << endl;
	//cout << "\nConstraint bounds (rows) Specified" << endl;
	//cout << "\nTotal number of rows: " << rowCount - 1 << endl;

	outPutFile << "\nCoefficient Matrix specified" << endl;
	clock_t end1 = clock(); // stop the timer
	double elapsed_secs1 = double(end1 - begin) / CLOCKS_PER_SEC; // Calculate the time required to populate the constraint matrix and objective coefficients
	outPutFile << "\nTotal time taken to define the rows, columns, objective and populate the coefficient matrix = " << elapsed_secs1 << " s " << endl;
	// RUN THE OPTIMIZATION SIMULATION ALGORITHM //
	cout << "\nSimulation in Progress. Wait !!! ....." << endl;
	modelMOMILP->optimize(); // Solves the optimization problem
	int stat = modelMOMILP->get(GRB_IntAttr_Status); // Outputs the solution status of the problem 

	// DISPLAY THE SOLUTION DETAILS //
	if (stat == GRB_INFEASIBLE){
		outPutFile << "\nThe solution to the problem is INFEASIBLE." << endl;
		cout << "\nThe solution to the problem is INFEASIBLE." << endl;
		delete modelMOMILP; // Free the memory of the GUROBI Problem Model
	} else if (stat == GRB_INF_OR_UNBD) {
		outPutFile << "\nNO FEASIBLE or BOUNDED solution to the problem exists." << endl;
		cout << "\nNO FEASIBLE or BOUNDED solution to the problem exists." << endl;
		delete modelMOMILP; // Free the memory of the GUROBI Problem Model
	} else if (stat == GRB_UNBOUNDED) {
		outPutFile << "\nThe solution to the problem is UNBOUNDED." << endl;
		cout << "\nThe solution to the problem is UNBOUNDED." << endl;
		delete modelMOMILP; // Free the memory of the GUROBI Problem Model
	} else if (stat == GRB_OPTIMAL) {
		outPutFile << "\nThe solution to the problem is OPTIMAL." << endl;
		cout << "\nThe solution to the problem is OPTIMAL." << endl;

		//Get the Optimal Objective Value results//
		z = modelMOMILP->get(GRB_DoubleAttr_ObjVal);

		// Open separate output files for writing results of different variables
		string outCandFlowFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/candidateLinesGUROBI/candFlowMWMOGUROBI.txt";
		string outCandDecFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/candidateLinesGUROBI/candLineDecisionMOGUROBI.txt";
		string outExtAngFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputAnglesGUROBI/externalAngleMOGUROBI.txt";
		string constraintSat = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/ConstraintSatGUROBI/constraintSatGUROBI.txt";
		ofstream candFlowMWOut(outCandFlowFileName, ios::out); //switchOnOut
		ofstream candLineDecisionOut(outCandDecFileName, ios::out); //switchOffOut
		ofstream externalAngleOut(outExtAngFileName, ios::out);
		ofstream constSatOut(constraintSat, ios::out);
		outPutFile << "\nThe Optimal Objective value (Line Building Decision cost) is: " << -z << endl;
		powerGenOut << "\nThe Optimal Objective value (Line Building Decision cost) is: " << -z << endl;
		cout << "\nThe Optimal Objective value (Line Building Decision cost) is: " << -z << endl;
		x.push_back(0); // Initialize the decision Variable vector

		// Display Shared Candidate lines' Power Flow variables
		candFlowMWOut << "\n****************** SHARED CANDIDATE LINES POWER FLOW VALUES *********************" << endl;
		candFlowMWOut << "CANDIDATE LINE ID" << "\t" << "MW FLOW VALUE" << "\n";
		int arrayInd = 1;
        	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
			int diffCandCounter2 = 0; // counter flag to indicate the first element of the candSerial list
			for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
				if (diffCandCounter2 > 0) { // Skip the first element, since it's a dummy "0"
					x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
					candFlowMWOut << (*candIterator) << "\t" << ((decvar[arrayInd]).get(GRB_DoubleAttr_X))*100 << " MW" << endl;
					++arrayInd;
				}
				++diffCandCounter2;
			}
		}
		candFlowMWOut << "Finished writing Shared Candidate lines' Power Flow variables" << endl;

		// Display Shared Candidate lines' Construction Decisions
		candLineDecisionOut << "\n****************** SHARED CANDIDATE LINES CONSTRUCTION DECISIONS *********************" << endl;
		candLineDecisionOut << "CANDIDATE LINE ID" << "\t" << "CONSTRUCTION DECISIONS" << "\n";
		int candInd1 = 0; // Initialize the counter for indexing the candidate lines
		int diffCandCounter1 = 0; // counter flag to indicate the first element of the candSerial list
		for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
			if (diffCandCounter1 > 0) { // Skip the first element, since it's a dummy "0"
				x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
				interimIntDecVar.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
				candLineDecisionOut << (*candIterator) << "\t" << ((decvar[arrayInd]).get(GRB_DoubleAttr_X)) << endl;
				++arrayInd;
				++candInd1;
			}
			++diffCandCounter1;
		}
		candLineDecisionOut << "Finished writing Shared Candidate lines' Construction decisions" << endl;

		// Display shared node angles
		externalAngleOut << "\n****************** OUTER ZONAL NODE VOLTAGE PHASE ANGLE VALUES *********************" << endl;
		externalAngleOut << "SHARED NODE GLOBAL RANK" << "\t" << "SHARED NODE LOCAL RANK" << "\t" <<"SHARED NODE ZONE ID" << "\t" <<"VOLTAGE PHASE ANGLE" << "\n";
	        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
			int otherInd1 = 0; // Initialize the counter for indexing the other-zone nodes
			int dummyStart = 0;
			for (diffZNIt = sharedGlobalList.begin(); diffZNIt != sharedGlobalList.end(); ++diffZNIt){
				if (dummyStart > 0) { // Skip the dummy element "0" at the beginning
					x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
					interimContDecVar.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
					externalAngleOut << (*diffZNIt) << "\t" << nodeList[otherInd1] << "\t" << zoneList[otherInd1] << "\t" << ((decvar[arrayInd]).get(GRB_DoubleAttr_X)) << endl;
					++arrayInd;
				}
				++dummyStart;
				++otherInd1;
			}
		}
	        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
			int otherInd1 = 0; // Initialize the counter for indexing the other-zone nodes
			int dummyStart = 0;
			for (diffZNIt = nodeList.begin(); diffZNIt != nodeList.end(); ++diffZNIt){
				//if (dummyStart > 0) { // Skip the dummy element "0" at the beginning
					//x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
					//interimContDecVar.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
					//*cout << "\nMessage from MO. The rank of shared node is " << (*diffZNIt) << endl;
					//++arrayInd;
					//++otherInd1;
				//}
				++dummyStart;
			}
			for (diffZNIt = sharedGlobalList.begin(); diffZNIt != sharedGlobalList.end(); ++diffZNIt){
				//if (dummyStart > 0) { // Skip the dummy element "0" at the beginning
					//x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
					//interimContDecVar.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
					//* cout << "\nMessage from MO. The rank of shared zone is " << (*diffZNIt) << endl;
					//++arrayInd;
					//++otherInd1;
				//}
				++dummyStart;
			}
		}
		externalAngleOut << "Finished writing shared node voltage phase angle values" << endl;

		constSatOut << "TESTING FOR CONSTRAINT SATISFACTION" << endl;
		double LHS[dimRow+1];
		constSatOut << "\nTesting of interzonal transmission limits" << endl; 
		// Dummy Constraint //
		LHS[0] = 0*((decvar[0]).get(GRB_DoubleAttr_X));
		//Non-Zero entries of A matrix (Constraint/Coefficient matrix entries)//
		int rCount = 1; // Initialize the row count	
		// Coefficients corresponding to shared existing Line Forward Flow Limit Constraints
       	 	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
			int jockey = 0;
			int diffSerCounter = 0;
			constSatOut << "\nCoefficients corresponding to shared existing Line Forward Flow Limit Constraints" << endl;
			for (exsharedIterator = SESerial.begin(); exsharedIterator != SESerial.end(); ++exsharedIterator){
				LHS[rCount] = 0;
				if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
					constSatOut << "\nRow: " << "\tColumn(From Node): " << "\tFrom Reciprocal Reactance" << "\tFrom Angle" << "\tLHS value" << endl;
					LHS[rCount] += (1/(SEReactance.at(jockey)))*((decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEFromRank.at(jockey))]).get(GRB_DoubleAttr_X));
					constSatOut << "\n" << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEFromRank.at(jockey)) << "\t" << 1/(SEReactance.at(jockey)) << "\t" << ((decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEFromRank.at(jockey))]).get(GRB_DoubleAttr_X)) << "\t" << LHS[rCount] << endl;
					constSatOut << "\nRow: " << "\tColumn(To Node): " << "\tTo Reciprocal Reactance" << "\tTo Angle" << "\tLHS value" << endl;
					LHS[rCount] += (-1/(SEReactance.at(jockey)))*((decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEToRank.at(jockey))]).get(GRB_DoubleAttr_X));
					constSatOut << "\n" << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEToRank.at(jockey)) << "\t" << 1/(SEReactance.at(jockey)) << "\t" << ((decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEToRank.at(jockey))]).get(GRB_DoubleAttr_X)) << "\t" << LHS[rCount] << endl;
					constSatOut << LHS[rCount] << "<=" << (SECapacity.at(jockey)) << endl;
					++rCount; // Increment the row count to point to the next transmission line object
					++jockey;
				}
				++diffSerCounter;	
			}
		}
		// Coefficients corresponding to shared existing Line Reverse Flow Limit Constraints
        	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
			int jockey = 0;
			int diffSerCounter = 0;
			constSatOut << "\nCoefficients corresponding to shared existing Line Reverse Flow Limit Constraints" << endl;
			for (exsharedIterator = SESerial.begin(); exsharedIterator != SESerial.end(); ++exsharedIterator){
				LHS[rCount] = 0;
				if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
					constSatOut << "\nRow: " << "\tColumn(From Node): " << "\tFrom Reciprocal Reactance" << "\tFrom Angle" << "\tLHS value" << endl;
					LHS[rCount] += (1/(SEReactance.at(jockey)))*((decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEFromRank.at(jockey))]).get(GRB_DoubleAttr_X));
					constSatOut << "\n" << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEFromRank.at(jockey)) << "\t" << 1/(SEReactance.at(jockey)) << "\t" << ((decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEFromRank.at(jockey))]).get(GRB_DoubleAttr_X)) << "\t" << LHS[rCount] << endl;
					constSatOut << "\nRow: " << "\tColumn(To Node): " << "\tTo Reciprocal Reactance" << "\tTo Angle" << "\tLHS value" << endl;
					LHS[rCount] += (-1/(SEReactance.at(jockey)))*((decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEToRank.at(jockey))]).get(GRB_DoubleAttr_X));
					constSatOut << "\n" << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEToRank.at(jockey)) << "\t" << 1/(SEReactance.at(jockey)) << "\t" << ((decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEToRank.at(jockey))]).get(GRB_DoubleAttr_X)) << "\t" << LHS[rCount] << endl;
					constSatOut << LHS[rCount] << ">=" << (SECapacity.at(jockey)) << endl;
					++rCount; // Increment the row count to point to the next transmission line object
					++jockey;
					//cout << "\nTest after shared existing Line Reverse Flow Limit Constraint " << diffSerCounter << endl;
				}
				++diffSerCounter;
			}
		}
		constSatOut << "\nTotal number of rows after accounting for Line Flow limits for shared existing lines: " << rCount << endl;
	
		// Coefficients corresponding to shared candidate Line Forward Flow Limit Constraints
        	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
			int jockey = 0;
			int diffSerCounter = 0;
			constSatOut << "\nCoefficients corresponding to shared candidate Line Forward Flow Limit Constraints" << endl;
			for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
				LHS[rCount] = 0;
				if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
					constSatOut << "\nRow: " << "\tColumn(flow continous variable): " << "\tFlow Value" << "\tLHS value" << endl;
					LHS[rCount] += (1)*((decvar[(rCount-2*countOfScenarios*sharedELineNumber)]).get(GRB_DoubleAttr_X));
					constSatOut << "\n" << rCount << "\t" << (rCount-2*countOfScenarios*sharedELineNumber) << "\t" << ((decvar[(rCount-2*countOfScenarios*sharedELineNumber)]).get(GRB_DoubleAttr_X)) << LHS[rCount] << endl;
					constSatOut << "\nRow: " << "\tColumn(Integer variable): " << "\tLine Building decision" << "\tLHS value" << endl;
					LHS[rCount] += (-candCapacity.at(jockey))*((decvar[(countOfScenarios*candLineNumber+rCount-2*countOfScenarios*sharedELineNumber-scenCounter*candLineNumber)]).get(GRB_DoubleAttr_X));
					constSatOut << "\n" << rCount << "\t" << (countOfScenarios*candLineNumber+rCount-2*countOfScenarios*sharedELineNumber-scenCounter*candLineNumber) << "\t" << ((decvar[(countOfScenarios*candLineNumber+rCount-2*countOfScenarios*sharedELineNumber-scenCounter*candLineNumber)]).get(GRB_DoubleAttr_X)) << "\t" << LHS[rCount] << endl;
					constSatOut << LHS[rCount] << "<=" << 0 << endl;
					++rCount;
					++jockey;
					//cout << "\nTest after shared candidate Line Forward Flow Limit Constraint " << diffSerCounter << endl;
				}
				++diffSerCounter;
			}
		}
		// Coefficients corresponding to shared candidate Line Reverse Flow Limit Constraints
        	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
			int jockey = 0;
			int diffSerCounter = 0;
			constSatOut << "\nCoefficients corresponding to shared candidate Line Reverse Flow Limit Constraints" << endl;
			for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
				LHS[rCount] = 0;
				if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
					constSatOut << "\nRow: " << "\tColumn(flow continous variable): " << "\tFlow Value" << "\tLHS value" << endl;
					LHS[rCount] += (1)*((decvar[(rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber))]).get(GRB_DoubleAttr_X));
					constSatOut << "\n" << rCount << "\t" << (rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber)) << "\t" << ((decvar[(rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber))]).get(GRB_DoubleAttr_X)) << LHS[rCount] << endl;
					constSatOut << "\nRow: " << "\tColumn(Integer variable): " << "\tLine Building decision" << "\tLHS value" << endl;
					LHS[rCount] += (candCapacity.at(jockey))*((decvar[(rCount-countOfScenarios*(2*sharedELineNumber)-scenCounter*candLineNumber)]).get(GRB_DoubleAttr_X));
					constSatOut << "\n" << rCount << "\t" << (rCount-countOfScenarios*(2*sharedELineNumber)-scenCounter*candLineNumber) << "\t" << ((decvar[(rCount-countOfScenarios*(2*sharedELineNumber)-scenCounter*candLineNumber)]).get(GRB_DoubleAttr_X)) << "\t" << LHS[rCount] << endl;
					constSatOut << LHS[rCount] << ">=" << 0 << endl;
					++rCount;
					++jockey;
					//cout << "\nTest after shared candidate Line Reverse Flow Limit Constraint " << diffSerCounter << endl;
				}
				++diffSerCounter;
			}
		}
		constSatOut << "\nTotal number of rows after accounting for Line Flow limits for shared candidate lines: " << rCount << endl;
		constSatOut << " Coefficients corresponding to shared candidate Line Definition upper bound" << endl;
        	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
			int jockey = 0;
			int diffSerCounter = 0;
			constSatOut << "\nCoefficients corresponding to shared candidate Line Definition upper bound" << endl;
			for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
				LHS[rCount] = 0;
				if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
					constSatOut << "\nRow: " << "\tColumn(flow continous variable): " << "\tFlow Value" << "\tLHS value" << endl;
					LHS[rCount] += (1)*((decvar[(rCount-countOfScenarios*(2*sharedELineNumber+2*candLineNumber))]).get(GRB_DoubleAttr_X));
					constSatOut << "\n" << rCount << "\t" << (rCount-countOfScenarios*(2*sharedELineNumber+2*candLineNumber)) << "\t" << ((decvar[(rCount-countOfScenarios*(2*sharedELineNumber+2*candLineNumber))]).get(GRB_DoubleAttr_X)) << "\t" << LHS[rCount] << endl;
					constSatOut << "\nRow: " << "\tColumn(From Node): " << "\tFrom Reciprocal Reactance" << "\tFrom Angle" << "\tLHS value" << endl;
					LHS[rCount] += (-1/(candReactance.at(jockey)))*((decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candFromRank.at(jockey))]).get(GRB_DoubleAttr_X));
					constSatOut << "\n" << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candFromRank.at(jockey)) << "\t" << 1/(candReactance.at(jockey)) << "\t" << ((decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candFromRank.at(jockey))]).get(GRB_DoubleAttr_X)) << "\t" << LHS[rCount] << endl;
					constSatOut << "\nRow: " << "\tColumn(To Node): " << "\tTo Reciprocal Reactance" << "\tTo Angle" << "\tLHS value" << endl;
					LHS[rCount] += (1/(candReactance.at(jockey)))*((decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candToRank.at(jockey))]).get(GRB_DoubleAttr_X));
					constSatOut << "\n" << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candToRank.at(jockey)) << "\t" << 1/(candReactance.at(jockey)) << "\t" << ((decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candToRank.at(jockey))]).get(GRB_DoubleAttr_X)) << "\t" << LHS[rCount] << endl;
					constSatOut << "\nRow: " << "\tColumn(Integer variable): " << "\tLine Building decision" << "\tLHS value" << endl;
					LHS[rCount] += (2.5*(candCapacity.at(jockey)))*((decvar[(rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber)-scenCounter*candLineNumber)]).get(GRB_DoubleAttr_X));
					constSatOut << rCount << "\t" << (rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber)-scenCounter*candLineNumber) << "\t" << ((decvar[(rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber)-scenCounter*candLineNumber)]).get(GRB_DoubleAttr_X)) << "\t" << LHS[rCount] << endl;
					constSatOut << LHS[rCount] << "<=" << 2.5*(candCapacity.at(jockey)) << endl;
					++rCount; // Increment the row count to point to the next transmission line object
					++jockey;
					//cout << "\nTest after shared candidate Line Definition upper bound " << diffSerCounter << endl;
				}
				++diffSerCounter;
			}
		}
		// Coefficients corresponding to shared candidate Line Definition lower bound
        	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
			int jockey = 0;
			int diffSerCounter = 0;
			constSatOut << "\nCoefficients corresponding to shared candidate Line Definition lower bound" << endl;
			for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
				LHS[rCount] = 0;
				if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
					constSatOut << "\nRow: " << "\tColumn(flow continous variable): " << "\tFlow Value" << "\tLHS value" << endl;
					LHS[rCount] += (1)*((decvar[(rCount-countOfScenarios*(2*sharedELineNumber+3*candLineNumber))]).get(GRB_DoubleAttr_X));
					constSatOut << "\n" << rCount << "\t" << (rCount-countOfScenarios*(2*sharedELineNumber+3*candLineNumber)) << "\t" << ((decvar[(rCount-countOfScenarios*(2*sharedELineNumber+3*candLineNumber))]).get(GRB_DoubleAttr_X)) << "\t" << LHS[rCount] << endl;
					constSatOut << "\nRow: " << "\tColumn(From Node): " << "\tFrom Reciprocal Reactance" << "\tFrom Angle" << "\tLHS value" << endl;
					LHS[rCount] += (-1/(candReactance.at(jockey)))*((decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candFromRank.at(jockey))]).get(GRB_DoubleAttr_X));
					constSatOut << "\n" << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candFromRank.at(jockey)) << "\t" << 1/(candReactance.at(jockey)) << "\t" << ((decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candFromRank.at(jockey))]).get(GRB_DoubleAttr_X)) << "\t" << LHS[rCount] << endl;
					constSatOut << "\nRow: " << "\tColumn(To Node): " << "\tTo Reciprocal Reactance" << "\tTo Angle" << "\tLHS value" << endl;
					LHS[rCount] += (1/(candReactance.at(jockey)))*((decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candToRank.at(jockey))]).get(GRB_DoubleAttr_X));
					constSatOut << "\n" << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candToRank.at(jockey)) << "\t" << 1/(candReactance.at(jockey)) << "\t" << ((decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candToRank.at(jockey))]).get(GRB_DoubleAttr_X)) << "\t" << LHS[rCount] << endl;
					constSatOut << "\nRow: " << "\tColumn(Integer variable): " << "\tLine Building decision" << "\tLHS value" << endl;
					LHS[rCount] += (-(2.5*(candCapacity.at(jockey))))*((decvar[(rCount-countOfScenarios*(2*sharedELineNumber+2*candLineNumber)-scenCounter*candLineNumber)]).get(GRB_DoubleAttr_X));
					constSatOut << rCount << "\t" << (rCount-countOfScenarios*(2*sharedELineNumber+2*candLineNumber)-scenCounter*candLineNumber) << "\t" << ((decvar[(rCount-countOfScenarios*(2*sharedELineNumber+2*candLineNumber)-scenCounter*candLineNumber)]).get(GRB_DoubleAttr_X)) << "\t" << LHS[rCount] << endl;
					constSatOut << LHS[rCount] << ">=" << -2.5*(candCapacity.at(jockey)) << endl;
					++rCount; // Increment the row count to point to the next transmission line object
					++jockey;
					//cout << "\nTest after shared candidate Line Definition lower bound " << diffSerCounter << endl;
				}
				++diffSerCounter;
			}
		}
		clock_t end2 = clock(); // stop the timer
		double elapsed_secs2 = double(end2 - begin) / CLOCKS_PER_SEC; // Calculate the Total Time
		outPutFile << "\nTotal time taken to solve the MILP Line Construction Decision Making Problem instance and retrieve the results = " << elapsed_secs2 << " s " << endl;
		cout << "\nTotal time taken to solve the MILP Line Construction Decision Making Problem instance and retrieve the results = " << elapsed_secs2 << " s " << endl;

		// Close the different output files
		candFlowMWOut.close();
		candLineDecisionOut.close();
		externalAngleOut.close();
		constSatOut.close();
		cout << "\nSimulation Completed.\nResults written on the different output files" << endl;
		delete modelMOMILP; // Free the memory of the GUROBI Problem Model
	}
	outPutFile.close();
	powerGenOut.close();
	//%%calcMILPBounds(); // Calculate the bounds
} // Function MILP() ends

double Marketover::LBMarketoverGUROBI(double LagMultXi[], double LagMultPi[], int totalCandLineNum, int totalSharedNodeNum, GRBEnv* environmentGUROBI) // Function LBMarketover() calculates the lower bound of the Mixed Integer Linear Programming Solver routine by calling GLPK routines for average heat rate objective
{
	// CREATION OF THE MIP SOLVER INSTANCE //
	clock_t begin = clock(); // start the timer
	vector<int>::iterator diffZNIt; // Iterator for diffZoneNodeID
	vector<Node*>::iterator nodeIterator; // Iterator for node objects
	vector<double>::iterator candCapIterator; // Iterator for candidate lines for capacity
	vector<int>::iterator candIterator; // Iterator for candidate lines	
	vector<double>::iterator exsharedCapIterator; // Iterator for shared existing lines for capacity
	vector<int>::iterator exsharedIterator;

	string outSummaryFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputSummaryGUROBIBounds/OutSumBoundMOGUROBI.txt";
	ofstream outPutFile(outSummaryFileName, ios::out); // Create Output File to output the Summary of Results
	if (!outPutFile){
		cerr << "\nCouldn't open the file" << endl;
		exit(1);
	}

        int dimRow = countOfScenarios*(4 * candLineNumber + 2 * sharedELineNumber); // Total number of rows of the A matrix (number of structural constraints of the LP) first term to account for lower and upper line limits & lower and upper definition limits of candidate shared lines, and second term for lower and upper line limits for shared existing lines
        int dimCol = countOfScenarios*nodeNumber+(countOfScenarios+1)*candLineNumber; // Total number of columns of the LP (number of Decision Variables) first term to account for voltage phase angles for inter-zonal lines' nodes, and second term for power flow values and binary integer decision variable values for shared candidate lines
	outPutFile << "\nTotal Number of Structural Constraints (Rows) is: " << dimRow << endl;
	outPutFile << "\nTotal Number of Decision Variables (Columns) is: " << dimCol << endl;
	//cout << "\nTotal Number of Structural Constraints (Rows) is: " << dimRow << endl;
	//cout << "\nTotal Number of Decision Variables (Columns) is: " << dimCol << endl;
	//cout << "\nTotal Number of Shared Nodes is: " << nodeNumber << endl;	
	//cout << "\nTotal Number of shared candidate lines is: " << candLineNumber << endl;
	//cout << "\nTotal Number of shared existing lines is: " << sharedELineNumber << endl;
	GRBModel *modelMOMILP = new GRBModel(*environmentGUROBI);
    	modelMOMILP->set(GRB_StringAttr_ModelName, "assignmentBoundMO");
	GRBVar decvar[dimCol+1];
	double z; // variable to store the objective value
	vector<double> x; // Coefficient matrix entires, objective function, and decision variables

	// SPECIFICATION OF PROBLEM PARAMETERS //
	// Dummy Decision Variable //
	decvar[0] = modelMOMILP->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);

	//Decision Variable Definitions, Bounds, and Objective Function Co-efficients//
	int colCount = 1;

	//Columns corresponding to Shared Candidate Line Flows continuous variables//
	outPutFile << "\nColumns corresponding to Shared Candidate Line Flows continuous variables" << endl;
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candCapIterator = candCapacity.begin(); candCapIterator != candCapacity.end(); ++candCapIterator){
			decvar[colCount] = modelMOMILP->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
			outPutFile << colCount << "\t";
			outPutFile << 0 << endl;
			++colCount;
		}
	}
	outPutFile << "\nTotal number of columns after accounting for Shared Candidate Line Flows continuous variables: " << colCount << endl;

	//Columns corresponding to Shared Candidate Line Construction Decision Relaxed Binary Integer variables//
	outPutFile << "\nColumns corresponding to Shared Candidate Line Construction Decision Relaxed Binary Integer variables" << endl;
	int candInd = 0; // Initialize the counter for indexing the candidate lines
	int diffCandCounter = 0; // counter flag to indicate the first element of the candSerial list
	for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
		if (diffCandCounter > 0) { // Skip the first element, since it's a dummy "0"
			decvar[colCount] = modelMOMILP->addVar(0, 1, 0.0, GRB_CONTINUOUS);
			outPutFile << colCount << "\t";
			outPutFile << LagMultPi[candInd] << "\t" << candInd << endl;
			++candInd;
			++colCount;
		}
		++diffCandCounter;
	}
	outPutFile << "\nTotal number of columns after accounting for Shared Candidate Line Construction Decision Relaxed Binary Integer variables: " << colCount << endl;

	//Columns corresponding to Voltage Phase Angles continuous variables for other zone nodes for shared lines//
	outPutFile << "\nColumns corresponding to Voltage Phase Angles continuous variables for other zone nodes for shared lines" << endl;
	int otherInd = 0; // Initialize the counter for indexing the other-zone nodes
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int columnCount = 1;
		int diffNodeCounter = 0; // counter flag to indicate the first element of the diffZoneNodeID list
		//for (diffZNIt = diffZoneNodeID.begin(); diffZNIt != diffZoneNodeID.end(); ++diffZNIt){
			//cout << " Column count after the shared node is " << columnCount << endl;	
			//++columnCount;
		//}//
		for (diffZNIt = nodeList.begin(); diffZNIt != nodeList.end(); ++diffZNIt){
			if (diffNodeCounter > 0) { // Skip the first element, since it's a dummy "0"
				decvar[colCount] = modelMOMILP->addVar(0, (44/7), 0.0, GRB_CONTINUOUS);
				outPutFile << colCount << "\t";
				outPutFile << LagMultXi[otherInd] << "\t" << otherInd << endl;
				++otherInd;	
				++colCount;
			}
			++diffNodeCounter;
		}
	}
	outPutFile << "\nDecision Variables and Objective Function defined" << endl;
	outPutFile << "\nTotal Number of columns: " << colCount - 1 << endl;
	//cout << "\nDecision Variables and Objective Function defined" << endl;
	//cout << "\nTotal Number of columns: " << colCount - 1 << endl;
	//Setting Objective//
	GRBLinExpr obj = 0.0;
	// Objective Contribution from Dummy Decision Variable //
	obj += 0*(decvar[0]);
	colCount = 1;

	//Columns corresponding to Shared Candidate Line Flows continuous variables//
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candCapIterator = candCapacity.begin(); candCapIterator != candCapacity.end(); ++candCapIterator){
			obj += 0*(decvar[colCount]);
			++colCount;
		}
	}
	//Columns corresponding to Shared Candidate Line Construction Decision Binary Integer variables//
	candInd = 0; // Initialize the counter for indexing the candidate lines
	diffCandCounter = 0; // counter flag to indicate the first element of the candSerial list
	for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
		if (diffCandCounter > 0) { // Skip the first element, since it's a dummy "0"
			obj += (LagMultPi[candInd])*(decvar[colCount]);
			++candInd;
			++colCount;
		}
		++diffCandCounter;
	}

	//Columns corresponding to Voltage Phase Angles continuous variables for other zone nodes for shared lines//
	otherInd = 0; // Initialize the counter for indexing the other-zone nodes
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int columnCount = 1;
		int diffNodeCounter = 0; // counter flag to indicate the first element of the diffZoneNodeID list
		for (diffZNIt = nodeList.begin(); diffZNIt != nodeList.end(); ++diffZNIt){
			if (diffNodeCounter > 0) { // Skip the first element, since it's a dummy "0"
				obj += (LagMultXi[otherInd])*(decvar[colCount]);
				++otherInd;	
				++colCount;
			}
			++diffNodeCounter;
		}
	}
	//cout << "\nDecision Variables and Objective Function defined" << endl;
	//cout << "\nTotal Number of columns: " << colCount - 1 << endl;

	modelMOMILP->setObjective(obj, GRB_MAXIMIZE);
	//cout << " Objective Function and Decision Variables have been defined and the colCount is " << colCount-1 << endl;

	//Row Definitions: Specification of RHS or b vector of b<=Ax<=b//
	GRBLinExpr lhs[dimRow+1];
	//Row Definitions and Bounds Corresponding to Constraints/
	string outPGenFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputPowerGUROBIBounds/OutPowerGenMOGUROBI.txt"; 
	ofstream powerGenOut(outPGenFileName, ios::out);
	if (!powerGenOut){
		cerr << "\nCouldn't open the file" << endl;
		exit(1);
	}
	outPutFile << "\nTesting of interzonal transmission limits" << endl; 
	// Dummy Constraint //
	lhs[0] = 0*(decvar[0]);
	modelMOMILP->addConstr(lhs[0], GRB_EQUAL, 0);


	//Non-Zero entries of A matrix (Constraint/Coefficient matrix entries)//
	int rCount = 1; // Initialize the row count
	
	// Coefficients corresponding to shared existing Line Forward Flow Limit Constraints
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared existing Line Forward Flow Limit Constraints" << endl;
		for (exsharedIterator = SESerial.begin(); exsharedIterator != SESerial.end(); ++exsharedIterator){
			lhs[rCount] = 0;
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				//cout << "Test Message 2" << endl;
				//cout << "SE From Rank" << SEFromRank.at(jockey) << endl;
				//cout << "SE Reactance" << 1/(SEReactance.at(jockey)) << endl;
				lhs[rCount] += (1/(SEReactance.at(jockey)))*(decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEFromRank.at(jockey))]);
				outPutFile << "\n" << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEFromRank.at(jockey)) << "\t" << 1/(SEReactance.at(jockey)) << endl;
				lhs[rCount] += (-1/(SEReactance.at(jockey)))*(decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEToRank.at(jockey))]);
				outPutFile << "\n" << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEToRank.at(jockey)) << "\t" << -1/(SEReactance.at(jockey)) << endl;
				modelMOMILP->addConstr(lhs[rCount] <= (SECapacity.at(jockey)));
				outPutFile << rCount << "\t";
				outPutFile << (SECapacity.at(jockey)) << endl;
				++rCount; // Increment the row count to point to the next transmission line object
				++jockey;
				//cout << "\nTest after shared existing Line Forward Flow Limit Constraint " << diffSerCounter << endl;
			}
			++diffSerCounter;	
		}
	}
	// Coefficients corresponding to shared existing Line Reverse Flow Limit Constraints
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared existing Line Reverse Flow Limit Constraints" << endl;
		for (exsharedIterator = SESerial.begin(); exsharedIterator != SESerial.end(); ++exsharedIterator){
			lhs[rCount] = 0;
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				lhs[rCount] += (1/(SEReactance.at(jockey)))*(decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEFromRank.at(jockey))]);
				outPutFile << "\n" << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEFromRank.at(jockey)) << "\t" << 1/(SEReactance.at(jockey)) << endl;
				lhs[rCount] += (-1/(SEReactance.at(jockey)))*(decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEToRank.at(jockey))]);
				outPutFile << "\n" << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+SEToRank.at(jockey)) << "\t" << -1/(SEReactance.at(jockey)) << endl;
				modelMOMILP->addConstr(lhs[rCount] >= -(SECapacity.at(jockey)));
				outPutFile << rCount << "\t";
				outPutFile << -(SECapacity.at(jockey)) << endl;
				++rCount; // Increment the row count to point to the next transmission line object
				++jockey;
				//cout << "\nTest after shared existing Line Reverse Flow Limit Constraint " << diffSerCounter << endl;
			}
			++diffSerCounter;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for Line Flow limits for shared existing lines: " << rCount << endl;
	
	// Coefficients corresponding to shared candidate Line Forward Flow Limit Constraints
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared candidate Line Forward Flow Limit Constraints" << endl;
		for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
			lhs[rCount] = 0;
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				lhs[rCount] += (1)*(decvar[(rCount-2*countOfScenarios*sharedELineNumber)]);
				outPutFile << "\n" << rCount << "\t" << (rCount-2*countOfScenarios*sharedELineNumber) << "\t" << 1 << endl;
				lhs[rCount] += (-candCapacity.at(jockey))*(decvar[(countOfScenarios*candLineNumber+rCount-2*countOfScenarios*sharedELineNumber-scenCounter*candLineNumber)]);
				outPutFile << "\n" << rCount << "\t" << (countOfScenarios*candLineNumber+rCount-2*countOfScenarios*sharedELineNumber-scenCounter*candLineNumber) << "\t" << -candCapacity.at(jockey) << endl;
				modelMOMILP->addConstr(lhs[rCount] <= 0.0);
				outPutFile << rCount << "\t";
				outPutFile << 0.0 << endl;
				++rCount;
				++jockey;
				//cout << "\nTest after shared candidate Line Forward Flow Limit Constraint " << diffSerCounter << endl;
			}
			++diffSerCounter;
		}
	}
	// Coefficients corresponding to shared candidate Line Reverse Flow Limit Constraints
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared candidate Line Reverse Flow Limit Constraints" << endl;
		for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
			lhs[rCount] = 0;
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				lhs[rCount] += (1)*(decvar[(rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber))]);
				outPutFile << "\n" << rCount << "\t" << (rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber)) << "\t" << 1 << endl;
				lhs[rCount] += (candCapacity.at(jockey))*(decvar[(rCount-countOfScenarios*(2*sharedELineNumber)-scenCounter*candLineNumber)]);
				outPutFile << "\n" << rCount << "\t" << (rCount-countOfScenarios*(2*sharedELineNumber)-scenCounter*candLineNumber) << "\t" << candCapacity.at(jockey) << endl;
				modelMOMILP->addConstr(lhs[rCount] >= 0.0);
				outPutFile << rCount << "\t";
				outPutFile << 0.0 << endl;
				++rCount;
				++jockey;
				//cout << "\nTest after shared candidate Line Reverse Flow Limit Constraint " << diffSerCounter << endl;
			}
			++diffSerCounter;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for Line Flow limits for shared candidate lines: " << rCount << endl;
	
	// Coefficients corresponding to shared candidate Line Definition upper bound
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared candidate Line Definition upper bound" << endl;
		for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
			lhs[rCount] = 0;
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				lhs[rCount] += (1)*(decvar[(rCount-countOfScenarios*(2*sharedELineNumber+2*candLineNumber))]);
				outPutFile << rCount << "\t" << (rCount-countOfScenarios*(2*sharedELineNumber+2*candLineNumber)) << " Power term for candidate line " << diffSerCounter << "-th candidate line" << endl;
				lhs[rCount] += (-1/(candReactance.at(jockey)))*(decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candFromRank.at(jockey))]);
				outPutFile << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candFromRank.at(jockey)) << "\t" << -1/(candReactance.at(jockey)) << endl;
				lhs[rCount] += (1/(candReactance.at(jockey)))*(decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candToRank.at(jockey))]);
				outPutFile << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candToRank.at(jockey)) << "\t" << 1/(candReactance.at(jockey)) << endl;
				lhs[rCount] += (2.5*(candCapacity.at(jockey)))*(decvar[(rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber)-scenCounter*candLineNumber)]);
				outPutFile << rCount << "\t" << (rCount-countOfScenarios*(2*sharedELineNumber+candLineNumber)-scenCounter*candLineNumber) << "\t" << BIGM << endl;
				modelMOMILP->addConstr(lhs[rCount] <= 2.5*(candCapacity.at(jockey)));
				outPutFile << rCount << "\t";
				outPutFile << BIGM << endl;
				++rCount; // Increment the row count to point to the next transmission line object
				++jockey;
				//cout << "\nTest after shared candidate Line Definition upper bound " << diffSerCounter << endl;
			}
			++diffSerCounter;
		}
	}
	// Coefficients corresponding to shared candidate Line Definition lower bound
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int jockey = 0;
		int diffSerCounter = 0;
		outPutFile << "\nCoefficients corresponding to shared candidate Line Definition lower bound" << endl;
		for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
			lhs[rCount] = 0;
			if (diffSerCounter > 0) { // Skip the first element, since it's a dummy "0"
				lhs[rCount] += (1)*(decvar[(rCount-countOfScenarios*(2*sharedELineNumber+3*candLineNumber))]);
				outPutFile << rCount << "\t" << (rCount-countOfScenarios*(2*sharedELineNumber+3*candLineNumber)) << " Power term for candidate line " << diffSerCounter << "-th candidate line" << endl;
				lhs[rCount] += (-1/(candReactance.at(jockey)))*(decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candFromRank.at(jockey))]);
				outPutFile << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candFromRank.at(jockey)) << "\t" << -1/(candReactance.at(jockey)) << endl;
				lhs[rCount] += (1/(candReactance.at(jockey)))*(decvar[((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candToRank.at(jockey))]);
				outPutFile << rCount << "\t" << ((countOfScenarios+1)*candLineNumber + scenCounter*nodeNumber+candToRank.at(jockey)) << "\t" << 1/(candReactance.at(jockey)) << endl;
				lhs[rCount] += (-(2.5*(candCapacity.at(jockey))))*(decvar[(rCount-countOfScenarios*(2*sharedELineNumber+2*candLineNumber)-scenCounter*candLineNumber)]);
				outPutFile << rCount << "\t" << (rCount-countOfScenarios*(2*sharedELineNumber+2*candLineNumber)-scenCounter*candLineNumber) << "\t" << -BIGM << endl;
				modelMOMILP->addConstr(lhs[rCount] >= -(2.5*(candCapacity.at(jockey))));
				outPutFile << rCount << "\t";
				outPutFile << -BIGM << endl;
				++rCount; // Increment the row count to point to the next transmission line object
				++jockey;
				//cout << "\nTest after shared candidate Line Definition lower bound " << diffSerCounter << endl;
			}
			++diffSerCounter;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for shared candidate lines flow definition: " << rCount << endl;
	outPutFile << "\nConstraint bounds (rows) Specified" << endl;
	outPutFile << "\nTotal number of rows: " << rCount - 1 << endl;
	//cout << "\nConstraint bounds (rows) Specified" << endl;
	//cout << "\nTotal number of rows: " << rowCount - 1 << endl;

	outPutFile << "\nCoefficient Matrix specified" << endl;
	clock_t end1 = clock(); // stop the timer
	double elapsed_secs1 = double(end1 - begin) / CLOCKS_PER_SEC; // Calculate the time required to populate the constraint matrix and objective coefficients
	outPutFile << "\nTotal time taken to define the rows, columns, objective and populate the coefficient matrix = " << elapsed_secs1 << " s " << endl;
	// RUN THE OPTIMIZATION SIMULATION ALGORITHM //
	cout << "\nSimulation in Progress. Wait !!! ....." << endl;
	modelMOMILP->optimize(); // Solves the optimization problem
	int stat = modelMOMILP->get(GRB_IntAttr_Status); // Outputs the solution status of the problem 

	// DISPLAY THE SOLUTION DETAILS //
	if (stat == GRB_INFEASIBLE){
		outPutFile << "\nThe solution to the problem is INFEASIBLE." << endl;
		cout << "\nThe solution to the problem is INFEASIBLE." << endl;
		delete modelMOMILP; // Free the memory of the GUROBI Problem Model
		return 0;
	} else if (stat == GRB_INF_OR_UNBD) {
		outPutFile << "\nNO FEASIBLE or BOUNDED solution to the problem exists." << endl;
		cout << "\nNO FEASIBLE or BOUNDED solution to the problem exists." << endl;
		delete modelMOMILP; // Free the memory of the GUROBI Problem Model
		return 0;
	} else if (stat == GRB_UNBOUNDED) {
		outPutFile << "\nThe solution to the problem is UNBOUNDED." << endl;
		cout << "\nThe solution to the problem is UNBOUNDED." << endl;
		delete modelMOMILP; // Free the memory of the GUROBI Problem Model
		return 0;
	} else if (stat == GRB_OPTIMAL) {
		outPutFile << "\nThe solution to the problem is OPTIMAL." << endl;
		cout << "\nThe solution to the problem is OPTIMAL." << endl;
	}

	//Get the Optimal Objective Value results//
	z = modelMOMILP->get(GRB_DoubleAttr_ObjVal);

	// Open separate output files for writing results of different variables
	string outCandFlowFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/candidateLinesGUROBIBounds/candFlowMWBoundMOGUROBI.txt";
	string outCandDecFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/candidateLinesGUROBIBounds/candLineDecisionBoundMOGUROBI.txt";
	string outExtAngFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputAnglesGUROBIBounds/externalAngleBoundMOGUROBI.txt";
	ofstream candFlowMWOut(outCandFlowFileName, ios::out); //switchOnOut
	ofstream candLineDecisionOut(outCandDecFileName, ios::out); //switchOffOut
	ofstream externalAngleOut(outExtAngFileName, ios::out);
	outPutFile << "\nThe Optimal Objective value (Line Building Decision cost) is: " << -z << endl;
	powerGenOut << "\nThe Optimal Objective value (Line Building Decision cost) is: " << -z << endl;
	cout << "\nThe Optimal Objective value (Line Building Decision cost) is: " << -z << endl;
	x.push_back(0); // Initialize the decision Variable vector

	// Display Shared Candidate lines' Power Flow variables
	candFlowMWOut << "\n****************** SHARED CANDIDATE LINES POWER FLOW VALUES *********************" << endl;
	candFlowMWOut << "CANDIDATE LINE ID" << "\t" << "MW FLOW VALUE" << "\n";
	int arrayInd = 1;
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int diffCandCounter2 = 0; // counter flag to indicate the first element of the candSerial list
		for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
			if (diffCandCounter2 > 0) { // Skip the first element, since it's a dummy "0"
				x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
				candFlowMWOut << (*candIterator) << "\t" << ((decvar[arrayInd]).get(GRB_DoubleAttr_X))*100 << " MW" << endl;
				++arrayInd;
			}
			++diffCandCounter2;
		}
	}
	candFlowMWOut << "Finished writing Shared Candidate lines' Power Flow variables" << endl;

	// Display Shared Candidate lines' Construction Decisions
	candLineDecisionOut << "\n****************** SHARED CANDIDATE LINES CONSTRUCTION DECISIONS *********************" << endl;
	candLineDecisionOut << "CANDIDATE LINE ID" << "\t" << "CONSTRUCTION DECISIONS" << "\n";
	int candInd1 = 0; // Initialize the counter for indexing the candidate lines
	int diffCandCounter1 = 0; // counter flag to indicate the first element of the candSerial list
	for (candIterator = candSerial.begin(); candIterator != candSerial.end(); ++candIterator){
		if (diffCandCounter1 > 0) { // Skip the first element, since it's a dummy "0"
			x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			candLineDecisionOut << (*candIterator) << "\t" << ((decvar[arrayInd]).get(GRB_DoubleAttr_X)) << endl;
			++arrayInd;
			++candInd1;
		}
		++diffCandCounter1;
	}
	candLineDecisionOut << "Finished writing Shared Candidate lines' Construction decisions" << endl;

	// Display shared node angles
	externalAngleOut << "\n****************** OUTER ZONAL NODE VOLTAGE PHASE ANGLE VALUES *********************" << endl;
	externalAngleOut << "EXTERNAL NODE ID" << "\t" << "VOLTAGE PHASE ANGLE" << "\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int otherInd1 = 0; // Initialize the counter for indexing the other-zone nodes
		int dummyStart = 0;
		for (diffZNIt = nodeList.begin(); diffZNIt != nodeList.end(); ++diffZNIt){
			if (dummyStart > 0) { // Skip the dummy element "0" at the beginning
				x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
				externalAngleOut << (*diffZNIt) << "\t" << ((decvar[arrayInd]).get(GRB_DoubleAttr_X)) << endl;
				++arrayInd;
				++otherInd1;
			}
			++dummyStart;
		}
	}
	externalAngleOut << "Finished writing shared node voltage phase angle values" << endl;

	delete modelMOMILP; // Free the memory of the GUROBI Problem Model
	clock_t end2 = clock(); // stop the timer
	double elapsed_secs2 = double(end2 - begin) / CLOCKS_PER_SEC; // Calculate the Total Time
	outPutFile << "\nTotal time taken to solve the MILP Line Construction Decision Making Problem instance and retrieve the results = " << elapsed_secs2 << " s " << endl;
	cout << "\nTotal time taken to solve the MILP Line Construction Decision Making Problem instance and retrieve the results = " << elapsed_secs2 << " s " << endl;

	// Close the different output files
	outPutFile.close();
	powerGenOut.close();
	candFlowMWOut.close();
	candLineDecisionOut.close();
	externalAngleOut.close();
	cout << "\nSimulation Completed.\nResults written on the different output files" << endl;
	//%%calcMILPBounds(); // Calculate the bounds
} // Function MILP() ends

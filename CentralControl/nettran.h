#ifndef NETTRAN_H
#define NETTRAN_H
#include <vector>
#include "gurobi_c++.h" // includes definition of the GUROBI solver header file
#include "powergenerator.h" // includes defintion of Generator class
#include "transl.h" // includes definition of Transmission line class
#include "load.h" // includes definition of load class
#include "node.h" // includes definition of node class
#include "sharedLine.h" // includes definition of existing shared lines between two different zones
#include "candidateLine.h" // includes definition of candidate lines shared between two different zones
#include "intcandidateLine.h" // includes the definition of intra zonal candidate lines

using namespace std;

// Nettran class definition
class Nettran {
private:
	int nodeNumber; // Number of Nodes
	vector<int> nodeNumVector; // Vector of subnet number of nodes
	vector<int> SELineSerList; // List of the global ID of SE lines
	vector<int> candLineSerList; // List of the global ID of shared candidate lines
	int univNodeNum = 0; // Initialize the total number of nodes of all the zones
	int univGenNum = 0; // Initialize the universal generator number
	int univTranNum = 0; // Initialize the universal transmission line number
	int univLoadNum = 0; // Initialize the universal load number
	int univIntCandNum = 0; // Initialize the universal intra zonal candidate line number
	int univSELineNum = 0; // Initialize the universal shared existing line number
	int univCandLineNum = 0; // Initialize the universal shared candidate line number
	int genNumber; // Number of Generators
	int tranNumber; // Number of Transmission lines
	int loadNumber; // Number of Loads
	int zonalCount; // Total number of different zones
	int countOfScenarios; // Number of total random scenarios
	int sharedELines; // Number of shared existing transmission lines
	int sharedCLines; // Number of shared candidate transmission lines
	int internalCLines; // NUmber of internal candidate transmission lines
	int simMode; // Mode of the simulation: Average Heat Rate, Piece-wise Linear, Polynomial (Cubic) Convex cost curve
	char genFile[ 100 ]; // String for storing the name of the generator file
	char tranFile[ 100 ]; // String for storing the name of the transmission line file
	char loadFile[ 100 ]; // String for storing the name of the load file
	char netFile[ 100 ]; // String for storing the name of the network file
	char sharedLineFile[ 100 ]; // String for storing the name of the shared existing lines file
	char candLineFile[ 100 ]; // String for storing the name of the candidate lines file	
	char intCandLineFile[ 100 ]; // String for storing the name of the internal/intra-zonal candidate lines file
	vector<Powergenerator*> genObject; // Set of Powergenerator objects
	vector<Load*> loadObject; // set of load objects
	vector<transmissionLine*> translObject; // set of transmission line objects
	vector<Node*> nodeObject; // set of node objects 
	vector<SELine*> SELineObject; // set of existing shared lines objects 
	vector<candLine*> candLineObject; // set of shared candidate lines objects
	vector<intCandLine*> intCandLineObject; // set of internal candidate line objects	
	vector<double> probability; // Probability of the different scenarios

public:
	// Functions for building the Network
	Nettran(string[], int, int); // Constructor
	~Nettran(); // destructor
	//void setNetwork(int); // Set function to set up network variables and validate the range of values
	// Functions for carrying out the Actual Simulation
	double MILPAvgHRGUROBI(GRBEnv*); // Calls GLPK solver to solve the MILP problem for Horizontal Coordination for Average Heat Rate Objective function
	void MILPPiecewiseLin(void); // Calls GLPK solver to solve the MILP problem for Horizontal Coordination for Piecewise Linear Objective function 
	void MILPPolynomial(); // Calls GLPK solver to solve the MILP problem for Horizontal Coordination for Polynomial Convex Cubic Objective function 
	void assignProb(); 
}; // end of class definition

#endif
//NETTRAN_H

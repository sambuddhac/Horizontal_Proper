#ifndef MARKETOVER_H
#define MARKETOVER_H
#include "gurobi_c++.h" // includes definition of the GUROBI solver header file
#include <vector>
#include "nettran.h"

using namespace std;

class Nettran; // Forward declaration of the Nettran class

// Marketover class definition
class Marketover {
private:
	int candLineNumber; // Number of shared candidate lines
	int sharedELineNumber; // Number of shared existing Transmission lines
	int zoneNumber; // Number of zones
	int nodeNumber; // Number of shared lines' nodes 
	int lpSolveAlgo; // Simplex for 1 and IPM for 2
	int countOfScenarios; // Number of total random scenarios
	vector<Nettran*> subNetVector; // Vector of zonal subnetworks
	vector<int> zoneList; // List of zones between which exsiting and candidate shared lines exist
	vector<int> nodeList; // List of nodes at the ends of the existing and candidate shared lines
	vector<int> sharedGlobalList; // List of global rankings of shared nodes
	vector<int> SESerial; // Serial list of shared existing lines
	vector<int> SEFromRank; // List of rank of from nodes of shared existing lines
	vector<int> SEToRank; // List of rank of to nodes of shared existing lines
	vector<double> SEReactance; // List of reactances of shared existing lines
	vector<double> SECapacity; // List of line flow limits of shared existing lines
	vector<int> candSerial; // Serial list of shared candidate lines
	vector<int> candFromRank; // List of rank of from nodes of shared candidate lines
	vector<int> candToRank; // List of rank of to nodes of shared candidate lines
	vector<double> candReactance; // List of reactances of shared candidate lines
	vector<double> candCapacity; // List of line flow limits of shared candidate lines
	vector<int> lineInterDecision; // Interim decisions of the different zones for building lines
	vector<int> lineInterDecVec; // Interim decisions of the different zones for building lines
	vector<int> lineDecIndex; // Global rank of the variable for integer decision variable or global rank of the candidate line
	vector<int> zonalIndVectorInt; // vector of zonal indices for tracking zonal decisions of integer variables
	vector<double> phaseAngleDecision[100]; // Interim decisions of the different zones for the phase angle values
	vector<double> phAngDecVec; // Interim decisions of the different zones for the phase angle values
	vector<int> scenVector; // Vector of scenario counts or scenario indices
	vector<int> angleDecIndex[100]; // Global rank of the variable for angle decision variable or global rank of the shared nodes
	vector<int> zonalIndVectorCont[100]; // vector of zonal indices for tracking zonal decisions of continuous angle variables 
	vector<double> interimContDecVar; // Interim continuous decision variable of the marketoverseer
	vector<double> interimContDecVarPrev; // Interim continuous decision variable of the marketoverseer from previous iteration
	vector<int> interimIntDecVar; // Interim integer decision variables of the marketOverseer
	vector<int> interimIntDecVarPrev; // Interim integer decision variables of the marketOverseer from previous iteration
	vector<int> constructedRanks; // Vector of ranks of constructed lines
	int compareBasis; // While updating the Lagrange multipliers, 0 when the current zonal updates are compared to previous iteration MO update;else, 1

public:
	// Functions for building the Market Overseer
	Marketover(int, int, int, int, vector <Nettran*>*, vector <int>*, vector <int>*, vector <int>*, vector <int>*, vector <int>*, vector <int>*, vector <double>*, vector <double>*, vector <int>*, vector <int>*, vector <int>*, vector <double>*, vector <double>*, int, int, int); // Constructor
	~Marketover(); // destructor
	void populateLineDec(int, int, int); // Method to pass the intermediate message for the zonal line building decision to the MO
	void populateAngleDec(double, int, int, int); // Method to pass the intermediate message for the zonal node angle decision to the MO
	// Functions for carrying out the Actual Simulation
	void MILPMarketover(double [], double [], int, int); // Calls GLPK solver to solve the MILP problem for attaining consensus
	void MILPMarketoverGUROBI(double [], double [], int, int, GRBEnv*); // Calls GLPK solver to solve the MILP problem for attaining consensus
	double LBMarketover(double [], double [], int, int); // Function LBMarketover() calculates the lower bound of the Mixed Integer Linear Programming Solver routine by calling GLPK routines for average heat rate objective
	double LBMarketoverGUROBI(double [], double [], int, int, GRBEnv*); // Function LBMarketover() calculates the lower bound of the Mixed Integer Linear Programming Solver routine by calling GLPK routines for average heat rate objective
	double rewardPenaltyCont(double, int, int); // Updates the rewards/penalties for bus voltage phase angle consensus
	double rewardPenaltyInteger(double, int, int); // Updates the rewards/penalties for line building decision consensus
	double getGlobalUpper(double [], double [], double [], int); // Gets the global upper bound for the investment coordination problem
	double getGlobalLower(double [], int); // Gets the global lower bound for the investment coordination problem
	void clearVectors(); // Clears the different interim vectors for making them ready for the next iteration
	double getGlobalConsensus(); // Gets the global consensus for the investment coordination problem
	void finDecLineConstr(); // Final decisions on the construction stauses of candidate lines taking into account the decisions of different zones
	int scanBuiltLinesList(int); // Scans the list of built candidate lines for the input argument of global rank to check if this line is built
	void bufferintermediateDecision(int); // Buffering the previous iterations' MO decision values, (for comparison basis equal to 0)
	void clearDelayedVectors(); // Clears the different interim vectors only buffer vectors
	void rewardPenaltyUpdate(double [], double [], int, int); // Update the Lagrange Multipliers
}; // end of class definition

#endif
//MARKETOVER_H

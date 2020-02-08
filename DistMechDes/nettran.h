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
#include "marketOverseer.h" // includes definition of the market overseer class for passing messages

using namespace std;

class Marketover; // Forward declaration of the Marketover class

// Nettran class definition
class Nettran {
private:
	int nodeNumber; // Number of Nodes
	int otherNodeCount; // Number of the other-zone nodes for shared existing and shared candidate lines
	int existingOtherZoneNodeNum=0; // The number of other zone nodes that are the ends of the SE lines and the built shared cand lines
	int genNumber; // Number of Generators
	int tranNumber; // Number of Transmission lines
	int loadNumber; // Number of Loads
	int zonalCount; // Total number of different zones
	int countOfScenarios; // Number of total random scenarios
	int sharedELines; // Number of shared existing transmission lines
	int sharedCLines; // Number of shared candidate transmission lines
	int internalCLines; // NUmber of internal candidate transmission lines
	int realizedCLines; // Number of shared candidate transmission lines that are actually built
	int realizedIntCLines; // Number of internal candidate transmission lines that are actually built
	int simMode; // Mode of the simulation: Average Heat Rate, Piece-wise Linear, Polynomial (Cubic) Convex cost curve
	int zonalIndex; // Zone Index
	int lpSolveAlgo; // Simplex for 1 and IPM for 2
	int containsFlag = 0; // Indicates whether a particular node in a particular connected zone has been accounted for or not; 0 for no, 1 for yes
	int containsFlagGlob = 0; // Indicates whether a particular rank has been accounted for or not; 0 for no, 1 for yes
	char genFile[ 100 ]; // String for storing the name of the generator file
	char tranFile[ 100 ]; // String for storing the name of the transmission line file
	char loadFile[ 100 ]; // String for storing the name of the load file
	char netFile[ 100 ]; // String for storing the name of the network file
	char sharedLineFile[ 100 ]; // String for storing the name of the shared existing lines file
	char candLineFile[ 100 ]; // String for storing the name of the candidate lines file	
	char intCandLineFile[ 100 ]; // String for storing the name of the internal/intra-zonal candidate lines file
	vector<int> diffZoneNodeID; // List of node IDs of different zones that are connected via existing or candidate shared lines
	vector<int> diffZoneID; // List of zone IDs that are connected via existing or candidate shared lines for each line
	vector<int> globalRankDiffNode; // List of global rankings of the other zone connected nodes
	vector<int> diffZoneNodeExistingID; // list of external-zone existing node ID's for only SE lines and constructed shared candidate lines
	vector<int> diffZoneExistingID; // list of external-zone existing zone ID's for only SE lines and constructed shared candidate lines
	vector<int> globalExistingRank; // list of global ranking of external zone nodes for only SE lines and shared candidate lines
	vector<int>::iterator globalIterator; // Iterator for the vector globalRankDiffNode
	vector<Powergenerator*> genObject; // Set of Powergenerator objects
	vector<Load*> loadObject; // set of load objects
	vector<transmissionLine*> translObject; // set of transmission line objects
	vector<Node*> nodeObject; // set of node objects 
	vector<SELine*> SELineObject; // set of existing shared lines objects 
	vector<candLine*> candLineObject; // set of shared candidate lines objects
	vector<intCandLine*> intCandLineObject; // set of internal candidate line objects	
	vector<candLine*> realCandLine; // set of realized shared candidate lines objects
	vector<intCandLine*> realIntCandLineObject; // set of realized internal candidate line objects
	vector<int>::iterator otherZoneIter; // Iterator for sharedZoneList
	vector<int>::iterator otherZoneNodeIter; // Iterator for sharedZoneNodeList
	vector<SELine*>::iterator sharedELineIt; // Iterator for sharedExisting line List
	vector<SELine*>::iterator sharedELineFromIt; // Iterator for sharedExisting line List for From node
	vector<SELine*>::iterator sharedELineFZoneIt; // Iterator for sharedExisting line List for From zone
	vector<SELine*>::iterator sharedELineToIt; // Iterator for sharedExisting line List for To node
	vector<SELine*>::iterator sharedELineTZoneIt; // Iterator for sharedExisting line List for To zone
	vector<SELine*>::iterator sharedELineReactIt; // Iterator for sharedExisting line List for Reactance
	vector<SELine*>::iterator sharedELineCapIt; // Iterator for sharedExisting line List for MW flow limit
	vector<candLine*>::iterator sharedCandLineIt; // Iterator for sharedCandidate line List
	vector<candLine*>::iterator sharedCandLineFromIt; // Iterator for sharedCandidate line List for From node
	vector<candLine*>::iterator sharedCandLineFZoneIt; // Iterator for sharedCandidate line List for From zone
	vector<candLine*>::iterator sharedCandLineToIt; // Iterator for sharedCandidate line List for To node
	vector<candLine*>::iterator sharedCandLineTZoneIt; // Iterator for sharedCandidate line List for To zone
	vector<candLine*>::iterator sharedCandLineReactIt; // Iterator for sharedCandidate line List for Reactance
	vector<candLine*>::iterator sharedCandLineCapIt; // Iterator for sharedCandidate line List for MW flow limit
	vector<double> thetaBuffer; // List of intermediate decision variable values of the APP iterate for shared nodes
	vector<double> diffBuffer; // List of the intermediate disagreements between the shared variable values among different zones
	vector<int> globRankBuffer; // List of global ranks of the shared nodes
	vector<double> probability; // Probability of the different scenarios

public:
	// Functions for building the Network
	Nettran(string[], int, int, int, int); // Constructor
	~Nettran(); // destructor
	//void setNetwork(int); // Set function to set up network variables and validate the range of values
	// Functions for carrying out the Actual Simulation
	double MILPAvgHR(Marketover &coordInstance, double [], double [], int, int); // Calls GLPK solver to solve the MILP problem for Horizontal Coordination for Average Heat Rate Objective function
	double MILPAvgHRGUROBI(Marketover &coordInstance, double [], double [], int, int, GRBEnv*); // Calls GLPK solver to solve the MILP problem for Horizontal Coordination for Average Heat Rate Objective function
	double APPQPAvgHR(Marketover &coordInstance, double [], int, int, GRBEnv*, int); // Calls the CVXGEN custom solver object and solver method to solve the problem of determining the values of the continuous variables
	void MILPPiecewiseLin(void); // Calls GLPK solver to solve the MILP problem for Horizontal Coordination for Piecewise Linear Objective function 
	void MILPPolynomial(); // Calls GLPK solver to solve the MILP problem for Horizontal Coordination for Polynomial Convex Cubic Objective function 
	double calcMILPBounds(double [], double [], int, int); // Calls GLPK solver to calculate the bounds on optimal value, fixing the RHS of non-anticipativity after each iteration//%%
	double calcMILPBoundsGUROBI(double [], double [], int, int, GRBEnv*); // Calls GLPK solver to calculate the bounds on optimal value, fixing the RHS of non-anticipativity after each iteration//%%
	int getConnZone(int); // returns the zone indices from the vector, diffZoneID
	int getConnNode(int); // returns the node indices from the vector, diffZoneNodeID
	int getSESerial(int); // returns the pointer to the base of the vector, SELineObject
	int getSEFromNode(int); // returns the from node indices from the vector, SELineObject
	int getSEFromZone(int); // returns the from zone indices from the vector, SELineObject
	int getSEToNode(int); // returns the to node indices from the vector, SELineObject
	int getSEToZone(int); // returns the to zone indices from the vector, SELineObject
	double getSEReactance(int); // returns the line reactance from the vector, SELineObject
	double getSECapacity(int); // returns the line flow capacity from the vector, SELineObject
	void setSEFromRank(int, int); // sets the from rank for the internal zone node end of the shared existing line
	void setSEToRank(int, int); // sets the to rank for the internal zone node end of the shared existing line
	void setSEFromRankConn(int, int); // populates the vector of the global ranks of all the external from nodes of SE lines
	void setSEToRankConn(int, int); // populates the vector of the global ranks of all the external to nodes of SE lines
	int getCandSerial(int); // returns the pointer to the base of the vector, candLineObject
	int getCandFromNode(int); // returns the from node indices from the vector, candLineObject
	int getCandFromZone(int); // returns the from zone indices from the vector, candLineObject
	int getCandToNode(int); // returns the to node indices from the vector, candLineObject
	int getCandToZone(int); // returns the to zone indices from the vector, candLineObject
	double getCandReactance(int); // returns the line reactance from the vector, candLineObject
	double getCandCapacity(int); // returns the line flow capacity from the vector, candLineObject
	void setCandFromRank(int, int); // sets the from rank for the internal zone node end of the shared candidate line
	void setCandToRank(int, int); // sets the to rank for the internal zone node end of the shared candidate line
	void setCandFromRankConn(int, int); // populates the vector of the global ranks of all the external from nodes of cand lines
	void setCandToRankConn(int, int); // populates the vector of the global ranks of all the external to nodes of cand lines
	void assignCandGlobalRank(int, int); // assigns the global rank of the shared candidate line
	void setRealizedCLines(Marketover &coordInstance); // Assigns the realizedCLine variable, the number of candidate lines that are actually built
	int returnMultiplicity(); // Returns the total multiplicity of the shared nodes
	int getNumberOfScenarios(); // Returns the number of scenarios
	vector<double> getZonalDecision(); // Returns the intermediate decision variable values from APP
	vector<int> getZonalRanks(); // Returns the global ranks of the shared decision variables from APP
	void TestBuiltExternalNodes();
}; // end of class definition

#endif
//NETTRAN_H

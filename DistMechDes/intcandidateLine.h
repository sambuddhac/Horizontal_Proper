// candLine class definition.
// Member functions defined in intcandidateLine.cpp
#ifndef INTCANDLINE_H
#define INTCANDLINE_H

class Node; // forward declaration

class intCandLine {

public:
	intCandLine( int, Node *, Node *, double, double, double, int, double, int ); // constructor
	~intCandLine(); // destructor
	int getTranslID(); // returns the ID of the transmission line
	int getTranslNodeID1(); // returns ID number of from node end to which the transmission line is connected
	int getTranslNodeID2(); // returns ID number of to node end to which the transmission line is connected
	double getFlowLimit(); // Gets the value of power flow line limit
	double getReactance(); // Returns the values of the reactance
	void setTranData(double, int, double); 
	double getInvestCost(); // Returns the total investment cost
	int returnPresAbsStatus(); // Returns the construction status of the candidate line
	void setPresAbsStatus(); // Sets the construction status of the candidate line
	void modifyNodeReact(); // function to modify the nodal connected reactance, if the candidate line is actually built

private:
	int translID; // transmissionLine object id number
	double ptMax; // Maximum MW transfer capability
	double reacT; // Resistance and Reactance of the transmission line
	Node *connNodetPtr1, *connNodetPtr2; // pointers to the from and to node objects
	double capitalCost; // Capital cost for the construction of the line
	int lifeYears; // Type of device; 1 for Generating device, 0 otherwise
	int statusOfConstruction; // Default value is 0; 1 if constructed
	double rateInterest; // Rate of Interest
	int resFromStageI; // Result from Stage-I regarding the presence or absence of the candidate transmission line

}; //end class intCandLine

#endif // INTCANDLINE_H
	

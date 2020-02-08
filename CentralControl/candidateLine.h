// candLine class definition.
// Member functions defined in candidateLine.cpp
#ifndef CANDLINE_H
#define CANDLINE_H

class Node; // forward declaration

class candLine {

public:
	candLine( int, int, Node *, Node *, double, double, double, int, double, int ); // constructor
	~candLine(); // destructor
	int getTranslID(); // returns the ID of the transmission line
	int getFromNodeID(); // returns ID number of from node
	int getToNodeID(); // returns ID number of to node
	double getFlowLimit(); // Gets the value of power flow line limit
	double getReactance(); // Returns the values of the reactance
	void setTranData(double, int, double); 
	double getInvestCost(); // Returns the total investment cost
	int returnPresAbsStatus(); // Returns the construction status of the candidate line

private:
	int translID; // transmissionLine object id number
	int localIndex; // local Serial number in the present zone
	double ptMax; // Maximum MW transfer capability
	double reacT; // Resistance and Reactance of the transmission line
	Node *connNodetPtr1; // pointers to the from node object of the shared candidate transmission line
	Node *connNodetPtr2; // pointers to the to node object of the shared candidate transmission line
	double capitalCost; // Capital cost for the construction of the line
	int lifeYears; // Type of device; 1 for Generating device, 0 otherwise
	double rateInterest; // Rate of Interest
	int status; // Presence or absence

}; //end class candLine

#endif // CANDLINE_H
	

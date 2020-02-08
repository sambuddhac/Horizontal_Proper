// transmissionLine class definition.
// Member functions defined in transl.cpp
#ifndef TRANSL_H
#define TRANSL_H

class Node; // forward declaration

class transmissionLine {

public:
	transmissionLine( int, Node *, Node *, double, double ); // constructor
	~transmissionLine(); // destructor
	int getTranslID(); // returns the ID of the transmission line
	int getTranslNodeID1(); // returns ID number of node at end-1 to which the transmission line is connected
	int getTranslNodeID2(); // returns ID number of node at end-2 to which the transmission line is connected
	double getFlowLimit(); // Gets the value of power flow line limit
	double getReactance(); // Gets the reactance of the transmission line

private:
	int translID; // transmissionLine object id number
	double ptMax; // Maximum MW transfer capability
	int deviceNature; // Type of device; 1 for Generating device, 0 otherwise
	double reacT; // Resistance and Reactance of the transmission line
	Node *connNodet1Ptr, *connNodet2Ptr; // pointers to the node objects at the two ends of the transmission line

}; //end class transmissionLine

#endif // TRANSL_H
	

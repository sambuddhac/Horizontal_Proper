// SELine class definition.
// Member functions defined in sharedLine.cpp
#ifndef SELINE_H
#define SELINE_H

class Node; // forward declaration

class SELine {

public:
	SELine( int, int, Node *, Node *, double, double ); // constructor
	~SELine(); // destructor
	int getTranslID(); // returns the ID of the transmission line
	int getFromNodeID(); // returns ID number of from node
	int getToNodeID(); // returns ID number of to node
	double getFlowLimit(); // Gets the value of power flow line limit
	double getReactance(); // returns the reactance value

private:
	int translID; // transmissionLine object id number or global serial number
	int localIndex; // local Serial number in the zone
	double ptMax; // Maximum MW transfer capability
	double reacT; // Resistance and Reactance of the transmission line
	Node *connNodetPtr1; // pointers to the node objects at the from end of the transmission line
	Node *connNodetPtr2; // pointers to the node objects at the to end of the transmission line
	
}; //end class SELine

#endif // SELINE_H
	

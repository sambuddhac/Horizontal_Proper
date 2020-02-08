// SELine class definition.
// Member functions defined in sharedLine.cpp
#ifndef SELINE_H
#define SELINE_H

class Node; // forward declaration

class SELine {

public:
	SELine( int, int, Node *, int, int, int, int, int, double, double ); // constructor
	~SELine(); // destructor
	int getTranslID(); // returns the ID of the transmission line
	int getIntlNodeID(); // returns ID number of intra-zonal node end to which the transmission line is connected
	int getIntlZoneID(); // returns ID number of intra-zonal node end zone to which the transmission line is connected
	int getExtNodeID(); // returns ID number of outer-zonal node end to which the transmission line is connected
	int getExtZoneID(); // returns ID number of outer-zonal node end zone to which the transmission line is connected
	int getExtNodeRank(); // returns serial rank in the list of nodes of the outer end node, to which the transmission line is connected
	int getExtNodeGlobalRank(); // returns global rank in the list of nodes of the outer end node, to which the transmission line is connected
	double getFlowLimit(); // Gets the value of power flow line limit
	int getFlowDir(); // returns the value of the direction flag indicating whether the intra-zonal node end of the line is from (+1) or to (-1) end
	double getReactance(); // returns the reactance value
	void outerNodeIndex(int, int); // creates the outer-zone node
	void assignRank(int); // assigns rank to the from/to node
	void connectRank(int); // assigns rank to otherNodeGlobal

private:
	int fromToFlag; // Flag to indicate, whether the intra-zonal node end of the line is from (+1) or to (-1) end
	int fromToOuter; // Flag to indicate, whether the outer-zonal node end of the line is from (-1) or to (+1) end 
	int translID; // transmissionLine object id number
	int sharedIndex; // Serial number in the list of shared lines
	double ptMax; // Maximum MW transfer capability
	double reacT; // Resistance and Reactance of the transmission line
	Node *connNodetPtr; // pointers to the node objects at the two ends of the transmission line
	int otherNodeRank; // other-zone nodes' ranking in the list of other-zone nodes
	int otherNodeGlobal; // other-zone nodes' global ranking
	int fromZone; // From zone number
	int toZone; // To zone number
	int fromNode; // From node index
	int toNode; // To node index
	
}; //end class SELine

#endif // SELINE_H
	

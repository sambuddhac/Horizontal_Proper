// candLine class definition.
// Member functions defined in candidateLine.cpp
#ifndef CANDLINE_H
#define CANDLINE_H

class Node; // forward declaration

class candLine {

public:
	candLine( int, int, Node *, int, int, int, int, int, double, double, double, int, double, int, int ); // constructor
	~candLine(); // destructor
	int getTranslID(); // returns the ID of the transmission line
	int getIntlNodeID(); // returns ID number of intra-zonal node end to which the transmission line is connected
	int getIntlZoneID(); // returns ID number of intra-zonal node end zone to which the transmission line is connected
	int getExtNodeID(); // returns ID number of outer-zonal node end to which the transmission line is connected
	int getExtZoneID(); // returns ID number of outer-zonal node end zone to which the transmission line is connected
	int getExtNodeRank(); // returns serial rank in the list of nodes of the outer end node, to which the transmission line is connected
	int getExtNodeGlobalRank(); // returns global rank in the list of nodes of the outer end node, to which the transmission line is connected
	double getFlowLimit(); // Gets the value of power flow line limit
	int getFlowDir(); // returns the value of the direction flag indicating whether the intra-zonal node end of the line is from (+1) or to (-1) end
	double getReactance(); // Returns the values of the reactance
	void outerNodeIndex(int, int); // creates the outer-zone node
	void assignRank(int); // assigns rank to the from/to node
	void connectRank(int); // assigns rank to otherNodeGlobal
	void setTranData(double, int, double); 
	void assignLineRank(int); // Assigns global rank to the candidate line
	int getGlobalRank(); // Returns the global rank of the candidate line
	int getOtherZone(); // returns the ID number of the external zone 
	double getInvestCost(); // Returns the total investment cost
	int returnPresAbsStatus(); // Returns the construction status of the candidate line
	void setPresAbsStatus(); // Sets the construction status of the candidate line
	void modifyNodeReact(); // function to modify the nodal connected reactance, if the candidate line is actually built
	int returnOwnership(); // Returns the value of ownership

private:
	int fromToFlag; // Flag to indicate, whether the intra-zonal node end of the line is from (+1) or to (-1) end
	int fromToOuter; // Flag to indicate, whether the outer-zonal node end of the line is from (-1) or to (+1) end 
	int ownership; // Indicates the cost bearing for this line, of this zone. 2, 1, 0 for full, half, and no ownership, respectively;
	int translID; // transmissionLine object id number
	int sharedIndex; // Serial number in the list of shared candidate lines
	double ptMax; // Maximum MW transfer capability
	double reacT; // Resistance and Reactance of the transmission line
	Node *connNodetPtr; // pointers to the intra-zone node objects of the shared candidate transmission line
	int otherNodeRank; // other-zone nodes' ranking in the list of other-zone nodes
	int otherNodeGlobal; // other-zone nodes' global ranking
	int fromZone; // From zone number
	int toZone; // To zone number
	int fromNode; // From node index
	int toNode; // To node index
	double capitalCost; // Capital cost for the construction of the line
	int lifeYears; // Type of device; 1 for Generating device, 0 otherwise
	double rateInterest; // Rate of Interest
	int globalRank; // Global rank of the candidate line
	int resFromStageI; // Result from Stage-I regarding the presence or absence of the candidate transmission line

}; //end class candLine

#endif // CANDLINE_H
	

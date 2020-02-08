// Node class definition.
// Member functions defined in node.cpp
#ifndef NODE_H
#define NODE_H
#include <vector>

using namespace std;

class Node {

public:
	Node( int, int ); // constructor
	~Node(); // destructor
	int getNodeID(); // returns ID of the node to the caller
	int getGenSer(int); // returns the serial number of the generators connected to the node
	int getGenLength(); // Returns the number of generators connected to the node
	int getConNodeLength(); // Returns the number of intra zonal nodes connected to the node
	int getConnSer(int); // Returns the serial number of the intra zonal nodes connected to this node
	double getConnReact(int); // Returns the reactance of the connected internal node
	int getExtraNodeLength(); // Returns the number of outer zonal nodes connected to this node
	int getExtConnSer(int); // Returns the serial number of the variable listing of the outer zonal node
	double getExtConnReact(int); // Returns the reactance of the connected external node
	int getCandLineLengthF(); // Returns the number of candidate shared lines connected to this from node
	int getCandSerF(int); // Returns the serial number of the shared candidate line connected to this from node
	int getCandLineLengthT(); // Returns the number of candidate shared lines connected to this to node
	int getCandSerT(int); // Returns the serial number of the shared candidate line connected to this to node
	int getIntCandLineLengthF(); // Returns the number of candidate internal lines connected to this from node
	int getIntCandSerF(int); // Returns the serial number of the internal candidate line connected to this from node
	int getIntCandLineLengthT(); // Returns the number of candidate internal lines connected to this to node
	int getIntCandSerT(int); // Returns the serial number of the internal candidate line connected to this to node
	double devpinitMessage(int) const; // returns the total value of connected load for each scenario
	void setgConn( int ); // Function to set number of generators connected
	void settConn( int, int, double, int ); // Function to set number of transmission lines connected
	void setSEConn( int, int, double, int ); // Function to set the number of shared existing transmission lines connected 
	void setlConn( int, vector<double>* ); // Function to set number of loads connected
	double getToReact(); // Sum of reciprocal of reactances of all intra and shared exiisting lines for which this is to node
	double getFromReact(); // Sum of reciprocal of reactances of all intra and shared exiisting lines for which this is from node
	int getSharedFlag(); // returns the status if this node is connected to a shared existing transmission line
	int getCandFlag(); // returns the status if this node i connected to a shared cand transmission line
	void sendExtNodeInfo( int, int, double, int ); // Function to populate the connected outer-node list
	void setCandConn( int, int, double, int ); // Function to set the number of shared cand transmission lines connected 
	void setIntCandConn( int, int, double, int, int ); // Function to set the number of internal cand transmission lines connected
	void assignGlobalRank( int ); // Assigns the global rank to the nodes that are ends of shared lines
	void populateGlobalConn( int ); // Populates the extNodeGlobalRank vector with the global ranks
	int getGlobalRank(); // Returns the global rank of this node
	void modifyReactAPP( int, int, double, int, int ); // Modifies the to and from reactances of lines connected to this node, to account for the newly constructed lines
	void initLoad( int ); // Initialize the default loads on all nodes to zero
	int getNodeMultiplicity();// get the multiplicity of the node i.e: the number of different zones (other than the one where it belongs) to which it is connected
	int getBuiltCandFlag(); 

private:
	int nodeID; // Node object id number
	int zoneID; // Zone to which the node belongs
	int globalRank; // Global rank in the stitched list of shared line end nodes
	int gConnNumber, tConnNumber, lConnNumber, sharedExConnNumber, candConnNumber, builtCandConnNumber, builtIntCandConnNumber, intCandConnNumber; // number of generators, transmission lines, loads, shared existing, shared candidate lines, internal candidate lines, and built candidate lines connected to the node
	int PDevCount; // Number of devices connected to the particular node object
	int sharedFlag; // node flag to indicate whether a shared existing or candidate line has been connected to a node
	double fromReact; // Sum of reciprocals of reactances of lines for which this is the from node
	double toReact; // Sum of reciprocals of reactances of lines for which this is the to node
	int connSharedPoint=0; // flag to indicate if this node is either the from or to end of any SE line. Default value 0 indicates it isn't
	int connCandPoint=0; // flag to indicate if this node is either the from or to end of any Cand line. Default value 0 indicates it isn't
	int connBuiltCandPoint=0; // flag to indicate if this node is either the from or to end of any constructed Cand line. Default value 0 indicates it isn't
	int connIntCandPoint=0; // flag to indicate if this node is either the from or to end of any internal Cand line. Default value 0 indicates it isn't 
	vector<double> connLoadVal; // Scenario values of connected load to this node
	int multiplicity=0; // Count of number of different zones to which the node shares existing or candidate lines, default value is 0
	vector< int > genSerialNum; // vector consisting of the serial numbers of generators connected to a particular node
	vector< int > loadSerialNum; // vector consisting of the serial numbers of loads connected to a particular node
	vector< int > tranFromSerial; // vector consisting of the transmission lines for which the node is from node
	vector< int > tranToSerial; // vector consisting of the transmission lines for which the node is to node
	vector< int > connNodeList; // List of intra-zonal nodes that are directly connected to this node via transmission lines
	vector< double > connReactRec; // List of reciprocals of reactances of the intra zone lines connected to the node
	vector< int > SEFromSerial; // vector consisting of the SE lines for which the node is from node
	vector< int > SEToSerial; // vector consisting of the SE lines for which the node is to node
	vector< int > shareNodeList; // List of outer-zone nodes connected to the node via shared existing transmission lines
	vector< double > shareReactRec; // List of reciprocals of reactances of the shared existing lines connected to the node
	vector< int > CandFromSerial; // vector consisting of the cand lines for which the node is from node
	vector< int > CandToSerial; // vector consisting of the cand lines for which the node is to node
	vector< int > IntCandFromSerial; // vector consisting of the internal cand lines for which the node is from node
	vector< int > IntCandToSerial; // vector consisting of the internal cand lines for which the node is to node
	vector< int > builtCandFromSerial; // vector consisting of the built cand lines for which the node is from node
	vector< int > builtCandToSerial; // vector consisting of the built cand lines for which the node is to node
	vector< int > builtIntCandFromSerial; // vector consisting of the built internal cand lines for which the node is from node
	vector< int > builtIntCandToSerial; // vector consisting of the built internal cand lines for which the node is to node	
	vector< int > CandNodeList; // List of outer-zone nodes connected to the node via shared cand transmission lines
	vector< int > extNodeGlobalRank; // Global ranking of external zone nodes connected to this node via shared existing or candidate lines
	vector< double > CandReactRec; // List of reciprocals of reactances of the shared cand lines connected to the node
	vector< int > connectedZoneList; // List of connected zones for nodes connected to shared existing and/or shared candidate lines
}; //end class Node

#endif // NODE_H
	

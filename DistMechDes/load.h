// Load class definition.
// Member functions defined in load.cpp
#ifndef LOAD_H
#define LOAD_H
#include <vector>
using namespace std;
class Node; // forward declaration

class Load {

public:
	Load( int, Node *, int, double[] ); // constructor
	~Load(); // destructor
	void setLoadValue(double[]); // Sets the load for each scenario
	int getLoadID(); // returns the ID of the Load
	int getLoadNodeID(); // returns the ID number of the node to which the load is connected

private:
	int loadID; // Load object id number
	int numberOfScenarios; // Number of random scenarios
	int deviceNature; // 1 for Generating device, 0 otherwise
	vector<double> Pl; // MW consumption
	Node *connNodelPtr; // connection node object

}; //end class Load

#endif // LOAD_H
	

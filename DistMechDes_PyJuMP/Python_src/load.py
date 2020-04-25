// Member functions for class Load.
#include <iostream>
#include <vector>
// include Load class definition from load.h and node.h
#include "load.h"
#include "node.h"

using namespace std;

Load::Load( int idOfLoad, Node *nodeConnl, int scenarioCount, double Load_P[] ) // constructor begins
	: loadID( idOfLoad ),
	  connNodelPtr( nodeConnl ),
	  numberOfScenarios(scenarioCount),
	  deviceNature( 0 )
{
	setLoadValue(Load_P); // Sets the load for each scenario
	connNodelPtr->setlConn( idOfLoad, &Pl ); // increments the load connection variable to node
} // constructor ends

Load::~Load() // destructor
{
	//cout << "\nThe load object having ID " << loadID << " have been destroyed.\n";

} // end of destructor

void Load::setLoadValue(double Load_P[]) // Sets the load for each scenario
{
	for (int i=0; i<numberOfScenarios; ++i) {
		Pl.push_back(Load_P[i]);
	}
}

int Load::getLoadID() // function getLoadID begins
{
	return loadID; // returns the ID of the load object
} // end of getLoadID function

int Load::getLoadNodeID() // function getLoadNodeID begins
{
	return connNodelPtr->getNodeID(); // returns the ID number of the node to which the load object is connected
} // end of getLoadNodeID function
	

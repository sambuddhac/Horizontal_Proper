// Member functions for class transmissionLine.
#include <iostream>
#include <iomanip>
#include <cmath>
// include transmissionLine class definition from transl.h, Node class definition from node.h
#include "transl.h"
#include "node.h"

using namespace std;

transmissionLine::transmissionLine( int idOfTransl, Node *nodeConnt1, Node *nodeConnt2, double PowertMax, double Reactance) // constructor begins
	: translID( idOfTransl ),
	  connNodet1Ptr( nodeConnt1 ),
	  connNodet2Ptr( nodeConnt2 ),
	  ptMax( PowertMax ),
	  reacT( Reactance ),
	  deviceNature( 0 )
{
	//*cout << "\nInitializing the parameters of the transmission line with ID: " << translID << endl;
	int fromNode=connNodet1Ptr->getNodeID();
	int toNode=connNodet2Ptr->getNodeID();
	connNodet1Ptr->settConn( idOfTransl, 1, reacT, toNode ); // increments the txr line connection variable to node 1
	connNodet2Ptr->settConn( idOfTransl, -1, reacT, fromNode ); // increments the txr line connection variable to node 2
} // constructor ends

transmissionLine::~transmissionLine() // destructor
{
	//cout << "\nThe transmission line object having ID " << translID << " have been destroyed.\n";

} // end of destructor

int transmissionLine::getTranslID() // function gettranslID begins
{
	return translID; // returns the ID of the generator object
} // end of gettranslID function

double transmissionLine::getFlowLimit() // function getFlowLimit begins
{
	return ptMax; // returns the Maximum power flow limit
} // end of getFlowLimit function

int transmissionLine::getTranslNodeID1() // function getGenNodeID begins
{
	return connNodet1Ptr->getNodeID(); // returns the ID number of the node to which the generator object is connected
} // end of getGenNodeID function

int transmissionLine::getTranslNodeID2() // function getGenNodeID begins
{
	return connNodet2Ptr->getNodeID(); // returns the ID number of the node to which the generator object is connected
} // end of getGenNodeID function

double transmissionLine::getReactance()
{
	return reacT;
}
	


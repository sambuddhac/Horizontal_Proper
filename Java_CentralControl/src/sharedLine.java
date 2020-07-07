// Member functions for class SELine.
#include <iostream>
#include <iomanip>
#include <cmath>
// include SELine class definition from sharedLine.h, Node class definition from node.h
#include "sharedLine.h"
#include "node.h"

using namespace std;

SELine::SELine( int localRank, int idOfTransl, Node *nodeConnt1, Node *nodeConnt2, double PowertMax, double Reactance ) // constructor begins
	: translID( idOfTransl ),
	  localIndex( localRank ),
	  connNodetPtr1( nodeConnt1 ),
	  connNodetPtr2( nodeConnt2 ),
	  ptMax( PowertMax ),
	  reacT( Reactance )	  
{
	int fromNode=connNodetPtr1->getNodeID();
	int toNode=connNodetPtr2->getNodeID();
	cout << "\nInitializing the parameters of the shared transmission line with ID: " << translID << endl;
	cout << "from node: " << fromNode << " To node: " << toNode << endl;
	connNodetPtr1->setSEConn( idOfTransl, 1, reacT, toNode ); // increments the txr line connection variable to node 1
	connNodetPtr2->setSEConn( idOfTransl, -1, reacT, fromNode ); // increments the txr line connection variable to node 1
} // constructor ends

SELine::~SELine() // destructor
{
	//cout << "\nThe transmission line object having ID " << translID << " have been destroyed.\n";

} // end of destructor

int SELine::getTranslID() // function gettranslID begins
{
	return translID; // returns the ID of the generator object
} // end of gettranslID function

int SELine::getFromNodeID() // function getFromNodeID begins
{
	return connNodetPtr1->getNodeID(); // returns the ID number of the from node
} // end of getFromNodeID function

int SELine::getToNodeID() // function getToNodeID begins
{
	return connNodetPtr2->getNodeID(); // returns the ID number of the to node
} // end of getToNodeID function


double SELine::getFlowLimit() // Function getFlowLimit gets the value of power flow line limit	
{
	return ptMax;
} // Function getFlowLimit ends

double SELine::getReactance()
{
	return reacT;
}

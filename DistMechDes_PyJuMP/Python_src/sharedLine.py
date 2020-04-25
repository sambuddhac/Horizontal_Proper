// Member functions for class SELine.
#include <iostream>
#include <iomanip>
#include <cmath>
// include SELine class definition from sharedLine.h, Node class definition from node.h
#include "sharedLine.h"
#include "node.h"

using namespace std;

SELine::SELine( int sharedRank, int idOfTransl, Node *nodeConnt, int fromN, int fromZ, int toN, int toZ, int zonalIDNum, double PowertMax, double Reactance ) // constructor begins
	: translID( idOfTransl ),
	  sharedIndex( sharedRank ),
	  connNodetPtr( nodeConnt ),
	  ptMax( PowertMax ),
	  reacT( Reactance ),
	  fromZone(fromZ),
	  toZone (toZ),
	  fromNode (fromN),
	  toNode (toN),
	  otherNodeGlobal(0)
{
	//*cout << "\nInitializing the parameters of the transmission line with ID: " << translID << endl;
	if (fromZone==zonalIDNum) {
		connNodetPtr->setSEConn( idOfTransl, 1, reacT, toZone ); // increments the txr line connection variable to node 1
		fromToFlag=1;
	}
	else {
		connNodetPtr->setSEConn( idOfTransl, -1, reacT, fromZone ); // increments the txr line connection variable to node 1
		fromToFlag=-1;
	}

} // constructor ends

SELine::~SELine() // destructor
{
	//cout << "\nThe transmission line object having ID " << translID << " have been destroyed.\n";

} // end of destructor

int SELine::getTranslID() // function gettranslID begins
{
	return translID; // returns the ID of the generator object
} // end of gettranslID function

int SELine::getIntlNodeID() // function getGenNodeID begins
{
	return connNodetPtr->getNodeID(); // returns the ID number of the node to which the generator object is connected
} // end of getGenNodeID function

int SELine::getIntlZoneID() // returns ID number of intra-zonal node end zone to which the transmission line is connected
{
	if (fromToFlag==1)
		return fromZone; // returns the ID number of the from zone if the intra-zonal node is the from node
	else 
		return toZone; // returns the ID number of the to zone if the intra-zonal node is the to node 
} // end of getIntlZoneID() function

int SELine::getExtNodeID() // returns ID number of outer-zonal node end to which the transmission line is connected
{
	if (fromToFlag==1)
		return toNode; // returns the ID number of the from zone if the intra-zonal node is the from node
	else 
		return fromNode; // returns the ID number of the to zone if the intra-zonal node is the to node
} // end of getGenNodeID function

int SELine::getExtZoneID() // returns ID number of outer-zonal node end zone to which the transmission line is connected
{
	if (fromToFlag==1)
		return toZone; // returns the ID number of the from zone if the intra-zonal node is the from node
	else 
		return fromZone; // returns the ID number of the to zone if the intra-zonal node is the to node
} // end of getGenNodeID function

int SELine::getExtNodeRank() // function getGenNodeID begins
{
	return otherNodeRank; // returns the ID number of the node to which the generator object is connected
} // end of getGenNodeID function

int SELine::getExtNodeGlobalRank() // function getExtNodeGlobalRank for the outside-zone node ID begins
{
	return otherNodeGlobal; // returns the global rank of the outside-zone node to which the SE line object is connected
} // end of getExtNodeGlobalRank function

double SELine::getFlowLimit() // Function getFlowLimit gets the value of power flow line limit	
{
	return ptMax;
} // Function getFlowLimit ends

int SELine::getFlowDir() // returns the value of the direction flag indicating whether the intra-zonal node end of the line is from (+1) or to (-1) end
{
	return fromToFlag;
} // Function getFlowDir ends

void SELine::outerNodeIndex(int rankOfOuterNode, int dirFlag)
{
	otherNodeRank=rankOfOuterNode;
	fromToOuter=dirFlag;
	connNodetPtr->sendExtNodeInfo( otherNodeRank, fromToOuter, reacT, 0 );
}

double SELine::getReactance()
{
	return reacT;
}

void SELine::assignRank(int ranking) // assigns rank to the from/to node
{
	connNodetPtr->assignGlobalRank(ranking);		
}

void SELine::connectRank(int ranking) // assigns rank to otherNodeGlobal
{
	otherNodeGlobal = ranking;
	connNodetPtr->populateGlobalConn(ranking);		
}

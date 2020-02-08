// Member functions for class Node.
#include <iostream>
// include Node class definition from node.h
#include <vector>
#include <algorithm>
#include "node.h"

using namespace std;

Node::Node( int universalID, int idOfNode, int zoneIndex ) // constructor begins
	: nodeID( idOfNode ),
	  zoneID( zoneIndex )
{
	//cout << "\nInitializing the parameters of the node with ID: " << nodeID << endl;

	// initialize the connected devices to zero for node
	gConnNumber = 0; // number of generators connected to a particular node
	tConnNumber = 0; // number of transmission lines connected to a particular node
	lConnNumber = 0; // number of loads connected to a particular node
	sharedExConnNumber = 0; // number of shared existing transmission lines connected to this node
	builtCandConnNumber = 0; // number of constructed candidate line connected to this node
	candConnNumber = 0; // number of shared candidate transmission lines connected to this node
	intCandConnNumber = 0; // number of internal candidate transmission lines connected to this node 
	sharedFlag = 0; // node flag to indicate whether a shared existing or candidate line has been connected to a node
	PDevCount = 0; // initialize number of devices connectedto a node to zero
	fromReact = 0.0; // Initialize the from reactance
	toReact = 0.0; // Initialize the to reactance
	globalRank = universalID; // sets the globalRank to universalID 
	cout << "Node number: " << globalRank << endl;
	cout << "Zone number: " << zoneID << endl;

} // constructor ends

Node::~Node() // destructor
{
	//cout << "\nThe node object having ID " << nodeID << " have been destroyed.\n";

} // end of destructor

int Node::getNodeID() // function getNodeID begins
{
	return globalRank; // returns node ID to the caller
} // end of function getNodeID

void Node::setgConn( int serialOfGen )
{
	++gConnNumber; // increment the number of generators connected by one whenever a generator is connected to the node
	genSerialNum.push_back( serialOfGen ); // records the serial number of the generator connected to the node 
}

void Node::settConn( int tranID, int dir, double react, int rankOfOther )
{
	++tConnNumber; // increment the number of txr lines connected by one whenever a txr line is connected to the node
	if ( dir == 1 ) {
		tranFromSerial.push_back(tranID);
		fromReact += (1/react);	
		if (std::find(connNodeList.begin(), connNodeList.end(), rankOfOther) != connNodeList.end()) { // If predecided Gen value is given for this particular Powergenerator
			auto pos = std::find(connNodeList.begin(), connNodeList.end(), rankOfOther) - connNodeList.begin(); // find the position of the Powergenerator in the chart of predecided values
			connReactRec[pos] -= 1/react;
		}
		else {
			connNodeList.push_back(rankOfOther);
			connReactRec.push_back(-1/react);
		}
	
	}
	else {
		tranToSerial.push_back(tranID);
		toReact -= (1/react);
		if (std::find(connNodeList.begin(), connNodeList.end(), rankOfOther) != connNodeList.end()) { // If predecided Gen value is given for this particular Powergenerator
			auto pos = std::find(connNodeList.begin(), connNodeList.end(), rankOfOther) - connNodeList.begin(); // find the position of the Powergenerator in the chart of predecided values
			connReactRec[pos] += 1/react;
		}
		else {
			connNodeList.push_back(rankOfOther);
			connReactRec.push_back(1/react);
		}
	}
}

void Node::setSEConn( int tranID, int dir, double react, int rankOfOther )
{
	++sharedExConnNumber; // increment the number of shared existing txr lines connected by one whenever a txr line is connected to the node
	connSharedPoint = 1; // set the flag to indicate this is the from or to end of a shared existing transmission line
	if ( dir == 1 ) {
		SEFromSerial.push_back(tranID);
		fromReact += (1/react);	
		if (std::find(connNodeList.begin(), connNodeList.end(), rankOfOther) != connNodeList.end()) { // If predecided Gen value is given for this particular Powergenerator
			auto pos = std::find(connNodeList.begin(), connNodeList.end(), rankOfOther) - connNodeList.begin(); // find the position of the Powergenerator in the chart of predecided values
			connReactRec[pos] -= 1/react;
		}
		else {
			connNodeList.push_back(rankOfOther);
			connReactRec.push_back(-1/react);
		}	
	}
	else {
		SEToSerial.push_back(tranID);
		toReact -= (1/react);
		if (std::find(connNodeList.begin(), connNodeList.end(), rankOfOther) != connNodeList.end()) { // If predecided Gen value is given for this particular Powergenerator
			auto pos = std::find(connNodeList.begin(), connNodeList.end(), rankOfOther) - connNodeList.begin(); // find the position of the Powergenerator in the chart of predecided values
			connReactRec[pos] += 1/react;
		}
		else {
			connNodeList.push_back(rankOfOther);
			connReactRec.push_back(1/react);
		}
	}
	cout << "The node " << globalRank << " in zone " << zoneID << " accessed by the SE line " << tranID << " with reactance " << react << " and dir " << dir << endl;
}

void Node::setCandConn( int tranID, int dir, double react, int rankOfOther )
{
	++candConnNumber; // increment the number of shared cand txr lines connected by one whenever a txr line is connected to the node
	if ( dir == 1 ) {
		CandFromSerial.push_back(tranID);	
	}
	else {
		CandToSerial.push_back(tranID);
	}
	connCandPoint=1; // Flag set to indicate that this node is connected to a cand line
	cout << "The node " << nodeID << " in zone " << zoneID << " accessed by the shared candidate line " << tranID << " with reactance " << react << endl;
}

void Node::setIntCandConn( int tranID, int dir, double react, int rankOfOther, int constStat )
{
	++intCandConnNumber; // increment the number of shared cand txr lines connected by one whenever a txr line is connected to the node
	if ( dir == 1 ) {
		IntCandFromSerial.push_back(tranID);
		fromReact += constStat*(1/react);	
		if (std::find(connNodeList.begin(), connNodeList.end(), rankOfOther) != connNodeList.end()) { // If predecided Gen value is given for this particular Powergenerator
			auto pos = std::find(connNodeList.begin(), connNodeList.end(), rankOfOther) - connNodeList.begin(); // find the position of the Powergenerator in the chart of predecided values
			connReactRec[pos] -= constStat*(1/react);
		}
		else {
			connNodeList.push_back(rankOfOther);
			connReactRec.push_back(constStat*(-1/react));
		}	
	}
	else {
		IntCandToSerial.push_back(tranID);
		toReact -= constStat*(1/react);
		if (std::find(connNodeList.begin(), connNodeList.end(), rankOfOther) != connNodeList.end()) { // If predecided Gen value is given for this particular Powergenerator
			auto pos = std::find(connNodeList.begin(), connNodeList.end(), rankOfOther) - connNodeList.begin(); // find the position of the Powergenerator in the chart of predecided values
			connReactRec[pos] += constStat*(1/react);
		}
		else {
			connNodeList.push_back(rankOfOther);
			connReactRec.push_back(constStat*(1/react));
		}
	}
	connIntCandPoint=1; // Flag set to indicate that this node is connected to an internal cand line
}

void Node::setlConn( int lID, vector<double>* loadVal )
{
	++lConnNumber; // increment the number of loads connected by one whenever a load is connected to the node
	loadSerialNum.push_back(lID);
	connLoadVal.clear();
	connLoadVal = *loadVal; // total connected load

}

int Node::getGenLength() // function getNodeID begins
{
	return genSerialNum.size(); // returns node ID to the caller
} // end of function getNodeID

int Node::getGenSer(int colCount)
{
	return genSerialNum.at(colCount-1);
}

// function redContNodeCount begins

void Node::initLoad( int scenNum ) // Initialize the default loads on all nodes to zero
{
	for (int i = 0; i < scenNum; ++i)
		connLoadVal.push_back(0);
}

double Node::devpinitMessage(int scenC) const// function devpinitMessage begins
{
	return connLoadVal.at(scenC); // return the total connected load
} // function devpinitMessage ends

int Node::getSharedFlag()
{
	return connSharedPoint; // return the status if this node is connected to a shared existing line
}		

int Node::getCandFlag()
{
	return connCandPoint; // return the status if this node is connected to a shared cand line
}

double Node::getToReact()
{
	return toReact; // return the total reciprocal of reactances for which this is the to node
}

double Node::getFromReact()
{
	return fromReact; // return the total reciprocal of reactances for which this is the from node
}

int Node::getConNodeLength()
{
	return connNodeList.size(); // returns the length of the vector containing the connected intra-zonal nodes
}

int Node::getConnSer(int colCount)
{
	return connNodeList.at(colCount-1); // returns the serial number of the connected internal node at this position
}

double Node::getConnReact(int colCount)
{
	return connReactRec.at(colCount-1); // returns the serial number of the connected internal node at this position
}

int Node::getCandLineLengthF()
{
	return CandFromSerial.size(); // returns the number of cand lines connected to this from node
}

int Node::getCandLineLengthT()
{
	return CandToSerial.size(); // returns the number of cand lines connected to this to node
}

int Node::getCandSerF(int colCount)
{
	return CandFromSerial.at(colCount-1); // returns the serial number of the cand line at this position
}

int Node::getCandSerT(int colCount)
{
	return CandToSerial.at(colCount-1); // returns the serial number of the cand line at this position
}

int Node::getIntCandLineLengthF()
{
	return IntCandFromSerial.size(); // returns the number of cand lines connected to this from node
}

int Node::getIntCandLineLengthT()
{
	return IntCandToSerial.size(); // returns the number of cand lines connected to this to node
}

int Node::getIntCandSerF(int colCount)
{
	return IntCandFromSerial.at(colCount-1); // returns the serial number of the cand line at this position
}

int Node::getIntCandSerT(int colCount)
{
	return IntCandToSerial.at(colCount-1); // returns the serial number of the cand line at this position
}


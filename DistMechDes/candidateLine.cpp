// Member functions for class candLine.
#include <iostream>
#include <iomanip>
#include <cmath>
// include candLine class definition from candidateLine.h, Node class definition from node.h
#include "candidateLine.h"
#include "node.h"

using namespace std;

candLine::candLine( int sharedRank, int idOfTransl, Node *nodeConnt, int fromN, int fromZ, int toN, int toZ, int zonalIDNum, double PowertMax, double Reactance, double ROI, int life, double cap, int absPres, int owner ) // constructor begins
	: translID( idOfTransl ),
	  sharedIndex( sharedRank ),
	  connNodetPtr( nodeConnt ),
	  ptMax( PowertMax ),
	  reacT( Reactance ),
	  fromZone(fromZ),
	  toZone (toZ),
	  fromNode (fromN),
	  toNode (toN),
	  otherNodeGlobal(0),
	  resFromStageI(absPres),
	  globalRank(0),
	  ownership(owner)
{
	//*cout << "\nInitializing the parameters of the transmission line with ID: " << translID << endl;
	if (fromZone==zonalIDNum) {
		connNodetPtr->setCandConn( idOfTransl, 1, reacT, toZone ); // increments the txr line connection variable to node 1
		fromToFlag=1;
	}
	else {
		connNodetPtr->setCandConn( idOfTransl, -1, reacT, fromZone ); // increments the txr line connection variable to node 1
		fromToFlag=-1;
	}
	setTranData(cap, life, ROI); // calls setTranData member function to set the parameter values

} // constructor ends

candLine::~candLine() // destructor
{
	//cout << "\nThe transmission line object having ID " << translID << " have been destroyed.\n";

} // end of destructor

int candLine::returnOwnership() // Returns the value of ownership
{
	return ownership;
}

void candLine::modifyNodeReact() // function to modify the nodal connected reactance, if the candidate line is actually built
{
	// If the connected node is the from node
	if (fromToFlag==1) {
		connNodetPtr->modifyReactAPP( translID, 1, reacT, otherNodeRank, 1 ); // increments the txr line connection variable to node 1
		connNodetPtr->sendExtNodeInfo( otherNodeRank, fromToOuter, reacT, 0 );
	}
	// If the connected node is the to node
	else {
		connNodetPtr->modifyReactAPP( translID, -1, reacT, otherNodeRank, 1 ); // increments the txr line connection variable to node 1
		connNodetPtr->sendExtNodeInfo( otherNodeRank, fromToOuter, reacT, 0 );
	}
} // function ends

int candLine::getTranslID() // function gettranslID begins
{
	return translID; // returns the ID of the generator object
} // end of gettranslID function

int candLine::getIntlNodeID() // function getGenNodeID for the intra-zone node ID begins
{
	return connNodetPtr->getNodeID(); // returns the ID number of the intra-zone node to which the candidate line object is connected
} // end of getGenNodeID function

int candLine::getIntlZoneID() // returns ID number of intra-zonal node end zone to which the transmission line is connected
{
	if (fromToFlag==1)
		return fromZone; // returns the ID number of the from zone if the intra-zonal node is the from node
	else 
		return toZone; // returns the ID number of the to zone if the intra-zonal node is the to node 
} // end of getIntlZoneID() function

int candLine::getExtNodeID() // returns ID number of outer-zonal node end to which the transmission line is connected
{
	if (fromToFlag==1)
		return toNode; // returns the ID number of the from zone if the intra-zonal node is the from node
	else 
		return fromNode; // returns the ID number of the to zone if the intra-zonal node is the to node
} // end of getGenNodeID function

int candLine::getExtZoneID() // returns ID number of outer-zonal node end zone to which the transmission line is connected
{
	if (fromToFlag==1)
		return toZone; // returns the ID number of the from zone if the intra-zonal node is the from node
	else 
		return fromZone; // returns the ID number of the to zone if the intra-zonal node is the to node
} // end of getGenNodeID function

int candLine::getExtNodeRank() // function getGenNodeID for the outside-zone node ID begins
{
	return otherNodeRank; // returns the ID number of the outside-zone node to which the candidate line object is connected
} // end of getGenNodeID function

int candLine::getExtNodeGlobalRank() // function getExtNodeGlobalRank for the outside-zone node ID begins
{
	return otherNodeGlobal; // returns the global rank of the outside-zone node to which the candidate line object is connected
} // end of getExtNodeGlobalRank function

double candLine::getFlowLimit() // Function getFlowLimit gets the value of power flow line limit	
{
	return ptMax;
} // Function getFlowLimit ends

int candLine::getFlowDir() // returns the value of the direction flag indicating whether the intra-zonal node end of the line is from (+1) or to (-1) end
{
	return fromToFlag;
} // Function getFlowDir ends

void candLine::outerNodeIndex(int rankOfOuterNode, int dirFlag)
{
	otherNodeRank=rankOfOuterNode;
	fromToOuter=dirFlag;
	connNodetPtr->sendExtNodeInfo( otherNodeRank, fromToOuter, reacT, 1 );
}

double candLine::getReactance()
{
	return reacT;
}

void candLine::assignRank(int ranking) // assigns rank to the from/to node
{
	connNodetPtr->assignGlobalRank(ranking);		
}

void candLine::connectRank(int ranking) // assigns rank to otherNodeGlobal
{
	otherNodeGlobal = ranking;
	connNodetPtr->populateGlobalConn(ranking);		
}

int candLine::getOtherZone() // function getOtherZone returns the ID number of the outside zone to which the other node is connected
{
	if (fromToFlag==1)
		return toZone; // returns the ID number of the outside-zone to which the candidate line object is connected
	else
		return fromZone;
} // end of getOtherZone function

void candLine::setTranData(double capC, int lifeTime, double interRate) // member function to set parameter values of transmission lines
{
	capitalCost=capC;
	lifeYears=lifeTime;
	rateInterest=interRate;
} // end function for setting parameter values

double candLine::getInvestCost() //member function getInvestCost begins
{
	return (capitalCost*rateInterest*(pow((1+rateInterest), lifeYears)))/(pow((1+rateInterest), lifeYears)-1);//(1+rateInterest);capitalCost/100;//
	//return capitalCost;
}

void candLine::assignLineRank(int globRank) // Assigns global rank to the candidate line
{
	globalRank = globRank; // Global rank of the candidate line
}

int candLine::getGlobalRank() // Returns the global rank of the candidate line
{
	return globalRank; // Global rank of the candidate line
}

int candLine::returnPresAbsStatus() // Returns the construction status of the candidate line
{
	return resFromStageI; 
}

void candLine::setPresAbsStatus() // Sets the construction status of the candidate line
{
	resFromStageI=1; 
}


// Member functions for class candLine.
#include <iostream>
#include <iomanip>
#include <cmath>
// include candLine class definition from candidateLine.h, Node class definition from node.h
#include "candidateLine.h"
#include "node.h"

using namespace std;

candLine::candLine( int localRank, int idOfTransl, Node *nodeConnt1, Node *nodeConnt2, double PowertMax, double Reactance, double ROI, int life, double cap, int absPres ) // constructor begins
	: translID( idOfTransl ),
	  localIndex( localRank ),
	  connNodetPtr1( nodeConnt1 ),
	  connNodetPtr2( nodeConnt2 ),
	  ptMax( PowertMax ),
	  reacT( Reactance ),
	  status(absPres)
{
	int fromNode=connNodetPtr1->getNodeID();
	int toNode=connNodetPtr2->getNodeID();
	//*cout << "\nInitializing the parameters of the transmission line with ID: " << translID << endl;
	connNodetPtr1->setCandConn( idOfTransl, 1, reacT, toNode ); // increments the txr line connection variable to node 1
	connNodetPtr2->setCandConn( idOfTransl, -1, reacT, fromNode ); // increments the txr line connection variable to node 1
	setTranData(cap, life, ROI); // calls setTranData member function to set the parameter values

} // constructor ends

candLine::~candLine() // destructor
{
	//cout << "\nThe transmission line object having ID " << translID << " have been destroyed.\n";

} // end of destructor

int candLine::getTranslID() // function gettranslID begins
{
	return translID; // returns the ID of the generator object
} // end of gettranslID function

int candLine::getFromNodeID() // function getFromNodeID for the from node ID begins
{
	return connNodetPtr1->getNodeID(); // returns the ID number of the from node
} // end of getFromNodeID function

int candLine::getToNodeID() // function getToNodeID for the to node ID begins
{
	return connNodetPtr2->getNodeID(); // returns the ID number of the to node
} // end of getToNodeID function

double candLine::getFlowLimit() // Function getFlowLimit gets the value of power flow line limit	
{
	return ptMax;
} // Function getFlowLimit ends

double candLine::getReactance()
{
	return reacT;
}

void candLine::setTranData(double capC, int lifeTime, double interRate) // member function to set parameter values of transmission lines
{
	capitalCost=capC;
	lifeYears=lifeTime;
	rateInterest=interRate;
} // end function for setting parameter values

double candLine::getInvestCost() //member function getInvestCost begins
{
	//return (capitalCost*rateInterest*(pow((1+rateInterest), lifeYears)))/(pow((1+rateInterest), lifeYears)-1);//(1+rateInterest);capitalCost/100;
	return capitalCost;
}

int candLine::returnPresAbsStatus() // Returns the construction status of the candidate line
{
	return status;
}


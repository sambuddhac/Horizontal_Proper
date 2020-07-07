// Member functions for class intCandLine.
#include <iostream>
#include <iomanip>
#include <cmath>
// include candLine class definition from intcandidateLine.h, Node class definition from node.h
#include "intcandidateLine.h"
#include "node.h"

using namespace std;

intCandLine::intCandLine( int idOfTransl, Node *nodeConnt1, Node *nodeConnt2, double PowertMax, double Reactance, double ROI, int life, double cap, int absPres ) // constructor begins
	: translID( idOfTransl ),
	  connNodetPtr1( nodeConnt1 ),
	  connNodetPtr2( nodeConnt2 ),
	  ptMax( PowertMax ),
	  reacT( Reactance ),
	  statusOfConstruction(absPres)
{
	//*cout << "\nInitializing the parameters of the transmission line with ID: " << translID << endl;
	int fromNode=connNodetPtr1->getNodeID();
	int toNode=connNodetPtr2->getNodeID();
	connNodetPtr1->setIntCandConn( idOfTransl, 1, reacT, toNode, statusOfConstruction ); // increments the txr line connection variable to node 1
	connNodetPtr2->setIntCandConn( idOfTransl, -1, reacT, fromNode, statusOfConstruction ); // increments the txr line connection variable to node 2
	setTranData(cap, life, ROI); // calls setTranData member function to set the parameter values

} // constructor ends

intCandLine::~intCandLine() // destructor
{
	//cout << "\nThe transmission line object having ID " << translID << " have been destroyed.\n";

} // end of destructor

int intCandLine::getTranslID() // function gettranslID begins
{
	return translID; // returns the ID of the generator object
} // end of gettranslID function


int intCandLine::getTranslNodeID1() // function getGenNodeID begins
{
	return connNodetPtr1->getNodeID(); // returns the ID number of the node to which the generator object is connected
} // end of getGenNodeID function

int intCandLine::getTranslNodeID2() // function getGenNodeID begins
{
	return connNodetPtr2->getNodeID(); // returns the ID number of the node to which the generator object is connected
} // end of getGenNodeID function

double intCandLine::getFlowLimit() // Function getFlowLimit gets the value of power flow line limit	
{
	return ptMax;
} // Function getFlowLimit ends

double intCandLine::getReactance()
{
	return reacT;
}

void intCandLine::setTranData(double capC, int lifeTime, double interRate) // member function to set parameter values of transmission lines
{
	capitalCost=capC;
	lifeYears=lifeTime;
	rateInterest=interRate;
} // end function for setting parameter values

double intCandLine::getInvestCost() //member function getInvestCost begins
{
	//return (capitalCost*rateInterest*(pow((1+rateInterest), lifeYears)))/(pow((1+rateInterest), lifeYears)-1);//(1+rateInterest);capitalCost/100;//
	return capitalCost;
}

int intCandLine::returnPresAbsStatus() // Returns the construction status of the candidate line
{
	return statusOfConstruction; 
}


#include <vector>
#include "powergenerator.h" // Include the definition of powergenerator class
// Include definition of Node class 
#include "node.h"
#include <iostream>
using namespace std;

// constructor definition for the piecewise linear objective function
Powergenerator::Powergenerator(int ID, Node *nodeConng, double incCost, double noLoad, double max, double Min)
: genID(ID), // Initializer list to initialize data members that don't need validity check
connNodegPtr( nodeConng ),
noLoadCost(noLoad)
{
	incrementalCost = incCost; // Initialize this data member to zero (unused for piecewise linear objective)
	//cout << "Incremental cost " << incrementalCost << endl;
	setGenParamsSimple(max, Min); // call the set function to perform validity check on parameter value ranges and assign the values 
	//cout << "Limits defined" << endl;
	connNodegPtr->setgConn( genID ); // increments the generation connection variable to node
	//cout << "Node connected" << endl;
} // end of constructor


void Powergenerator::setGenParamsSimple(double max, double Min) // set function to set the Powergenerator class data members min max limits
{
		pMax = max;
		pMin = Min;
} // end of setGenParams function

Powergenerator::~Powergenerator() // destructor definition
{
	//cout << "\nGenerator object " << genNum << "  destroyed" << endl;
} // end of destructor

int Powergenerator::getGenID() // returns the Powergenerator ID number
{
	return genID;
} // end of getGenID function

double Powergenerator::getPMax() // function getPMax begins
{
	return pMax; 
} // getPMax ends

double Powergenerator::getPMin() // function getPMin begins
{
	return pMin;
} // getPMax ends

double Powergenerator::getLinCoeff() // Gets the linear coefficient (Incremental production cost
{
	return incrementalCost;
} // function getLinCoeff ends

double Powergenerator::getNLCost() // Gets the no load cost 
{
	return noLoadCost;
} // function getNLCost ends 


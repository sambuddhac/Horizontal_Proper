#include <vector>
from numpy import array
#include "powergenerator.h" // Include the definition of powergenerator class
import powergenerator #I am not sure abot this line
# Include definition of Node class 
#include "node.h"
from Python_src.node import node
#include <iostream>
from math import *
#using namespace std;

# constructor definition for the piecewise linear objective function
class Powergenerator(object):
	def __init__(self,int(ID), Node.nodeConng, float(incCost), float(noLoad), float(Max), float(Min)):
 		self.genID=ID
 # Initializer list to initialize data members that don't need validity check
 		self.connNodegPtr= nodeConng
 		self.noLoadCost=noLoad
 		self.incrementalCost = incCost # Initialize this data member to zero (unused for piecewise linear objective)
	#cout << "Incremental cost " << incrementalCost << endl;
		self.setGenParamsSimple(max, Min) # call the set function to perform validity check on parameter value ranges and assign the values 
	#cout << "Limits defined" << endl;
		self.connNodegPtr.setgConn(self.genID) # increments the generation connection variable to node
	#cout << "Node connected" << endl;
 # end of constructor


 def setGenParamsSimple(self, float(Max), float(Min)): # set function to set the Powergenerator class data members min max limits
	self.pMax = Max
	self.pMin = Min
 # end of setGenParams function
"""
Powergenerator::~Powergenerator() // destructor definition
{
	//cout << "\nGenerator object " << genNum << "  destroyed" << endl;
} // end of destructor
"""
def getGenID(self): # returns the Powergenerator ID number
	return self.genID
 # end of getGenID function

def PowerggetPMax(self): # function getPMax begins
	return self.pMax 
# getPMax ends

def getPMin(self): # function getPMin begins
	return self.pMin
# getPMax ends

def getLinCoeff(self): # Gets the linear coefficient (Incremental production cost
	return self.incrementalCost
# function getLinCoeff ends

def getNLCost(self): # Gets the no load cost 
	return self.noLoadCost
# function getNLCost ends 


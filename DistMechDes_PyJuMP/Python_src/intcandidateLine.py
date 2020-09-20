#Member functions for class intCandLine.
from math import *
import numpy as np
#include <iostream>
#include <iomanip>
#include <cmath>
// include candLine class definition from intcandidateLine.h, Node class definition from node.h
#include "intcandidateLine.h"
from Python_src.node import node
#include "node.h"
#using namespace std;

class intCandLine(object):
	def _init_(self, idOfTransl, nodeConnt1, nodeConnt2, PowertMax, Reactance, ROI, life, cap, absPres): # constructor begins
		self.translID=idOfTransl
	  	self.connNodetPtr1=nodeConnt1
	  	self.connNodetPtr2=nodeConnt2
	  	self.ptMax=PowertMax
	  	self.reacT=Reactance
	  	self.resFromStageI=absPres
	  	self.statusOfConstruction=absPres
		self.fromNode= self.connNodetPtr1.getNodeID(self) ### I am not sure about this definition
		self.toNode= self.connNodetPtr2.getNodeID(self)    ### same
		self.connNodetPtr1.setIntCandConn(self, idOfTransl, 1, reacT, toNode, statusOfConstruction) # increments the txr line connection variable to node 1
		self.connNodetPtr2.setIntCandConn(self,idOfTransl, -1, reacT, fromNode, statusOfConstruction) # increments the txr line connection variable to node 2
		self.setTranData(self, cap, life, ROI) # calls setTranData member function to set the parameter values

# constructor ends
"""
intCandLine::~intCandLine() // destructor
{
	//cout << "\nThe transmission line object having ID " << translID << " have been destroyed.\n";

} // end of destructor
"""
	def modifyNodeReact(self): # function to modify the nodal connected reactance, if the candidate line is actually built
		self.fromNode=self.connNodetPtr1.getNodeID(self)
		self.toNode=self.connNodetPtr2.getNodeID(self)
		self.connNodetPtr1.modifyReactAPP(self, translID, 1, reacT, toNode, 0) #increments the txr line connection variable to node 1
		self.connNodetPtr2.modifyReactAPP(self, translID, -1, reacT, fromNode, 0)  #increments the txr line connection variable to node 1
#	function ends

	def getTranslID(self): # function gettranslID begins
		return self.translID # returns the ID of the generator object
 # end of gettranslID function


	def getTranslNodeID1(self): # function getGenNodeID begins
		return self.connNodetPtr1.getNodeID(self) # returns the ID number of the node to which the generator object is connected
# end of getGenNodeID function

	def getTranslNodeID2(self): # function getGenNodeID begins
		return self.connNodetPtr2.getNodeID(self) # returns the ID number of the node to which the generator object is connected
# end of getGenNodeID function

	def getFlowLimit(self): # Function getFlowLimit gets the value of power flow line limit	
		return self.ptMax
# Function getFlowLimit ends

	def getReactance(self):
		return self.reacT	### I am not sure about self, actually I am a bit confused with "self"

	def setTranData(self, capC, lifeTime, interRate): # member function to set parameter values of transmission lines
		self.capitalCost=capC
		self.lifeYears=lifeTime
		self.rateInterest=interRate
# end function for setting parameter values

	def getInvestCost(self): #member function getInvestCost begins
		return (self.capitalCost*self.rateInterest*(pow((1+self.rateInterest), self.lifeYears)))/(pow((1+self.rateInterest), self.lifeYears)-1) #(1+rateInterest);capitalCost/100;//
	#return capitalCost;

	def returnPresAbsStatus(self): # Returns the construction status of the candidate line
		return self.resFromStageI

	def setPresAbsStatus(self): # Sets the construction status of the candidate line
		self.resFromStageI=1 


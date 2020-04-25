#Member functions for class intCandLine
#include class definition from intcandidateLine, Node class definition from node
from Python_src.node import Node
from Python_src.log import log

class intCandLine(object):
	def __init__(self, idOfTransl, nodeConnt1, nodeConnt2, PowertMax, Reactance, ROI, life, cap, absPres ): #constructor begins
		self.translID = idOfTransl
		self.connNodetPtr1 = nodeConnt1
		self.connNodetPtr2 = nodeConnt2
		self.ptMax = PowertMax
		self.reacT = Reactance
		self.statusOfConstruction = absPres
		self.fromNode = self.connNodetPtr1.getNodeID()
		self.toNode = self.connNodetPtr2.getNodeID()
		self.connNodetPtr1.setIntCandConn(idOfTransl, 1, self.reacT, self.toNode, self.statusOfConstruction) #increments the txr line connection variable to node 1
		self.connNodetPtr2.setIntCandConn(idOfTransl, -1, self.reacT, self.fromNode, self.statusOfConstruction) #increments the txr line connection variable to node 2
		self.setTranData(cap, life, ROI) #calls setTranData member function to set the parameter values
		#constructor ends

	def __del__(): #destructor
		#log.info("\nThe transmission line object having ID {} have been destroyed.".format(self.translID))
		#end of destructor

	def getTranslID(self): #function gettranslID begins
		return self.translID #returns the ID of the generator object
		#end of gettranslID function

	def getTranslNodeID1(self): #function getGenNodeID begins
		return self.connNodetPtr1.getNodeID() #returns the ID number of the node to which the generator object is connected
		#end of getGenNodeID function

	def getTranslNodeID2(self): #function getGenNodeID begins
		return self.connNodetPtr2.getNodeID() #returns the ID number of the node to which the generator object is connected
		#end of getGenNodeID function

	def getFlowLimit(self): #Function getFlowLimit gets the value of power flow line limit
		return self.ptMax
		#Function getFlowLimit ends

	def getReactance(self):
		return self.reacT

	def setTranData(self, capC, lifeTime, interRate): #member function to set parameter values of transmission lines
		self.capitalCost = capC
		self.lifeYears = lifeTime
		self.rateInterest = interRate
		#end function for setting parameter values

	def getInvestCost(self): #member function getInvestCost begins
		#return (self.capitalCost*self.rateInterest*(((1+self.rateInterest) ** self.lifeYears)))/(((1+self.rateInterest) ** self.lifeYears)-1) #1+self.rateInterest);self.capitalCost/100
		return self.capitalCost

	def returnPresAbsStatus(self): #Returns the construction status of the candidate line
		return self.statusOfConstruction


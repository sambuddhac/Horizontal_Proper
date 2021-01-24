#Member functions for class candLine
#include Node class definition
from Python_src.node import Node
from Python_src.log import log

class candLine(object):
	def __init__(self, localRank, idOfTransl, nodeConnt1, nodeConnt2, PowertMax, Reactance, ROI, life, cap, absPres): #constructor begins
		self.translID = idOfTransl
		self.localIndex = localRank
		self.connNodetPtr1 = nodeConnt1
		self.connNodetPtr2 = nodeConnt2
		self.ptMax = PowertMax
		self.reacT = Reactance
		self.status = absPres
		self.fromNode = connNodetPtr1.getNodeID()
		self.toNode = connNodetPtr2.getNodeID()
		#log.info("\nInitializing the parameters of the transmission line with ID: {}".format(self.translID))
		self.connNodetPtr1.setCandConn(idOfTransl, 1, self.reacT, self.toNode ) #increments the txr line connection variable to node 1
		self.connNodetPtr2.setCandConn(idOfTransl, -1, self.reacT, self.fromNode ) #increments the txr line connection variable to node 1
		self.setTranData(cap, life, ROI) #calls setTranData member function to set the parameter values
		#constructor ends

	#def __del__(): #destructor
		#log.info("\nThe transmission line object having ID {} have been destroyed.\n".format(self.translID))
		#end of destructor

	def getTranslID(self): #function gettranslID begins
		return self.translID #returns the ID of the generator object
		#end of gettranslID function

	def getFromNodeID(self): #function getFromNodeID for the from node ID begins
		return self.connNodetPtr1.getNodeID() #returns the ID number of the from  node 
		#end of getFromNodeID function

	def getToNodeID(self): #function getToNodeID for the to node ID begins
		return self.connNodetPtr2.getNodeID() #returns the ID number of the to node 
		#end of getToNodeID function

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
		#return (self.capitalCost*self.rateInterest*(((1+self.rateInterest) ** self.lifeYears)))/(((1+self.rateInterest) ** self.lifeYears)-1) #(1+self.rateInterest) #self.capitalCost/100
		return self.capitalCost

	def returnPresAbsStatus(self): #Returns the construction status of the candidate line
		return self.status

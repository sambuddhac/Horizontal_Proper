#Member functions for class SELine.
#include Node class definition from node.py
from Python_src.node import Node
from Python_src.log import log

class SELine(object):
	def __init__(self, localRank, idOfTransl, nodeConnt1, nodeConnt2, PowertMax, Reactance): #constructor begins
		self.translID = idOfTransl
		self.localIndex = localRank
		self.connNodetPtr1 = nodeConnt1
		self.connNodetPtr2 = nodeConnt2
		self.ptMax = PowertMax
		self.reacT = Reactance	
		self.fromNode = self.connNodetPtr1.getNodeID()
		self.toNode = self.connNodetPtr2.getNodeID()
		log.info("\nInitializing the parameters of the shared transmission line with ID: {}".format(self.translID))
		log.info("from node: {} To node: {}".format(self.fromNode, self.toNode))
		self.connNodetPtr1.setSEConn(self.translID, 1, self.reacT, self.toNode) #increments the txr line connection variable to node 1
		self.connNodetPtr2.setSEConn(self.translID, -1, self.reacT, self.fromNode) #increments the txr line connection variable to node 1 
		#constructor ends

	#def __del__(): #destructor
		#log.info("\nThe transmission line object having ID {} have been destroyed.\n".format(self.translID))
		# end of destructor

	def getTranslID(self): #function gettranslID begins
		return self.translID #returns the ID of the generator object 
		#end of gettranslID function

	def getFromNodeID(self): #function getFromNodeID begins
		return self.connNodetPtr1.getNodeID() #returns the ID number of the from node 
		#end of getFromNodeID function

	def getToNodeID(self): #function getToNodeID begins
		return self.connNodetPtr2.getNodeID() #returns the ID number of the to node
		# end of getToNodeID function


	def getFlowLimit(self): #Function getFlowLimit gets the value of power flow line limit
		return self.ptMax
		#Function getFlowLimit ends

	def getReactance(self):
		return self.reacT

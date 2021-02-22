#Member functions for class transmissionLine
#include transmissionLine class definition from transl.h, Node class definition from node.h
from Python_src.node import Node
from Python_src.log import log

class transmissionLine(object):
	def __init__(self, idOfTransl, nodeConnt1, nodeConnt2, PowertMax, Reactance): #constructor begins
		self.translID = idOfTransl
		self.connNodet1Ptr = nodeConnt1
		self.connNodet2Ptr = nodeConnt2
		self.ptMax = PowertMax
		self.reacT = Reactance
		self.deviceNature = 0
		#log.info("\nInitializing the parameters of the transmission line with ID: {}".format(self.translID))
		self.fromNode = self.connNodet1Ptr.getNodeID()
		self.toNode = self.connNodet2Ptr.getNodeID()
		self.connNodet1Ptr.settConn(self.idOfTransl, 1, self.reacT, self.toNode) #increments the txr line connection variable to node 1
		self.connNodet2Ptr.settConn(self.idOfTransl, -1, self.reacT, self.fromNode) #increments the txr line connection variable to node 2
		# constructor ends

	#def __del__(): #destructor
		#log.info("\nThe transmission line object having ID {} have been destroyed.\n".format(self.translID))
		#end of destructor

	def getTranslID(self): #function gettranslID begins
		return self.translID #returns the ID of the generator object
		#end of gettranslID function

	def getFlowLimit(self): #function getFlowLimit begins
		return self.ptMax #returns the Maximum power flow limit
		# end of getFlowLimit function

	def getTranslNodeID1(self): #function getGenNodeID begins
		return self.connNodet1Ptr.getNodeID() #returns the ID number of the node to which the generator object is connected 
		#end of getGenNodeID function

	def getTranslNodeID2(self): #function getGenNodeID begins
		return self.connNodet2Ptr.getNodeID() #returns the ID number of the node to which the generator object is connected 
		#end of getGenNodeID function

	def getReactance(self):
		return self.reacT
	


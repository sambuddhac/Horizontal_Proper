#Member functions for class SELine.
#include Node class definition from node.py
from Python_src.node import Node
from Python_src.log import log


class SELine(object):
	def __init__(self, sharedRank, idOfTransl, nodeConnt, fromN, fromZ, toN, toZ, zonalIDNum, PowertMax, Reactance): #constructor begins
		self.translID = idOfTransl
		self.sharedIndex = sharedRank
		self.connNodetPtr = nodeConnt
		self.ptMax = PowertMax
		self.reacT = Reactance
		self.fromZone = fromZ
		self.toZone = toZ
		self.fromNode = fromN
		self.toNode = toN
		self.otherNodeGlobal = 0
		log.info("\nInitializing the parameters of the transmission line with ID: {}".format(self.translID))
		if self.fromZone==zonalIDNum:
			self.connNodetPtr.setSEConn(idOfTransl, 1, self.reacT, self.toZone ) #increments the txr line connection variable to node 1
			self.fromToFlag=1
		else:
			self.connNodetPtr.setSEConn(idOfTransl, -1, self.reacT, self.fromZone ) #increments the txr line connection variable to node 1
			self.fromToFlag=-1
		#constructor ends

	def getTranslID(self): #function gettranslID begins
		return self.translID #returns the ID of the generator object 
		#end of gettranslID function

	def getIntlNodeID(self): #function getGenNodeID begins
		return self.connNodetPtr.getNodeID() #returns the ID number of the node to which the generator object is connected
		#end of getGenNodeID function

	def getIntlZoneID(self): #returns ID number of intra-zonal node end zone to which the transmission line is connected
		if self.fromToFlag==1:
			return self.fromZone #returns the ID number of the from zone if the intra-zonal node is the from node
		else: 
			return self.toZone #returns the ID number of the to zone if the intra-zonal node is the to node 
		#end of getIntlZoneID() function

	def getExtNodeID(): #returns ID number of outer-zonal node end to which the transmission line is connected
{
	if (fromToFlag==1)
		return toNode; // returns the ID number of the from zone if the intra-zonal node is the from node
	else 
		return fromNode; // returns the ID number of the to zone if the intra-zonal node is the to node
} // end of getGenNodeID function

	def getExtZoneID(): #returns ID number of outer-zonal node end zone to which the transmission line is connected
{
	if (fromToFlag==1)
		return toZone; // returns the ID number of the from zone if the intra-zonal node is the from node
	else 
		return fromZone; // returns the ID number of the to zone if the intra-zonal node is the to node
} // end of getGenNodeID function

	def getExtNodeRank(): #function getGenNodeID begins
{
	return otherNodeRank; // returns the ID number of the node to which the generator object is connected
} // end of getGenNodeID function

	def getExtNodeGlobalRank(): #function getExtNodeGlobalRank for the outside-zone node ID begins
{
	return otherNodeGlobal; // returns the global rank of the outside-zone node to which the SE line object is connected
} // end of getExtNodeGlobalRank function

	def getFlowLimit(): #Function getFlowLimit gets the value of power flow line limit	
{
	return ptMax;
} // Function getFlowLimit ends

	def getFlowDir(): #returns the value of the direction flag indicating whether the intra-zonal node end of the line is from (+1) or to (-1) end
{
	return fromToFlag;
} // Function getFlowDir ends

	def outerNodeIndex(int rankOfOuterNode, int dirFlag):
{
	otherNodeRank=rankOfOuterNode;
	fromToOuter=dirFlag;
	connNodetPtr->sendExtNodeInfo( otherNodeRank, fromToOuter, reacT, 0 );
}

	def getReactance():
{
	return reacT;
}
	def assignRank(int ranking): #assigns rank to the from/to node
{
	connNodetPtr->assignGlobalRank(ranking);		
}

	def connectRank(int ranking): #assigns rank to otherNodeGlobal
{
	otherNodeGlobal = ranking;
	connNodetPtr->populateGlobalConn(ranking);		
}

#Member functions for class candLine
#include Node class definition
from Python_src.node import Node
from Python_src.log import log

class candLine(object):
	def __init__(self, sharedRank, idOfTransl, nodeConnt, fromN, fromZ, toN, toZ, zonalIDNum, PowertMax, Reactance, ROI, life, cap, absPres, owner): #constructor begins
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
		self.resFromStageI = absPres
		self.globalRank = 0
		self.ownership = owner
		#log.info("\nInitializing the parameters of the transmission line with ID: {}".format(self.translID))
		if self.fromZone==zonalIDNum:
			self.connNodetPtr.setCandConn(idOfTransl, 1, self.reacT, self.toZone) #increments the txr line connection variable to node 1
			self.fromToFlag=1
		else:
			self.connNodetPtr.setCandConn(idOfTransl, -1, self.reacT, self.fromZone) #increments the txr line connection variable to node 1
			self.fromToFlag=-1
		self.setTranData(cap, life, ROI) #calls setTranData member function to set the parameter values

		#constructor ends

	def returnOwnership(self): #Returns the value of ownership
		return self.ownership

	def modifyNodeReact(self): #function to modify the nodal connected reactance, if the candidate line is actually built
		#If the connected node is the from node
		if self.fromToFlag==1:
			self.connNodetPtr.modifyReactAPP(self.translID, 1, self.reacT, self.otherNodeRank, 1) #increments the txr line connection variable to node 1
			self.connNodetPtr.sendExtNodeInfo(self.otherNodeRank, self.fromToOuter, self.reacT, 0)
		#If the connected node is the to node
		else:
			self.connNodetPtr.modifyReactAPP(self.translID, -1, self.reacT, self.otherNodeRank, 1) #increments the txr line connection variable to node 1
			self.connNodetPtr.sendExtNodeInfo(self.otherNodeRank, self.fromToOuter, self.reacT, 0)
		#function ends

	def getTranslID(self): #function gettranslID begins
		return self.translID #returns the ID of the generator object
		#end of gettranslID function

	def getIntlNodeID(self): #function getGenNodeID for the intra-zone node ID begins
		return self.connNodetPtr.getNodeID() #returns the ID number of the intra-zone node to which the candidate line object is connected
		#end of getGenNodeID function

	def getIntlZoneID(self): #returns ID number of intra-zonal node end zone to which the transmission line is connected
		if self.fromToFlag==1:
			return self.fromZone #returns the ID number of the from zone if the intra-zonal node is the from node
		else:
			return self.toZone #returns the ID number of the to zone if the intra-zonal node is the to node 
	#end of getIntlZoneID() function

	def getExtNodeID(self): #returns ID number of outer-zonal node end to which the transmission line is connected
		if self.fromToFlag==1:
			return self.toNode #returns the ID number of the from zone if the intra-zonal node is the from node
		else: 
			return self.fromNode #returns the ID number of the to zone if the intra-zonal node is the to node
	#end of getGenNodeID function

	def getExtZoneID(self): #returns ID number of outer-zonal node end zone to which the transmission line is connected
		if self.fromToFlag==1:
			return self.toZone #returns the ID number of the from zone if the intra-zonal node is the from node
		else:
			return self.fromZone #returns the ID number of the to zone if the intra-zonal node is the to node
	#end of getGenNodeID function

	def getExtNodeRank(self): #function getGenNodeID for the outside-zone node ID begins
		return self.otherNodeRank #returns the ID number of the outside-zone node to which the candidate line object is connected
	#end of getGenNodeID function

	def getExtNodeGlobalRank(self): #function getExtNodeGlobalRank for the outside-zone node ID begins
		return self.otherNodeGlobal #returns the global rank of the outside-zone node to which the candidate line object is connected
	#end of getExtNodeGlobalRank function

	def getFlowLimit(self): #Function getFlowLimit gets the value of power flow line limit
		return self.ptMax
	#Function getFlowLimit ends

	def getFlowDir(self): #returns the value of the direction flag indicating whether the intra-zonal node end of the line is from (+1) or to (-1) end
		return self.fromToFlag
	#Function getFlowDir ends

	def outerNodeIndex(self, rankOfOuterNode, dirFlag):
		self.otherNodeRank=rankOfOuterNode
		self.fromToOuter=dirFlag
		self.connNodetPtr.sendExtNodeInfo(self.otherNodeRank, self.fromToOuter, self.reacT, 1)

	def getReactance(self):
		return self.reacT

	def assignRank(self, ranking): #assigns rank to the from/to node
		self.connNodetPtr.assignGlobalRank(ranking)

	def connectRank(self, ranking): #assigns rank to otherNodeGlobal
		self.otherNodeGlobal = ranking
		self.connNodetPtr.populateGlobalConn(ranking)

	def getOtherZone(self): #function getOtherZone returns the ID number of the outside zone to which the other node is connected
		if self.fromToFlag==1:
			return self.toZone #returns the ID number of the outside-zone to which the candidate line object is connected
		else:
			return self.fromZone
	#end of getOtherZone function

	def setTranData(self, capC, lifeTime, interRate): #member function to set parameter values of transmission lines
		self.capitalCost=capC
		self.lifeYears=lifeTime
		self.rateInterest=interRate
	#end function for setting parameter values

	def getInvestCost(self): #member function getInvestCost begins
		return (self.capitalCost*self.rateInterest*((1+self.rateInterest)**self.lifeYears))/(((1+self.rateInterest)**self.lifeYears)-1) #(1+rateInterest);capitalCost/100#
		#return capitalCost

	def assignLineRank(self, globRank): #Assigns global rank to the candidate line
		self.globalRank = globRank #Global rank of the candidate line

	def getGlobalRank(self): #Returns the global rank of the candidate line
		return self.globalRank #Global rank of the candidate line

	def returnPresAbsStatus(self): #Returns the construction status of the candidate line
		return self.resFromStageI

	def setPresAbsStatus(self): #Sets the construction status of the candidate line
		self.resFromStageI=1
	#end

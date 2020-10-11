#Member functions for class Node.
from math import *
import numpy as np
#include <iostream>
# include Node class definition from node.h
#include <vector>
#include <algorithm>
#include "node.h"
#using namespace std;

class Node(object):
	def _init_(self, idOfNode, zoneIndex): # constructor begins
		self.nodeID=idOfNode
		self.zoneID=zoneIndex
#cout << "\nInitializing the parameters of the node with ID: " << nodeID << endl;
#initialize the connected devices to zero for node
		self.gConnNumber = 0 # number of generators connected to a particular node
		self.tConnNumber = 0 # number of transmission lines connected to a particular node
		self.lConnNumber = 0 # number of loads connected to a particular node
		self.sharedExConnNumber = 0  # number of shared existing transmission lines connected to this node
		self.builtCandConnNumber = 0 # number of constructed candidate line connected to this node
		self.builtIntCandConnNumber = 0 # number of constructed candidate line connected to this node
		self.candConnNumber = 0 # number of shared candidate transmission lines connected to this node
		self.intCandConnNumber = 0 # number of internal candidate transmission lines connected to this node 
		self.sharedFlag = 0 # node flag to indicate whether a shared existing or candidate line has been connected to a node
		self.PDevCount = 0 # initialize number of devices connectedto a node to zero
		self.fromReact = 0.0 # Initialize the from reactance
		self.toReact = 0.0  # Initialize the to reactance
		self.globalRank = 0 # sets the globalRank to default value of 0 

	# constructor ends

	def getNodeID(self): # function getNodeID begins
		return self.nodeID #returns node ID to the caller
 # end of function getNodeID

	def setgConn(self, serialOfGen):
		self.gConnNumber+=1 # increment the number of generators connected by one whenever a generator is connected to the node
		self.genSerialNum.append(serialOfGen) # records the serial number of the generator connected to the node 
### 

	def settConn(self, tranID, dir, react, rankOfOther):
		self.tConnNumber+=1 # increment the number of txr lines connected by one whenever a txr line is connected to the node
		if dir == 1:
			self.tranFromSerial.append(tranID)
			self.fromReact += 1/react	
			if  rankOfOther in self.connNodeList:  # If predecided Gen value is given for this particular Powergenerator
				pos= self.connNodeList.index(rankOfOther) # find the position of the Powergenerator in the chart of predecided values
				self.connReactRec[pos] -= 1/react
			else: 
				self.connNodeList.append(rankOfOther)
				self.connReactRec.append(-1/react)
	
		else:
			self.tranToSerial.append(tranID)
			self.toReact -= (1/react)
			if rankOfOther in self.connNodeList: # If predecided Gen value is given for this particular Powergenerator
				pos = self.connNodeList.index(rankOfOther) # find the position of the Powergenerator in the chart of predecided values
				self.connReactRec[pos] += 1/react
			else:
				self.connNodeList.append(rankOfOther)
				self.connReactRec.append(1/react)

	def setSEConn(self, tranID, dir, react, connectZone):
		self.sharedExConnNumber +=1  # increment the number of shared existing txr lines connected by one whenever a txr line is connected to the node
		if  dir == 1:
			self.SEFromSerial.append(tranID)
			self.fromReact += (1/react)		
		else:
			self.SEToSerial.append(tranID)
			self.toReact -= (1/react)
		self.connSharedPoint=1 #Flag set to indicate that this node is connected to an SE line
		if connectZone not in self.connectedZoneList: # If the connected zone isn't in the list
			self.connectedZoneList.append(connectZone) # Put it on the list
			self.multiplicity +=1 # increase the multiplicity by 1

	def modifyReactAPP(self, tranID, dir, react, rankOfOther, deviceType): # Modifies the to and from reactances of lines connected to this node, to account for the newly constructed lines
		if deviceType== 1:  # If shared candidate line
			self.builtCandConnNumber+=1 # increment the number of shared constructed candidate lines connected by one whenever the line is connected to the node
			if  dir == 1:  
				self.builtCandFromSerial.append(tranID)
				self.fromReact += (1/react)		
			else:
				self.builtCandToSerial.append(tranID)
				self.toReact -= (1/react)
	
		else: # If internal candidate line
			self.builtIntCandConnNumber+=1 # increment the number of shared constructed candidate lines connected by one whenever the line is connected to the node
			if  dir == 1:
				self.builtIntCandFromSerial.append(tranID)
				self.fromReact += (1/react)	
				if  rankOfOther in  self.connNodeList: # If predecided Gen value is given for this particular Powergenerator
					pos = self.connNodeList.index(rankOfOther)  # find the position of the Powergenerator in the chart of predecided values
					self.connReactRec[pos] -= (1/react)
				else:
					self.connNodeList.append(rankOfOther)
					self.connReactRec.append((-1/react))
			else:
				self.builtIntCandToSerial.append(tranID)
				self.toReact -= (1/react)
				if  rankOfOther in self.connNodeList: # If predecided Gen value is given for this particular Powergenerator
					pos = self.connNodeList.index(rankOfOther) # find the position of the Powergenerator in the chart of predecided values
					self.connReactRec[pos] += (1/react)
				else:
					self.connNodeList.append(rankOfOther)
					self.connReactRec.append((1/react))
		self.connBuiltCandPoint=1

	def setCandConn(self, tranID, dir, react, connectZone):
		self.candConnNumber+=1 # increment the number of shared cand txr lines connected by one whenever a txr line is connected to the node
		if dir == 1:
			self.CandFromSerial.append(tranID)
		else:
			self.CandToSerial.append(tranID)
		self.connCandPoint=1 # Flag set to indicate that this node is connected to a cand line
		if  connectZone not in self.connectedZoneList: # If the connected zone isn't in the list
			self.connectedZoneList.append(connectZone) # Put it on the list
			self.multiplicity+=1 # increase the multiplicity by 1

	def getNodeMultiplicity(self): # get the multiplicity of the node i.e: the number of different zones (other than the one where it belongs) to which it is connected
		return self.multiplicity

	def setIntCandConn(self, tranID, dir, react, rankOfOther, constStat):
		self.intCandConnNumber+=1 # increment the number of shared cand txr lines connected by one whenever a txr line is connected to the node
		if dir == 1:
			self.IntCandFromSerial.append(tranID)
			self.fromReact += constStat*(1/react)	
			if  rankOfOther in self.connNodeList: # If predecided Gen value is given for this particular Powergenerator
				pos = self.connNodeList.index(rankOfOther) # find the position of the Powergenerator in the chart of predecided values
				self.connReactRec[pos] -= constStat*(1/react)
			else:
				self.connNodeList.append(rankOfOther)
				self.connReactRec.append(constStat*(-1/react))
		else:
			self.IntCandToSerial.append(tranID)
			self.toReact -= constStat*(1/react)
			if rankOfOther in self.connNodeList: # If predecided Gen value is given for this particular Powergenerator
				pos = self.connNodeList.index(rankOfOther)  # find the position of the Powergenerator in the chart of predecided values
				self.connReactRec[pos] += constStat*(1/react)
			else:
				self.connNodeList.append(rankOfOther)
				self.connReactRec.append(constStat*(1/react))
		self.connIntCandPoint=1	 # Flag set to indicate that this node is connected to an internal cand line

	def setlConn(self, lID, loadVal):
		self.lConnNumber +=1 # increment the number of loads connected by one whenever a load is connected to the node
		self.loadSerialNum.append(lID)
		self.connLoadVal = []
		self.connLoadVal = loadVal ####

	def getGenLength(self): # function getNodeID begins
		return self.genSerialNum.len() # returns node ID to the caller 

	def getGenSer(self, colCount):
		return self.genSerialNum[colCount-1] ###not sure at is

# function redContNodeCount begins

	def initLoad(self,scenNum): # Initialize the default loads on all nodes to zero
		i = 0 ### for (int i = 0; i < scenNum; ++i)
		for i in range(scenNum):
			self.connLoadVal.append(0) ###Not sure

	def devpinitMessage(self,scenC): # function devpinitMessage begins
		return self.connLoadVal[scenC] # return the total connected load ###Not sure
# function devpinitMessage ends

	def sendExtNodeInfo(self, rankOfOuter, direction,reactance, indicatorSECand): # Function to populate the connected outer-node list
		if rankOfOuter in self.shareNodeList: # If outer node rank is present in the list
			pos = self.shareNodeList.index(rankOfOuter)  # find the position of the outer node
			if direction == 1:
				self.shareReactRec[pos] -= (1-indicatorSECand)*(1/reactance)
			else:
				self.shareReactRec[pos] += (1-indicatorSECand)*(1/reactance)
		else:
			if direction == 1:
				self.shareNodeList.append(rankOfOuter)
				self.shareReactRec.append((1-indicatorSECand)*(-1/reactance))
			else:
				self.shareNodeList.append(rankOfOuter)
				self.shareReactRec.append((1-indicatorSECand)*(1/reactance))

	def getSharedFlag(self):
		return self.connSharedPoint # return the status if this node is connected to a shared existing line		

	def getCandFlag(self):
		return self.connCandPoint  # return the status if this node is connected to a shared cand line

	def getBuiltCandFlag(self):
		return self.connBuiltCandPoint #return the status if this node is connected to a shared cand line that is built

	def getToReact(self):
		return self.toReact  #return the total reciprocal of reactances for which this is the to node

	def getFromReact(self):
		return self.fromReact # return the total reciprocal of reactances for which this is the from node

	def getConNodeLength(self):
		return self.connNodeList.len() # returns the length of the vector containing the connected intra-zonal nodes

	def getConnSer(self,colCount):
		return self.connNodeList[colCount-1] # returns the serial number of the connected internal node at this position

	def getConnReact(self,colCount):
		return self.connReactRec[colCount-1] # returns the serial number of the connected internal node at this position

	def getExtraNodeLength(self):
		return self.shareNodeList.len() # returns the length of the vector containing the connected outer-zonal nodes

	def getExtConnSer(self,colCount):
		return self.shareNodeList[colCount-1] # returns the serial number of the connected external node at this position

	def getExtConnReact(self,colCount):
		return self.shareReactRec[colCount-1] # returns the serial number of the connected internal node at this position

	def getCandLineLengthF(self):
		return self.CandFromSerial.len() # returns the number of cand lines connected to this from node

	def getCandLineLengthT(self):
		return self.CandToSerial.len() # returns the number of cand lines connected to this to node

	def getCandSerF(self,colCount):
		return self.CandFromSerial[colCount-1] #returns the serial number of the cand line at this position

	def getCandSerT(self,colCount):
		return self.CandToSerial[colCount-1]  # returns the serial number of the cand line at this position

	def getIntCandLineLengthF(self):
		return self.IntCandFromSerial.len() # returns the number of cand lines connected to this from node

	def getIntCandLineLengthT(self):
		return self.IntCandToSerial.len()	# returns the number of cand lines connected to this to node

	def getIntCandSerF(self,colCount):
		return self.IntCandFromSerial[colCount-1] # returns the serial number of the cand line at this position

	def getIntCandSerT(self,colCount):
		return self.IntCandToSerial[colCount-1]	# returns the serial number of the cand line at this position

	def assignGlobalRank(self,rank): # Assigns the global rank to the nodes that are ends of shared lines
		self.globalRank = rank # sets the rank 
	def populateGlobalConn(self,rank): # Populates the extNodeGlobalRank vector with the global ranks
		if rank in self.shareNodeList:	#If outer node rank is not present in the list
			self.extNodeGlobalRank.append(rank)
	def getGlobalRank(self): # Returns the global rank of this node
		return self.globalRank # Global rank in the stitched list of shared line end nodes

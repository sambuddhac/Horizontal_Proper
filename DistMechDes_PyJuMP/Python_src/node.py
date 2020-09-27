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
		self.sharedExConnNumber = 0 # number of shared existing transmission lines connected to this node
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
"""
Node::~Node() // destructor
{
	//cout << "\nThe node object having ID " << nodeID << " have been destroyed.\n";

} // end of destructor
"""
	def getNodeID(self): # function getNodeID begins
		return self.nodeID #returns node ID to the caller
 # end of function getNodeID

	def setgConn(self, serialOfGen):
		self.gConnNumber+=1 # increment the number of generators connected by one whenever a generator is connected to the node
		self.genSerialNum.append(serialOfGen) # records the serial number of the generator connected to the node 
### I am not sure about the above codes especially the append command 

	def settConn(self, tranID, dir, react, rankOfOther):
		self.tConnNumber+=1 # increment the number of txr lines connected by one whenever a txr line is connected to the node
		if dir == 1:
			self.tranFromSerial.append(tranID)
			self.fromReact += 1/react
			### I am not sure about the following line	
			if  rankOfOther in self.connNodeList:  # If predecided Gen value is given for this particular Powergenerator
				self.pos= self.connNodeList.index(rankOfOther) # find the position of the Powergenerator in the chart of predecided values
				self.connReactRec[self.pos] -= 1/react
			else: 
				self.connNodeList.append(rankOfOther)
				self.connReactRec.append(-1/react)
	
		else:
			self.tranToSerial.append(tranID)
			self.toReact -= (1/react)
			if rankOfOther in self.connNodeList: # If predecided Gen value is given for this particular Powergenerator
				self.pos= self.connNodeList.index(rankOfOther) # find the position of the Powergenerator in the chart of predecided values
				self.connReactRec[self.pos] += 1/react
			else:
				self.connNodeList.append(rankOfOther)
				self.connReactRec.append(1/react)

	def self.setSEConn(self, tranID, dir, react, connectZone):
		self.sharedExConnNumber += 1  # increment the number of shared existing txr lines connected by one whenever a txr line is connected to the node
		if  dir == 1:
			self.SEFromSerial.append(tranID)
			self.fromReact += (1/react)		
		else:
			self.SEToSerial.append(tranID)
			self.toReact -= (1/react)
		self.connSharedPoint=1 #Flag set to indicate that this node is connected to an SE line
		if connectZone not in self.connectedZoneList: # If the connected zone isn't in the list
			connectedZoneList.append(connectZone) # Put it on the list
			multiplicity +=1 # increase the multiplicity by 1

	def self.modifyReactAPP(self, tranID, dir, react, rankOfOther, deviceType): # Modifies the to and from reactances of lines connected to this node, to account for the newly constructed lines
		if deviceType== 1:  # If shared candidate line
			builtCandConnNumber+=1 # increment the number of shared constructed candidate lines connected by one whenever the line is connected to the node
			if  dir == 1:  
				self.builtCandFromSerial.append(tranID)
				self.fromReact += (1/react)		
			else:
				self.builtCandToSerial.append(tranID)
				self.toReact -= (1/react)
	
		else: # If internal candidate line
			builtIntCandConnNumber+=1 # increment the number of shared constructed candidate lines connected by one whenever the line is connected to the node
			if  dir == 1:
				self.builtIntCandFromSerial.append(tranID)
				self.fromReact += (1/react)	
				if  rankOfOther in  self.connNodeList: # If predecided Gen value is given for this particular Powergenerator
					self.pos = self.connNodeList.index(rankOfOther)  # find the position of the Powergenerator in the chart of predecided values
					self.connReactRec[self.pos] -= (1/react)
				else:
					self.connNodeList.append(rankOfOther)
					self.connReactRec.append((-1/react))
			else:
				self.builtIntCandToSerial.append(tranID)
				self.toReact -= (1/react)
				if  rankOfOther in self.connNodeList: # If predecided Gen value is given for this particular Powergenerator
					self.pos = self.connNodeList.index(rankOfOther) # find the position of the Powergenerator in the chart of predecided values
					self.connReactRec[self.pos] += (1/react)
				else:
					self.connNodeList.append(rankOfOther)
					self.connReactRec.append((1/react))
		connBuiltCandPoint=1

	def setCandConn(self, tranID, dir, react, connectZone):
		self.candConnNumber+=1 # increment the number of shared cand txr lines connected by one whenever a txr line is connected to the node
		if dir == 1:
			self.CandFromSerial.append(tranID)
		else:
			self.CandToSerial.append(tranID)
		self.connCandPoint=1 # Flag set to indicate that this node is connected to a cand line
		if  connectZone not in self.connectedZoneList # If the connected zone isn't in the list
			self.connectedZoneList.append(connectZone) # Put it on the list
			multiplicity+=1 # increase the multiplicity by 1

	def getNodeMultiplicity(self): # get the multiplicity of the node i.e: the number of different zones (other than the one where it belongs) to which it is connected
		return self.multiplicity

	def setIntCandConn(self, tranID, dir, react, rankOfOther, constStat):
		self.intCandConnNumber+=1 # increment the number of shared cand txr lines connected by one whenever a txr line is connected to the node
		if dir == 1:
			self.IntCandFromSerial.append(tranID)
			self.fromReact += constStat*(1/react)	
			if  rankOfOther in self.connNodeList: # If predecided Gen value is given for this particular Powergenerator
				self.pos = self.connNodeList.index(rankOfOther) # find the position of the Powergenerator in the chart of predecided values
				self.connReactRec[self.pos] -= constStat*(1/react)
			else:
				self.connNodeList.append(rankOfOther)
				self.connReactRec.append(constStat*(-1/react))
		else:
			self.IntCandToSerial.append(tranID)
			self.toReact -= constStat*(1/react)
			if rankOfOther in self.connNodeList: # If predecided Gen value is given for this particular Powergenerator
				self.pos = self.connNodeList.index(rankOfOther)  # find the position of the Powergenerator in the chart of predecided values
				self.connReactRec[self.pos] += constStat*(1/react)
			else:
				self.connNodeList.append(rankOfOther)
				self.connReactRec.append(constStat*(1/react))
		self.connIntCandPoint=1; # Flag set to indicate that this node is connected to an internal cand line

	def setlConn(self, lID, loadVal):
		self.lConnNumber+=1 # increment the number of loads connected by one whenever a load is connected to the node
		self.loadSerialNum.append(lID)
		self.connLoadVal.clear(self) ### I am not sure
		self.connLoadVal = *loadVal # total connected load

	def getGenLength(self): # function getNodeID begins
		return self.genSerialNum.len(self) # returns node ID to the caller ### Not sure

	def getGenSer(self, colCount):
		return self.genSerialNum.index(colCount-1) ###not sure

# function redContNodeCount begins

	def initLoad(self,scenNum): # Initialize the default loads on all nodes to zero
		i = 0
		for i in range (scenNum):
			self.connLoadVal.append(0) ###Not sure

	def devpinitMessage(self,scenC): # function devpinitMessage begins
		return self.connLoadVal.index(scenC) # return the total connected load ###Not sure
# function devpinitMessage ends

	def sendExtNodeInfo(self, rankOfOuter, direction,reactance, indicatorSECand): # Function to populate the connected outer-node list
		if rankOfOuter in self.shareNodeList: # If outer node rank is present in the list
			auto pos = std::find(shareNodeList.begin(), shareNodeList.end(), rankOfOuter) - shareNodeList.begin(); // find the position of the outer node
			if (direction == 1) {
				shareReactRec[pos] -= (1-indicatorSECand)*(1/reactance);
		}
			else {
				shareReactRec[pos] += (1-indicatorSECand)*(1/reactance);
		}
	}
		else {
			if (direction == 1) {
				shareNodeList.push_back(rankOfOuter);
				shareReactRec.push_back((1-indicatorSECand)*(-1/reactance));
		}
			else {
				shareNodeList.push_back(rankOfOuter);
				shareReactRec.push_back((1-indicatorSECand)*(1/reactance));
		}

	}	
}

int Node::getSharedFlag()
{
	return connSharedPoint; // return the status if this node is connected to a shared existing line
}		

int Node::getCandFlag()
{
	return connCandPoint; // return the status if this node is connected to a shared cand line
}

int Node::getBuiltCandFlag()
{
	return connBuiltCandPoint; // return the status if this node is connected to a shared cand line that is built
}

double Node::getToReact()
{
	return toReact; // return the total reciprocal of reactances for which this is the to node
}

double Node::getFromReact()
{
	return fromReact; // return the total reciprocal of reactances for which this is the from node
}

int Node::getConNodeLength()
{
	return connNodeList.size(); // returns the length of the vector containing the connected intra-zonal nodes
}

int Node::getConnSer(int colCount)
{
	return connNodeList.at(colCount-1); // returns the serial number of the connected internal node at this position
}

double Node::getConnReact(int colCount)
{
	return connReactRec.at(colCount-1); // returns the serial number of the connected internal node at this position
}

int Node::getExtraNodeLength()
{
	return shareNodeList.size(); // returns the length of the vector containing the connected outer-zonal nodes
}

int Node::getExtConnSer(int colCount)
{
	return shareNodeList.at(colCount-1); // returns the serial number of the connected external node at this position
}

double Node::getExtConnReact(int colCount)
{
	return shareReactRec.at(colCount-1); // returns the serial number of the connected internal node at this position
}

int Node::getCandLineLengthF()
{
	return CandFromSerial.size(); // returns the number of cand lines connected to this from node
}

int Node::getCandLineLengthT()
{
	return CandToSerial.size(); // returns the number of cand lines connected to this to node
}

int Node::getCandSerF(int colCount)
{
	return CandFromSerial.at(colCount-1); // returns the serial number of the cand line at this position
}

int Node::getCandSerT(int colCount)
{
	return CandToSerial.at(colCount-1); // returns the serial number of the cand line at this position
}

int Node::getIntCandLineLengthF()
{
	return IntCandFromSerial.size(); // returns the number of cand lines connected to this from node
}

int Node::getIntCandLineLengthT()
{
	return IntCandToSerial.size(); // returns the number of cand lines connected to this to node
}

int Node::getIntCandSerF(int colCount)
{
	return IntCandFromSerial.at(colCount-1); // returns the serial number of the cand line at this position
}

int Node::getIntCandSerT(int colCount)
{
	return IntCandToSerial.at(colCount-1); // returns the serial number of the cand line at this position
}

void Node::assignGlobalRank( int rank ) // Assigns the global rank to the nodes that are ends of shared lines
{
	globalRank = rank; // sets the rank 
}

void Node::populateGlobalConn( int rank ) // Populates the extNodeGlobalRank vector with the global ranks
{
	if (std::find(extNodeGlobalRank.begin(), extNodeGlobalRank.end(), rank) == shareNodeList.end()) // If outer node rank is not present in the list
		extNodeGlobalRank.push_back(rank);
}

int Node::getGlobalRank() // Returns the global rank of this node
{
	return globalRank; // Global rank in the stitched list of shared line end nodes
}

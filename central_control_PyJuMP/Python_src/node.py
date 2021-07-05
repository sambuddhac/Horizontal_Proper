#Member functions for class Node.
from Python_src.log import log

class Node(object):
	def __init__(self, universalID, idOfNode, zoneIndex): #constructor begins
		self.nodeID = idOfNode #Node object id number
		self.zoneID = zoneIndex #Zone to which the node belongs
		#log.info("\nInitializing the parameters of the node with ID: {}".format(nodeID))
		#initialize the connected devices to zero for node
		self.gConnNumber = 0 #number of generators connected to a particular node
		self.tConnNumber = 0 #number of transmission lines connected to a particular node
		self.lConnNumber = 0 #number of loads connected to a particular node
		self.sharedExConnNumber = 0 #number of shared existing transmission lines connected to this node
		self.builtCandConnNumber = 0 # number of constructed candidate line connected to this node
		self.candConnNumber = 0 # number of shared candidate transmission lines connected to this node
		self.intCandConnNumber = 0 #number of internal candidate transmission lines connected to this node 
		self.sharedFlag = 0 #node flag to indicate whether a shared existing or candidate line has been connected to a node #node flag to indicate whether a shared existing or candidate line has been connected to a node
		self.PDevCount = 0 #initialize number of devices connectedto a node to zero #Number of devices connected to the particular node object
		self.fromReact = 0.0 #Initialize the from reactance #Sum of reciprocals of reactances of lines for which this is the from node
		self.toReact = 0.0 #Initialize the to reactance #Sum of reciprocals of reactances of lines for which this is the to node
		self.globalRank = universalID #sets the globalRank to universalID #Global rank in the stitched list of shared line end nodes
		log.info("Node number: {}".format(globalRank))
		log.info("Zone number: {}".format(zoneID)) 
		self.genSerialNum = [] #vector consisting of the serial numbers of generators connected to a particular node
		self.tranFromSerial = [] #vector consisting of the transmission lines for which the node is from node
		self.connNodeList = [] #List of intra-zonal nodes that are directly connected to this node via transmission lines
		self.connReactRec = [] #List of reciprocals of reactances of the intra zone lines connected to the node
		self.tranToSerial = []	#vector consisting of the transmission lines for which the node is to node
		self.connLoadVal = [] #Scenario values of connected load to this node
		self.loadSerialNum = [] #vector consisting of the serial numbers of loads connected to a particular node
		self.SEFromSerial = [] #vector consisting of the SE lines for which the node is from node
		self.SEToSerial = [] #vector consisting of the SE lines for which the node is to node
		self.shareNodeList = [] #List of outer-zone nodes connected to the node via shared existing transmission lines
		self.shareReactRec = [] #List of reciprocals of reactances of the shared existing lines connected to the node
		self.CandFromSerial = [] #vector consisting of the cand lines for which the node is from node
		self.CandToSerial = [] #vector consisting of the cand lines for which the node is to node
		self.IntCandFromSerial = [] #vector consisting of the internal cand lines for which the node is from node
		self.IntCandToSerial = [] #vector consisting of the internal cand lines for which the node is to node
		self.builtCandFromSerial = [] #vector consisting of the built cand lines for which the node is from node
		self.builtCandToSerial = [] #vector consisting of the built cand lines for which the node is to node	
		self.CandNodeList = [] #List of outer-zone nodes connected to the node via shared cand transmission lines
		self.CandReactRec = [] #List of reciprocals of reactances of the shared cand lines connected to the node
		self.connSharedPoint = 0 #flag to indicate if this node is either the from or to end of any SE line. Default value 0 indicates it isn't
		self.connCandPoint = 0 #flag to indicate if this node is either the from or to end of any Cand line. Default value 0 indicates it isn't
		self.connIntCandPoint = 0 #flag to indicate if this node is either the from or to end of any internal Cand line. Default value 0 indicates it isn't 
		#constructor ends

	def __del__(self): #destructor
		#log.info("\nThe node object having ID {" <<  << "} have been destroyed.\n".format(self.nodeID))
		#end of destructor

	def getNodeID(self): #function getNodeID begins
		return self.globalRank #returns node ID to the caller
		# end of function getNodeID

	def setgConn(self, serialOfGen):
		self.gConnNumber += 1 #increment the number of generators connected by one whenever a generator is connected to the node
		self.genSerialNum.append(serialOfGen) #records the serial number of the generator connected to the node

	def settConn(self, tranID, dir, react, rankOfOther):
		self.tConnNumber += 1 #increment the number of txr lines connected by one whenever a txr line is connected to the node
		if dir == 1:
			self.tranFromSerial.appendk(tranID)
			self.fromReact += 1/react	
			if rankOfOther in self.connNodeList: #If predecided Gen value is given for this particular Powergenerator
				pos = self.connNodeList.index(rankOfOther) #find the position of the Powergenerator in the chart of predecided values
				self.connReactRec[pos] -= 1/react
			else:
				self.connNodeList.append(rankOfOther)
				self.connReactRec.append(-1/react)
		else:
			self.tranToSerial.append(tranID)
			self.toReact -= 1/react
			if rankOfOther in self.connNodeList: #If predecided Gen value is given for this particular Powergenerator
				pos = self.connNodeList.index(rankOfOther) #find the position of the Powergenerator in the chart of predecided values
				self.connReactRec[pos] += 1/react
			else:
				self.connNodeList.append(rankOfOther)
				self.connReactRec.append(1/react)

	def setSEConn(self, tranID, dir, react, rankOfOther):
		self.sharedExConnNumber += 1 #increment the number of shared existing txr lines connected by one whenever a txr line is connected to the node
		self.connSharedPoint = 1 #set the flag to indicate this is the from or to end of a shared existing transmission line
		if dir == 1:
			self.SEFromSerial.append(tranID)
			self.fromReact += 1/react
			if rankOfOther in self.connNodeList: #If predecided Gen value is given for this particular Powergenerator
				pos = self.connNodeList.index(rankOfOther) #find the position of the Powergenerator in the chart of predecided values
				self.connReactRec[pos] -= 1/react
			else:
				self.connNodeList.append(rankOfOther)
				self.connReactRec.append(-1/react)
		else:
			self.SEToSerial.append(tranID)
			self.toReact -= 1/react
			if rankOfOther in self.connNodeList: #If predecided Gen value is given for this particular Powergenerator
				pos = connNodeList.index(rankOfOther) #find the position of the Powergenerator in the chart of predecided values
				self.connReactRec[pos] += 1/react
			else:
				self.connNodeList.append(rankOfOther)
				self.connReactRec.append(1/react)
		log.info("The node {} in zone {} accessed by the SE line {} with reactance {} and dir {}".format(self.globalRank, self.zoneID, tranID, react, dir))

	def setCandConn(self, tranID, dir, react, rankOfOther):
		self.candConnNumber += 1 #increment the number of shared cand txr lines connected by one whenever a txr line is connected to the node
		if dir == 1:
			self.CandFromSerial.append(tranID)
		else:
			self.CandToSerial.append(tranID)
		self.connCandPoint = 1 #Flag set to indicate that this node is connected to a cand line
		log.info("The node {} in zone {} accessed by the shared candidate line {} with reactance {}".format(self.nodeID, self.zoneID, tranID, react))

	def setIntCandConn(self, tranID, dir, react, rankOfOther, constStat):
		self.intCandConnNumber += 1 #increment the number of shared cand txr lines connected by one whenever a txr line is connected to the node
		if dir == 1:
			self.IntCandFromSerial.append(tranID)
			self.fromReact += constStat*(1/react)	
			if rankOfOther in self.connNodeList: #If predecided Gen value is given for this particular Powergenerator
				pos = self.connNodeList.index(rankOfOther) #find the position of the Powergenerator in the chart of predecided values
				self.connReactRec[pos] -= constStat*(1/react)
			else:
				self.connNodeList.append(rankOfOther)
				self.connReactRec.append(constStat*(-1/react))
		else:
			self.IntCandToSerial.append(tranID)
			self.toReact -= constStat*(1/react)
			if rankOfOther in self.connNodeList: #If predecided Gen value is given for this particular Powergenerator
				pos = self.connNodeList.index(rankOfOther) #find the position of the Powergenerator in the chart of predecided values
				self.connReactRec[pos] += constStat*(1/react)
			else:
				self.connNodeList.append(rankOfOther)
				self.connReactRec.append(constStat*(1/react))
		self.connIntCandPoint = 1 #Flag set to indicate that this node is connected to an internal cand line

	def setlConn(self, lID, loadVal):
		self.lConnNumber += 1 #increment the number of loads connected by one whenever a load is connected to the node
		self.loadSerialNum.append(lID)
		self.connLoadVal = []
		self.connLoadVal = loadVal #total connected load

	def getGenLength(self): #function getNodeID begins
		return len(self.genSerialNum) #returns node ID to the caller 
		#end of function getNodeID

	def getGenSer(colCount):
		return self.genSerialNum[colCount-1]

	#function redContNodeCount begins

	def initLoad(self, scenNum): #Initialize the default loads on all nodes to zero
		for i in range(scenNum):
			self.connLoadVal.append(0)

	def devpinitMessage(self, scenC): #function devpinitMessage begins
		return self.connLoadVal[scenC] #return the total connected load 
		#function devpinitMessage ends

	def getSharedFlag(self):
		return self.connSharedPoint #return the status if this node is connected to a shared existing line

	def getCandFlag(self):
		return self.connCandPoint #return the status if this node is connected to a shared cand line

	def getToReact(self):
		return self.toReact #return the total reciprocal of reactances for which this is the to node

	def getFromReact(self):
		return self.fromReact #return the total reciprocal of reactances for which this is the from node

	def getConNodeLength(self):
		return len(self.connNodeList) #returns the length of the vector containing the connected intra-zonal nodes

	def getConnSer(self, colCount):
		return self.connNodeList[colCount-1] #returns the serial number of the connected internal node at this position

	def getConnReact(self, colCount):
		return self.connReactRec[colCount-1] #returns the serial number of the connected internal node at this position

	def getCandLineLengthF(self):
		return len(self.CandFromSerial) #returns the number of cand lines connected to this from node

	def getCandLineLengthT(self):
		return len(self.CandToSerial) #returns the number of cand lines connected to this to node

	def getCandSerF(self, colCount):
		return self.CandFromSerial[colCount-1] #returns the serial number of the cand line at this position

	def getCandSerT(self, colCount):
		return self.CandToSerial[colCount-1] #returns the serial number of the cand line at this position

	def getIntCandLineLengthF(self):
		return len(self.IntCandFromSerial) #returns the number of cand lines connected to this from node

	def getIntCandLineLengthT(self):
		return len(self.IntCandToSerial) #returns the number of cand lines connected to this to node

	def getIntCandSerF(self, colCount):
		return self.IntCandFromSerial[colCount-1] #returns the serial number of the cand line at this position

	def getIntCandSerT(self, colCount):
		return self.IntCandToSerial[colCount-1] #returns the serial number of the cand line at this position


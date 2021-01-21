countOfScenarios = 20
angleDecIndex = []
for scenPos in range(countOfScenarios):
        angleDecIndex.append([])
for scenPos in range(countOfScenarios):
        for i in range(scenPos+10):
                (angleDecIndex[scenPos]).append((i+7)%6)
revisedUpperBound = 0 #total upper bound initialized
#Calculate the MO upper bound, which is the minimum value of the MO objective, calculated with a given mix of the intermediate, shared continuous variable values from different regions/zones
for scenPos in range(countOfScenarios):
	rankPresChecker = [] #Vector to check if a particular rank has been accounted for 
	length = len(angleDecIndex[scenPos]) #length of the rank vector
	print("The Length of the angleDecIndex vector corresponding to scenario index {} is {}".format(scenPos, length))
        #print("The Length of the angleDecIndex vector corresponding to scenario index {} is {}".format(scenPos, length))
	#"""
	for i in range(length): #Put as many zeroes on the rankPresChecker vector as is the length of the rank vector
		rankPresChecker.append(0)
	interAngleTermVec = [] #Intermediate vector for storing the costs of the angle terms from each sub-network
	tracker = 0 #tracker to track the position of the angleDecIndexIterator
	for angleDecIndexIterator in angleDecIndex[scenPos]: #Iterate through rank vector
		interNodeRank = angleDecIndexIterator #Store the value of the rank of the present node iterate in the list
		print("The value of the internode rank, which is the element of the vector is {}".format(interNodeRank))
		if rankPresChecker[tracker] == 0: #If this node rank hasn't been already accounted for
			indexList = [pos1 for pos1 in range(len(angleDecIndex[scenPos])) if angleDecIndex[scenPos][pos1] == interNodeRank] # Get all the indices of interNodeRank in the angleDecIndex vector
			for indexListIterator in indexList: # while all the different positions of this rank hasn't been accounted for
				rankPresChecker[indexListIterator] = 1 #Mark the element in the corresponding position of checker vector 1, to indicate this rank has been accounted for
				print("pos1 is equal to {}".format(indexListIterator))
			indexList = []
		tracker += 1 #Increment the tracker
	print("The length of the rankPresChecker is {}".format(len(rankPresChecker)))
	for rankIterator in rankPresChecker:
		print("The {}-th position element of the rankPresChecker vector is {}".format(rankPresChecker.index(rankIterator), rankIterator))
		#"""

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
using namespace std;
int main()
{

        int countOfScenarios = 20;
        vector<int> angleDecIndex[100];
        for (int scenPos = 0; scenPos < countOfScenarios; ++scenPos) {
                for (int i = 0; i < (scenPos+10); ++i) {
                        (angleDecIndex[scenPos]).push_back((i+7)%6);
                }
        }
	// Sum up the regional upper bounds, which are the tentative regional minima, at the end of every iteration
	double revisedUpperBound = 0; // total upper bound initialized
        /*
	for ( int i = 0; i < numberOfZones; ++i ) {
		revisedUpperBound += regionalUpper[i];
	}*/
	// Calculate the MO upper bound, which is the minimum value of the MO objective, calculated with a given mix of the intermediate, shared continuous variable values from different regions/zones
	for (int scenPos = 0; scenPos < countOfScenarios; ++scenPos) {
		int interNodeRank; // Intermediate variable for storing the rank of the node 
		vector<int> rankPresChecker; // Vector to check if a particular rank has been accounted for 
		int length = (angleDecIndex[scenPos]).size(); // length of the rank vector
                cout << "The Length of the angleDecIndex vector corresponding to scenario index " << scenPos << " is " << length << endl;
		vector<int>::iterator angleDecIndexIterator; // Iterator for the rank vector
		for (int i = 0; i < length; ++i) // Put as many zeroes on the rankPresChecker vector as is the length of the rank vector 
			rankPresChecker.push_back(0);
		vector<double> interAngleTermVec; // Intermediate vector for storing the costs of the angle terms from each sub-network
		int tracker = 0; // tracker to track the position of the angleDecIndexIterator
		for (angleDecIndexIterator=(angleDecIndex[scenPos]).begin(); angleDecIndexIterator!=(angleDecIndex[scenPos]).end(); ++angleDecIndexIterator) { // Iterate through rank vector
			interNodeRank = (*angleDecIndexIterator); // Store the value of the rank of the present node iterate in the list 
                        cout << "The value of the internode rank, which is the element of the vector is " << interNodeRank << endl;
			if ( rankPresChecker.at(tracker) == 0 ) { // If this node rank hasn't been already accounted for
				auto pos = std::find((angleDecIndex[scenPos]).begin(), (angleDecIndex[scenPos]).end(), interNodeRank); // find the first position of this rank in the vector
				while(pos != (angleDecIndex[scenPos]).end()) // while all the different positions of this rank hasn't been accounted for 
    				{
      					auto pos1 = std::distance((angleDecIndex[scenPos]).begin(), pos); // get the location in the vector, of this rank
					rankPresChecker.at(pos1) = 1; // Mark the element in the corresponding position of checker vector 1, to indicate this rank has been accounted for 
                                        cout << "pos1 is equal to " << pos1 << endl;
					//interAngleTermVec.push_back(-((phaseAngleDecision[scenPos]).at(pos1))*(LagMultXi[scenPos*nodeNumber+(*angleDecIndexIterator-1)])); // Calculate cost term
      					pos = std::find(pos + 1, (angleDecIndex[scenPos]).end(), interNodeRank); // Find position of the next occurence of this rank
   				}
				//double smallest_element = *min_element(interAngleTermVec.begin(), interAngleTermVec.end());
				//revisedUpperBound += smallest_element;
				//interAngleTermVec.clear();					
			}
			++tracker; // Increment the tracker
		}
                vector<int>::iterator rankIterator; // Iterator for the rank vector
                for (rankIterator=rankPresChecker.begin(); rankIterator!=rankPresChecker.end(); ++rankIterator) {
                        cout << "The " << std::distance(rankPresChecker.begin(), rankIterator) << "-th position element of the rankPresChecker vector is " << *rankIterator << endl;

                }
	}

        return 0;
}
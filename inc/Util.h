#ifndef _UTIL_H_
#define _UTIL_H_

#define DEBUG_MSG(message) cout << "[DEBUG: " << message << "]" << endl;
#define UNSET_IDX -1
#define P_PER_NODE 16

#include "Coor.h"
// #include "Pin.h"

// #include <iostream>
// #include <unordered_map>
// #include <unordered_set>
// #include <string>
// #include <cassert>
// #include <vector>
// #include <fstream>
#include <cmath>
// #include <climits>
// #include <queue>
// #include <limits.h>
// using namespace std;

double EuclideanDistance(const Coor &p1, const Coor &p2);
double MangattanDistance(const Coor &p1, const Coor &p2);
// string getNewFFName(const string&, const int&); // given prefix string and counter return the concate name = (prefix + string(int))
// double HPWL(const Coor&, const Coor&);
// // double HPWL(const Net& n, vector<Coor>& c, Manager& mgr);
#endif
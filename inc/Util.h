#ifndef _UTIL_H_
#define _UTIL_H_

#define DEBUG_MSG(message) std::cout << "[DEBUG: " << message << "]" << std::endl;
#define UNSET_IDX -1
#define P_PER_NODE 16

// Parameter for KNN
#define MAX_NEIGHBORS 20
// #define MAX_SQUARE_DISPLACEMENT 100000000
// #define MAX_BANDWIDTH 10000
#define BANDWIDTH_SELECTION_NEIGHBOR 4
#define SHIFT_TOLERANCE 0.1
#define SQUARE_EPSILON 25000000

// OpenMP parallel CPU core
#define MAX_THREADS 8

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

class Param{

    public:
        Param();
        ~Param();
        double MAX_SQUARE_DISPLACEMENT;
        double MAX_BANDWIDTH;
};

double SquareEuclideanDistance(const Coor &p1, const Coor &p2);
double MangattanDistance(const Coor &p1, const Coor &p2);
double GaussianKernel(const Coor &p1, const Coor &p2, double bandwidth);
// string getNewFFName(const string&, const int&); // given prefix string and counter return the concate name = (prefix + string(int))
double HPWL(const Coor&, const Coor&);
// // double HPWL(const Net& n, vector<Coor>& c, Manager& mgr);
#endif
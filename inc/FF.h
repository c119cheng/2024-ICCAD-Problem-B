#ifndef _FF_H_
#define _FF_H_

#include <iostream>
#include <string>
#include <unordered_map>
#include <algorithm>
#include "Instance.h"
#include "Manager.h"
#include "Util.h"

class Manager;
class FF;
class Gate;
struct PrevStage
{
    FF* ff; // start point of critical path (FF)
    Gate* outputGate; // start ff's output gate
    std::string pinName; // input pin of outputGate
};

typedef PrevStage NextStage; // ff -> the end of critical path
                             // outputGate -> your output gate for this critical path
                             // outputGate's pinName

class FF : public Instance{
private:
    std::unordered_map<std::string, double> TimingSlack;
    
    // ######################################### used in cluster ########################################################
    int ffIdx;
    int clusterIdx;
    Coor newCoor;
    double bandwidth;   // used in gaussian kernel function
    bool isShifting;
    // pair<other ffId, euclidean distance>, store the neighbor ff with their Id and the distance to this FF
    std::vector<std::pair<int, double>> NeighborFFs;

    // ######################################### used in Preprocessing ########################################################
    PrevStage prevStage; // {prev stage FF/INPUT, {prevFF's output cell on critical path, output cell pin}}
                                                                    // if prev stage FF is nullptr, cur(this) FF is directly connect with prev stage or is IO
                                                                    // use prevInstance
    std::pair<Instance*, std::string> prevInstance; // prev instance on critical path and its putput pin
    std::vector<NextStage> nextStage;
    Coor originalD, originalQ; // initial location for FF list, only can be set in mgr.Debank
    double originalQpinDelay;
public:
    FF();
    ~FF();

    // Setters
    void setTimingSlack(const std::string &pinName, double slack);
    void setFFIdx(int ffIdx);
    void setClusterIdx(int clusterIdx);
    void setNewCoor(const Coor &coor);
    void setBandwidth(const Manager &mgr);
    void addNeighbor(int ffIdx, double euclidean_distance);
    void setIsShifting(bool shift);
    void setPrevStage(PrevStage);
    void setPrevInstance(std::pair<Instance*, std::string>);
    void addNextStage(NextStage);
    void setOriginalCoor(const Coor& coorD, const Coor& coorQ);
    void setOriginalQpinDelay(double);
    
    // Getter
    double getTimingSlack(const std::string &pinName)const;
    int getFFIdx()const;
    bool getIsCluster()const;
    int getClusterIdx()const;
    Coor getNewCoor()const;
    double getBandwidth()const;
    std::pair<int, double> getNeighbor(int idx)const;
    int getNeighborSize()const;
    bool getIsShifting()const;
    PrevStage getPrevStage()const;
    std::pair<Instance*, std::string> getPrevInstance()const;
    std::vector<NextStage> getNextStage()const;
    Coor getOriginalD()const;
    Coor getOriginalQ()const;
    double getOriginalQpinDelay()const;
    
    // ######################################### used in cluster ########################################################
    void sortNeighbors();
    double shift(const Manager &mgr);     // shift the ff and return the euclidean distance from origin coordinate
    // ######################################### used in cluster ########################################################

    friend std::ostream &operator<<(std::ostream &os, const FF &ff);
};


#endif

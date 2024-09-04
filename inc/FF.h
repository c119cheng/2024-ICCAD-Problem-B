#ifndef _FF_H_
#define _FF_H_

#include <iostream>
#include <string>
#include <unordered_map>
#include <queue>
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
    Coor outputGateCoor;
    // Gate* outputGate; // start ff's output gate
    // std::string pinName; // input pin of outputGate
    double originalCost; // of all path, without ff q pin delay and ff->gate hpwl and gate->curFF
    double curCost;
    // double originalHPWL; // from FF -> its output GATE
    // double originalQ;    // ff's original q pin delay
};

struct NextStage
{
    FF* ff; // ff -> the end of critical path
    Gate* outputGate; // outputGate -> your output gate for this critical path
    std::string pinName; // outputGate's pinName
    size_t pathID; // path id in next stage prevFF
};
enum class CellType{
    IO = 0,
    FF = 1,
    GATE = 2
};

struct PrevInstance{
    Instance* instance;
    CellType cellType;
    std::string pinName;
};

class pathComp
{
public:
    bool operator()(const PrevStage* value1, const PrevStage* value2)
    {
        return value1->curCost < value2->curCost;
    }
};
class pathPriorityQueue : public std::priority_queue<PrevStage*, std::vector<PrevStage*>, pathComp>
{
public:
    void update(){ // make heap again
        std::make_heap(this->c.begin(), this->c.end(), this->comp);
    }
};

class FF : public Instance{
public:
    std::unordered_map<std::string, double> TimingSlack;
    std::vector<FF*> clusterFF;
    // ######################################### used in cluster ########################################################
    int ffIdx;
    int clusterIdx;
    Coor newCoor;
    double bandwidth;   // used in gaussian kernel function
    bool isShifting;
    // pair<other ffId, euclidean distance>, store the neighbor ff with their Id and the distance to this FF
    std::vector<std::pair<int, double>> NeighborFFs;
    int clkIdx;
    bool isLegalize;

    // ######################################### used in Preprocessing ########################################################
    std::vector<PrevStage*> prevFF; // if shift regist, it will be empty, other wise need to find cur critical path from here
    PrevInstance prevInstance;     // prev instance on critical path and its output pin
    std::vector<NextStage> nextStage; // include all the path, should check which one is current critical path
    std::vector<NextStage> nextStageCriticalPath;
    bool nextStageCriticalPathDirty;
    Coor originalD, originalQ; // initial location for FF list, only can be set in mgr.Debank
    double originalQpinDelay, originalSlack;
    double originalCriticalPathCost; // largest critical path in input design, for calculate slack

    FF* physicalFF;
    int slot;

    // ######################################### Fixed flag ########################################################
    bool fixed;

    // ######################################### for timing ########################################################
    double curCriticalPathCost;
    PrevStage curCriticalPath;
    pathPriorityQueue pathPQ;

public:
    FF();
    explicit FF(int size);
    ~FF();

    // Setters
    void setTimingSlack(const std::string &pinName, double slack);
    void addClusterFF(FF* inputFF, int slot);
    void setFFIdx(int ffIdx);
    void setClusterIdx(int clusterIdx);
    void setClkIdx(int clkIdx);
    void setNewCoor(const Coor &coor);
    void setBandwidth(const Manager &mgr);
    void addNeighbor(int ffIdx, double euclidean_distance);
    void setIsShifting(bool shift);
    void addPrevStage(PrevStage);
    void setPrevInstance(PrevInstance);
    void addNextStage(NextStage);
    void setOriginalCoor(const Coor& coorD, const Coor& coorQ);
    void setOriginalQpinDelay(double);
    void setOriginalSlack(double slack);
    void setOriginalPathCost(double cost);
    void setPhysicalFF(FF* targetFF, int slot);
    void setClusterSize(int);
    void setFixed(bool fixed);
    void setIsLegalize(bool isLegalize);
    // Getter
    double getTimingSlack(const std::string &pinName)const;
    std::vector<FF*>& getClusterFF();
    int getFFIdx()const;
    bool getIsCluster()const;
    int getClusterIdx()const;
    int getClkIdx()const;
    Coor getNewCoor()const;
    double getBandwidth()const;
    std::pair<int, double> getNeighbor(int idx)const;
    int getNeighborSize()const;
    bool getIsShifting()const;
    const std::vector<PrevStage*>& getPrevStage()const;
    PrevInstance getPrevInstance()const;
    std::vector<NextStage> getNextStage()const;
    std::vector<NextStage> getNextStageCriticalPath();
    Coor getOriginalD()const;
    Coor getOriginalQ()const;
    double getOriginalQpinDelay()const;
    double getOriginalSlack()const;
    double getOriginalPathCost()const;
    FF* getPhysicalFF()const;
    int getSlot()const;
    bool getFixed()const;
    bool getIsLegalize()const;
    std::string getPhysicalPinName();
    std::vector<std::pair<Coor, double>> getCriticalCoor(); // return the relative coor on critical path
    size_t getCriticalSize(); // return the size of all critical path both Q and D pin
    double getAllSlack(); // return the slack of all critical path both Q and D pin (can be positive)
    double getCost(); // return overall cost of MBFF(can be 1 bit), include TNS of Q pin (next stage FFs)
    // ######################################### used in cluster ########################################################
    void sortNeighbors();
    double shift(const std::vector<FF *> &FFs);     // shift the ff and return the euclidean distance from origin coordinate
    // ######################################### used in cluster ########################################################


    // --------------------------- for timing-----------------------
    void getNS(double& TNS, double& WNS); // getNS, getTNS, getWNS will call updateSlack
    double getTNS();
    double getWNS();
    void updateSlack();
    void updateAllCriticalPath(); // update critical path of all FF in cluster
    void updateAllNextCriticalPath();
    bool isCriticalPath(size_t id);
    double getMoveTNS(size_t id, Coor newQCoor); // given the critical path id and it new qpin coor, get the TNS;
    // --------------------------------------------------------------
    
    void clear(); // clear all the data

    friend std::ostream &operator<<(std::ostream &os, const FF &ff);
    friend class postBankingObjFunction;
    friend class DetailPlacement;
    static double DisplacementDelay;
    static double alpha;
    static double beta;
    static double gamma;


    // --------------------------- function below can only be call by FF in clusterFF or FF_list in preprocess-----------------------
    double getSlack(); // don't touch is only for FF in FF_list

private:
    PrevStage getCriticalPath();
    void updateSelfCriticalPath();
    void updateNextStageCriticalPath(); // update all nextStage path of clusterFF
    void updatePathCost(size_t idx); // get the new cost and maintain pq  
};


#endif

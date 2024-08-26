#ifndef _GATE_H_
#define _GATE_H_

#include <iostream>
#include <unordered_map>
#include "Instance.h"

class Gate;
class FF;
struct Path
{
    FF* ff; // start point of path INPUT/FF
    Gate* outputGate; // 
    std::string pinName; // input pin of outputGate
    double cost; // QpinDelay + accumulate HPWL
    double originalHPWL; // original HPWL from ff -> outputGATE
};

class Gate : public Instance{
private:
    // for delay propagation
    int visitedTime;
    std::vector<Path> PathList; // FF coming to this gate
public:
    Gate();
    ~Gate();
    // setters
    void updateVisitedTime();
    void addPath(Path ff);

    // getters
    int getVisitedTime();
    const std::vector<Path>& getPath();

    void removeDuplicatePath(); // remove the path with same start point but smaller cost
    friend std::ostream &operator<<(std::ostream &os, const Gate &gate);
};

#endif

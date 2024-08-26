#include "Gate.h"
#include "FF.h"
Gate::Gate() : visitedTime(0){
}

Gate::~Gate(){}

// setters
void Gate::updateVisitedTime(){
    this->visitedTime++;
}

void Gate::addPath(Path ff){
    PathList.push_back(ff);
}

// getters 
int Gate::getVisitedTime(){
    return this->visitedTime;
}

const std::vector<Path>& Gate::getPath(){
    return this->PathList;
}

void Gate::removeDuplicatePath(){
    std::unordered_map<std::string, Path> m;
    for(size_t i=0;i<PathList.size();i++){
        std::string startPoint;
        if(PathList[i].ff)
            startPoint = PathList[i].ff->getInstanceName() + PathList[i].outputGate->getInstanceName() + PathList[i].pinName;
        else
            startPoint = "IO";
        if(m.count(startPoint)){
            if(PathList[i].cost > m[startPoint].cost)
                m[startPoint] = PathList[i];
        }
        else{
            m[startPoint] = PathList[i];
        }
    }
    PathList = std::vector<Path>(m.size());
    size_t i=0;
    for(auto p : m)
        PathList[i++] = p.second;
}

std::ostream &operator<<(std::ostream &os, const Gate &gate){
    os << "Instance Name: " << gate.instanceName << std::endl;
    os << "Coor: " << gate.coor << std::endl;
    os << "CellName: " << gate.getCell()->getCellName() << std::endl;
    return os;
}
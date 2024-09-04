#include "FF.h"
double FF::DisplacementDelay = 0;
double FF::alpha = 0;
double FF::beta = 0;
double FF::gamma = 0;
FF::FF() : 
    Instance(),
    ffIdx(UNSET_IDX),
    clusterIdx(UNSET_IDX),
    newCoor({0, 0}),
    bandwidth(0),
    isShifting(true),
    clkIdx(UNSET_IDX),
    isLegalize(false),
    prevStage({nullptr, nullptr, ""}),
    prevInstance({nullptr, CellType::IO, ""}),
    originalD({0, 0}),
    originalQ({0, 0}),
    originalQpinDelay(0),
    physicalFF(nullptr),
    slot(0),
    fixed(true){
}

FF::FF(int size) : Instance(), clusterFF(size, nullptr){
    TimingSlack.clear();
    clusterFF.clear();
    ffIdx = UNSET_IDX;
    clusterIdx = UNSET_IDX;
    newCoor = {0, 0};
    coor = {0, 0};
    bandwidth = 0;
    isShifting = true;
    NeighborFFs.clear();
    clkIdx = UNSET_IDX;
    isLegalize = false;
    prevInstance = {nullptr, CellType::IO, ""};
    nextStage.clear();
    originalD = {0, 0};
    originalQ = {0, 0};
    originalQpinDelay = 0;
    physicalFF = nullptr;
    slot = -1;
    fixed = true;
}

FF::~FF(){}

// Setter
void FF::setTimingSlack(const std::string &pinName, double slack){
    TimingSlack[pinName] = slack;
}

void FF::addClusterFF(FF* inputFF, int slot){
    clusterFF[slot] = inputFF;
}

void FF::setFFIdx(int ffIdx){
    this->ffIdx = ffIdx;
}

void FF::setClusterIdx(int clusterIdx){
    this->clusterIdx = clusterIdx;
}

void FF::setClkIdx(int clkIdx){
    this->clkIdx = clkIdx;
}

void FF::setIsLegalize(bool isLegalize){
    this->isLegalize = isLegalize;
}

void FF::setNewCoor(const Coor &coor){
    this->newCoor = coor;
}

void FF::setBandwidth(const Manager &mgr){
    if(getNeighborSize() > BANDWIDTH_SELECTION_NEIGHBOR){
        bandwidth = NeighborFFs[BANDWIDTH_SELECTION_NEIGHBOR].second;
    }
    else{
        bandwidth = NeighborFFs.back().second;
    }

    if(bandwidth > mgr.param.MAX_SQUARE_DISPLACEMENT){
        bandwidth = mgr.param.MAX_BANDWIDTH;
    }
    else{
        bandwidth = std::sqrt(bandwidth);
    }
}

void FF::addNeighbor(int ffIdx, double euclidean_distance){
    NeighborFFs.push_back({ffIdx, euclidean_distance});
}

void FF::setIsShifting(bool shift){
    this->isShifting = shift;
}

void FF::addPrevStage(PrevStage inputStage){
    this->prevFF.push_back(new PrevStage(inputStage));
    pathPQ.push(prevFF[prevFF.size()-1]);
}

void FF::setPrevInstance(const PrevInstance& inputInstance){
    this->prevInstance = inputInstance;
}

void FF::addNextStage(const NextStage& input){
    this->nextStage.push_back(input);
}

/**
 * @brief 
 * 
 * @param targetFF 
 * @param slot 
 */
void FF::setPhysicalFF(FF* targetFF, int slot){
    this->physicalFF = targetFF;
    this->slot = slot;
}

void FF::setClusterSize(int n){
    assert(n == cell->getBits() && "clusterFF size should be equal to cell bit");
    clusterFF = vector<FF*>(n, nullptr);
}

void FF::setFixed(bool fixed){
    this->fixed = fixed;
}

void FF::setOriginalCoor(const Coor& coorD, const Coor& coorQ){
    this->originalD = coorD;
    this->originalQ = coorQ;
}

void FF::setOriginalQpinDelay(double in){
    this->originalQpinDelay = in;
}

void FF::setOriginalSlack(double slack){
    this->originalSlack = slack;
}

void FF::setOriginalPathCost(double cost){
    this->originalCriticalPathCost = cost;
}


// Getter
double FF::getTimingSlack(const std::string &pinName)const{
    auto it = TimingSlack.find(pinName);
    if (it == TimingSlack.end()) {
        throw std::out_of_range("Pin name not found");
    }
    return it->second;
}

std::vector<FF*>& FF::getClusterFF(){
    return this->clusterFF;
}

int FF::getFFIdx()const{
    return ffIdx;
}

bool FF::getIsCluster()const{
    return clusterIdx != UNSET_IDX;
}

int FF::getClusterIdx()const{
    return clusterIdx;
}

int FF::getClkIdx()const{
    return clkIdx;
}

bool FF::getIsLegalize()const{
    return isLegalize;
}

Coor FF::getNewCoor()const{
    return newCoor;
}

double FF::getBandwidth()const{
    return bandwidth;
}

std::pair<int, double> FF::getNeighbor(int idx)const{
    return NeighborFFs[idx];
}

int FF::getNeighborSize()const{
    return NeighborFFs.size();
}

bool FF::getIsShifting()const{
    return isShifting;
}

const vector<PrevStage*>& FF::getPrevStage()const{
    return prevFF;
}

PrevInstance FF::getPrevInstance()const{
    return this->prevInstance;
}

std::vector<NextStage> FF::getNextStage()const{
    return this->nextStage;
}

std::vector<NextStage> FF::getNextStageCriticalPath(){
    if(nextStageCriticalPathDirty){
        nextStageCriticalPath.clear();
        nextStageCriticalPathDirty = false;
        for(size_t i=0;i<nextStage.size();i++){
            NextStage next = nextStage[i];
            if(next.outputGate == nullptr || next.ff->isCriticalPath(next.pathID)){
                nextStageCriticalPath.push_back(next);
            }
        }
    }
    return nextStageCriticalPath;
}

Coor FF::getOriginalD()const{
    return originalD;
}

Coor FF::getOriginalQ()const{
    return originalQ;
}

double FF::getOriginalQpinDelay()const{
    return originalQpinDelay;
}

double FF::getOriginalSlack()const{
    return originalSlack;
}

double FF::getOriginalPathCost()const{
    return originalCriticalPathCost;
}

FF* FF::getPhysicalFF()const{
    return physicalFF;
}

int FF::getSlot()const{
    return slot;
}

bool FF::getFixed()const{
    return fixed;
}

void FF::sortNeighbors(){
    auto FFcmp = [](const std::pair<int, double> &neighbor1, const std::pair<int, double> &neighbor2){
        return neighbor1.second < neighbor2.second;
    };
    std::sort(NeighborFFs.begin(), NeighborFFs.end(), FFcmp);
}

double FF::shift(const std::vector<FF *> &FFs){
    double x_shift = 0;
    double y_shift = 0;
    double scale_factor = 0;
    for(size_t i = 0; i < NeighborFFs.size(); i++){
        FF *ffneighbor = FFs[NeighborFFs[i].first];
        double bandwidth_i = ffneighbor->getBandwidth();
        double weight = GaussianKernel(newCoor, ffneighbor->getCoor(), bandwidth_i) / std::pow(bandwidth_i, 4);
        // DEBUG_MSG(weight)
        x_shift += ffneighbor->getCoor().x * weight;
        y_shift += ffneighbor->getCoor().y * weight;
        scale_factor += weight;
    }
    assert(std::isnormal(scale_factor));
    x_shift /= scale_factor;
    y_shift /= scale_factor;
    double euclidean_distance = std::sqrt(SquareEuclideanDistance({x_shift, y_shift}, newCoor));
    setNewCoor({x_shift, y_shift});
    return euclidean_distance;
}

void FF::getNS(double& TNS, double& WNS){
    this->updateSlack();
    TNS = 0;
    WNS = 0;
    for(const auto& slack_m : TimingSlack){
        if(slack_m.second < 0){
            TNS += slack_m.second;
            WNS = std::min(slack_m.second, WNS);
        }
    }
    TNS = std::abs(TNS);
    WNS = std::abs(WNS);
}

double FF::getTNS(){
    this->updateSlack();
    double TNS = 0;
    for(const auto& slack_m : TimingSlack){
        if(slack_m.second < 0){
            TNS += slack_m.second;
        }
    }
    TNS = std::abs(TNS);
    return TNS;
}

double FF::getWNS(){
    this->updateSlack();
    double WNS = 0;
    for(const auto& slack_m : TimingSlack){
        if(slack_m.second < 0){
            WNS = std::min(slack_m.second, WNS);
        }
    }
    WNS = std::abs(WNS);
    return WNS;
}

void FF::updateSlack(){
    int curSlot = 0;
    for(auto& FF : clusterFF){
        double slack = FF->getSlack();
        if(clusterFF.size() == 1){
            TimingSlack["D"] = slack;
        }
        else{
            TimingSlack["D" + std::to_string(curSlot)] = slack;
        }
        curSlot++;
    }
}

void FF::clear(){
    TimingSlack.clear();
    clusterFF.clear();
    NeighborFFs.clear();
    prevFF.clear();
    prevInstance = {nullptr, CellType::IO, ""};
    nextStage.clear();
    physicalFF = nullptr;
}

std::ostream &operator<<(std::ostream &os, const FF &ff){
    os << "Instance Name: " << ff.instanceName << std::endl;
    os << "Coor: " << ff.coor << std::endl;
    os << "CellName: " << ff.getCell()->getCellName() << std::endl;
    for(auto &pair : ff.TimingSlack)
        os << "Pin[" << pair.first << "] Slack: " << pair.second << std::endl;
    return os;
}


/**
 * @brief 
 * 
 * @return the current critical path, also set in FF
*/
PrevStage FF::getCriticalPath(){
    curCriticalPath = *pathPQ.top();
    curCriticalPathCost = curCriticalPath.curCost;
    return curCriticalPath;
}

/**
 * @brief 
 * 
 * @return double slack, the value of the slack, positive is positive, negative is negative
 */
double FF::getSlack(){
    FF* cur_ff = this;
    // update slack for new location
    Coor newCoorD = this->physicalFF->getNewCoor() + this->physicalFF->getPinCoor("D" + this->getPhysicalPinName());
    double delta_hpwl = 0;
    Coor inputCoor;
    // D pin delta HPWL
    PrevInstance prevInstance = cur_ff->getPrevInstance();
    if(prevInstance.instance){
        if(prevInstance.cellType == CellType::IO){
            inputCoor = prevInstance.instance->getCoor();
            double old_hpwl = HPWL(inputCoor, cur_ff->getOriginalD());
            double new_hpwl = HPWL(inputCoor, newCoorD);
            delta_hpwl += old_hpwl - new_hpwl;
        }
        else if(prevInstance.cellType == CellType::GATE){
            inputCoor = prevInstance.instance->getCoor() + prevInstance.instance->getPinCoor(prevInstance.pinName);
            double old_hpwl = HPWL(inputCoor, cur_ff->getOriginalD());
            double new_hpwl = HPWL(inputCoor, newCoorD);
            delta_hpwl += old_hpwl - new_hpwl;
        }
        else{
            FF* inputFF = dynamic_cast<FF*>(prevInstance.instance);
            inputCoor = inputFF->getOriginalQ();
            Coor newCoorQ = inputFF->physicalFF->getNewCoor() + inputFF->physicalFF->getPinCoor("Q" + inputFF->getPhysicalPinName());
            double old_hpwl = HPWL(inputCoor, cur_ff->getOriginalD());
            double new_hpwl = HPWL(newCoorQ, newCoorD);
            delta_hpwl += old_hpwl - new_hpwl;
        }
    }

    // Q pin delta HPWL (prev stage FFs Qpin)
    double prevCost = 0;
    if(prevFF.size()!=0){
        this->getCriticalPath();
        prevCost = originalCriticalPathCost - curCriticalPathCost;
    }
    // get new slack
    double newSlack = originalSlack + prevCost + FF::DisplacementDelay * delta_hpwl;
    return newSlack;
}

std::string FF::getPhysicalPinName(){
    if(this->physicalFF->getCell()->getBits() == 1)
        return "";
    else
        return std::to_string(this->slot); 
}

std::vector<std::pair<Coor, double>> FF::getCriticalCoor(){
    std::vector<std::pair<Coor, double>> coorList;
    for(auto& curFF : clusterFF){
        // D pin
        Coor inputCoor;
        PrevInstance prev = curFF->getPrevInstance();
        if(prev.instance){
            if(prev.cellType == CellType::IO){
                inputCoor = prev.instance->getCoor();
            }
            else if(prev.cellType == CellType::GATE){
                inputCoor = prev.instance->getCoor() + prev.instance->getPinCoor(prev.pinName);
            }
            else{
                FF* inputFF = dynamic_cast<FF*>(prev.instance);
                inputCoor = inputFF->physicalFF->getNewCoor() + inputFF->physicalFF->getPinCoor("Q" + inputFF->getPhysicalPinName());
            }
        }
        coorList.push_back({inputCoor, curFF->getSlack()});

        // Q pin
        for(auto& next : curFF->getNextStageCriticalPath()){
            Coor outputCoor;
            if(next.outputGate){
                outputCoor = next.outputGate->getCoor() + next.outputGate->getPinCoor(next.pinName);
            }
            else{
                outputCoor = next.ff->physicalFF->getNewCoor() + next.ff->physicalFF->getPinCoor("D" + next.ff->getPhysicalPinName());
            }
            coorList.push_back({outputCoor, next.ff->getSlack()});
        }
    }
    return coorList;
}

size_t FF::getCriticalSize(){
    size_t num = 0;
    for(auto& curFF : clusterFF)
        num += 1 + curFF->getNextStageCriticalPath().size();
    return num;
}

double FF::getAllSlack(){
    double slack = 0;
    for(auto& curFF : clusterFF){
        slack += curFF->getSlack();
        for(auto& next : curFF->getNextStageCriticalPath())
            slack += next.ff->getSlack();
    }
    return slack;
}

double FF::getCost(){
    double timingCost = 0;
    double areaCost = cell->getArea();
    double powerCost = cell->getGatePower();

    for(auto& curFF : clusterFF){
        double slack = curFF->getSlack();
        timingCost += slack < 0 ? -slack : 0;
        for(auto& next : curFF->getNextStageCriticalPath()){
            if(next.outputGate == nullptr || next.ff->isCriticalPath(next.pathID)){ // only consider slack of shift register or critical path
                slack = next.ff->getSlack();
                timingCost += slack < 0 ? -slack : 0;
            }
        }
    }
    return FF::alpha * timingCost + FF::beta * powerCost + FF::gamma * areaCost;
}

bool FF::isCriticalPath(size_t id){
    return prevFF[id] == pathPQ.top();
}

void FF::updateAllCriticalPath(){
    #pragma omp parallel for num_threads(MAX_THREADS)
    for(size_t i=0;i<clusterFF.size();i++){
        FF* curFF = clusterFF[i];
        curFF->updateSelfCriticalPath();
    }
    nextStageCriticalPathDirty = true;
}

void FF::updateAllNextCriticalPath(){
    #pragma omp parallel for num_threads(MAX_THREADS)
    for(size_t i=0;i<clusterFF.size();i++){
        FF* curFF = clusterFF[i];
        curFF->updateNextStageCriticalPath();
    }
    nextStageCriticalPathDirty = true;
}

void FF::updateSelfCriticalPath(){
    FF* curFF = this;
    bool dirty = false;
    #pragma omp parallel for num_threads(MAX_THREADS)
    for(size_t idx=0;idx<curFF->prevFF.size();idx++){
        PrevStage* p = prevFF[idx];
        if(p->ff == nullptr){ // IO -> FF, cost never change
            continue;
        }
        double newCost = p->originalCost 
            + FF::DisplacementDelay * HPWL(p->ff->getPhysicalFF()->getNewCoor() + p->ff->getPhysicalFF()->getCell()->getPinCoor("Q" + p->ff->getPhysicalPinName()),
                        p->outputGateCoor)
            + p->ff->getPhysicalFF()->getCell()->getQpinDelay();
        if(newCost != p->curCost){
            p->curCost = newCost;
            if(newCost > pathPQ.top()->curCost || p == pathPQ.top()) // only update pq when is critical path or its cost is larger
                dirty = true;
        }
    }
    if(dirty)
        curFF->pathPQ.update();
}

void FF::updateNextStageCriticalPath(){
    #pragma omp parallel for num_threads(MAX_THREADS)
    for(size_t i=0;i<clusterFF.size();i++){
        FF* curFF = clusterFF[i];
        for(size_t j=0;j<curFF->nextStage.size();j++){
            if(curFF->nextStage[j].outputGate != nullptr) // is not shift register
                curFF->nextStage[j].ff->updatePathCost(curFF->nextStage[j].pathID);
        }
    }
}

void FF::updatePathCost(size_t idx){
    PrevStage* p = prevFF[idx];
    double newCost = p->originalCost 
        + FF::DisplacementDelay * HPWL(p->ff->getPhysicalFF()->getNewCoor() + p->ff->getPhysicalFF()->getCell()->getPinCoor("Q" + p->ff->getPhysicalPinName()),
                    p->outputGateCoor)
        + p->ff->getPhysicalFF()->getCell()->getQpinDelay();
    if(newCost != p->curCost){
        p->curCost = newCost;
        if(newCost > pathPQ.top()->curCost || p == pathPQ.top()) // only update pq when is critical path or its cost is larger
            pathPQ.update();
    }
}

double FF::getMoveTNS(size_t id, Coor newQCoor){
    PrevStage* p = prevFF[id];
    double newCost = p->originalCost 
        + FF::DisplacementDelay * HPWL(newQCoor, p->outputGateCoor)
        + p->ff->getPhysicalFF()->getCell()->getQpinDelay();
    double slack;

    FF* cur_ff = this;
    // update slack for new location
    Coor newCoorD = this->physicalFF->getNewCoor() + this->physicalFF->getPinCoor("D" + this->getPhysicalPinName());
    double delta_hpwl = 0;
    Coor inputCoor;
    // D pin delta HPWL
    PrevInstance prevInstance = cur_ff->getPrevInstance();
    if(prevInstance.instance){
        if(prevInstance.cellType == CellType::IO){
            inputCoor = prevInstance.instance->getCoor();
            double old_hpwl = HPWL(inputCoor, cur_ff->getOriginalD());
            double new_hpwl = HPWL(inputCoor, newCoorD);
            delta_hpwl += old_hpwl - new_hpwl;
        }
        else if(prevInstance.cellType == CellType::GATE){
            inputCoor = prevInstance.instance->getCoor() + prevInstance.instance->getPinCoor(prevInstance.pinName);
            double old_hpwl = HPWL(inputCoor, cur_ff->getOriginalD());
            double new_hpwl = HPWL(inputCoor, newCoorD);
            delta_hpwl += old_hpwl - new_hpwl;
        }
        else{
            FF* inputFF = dynamic_cast<FF*>(prevInstance.instance);
            inputCoor = inputFF->getOriginalQ();
            Coor newCoorQ = inputFF->physicalFF->getNewCoor() + inputFF->physicalFF->getPinCoor("Q" + inputFF->getPhysicalPinName());
            double old_hpwl = HPWL(inputCoor, cur_ff->getOriginalD());
            double new_hpwl = HPWL(newCoorQ, newCoorD);
            delta_hpwl += old_hpwl - new_hpwl;
        }
    }

    if(newCost > pathPQ.top()->curCost)
        slack = originalSlack + originalCriticalPathCost - newCost + FF::DisplacementDelay * delta_hpwl;
    else
        slack = 0;
    return slack < 0 ? -slack : 0;
}
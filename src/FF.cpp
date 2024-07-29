#include "FF.h"
double FF::DisplacementDelay = 0;
FF::FF() : Instance(){
    ffIdx = UNSET_IDX;
    clusterIdx = UNSET_IDX;
    coor = {0, 0};
    bandwidth = 0;
    isShifting = true;
    clkIdx = UNSET_IDX;
    prevStage = {nullptr, nullptr, ""};
    prevInstance = {nullptr, CellType::IO, ""};
}

FF::FF(int size) : Instance(), clusterFF(size, nullptr){
    ffIdx = UNSET_IDX;
    clusterIdx = UNSET_IDX;
    coor = {0, 0};
    bandwidth = 0;
    isShifting = true;
    clkIdx = UNSET_IDX;
    prevStage = {nullptr, nullptr, ""};
    prevInstance = {nullptr, CellType::IO, ""};
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

void FF::setPrevStage(PrevStage inputStage){
    this->prevStage = inputStage;
}

void FF::setPrevInstance(PrevInstance inputInstance){
    this->prevInstance = inputInstance;
}

void FF::addNextStage(NextStage input){
    this->nextStage.push_back(input);
}

void FF::setPhysicalFF(FF* targetFF, int slot){
    this->physicalFF = targetFF;
    this->slot = slot;
}

void FF::setClusterSize(int n){
    assert(n == cell->getBits() && "clusterFF size should be equal to cell bit");
    clusterFF = vector<FF*>(n, nullptr);
}

void FF::setOriginalCoor(const Coor& coorD, const Coor& coorQ){
    this->originalD = coorD;
    this->originalQ = coorQ;
}

void FF::setOriginalQpinDelay(double in){
    this->originalQpinDelay = in;
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

PrevStage FF::getPrevStage()const{
    return this->prevStage;
}

PrevInstance FF::getPrevInstance()const{
    return this->prevInstance;
}

std::vector<NextStage> FF::getNextStage()const{
    return this->nextStage;
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

FF* FF::getPhysicalFF()const{
    return physicalFF;
}

int FF::getSlot()const{
    return slot;
}

void FF::sortNeighbors(){
    auto FFcmp = [](const std::pair<int, double> &neighbor1, const std::pair<int, double> &neighbor2){
        return neighbor1.second < neighbor2.second;
    };
    std::sort(NeighborFFs.begin(), NeighborFFs.end(), FFcmp);
}

double FF::shift(std::vector<FF *> &FFs){
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
    prevStage = {nullptr, nullptr, ""};
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


double FF::getSlack(){
    FF* cur_ff = this;
    // update slack for new location
    Coor newCoorD = this->physicalFF->getNewCoor() + this->physicalFF->getPinCoor("D" + this->getPhysicalPinName());
    double delta_hpwl = 0;
    double delta_q = 0; // delta q pin delay
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
    const PrevStage& prev = cur_ff->getPrevStage();
    if(prev.ff){
        Coor originalInput, newInput;

        FF* inputFF = prev.ff;
        originalInput = inputFF->getOriginalQ();
        newInput = inputFF->physicalFF->getNewCoor() + inputFF->physicalFF->getPinCoor("Q" + inputFF->getPhysicalPinName());
        delta_q = inputFF->getOriginalQpinDelay() - inputFF->physicalFF->getCell()->getQpinDelay();
        
        inputCoor = prev.outputGate->getCoor() + prev.outputGate->getPinCoor(prev.pinName);
        double old_hpwl = HPWL(inputCoor, originalInput);
        double new_hpwl = HPWL(inputCoor, newInput);
        delta_hpwl += old_hpwl - new_hpwl;
    }
    // get new slack
    double newSlack = cur_ff->getTimingSlack("D") + (delta_q) + FF::DisplacementDelay * delta_hpwl;
    return newSlack;
}

std::string FF::getPhysicalPinName(){
    if(this->physicalFF->getCell()->getBits() == 1)
        return "";
    else
        return std::to_string(this->slot); 
}

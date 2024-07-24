#include "Legalizer.h"

Legalizer::Legalizer(Manager& mgr) : mgr(mgr){
    timer.start();
}

Legalizer::~Legalizer(){
    for(auto &ff : ffs){
        delete ff;
    }
    for(auto &gate : gates){
        delete gate;
    }
    for(auto &row : rows){
        delete row;
    }
}

void Legalizer::run(){
    LoadFF();
    LoadGate();
    LoadPlacementRow();
    SliceRowsByRows();
    SliceRowsByGate();
    Tetris();
    LegalizeWriteBack();
    timer.stop();
}

void Legalizer::LoadFF(){
    DEBUG_LGZ("Load FF to Databse");
    for(const auto &pair : mgr.FF_Map){
        Node *ff = new Node();
        ff->setName(pair.second->getInstanceName());
        ff->setGPCoor(pair.second->getNewCoor());
        ff->setLGCoor(Coor(DBL_MAX, DBL_MAX));
        ff->setCell(pair.second->getCell());
        ff->setW(pair.second->getW());
        ff->setH(pair.second->getH());
        ff->setIsPlace(false);
        ffs.emplace_back(ff);
    }
}

void Legalizer::LoadGate(){
    DEBUG_LGZ("Load Gate to Databse");
    for(const auto &pair: mgr.Gate_Map){
        Node *gate = new Node();
        gate->setName(pair.second->getInstanceName());
        gate->setGPCoor(pair.second->getCoor());
        gate->setLGCoor(pair.second->getCoor());
        gate->setCell(nullptr);
        gate->setW(pair.second->getW());
        gate->setH(pair.second->getH());
        gate->setIsPlace(false);
        gates.emplace_back(gate);
    }
}

void Legalizer::LoadPlacementRow(){
    DEBUG_LGZ("Load Placement Row to Databse");
    std::vector<PlacementRow> PlacementRows = mgr.die.getPlacementRows();
    for(size_t i = 0; i < PlacementRows.size(); i++){
        Row *row = new Row();
        row->setStartCoor(PlacementRows[i].startCoor);
        row->setSiteHeight(mgr.die.getDieBorder().y - PlacementRows[i].startCoor.y);
        row->setSiteWidth(PlacementRows[i].siteWidth);
        row->setNumOfSite(PlacementRows[i].NumOfSites);
        row->setEndX(PlacementRows[i].startCoor.x + PlacementRows[i].siteWidth * PlacementRows[i].NumOfSites);

        // init Row::subrows
        Subrow *subrow = new Subrow();
        subrow->setStartX(PlacementRows[i].startCoor.x);
        subrow->setEndX(PlacementRows[i].startCoor.x + PlacementRows[i].siteWidth * PlacementRows[i].NumOfSites);
        subrow->setFreeWidth(PlacementRows[i].siteWidth * PlacementRows[i].NumOfSites);
        subrow->setHeight(mgr.die.getDieBorder().y - PlacementRows[i].startCoor.y);
        row->addSubrows(subrow);
        rows.emplace_back(row);
    }

    // sort row by the y coordinate in ascending order, if tie, sort by x in ascending order
    std::sort(rows.begin(), rows.end(), [](const Row *a, const Row *b){
        return *a < *b;
    });
}


void Legalizer::SliceRowsByRows(){
    for(size_t i = 0; i < rows.size() - 1; i++){
        double startX = rows[i]->getStartCoor().x;
        double endX = rows[i]->getEndX();
        double siteHeight = rows[i]->getSiteHeight();
        std::list<XTour> xList;
        for(size_t j = i + 1; j < rows.size(); j++){
            if(rows[j]->getEndX() <= startX || rows[j]->getStartCoor().x >= endX) continue;
            Node * upRow = new Node();
            upRow->setGPCoor(rows[j]->getStartCoor());
            upRow->setLGCoor(rows[j]->getStartCoor());
            upRow->setW(rows[j]->getEndX() - rows[j]->getStartCoor().x);
            upRow->setH(rows[j]->getSiteHeight());
            rows[i]->slicing(upRow);
            UpdateXList(rows[j]->getStartCoor().x, rows[j]->getEndX(), xList);

            if((*xList.begin()).startX <= startX && (*xList.begin()).endX >= endX){
                siteHeight = rows[j]->getStartCoor().y - rows[i]->getStartCoor().y;
                break;
            }
        }
        rows[i]->setSiteHeight(siteHeight);
    }
}

void Legalizer::SliceRowsByGate(){
    DEBUG_LGZ("Seperate PlacementRows by Gate Cell");
    for(const auto &gate : gates){
        for(auto &row : rows){
            if(row->getStartCoor().y > gate->getGPCoor().y + gate->getH()) break;
            row->slicing(gate);
        }
        gate->setIsPlace(true);
    }
}


void Legalizer::Tetris(){
    DEBUG_LGZ("Start Legalize FF");
    std::sort(ffs.begin(), ffs.end(), [](const Node *a, const Node *b){
        double costA = a->getH() * a->getW();
        double costB = b->getH() * b->getW();
        return costA > costB;
    });

    for(const auto &ff : ffs){
        size_t closest_row_idx = FindClosestRow(ff);
        double minDisplacement = PlaceFF(ff, closest_row_idx);
        int down_row_idx = closest_row_idx - 1;
        int up_row_idx = closest_row_idx + 1;

        // local search down
        while(down_row_idx >= 0 && std::abs(ff->getGPCoor().y - rows[down_row_idx]->getStartCoor().y) < minDisplacement){
            double downDisplacement = PlaceFF(ff, down_row_idx);
            minDisplacement = minDisplacement < downDisplacement ? minDisplacement : downDisplacement;
            down_row_idx--;
        }

        // local search up
        while(up_row_idx < (int)rows.size() && std::abs(ff->getGPCoor().y - rows[up_row_idx]->getStartCoor().y) < minDisplacement){
            double upDisplacement = PlaceFF(ff, up_row_idx);
            minDisplacement = minDisplacement < upDisplacement ? minDisplacement : upDisplacement;
            up_row_idx++;
        }

        if(ff->getIsPlace()){
            for(auto &row : rows){
                if(row->getStartCoor().y > ff->getLGCoor().y + ff->getH()) break;
                row->slicing(ff);
            }
            continue;
        }

        // change cell type (in top 3) and local search
        DEBUG_LGZ("Legalize Failed");
    }
}

void Legalizer::LegalizeWriteBack(){
    DEBUG_LGZ("Write Back Legalize Coordinate");
    for(const auto &ff : ffs){
        if(ff->getIsPlace()){
            mgr.FF_Map[ff->getName()]->setNewCoor(ff->getLGCoor());
            //std::cout << *ff << std::endl;
        }
        else{
            mgr.FF_Map[ff->getName()]->setNewCoor(Coor(0, 0));
        }
    }
}

// Helper Function
void Legalizer::UpdateXList(double start, double end, std::list<XTour> & xList){
    if(xList.empty()){
        XTour insertRow;
        insertRow.startX = start;
        insertRow.endX = end;
        xList.push_back(insertRow);
        return;
    }

    auto iter1 = xList.begin();
    bool startIsFound = false;
    while(iter1 != xList.end() && start >= (*iter1).startX){
        if(start >= (*iter1).startX && start <= (*iter1).endX){
            startIsFound = true;
            break;
        } 
        iter1++;
    }
    auto iter2 = iter1;
    bool endIsFound = false;
    while(iter2 != xList.end() && end >= (*iter2).startX){
        if(end >= (*iter2).startX && end <= (*iter2).endX){
            endIsFound = true;
            break;
        }
        iter2++;
    }
    if(!startIsFound && !endIsFound){
        std::cout << "All Not Found!" << std::endl;
        XTour insertRow = {start, end};
        if(iter1 != iter2){
            auto iter = xList.erase(iter1, iter2);
            xList.insert(iter, insertRow);
        }else{
            xList.insert(iter1, insertRow);
        }
    }else{
        if(iter1 == iter2){
            (*iter1).endX = ((*iter1).endX >= end) ? (*iter1).endX : end;
        }else{
            XTour insertRow = {start, end};
            std::list<XTour>::iterator eraseEnd;
            if(endIsFound){
                eraseEnd = next(iter2);
            }else{
                eraseEnd = iter2;
            }
            auto iter = xList.erase(iter1, eraseEnd);
            xList.insert(iter, insertRow);
        }
    }
}

size_t Legalizer::FindClosestRow(Node *ff){
    double min_distance = std::abs(ff->getGPCoor().y - rows[0]->getStartCoor().y);
    for(size_t i = 1; i < rows.size(); i++){
        double curr_distance = std::abs(ff->getGPCoor().y - rows[i]->getStartCoor().y);
        if(curr_distance < min_distance){
            min_distance = curr_distance;
        }
        else{
            return i - 1;
        }
    }
    return rows.size() - 1;
}

double Legalizer::PlaceFF(Node *ff, size_t row_idx){
    double minDisplacement = ff->getDisplacement();
    const auto &subrows = rows[row_idx]->getSubrows();
    for(size_t i = 0; i < subrows.size(); i++){
        const auto &subrow = subrows[i];
        double alignedStartX = rows[row_idx]->getStartCoor().x + std::ceil((int)(subrow->getStartX() - rows[row_idx]->getStartCoor().x) / rows[row_idx]->getSiteWidth()) * rows[row_idx]->getSiteWidth();
        for(int x = alignedStartX; x <= subrow->getEndX(); x += rows[row_idx]->getSiteWidth()){
            Coor currCoor = Coor(x, rows[row_idx]->getStartCoor().y);
            if(ff->getDisplacement(currCoor) > minDisplacement){
                Coor subrowEndCoor = Coor(subrow->getEndX(), rows[row_idx]->getStartCoor().y);
                // If current subrow can't find better solution
                if(ff->getDisplacement(subrowEndCoor) > minDisplacement) break;
                continue;
            } 
            bool placeable = ContinousAndEmpty(x, rows[row_idx]->getStartCoor().y, ff->getW(), ff->getH(), row_idx);
            double displacement = ff->getDisplacement(currCoor);
            if(placeable && displacement < minDisplacement){
                minDisplacement = displacement;
                ff->setLGCoor(currCoor);
                ff->setIsPlace(true);
            }
        }
    }
    return minDisplacement;
}

bool Legalizer::ContinousAndEmpty(double startX, double startY, double w, double h, int row_idx){
    double endY = startY + h;
    double currentY = startY;

    for(size_t i = row_idx; i < rows.size(); i++){
        double rowStartY = rows[i]->getStartCoor().y;
        double rowEndY = rowStartY + rows[i]->getSiteHeight();

        // check if the current row is complete beyond the target endY
        if(rowStartY >= endY) break;
        double placeH;
        if(h >= rows[i]->getSiteHeight()){
            placeH = rows[i]->getSiteHeight();
        }else{
            placeH = h;
        }
        // check if this row can place the cell from startX to startX + w
        if(rows[i]->canPlace(startX, startX + w, placeH)){
            // Ensure the row covers the currentY to at least part of the target range
            if(rowStartY <= currentY && rowEndY > currentY){
                currentY = rowEndY;
            }

            // Check if we reached the target height
            if(currentY >= endY) return true;
        }
    }
    return false;
}

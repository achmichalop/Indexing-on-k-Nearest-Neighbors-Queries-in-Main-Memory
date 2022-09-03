#include "process.h"

int myQuotient(double numer, double denom) {
        return int(numer/denom + EPS);
    };


double myRemainder(double numer, double denom, int q) {
    double rem = double(numer - q*denom);
    return ((abs(rem) < EPS) ? 0: rem);
};


int getCellId(int x, int y, int numCellsPerDimension) {
    return (y * numCellsPerDimension + x);
};


Coord minDist(Coordinates &point, Coord xStart, Coord yStart, Coord xEnd, Coord yEnd){
    Coord sum = 0.0;
    Coord diff;


    diff = 0.0;
    if ( point.x < xStart) {
        diff = xStart - point.x;
    }
    else if ( point.x > xEnd){
        diff = point.x- xEnd;
    }
    sum = diff*diff;
    diff = 0.0;

    if ( point.y < yStart ) {
        diff = yStart- point.y;
    }
    else if ( point.y > yEnd ){
        diff = point.y - yEnd;
    }
    sum += diff*diff;
            
    return sum;
}

void init_pq_Cell(priority_queue<CellNode, vector<CellNode>, CellComparator>*pqueue, Coordinates point, int runNumPartitionsPerRelation)
{
    int xStartCell, yStartCell, cellNumber;
    Coord partitionExtent = 1.0/runNumPartitionsPerRelation;
    Coord dist = 0.0;
    Coord diff = 0.0;

    xStartCell = myQuotient(point.x + EPS, partitionExtent);
    yStartCell = myQuotient(point.y + EPS, partitionExtent);
    
    int Cq = getCellId(xStartCell, yStartCell, runNumPartitionsPerRelation);
    
    CellNode pqNode;
    pqNode.dist = -1.0;
    pqNode.id = Cq;
    pqNode.flag = 4;
    pqNode.level = 0;

    pqueue->push(pqNode);    
    
    if (yStartCell-1 >= 0) {
        cellNumber = Cq - runNumPartitionsPerRelation;
        diff = point.y - yStartCell * partitionExtent;
        dist = diff * diff;

        CellNode pqNode;
        pqNode.dist = dist;
        pqNode.id = cellNumber;
        pqNode.flag = 1;
        pqNode.level = 2;

        pqueue->push(pqNode);
    }

    if (xStartCell-1 >= 0) {
        cellNumber = Cq - 1;
        diff = point.x - xStartCell * partitionExtent;
        dist = diff * diff;

        CellNode pqNode;
        pqNode.dist = dist;
        pqNode.id = cellNumber;
        pqNode.flag = 3;
        pqNode.level = 2;

        pqueue->push(pqNode);
 
    }

    if (yStartCell+1 <= runNumPartitionsPerRelation - 1) {
        cellNumber = Cq + runNumPartitionsPerRelation;
        diff = partitionExtent * (yStartCell+1) - point.y;
        dist = diff * diff;
        
        CellNode pqNode;
        pqNode.dist = dist;
        pqNode.id = cellNumber;
        pqNode.flag = 7;
        pqNode.level = 2;

        pqueue->push(pqNode);

    }

    if (xStartCell+1 <= runNumPartitionsPerRelation - 1) {
        cellNumber = Cq + 1;
        diff = partitionExtent * (xStartCell+1) - point.x;
        dist = diff * diff;
        
        CellNode pqNode;
        pqNode.dist = dist;
        pqNode.id = cellNumber;
        pqNode.flag = 5;
        pqNode.level = 2;

        pqueue->push(pqNode);
    }    
}


void push_next_level_cells(priority_queue<CellNode, vector<CellNode>, CellComparator> *pqueue, CellNode node, Coordinates point, int runNumPartitionsPerRelation, Coord bound)
{
    Coord dist;
    Coord partitionExtent = 1.0/runNumPartitionsPerRelation; 
    int level = node.level/2;
    int xStartCell = node.id % runNumPartitionsPerRelation;
    int yStartCell = node.id / runNumPartitionsPerRelation;
    Coord crossDist = node.dist;
    
    CellNode pqNode;
    pqNode.dist = crossDist;
    pqNode.id = node.id;
    pqNode.flag = node.flag;
    pqNode.level = 0;
    
    pqueue->push(pqNode);
  
    if (node.flag == 1) 
    {   
        Coord diffLeft = point.x - (xStartCell) * partitionExtent;

        for (int i = 1; i <= level; i++) {
            if (xStartCell-i < 0)
                break;

            int cellNumber = node.id - i;
            
            Coord diff = diffLeft * diffLeft;
            dist = crossDist + diff;

            if (dist >= bound)
                break;
            
            CellNode pqNode;
            pqNode.dist = dist;
            pqNode.id = cellNumber;
            pqNode.flag = node.flag - 1;
            pqNode.level = 0;

            pqueue->push(pqNode);

            diffLeft += partitionExtent;
        }

        Coord diffRight = partitionExtent * (xStartCell+1) - point.x;

        for (int i = 1; i <= level-1; i++) {
            if (xStartCell+i > runNumPartitionsPerRelation - 1)
                break;
            
            int cellNumber = node.id + i;

            Coord diff = diffRight * diffRight;
            dist = crossDist + diff;

            if (dist >= bound)
                break;
            
            CellNode pqNode;
            pqNode.dist = dist;
            pqNode.id = cellNumber;
            pqNode.flag = node.flag + 1;
            pqNode.level = 0;

            pqueue->push(pqNode);

            diffRight += partitionExtent;
        }

        if (yStartCell-1 >= 0) {
            Coord diff = point.y - yStartCell * partitionExtent;
            dist = diff * diff;

            if (dist < bound) {
                int cellNumber = node.id - runNumPartitionsPerRelation;
                
                CellNode pqNode;
                pqNode.dist = dist;
                pqNode.id = cellNumber;
                pqNode.flag = node.flag;
                pqNode.level = node.level + 2;

                pqueue->push(pqNode);
            }
        } 
    }
    else if (node.flag == 7)
    {   
        Coord diffRight = partitionExtent * (xStartCell+1) - point.x;

        for (int i = 1; i <= level; i++) {
            if (xStartCell+i > runNumPartitionsPerRelation - 1)
                break;
            
            int cellNumber = node.id + i;

            Coord diff = diffRight * diffRight;
            dist = crossDist + diff;

            if (dist >= bound)
                break;

            CellNode pqNode;
            pqNode.dist = dist;
            pqNode.id = cellNumber;
            pqNode.flag = node.flag + 1;
            pqNode.level = 0;

            pqueue->push(pqNode);

            diffRight += partitionExtent;
        }

        Coord diffLeft = point.x - (xStartCell) * partitionExtent;

        for (int i = 1; i <= level-1; i++) {
            if (xStartCell-i < 0)
                break;

            int cellNumber = node.id - i;
            
            Coord diff = diffLeft * diffLeft;
            dist = crossDist + diff;

            if (dist >= bound)
                break;
            
            CellNode pqNode;
            pqNode.dist = dist;
            pqNode.id = cellNumber;
            pqNode.flag = node.flag - 1;
            pqNode.level = 0;

            pqueue->push(pqNode);

            diffLeft += partitionExtent;
        }

        if (yStartCell+1 <= runNumPartitionsPerRelation - 1) {
            Coord diff = partitionExtent * (yStartCell+1) - point.y;
            dist = diff * diff;

            if (dist < bound) {
                int cellNumber = node.id + runNumPartitionsPerRelation;
                
                CellNode pqNode;
                pqNode.dist = dist;
                pqNode.id = cellNumber;
                pqNode.flag = node.flag;
                pqNode.level = node.level + 2;

                pqueue->push(pqNode);

            }
        }
    }
    else if (node.flag == 3) {
        Coord diffTop = partitionExtent * (yStartCell+1) - point.y;

        for (int i = 1; i <= level; i++) {
            if (yStartCell+i > runNumPartitionsPerRelation - 1)
                break;
            
            int cellNumber = node.id + (i * runNumPartitionsPerRelation);

            Coord diff = diffTop * diffTop;
            dist = crossDist + diff;

            if (dist >= bound)
                break;


            CellNode pqNode;
            pqNode.dist = dist;
            pqNode.id = cellNumber;
            pqNode.flag = node.flag + 3;
            pqNode.level = 0;

            pqueue->push(pqNode);

            diffTop += partitionExtent;
        }

        Coord diffBottom = point.y - yStartCell * partitionExtent;

        for (int i = 1; i <= level-1; i++) {
            if (yStartCell-i < 0)
                break;
            
            int cellNumber = node.id - (i * runNumPartitionsPerRelation);

            Coord diff = diffBottom * diffBottom;
            dist = crossDist + diff;

            if (dist >= bound)
                break;

            CellNode pqNode;
            pqNode.dist = dist;
            pqNode.id = cellNumber;
            pqNode.flag = node.flag - 3;
            pqNode.level = 0;

            pqueue->push(pqNode);
            
            diffBottom += partitionExtent;
        }

        if (xStartCell-1 >= 0) {
            Coord diff = point.x - xStartCell * partitionExtent;
            dist = diff * diff;

            if (dist < bound) {
                int cellNumber = node.id - 1;
                
                CellNode pqNode;
                pqNode.dist = dist;
                pqNode.id = cellNumber;
                pqNode.flag = node.flag;
                pqNode.level = node.level + 2;

                pqueue->push(pqNode);
            }
        }
    }
    else // if (node.flag == 5) => right direction!!
    {
        Coord diffBottom = point.y - yStartCell * partitionExtent;

        for (int i = 1; i <= level; i++) {
            if (yStartCell-i < 0)
                break;
            
            int cellNumber = node.id - (i * runNumPartitionsPerRelation);

            Coord diff = diffBottom * diffBottom;
            dist = crossDist + diff;

            if (dist >= bound)
                break;


            CellNode pqNode;
            pqNode.dist = dist;
            pqNode.id = cellNumber;
            pqNode.flag = node.flag - 3;
            pqNode.level = 0;

            pqueue->push(pqNode);
            
            diffBottom += partitionExtent;
        }

        Coord diffTop = partitionExtent * (yStartCell+1) - point.y;

        for (int i = 1; i <= level-1; i++) {
            if (yStartCell+i > runNumPartitionsPerRelation - 1)
                break;
            
            int cellNumber = node.id + (i * runNumPartitionsPerRelation);

            Coord diff = diffTop * diffTop;
            dist = crossDist + diff;

            if (dist >= bound)
                break;
            

            CellNode pqNode;
            pqNode.dist = dist;
            pqNode.id = cellNumber;
            pqNode.flag = node.flag + 3;
            pqNode.level = 0;

            pqueue->push(pqNode);

            diffTop += partitionExtent;
        }

        if (xStartCell+1 <= runNumPartitionsPerRelation - 1) {
            Coord diff = partitionExtent * (xStartCell+1) - point.x;
            dist = diff * diff;

            if (dist < bound) {
                int cellNumber = node.id + 1;
                
                CellNode pqNode;
                pqNode.dist = dist;
                pqNode.id = cellNumber;
                pqNode.flag = node.flag;
                pqNode.level = node.level + 2;

                pqueue->push(pqNode);

            }
        }
    }
}

void pushObjects(PriorityQueue *pqueue, Coordinates point, CellNode node, Relation *pR, size_t *pR_sizes, Coord *bound)
{
    Coord dist;
    int cellNumber = node.id;

    if (node.flag == 0)
    {
        for (int j = 0; j < 16; j+=8) {
            for (int i = pR_sizes[cellNumber * 16 + j]; i < pR_sizes[cellNumber * 16 + j + 1]; i++) {
                
                dist = minDist(point, pR[cellNumber][i].xStart, pR[cellNumber][i].yStart, pR[cellNumber][i].xEnd, pR[cellNumber][i].yEnd);

                PriorityQueueNode ObjNode = PriorityQueueNode(dist, pR[cellNumber][i].id);
                pqueue->insert(ObjNode);
            }

            for (int i = pR_sizes[cellNumber * 16 + j + 2]; i < pR_sizes[cellNumber * 16 + j + 3]; i++) {
                dist = minDist(point, pR[cellNumber][i].xStart, pR[cellNumber][i].yStart, pR[cellNumber][i].xEnd, pR[cellNumber][i].yEnd);

                PriorityQueueNode ObjNode = PriorityQueueNode(dist, pR[cellNumber][i].id);
                pqueue->insert(ObjNode);

            }
        }
    }
    else if (node.flag == 1)
    {
        for (int j = 0; j < 16; j+=2) {
            for (int i = pR_sizes[cellNumber * 16 + j]; i < pR_sizes[cellNumber * 16 + j + 1]; i++) {
                dist = minDist(point, pR[cellNumber][i].xStart, pR[cellNumber][i].yStart, pR[cellNumber][i].xEnd, pR[cellNumber][i].yEnd);

                PriorityQueueNode ObjNode = PriorityQueueNode(dist, pR[cellNumber][i].id);
                pqueue->insert(ObjNode);
            }
        }
    }
    else if (node.flag == 2)
    {
        for (int j = 0; j < 8; j+=2) {
            for (int i = pR_sizes[cellNumber * 16 + j]; i < pR_sizes[cellNumber * 16 + j + 1]; i++) {
                dist = minDist(point, pR[cellNumber][i].xStart, pR[cellNumber][i].yStart, pR[cellNumber][i].xEnd, pR[cellNumber][i].yEnd);

                PriorityQueueNode ObjNode = PriorityQueueNode(dist, pR[cellNumber][i].id);
                pqueue->insert(ObjNode);
            }
        }
    }
    else if (node.flag == 3)
    {
        for (int j = 0; j < 16; j+=8) {
            for (int i = pR_sizes[cellNumber * 16 + j]; i < pR_sizes[cellNumber * 16 + j + 4]; i++) {
                dist = minDist(point, pR[cellNumber][i].xStart, pR[cellNumber][i].yStart, pR[cellNumber][i].xEnd, pR[cellNumber][i].yEnd);

                PriorityQueueNode ObjNode = PriorityQueueNode(dist, pR[cellNumber][i].id);
                pqueue->insert(ObjNode);
            }
        }
    }
    else if (node.flag == 4)
    {

        for (int i = 0; i < pR[cellNumber].size(); i++) {
            dist = minDist(point, pR[cellNumber][i].xStart, pR[cellNumber][i].yStart, pR[cellNumber][i].xEnd, pR[cellNumber][i].yEnd);

            PriorityQueueNode ObjNode = PriorityQueueNode(dist, pR[cellNumber][i].id);
            pqueue->insert(ObjNode);
        }
        

    }
    else if (node.flag == 5)
    {
        for (int i = pR_sizes[cellNumber * 16]; i < pR_sizes[cellNumber * 16 + 8]; i++) {
            dist = minDist(point, pR[cellNumber][i].xStart, pR[cellNumber][i].yStart, pR[cellNumber][i].xEnd, pR[cellNumber][i].yEnd);

            PriorityQueueNode ObjNode = PriorityQueueNode(dist, pR[cellNumber][i].id);
            pqueue->insert(ObjNode);
        }
    }
    else if (node.flag == 6)
    {
        for (int j = 0; j < 16; j+=8) {
            for (int i = pR_sizes[cellNumber * 16 + j]; i < pR_sizes[cellNumber * 16 + j + 2]; i++) {
                dist = minDist(point, pR[cellNumber][i].xStart, pR[cellNumber][i].yStart, pR[cellNumber][i].xEnd, pR[cellNumber][i].yEnd);
                
                PriorityQueueNode ObjNode = PriorityQueueNode(dist, pR[cellNumber][i].id);
                pqueue->insert(ObjNode);
            }
        }
    }
    else if (node.flag == 7)
    {
        for (int j = 0; j < 16; j+=4) {
            for (int i = pR_sizes[cellNumber * 16 + j]; i < pR_sizes[cellNumber * 16 + j + 2]; i++) {
                dist = minDist(point, pR[cellNumber][i].xStart, pR[cellNumber][i].yStart, pR[cellNumber][i].xEnd, pR[cellNumber][i].yEnd);

                PriorityQueueNode ObjNode = PriorityQueueNode(dist, pR[cellNumber][i].id);
                pqueue->insert(ObjNode);
            }
        }
    }
    else 
    {
        for (int j = 0; j < 8; j+=4) {
            for (int i = pR_sizes[cellNumber * 16 + j]; i < pR_sizes[cellNumber * 16 + j + 2]; i++) {
                dist = minDist(point, pR[cellNumber][i].xStart, pR[cellNumber][i].yStart, pR[cellNumber][i].xEnd, pR[cellNumber][i].yEnd);
                
                PriorityQueueNode ObjNode = PriorityQueueNode(dist, pR[cellNumber][i].id);
                pqueue->insert(ObjNode);
            }
        }
    }

    if (pqueue->size() == pqueue->capacity) {
        *bound = pqueue->at(0).dist;
    }
}

void nearest(int knn_size, Coordinates point, Relation *pR, size_t *pR_sizes, int runNumPartitionsPerRelation) 
{

    Coord bound = numeric_limits<double>::max();

    PriorityQueue pqueueObjs(knn_size);
    priority_queue <CellNode, vector<CellNode>, CellComparator> pqueueCells;

    CellNode cellPop;

    init_pq_Cell(&pqueueCells, point, runNumPartitionsPerRelation);

    cellPop = pqueueCells.top();
    pqueueCells.pop();                 
    pushObjects(&pqueueObjs, point, cellPop, pR, pR_sizes, &bound);
    
    while (!pqueueCells.empty())
    {   

        cellPop = pqueueCells.top();
        pqueueCells.pop();

        if (pqueueObjs.currCapacity == knn_size && cellPop.dist >= bound)
            break;

        if (cellPop.level > 0) {
            push_next_level_cells(&pqueueCells, cellPop, point, runNumPartitionsPerRelation, bound);
        }
        else {
            pushObjects(&pqueueObjs, point, cellPop, pR, pR_sizes, &bound);
        }
    }

}

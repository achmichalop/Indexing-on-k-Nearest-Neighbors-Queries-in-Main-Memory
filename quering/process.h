#ifndef _QUERY_H_
#define _QUERY_H_

#include "../def.h"
#include "../containers/relation.h"
#include "../PriorityQueue/priorityQueue.h"

class CellNode {
    public:
        Coord dist;
        int id, flag, level;

};

struct CellComparator
{
    bool operator()(const CellNode& lhs, const CellNode& rhs)
    {
        return lhs.dist > rhs.dist;
    }
};
int myQuotient(double numer, double denom);

double myRemainder(double numer, double denom, int q);

int getCellId(int x, int y, int numCellsPerDimension);


Coord minDist(Coordinates &point, Coord xStart, Coord yStart, Coord xEnd, Coord yEnd);

void init_pq_Cell(priority_queue<CellNode, vector<CellNode>, CellComparator>*pqueue, Coordinates point, int runNumPartitionsPerRelation);

void push_next_level_cells(priority_queue<CellNode, vector<CellNode>, CellComparator> *pqueue, CellNode node, Coordinates point, int runNumPartitionsPerRelation, Coord bound);

void pushObjects(PriorityQueue *pqueue, Coordinates point, CellNode node, Relation *pR, size_t *pR_sizes, Coord *bound);

void nearest(int knn_size, Coordinates point, Relation *pR, size_t *pR_sizes, int runNumPartitionsPerRelation);

#endif //_QUERY_H_

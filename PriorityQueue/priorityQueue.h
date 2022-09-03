#ifndef _PRIORITYQUEUE_H_
#define _PRIORITYQUEUE_H_

#include "../def.h"


class PriorityQueueNode{
    public :
        Coord dist;
        RecordId id;

        PriorityQueueNode();
        PriorityQueueNode(Coord dist, RecordId id);
        ~PriorityQueueNode();
        bool operator < (const PriorityQueueNode& rhs) const;


};


class PriorityQueue : public vector<PriorityQueueNode>{
    public:
        int currCapacity, capacity;

        PriorityQueue(int k);
        ~PriorityQueue();
        void compute();
        void enqueue(PriorityQueueNode pqNode);
        int parent (int position);
        void movedown();
        void insert(PriorityQueueNode pqNode);
        PriorityQueueNode dequeue();


};



#endif //_PRIORITYQUEUE_H_
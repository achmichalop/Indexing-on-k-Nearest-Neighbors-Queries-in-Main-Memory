#include "priorityQueue.h"
#include <random>



PriorityQueueNode :: PriorityQueueNode(){

}


PriorityQueueNode :: PriorityQueueNode(Coord dist, RecordId id){
    this->dist = dist;
    this->id = id;
} 


PriorityQueueNode :: ~PriorityQueueNode(){

}

int PriorityQueue :: parent (int position){
    if (position%2){
        return position/2;
    }
    else{
        return (position-1)/2;
    }
}


void PriorityQueue :: insert(PriorityQueueNode pqNode){
    
    if (currCapacity < capacity){     
        enqueue(pqNode);
        currCapacity++;
    }   
    else{
        if (pqNode.dist > this->at(0).dist){
        //if (num > this->at(0).dist){
            return;
        }
        else{
            this->at(0) = pqNode;

            movedown();
        }
    }
}


void PriorityQueue :: compute(){
    Coord num;
    Coord a[20];

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(1.0, 2.0);

    for ( int i = 0 ; i < 20 ; i++ ){
        num = dis(gen);
        a[i] = num;
        PriorityQueueNode pqNode = PriorityQueueNode(num,i);

        if (currCapacity < capacity){
            
            enqueue(pqNode);
            currCapacity++;
        }   
        else{
            if (num < this->at(0).dist){
            //if (num > this->at(0).dist){
                this->at(0) = pqNode;

                movedown();
            }
        }
    }
}


void PriorityQueue :: movedown(){

    int pos = 0;
    int swap;
    PriorityQueueNode tmp;

    while (pos*2 + 1 < this->size()){
        if (pos*2 + 2 < this->size()){
            //if (this->at(pos+1).dist < this->at(pos+2).dist ){
            if (this->at(pos*2+1).dist < this->at(pos*2+2).dist ){
                swap = pos*2+2;
            }
            else{
                
                swap = pos*2+1;
            }
        }
        else{
            swap = pos*2+1;
        }

        //if (this->at(pos).dist > this->at(swap).dist){
        if (this->at(pos).dist > this->at(swap).dist){
           break;
        }
        else{
            tmp = this->at(swap);
            this->at(swap) = this->at(pos);
            this->at(pos) = tmp;

            pos = swap;
        }
    }


}

void PriorityQueue :: enqueue(PriorityQueueNode pqNode){

    int lastPos = this->size();
    int p;
    PriorityQueueNode tmp;

    this->push_back(pqNode);
    
    while (lastPos > 0){
        p = parent(lastPos);

        //if (pqNode.dist < this->at(p).dist){
        if (pqNode.dist < this->at(p).dist){
            break;
        }
        else{
            tmp = this->at(p);
            this->at(p) = pqNode;
            this->at(lastPos) = tmp;

            lastPos = parent(lastPos);
        }

    }
}

PriorityQueueNode PriorityQueue :: dequeue() {
    PriorityQueueNode pqNode;
    pqNode = this->at(this->currCapacity-1);
    this->pop_back();
    this->currCapacity--;

    return pqNode;
} 


PriorityQueue :: PriorityQueue(int k){
    this->currCapacity = 0;
    this->capacity = k;
    this->reserve(this->capacity);
}


PriorityQueue :: ~PriorityQueue(){

}

bool PriorityQueueNode::operator < (const PriorityQueueNode& rhs) const
{
    return this->dist > rhs.dist;
}



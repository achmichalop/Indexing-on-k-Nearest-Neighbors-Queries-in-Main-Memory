#include <string>
#include <queue>
#include <limits>

#include "def.h"
#include "getopt.h"
#include "./containers/relation.h"
#include "./partitioning/partition.h"
#include "./grid/twoLevel.h"
#include "./PriorityQueue/priorityQueue.h"
#include "./quering/process.h"

void usage(){
    cerr << "NAME" << endl;
    cerr << "       ./twoLevel - range query using the 2-level algorithm" << endl << endl;
    cerr << "USAGE" << endl;
    cerr << "       ./twoLevel [OPTION]... [FILE1] [FILE2]" << endl << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "       Mandatory arguments" << endl << endl;
    cerr << "       -p" << endl;
    cerr << "              number of partions per dimension" << endl;
    cerr << "       -w" << endl;
    cerr << "              window query" << endl;
    cerr << "       Other arguments" << endl << endl;
    cerr << "       -i" << endl;
    cerr << "              number of iterations" << endl;
    cerr << "       -m" << endl;
    cerr << "              display this help message and exit" << endl << endl;
    cerr << "EXAMPLES" << endl;
    cerr << "       Window range query using the 2-level algorithm with 3000 partitions per dimension." << endl;
    cerr << "              ./twoLevel -p 3000 -w TIGER_ROADS_mbr.inp TIGER_c0.1%_n10000.qry" << endl;
    cerr << "\n" << endl;
    exit(1);
}

int main(int argc, char **argv)
{
    char c;
    Timer tim;
    Relation R, S, *pR;
    int knn_size = 0;
    int runNumPartitionsPerRelation = -1;
    double timeIndexingOrPartitioning = 0.0, timeKNNalgorithm = 0.0;
    size_t *pR_sizes;
    int runNumPartitions = -1;
    int queryMethod=-1;
    Coord epsilon = -1.0, x, y;
    Coord minX, maxX, minY, maxY, diffX, diffY, maxExtend;
    int NUM_ITERATIONS = 1;
    Coordinates point;

    while ((c = getopt(argc, argv, "p:f:wui:m")) != -1)
    {
        switch (c)
        {
            case 'p':
                runNumPartitionsPerRelation = atoi(optarg);
                break;
            case 'f':
                  knn_size = atoi(optarg);
                break;
            case 'w':
                queryMethod = NEAREST;
                break;
            case 'i':
                NUM_ITERATIONS = atoi(optarg);
                break;
            case 'm':
                usage();
                break;
            default:
                cerr << "Wrong arguments! ";
                usage();
                break;
        }
    }
    

    if(runNumPartitionsPerRelation == -1)
    {
        cerr << "Number of partitions is missing" << endl;
        usage();
    }

    if(queryMethod == -1)
    {
        cerr << "Query method is missing" << endl;
        usage();
    }

    if (queryMethod == DISK_QUERY){
        if (epsilon < 0){
            cout<<"Radius value is missing"<<endl;
            usage();
        }
    }
    
    // Load inputs
   #pragma omp parallel sections
    {
         #pragma omp section
         {
            R.load(argv[optind]);
         }
         #pragma omp section
         {
            S.load(argv[optind+1]);
         }
    }

    runNumPartitions = runNumPartitionsPerRelation * runNumPartitionsPerRelation;

    minX = min(R.minX, S.minX);
    maxX = max(R.maxX, S.maxX);
    minY = min(R.minY, S.minY);
    maxY = max(R.maxY, S.maxY);
    diffX = maxX - minX;
    diffY = maxY - minY;
    maxExtend = (diffX<diffY)?diffY:diffX;


    R.normalize(minX, maxX, minY, maxY, maxExtend);
    S.normalize(minX, maxX, minY, maxY, maxExtend);

    pR = new Relation[runNumPartitions];

    pR_sizes = new size_t[16*runNumPartitions];
    memset(pR_sizes, 0, 16*runNumPartitions*sizeof(size_t));
    
    tim.start();
    partition::twoLevel::single::PartitionTwoDimensional(R, pR, pR_sizes, runNumPartitionsPerRelation);
    timeIndexingOrPartitioning = tim.stop();
    

    switch (queryMethod)
    {
        case NEAREST:

            double resultTime = 0.0;

            for (size_t i = 0; i < S.size(); i++) {
                auto s = S[i];
                x = (s.xEnd+s.xStart)/2;
                y = (s.yEnd+s.yStart)/2;
                point = Coordinates(x,y);
                vector<Coord> times;
                
                for (size_t j = 0; j < NUM_ITERATIONS; j++) {

                    tim.start();

                    nearest(knn_size, point, pR, pR_sizes, runNumPartitionsPerRelation);

                    timeKNNalgorithm = tim.stop();
                    times.push_back(timeKNNalgorithm);
                    timeKNNalgorithm = 0.0;
                }

                sort(times.begin(), times.end());

                double sum = 0.0;

                for (size_t j = 1; j < NUM_ITERATIONS-1; j++) {
                    sum += times[j];
                }
                timeKNNalgorithm = sum / (NUM_ITERATIONS-2);
                
                // // // // cout << "Method\t\t" << "Indexing Time\t" << "k-NN Algorithm Time" << endl;
                // // // // cout << "Polygons\t" << timeIndexingOrPartitioning << "\t" << timeKNNalgorithm << endl<<endl;
                resultTime+= timeKNNalgorithm;
                timeKNNalgorithm = 0.0;
                
            }

            resultTime = resultTime / S.size();
            cout << "Method\t\t" << "Indexing Time\t" << "k-NN Average Time" << endl;
            cout << "Polygons\t\t" << timeIndexingOrPartitioning << " \t" << resultTime << endl<<endl;

            break;
    }

    delete[] pR_sizes;
    delete[] pR;
} 

    

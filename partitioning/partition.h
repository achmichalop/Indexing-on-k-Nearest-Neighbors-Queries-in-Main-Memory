
namespace partition
{
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


    int findReferenceCell(double x, double y, double cellExtent, int numCellsPerDimension) {
        int xInt,yInt;

        xInt = (x + EPS)/cellExtent;
        yInt = (y + EPS)/cellExtent;

        return (yInt * numCellsPerDimension + xInt);
    };

    int FindRelevantTiles(Record &s,vector<int> &cornerCells,vector<int> &insideCells,int runNumPartitionsPerRelation)
    {
        double xStartCell, yStartCell, xEndCell, yEndCell;

        int queryCase;
        int runNumPartitions = runNumPartitionsPerRelation*runNumPartitionsPerRelation;
        Coord partitionExtent = 1.0/runNumPartitionsPerRelation;
        xStartCell = myQuotient(s.xStart + EPS, partitionExtent);
        yStartCell = myQuotient(s.yStart + EPS, partitionExtent);
        auto xEnd = myRemainder(s.xEnd, partitionExtent, int(myQuotient(s.xEnd + EPS, partitionExtent)));
        auto yEnd = myRemainder(s.yEnd, partitionExtent, int(myQuotient(s.yEnd + EPS, partitionExtent)));
        

        if (s.xEnd + EPS >= 1) {
            xEndCell = runNumPartitionsPerRelation - 1;
        }
        else {
            xEndCell = myQuotient(s.xEnd + EPS, partitionExtent);
        }

        if (s.yEnd + EPS >= 1) {
            yEndCell = runNumPartitionsPerRelation - 1;
        }
        else {
            yEndCell = myQuotient(s.yEnd + EPS, partitionExtent);
        }

        int bottomLeft = getCellId(xStartCell, yStartCell, runNumPartitionsPerRelation);
        int bottomRight = getCellId(xEndCell, yStartCell, runNumPartitionsPerRelation);
        int topLeft = getCellId(xStartCell, yEndCell, runNumPartitionsPerRelation);
        int topRight = getCellId(xEndCell, yEndCell, runNumPartitionsPerRelation);
        
        
        if (bottomLeft == topLeft && bottomLeft == bottomRight ){
            cornerCells.push_back(bottomLeft);
            cornerCells.push_back(-1);
            cornerCells.push_back(-1);
            cornerCells.push_back(-1);
            queryCase = 0 ;
        }
        else if (bottomLeft == topLeft){
            cornerCells.push_back(bottomLeft);
            cornerCells.push_back(bottomRight);
            cornerCells.push_back(-1);
            cornerCells.push_back(-1);
            queryCase = 1;
        }
        else if ( bottomLeft == bottomRight){
            cornerCells.push_back(bottomLeft);
            cornerCells.push_back(-1);
            cornerCells.push_back(topLeft);
            cornerCells.push_back(-1);
            queryCase = 2;
        }
        else{
        
            cornerCells.push_back(bottomLeft);
            cornerCells.push_back(bottomRight);
            cornerCells.push_back(topLeft);
            cornerCells.push_back(topRight);
            queryCase = 3;
        }

        int rightLimit = bottomRight;
        int start = bottomLeft;

        int ii = start; 
        while( ii<=topRight)
        {
            if(ii <= rightLimit)
            {
                if(ii != start && ii != rightLimit && start != bottomLeft && rightLimit != topRight)  
                {
                    insideCells.push_back(ii);
                }
                ii++;
            }
            else{
                rightLimit +=runNumPartitionsPerRelation;
                start += runNumPartitionsPerRelation;
                ii=start;
            }
        }
        
        return queryCase;
    }

    namespace twoLevel
    {
        namespace single
        {

            void PartitionUniform(const Relation& R, Relation *pR, size_t *pR_sizes, int runNumPartitionsPerRelation)
            {
                int runNumPartitions = runNumPartitionsPerRelation*runNumPartitionsPerRelation;
                Coord partitionExtent = 1.0/runNumPartitionsPerRelation;
                double xStartCell, yStartCell, xEndCell, yEndCell;
                int firstCell, lastCell;
                Timer tim;
                double timepR = 0, timeDecomp = 0;

                for (size_t i = 0; i < R.numRecords; i++){
                    auto &r = R[i];
                    
                    // Determine cell for (rec.xStart, rec.yStart)
                    xStartCell = myQuotient(r.xStart + EPS, partitionExtent);
                    yStartCell = myQuotient(r.yStart + EPS, partitionExtent);
                    firstCell = getCellId(xStartCell, yStartCell, runNumPartitionsPerRelation);

                    // Determine cell for (rec.xEnd, rec.yEnd)
                    auto xEnd = myRemainder(r.xEnd, partitionExtent, int(myQuotient(r.xEnd + EPS, partitionExtent)));
                    auto yEnd = myRemainder(r.yEnd, partitionExtent, int(myQuotient(r.yEnd + EPS, partitionExtent)));

                    if (r.xEnd + EPS >= 1) {
                        xEndCell = runNumPartitionsPerRelation - 1;
                    }
                    else {
                        xEndCell = myQuotient(r.xEnd + EPS, partitionExtent);
                    }

                    if (r.yEnd + EPS >= 1) {
                        yEndCell = runNumPartitionsPerRelation - 1;
                    }
                    else {
                        yEndCell = myQuotient(r.yEnd + EPS, partitionExtent);
                    }

                    lastCell = getCellId(xEndCell, yEndCell, runNumPartitionsPerRelation);
                    
                    // cout << "For object " << (i+1) << ":\n";
                    // cout << "\txStart = " << r.xStart << ", yStart = " << r.yStart << endl;
                    // cout << "\txEnd = " << r.xEnd << ", yEnd = " << r.yEnd <<  "\n\n";
      
                    // cout << "\txStartC = " << xStartCell << ", yStartC = " << yStartCell << ", Fcell = " << firstCell << endl;
                    // cout << "\txEndC = " << xEndCell << ", yEndC = " << yEndCell << ", Lcell = " << lastCell << "\n\n";
                    
                    

                    // Put record in cells.
                    if (firstCell == lastCell) {
                        // cout << "\nCell " << firstCell << ": Inside-Inside" << endl;
                        pR_sizes[firstCell * 16]++;
                    }
                    else if (xStartCell == xEndCell) {
                        for (int j = yStartCell; j <= yEndCell; j++) {
                            if (j == yStartCell) {
                                // cout << "\nCell " << firstCell << ": Inside-Starts in" << endl;
                                pR_sizes[firstCell * 16 + 1]++;
                            }
                            else if (j == yEndCell) {
                                // cout << "\nCell " << lastCell << ": Inside-Ends in" << endl;
                                pR_sizes[lastCell * 16 + 2]++;
                            }
                            else {
                                int cellNumber = getCellId(xStartCell, j, runNumPartitionsPerRelation); 
                                // cout << "\nCell " << cellNumber << ": Inside-Covers" << endl;
                                pR_sizes[cellNumber * 16 + 3]++;
                            }
                        }
                    }
                    else {
                        int caseFlag = 0;
                        
                        if (yStartCell == yEndCell)
                        {

                            for (int i = xStartCell; i <= xEndCell; i++) {
                                int cellNumber = getCellId(i, yStartCell, runNumPartitionsPerRelation);

                                if (i == xStartCell) {
                                    // cout << "\nCell " << cellNumber << ": Starts in-Inside" << endl;
                                    pR_sizes[cellNumber * 16 + 4]++;
                                }
                                else if (i == xEndCell) {
                                    // cout << "\nCell " << cellNumber << ": Ends in-Inside" << endl;
                                    pR_sizes[cellNumber * 16 + 8]++;
                                }
                                else {
                                    // cout << "\nCell " << cellNumber << ": Covers-Inside" << endl;
                                    pR_sizes[cellNumber * 16 + 12]++;
                                }
                            }
                        }
                        else {
                            for (int j = yStartCell; j <= yEndCell; j++) {
                                if (j == yStartCell) {
                                    caseFlag = 1;
                                }
                                else if (j == yEndCell) {
                                    caseFlag = 2;
                                }
                                else {
                                    caseFlag = 3;
                                }
                                for (int i = xStartCell; i <= xEndCell; i++) {
                                    int cellNumber = getCellId(i, j, runNumPartitionsPerRelation);

                                    if (i == xStartCell) {
                                        // if (caseFlag == 1) {
                                        //     // cout << "\nCell " << cellNumber << ": Starts in-Starts in" << endl;
                                        // }
                                        // else if (caseFlag == 2) {
                                        //     // cout << "\nCell " << cellNumber << ": Starts in-Ends in" << endl;
                                        // }
                                        // else {
                                        //     // cout << "\nCell " << cellNumber << ": Starts in-Covers" << endl;
                                        // }
                                        pR_sizes[cellNumber * 16 + 4 + caseFlag]++;
                                    }
                                    else if (i == xEndCell) {
                                        // if (caseFlag == 1) {
                                        //     // cout << "\nCell " << cellNumber << ": Ends in-Starts in" << endl;
                                        // }
                                        // else if (caseFlag == 2) {
                                        //     // cout << "\nCell " << cellNumber << ": Ends in-Ends in" << endl;
                                        // }
                                        // else {
                                        //     // cout << "\nCell " << cellNumber << ": Ends in-Covers" << endl;
                                        // }
                                        pR_sizes[cellNumber * 16 + 8 + caseFlag]++;
                                    }
                                    else {
                                        // if (caseFlag == 1) {
                                        //     // cout << "\nCell " << cellNumber << ": Covers-Starts in" << endl;
                                        // }
                                        // else if (caseFlag == 2) {
                                        //     // cout << "\nCell " << cellNumber << ": Covers-Ends in" << endl;
                                        // }
                                        // else {
                                        //     // cout << "\nCell " << cellNumber << ": Covers-Covers" << endl;
                                        // }
                                        pR_sizes[cellNumber * 16 + 12 + caseFlag]++;
                                    }
                                }
                            }
                        }
                    }
                }

                for (int i = 0; i < runNumPartitions; i++) {
                    int counter = 0;
                    // cout << "For Cell " << i << ":\n";
                    for (int j = 0; j < 16; j += 4) {
                        counter += pR_sizes[i*16+j] + pR_sizes[i*16+j+1] + pR_sizes[i*16+j+2] + pR_sizes[i*16+j+3];
                        // if (j < 4) {
                        //     // cout << "\tInside - Inside: " << pR_sizes[i*16+j] << endl; 
                        //     // cout << "\tInside - Starts In: " << pR_sizes[i*16+j+1] << endl;
                        //     // cout << "\tInside - Ends In: " << pR_sizes[i*16+j+2] << endl;
                        //     // cout << "\tInside - Covers: " << pR_sizes[i*16+j+3] << endl;
                        // }
                        // else if (j < 8) {
                        //     // cout << "\tStarts In - Inside: " << pR_sizes[i*16+j] << endl;
                        //     // cout << "\tStarts In - Starts In: " << pR_sizes[i*16+j+1] << endl;
                        //     // cout << "\tStarts In - Ends In: " << pR_sizes[i*16+j+2] << endl;
                        //     // cout << "\tStarts In - Covers: " << pR_sizes[i*16+j+3] << endl;
                        // }
                        // else if (j < 12) {
                        //     // cout << "\tEnds In - Inside: " << pR_sizes[i*16+j] << endl;
                        //     // cout << "\tEnds In - Stars In: " << pR_sizes[i*16+j+1] << endl;
                        //     // cout << "\tEnds In - Ends In: " << pR_sizes[i*16+j+2] << endl;
                        //     // cout << "\tEnds In - Covers: " << pR_sizes[i*16+j+3] << endl;
                        // }
                        // else {
                        //     // cout << "\tCovers - Inside: " << pR_sizes[i*16+j] << endl;
                        //     // cout << "\tCovers - Starts In: " << pR_sizes[i*16+j+1] << endl;
                        //     // cout << "\tCovers - Ends In: " << pR_sizes[i*16+j+2] << endl;
                        //     // cout << "\tCovers - Covers: " << pR_sizes[i*16+j+3] << endl;
                        // }

                        if (j != 0)
                            pR_sizes[i*16+j] += pR_sizes[i*16+j-1];
                        
                        pR_sizes[i*16+j+1] += pR_sizes[i*16+j];
                        pR_sizes[i*16+j+2] += pR_sizes[i*16+j+1];
                        pR_sizes[i*16+j+3] += pR_sizes[i*16+j+2];
                    }
                    pR[i].resize(counter);
                    pR[i].numRecords = counter;
                    // cout << "NumRecords = " << pR[i].numRecords << "\n\n";
                }

                for (size_t i = 0; i < R.numRecords; i++){
                    auto &r = R[i];
                    
                    // Determine cell for (rec.xStart, rec.yStart)
                    xStartCell = myQuotient(r.xStart + EPS, partitionExtent);
                    yStartCell = myQuotient(r.yStart + EPS, partitionExtent);
                    firstCell = getCellId(xStartCell, yStartCell, runNumPartitionsPerRelation);

                    // Determine cell for (rec.xEnd, rec.yEnd)
                    auto xEnd = myRemainder(r.xEnd, partitionExtent, int(myQuotient(r.xEnd + EPS, partitionExtent)));
                    auto yEnd = myRemainder(r.yEnd, partitionExtent, int(myQuotient(r.yEnd + EPS, partitionExtent)));

                    if (r.xEnd + EPS >= 1) {
                        xEndCell = runNumPartitionsPerRelation - 1;
                    }
                    else {
                        xEndCell = myQuotient(r.xEnd + EPS, partitionExtent);
                    }

                    if (r.yEnd + EPS >= 1) {
                        yEndCell = runNumPartitionsPerRelation - 1;
                    }
                    else {
                        yEndCell = myQuotient(r.yEnd + EPS, partitionExtent);
                    }

                    lastCell = getCellId(xEndCell, yEndCell, runNumPartitionsPerRelation);
                    
                    // cout << "For object " << (i+1) << ":\n";
                    // cout << "\txStart = " << r.xStart << ", yStart = " << r.yStart << endl;
                    // cout << "\txEnd = " << r.xEnd << ", yEnd = " << r.yEnd <<  "\n\n";
        
                    // cout << "\txStartC = " << xStartCell << ", yStartC = " << yStartCell << ", Fcell = " << firstCell << endl;
                    // cout << "\txEndC = " << xEndCell << ", yEndC = " << yEndCell << ", Lcell = " << lastCell << "\n\n";
                    
                    

                    // Put record in cells.
                    if (firstCell == lastCell) {
                        pR[firstCell][--pR_sizes[firstCell*16]] = r;
                    }
                    else if (xStartCell == xEndCell) {
                        for (int j = yStartCell; j <= yEndCell; j++) {
                            if (j == yStartCell) {
                                pR[firstCell][--pR_sizes[firstCell*16+1]] = r;
                            }
                            else if (j == yEndCell) {
                                pR[lastCell][--pR_sizes[lastCell*16+2]] = r;
                            }
                            else {
                                int cellNumber = getCellId(xStartCell, j, runNumPartitionsPerRelation); 
                                pR[cellNumber][--pR_sizes[cellNumber*16+3]] = r;
                            }
                        }
                    }
                    else {
                        int caseFlag = 0;
                        
                        if (yStartCell == yEndCell)
                        {
                            for (int i = xStartCell; i <= xEndCell; i++) {
                                int cellNumber = getCellId(i, yStartCell, runNumPartitionsPerRelation);

                                if (i == xStartCell) {
                                    pR[cellNumber][--pR_sizes[cellNumber * 16 + 4]] = r;
                                }
                                else if (i == xEndCell) {
                                    pR[cellNumber][--pR_sizes[cellNumber * 16 + 8]] = r;
                                }
                                else {
                                    pR[cellNumber][--pR_sizes[cellNumber * 16 + 12]] = r;
                                }
                            }
                        }
                        else {
                            for (int j = yStartCell; j <= yEndCell; j++) {
                                if (j == yStartCell) {
                                    caseFlag = 1;
                                }
                                else if (j == yEndCell) {
                                    caseFlag = 2;
                                }
                                else {
                                    caseFlag = 3;
                                }

                                for (int i = xStartCell; i <= xEndCell; i++) {
                                    int cellNumber = getCellId(i, j, runNumPartitionsPerRelation);

                                    if (i == xStartCell) {
                                        pR[cellNumber][--pR_sizes[cellNumber * 16 + 4 + caseFlag]] = r;
                                    }
                                    else if (i == xEndCell) {
                                        pR[cellNumber][--pR_sizes[cellNumber * 16 + 8 + caseFlag]] = r;
                                    }
                                    else {
                                        pR[cellNumber][--pR_sizes[cellNumber * 16 + 12 + caseFlag]] = r;
                                    }
                                }
                            }
                        }
                    }
                }
            };


            void PartitionTwoDimensional(Relation& R, Relation *pR, size_t *pR_sizes, int runNumPartitionsPerRelation)
            {
                PartitionUniform(R, pR, pR_sizes, runNumPartitionsPerRelation);            
            };
            
        } 
        ;
    }
}


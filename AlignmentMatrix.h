//
// Created by tom on 30/08/18.
//

#ifndef ALIGNMENT_HPC_ALIGNMENTMATRIX_H
#define ALIGNMENT_HPC_ALIGNMENTMATRIX_H

#include <iostream>
#include <map>
#include <vector>

class AlignmentMatrix {

public:
    //variables
    int *scoring_matrix;
    int x, y;
    int mat_size;
    char *x_order, *y_order;
    //constructor/destructor
    AlignmentMatrix(char *x_in, char *y_in, int size, std::vector<int> &scores);
    ~AlignmentMatrix();
};





#endif //ALIGNMENT_HPC_ALIGNMENTMATRIX_H

//
// Created by tom on 30/08/18.
//

#include <iostream>

#include "AlignmentMatrix.h"

AlignmentMatrix::AlignmentMatrix(char *x_in, char *y_in, int size, std::vector<int> &scores){
    x = size;
    y = size;

    x_order = x_in;
    y_order = y_in;

    std::cout << "Making alignment matrix object of " << size << " x " << size << std::endl;
    std::cout << "Score array consists of " << scores.size() << " elements" << std::endl;

    mat_size = x * y;

    scoring_matrix = new int[mat_size];

    for (int i = 0; i < mat_size; i++) {
      scoring_matrix[i] = scores[i];
    }
}

AlignmentMatrix::~AlignmentMatrix(){
  std::cout << "Cleaning up alignment matrix" << std::endl;
  delete [] scoring_matrix;
}

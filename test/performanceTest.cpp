/* 
 *  Copyright (c) 2010 Daisuke Okanohara
 * 
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 * 
 *   1. Redistributions of source code must retain the above Copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above Copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   3. Neither the name of the authors nor the names of its contributors
 *      may be used to endorse or promote products derived from this
 *      software without specific prior written permission.
 */

#include <iostream>
#include <set>
#include <sys/time.h>
#include "../src/redsvd.hpp"

using namespace std;
using namespace Eigen;

const int MAX_SIZE        = 100000;
const int SQUARE_MAX_SIZE = 10000;
const int SVD_MAX_SIZE    = 1000;

double getSec(){
  timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + (double)tv.tv_usec*1e-6;
}

void runDenseSVD(int row, int col, int r){
  MatrixXf A = MatrixXf::Random(row, col);

  cout << row << "\t" << col << "\t";

  {
    double start = getSec();
    REDSVD::RedSVD REDSVDOfA(A, r);
    cout << getSec() - start << "\t";
  }

  if (row <= SVD_MAX_SIZE && col <= SVD_MAX_SIZE){
    double start = getSec();
    JacobiSVD<MatrixXf> SVDOfA(A);
    cout << getSec() - start << "\t";
  }

  cout << endl;
}

void runSparseSVD(int row, int col, float ratio){
  Eigen::SparseMatrix<float, Eigen::RowMajor> A(row, col);

  cout << row << "\t" << col << "\t" << ratio << "\t";
  int colValNum = (int)(col * ratio);
  for (int i = 0; i < row; ++i){
    A.startVec(i);
    set<int> v;
    for (;(int)v.size() < colValNum;){
      int k = rand() % col;
      if (v.find(k) != v.end()) continue;
      v.insert(k);
    }
    for (set<int>::iterator it = v.begin();
	 it != v.end(); ++it){
      A.insertBack(i, *it) = (float)rand();
    }
  }
  A.finalize();

  double start = getSec();
  for (int r = 10; r <= 80; r *= 2){
    REDSVD::RedSVD redSVD(A, r);
    cout << getSec() - start << "\t";  
  }
  cout << endl;

}


void denseMatrixTest(){
  cout << "Dense Matrix" << endl;
  cout << "row\tcol\tredSVD\tSVD" << endl;

  cout << "row fixed:" << endl;
  vector<int> ranks;
  ranks.push_back(10);
  ranks.push_back(100);
  ranks.push_back(1000);
  
  for (size_t ir = 0; ir < ranks.size(); ++ir){
    int r = ranks[ir];
    cout << "rank:\t" << r << endl;
    for (int i = 100; i <= MAX_SIZE; i *= 10){
      runDenseSVD(i, 500, r);
    }
  }

  cout << "col fixed:" << endl;
  for (size_t ir = 0; ir < ranks.size(); ++ir){
    int r = ranks[ir];
    cout << "rank:\t" << r << endl;
    for (int i = 100; i <= MAX_SIZE; i *= 10){
      runDenseSVD(500, i, r);
    }
  }

  cout << "square:" << endl;
  for (size_t ir = 0; ir < ranks.size(); ++ir){
    int r = ranks[ir];
    cout << "rank:\t" << r << endl;
    for (int i = 100; i <= MAX_SIZE; i *= 10){
      runDenseSVD(i, i, r);
    }
  }


  cout << "full rank:" << endl;
  for (int i = 100; i <= SQUARE_MAX_SIZE; i *= 10){
    runDenseSVD(i, i, i);
  }
}


void sparseMatrixTest(){
  cout << "Sparse Matrix" << endl;
  cout << "row\tcol\tratio\t10\t20\t40\t80" << endl;
  for (float ratio = 0.001f; ratio <= 0.1f; ratio *= 10){
    for (int size = 100; size <= 100000; size *= 10){
      runSparseSVD(size, size, ratio);
    }
  }
}

int main(int argc, char* argv[]){
  sparseMatrixTest();
  denseMatrixTest();
  return 0;
}


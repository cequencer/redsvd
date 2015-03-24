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
#include "../src/redsvd.hpp"

using namespace std;
using namespace Eigen;
using namespace REDSVD;

const float EPS = 0.001f;

void denseTest(int row, int col, int actualRank, int estimateRank){
  cout << "dense " << row << "\t" << col << "\t" << actualRank << "\t" << estimateRank << endl;
  MatrixXf U = MatrixXf::Random(row, actualRank);
  MatrixXf V=  MatrixXf::Random(col, actualRank);
  Util::processGramSchmidt(U);
  Util::processGramSchmidt(V);
  VectorXf S(actualRank);
  for (int i = 0; i < actualRank; ++i){
    S(i) = pow(0.9f, i);
  }
  MatrixXf A = U * S.asDiagonal() * V.transpose();
  RedSVD svdOfA(A, estimateRank);
  

  for (int i = 0; i < estimateRank; ++i){
    cout << i << "\t" << log(S(i)) << "\t" <<  log(svdOfA.singularValues()(i)) << "\t" 
	 << fabs(U.col(i).dot(svdOfA.matrixU().col(i))) << "\t"
	 << fabs(V.col(i).dot(svdOfA.matrixV().col(i))) << endl;
  }
}

void sparseTest(int row, int col){
  cout << "sparse " << row << "\t" << col << endl;
  SMatrixXf A(row, col);
  int colVal = (int)(col * 0.10);
  for (int i = 0; i < row; ++i){
    A.startVec(i);
    set<int> v;
    while ((int)v.size() < colVal){
      int k = rand() % col;
      if (v.find(k) != v.end()) continue;
      v.insert(k);
    }
    for (set<int>::iterator it = v.begin();
	 it != v.end(); ++it){
      A.insertBack(i, *it) = (float)rand() / RAND_MAX;
    }
  }
  A.finalize();

  RedSVD svdOfA(A);
  MatrixXf newA = svdOfA.matrixU() * svdOfA.singularValues().asDiagonal() * svdOfA.matrixV().transpose();

  MatrixXf oldA = MatrixXf::Zero(row, col);
  for (int k=0; k<A.outerSize(); ++k) {
    for (SMatrixXf::InnerIterator it(A, k); it; ++it) {
      oldA(it.row(), it.col()) = it.value();
    }
  }

  bool error = false;
  for (int i = 0; i < row; ++i){
    for (int j = 0; j < col; ++j){
      if (fabs(oldA(i, j) - newA(i, j)) <= EPS) continue;
      error = true;
    }
  }
}

int main(int argc, char* argv[]){
  denseTest(100,  100,  10, 10);
  denseTest(100,  100,  20, 10);
  denseTest(100,  100, 100, 10);

  denseTest(1000,  1000,   10, 10);
  denseTest(1000,  1000,   20, 10);
  denseTest(1000,  1000,  100, 10);
  denseTest(1000,  1000, 1000, 10);

  return 0;
}

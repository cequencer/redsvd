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

#include <gtest/gtest.h>
#include "../src/redsvd.hpp"
#include "../src/redsvdFile.hpp"

using namespace std;
using namespace Eigen;

const float EPS = 0.001f;

TEST(redsvd, trivial){
  MatrixXf A;
  REDSVD::RedSVD redsvd;
  redsvd.run(A, 0);
  MatrixXf U = redsvd.matrixU();
  VectorXf S = redsvd.singularValues();
  MatrixXf V = redsvd.matrixV();
  ASSERT_EQ(0, U.rows());
  ASSERT_EQ(0, U.cols());
  ASSERT_EQ(0, V.rows());
  ASSERT_EQ(0, V.cols());
  ASSERT_EQ(0, S.rows());
  ASSERT_EQ(1, S.cols());
}


TEST(redsvd, artificial){
  int rowN = 3;
  int colN = 5;
  MatrixXf A(rowN, colN);
  A << 1,  2,  3,  4,  5,
       6,  7,  8,  9, 10, 
      11, 12, 13, 14, 15;
  
  int r = 2;
  REDSVD::RedSVD redsvd;
  redsvd.run(A, r);
  MatrixXf U = redsvd.matrixU();
  VectorXf S = redsvd.singularValues();
  MatrixXf V = redsvd.matrixV();

  ASSERT_EQ(rowN, U.rows());
  ASSERT_EQ(r,    U.cols());
  ASSERT_EQ(colN, V.rows());
  ASSERT_EQ(r,    V.cols());
  ASSERT_EQ(r,    S.rows());
  ASSERT_EQ(1,    S.cols());

  for (int i = 0; i < U.cols(); ++i){
    ASSERT_NEAR(1.f, U.col(i).norm(), EPS);
    for (int j = 0; j < i; ++j){
      ASSERT_NEAR(0.f, U.col(i).dot(U.col(j)), EPS);
    }
  }

  for (int i = 0; i < V.cols(); ++i){
    ASSERT_NEAR(1.f, V.col(i).norm(), EPS);
    for (int j = 0; j < i; ++j){
      ASSERT_NEAR(0.f, V.col(i).dot(V.col(j)), EPS);
    }
  }

  ASSERT_NEAR(0.f, (A - U * S.asDiagonal() * V.transpose()).norm(), EPS);
}

TEST(redsvd, random){
  int rowN = 20;
  int colN = 70;
  MatrixXf A = MatrixXf::Random(rowN, colN);

  int r = 20;
  REDSVD::RedSVD redsvd;
  redsvd.run(A, r);
  MatrixXf U = redsvd.matrixU();
  VectorXf S = redsvd.singularValues();
  MatrixXf V = redsvd.matrixV();

  ASSERT_EQ(rowN, U.rows());
  ASSERT_EQ(r,    U.cols());
  ASSERT_EQ(colN, V.rows());
  ASSERT_EQ(r,    V.cols());
  ASSERT_EQ(r,    S.rows());
  ASSERT_EQ(1,    S.cols());

  for (int i = 0; i < U.cols(); ++i){
    ASSERT_NEAR(1.f, U.col(i).norm(), EPS);
    for (int j = 0; j < i; ++j){
      ASSERT_NEAR(0.f, U.col(i).dot(U.col(j)), EPS);
    }
  }

  for (int i = 0; i < V.cols(); ++i){
    ASSERT_NEAR(1.f, V.col(i).norm(), EPS);
    for (int j = 0; j < i; ++j){
      ASSERT_NEAR(0.f, V.col(i).dot(V.col(j)), EPS);
    }
  }

  ASSERT_NEAR(0.f, (A - U * S.asDiagonal() * V.transpose()).norm(), EPS);
}

TEST(redsvd, random_usv){
  const int row  = 50;
  const int col  = 30;
  const int rank = 10;
  MatrixXf U(row, rank);
  MatrixXf V(col, rank);
  REDSVD::sampleGaussianMat(U);
  REDSVD::processGramSchmidt(U);
  REDSVD::sampleGaussianMat(V);
  REDSVD::processGramSchmidt(V);
  VectorXf S(rank);
  for (int i = 0; i < rank; ++i){
    S(i) = rank-i;
  }
  MatrixXf A = U * S.asDiagonal() * V.transpose();
  REDSVD::RedSVD svdOfA(A, rank);
  

  for (int i = 0; i < rank; ++i){
    ASSERT_NEAR(1.f, fabs(U.col(i).dot(svdOfA.matrixU().col(i))), EPS);
    ASSERT_NEAR(1.f, fabs(V.col(i).dot(svdOfA.matrixV().col(i))), EPS);
    ASSERT_NEAR(0.f, S(i) - svdOfA.singularValues()(i), EPS);
  }
}

TEST(redsvd, file){
  int rowN = 20;
  int colN = 70;
  MatrixXf A = MatrixXf::Random(rowN, colN);

  int r = 20;
  REDSVD::RedSVD redsvd(A, r);
  REDSVD::RedPCA redpca(A, r);
  MatrixXf U = redsvd.matrixU();
  VectorXf S = redsvd.singularValues();
  MatrixXf V = redsvd.matrixV();
  REDSVD::writeMatrix("test", redsvd);
  REDSVD::writeMatrix("test", redpca);
}


TEST(redsvd, sparse){
  const int row = 10;
  const int col = 10;
  REDSVD::SMatrixXf A(row, col);
  REDSVD::RedSVD svdOfA(A);

  MatrixXf newA = svdOfA.matrixU() * svdOfA.singularValues().asDiagonal() * svdOfA.matrixV().transpose();
  
  for (int i = 0; i < row; ++i){
    for (int j = 0; j < col; ++j){
      ASSERT_FALSE(isnan(newA(i, j)));
    }
  }
}

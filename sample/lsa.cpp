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

#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sys/time.h>
#include "lsa.hpp"

using namespace std;
using namespace REDSVD;

LSA::LSA(){
}

LSA::~LSA(){
}

string normalize(const string& word){
  string ret;
  size_t p = 0;
  while (p < word.size() && !isalnum(word[p])){
    ++p;
  }
  for (; p < word.size(); ++p){
    if (!isalnum(word[p])) break;
    ret += tolower(word[p]);
  }
  return ret;
}

double getSec(){
  timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + (double)tv.tv_usec*1e-6;
}


void LSA::readFile(const char* fn, fv_t& fv){
  fv.clear();
  ifstream ifs(fn);
  if (!ifs){
    throw string("cannot open ") + fn;
  }

  string word;
  vector<string> words;
  while (ifs >> word){
    string nword = normalize(word);
    if (nword.size() == 0) continue;
    words.push_back(nword);
    if (words.size() == 100) break;
  }
  sort(words.begin(), words.end());
  words.erase(unique(words.begin(), words.end()), words.end());
  
  if (words.size() <= 5){
    return;
  }
  for (size_t i = 0; i <words.size(); ++i){
    fv.push_back(make_pair(getID(words[i]), 1.f));
  }
  sort(fv.begin(), fv.end());
}


void LSA::processFileList(const char* fn, const int r){
  vector<fv_t> fvs;
  ifstream ifs(fn);
  if (!ifs){
    throw string("cannot open ") + fn;
  }

  double start = getSec();
  cout << "read ... " << flush;
  for (string line; getline(ifs, line); ){
    fv_t fv;
    readFile(line.c_str(), fv);
    if (fv.size() >= 10){
      fvs.push_back(fv);
    }
  }
  cout << getSec() - start << " sec." << endl;


  SMatrixXf A;
  Util::convertFV2Mat(fvs, A);
  cout << "docN:" << A.rows() << " wordN:" << A.cols() << " wordTotalN:" << A.nonZeros() << endl;
  cout << "SVD ..." << flush;
  start = getSec();
  RedSVD svdOfA(A, r);
  cout << getSec() - start << "sec." << endl;

  svdOfA.matrixV();
}

int LSA::getID(const std::string& str){
  map<string, int>::const_iterator it = word2ID.find(str);
  if (it == word2ID.end()){
    int newID = (int)word2ID.size();
    word2ID[str] = newID;
    return newID;
  } else {
    return it->second;
  }
}

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

#ifndef LSA_HPP__
#define LSA_HPP__

#include <map>
#include <vector>
#include <string>
#include "../src/redsvd.hpp"

class LSA{
public:
  LSA();
  ~LSA();

  void processFileList(const char* fn, const int r);
private:
  void readFile(const char* fn, REDSVD::fv_t& fv);
  int getID(const std::string& str);

  std::map<std::string, int> word2ID;
};

#endif // LSA_HPP__

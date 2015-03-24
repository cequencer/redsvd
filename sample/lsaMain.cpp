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
#include "lsa.hpp"

using namespace std;

// use only for measuring the performance of RedSVD 

int main(int argc, char* argv[]){
  if (argc != 3){
    cerr << argv[0] << " filelist rank=INT" << endl;
    return -1;
  }

  int rank = atoi(argv[2]);
  if (rank <= 0){
    cerr << "rank=" << rank << endl
	 << "rank should be positive interger" << endl;
    return -1;
  }

  try {
    LSA lsa;
    lsa.processFileList(argv[1], atoi(argv[2]));
  } catch (string& errorMessage){
    cerr << errorMessage << endl;
    return -1;
  }
  return 0;
}

#ifndef SNOWUTILS_H
#define SNOWUTILS_H

#include <string>
#include <time.h>
#include <ctime>
#include <vector>
#include <unistd.h>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <stdio.h>

#include "SeqLib/SeqLibCommon.h"

#if __cplusplus > 199711L //indicates that I am using a compiler that is compatible with the standard defined in November of 1997. 
  #include <memory>
  #include <unordered_set>
  #include <unordered_map>
  #define SeqHashMap std::unordered_map
  #define SeqHashSet std::shared_ptr
  #define HAVE_C11 1 
#else

#ifdef __APPLE__
  #include <memory>
  #include <unordered_set>
  #include <unordered_map>
  #define SeqHashMap std::unordered_map
  #define SeqHashSet std::unordered_set
  #define SeqPointer std::shared_ptr
#else
  #include <tr1/memory>
  #include <tr1/unordered_set>
  #include <tr1/unordered_map>
  #define SeqHashMap std::tr1::unordered_map
  #define SeqHashSet std::tr1::unordered_set
  #define SeqPointer std::tr1::shared_ptr
#endif
#endif

namespace SeqLib {
  template<typename T>
    inline std::string tostring(T d){
  }
}


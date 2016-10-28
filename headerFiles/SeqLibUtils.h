#ifndef SNOWUTILS_H
#define SNOWUTILS_H

#include <string>
#include <time.h> // clock_gettime()
#include <ctime>
#include <vector>
#include <unistd.h>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <stdio.h>

#include "SeqLibCommon.h"

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
      std::stringstream ss;
      ss << d;
      return ss.str();
    }

    /** Check if a file is readable and exists
     * @param name Name of a file to test
     * @return File is readable and exists. 
     */
    inline bool read_access_test (const std::string& name){
      return(access(name.c_str(), R_OK) == 0);
    }

    /** Format an integer to include commas
     * @param data Number to format
     * @return String with formatted number containing commas
     */
      template<typename T>
	inline std::string AddCommas(T data) {
	std::stringstream ss;
	ss << data;
	std::string s = ss.str();
	if (s.length() > 3) // if the length of the input string is greater than 3, go ahead ...
	  for(int i = s.length()-3; i > 0; i -= 3)
	    s.insert(i, ",");
	return s;
      }

      /** Display the runtime (CPU and Wall)
       * 
       * @param start Running timer 
       * @return Time formatted as "CPU: XmYs Wall:XmYs"
       * @note Does nto work on OSX or Windows (returns "not configured")
       */
     
      /** Reverse complement in-place sequence containing upper/lower case ACTGN
       * @param a Sequence to be reverse complemented
       */
}
#endif 

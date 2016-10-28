#ifndef SEQLIB_COMMON_H
#define SEQLIB_COMMON_H

/*! \mainpage SeqLib 1.0                                                                                                         
 *                                                                                                                                                   
 * \section intro_sec Introduction                                                                                                              
 *                                                                                                                                                                                
 * SeqLib is a C++ package for querying BAM/SAM/CRAM files with HTSlib, performing                                                                                                 
 * BWA-MEM operations in memory, and performing sequence assembly with FermiKit.                                                                                                    
 * See https://github.com/walaj/SeqLib for                                                                                                                                          
 * full description                                                                                                                                                                  */


/** HTSlib/BWA-MEM/BLAT/Fermi operations */
#include <string>
#include <vector>
#include <iostream>
namespace SeqLib {
  static const char RCOMPLEMENT_TABLE[128] = {' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ', ' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ','T',' ', 'G',' ',' ',' ','C',' ',' ',' ',' ',' ',' ','N',' ',' ',' ',' ',' ','A',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ','t',' ','g',' ',' ',' ','c',' ',' ',' ',' ',' ',' ','n',' ',' ',' ',' ',' ','a',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '};
}
#endif     

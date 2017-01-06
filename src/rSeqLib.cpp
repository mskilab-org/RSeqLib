#include "SeqLib/RefGenome.h"
#include <stdlib.h>
#include "SeqLib/FastqReader.h"
#include "SeqLib/BWAWrapper.h"
#include <Rcpp.h>
#include "SeqLib/UnalignedSequence.h"
#include "SeqLib/FermiAssembler.h"

using namespace SeqLib;
using namespace Rcpp;

class BWA {
public:
  
  BWA() : bwa( new BWAWrapper) {}
  
  void make_string_index( std::string& querySeq, std::string& seqname ) {
    UnalignedSequenceVector usv = {{seqname, querySeq}};
    Rcpp::XPtr<BWAWrapper> bwa( new BWAWrapper);    
    bwa->ConstructIndex(usv);    
    this->bwa = bwa;
  }
  
  void load_fasta_index(std::string& fasta)  
  { 
    Rcpp::XPtr<BWAWrapper> bwa( new BWAWrapper);
    bwa->LoadIndex(fasta);
    this->bwa = bwa;
  }

  std::string query(std::string& query, std::string& qname, bool hardclip, 
		    double keep_sec_with_frac_of_primary_score, int max_secondary)
  {
    BamRecordVector results;
    
    //  hardclip=false, secondary score cutoff=0.9, max secondary alignments=10
    this->bwa->AlignSequence(query, qname, results, hardclip,
			     keep_sec_with_frac_of_primary_score, max_secondary);

    //    std::cout << "the size of the results object is: " << results.size() << std::endl;

    std::ostringstream stream;
    
    //  print results to stdout
    for (auto& i : results)
      stream << i << std::endl;

    std::cout << results.size() << std::endl;
    
    std::string str =  stream.str();  
    return(str);
  }

private:
  Rcpp::XPtr<BWAWrapper> bwa;
};

using namespace Rcpp ;

/** create an external pointer to a BWA object */
// [[Rcpp:>:export]]
RcppExport SEXP BWA__new(){
  return Rcpp::XPtr<BWA>( new BWA, true ) ;
}

/** make index from string*/
// [[Rcpp::export]]
void BWA__from_string( SEXP xp, std::string& querySeq, std::string& seqname) {
  Rcpp::XPtr<BWA> bwa(xp) ;
  bwa->make_string_index(querySeq, seqname);
}

/** grab index from FASTA */
// [[Rcpp::export]]
void BWA__from_fasta( SEXP xp, std::string& fasta) {
  Rcpp::XPtr<BWA> bwa(xp) ;
  bwa->load_fasta_index(fasta);
}

/** align sequence */
// [[Rcpp::export]]
std::string BWA__query( SEXP xp, std::string& qstring, std::string& qname,
		       bool hardclip, double keep_sec_with_frac_of_primary_score, int max_secondary) {  
  Rcpp::XPtr<BWA> bwa(xp) ;
  return(bwa->query(qstring, qname, hardclip,
		    keep_sec_with_frac_of_primary_score, max_secondary));
}


// ################## S
// ##################  E
// ##################   Q
// ##################    L
// ##################     I
// ##################      B


// [[Rcpp::export]]
bool fastqReader_Open(const std::string& f) {
  // initialize object of type FastqReader;
    SeqLib::FastqReader instance; 
  
  // check to see if file is available;
  bool available = instance.Open(f);

  return(available);
  //bool available FastqReader::Open(f);
}

#include "SeqLib/RefGenome.h"
#include <stdlib.h>
#include "SeqLib/FastqReader.h"
#include "SeqLib/BWAWrapper.h"
#include <Rcpp.h>
#include "SeqLib/UnalignedSequence.h"

#include "SeqLib/FermiAssembler.h"

using namespace SeqLib;
using namespace Rcpp;

// [[Rcpp::export]]

// class Uniform {
// public:
//   Uniform(double min_, double max_) : min(min_), max(max_) {}
//   NumericVector draw(int n) {
//     RNGScope scope;
//     return runif( n, min, max );
//   }
// private:
//   double min, max;
// };

// namespace Rcpp {
// template <> SEXP wrap( Uniform& );
// template <> Uniform as( SEXP ) throw(not_compatible);
// }

// namespace Rcpp {
// template <> SEXP wrap<Uniform>( Uniform& el ) {
//   int val = 10;
//   Rcpp::Language call( "new", Symbol( "Uniform" ), val ) ;
//   return call.eval() ;
// };

// template <> Uniform as<Uniform>( SEXP s ) throw(not_compatible) {
// try {
// if ( TYPEOF(s) != S4SXP ) {
// ::Rf_error( "supplied object is not of type Uniform." );
// }

// Rcpp::S4 s4obj( s );
// if ( !s4obj.is("Rcpp_Uniform" ) ) {
// ::Rf_error( "supplied object is not of type Uniform." );
// }

// Rcpp::Environment env( s4obj );
// Rcpp::XPtr<Uniform> xptr( env.get(".pointer") );


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

//bool fastqReader_GetNextSequence(struct UnalignedSequence f){
//return false; 
//}

// what I want to be able to do is create an instance of BamRecordVector and have it available in the global environment. 
// in order to create a BamRecordVector, I need a bunch of BamRecord objects. 


// [[Rcpp::export]]
void createBamVector(){
  //SeqLib::BamRecordVector instance; 
}

// [[Rcpp::export]]
void createBamRecord(){
  
}

// For some reason, the LoadIndex function is being converted to accept char rather than strings. 


// [[Rcpp::export]]
void alignSequence(const std::string& seq, const std::string& name, bool hardclip, double keep_sec_with_frac_of_primary_score, int max_secondary) {

  RefGenome ref;
  
  ref.LoadIndex("hg19.fasta");
  
  SeqLib::BWAWrapper instance;
  
  //instance.AlignSequence(name, vec, hardclip, keep_sec_with_frac_of_primary_score, max_secondary);
}

//const std::string& seq_name, const std::string& queryString, const bool& hardclip, const short& cutoff, const short& max_secondary
//void run_bwa()

// '@export
// [[Rcpp::export]]
RcppExport SEXP load_index_from_sequence(std::string& querySeq, std::string& seqname )
{
  UnalignedSequenceVector usv = {{seqname, seqname}};
  //  BWAWrapper bwa;
  // it looks like the idx is equal to NULL. 


  //  const void * address = static_cast<const void*>(&bwa);

  Rcpp::XPtr<BWAWrapper> bwa( new BWAWrapper);
  
  bwa->ConstructIndex(usv);
  
  //  std::stringstream ss;
  // ss << address;
  //return(ss.str());
  return(bwa);
}


// '@export
// [[Rcpp::export]]
std::string load_index_from_fasta(std::string& fasta)  
{
  
  BWAWrapper bwa;
  bwa.LoadIndex(fasta);

  const void * address = static_cast<const void*>(&bwa);
  std::stringstream ss;
  ss << address;
  return(ss.str());
}


// std::string run_bwa(Rcpp::XPtr<SeqLib::BWAWrapper> bwap, std::string& querySeq)  
// '@export
// [[Rcpp::export]]
std::string run_bwa(Rcpp::XPtr<int> bwap, std::string& querySeq)  
{
  // RefGenome ref;
   // ref.LoadIndex("/gpfs/commons/home/knagdimov/hg19.fasta");
   
   //  get sequence at given locus.
   // std::string seq = ref.QueryRegion("1", 1000000, 1001000);
   //   std::cout << seq << std::endl; 
   
   // Make an in-memory BWA-MEM index of region

   //   long int addy = strtol(address.c_str(), NULL, 16);
   
   //   SeqLib::BWAWrapper* bwap = reinterpret_cast<SeqLib::BWAWrapper*>(addy);
   // std::cout << addy;

   // std:: cout << bwap->NumSequences() << std::endl;
   
   //  &bwa = address;
   //   SeqLib::BWAWrapper bwa = address*
   // UnalignedSequenceVector usv = {{"chr_reg1", seq}};
   // it looks like the idx is equal to NULL.
   //   SeqLib::BWAWrapper bwa;
   //   bwap->ConstructIndex(usv);
   
   BamRecordVector results;
   // hardclip=false, secondary score cutoff=0.9, max secondary alignments=10
   //  std:: cout << bwa.NumSequences() << std::endl;
   //  bwap->AlignSequence(querySeq, "my_seq", results, false, 0.9, 10);
   
   //std::cout << "the size of the results object is: " << results.size() << std::endl;
   
   //   std::ostringstream stream;
   // print results to stdout
   //for (auto& i : results)
   // stream << i << std::endl;
   //    std::cout << i << std::endl;
   
   //std::string str =  stream.str();  
   //   return(str);
   return("FDASFDSFSD");
}


// Thursday, Dec 08, 2016 06:29:17 PM
// So what I want to be able to do 


// ################## F
// ##################  E
// ##################   R
// ##################    M
// ##################     I



// /// create an external pointer to a Uniform object
// // [[Rcpp::export]]
// RcppExport SEXP Uniform__new(SEXP min_, SEXP max_) {
//   // convert inputs to appropriate C++ types
//   double min = as<double>(min_), max = as<double>(max_);
//   // create a pointer to an Uniform object and wrap it
//   // as an external pointer
//   Rcpp::XPtr<Uniform> ptr( new Uniform( min, max ), true );
//   // return the external pointer to the R side
//   return ptr;
// }

// /// invoke the draw method
// // [[Rcpp::export]]
// RcppExport SEXP Uniform__draw( SEXP xp, SEXP n_ ) {
//   // grab the object as a XPtr (smart pointer) to Uniform
//   Rcpp::XPtr<Uniform> ptr(xp);
//   // convert the parameter to int
//   int n = as<int>(n_);
//   // invoke the function
//   NumericVector res = ptr->draw( n );
//   // return the result to R
//   return res;
// }


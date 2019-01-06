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

    //    std::cout << results.size() << std::endl;
    
    std::string str =  stream.str();  
    return(str);
  }

  void secondary_alignments(){
    
  }

private:
  Rcpp::XPtr<BWAWrapper> bwa;
};


// Fermi

class Fermi {
public:

  Fermi() : fermiAssembler(new FermiAssembler) {}

private:
  Rcpp::XPtr<FermiAssembler> fermiAssembler;
};


using namespace Rcpp ;


/** create an external pointer to a BWA object */
// [[Rcpp::export]]
RcppExport SEXP Fermi__new(){
  return Rcpp::XPtr<FermiAssembler>( new FermiAssembler, true ) ;
}

/** add reads to fermi assembler*/
// [[Rcpp::export]]
void Fermi__addReads(SEXP fermi,
                     Rcpp::StringVector qnames,
                     Rcpp::StringVector seqs,
                     Rcpp::StringVector quals
                     )
{
  Rcpp::XPtr<FermiAssembler> fermixp(fermi) ;

  UnalignedSequenceVector usv;

  for(int i=0; i < qnames.size(); i++)
    {
      UnalignedSequence read;
      read.Seq =  Rcpp::as<std::string>(seqs[i]);
      read.Name = Rcpp::as<std::string>(qnames[i]);
      read.Qual = Rcpp::as<std::string>(quals[i]);
      usv.push_back(read);
    }

  fermixp->AddReads(usv);
}


/** error correct reads*/
// [[Rcpp::export]]
void Fermi__correctReads(SEXP fermi)
{
  Rcpp::XPtr<FermiAssembler> fermixp(fermi) ;
  fermixp->CorrectReads();
}

/** perform Assembly*/
// [[Rcpp::export]]
void Fermi__performAssembly(SEXP fermi)
{
  Rcpp::XPtr<FermiAssembler> fermixp(fermi) ;
  fermixp->PerformAssembly();
}

/** get Contigs*/
// [[Rcpp::export]]
std::vector< std::string > Fermi__getContigs(SEXP fermi)
{
  const Rcpp::XPtr<FermiAssembler> fermixp(fermi) ;
  std::vector<std::string> contigs = fermixp->GetContigs();
  return(contigs);
}

/** create an external pointer to a BWA object */
// [[Rcpp::export]]
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

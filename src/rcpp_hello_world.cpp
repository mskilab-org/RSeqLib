#include <Rcpp.h>
#include "/Users/knagdimov/Documents/seqLib/headerFiles/SeqLib/SeqLibCommon.h"
#include <iostream>
#include <unistd.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List rcpp_hello_world() {

    CharacterVector x = CharacterVector::create( "foo", "bar2" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y ) ;
    std::cout << "hello" << endl;
    return z ;
}

// [[Rcpp::export]]
int bwa() {
    return 5; 
}

// [[Rcpp::export]]
int main()
{
  cout << SeqLib::RCOMPLEMENT_TABLE << endl;
  cout << __cplusplus << endl;
  cout << __APPLE__ << endl;
  cout << "hello" << endl;
  
  stringstream ss;
  ss << 100;
  cout << ss << endl;
  
  access();
}

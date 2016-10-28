#include <iostream>
#include <Rcpp.h>
#include "../headerFiles/SeqLibCommon.h"
#include "../headerFiles/SeqLibUtils.h"
using namespace Rcpp;

// [[Rcpp::export]]
List rcpp_hello_world() {

    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y ) ;
    std::cout << "hey, how are you";
    return z ;
}

// [[Rcpp::export]]
void whats_my_name(){
  std::cout << SeqLib::RCOMPLEMENT_TABLE;
  std::cout << "it's khagay" << std::endl; 
}

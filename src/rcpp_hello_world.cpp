#include <iostream>
#include <Rcpp.h>
#include "SeqLib/FermiAssembler.h"
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
  //std::cout << SeqLib::RCOMPLEMENT_TABLE;
  //std::cout << "it's khagay" << std::endl;
  SeqLib::FermiAssembler f;
}

// [[Rcpp::export]]
void start_fermi(){
  SeqLib::print_my_name();
  
}

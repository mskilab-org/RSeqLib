#include <iostream>
#include "SeqLib/FermiAssembler.h"

using namespace SeqLib;

// [[Rcpp::export]]
void whats_my_name(){
  //std::cout << SeqLib::RCOMPLEMENT_TABLE;
  std::cout << "it's khagay" << std::endl;
  //MAG_MIN_NSR_COEF .1
  //if(MAG_MIN_NSR_COEF == .1)
  //std::cout << "Hey, what's up!!!!!" << std::endl;
}

// [[Rcpp::export]]
void start_fermi(){
}

// [[Rcpp::export]]
void fermi() {
  SeqLib::FermiAssembler f;
}

int main(){
}

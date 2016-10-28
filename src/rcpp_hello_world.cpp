#include <Rcpp.h>

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

#include "Benchmarks.h"

Benchmarks::Benchmarks(){
  n_dim = 100;
  x_min = -100.0;
  x_max = +100.0;
}


Benchmarks::~Benchmarks(){
  /* emtpy */
}

double Benchmarks::getMin(){
  return x_min;
}

double Benchmarks::getMax(){
  return x_max;
}

uint Benchmarks::getID(){
  return ID;
}

void Benchmarks::setMin( const double _min ){
  x_min = _min;
}

void Benchmarks::setMax( const double _max ){
  x_max = _max;
}

void Benchmarks::setDim( const uint _dim ){
  n_dim = _dim;
}

#ifndef __jDE__
#define __jDE__

#include <tuple>
#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <random>
#include <algorithm>
#include <cassert>

#include "Benchmarks.h"

typedef std::vector<double> vDouble;

class jDE {

private:
  double f_lower;
  double f_upper;
  double T;

  vDouble F;
  vDouble CR;

  uint size;

  std::mt19937 rng;
public:
  jDE( uint );
  ~jDE();

  double getFLower();
  double getFUpper();
  double getT();

  void setFLower( const double );
  void setFUpper( const double );
  void setT( const double );

  /* update the F and CR values */
  void update();

  /* print all F values */
  void showF();

  /* print all CR values */
  void showCR();

  /*
   * Performs a iteration of DE/rand/1/bin using
   * jDE self-adaptation.
   * @params:
   *    - uint: num of dimensions
   *    - uint: pop size
   *    - vector<double> containing the genes(ps x ndim)
   *    - vector<double> to put the new generation
   *    - Benchmarks to check bounds
   */
  void runDE(uint, uint, const vDouble &, vDouble &, Benchmarks *);

  /*
   * Performs selection on the new offsprings
   * @params:
   *    - uint: num of dimensions
   *    - uint: pop size
   *    - vector<double> containing the genes(ps x ndim)
   *    - vector<double> containing the new genes(ps x ndim)
   *    - vector<double> containing the fitness of previous generation
   *    - vector<double> containing the fitness of new generation
   */
  void selection(uint, uint, vDouble &, const vDouble &, vDouble& fitness, const vDouble & );

};

#endif

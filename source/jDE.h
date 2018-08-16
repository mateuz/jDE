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
   */
  void runDE(uint, uint, const vDouble &, vDouble &);

};

#endif

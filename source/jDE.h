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

};

#endif

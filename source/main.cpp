#include <iostream>
#include <tuple>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <sys/time.h>
#include <cassert>

#include "jDE.h"

int main(int argc, char * argv[]){
  uint size = 20;
  uint ndim = 10;

  jDE * j = new jDE(size);
  j->update();

  //j->showF();
  //j->showCR();

  //j->update();

  //j->showF();
  //j->showCR();

  vDouble gen( size * ndim, 0.0);
  vDouble n_gen( size * ndim, 0.0);

  j->runDE(ndim, size, gen, n_gen);

  return 0;
}

#include "F2.h"
#include "IO.h"

#include <fstream>
#include <iostream>

F2::F2(uint _dim):Benchmarks()
{
  n_dim = _dim;

  ID = 2;
  x_min = -100.0;
  x_max = +100.0;

  /* ---------------------------------------------- */
  /* Load a shift vector to test the bench function */
  std::string file_name = "data-files/shift_rosenbrock.mat";
  std::string vec_name = "Shift - Rosenbrock [-100.0, +100.0]";
  IO * io = new IO();
  std::ifstream file(file_name);
  if( not file.is_open() ){
    std::cout << "Error opening file\n";
    exit(-1);
  }
  shift = io->load_vector<double>( vec_name, file ) ;
  file.close();
  /* ---------------------------------------------- */
}

F2::~F2(){
  /* empty */
}

double F2::compute(const vDouble gen, const uint ip){
  double s = 0.0, a, b, t;
  for( uint it = 0; it < (n_dim - 1); it++ ){
    a = gen[ip + it] - shift[it] + 1.00;
    b = gen[ip + it + 1] - shift[it + 1] + 1.00;

    t = (b - (a * a));
    s += (100.0 * t * t);
    t = (a - 1.0);
    s += (t * t);
  }
  return s;
}

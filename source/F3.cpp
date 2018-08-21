#include "F3.h"
#include "IO.h"

#include <fstream>
#include <iostream>

F3::F3(uint _dim):Benchmarks()
{
  n_dim = _dim;

  ID = 3;
  x_min = -600.0;
  x_max = +600.0;

  /* ---------------------------------------------- */
  /* Load a shift vector to test the bench function */
  std::string file_name = "data-files/shift_griewank.mat";
  std::string vec_name = "Shift - Griewank [-600.0, +600.0]";
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

F3::~F3(){
  /* empty */
}

double F3::compute(const vDouble gen, const uint ip){
  double s1 = 0.0, s2 = 1.0, z;
  for( uint it = 0; it < n_dim; it++ ){
    z = gen[ip + it] - shift[it];
    s1 += z * z;
    s2 *= cos(z/sqrt(it+1));
  }
  s1 /= 4000.0;
  return (s1 - s2 + 1.0);
}

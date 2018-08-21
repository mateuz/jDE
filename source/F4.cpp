#include "F4.h"
#include "IO.h"

#include <fstream>
#include <iostream>

F4::F4(uint _dim):Benchmarks()
{
  n_dim = _dim;

  ID = 4;
  x_min = -5.0;
  x_max = +5.0;

  /* ---------------------------------------------- */
  /* Load a shift vector to test the bench function */
  std::string file_name = "data-files/shift_rastrigin.mat";
  std::string vec_name = "Shift - Rastrigin [-5.0, +5.0]";
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

F4::~F4(){
  /* empty */
}

double F4::compute(const vDouble gen, const uint ip){
  double s = 0.0, z;
  for( uint it = 0; it < n_dim; it++ ){
    z = gen[ip + it] - shift[it];
    s += pow(z, 2.00) - 10.0 * cos( 2.0 * PI * z) + 10.0;
  }
  return s;
}

#include "F1.h"
#include "IO.h"

#include <fstream>
#include <iostream>

F1::F1(uint _dim):Benchmarks()
{
  n_dim = _dim;

  ID = 1;
  x_min = -100.0;
  x_max = +100.0;

  /* ---------------------------------------------- */
  /* Load a shift vector to test the bench function */
  std::string file_name = "data-files/shift_sphere.mat";
  std::string vec_name = "Shift - Sphere";
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

F1::~F1(){
  /* empty */
}

/*
 * vDouble represents the gene
 * uint start point of the individual (initial point)
 */
double F1::compute(const vDouble gen, const uint ip){
  double s = 0.0, z;
  for( uint it = 0; it < n_dim; it++ ){
    z = gen[ip + it] - shift[it];
    s += z * z;
  }
  return s;
}

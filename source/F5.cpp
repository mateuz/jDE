#include "F5.h"
#include "IO.h"

#include <fstream>
#include <iostream>
#include <cmath>

/*
 * Shifted and Rotated Rosenbrock's Function
 *
 * as defined in "Problem Definitions and Evaluation Criteria for the
 * CEC 2013 Special Session and Competition on Real-Parameter Optimization",
 * by Liang, J.J., Qu, B.-Y., Suganthan, P.N., Hernandez-Diaz, A.G.,
 * Computational Intelligence Laboratory, Zhengzhou University, Zhengzhou,
 * China and Nanyang Technological University, Singapore, Technical Report,
 * v. 2012, p. 3-18, 2013.
*/


F5::F5(uint _dim):Benchmarks()
{
  n_dim = _dim;

  ID = 5;
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

  /* ---------------------------------------------- */
  /* Load a rotate matrix                           */
  check(n_dim);
  file_name = "data-files/rot/M_D" + std::to_string(n_dim) + ".txt";
  vec_name = "M_D" + std::to_string(n_dim);

  file.open(file_name, std::ifstream::in);
  if( not file.is_open() ){
    std::cout << "Error opening rotation matrix file\n";
    exit(-1);
  }

  M = io->load_vector<double>( vec_name, file ) ;
  file.close();
}

F5::~F5(){
  /* empty */
}

double F5::compute(const vDouble gen, const uint ip){
  //double c = 2.048/100.0;
  uint i, j;
  std::vector<double> z(n_dim);

  //shift and multiply by 2.048/100.0
  for( i = 0; i < n_dim; i++ )
    z[i] = (gen[ip + i] - shift[i]);// * c;

  //rotate function
  std::vector<double> z_rot(n_dim);
  for( i = 0; i < n_dim; i++ ){
    z_rot[i] = 0.0;

    for( j = 0; j < n_dim; j++ )
			z_rot[i] += z[j] * M[i * n_dim + j];

    z_rot[i] += 1.0;
  }

  //calc fitness
  double s1 = 0.0, t1, t2;
  for( uint i = 0; i < (n_dim - 1); i++ ){
    t1 = (z_rot[i] * z_rot[i]) - z_rot[i + 1];
    t2 = z_rot[i] - 1.0;
    s1 += 100.0 * (t1 * t1) + (t2 * t2);
  }
  return (s1);
}

#include "F6.h"
#include "IO.h"

#include <fstream>
#include <iostream>
#include <cmath>

/*
 * Shifted and Rotated Griewank's Function
 *
 * as defined in "Problem Definitions and Evaluation Criteria for the
 * CEC 2013 Special Session and Competition on Real-Parameter Optimization",
 * by Liang, J.J., Qu, B.-Y., Suganthan, P.N., Hernandez-Diaz, A.G.,
 * Computational Intelligence Laboratory, Zhengzhou University, Zhengzhou,
 * China and Nanyang Technological University, Singapore, Technical Report,
 * v. 2012, p. 3-18, 2013.
*/


F6::F6(uint _dim):Benchmarks()
{
  n_dim = _dim;

  ID = 6;
  x_min = -100.0;
  x_max = +100.0;

  /* ---------------------------------------------- */
  /* Load a shift vector to test the bench function */
  std::string file_name = "data-files/shift_griewank.mat";
  std::string vec_name = "Shift - Griewank";
  IO * io = new IO();
  std::ifstream file(file_name);
  if( not file.is_open() ){
    std::cout << "Error opening data-files/shift_griewank.mat file\n";
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
  /*
  for (uint i=0, j=1; i< n_dim; i++, j++){
    printf("%+14.10lf ", shift[i]);
    if( j == 5 ){
      printf("\n");
      j = 0;
    }
  }
  */

  M = io->load_vector<double>( vec_name, file ) ;
  file.close();

  /*
  for (uint i=0, j=1; i< n_dim * n_dim; i++, j++){
    printf("%+14.10lf ", M[i]);
    if( j == 5 ){
      printf("\n");
      j = 0;
    }
  }
  */
}

F6::~F6(){
  /* empty */
}

double F6::compute(const vDouble gen, const uint ip){
  const double c = 600.0/100.0;

  std::vector<double> z(n_dim);
  for(uint i = 0; i < n_dim; i++ )
    z[i] = (gen[ip + i] - shift[i]) * c;

  //rotate function
  std::vector<double> z_rot(n_dim);
  for(uint i = 0; i < n_dim; i++ ){
    z_rot[i] = 0.0;
    for( uint j = 0; j < n_dim; j++ )
      z_rot[i] += z[j] * M[i*n_dim+j];
  }

  for( int i = 0; i < n_dim; i++ ){
    z[i] = z_rot[i] * pow(100.0,1.0*i/(n_dim-1)/2.0);
  }

  double s = 0.0, p = 1.0;
  for( uint i = 0; i < n_dim; i++ ){
    s += z[i] * z[i];
    p *= cos(z[i] / sqrt(i + 1.0));
  }
  return (1.0 + s/4000.0 - p);
}


double F6::compute(const double * gen, const uint ip){
  const double c = 600.0/100.0;

  std::vector<double> z(n_dim);
  for(uint i = 0; i < n_dim; i++ )
    z[i] = (gen[ip + i] - shift[i]) * c;

  //rotate function
  std::vector<double> z_rot(n_dim);
  for(uint i = 0; i < n_dim; i++ ){
    z_rot[i] = 0.0;
    for( uint j = 0; j < n_dim; j++ )
      z_rot[i] += z[j] * M[i*n_dim+j];
  }

  for( int i = 0; i < n_dim; i++ ){
    z[i] = z_rot[i] * pow(100.0,1.0*i/(n_dim-1)/2.0);
  }

  double s = 0.0, p = 1.0;
  for( int i = 0; i < n_dim; i++ ){
    s += z[i] * z[i];
    p *= cos(z[i] / sqrt(1.0+i));
  }
  return (1.0 + s/4000.0 - p);// - 500.0; //the shift point!!!
}

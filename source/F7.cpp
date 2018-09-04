#include "F7.h"
#include "IO.h"

#include <fstream>
#include <iostream>
#include <cmath>

/*
 * Shifted and Rotated Rastrigin's Function
 *
 * as defined in "Problem Definitions and Evaluation Criteria for the
 * CEC 2013 Special Session and Competition on Real-Parameter Optimization",
 * by Liang, J.J., Qu, B.-Y., Suganthan, P.N., Hernandez-Diaz, A.G.,
 * Computational Intelligence Laboratory, Zhengzhou University, Zhengzhou,
 * China and Nanyang Technological University, Singapore, Technical Report,
 * v. 2012, p. 3-18, 2013.
*/


F7::F7(uint _dim):Benchmarks()
{
  n_dim = _dim;

  ID = 7;
  x_min = -5.0;
  x_max = +5.0;

  /* ---------------------------------------------- */
  /* Load a shift vector to test the bench function */
  std::string file_name = "data-files/shift_rastrigin.mat";
  std::string vec_name = "Shift - Rastrigin [-5.0, +5.0]";
  IO * io = new IO();
  std::ifstream file(file_name);
  if( not file.is_open() ){
    std::cout << "Error opening data-files/shift_rastrigin.mat file\n";
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

F7::~F7(){
  /* empty */
}

double F7::compute(const vDouble gen, const uint ip){
  const uint d = (n_dim - 1);

  const double beta = 0.2;
  const double alpha = 10.0;
  const double c = 5.12/100.0;

  std::vector<double> z(n_dim);
  for(uint i = 0; i < n_dim; i++ )
    z[i] = (gen[ip + i] - shift[i]) * c;

  //first rotation
  std::vector<double> z_rot(n_dim);
  for(uint i = 0; i < n_dim; i++ ){
    z_rot[i] = 0.0;

    for( uint j = 0; j < n_dim; j++ )
      z_rot[i] += z[j] * M[i * n_dim + j];
  }

  //osz
  double xx;
  int sx;
  std::vector<double> xosz(n_dim);
  for( uint i = 0; i < n_dim; i++ ){
    if( i == 0 || i == (n_dim - 1) ){
      if( x[i] != 0 )
        xx = log(fabs( z[i] ));

      if( x[i] > 0 ){
        c1 = 10;
        c2 = 7.9;
      } else {
        c1 = 5.5;
        c2 = 3.1;
      }

      if( x[i] > 0 )
        sx = 1;
      else if( x[i] == 0 )
        sx = 0;
      else
        sx = -1;

      xosz[i] = sx * exp(xx + 0.049 * (sin(c1 * xx) + sin(c2 * xx)));
    } else {
      xosz[i] = x[i];
    }
  }

  //asy
  std::vector<double> xasy(n_dim);
  for( uint i = 0; i < n_dim; i++ ){
    if( x[i] > 0 )
      xasy[i] = pow(x[i], 1.0 + beta * i / d * pow(x[i], 0.5));
  }

  //second rotation

  //alpha
  for( uint i = 0; i < n_dim; i++ ){
    z[i] *= pow(alpha, (1.0 * i) / (d * 2.0));
  }

  //third rotation

  //evaluation
  double s = 0.0;
  for( uint i = 0; i < n_dim; i++ )
    s += (z[i] * z[i] - 10.0 * cos(2.0 * PI * z[i]) + 10.0);

  return s;
}

#include <iostream>
#include <tuple>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <sys/time.h>
#include <cassert>
#include <cmath>

#include "IO.h"
#include "jDE.h"

#define PI (3.141592653589793238462643383279)

/*
 * Shift - Sphere [-100.0, +100.0]
 */
double shift_sphere(const vDouble gen, const vDouble shift, const uint p, const uint ndim){
  double s = 0.0;
  for( uint it = 0; it < ndim; it++ )
    s += pow(gen[p + it] - shift[it], 2);
  return s;
}

/*
 * Shift - Rosenbrock [-100.0, +100.0]
 */
double shift_rosenbrock(const vDouble gen, const vDouble shift, const uint p, const uint ndim){
  double s = 0.0, a, b;
  for( uint it = 0; it < (ndim - 1); it++ ){
    a = gen[p + it] - shift[it] + 1.00;
    b = gen[p + it + 1] - shift[it + 1] + 1.00;
    s += 100.0*pow(pow(a, 2.00)-b,2.00)+pow(a-1.0,2.00);
  }
  return s;
}

/*
 * Shift - Griewank [-600.0, +600.0]
 */
double shift_griewank(const vDouble gen, const vDouble shift, const uint p, const uint ndim){
  double s1 = 0.0, s2 = 1.0, z;
  for( uint it = 0; it < ndim; it++ ){
    z = gen[p + it] - shift[it];
    s1 += z * z;
    s2 *= cos(z/sqrt(it+1));
  }
  s1 /= 4000.0;
  return (s1 - s2 + 1.0);
}

/*
 * Shift - Rastrigin [-5.0, +5.0]
 */
double shift_rastrigin(const vDouble gen, const vDouble shift, const uint p, const uint ndim){
    double s = 0.0, z;
    for( uint it = 0; it < ndim; it++ ){
      z = gen[p + it] - shift[it];
      s += pow(z, 2.00) - 10.0 * cos( 2 * PI * z) + 10.0;
    }
}

double compute(
  const vDouble gen,
  const vDouble shift,
  const uint ip,
  const uint ndim,
  double (*f)(const vDouble, const vDouble, const uint, const uint))
{
    return (*f)(gen, shift, ip, ndim);
}

int main(int argc, char * argv[]){
  uint NP = 20;
  uint ndim = 100;
  uint nEvals = 1000000;

  vDouble gen(NP * ndim);
  vDouble n_gen(NP * ndim);
  vDouble fitness(NP, 0.0);
  vDouble n_fitness(NP, 0.0);

  double x_min = -600.0;
  double x_max = +600.0;

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
  vDouble shift = io->load_vector<double>( vec_name, file ) ;
  file.close();
  //vDouble shift(ndim, 0.00);
/* ---------------------------------------------- */

  std::mt19937 rng;
  rng.seed(std::random_device{}());
  std::uniform_real_distribution<double> random(x_min, x_max);

  //randomly init genes
  for( auto it = gen.begin(); it != gen.end(); it++ )
    *it = random(rng);

  jDE * jde = new jDE(NP);

  //evaluate the first offspring
  for( uint i = 0; i < NP; i++ ){
    fitness[i] = compute(gen, shift, i * ndim, ndim, shift_griewank);
    //printf("%-3.2lf %-3.2lf = %-5.20E\n", gen[i * ndim], gen[i * ndim + 1], fitness[i]);
  }
  for( uint run = 0; run < nEvals; run += NP){
    jde->runDE(ndim, NP, gen, n_gen, x_min, x_max);

    for( uint i = 0; i < NP; i++ )
        n_fitness[i] = compute(n_gen, shift, i * ndim, ndim, shift_griewank);

    jde->selection(ndim, NP, gen, n_gen, fitness, n_fitness);
    //jde->update();
  }
  double a = fitness[0];
  uint pb = 0;
  for( uint i = 0; i < NP; i++ ){
    if( fitness[i] < a){
      pb = i;
      a = fitness[i];
    }
  }
  printf("%-3i | %-5.20E\n", pb, fitness[pb]);
  return 0;
}

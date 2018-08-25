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

#include "Benchmarks.h"
#include "F1.h"
#include "F2.h"
#include "F3.h"
#include "F4.h"

double stime(){
  struct timeval tv;
  struct timezone tz;
  gettimeofday(&tv, &tz);
  double mlsec = 1000.0 * ((double)tv.tv_sec + (double)tv.tv_usec/1000000.0);
  return mlsec/1000.0;
}

void show_params(
	uint n_runs,
	uint NP,
	uint n_evals,
	uint n_dim,
	std::string FuncObj
){
	printf(" | Number of Executions:                    %d\n", n_runs);
	printf(" | Population Size:                         %d\n", NP);
	printf(" | Number of Dimensions:                    %d\n", n_dim);
	printf(" | Number of Function Evaluations:          %d\n", n_evals);
	printf(" | Optimization Function:                   %s\n", FuncObj.c_str());
}

std::string toString(uint id){
  switch( id ){
    case 1:
      return "Shifted Sphere";
    case 2:
      return "Shifted Rosenbrock";
    case 3:
      return "Shifted Griewank";
    case 4:
      return "Shifted Rastringin";
    default:
      return "Unknown";
  }
}

int main(int argc, char * argv[]){
  uint NP = 50;
  uint ndim = 100;
  uint n_evals = 10000 * ndim;
  vDouble gen(NP * ndim);
  vDouble n_gen(NP * ndim);
  vDouble fitness(NP, 0.0);
  vDouble n_fitness(NP, 0.0);

  std::mt19937 rng;
  rng.seed(std::random_device{}());

  Benchmarks * B = new F2(ndim);

  uint n_runs = 10;
  uint f_id = 2;

  printf(" +==============================================================+ \n");
	printf(" |                      EXECUTION PARAMETERS                    | \n");
	printf(" +==============================================================+ \n");
	show_params(n_runs, NP, n_evals, ndim, toString(f_id));
	printf(" +==============================================================+ \n");

  std::uniform_real_distribution<double> random(B->getMin(), B->getMax());

  double tini, tend;
  for( int go = 1; go <= n_runs; go++ ){
    tini = stime();
    //randomly init genes
    for( auto it = gen.begin(); it != gen.end(); it++ )
      *it = random(rng);

    jDE * jde = new jDE(NP);

    for( uint i = 0; i < NP; i++ )
      fitness[i] = B->compute(gen, i * ndim);

    for( uint run = 0; run < n_evals; run += NP){
      jde->runDE(ndim, NP, gen, n_gen, B->getMin(), B->getMax());

      for( uint i = 0; i < NP; i++ )
        n_fitness[i] = B->compute(n_gen, i * ndim);

      jde->selection(ndim, NP, gen, n_gen, fitness, n_fitness);
      jde->update();
    }
    double a = fitness[0];
    uint pb = 0;
    for( uint i = 0; i < NP; i++ ){
      if( fitness[i] < a){
        pb = i;
        a = fitness[i];
      }
    }
    tend = stime();
    printf(" | Execution: %-2d Overall Best: %+.8lf Time(ms): %.8f\n", go, fitness[pb], tend-tini);
  }
  return 0;
}

#include <iostream>
#include <tuple>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <sys/time.h>
#include <cassert>
#include <cmath>
#include <chrono>

#include "IO.h"
#include "jDE.h"

#include "Benchmarks.h"
#include "F1.h"
#include "F2.h"
#include "F3.h"
#include "F4.h"
#include "F5.h"
#include "F6.h"
#include "F7.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

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
    case 5:
      return "Rotated Rosenbrock";
    case 6:
      return "Rotated Griewank";
    case 7:
      return "Rotated Rastrigin";
    default:
      return "Unknown";
  }
}

Benchmarks * getFunction(uint id, uint n_dim){
  Benchmarks * n;

  if( id == 1 ){
    n = new F1(n_dim);
    return n;
  }

  if( id == 2 ){
    n = new F2(n_dim);
    return n;
  }

  if( id == 3 ){
    n = new F3(n_dim);
    return n;
  }

  if( id == 4 ){
    n = new F4(n_dim);
    return n;
  }

  if( id == 5 ){
    n = new F5(n_dim);
    return n;
  }

  if( id == 6 ){
    n = new F6(n_dim);
    return n;
  }

  if( id == 7 ){
    n = new F7(n_dim);
    return n;
  }

  return NULL;
}

int main(int argc, char * argv[]){
  srand((unsigned)time(NULL));
  uint n_runs, NP, max_evals, n_dim, f_id;

  try{
    po::options_description config("Opções");
    config.add_options()
    ("runs,r"    , po::value<uint>(&n_runs)->default_value(1)    , "Number of Executions"          )
    ("pop_size,p", po::value<uint>(&NP)->default_value(20)       , "Population Size"               )
    ("dim,d"     , po::value<uint>(&n_dim)->default_value(10)    , "Number of Dimensions"          )
    ("func_obj,o", po::value<uint>(&f_id)->default_value(1)      , "Function to Optimize"          )
    ("max_eval,e", po::value<uint>(&max_evals)->default_value(10e5), "Number of Function Evaluations")
    ("help,h", "Mostrar texto de Ajuda");

    po::options_description cmdline_options;
    cmdline_options.add(config);
    po::variables_map vm;
    store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
    notify(vm);
    if( vm.count("help") ){
      std::cout << cmdline_options << "\n";
      return 0;
    }
  }catch(std::exception& e){
  std::cout << e.what() << "\n";
  return 1;
  }

	printf(" +==============================================================+ \n");
	printf(" |                      EXECUTION PARAMETERS                    | \n");
	printf(" +==============================================================+ \n");
	show_params(n_runs, NP, max_evals, n_dim, toString(f_id));
	printf(" +==============================================================+ \n");

  Benchmarks * B = NULL;
  B = getFunction(f_id, n_dim);

  if( B == NULL ){
    std::cout << "Unknown function! Exiting...\n";
    return 0;
  }

	std::vector< std::pair<double, double> > stats;

  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> random(0, 1);//[0, 1)

  const double x_min = B->getMin();
  const double x_max = B->getMax();

  jDE * jde = new jDE(NP);

  double tini, tend, best;
  for( int go = 1; go <= n_runs; go++ ){
    vDouble gen(NP * n_dim);
    vDouble n_gen(NP * n_dim);
    vDouble fitness(NP, 0.0);
    vDouble n_fitness(NP, 0.0);
    //randomly init genes
    for( uint i = 0; i < gen.size(); i++)
      gen[i] = (( x_max - x_min) * random(rng)) + x_min;

    //evaluate the initial population
    for( uint i = 0; i < NP; i++ )
      fitness[i] = B->compute(gen, i * n_dim);

    uint nfes = 0;
    tini = stime();
    //start the evolutive process;
    while(nfes < max_evals){
      jde->update();
      jde->runDE(n_dim, NP, gen, n_gen, x_min, x_max);
      for( int i = 0; i < NP; i++ ){
        n_fitness[i] = B->compute(n_gen, i * n_dim);
        nfes++;
        if( nfes >= max_evals ) break;
      }
      jde->selection(n_dim, NP, gen, n_gen, fitness, n_fitness);
    }
    best = *std::min_element(fitness.begin(),fitness.end());
    tend = stime();

    jde->reset();
    printf(" | Execution: %-2d Overall Best: %+.8lf Time(s): %.8f\n", go, best, tend-tini);
		stats.push_back(std::make_pair(best, tend-tini));
  }

	/* Statistics of the Runs */
	double FO_mean  = 0.0f, FO_std  = 0.0f;
	double T_mean   = 0.0f, T_std   = 0.0f;
	for( auto it = stats.begin(); it != stats.end(); it++){
		FO_mean += it->first;
		T_mean  += it->second;
	}
	FO_mean /= n_runs;
	T_mean  /= n_runs;
	for( auto it = stats.begin(); it != stats.end(); it++){
		FO_std += (( it->first - FO_mean )*( it->first  - FO_mean ));
		T_std  += (( it->second - T_mean )*( it->second - T_mean  ));
	}
	FO_std /= n_runs;
	FO_std = sqrt(FO_std);
	T_std /= n_runs;
	T_std = sqrt(T_std);
	printf(" +==============================================================+ \n");
	printf(" |                     EXPERIMENTS RESULTS                      | \n");
	printf(" +==============================================================+ \n");
	printf(" | Objective Function:\n");
	printf(" | \t mean:         %+.20E\n", FO_mean);
	printf(" | \t std:          %+.20E\n", FO_std);
	printf(" | Execution Time (ms): \n");
	printf(" | \t mean:         %+.3lf\n", T_mean);
	printf(" | \t std:          %+.3lf\n", T_std);
	printf(" +==============================================================+ \n");

  //delete jde;
  return 0;
}

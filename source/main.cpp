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

int main(int argc, char * argv[]){

  uint NP = 50;
  uint ndim = 100;
  uint nEvals = 10e4 * ndim;
  vDouble gen(NP * ndim);
  vDouble n_gen(NP * ndim);
  vDouble fitness(NP, 0.0);
  vDouble n_fitness(NP, 0.0);

  std::mt19937 rng;
  rng.seed(std::random_device{}());

  Benchmarks * B;
  for( uint go = 1; go <= 4; go++ ){
    switch (go) {
      case 1:
        B = new F1(ndim);
        break;
      case 2:
        B = new F2(ndim);
        break;
      case 3:
        B = new F3(ndim);
        break;
      case 4:
        B = new F4(ndim);
        break;
    }

    std::uniform_real_distribution<double> random(B->getMin(), B->getMax());

    //randomly init genes
    for( auto it = gen.begin(); it != gen.end(); it++ )
      *it = random(rng);

    jDE * jde = new jDE(NP);

    //evaluate the first offspring
    for( uint i = 0; i < NP; i++ ){
      fitness[i] = B->compute(gen, i * ndim);
      //printf("%-3.2lf %-3.2lf = %-5.20E\n", gen[i * ndim], gen[i * ndim + 1], fitness[i]);
    }
    for( uint run = 0; run < nEvals; run += NP){
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
    printf("F :: %d best individual %-3i with %-5.20E.\n", B->getID(), pb, fitness[pb]);
  }
  return 0;
}

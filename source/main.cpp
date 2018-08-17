#include <iostream>
#include <tuple>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <sys/time.h>
#include <cassert>

#include "Header.h"
#include "jDE.h"

int main(int argc, char * argv[]){
  uint NP = 20;
  uint ndim = 10;
  uint funID = 1;
  uint nEvals = 10000;

  Benchmarks * b = generateFuncObj( funID );
  b->setDimension(ndim);

  vDouble gen( NP * ndim);
  vDouble n_gen( NP * ndim);
  vDouble fitness(NP, 0.0);
  vDouble n_fitness(NP, 0.0);

  std::mt19937 rng;
  rng.seed(std::random_device{}());
  std::uniform_real_distribution<double> random(b->getMinX(), b->getMaxX());

  //randomly init genes
  for( auto it = gen.begin(); it != gen.end(); it++ )
    *it = random(rng);


  jDE * jde = new jDE(NP);

  //evaluate the first offspring
  for( uint i = 0; i < NP; i++ )
    fitness[i] = b->compute(gen, i * ndim);


  double a = 9999999999.0;
  uint pb = 0;
  for( uint run = 0; run < nEvals; run += NP){
    jde->runDE(ndim, NP, gen, n_gen, b);

    for( uint i = 0; i < NP; i++ )
        n_fitness[i] = b->compute(n_gen, i * ndim);

    jde->selection(ndim, NP, gen, n_gen, fitness, n_fitness);
    jde->update();

    for( uint i = 0; i < NP; i++ ){
      if( fitness[i] < a){
        pb = i;
        a = fitness[i];
      }
    }
    printf("%-3i | %-5.20E\n", run/NP, fitness[pb]);
  }
  return 0;
}

Benchmarks* generateFuncObj(int funcID){
	Benchmarks *fp;
	// run each of specified function in "configure.ini"
	if (funcID==1){
		fp = new F1();
	}else if (funcID==2){
		fp = new F2();
	}else if (funcID==3){
		fp = new F3();
	}else if (funcID==4){
		fp = new F4();
	}else if (funcID==5){
		fp = new F5();
	}else if (funcID==6){
		fp = new F6();
	}else if (funcID==7){
		fp = new F7();
	}else if (funcID==8){
		fp = new F8();
	}else if (funcID==9){
		fp = new F9();
	}else if (funcID==10){
		fp = new F10();
	}else if (funcID==11){
		fp = new F11();
	}else if (funcID==12){
		fp = new F12();
	}else if (funcID==13){
		fp = new F13();
	}else if (funcID==14){
		fp = new F14();
	}else if (funcID==15){
		fp = new F15();
	}else if (funcID==16){
		fp = new F16();
	}else if (funcID==17){
		fp = new F17();
	}else if (funcID==18){
		fp = new F18();
	}else if (funcID==19){
		fp = new F19();
	}else if (funcID==20){
		fp = new F20();
	}else{
		cerr<<"Fail to locate Specified Function Index"<<endl;
		exit(-1);
	}
	return fp;
}

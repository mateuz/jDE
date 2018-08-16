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
  uint size = 20;
  uint ndim = 10;
  uint funID = 1;

  jDE * j = new jDE(size);
  j->update();

  //j->showF();
  //j->showCR();

  //j->update();

  //j->showF();
  //j->showCR();

  vDouble gen( size * ndim, +100.0);
  vDouble n_gen( size * ndim, 0.0);

  vDouble fitness(size);
  vDouble n_fitness(size);

  Benchmarks * b = generateFuncObj( funID );
  b->setDimension(ndim);

  j->runDE(ndim, size, gen, n_gen, b);

  for( uint i = 0; i < size; i++ ){
    fitness[i] = b->compute(gen, i * ndim);
    printf("%.2lf\n", fitness[i]);
  }

  for( uint i = 0; i < size; i++ ){
    n_fitness[i] = b->compute(n_gen, i * ndim);
    printf("%.2lf\n", fitness[i]);
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

#include <iostream>
#include <tuple>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <sys/time.h>
#include <cassert>

#include "jDE.h"

int main(int argc, char * argv[]){
  jDE * j = new jDE(10);
  j->update();

  j->showF();
  j->showCR();

  j->update();

  j->showF();
  j->showCR();
  return 0;
}

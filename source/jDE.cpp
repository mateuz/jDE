#include "jDE.h"

jDE::jDE( uint _s ):
  size(_s)
{
  rng.seed(std::random_device{}());
  F.assign(size, 0.5);
  CR.assign(size, 0.9);

  f_lower = 0.1;
  f_upper = 0.9;
  T = 0.1;
}

jDE::~jDE(){
  /* empty */
}

double jDE::getFLower(){
  return f_lower;
}

double jDE::getFUpper(){
  return f_upper;
}

double jDE::getT(){
  return T;
}

void jDE::setFLower( const double FL ){
  f_lower = FL;
}

void jDE::setFUpper( const double FU ){
  f_upper = FU;
}

void jDE::setT( const double _T){
  T = _T;
}

void jDE::update(){
  assert( F.size() > 0 || CR.size() > 0);
  assert( F.size() == CR.size() );

  std::uniform_real_distribution<double> random(0.00, 1.00);

  double r1, r2, r3, r4;

  for( uint i = 0; i < size; i++ ){
    r1 = random(rng);
    r2 = random(rng);
    r3 = random(rng);
    r4 = random(rng);
    printf("r1 %.2lf r2 %.2lf r3 %.2lf r4 %.2lf\n", r1, r2, r3, r4);
    if( r2 < T )
      F[i] = f_lower + (r1 * f_upper);

    if( r4 < T )
      CR[i] = r3;
  }
}

void jDE::showF(){
  std::cout << "# F:" << std::endl;
  uint i = 0;
  for( auto it = F.begin(); it != F.end(); it++, i++ )
    std::cout << i << ": " << *it << std::endl;
}

void jDE::showCR(){
  std::cout << "# CR:" << std::endl;
  uint i = 0;
  for( auto it = CR.begin(); it != CR.end(); it++, i++ )
    std::cout << i << ": " << *it << std::endl;
}

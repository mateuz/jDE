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
    //printf("r1 %.2lf r2 %.2lf r3 %.2lf r4 %.2lf\n", r1, r2, r3, r4);
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

void jDE::runDE(uint ndim, uint ps, const vDouble& genes, vDouble& n_gen, Benchmarks * f){
  assert( ndim > 0 );
  assert( (ndim * ps) == genes.size() );
  assert( ps == size );

  std::uniform_real_distribution<double> random(0.00, 1.00);
  std::uniform_int_distribution<int> random_i(0, ps-1);

  int n1, n2, n3;
  double myCR, myF;

  uint i = 0, j = 0;
  for( ; i < ps; i++ ){

    /* Get three mutually different indexs */
    do n1 = random_i(rng); while (n1 == i);
    do n2 = random_i(rng); while (n2 == i || n2 == n1 );
    do n3 = random_i(rng); while (n3 == i || n3 == n1 || n3 == n2);
    //printf("%d %d %d\n", n1, n2, n3);

    myCR = CR[i];
    myF  = F[i];

    for( j = 0; j < ndim; j++ ){
      if( random(rng) <= myCR or ( j == ndim-1 ) ){
        n_gen[i * ndim + j] = genes[n1 * ndim + j] + myF * (genes[n2 * ndim + j] - genes[n3 * ndim + j] );

        if( n_gen[i * ndim + j] > f->getMaxX() )
          n_gen[i * ndim + j] = f->getMaxX();

        if( n_gen[i * ndim + j] < f->getMinX() )
          n_gen[i * ndim + j] = f->getMinX();

      } else
        n_gen[i * ndim + j] = genes[i * ndim + j];
    }
  }
}

void jDE::selection(
  uint ndim, uint ps,
  vDouble & genes, const vDouble & n_gen,
  vDouble& fitness, const vDouble & n_fitness
){
  assert( ndim > 0 );
  assert( ps > 0 );
  assert( genes.size() == n_gen.size() );
  assert( fitness.size() == n_fitness.size() );
  assert( genes.size() == (ndim * ps) );
  assert( fitness.size() ==  ps );

  for( uint i = 0; i < ps; i++ ){

    if( n_fitness[i] < fitness[i] ){
      for( uint j = 0; j < ndim; j++ )
        genes[i * ndim + j] = n_gen[ i * ndim + j];

      fitness[i] = n_fitness[i];
    }
  }
}

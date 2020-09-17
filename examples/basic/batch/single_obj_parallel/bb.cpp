#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

int main ( int argc , char ** argv ) {

  double f = 1e20, c1 = 1e20 , c2 = 1e20;
  double x[5];

  if ( argc >= 2 ) {
    c1 = 0.0 , c2 = 0.0;
    ifstream in ( argv[1] );
    for ( int i = 0 ; i < 5 ; i++ ) {
      in >> x[i];
      c1 += pow ( x[i]-1 , 2 );
      c2 += pow ( x[i]+1 , 2 );
    }
    f = x[4];
    if ( in.fail() )
      f = c1 = c2 = 1e20;
    else {
      c1 = c1 - 25;
      c2 = 25 - c2;
    }
    in.close();
  }
  cout << f << " " << c1 << " " << c2 << endl;
  return 0;
}

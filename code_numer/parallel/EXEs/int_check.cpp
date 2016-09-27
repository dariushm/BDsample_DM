#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>

using namespace std;
int main(int argc, char *argv[])
{

  double d1=30.0001;
  double d2=90000;
  int sw=10;
  double r1=d1/sw;
  double r2=d2/sw;
  double tol=1E-12;
  if(r1-int(r1)<tol)cout <<"integer! "<<r1<<endl;
  if(r2-int(r2)<tol)cout <<"integer! "<<r2<<endl;


} 

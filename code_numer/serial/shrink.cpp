/*
write out only some entries
maggie johnson
*/


#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>


using namespace std;


int main(int argc, char *argv[])
{
  ifstream infile(argv[1]);
  ofstream outfile(argv[2]);
  int t=0;
  double time;
  double ps, ps_ex;
  int write=50000;
  infile >>time >>ps>>ps_ex;
  outfile<<time<<'\t'<<ps<<'\t'<<ps_ex<<endl;
  
  while( !infile.eof() ){
    t++;
    infile >>time >>ps>>ps_ex;
    if(t%write==0){
      outfile<<time<<'\t'<<ps<<'\t'<<ps_ex<<endl;
      cout <<"time: "<<time<<endl;
    }
  }
  cout <<"total lines read: "<<t<<endl;
}


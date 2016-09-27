/*

calculate the solution to the diffusion equation with radiation
boundary conditions at r=sigma numerically, using
different integration schemes. explicit, and crank nicolson. First just solve the equation
that is spherically symmetric. 
maggie johnson
*/

#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include "mkl.h"
#include "mkl_cblas.h"


using namespace std;

void integrate_simple(int tbins, int nrbins, double delR, double deltat, double sigma, double ka, double Dtot, double BC, double **prt, double *passoc);
void integrate_simple_v2(int tbins, int nrbins, double delR, double deltat, double sigma, double ka, double Dtot, double BC, double **prt, double *passoc);
void crank_nicholson_out(int t1, int tbins, int nrbins, double delR, double deltat, double sigma, double ka, double Dtot, double BC, double *prt, double *passoc, double c, double &passoc1, ofstream &outfile);
void build_A_B_matrix(int N, double *A, double *B, double Dtot, double deltat, double c, double delR, double BC, double sigma);
void build_triD(int N, double *diag, double *sub, double *super, double Dtot, double deltat, double c, double delR, double BC, double sigma);
void build_Bband(int N, double *Bband, double Dtot, double deltat, double c, double delR, double BC, double sigma);
void build_triD_semi(int N, double *diag, double *sub, double *super, double Dtot, double deltat, double c, double delR, double BC, double sigma);

int main(int argc, char *argv[])
{
  int i,j;

  srand(static_cast<unsigned>(time(0)));

  /*1/D dp/dt=1/r^2*ddr(r^2 dpdr) = 2/r*dpdr+d2pdr2
    calculate change in time and space with delta function initial condition.
    p(r,0)=deltat(r-r0)/(4pir0^2)
  */
  
  /*first just use an explicit scheme (not Crank-Nicolson) because it will be
    easier to extend to 4 dimensions. Otherwise the matrix problem might become
    impossible.
  */
  double ka=10.0;
  double Dtot=10.0;
  int t;
   
  int nrbins=400;
  double maxR=20;
  double sigma=1.0;
  double delR=atof(argv[1]);//(maxR-sigma)/(1.0*nrbins);
  nrbins=int((maxR-sigma)/delR);
  //double *r=new double[nrbins];
  double r;
  int tbins=1000;
  double deltat=atof(argv[2]);//1E-4;//us
  
  double *prt=new double[nrbins];

  double *passoc=new double[tbins];
  passoc[0]=0.0;
  /*create initial conditions
    deltat function
  */
  
  double r0=atof(argv[3]);//6.05;
  double range=r0-sigma;
  cout <<"range: "<<range<<" index: "<<range/delR<<endl;
  double ir=range/delR;
  int index=int(round(ir));
  r0=index*delR+sigma;//so it's on a grid point
  cout <<"Delr: "<<delR<<" sigma: "<<sigma<<" r0: "<<r0<<" index: "<<index<<endl;
  for(i=0;i<nrbins;i++){
    r=i*delR+sigma;
    if(i==index)
      prt[i]=1.0/(4.0*M_PI*r*r*delR);
    else
      prt[i]=0.0;
  }
   /*Write out the distribution at each time step.*/
  char fname[300];
  sprintf(fname,"prt_crank_r0_%g_dt_%g.dat",r0, deltat);
  ofstream outfile;
  outfile.open(fname);
  outfile<<-0.0<<'\t';
  for(j=0;j<nrbins;j++){
    r=j*delR+sigma;
    outfile<<r<<'\t';
  }
  outfile<<endl;
  outfile<<0.0<<'\t';//time 0
  for(j=0;j<nrbins;j++){

    outfile<<prt[j]<<'\t';
  }
  outfile<<endl;
  /*Now iterate over time steps.*/
  double time;
  
  double BC=ka/(4.0*M_PI*Dtot*sigma*sigma);
  double c;
  double passoc1=0.0;
  c=0.5;//fully implicit if c=1, fully explicit if c=0.
  crank_nicholson_out(1, tbins,nrbins,  delR,  deltat,  sigma,  ka,  Dtot,  BC,  prt,  passoc,  c,  passoc1,  outfile);
  //c=0;
  //crank_nicholson_out(15, tbins, nrbins,  delR,  deltat,  sigma,  ka,  Dtot,  BC,  prt,  passoc,  c,  passoc1,  outfile);

  
  
  outfile.close();
  cout <<"finished. "<<endl;


       
}/*End main*/


void integrate_simple(int tbins, int nrbins, double delR, double deltat, double sigma, double ka, double Dtot, double BC, double **prt, double *passoc)
{
  int t, i;
  double time, r;
  double rhp, rhm;
  double passoc1=0;
  double psurvive=0;
  double psreal;
  double dr2=delR*delR;
  
  for(t=1;t<tbins;t++){
    time=t*deltat;
    r=sigma;
    rhp=r+delR*0.5;
    rhm=r-delR*0.5;
    i=0;
    //BC
    psurvive=0.0;
    prt[t][i]=prt[t-1][i]+Dtot*deltat/(r*r*dr2)*(rhp*rhp*prt[t-1][i+1]-2*r*r*prt[t-1][i]+rhm*rhm*(prt[t-1][i+1]-BC*2.0*delR*prt[t-1][i]));
    passoc1+=prt[t][0]*ka*deltat;
    psreal=1.0-passoc1;
    psurvive+=r*r*delR*prt[t][i]*4.0*M_PI;
    for(i=1;i<nrbins;i++){
      r=i*delR+sigma;
      rhp=r+delR*0.5;
      rhm=r-delR*0.5;
      prt[t][i]=prt[t-1][i]+Dtot*deltat/(r*r*dr2)*(rhp*rhp*prt[t-1][i+1]-2.0*r*r*prt[t-1][i]+rhm*rhm*prt[t-1][i-1]);
      psurvive+=r*r*delR*prt[t][i]*4.0*M_PI;
    }//done looping over space
    /*for each time step, integrate over all space to get survival prob*/
    passoc[t]=passoc1;
    for(i=0;i<nrbins;i++)
      prt[t][i]=prt[t][i]/psurvive*psreal;
  }//end looping over timesteps.
  
  
}

void integrate_simple_v2(int tbins, int nrbins, double delR, double deltat, double sigma, double ka, double Dtot, double BC, double **prt, double *passoc)
{
  int t, i;
  double time, r;
  double rhp, rhm;
  double passoc1=0;
  double psurvive=0;
  double psreal;
  double dr2=delR*delR;
  
  for(t=1;t<tbins;t++){
    time=t*deltat;
    r=sigma;
    rhp=r+delR*0.5;
    rhm=r-delR*0.5;
    i=0;
    //BC
    psurvive=0.0;
    prt[t][i]=prt[t-1][i]+Dtot*deltat/(r*delR)*(prt[t-1][i+1]-(prt[t-1][i+1]-BC*2.0*delR*prt[t-1][i]))+Dtot*deltat/dr2*(prt[t-1][i+1]-2.0*prt[t-1][i]+(prt[t-1][i+1]-2.0*delR*BC*prt[t-1][i]));
    passoc1+=prt[t][0]*ka*deltat;
    psreal=1.0-passoc1;
    psurvive+=r*r*delR*prt[t][i]*4.0*M_PI;
    for(i=1;i<nrbins;i++){
      r=i*delR+sigma;
      rhp=r+delR*0.5;
      rhm=r-delR*0.5;
      prt[t][i]=prt[t-1][i]+Dtot*deltat/(r*delR)*(prt[t-1][i+1]-prt[t-1][i-1])+Dtot*deltat/dr2*(prt[t-1][i+1]-2.0*prt[t-1][i]+prt[t-1][i-1]);
      psurvive+=r*r*delR*prt[t][i]*4.0*M_PI;
    }//done looping over space
    /*for each time step, integrate over all space to get survival prob*/
    passoc[t]=passoc1;
    for(i=0;i<nrbins;i++)
      prt[t][i]=prt[t][i]/psurvive*psreal;
  }//end looping over timesteps.
  
  
}
void crank_nicholson_out(int t1, int tbins, int nrbins, double delR, double deltat, double sigma, double ka, double Dtot, double BC, double *prt, double *passoc, double c, double &passoc1, ofstream &outfile)
{
  int t, i;
  double time, r;
  double rhp, rhm;
  
  double psurvive=0;
  double psreal;
  double dr2=delR*delR;
  int N=nrbins;
  /*Build LHS matrix*/
  //double c=0.5;//Crank Nicholson. If c=1.0, then fully implicit
  i=0;
  double *diag=new double[N];
  double *super=new double[N];
  double *sub=new double[N];
  
  build_triD( N,  diag,  sub,  super, Dtot, deltat, c,  delR, BC, sigma);
  

  /*build RHS B matrix, in full form and banded form*/
   double *B=new double[N*N];
   double *A=new double[N*N];
   build_A_B_matrix(N, A, B, Dtot, deltat, c, delR, BC, sigma);
    
  /*Build compact banded B matrix*/
  double *Bband=new double[N*3];
  build_Bband( N,  Bband,  Dtot, deltat,  c,  delR, BC,  sigma);
  double *bcol=new double[N];
  
  int nrhs=1;
  int info;
  char trans='N';
  int kl=1;
  double alpha=1;
  double beta=0;
  int lda=3;//or N?
  int incx=1;
  /*Would it be faster to invert the LHS matrix (that doesn't change), 
    Au(t+1)=Bu(t), so you end up with u(t+1)=(A^-1*B)u(t), where 
    the matrix M=(A^-1*B) is fixed. 
  */

  int *ipiv=new int[N];
//   ofstream af("Amat.dat");
//   ofstream bf("Bmat.dat");
//   int j;
//   for(i=0;i<N;i++){
//     for(j=0;j<N;j++){
//       af<<A[j*N+i]<<'\t';
//       bf<<B[j*N+i]<<'\t';
//     }
//     bf<<endl;
//     af<<endl;
//   }
  ofstream rhs("Rhs.dat");
  //passoc1=prt[0]*ka*deltat; 
  for(t=t1;t<tbins;t++){
    time=t*deltat;
    psurvive=0.0;
    cblas_dgbmv( CblasColMajor, CblasNoTrans, N, N, kl, kl,alpha, Bband, lda, prt, incx, beta, bcol, incx);
    
    //cblas_dgemv( CblasColMajor, CblasNoTrans, N, N, alpha, B, N, prt, incx, beta, bcol, incx);
    for(i=0;i<N;i++)
      rhs<<bcol[i]<<'\t';
    rhs<<endl;

    dgtsv(&N, &nrhs, sub, diag, super, bcol, &N, &info);
    //dgesv(&N, &nrhs, A, &N, ipiv, bcol2, &N, &info);
    cout <<"Solved linear equations "<<endl;
    for(i=0;i<N;i++)
      prt[i]=bcol[i];//new solution
    
    passoc1+=prt[0]*ka*deltat;
    psreal=1.0-passoc1;
    for(i=0;i<nrbins;i++){
      r=i*delR+sigma;
      psurvive+=r*r*delR*prt[i]*4.0*M_PI;
    }
    
    passoc[t]=passoc1;
    outfile<<passoc1<<'\t';
    for(i=0;i<nrbins;i++){
      prt[i]=prt[i]/psurvive*psreal;
      outfile<<prt[i]<<'\t';
    }
    outfile<<endl;
    build_triD( N,  diag,  sub,  super, Dtot, deltat, c,  delR, BC, sigma);
    build_Bband( N,  Bband,  Dtot, deltat,  c,  delR, BC,  sigma);
    //build_A_B_matrix(N, A, B, Dtot, deltat, c, delR, BC, sigma);
  }//end looping over timesteps.
  
  
}
void build_triD(int N, double *diag, double *sub, double *super, double Dtot, double deltat, double c, double delR, double BC, double sigma)
{
  int i;
  double r;
  double dr2=delR*delR;
  double Dtc=Dtot*deltat*c;
  i=0;
  r=sigma;
  diag[i]=1.0+2.0*Dtc/dr2-2.0*delR*BC*(Dtc/(r*delR)-Dtc/dr2);
  super[i]=-2.0*Dtc/dr2;
  //diag[i]=1.0+2.0*Dtc/dr2;
  //super[i]=-(Dtc/(r*delR)+Dtc/dr2);
    
  for(i=1;i<N-1;i++){
    r=i*delR+sigma;
    diag[i]=1.0+2.0*Dtc/dr2;
    super[i]=-(Dtc/(r*delR)+Dtc/dr2);
    sub[i-1]=Dtc/(r*delR)-Dtc/dr2;
  }
  i=N-1;
  r=i*delR+sigma;
  diag[i]=1.0+2.0*Dtc/dr2;
  sub[i-1]=Dtc/(r*delR)-Dtc/dr2;
}
void build_triD_semi(int N, double *diag, double *sub, double *super, double Dtot, double deltat, double c, double delR, double BC, double sigma)
{
  int i;
  double r;
  double dr2=delR*delR;
  double Dtc=Dtot*deltat*c;
  i=0;
  r=sigma;
  /*Make future boundary point fully explicit*/
  diag[i]=1.0;//+2.0*Dtc/dr2-2.0*delR*BC*(Dtc/(r*delR)-Dtc/dr2);
  super[i]=0;//-2.0*Dtc/dr2;
  i=1;
  r=i*delR+sigma;
  diag[i]=1.0;//+2.0*Dtc/dr2;
  super[i]=0;//-(Dtc/(r*delR)+Dtc/dr2);
  sub[i-1]=0;//Dtc/(r*delR)-Dtc/dr2;
  int start=1;//if 2, first two rows diag, if 1, only first row is diagonal(overwrite above i=1).
  for(i=start;i<N-1;i++){
    r=i*delR+sigma;
    diag[i]=1.0+2.0*Dtc/dr2;
    super[i]=-(Dtc/(r*delR)+Dtc/dr2);
    sub[i-1]=Dtc/(r*delR)-Dtc/dr2;
  }
  i=N-1;
  r=i*delR+sigma;
  diag[i]=1.0+2.0*Dtc/dr2;
  sub[i-1]=Dtc/(r*delR)-Dtc/dr2;
}
void build_Bband(int N, double *Bband, double Dtot, double deltat, double c, double delR, double BC, double sigma)
{
  //col=0
  //Bband[0*3+0]=0 =>outside matrix, col0, row0
//   Bband[0*3+1]=B[0];//col0, row1
//   Bband[0*3+2]=B[1];//col0, row2
//   for(i=1;i<N-1;i++){
//     //i is the column
//     Bband[i*3+0]=B[i*N+(i-1)];
//     Bband[i*3+1]=B[i*N+i];
//     Bband[i*3+2]=B[i*N+(i+1)];
//   }
//   i=N-1;
//   Bband[i*3+0]=B[i*N+(i-1)];
//   Bband[i*3+1]=B[i*N+i];
  //last entry is outside matrix on the bottom.

  //col=0
  //Bband[0*3+0]=0 =>outside matrix, col0, row0
  int i;
  double r;
  double dr2=delR*delR;
  double Dtmc=Dtot*deltat*(1.0-c);
  r=sigma;//i=0
  Bband[0*3+1]=1.0-Dtmc*2.0/dr2-2.0*delR*BC*(-Dtmc/(r*delR)+Dtmc/dr2);//diag;//row0, c1
  r=delR+sigma;//i=1
  Bband[0*3+2]=-Dtmc/(r*delR)+Dtmc/dr2;//sub//col0, row2

  for(i=1;i<N-1;i++){
    //i is the column
    r=(i-1)*delR+sigma;
    Bband[i*3+0]=Dtmc/(r*delR)+Dtmc/dr2;//super;
    r=i*delR+sigma;
    Bband[i*3+1]=1.0-Dtmc*2.0/dr2;//diag;
    r=(i+1)*delR+sigma;
    Bband[i*3+2]=-Dtmc/(r*delR)+Dtmc/dr2;//sub
  }
  Bband[1*3+0]=2.0*Dtmc/dr2;//row0, c2, other boundary node at sigma
  i=N-1;
  r=(i-1)*delR+sigma;
  Bband[i*3+0]=Dtmc/(r*delR)+Dtmc/dr2;//super;
  r=i*delR+sigma;
  Bband[i*3+1]=1.0-Dtmc*2.0/dr2;//diag
  //last entry is outside matrix on the bottom.


}
void build_A_B_matrix(int N, double *A, double *B, double Dtot, double deltat, double c, double delR, double BC, double sigma)
{

  int i;
  for(i=0;i<N*N;i++){
    B[i]=0.0;
    A[i]=0.0;
  }
  double r;
  double dr2=delR*delR;
  
  r=sigma;
  double Dtc=Dtot*deltat*c;
  double Dtmc=Dtot*deltat*(1.0-c);
  B[0]=1.0-Dtmc*2.0/dr2-2.0*delR*BC*(-Dtmc/(r*delR)+Dtmc/dr2);//diag;//row0, c1
  B[1*N]=2.0*Dtmc/dr2;//row0, c2
  A[0]=1.0+Dtc*2.0/dr2-2.0*delR*BC*(Dtc/(r*delR)-Dtc/dr2);//diag;//row0, c1
  A[1*N]=-2.0*Dtc/dr2;//row0, c2
  for(i=1;i<N-1;i++){
    r=i*delR+sigma;
    B[(i-1)*N+i]=-Dtmc/(r*delR)+Dtmc/dr2;//sub
    B[i*N+i]=1.0-Dtmc*2.0/dr2;//diag
    B[(i+1)*N+i]=Dtmc/(r*delR)+Dtmc/dr2;//super
    A[(i-1)*N+i]=Dtc/(r*delR)-Dtc/dr2;//sub
    A[i*N+i]=1.0+Dtc*2.0/dr2;//diag
    A[(i+1)*N+i]=-Dtc/(r*delR)-Dtc/dr2;//super
  }
  i=N-1;
  r=i*delR+sigma;
  B[i*N+i]=1.0-Dtmc*2.0/dr2;//diag
  B[(i-1)*N+i]=-Dtmc/(r*delR)+Dtmc/dr2;//sub
  A[i*N+i]=1.0+Dtc*2.0/dr2;//diag
  A[(i-1)*N+i]=Dtc/(r*delR)-Dtc/dr2;//sub
}

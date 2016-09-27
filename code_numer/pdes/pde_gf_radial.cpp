/*

calculate the solution to the diffusion equation with radiation
boundary conditions at r=sigma numerically, using
explicit integrator. First just solve the equation
that is spherically symmetric. 

v2, which discretizes 1/r dp/dr + d2p/dr2, is much
more accurate than v1, which discretizes 1/r^2(r^2dp/dr).
Both versions need to be renormalized to keep probability distibution
summing to survival probability.

Survival probability is calculated as 1-passoc, where passoc(t)=int[dt k*p(sigma, t)].
has to be integrated over all time starting from when p first hits the boundary

Trying to integrate over distribution is hopelessly inaccurate.
maggie johnson
*/

#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>



using namespace std;

void integrate_simple(int tbins, int nrbins, double delR, double deltat, double sigma, double ka, double Dtot, double BC, double **prt, double *passoc);
void integrate_simple_v2(int tbins, int nrbins, double delR, double deltat, double sigma, double ka, double Dtot, double BC, double **prt, double *passoc);



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
  double maxR=15;
  double sigma=1.0;
  double delR=atof(argv[1]);//(maxR-sigma)/(1.0*nrbins);
  nrbins=int((maxR-sigma)/delR);
  //double *r=new double[nrbins];
  double r;
  int tbins=1000;
  double deltat=atof(argv[2]);//1E-4;//us
  
  double **prt=new double*[tbins];
  for(i=0;i<tbins;i++)
    prt[i]=new double[nrbins];
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
      prt[0][i]=1.0/(4.0*M_PI*r*r*delR);
    else
      prt[0][i]=0.0;
  }
 
  /*Now iterate over time steps.*/
  double time;
  
  double BC=ka/(4.0*M_PI*Dtot*sigma*sigma);
  
  integrate_simple_v2(tbins, nrbins, delR, deltat, sigma, ka, Dtot, BC, prt, passoc);
  
  
  /*Write out the distribution at each time step.*/
  char fname[300];
  sprintf(fname,"prt_v2numer_r0_%g_dt_%g.dat",r0, deltat);
  ofstream outfile;
  outfile.open(fname);
  outfile<<-0.0<<'\t';
  for(t=0;t<tbins;t++)
    outfile<<passoc[t]<<'\t';
  outfile<<endl;
  /*now print probability distiburion given r0 and time */
  for(j=0;j<nrbins;j++){
    r=j*delR+sigma;
    outfile<<r<<'\t';
    for(t=0;t<tbins;t++){
      
      outfile<<prt[t][j]<<'\t';
    }
    outfile<<endl;
  }
  
  
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

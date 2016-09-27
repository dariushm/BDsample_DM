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


void integrate_simple_v2(int tbins, int nrbins, double delR, double deltat, double sigma, double ka, double Dtot, double BC, double **prt, double *passoc);
void integrate_spherical_v2(int tbins, int nrbins, double delR, double deltat, double sigma, double ka, double Dtot, double BC, double ***prt, double ***prt2, double *passoc, int polbins, double delpol, int azbins, double delaz, double r0);
//void integrate_spherical_capgrid(int tbins, int nrbins, double delR, double deltat, double sigma, double ka, double Dtot, double BC, double ***prt, double ***prt2, double *passoc, int polbins, double delpol, int azbins, double delaz, double r0, double mindel);
void integrate_spherical_capgrid(int tbins, int nrbins, double delR, double deltat, double sigma, double ka, double Dtot, double BC, double ***prt, double ***prt2, double *passoc, int polbins, double delpol, int azbins, double delaz, double r0, double mindel, double ***pfree, double tstart, double az0, double pol0);
void fullk_bc_npole(int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2);
void fullk_bc_spole(int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2);
void fullk_bc(int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2);
void fullk(int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2);
void fullk_npole(int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2);
void fullk_spole(int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2);
void gridk_bc_npole(int nbins, int nsmooth, int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2);
void gridk_bc_spole(int nbins, int nsmooth, int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2);
void gridk_bc(int nbins, int nsmooth, int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2);
void gridk_npole(int nbins, int nsmooth, int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2);
void gridk_spole(int nbins, int nsmooth, int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2);
void gridk(int nbins, int nsmooth, int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2);
void free_prop_sphere(double sigma, int nrbins, int polbins, int azbins, double delR, double delpol, double delaz, double Dtot, double time, double r0, double az0, double pol0, double ***pfree, double &norm);
void free_prop_cartesian(double sigma, int xbins, int ybins, int zbins, double delx, double dely, double delz, double Dtot, double time, double r0, double az0, double pol0, double ***pfree);

int main(int argc, char *argv[])
{
  int i,j, k;

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
  double maxR=10;
  double sigma=1.0;
  double delR=atof(argv[1]);//(maxR-sigma)/(1.0*nrbins);
  nrbins=int((maxR-sigma)/delR);
  //double *r=new double[nrbins];
  double r;
  int tbins=100;
  double deltat=atof(argv[2]);//1E-4;//us
  int polbins=50;
  int azbins=100;
  double delpol=M_PI/(1.0*polbins);
  double delaz=2.0*M_PI/(1.0*azbins);
  double ***prt=new double**[nrbins+1];
  double ***pfree=new double**[nrbins+1];
  double ***prt2=new double**[nrbins+1];
  for(i=0;i<nrbins+1;i++){
    prt[i]=new double*[polbins];
    pfree[i]=new double*[polbins];
    prt2[i]=new double*[polbins];
  }
  for(i=0;i<nrbins+1;i++){
    for(j=0;j<polbins;j++){
      prt[i][j]=new double[azbins];
      pfree[i][j]=new double[azbins];
      prt2[i][j]=new double[azbins];
    }
  }
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
  //  r0=index*delR+sigma;//so it's on a grid point
  double pol0=M_PI/2.0;
  int pindex=int(pol0/delpol);
  double az0=0;
  int aindex=int(az0/delaz);
  cout <<"Delr: "<<delR<<" sigma: "<<sigma<<" r0: "<<r0<<" index: "<<index<<endl;
  cout <<"Polar: "<<pol0<<" index; "<<pindex<<" az: "<<az0<<" index: "<<aindex<<endl;
  /*Deltat function in r, or in r, polar and az angles
    for the later, association rate should be the same, but wavefronts will be different.
  */
  double psurvive=0;
  double theta, st;
  double time=atof(argv[4]);
  double psurvivef1=0;
  free_prop_sphere(sigma, nrbins, polbins, azbins, delR, delpol, delaz, Dtot,  time,  r0,  az0,pol0, pfree, psurvivef1);  
  ofstream ffile;
  char cfname[300];
  sprintf(cfname, "free_sol_%g_%g.dat",r0, time);
  ffile.open(cfname);
  double tstart=time;
  for(i=0;i<nrbins+1;i++){
    r=i*delR+sigma;
    for(j=0;j<polbins;j++){
      theta=(j+0.5)*delpol;
      st=sin(theta);
      for(k=0;k<azbins;k++){
	ffile<<pfree[i][j][k]/psurvivef1<<'\t';
	//psurvive+=pfree[i][j][k]*delaz*delpol*st*r*r*delR;
      }
    }
    ffile<<endl;
  }
  cout <<"Total pfree 1?: "<<psurvivef1<<endl;

  
  
  psurvive=0;
  for(i=0;i<nrbins+1;i++){
    r=i*delR+sigma;
    for(j=0;j<polbins;j++){
      theta=(j+0.5)*delpol;
      st=sin(theta);
       
      for(k=0;k<azbins;k++){
	if(i==index){
	  if(j>45 &&j<55){
	    // if(k>5&&k<95)
	      prt[i][j][k]=2.0/(r*r*delR*4.0*M_PI);
	      //else
	      //prt[i][j][k]=0.0;
	  }else
	    prt[i][j][k]=0.0;
	}else{
	  prt[i][j][k]=0.0;
	}
	
	psurvive+=prt[i][j][k]*r*r*delR*st*delpol*delaz;
	// if(i==index && j==pindex &&k==aindex)
//  	  prt[i][j][k]=1.0/(r*r*delR*sin(pol0)*delpol*delaz);
// 	else
//  	  prt[i][j][k]=0.0;
	
      }
    }
  }
  for(i=0;i<nrbins+1;i++){
    for(j=0;j<polbins;j++){
      for(k=0;k<azbins;k++){
	prt[i][j][k]=pfree[i][j][k]/psurvivef1;
      }
    }
  }
  
  
  int xbins=100;
  int ybins=100;
  int zbins=100;
  double ***pcart=new double**[xbins];
  for(i=0;i<xbins;i++){
    pcart[i]=new double*[ybins];
  }
  for(i=0;i<xbins;i++){
    for(j=0;j<ybins;j++){
      pcart[i][j]=new double[zbins];
    }
  }
  double delx=0.01;
  double dely=0.01;
  double delz=0.01;
  free_prop_cartesian(sigma, xbins, ybins, zbins, delx, dely, delz, Dtot,  time,  r0,  az0,pol0, pcart);  
  
  /*Write out the distribution at each time step.*/
  
  /*Now iterate over time steps.*/

  
  double BC=ka/(4.0*M_PI*Dtot*sigma*sigma);
  ofstream tfile("polars.dat");
  for(j=0;j<polbins;j++)
    tfile<<(j+0.5)*delpol<<endl;
  ofstream sfile("azs.dat");
  for(j=0;j<azbins;j++)
    sfile<<j*delaz<<endl;
  ofstream rfile("rads.dat");
  for(j=0;j<nrbins;j++)
    rfile<<j*delR+sigma<<endl;
  
  double mindel=delpol;
  //integrate_spherical_capgrid(tbins, nrbins, delR, deltat, sigma, ka, Dtot, BC, prt, prt2, passoc, polbins, delpol, azbins, delaz, r0, mindel);
  integrate_spherical_capgrid(tbins, nrbins, delR, deltat, sigma, ka, Dtot, BC, prt, prt2, passoc, polbins, delpol, azbins, delaz, r0, mindel, pfree, tstart, az0, pol0);
  
  double **psimp=new double*[tbins];
  for(i=0;i<tbins;i++)
    psimp[i]=new double[nrbins+1];
  for(i=0;i<nrbins+1;i++){
    r=i*delR+sigma;
    if(i==index)
      psimp[0][i]=1.0/(4.0*M_PI*r*r*delR);
    else
      psimp[0][i]=0.0;
  }
  double *passoc_simp=new double[tbins];
  integrate_simple_v2( tbins,  nrbins,  delR,  deltat, sigma, ka,  Dtot,  BC, psimp, passoc_simp);
  char fname[300];
  sprintf(fname,"passoc_r0_%g_deltat_%g.dat", r0, deltat);
  ofstream pfile;
  pfile.open(fname);
  for(t=1;t<tbins;t++)
    pfile <<t*deltat<<' '<<passoc[t]<<' '<<passoc_simp[t]<<endl;

  
  
  
  cout <<"finished. "<<endl;


       
}/*End main*/
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

void integrate_spherical_v2(int tbins, int nrbins, double delR, double deltat, double sigma, double ka, double Dtot, double BC, double ***prt, double ***prt2, double *passoc, int polbins, double delpol, int azbins, double delaz, double r0)
{
  int t, i;
  double time, r;
  double passoc1=0;
  double psurvive=0;
  double psreal;
  double st;
  double theta;
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  int j, k;
  double dr2=delR*delR;
  double r2;
  char fname[300];
  ofstream outfile;
  int twrite=5;
  ofstream ffile;
  
  for(t=1;t<tbins;t++){
    time=t*deltat;
    if(t%twrite==0){
      
      sprintf(fname,"pspheret_v2numer_r0_%g_dt_%g.dat",r0, time);
      cout <<"new file: "<<fname<<" time: "<<time<<endl;
      outfile.open(fname);
      
    }
    psurvive=0.0;
    r=sigma;
    r2=r*r;
    i=0;//BC, r=sigma.
    j=0;//BC, theta=0+delta
    theta=(j+0.5)*delpol;//need to skip the poles.
    st=sin(theta);
    k=0;//BC, az=0.
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-2.0*delR*BC*prt[i][j][k]))+Dtot*deltat/r2*( cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j+1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j+1][k])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][azbins-1]) );//for i=0, j=0, k=0
    passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    
    for(k=1;k<azbins-1;k++){
      prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-2.0*delR*BC*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j+1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j+1][k])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][k-1]) );//for i=0, j=0, k interior
      passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
      psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    }
    k=azbins-1;//BC, az=2pi-delta.
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-2.0*delR*BC*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j+1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j+1][k])+1.0/(st*st*delaz2)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][k-1]) );//for i=0, j=0, k=end
    passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    
    j=polbins-1;//BC, theta=pi-delta
    theta=(j+0.5)*delpol;//need to skip the poles.
    st=sin(theta);
    k=0;//BC, az=0.
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-2.0*delR*BC*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j-1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j-1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][azbins-1]) );//for i=0, j=end, k=0
    passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    
    for(k=1;k<azbins-1;k++){
      prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-2.0*delR*BC*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j-1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j-1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][k-1]) );//for i=0, j=0, k interior
      passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
      psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    }
    k=azbins-1;//BC, az=2pi-delta.
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-2.0*delR*BC*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j-1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j-1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][k-1]) );//for i=0, j=end, k=end
    passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    
    k=0;
    for(j=1;j<polbins-1;j++){
      theta=(j+0.5)*delpol;
      st=sin(theta);
      prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-2.0*delR*BC*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][azbins-1]) );//for i=0, j=interior, k=0
      passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
      psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    }
    k=azbins-1;
    for(j=1;j<polbins-1;j++){
      theta=(j+0.5)*delpol;
      st=sin(theta);
      prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-2.0*delR*BC*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][k-1]) );//for i=0, j=interior, k=end
      passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
      psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    }
    
    /*below, all interior points at r=sigma*/
    for(j=1;j<polbins-1;j++){
      theta=(j+0.5)*delpol;//need to skip the poles.
      st=sin(theta);
      for(k=1;k<azbins-1;k++){
	//azangle=(k+0.5)*delaz;
	prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-2.0*delR*BC*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][k-1]) );
	passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
	psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
      }
    }	
    //passoc1/=(4.0*M_PI);
    psreal=1.0-passoc1/(4.0*M_PI);
    passoc[t]=passoc1/(4.0*M_PI);
    /************R interior**************/
    /*interior r points, at end i+1, p=0 (r->inf)*/
    for(i=1;i<nrbins;i++){
      r=i*delR+sigma;
      r2=r*r;
      /*polar and az interior*/
      for(j=1;j<polbins-1;j++){
	theta=(j+0.5)*delpol;//need to skip the poles.
	st=sin(theta);
	for(k=1;k<azbins-1;k++){
	  //azangle=(k+0.5)*delaz;
	  

	  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][k-1]) );
	  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
	}//done looping over az
      }//theta
      
      /*For interior r, also do boundaries on angles, 8 additional regions*/
      j=0;//BC, theta=0+delta
      theta=(j+0.5)*delpol;//need to skip the poles.
      st=sin(theta);
      k=0;//BC, az=0.
      prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j+1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j+1][k])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][azbins-1]) );//for i=inter, j=0, k=0
      
      psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
      
      for(k=1;k<azbins-1;k++){
	prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j+1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j+1][k])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][k-1]) );//for i=inter, j=0, k interior
	
	psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
      }
      k=azbins-1;//BC, az=2pi-delta.
      prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j+1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j+1][k])+1.0/(st*st*delaz2)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][k-1]) );//for i=inter, j=0, k=end
      
      psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
      
      j=polbins-1;//BC, theta=pi-delta
      theta=(j+0.5)*delpol;//need to skip the poles.
      st=sin(theta);
      k=0;//BC, az=0.
      prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j-1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j-1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][azbins-1]) );//for i=inter, j=end, k=0
      
      psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
      
      for(k=1;k<azbins-1;k++){
	prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j-1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j-1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][k-1]) );//for i=inter, j=end, k interior
	
	psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
      }
      k=azbins-1;//BC, az=2pi-delta.
      prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j-1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j-1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][k-1]) );//for i=inter, j=end, k=end
      
      psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
      
      k=0;
      for(j=1;j<polbins-1;j++){
	theta=(j+0.5)*delpol;
	st=sin(theta);
	prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][azbins-1]) );//for i=inter, j=interior, k=0
	
	psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
      }
      k=azbins-1;
      for(j=1;j<polbins-1;j++){
	theta=(j+0.5)*delpol;
	st=sin(theta);
	prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][k-1]) );//for i=inter, j=interior, k=end
	
	psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
      }
      
      
    }//space
    cout <<"time: "<<t<<" dt: "<<t*deltat<<" normalization: "<<psurvive<<" passoc: "<<passoc1<<" 1-passoc: "<<psreal<<endl;
    /*for each time step, integrate over all space to get survival prob*/
    for(i=0;i<nrbins;i++){
      for(j=0;j<polbins;j++){
	for(k=0;k<azbins;k++){
	  prt[i][j][k]=prt2[i][j][k]/psurvive*psreal;
	}
      }
    }
    
    if(t%twrite==0){
      /*for each time step, integrate over all space to get survival prob*/
      for(i=0;i<nrbins;i++){
	for(j=0;j<polbins;j++){
	  for(k=0;k<azbins;k++){
	    outfile<<prt[i][j][k]<<'\t';
	  }
	}
	outfile<<endl;
      }
      outfile.close();
      //solve free
    }
  }//end looping time steps.
    
}

//void integrate_spherical_capgrid(int tbins, int nrbins, double delR, double deltat, double sigma, double ka, double Dtot, double BC, double ***prt, double ***prt2, double *passoc, int polbins, double delpol, int azbins, double delaz, double r0, double mindel)
void integrate_spherical_capgrid(int tbins, int nrbins, double delR, double deltat, double sigma, double ka, double Dtot, double BC, double ***prt, double ***prt2, double *passoc, int polbins, double delpol, int azbins, double delaz, double r0, double mindel, double ***pfree, double tstart, double az0, double pol0)
{
  /*Change evaluations at grid points at poles, by only updating some grid points and 
    interpolating between them for the rest. 
  */
  int t, i;
  double time, r;
  double passoc1=0;
  double psurvive=0;
  double psreal;
  double st;
  double theta;
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  int j, k;
  double dr2=delR*delR;
  double r2;
  char fname[300];
  ofstream outfile;
  ofstream ffile;
  int twrite=1;
  int nsmooth;
  double arc;
  int gbins;
  double norm;
  for(t=1;t<tbins;t++){
    time=t*deltat+tstart;
    if(t%twrite==0){

      sprintf(fname,"pspheret_v2numer_r0_%g_dt_%g.dat",r0, time);
      cout <<"new file: "<<fname<<" time: "<<time<<endl;
      outfile.open(fname);
      sprintf(fname,"pfree_v2numer_r0_%g_dt_%g.dat",r0, time);
      cout <<"new file: "<<fname<<" time: "<<time<<endl;
      ffile.open(fname);

    }
    /*First do r=sigma with BC*/
    psurvive=0.0;
    r=sigma;
    r2=r*r;
    i=0;//BC, r=sigma.
    j=0;//BC, theta=0+delta
    theta=(j+0.5)*delpol;//need to skip the poles.
    st=sin(theta);

    //for each grid point, decide which azimuthals to include 
    arc=r*st*delaz;
    if(arc<mindel){
      //choose a subset of grid points in k to update, interpolate the rest
      nsmooth=int(ceil(mindel/arc));
      
      gbins=azbins/nsmooth;
      if(gbins%2==1){
	gbins-=1;//even number of bins

      }
      nsmooth=int(azbins/gbins);
      //cout <<" at r=: "<<r<<" theta=: "<<theta<<" too small arc: "<<arc<<" nsmooth: "<<nsmooth<<endl;
      //so every nsmooth steps only update the az
      //which bins to use? just start at zero and loop around.

      gridk_bc_npole(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      
      
    }else{
      fullk_bc_npole( i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);
    }

    for(j=1;j<polbins-1;j++){
      theta=(j+0.5)*delpol;//need to skip the poles.
      st=sin(theta);
      
      //for each grid point, decide which azimuthals to include 
      arc=r*st*delaz;
      if(arc<mindel){
	//choose a subset of grid points in k to update, interpolate the rest
	nsmooth=int(ceil(mindel/arc));
	gbins=azbins/nsmooth;
	if(gbins%2==1){
	  gbins-=1;//even number of bins
	  
	}
	nsmooth=int(azbins/gbins);
	//cout <<" at r=: "<<r<<" theta=: "<<theta<<" too small arc: "<<arc<<" nsmooth: "<<nsmooth<<endl;
	//so every nsmooth steps only update the az
	//which bins to use? just start at zero and loop around.
	gridk_bc(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      
	
      }else{
	fullk_bc( i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);
      }
    }
    j=polbins-1;
    theta=(j+0.5)*delpol;//need to skip the poles.
    st=sin(theta);
    
    //for each grid point, decide which azimuthals to include 
    arc=r*st*delaz;
    if(arc<mindel){
      //choose a subset of grid points in k to update, interpolate the rest
      nsmooth=int(ceil(mindel/arc));
      
      gbins=azbins/nsmooth;
      if(gbins%2==1){
	gbins-=1;//even number of bins
	
      }
      nsmooth=int(azbins/gbins);
      //cout <<" at r=: "<<r<<" theta=: "<<theta<<" too small arc: "<<arc<<" nsmooth: "<<nsmooth<<endl;
      //so every nsmooth steps only update the az
      //which bins to use? just start at zero and loop around.
      gridk_bc_spole(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      
      
    }else{
      fullk_bc_spole( i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);
    }
    /*End over all angles at the r=sigma boundary*/
    psreal=1.0-passoc1/(4.0*M_PI);
    passoc[t]=passoc1/(4.0*M_PI);

    
    /*************All r, non BC****************/
    for(i=1;i<nrbins;i++){
      r=i*delR+sigma;
      r2=r*r;
      j=0;
      theta=(j+0.5)*delpol;//need to skip the poles.
      st=sin(theta);
      //for each grid point, decide which azimuthals to include 
      arc=r*st*delaz;
      if(arc<mindel){
	//choose a subset of grid points in k to update, interpolate the rest
	nsmooth=int(ceil(mindel/arc));
	gbins=azbins/nsmooth;
	if(gbins%2==1){
	  gbins-=1;//even number of bins
	}
	nsmooth=int(azbins/gbins);
	//cout <<" at r=: "<<r<<" theta=: "<<theta<<" too small arc: "<<arc<<" nsmooth: "<<nsmooth<<endl;
	//so every nsmooth steps only update the az
	//which bins to use? just start at zero and loop around.
	
	gridk_npole(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      
      }else{
	fullk_npole( i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);
      }

      for(j=1;j<polbins-1;j++){
	theta=(j+0.5)*delpol;//need to skip the poles.
	st=sin(theta);
	//for each grid point, decide which azimuthals to include 
	arc=r*st*delaz;
	if(arc<mindel){
	  //choose a subset of grid points in k to update, interpolate the rest
	  
	  nsmooth=int(ceil(mindel/arc));
	  gbins=azbins/nsmooth;
	  if(gbins%2==1){
	    gbins-=1;//even number of bins
	  }
	  nsmooth=int(azbins/gbins);
	  //cout <<" at r=: "<<r<<" theta=: "<<theta<<" too small arc: "<<arc<<" nsmooth: "<<nsmooth<<endl;
	  //so every nsmooth steps only update the az
	  //which bins to use? just start at zero and loop around.
	  gridk(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      

	  
	}else{
	  fullk( i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);
	}
	
      }
      j=polbins-1;
      theta=(j+0.5)*delpol;//need to skip the poles.
      st=sin(theta);
      //for each grid point, decide which azimuthals to include 
      arc=r*st*delaz;
      if(arc<mindel){
	//choose a subset of grid points in k to update, interpolate the rest
	nsmooth=int(ceil(mindel/arc));
	gbins=azbins/nsmooth;
	if(gbins%2==1){
	  gbins-=1;//even number of bins
	}
	nsmooth=int(azbins/gbins);
	//cout <<" at r=: "<<r<<" theta=: "<<theta<<" too small arc: "<<arc<<" nsmooth: "<<nsmooth<<endl;
	//so every nsmooth steps only update the az
	//which bins to use? just start at zero and loop around.
	gridk_spole(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      
	
      }else{
	fullk_spole( i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);
      }
    }//end R
    /*renormalize*/
    cout <<"time: "<<t<<" dt: "<<time<<" normalization: "<<psurvive<<" passoc: "<<passoc[t]<<" 1-passoc: "<<psreal<<endl;
    /*for each time step, integrate over all space to get survival prob*/
    for(i=0;i<nrbins;i++){
      for(j=0;j<polbins;j++){
	for(k=0;k<azbins;k++){
	  prt[i][j][k]=prt2[i][j][k]/psurvive*psreal;
	}
      }
    }
    
    if(t%twrite==0){
      free_prop_sphere(sigma, nrbins, polbins, azbins, delR, delpol, delaz, Dtot,  time,  r0,  az0,pol0, pfree, norm);  
      cout <<"free prop normalization at time: "<<time<<" norm: "<<norm<<endl;
      /*for each time step, integrate over all space to get survival prob*/
      for(i=0;i<nrbins;i++){
	for(j=0;j<polbins;j++){
	  for(k=0;k<azbins;k++){
	    outfile<<prt[i][j][k]<<'\t';
	    ffile<<pfree[i][j][k]/norm<<'\t';
	  }
	}
	outfile<<endl;
	ffile<<endl;
      }
      outfile.close();
      ffile.close();
      
    }
  }
}


void fullk_bc_npole(int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2)
{
  /*For npole, instead of j-1, use j, kcross (the azimuth angle across the pole, same polar angle)*/
  int k=0;//BC, az=0.
  int kcross=k+azbins/2;//for k
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-2.0*delR*BC*prt[i][j][k]))+Dtot*deltat/r2*( cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j][kcross]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j][kcross])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][azbins-1]) );//for i=0, j=0, k=0
  passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  
  for(k=1;k<azbins-1;k++){
    if(k<azbins/2)
      kcross=k+azbins/2;//for k
    else
      kcross=k-(azbins/2);//for k
    
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-2.0*delR*BC*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j][kcross]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j][kcross])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][k-1]) );//for i=0, j=0, k interior
    passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  k=azbins-1;//BC, az=2pi-delta.
  kcross=k-azbins/2;
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-2.0*delR*BC*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j][kcross]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j][kcross])+1.0/(st*st*delaz2)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][k-1]) );//for i=0, j=0, k=end
  passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;

}  
void gridk_bc_npole(int nbins, int nsmooth, int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2)
{
  int k=0;//BC, az=0.
  int ncross=0+nbins/2;
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;

  int bleft=0;
  int eleft=0;
  int ex1=0;
  int ex2=0;
  if(azbins%nbins!=0){
    bleft=azbins-(nbins-1)*nsmooth;
    eleft=bleft-nsmooth;
    ex1=eleft/2;
    ex2=eleft-ex1;
  }
  //cout <<"In gridk BC Npole: nbins, "<<nbins<<" nsmooth: "<<nsmooth<<" last bin: "<<bleft<<" leftover: "<<eleft<<" ex1, ex2: "<<ex1<<' '<<ex2<<endl;
  double dpda;
  /*First half of azimuths*/
  int n=0;
  int nnext=0+nsmooth;
  int nprev=azbins-nsmooth-ex2;
  double delnext=nsmooth*delaz;
  int kcross=azbins/2+n*nsmooth;
  //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<endl;
  
  /*Two of the bins can be slightly larger than nsmooth*delaz, affects 4 second derivs*/
  //for unequal spacing
  double h2=delnext;
  double h1=(nsmooth+ex2)*delaz;
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j][kcross]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j][kcross])+1.0/(st*st*h1*h2*(h1+h2))*2.0*(h1*prt[i][j][nnext]-(h1+h2)*prt[i][j][k]+h2*prt[i][j][nprev]) );//for i=inter, j=0, k=0
      
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;

  int nhalf=nbins/2;
  double del1=delnext;
  //for(k=nsmooth;k<azbins-nsmooth;k+=nsmooth){
  for(n=1;n<nhalf-1;n++){
    k=n*nsmooth;
    nnext=(n+1)*nsmooth;
    nprev=(n-1)*nsmooth;
    
    //if(n<nbins/2)
    ncross=n+nbins/2;
    kcross=azbins/2+n*nsmooth;
    //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
    //else
    //ncross=n-nbins/2;
    
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j][kcross]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j][kcross])+1.0/(st*st*del1*del1)*(prt[i][j][nnext]-2.0*prt[i][j][k]+prt[i][j][nprev]) );//for i=inter, j=0, k interior
    
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  if(nbins>2){
    //at nhalf-1, the next bin is possibly larger by ex1 bins
    n=nhalf-1;
    k=n*nsmooth;
    ncross=n+nbins/2;
    kcross=azbins/2+n*nsmooth;
    nnext=(n+1)*nsmooth+ex1;
    nprev=(n-1)*nsmooth;
    //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
    delnext=(nsmooth+ex1)*delaz;
    h2=delnext;
    h1=del1;
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j][kcross]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j][kcross])+1.0/(st*st*h1*h2*(h1+h2))*2.0*(h1*prt[i][j][nnext]-(h1+h2)*prt[i][j][k]+h2*prt[i][j][nprev]) );//for i=inter, j=0, k interior
    
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    
    /*last in half*/
    n=nhalf;
    k=nsmooth*nhalf+ex1;
    nnext=(nhalf+1)*nsmooth+ex1;
    nprev=(nhalf-1)*nsmooth;
    ncross=n-nbins/2;
    kcross=0;
    //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
    h1=h2;
    h2=del1;
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j][kcross]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j][kcross])+1.0/(st*st*h1*h2*(h1+h2))*2.0*(h1*prt[i][j][nnext]-(h1+h2)*prt[i][j][k]+h2*prt[i][j][nprev]) );//for i=inter, j=0, k interior
    
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  
  delnext=(nsmooth)*delaz;//back to normal grid spacing
  /*Now do second half*/
  for(n=nhalf+1;n<nbins-1;n++){
    k=n*nsmooth+ex1;//add in extra nodes added to bind nhalf

    nnext=(n+1)*nsmooth+ex1;
    nprev=(n-1)*nsmooth+ex1;
    ncross=n-nbins/2;
    kcross=ncross*nsmooth;
    //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j][kcross]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j][kcross])+1.0/(st*st*del1*del1)*(prt[i][j][nnext]-2.0*prt[i][j][k]+prt[i][j][nprev]) );//for i=inter, j=0, k interior
    
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  n=nbins-1;
  k=azbins-nsmooth-ex2;//size of last bin is nsmooth+ex2
  nprev=k-nsmooth;
  nnext=0;
  ncross=n-nbins/2;
  kcross=ncross*nsmooth;
  //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
  delnext=(nsmooth+ex2)*delaz;
  h2=delnext;
  h1=del1;
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j][kcross]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j][kcross])+1.0/(st*st*h1*h2*(h1+h2))*2.0*(h1*prt[i][j][0]-(h1+h2)*prt[i][j][k]+h2*prt[i][j][nprev]) );//for i=inter, j=0, k=end
  
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  
  
  //interpolate all the nodes between
  delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int t;
  int kprev;
  for(n=1;n<nhalf;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
    for(t=0;t<nsmooth-1;t++){
      prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
	passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
	psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    }
  }
  n=nhalf;
  k=n*nsmooth+ex1;
  kprev=(n-1)*nsmooth;
  int ns=nsmooth+ex1;
  delnext=(nsmooth+ex1)*delaz;
  dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
  for(t=0;t<ns-1;t++){
    prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
    passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  delnext=(nsmooth)*delaz;//back to normal grid spacing  

  for(n=nhalf+1;n<nbins;n++){
    k=n*nsmooth+ex1;
    kprev=k-nsmooth;
    dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
    for(t=0;t<nsmooth-1;t++){
      prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
	passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
	psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    }
  }
  kprev=azbins-nsmooth-ex2;
  delnext=(nsmooth+ex2)*delaz;
  k=0;//last bin
  ns=azbins-kprev;
  dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
  for(t=0;t<ns-1;t++){
    prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
    passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  
  
}


void fullk_bc_spole(int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2)
{
  /*For south pole, instead of j+1, use j, kcross (the azimuth across the pole, at the same polar angle)*/
  int k=0;//BC, az=0.
  int kcross=k+azbins/2;
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-2.0*delR*BC*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j][kcross]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j][kcross]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][azbins-1]) );//for i=0, j=end, k=0
  passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  
  for(k=1;k<azbins-1;k++){
    if(k<azbins/2)
      kcross=k+azbins/2;
    else
      kcross=k-azbins/2;
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-2.0*delR*BC*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j][kcross]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j][kcross]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][k-1]) );//for i=0, j=end, k interior
    passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  k=azbins-1;//BC, az=2pi-delta.
  kcross=k-azbins/2;
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-2.0*delR*BC*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j][kcross]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j][kcross]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][k-1]) );//for i=0, j=end, k=end
  passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
}
void gridk_bc_spole(int nbins, int nsmooth, int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2)
{
  int k=0;//BC, az=0.
  int ncross=0+nbins/2;
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;

  int bleft, eleft;
  int ex1=0;
  int ex2=0;
  if(azbins%nbins!=0){
    bleft=azbins-(nbins-1)*nsmooth;
    eleft=bleft-nsmooth;
    ex1=eleft/2;
    ex2=eleft-ex1;
  }
  //cout <<"In gridk BC Spole: nbins, "<<nbins<<" nsmooth: "<<nsmooth<<" last bin: "<<bleft<<" leftover: "<<eleft<<" ex1, ex2: "<<ex1<<' '<<ex2<<endl;
  double dpda;
  /*First half of azimuths*/
  int n=0;
  int nnext=0+nsmooth;
  int nprev=azbins-nsmooth-ex2;
  double delnext=nsmooth*delaz;
  int kcross=azbins/2+n*nsmooth;
  //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<endl;
  
  /*Two of the bins can be slightly larger than nsmooth*delaz, affects 4 second derivs*/
  //for unequal spacing
  double h2=delnext;
  double h1=(nsmooth+ex2)*delaz;
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j][kcross]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j][kcross]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*h1*h2*(h1+h2))*2.0*(h1*prt[i][j][nnext]-(h1+h2)*prt[i][j][k]+h2*prt[i][j][nprev]) );//for i=inter, j=0, k=0
      
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;

  int nhalf=nbins/2;
  double del1=delnext;
  //for(k=nsmooth;k<azbins-nsmooth;k+=nsmooth){
  for(n=1;n<nhalf-1;n++){
    k=n*nsmooth;
    nnext=(n+1)*nsmooth;
    nprev=(n-1)*nsmooth;
    
    //if(n<nbins/2)
    ncross=n+nbins/2;
    kcross=azbins/2+n*nsmooth;
    //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
    //else
    //ncross=n-nbins/2;
    
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j][kcross]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j][kcross]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][nnext]-2.0*prt[i][j][k]+prt[i][j][nprev]) );//for i=inter, j=0, k interior
    
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  if(nbins>2){
    //at nhalf-1, the next bin is possibly larger by ex1 bins
    n=nhalf-1;
    k=n*nsmooth;
    ncross=n+nbins/2;
    kcross=azbins/2+n*nsmooth;
    nnext=(n+1)*nsmooth+ex1;
    nprev=(n-1)*nsmooth;
    //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
    delnext=(nsmooth+ex1)*delaz;
    h2=delnext;
    h1=del1;
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j][kcross]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j][kcross]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*h1*h2*(h1+h2))*2.0*(h1*prt[i][j][nnext]-(h1+h2)*prt[i][j][k]+h2*prt[i][j][nprev]) );//for i=inter, j=0, k interior
    
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    
    /*last in half*/
    n=nhalf;
    k=nsmooth*nhalf+ex1;
    nnext=(nhalf+1)*nsmooth+ex1;
    nprev=(nhalf-1)*nsmooth;
    ncross=n-nbins/2;
    kcross=0;
    //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
    h1=h2;
    h2=del1;
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j][kcross]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j][kcross]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*h1*h2*(h1+h2))*2.0*(h1*prt[i][j][nnext]-(h1+h2)*prt[i][j][k]+h2*prt[i][j][nprev]) );//for i=inter, j=0, k interior
    
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  
  delnext=(nsmooth)*delaz;//back to normal grid spacing
  /*Now do second half*/
  for(n=nhalf+1;n<nbins-1;n++){
    k=n*nsmooth+ex1;//add in extra nodes added to bind nhalf

    nnext=(n+1)*nsmooth+ex1;
    nprev=(n-1)*nsmooth+ex1;
    ncross=n-nbins/2;
    kcross=ncross*nsmooth;
    //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j][kcross]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j][kcross]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][nnext]-2.0*prt[i][j][k]+prt[i][j][nprev]) );//for i=inter, j=0, k interior
    
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  n=nbins-1;
  k=azbins-nsmooth-ex2;//size of last bin is nsmooth+ex2
  nprev=k-nsmooth;
  nnext=0;
  ncross=n-nbins/2;
  kcross=ncross*nsmooth;
  //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
  delnext=(nsmooth+ex2)*delaz;
  h2=delnext;
  h1=del1;
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j][kcross]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j][kcross]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*h1*h2*(h1+h2))*2.0*(h1*prt[i][j][0]-(h1+h2)*prt[i][j][k]+h2*prt[i][j][nprev]) );//for i=inter, j=0, k=end
  
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  
  
  //interpolate all the nodes between
  delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int t;
  int kprev;
  for(n=1;n<nhalf;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
    for(t=0;t<nsmooth-1;t++){
      prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
	passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
	psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    }
  }
  n=nhalf;
  k=n*nsmooth+ex1;
  kprev=(n-1)*nsmooth;
  int ns=nsmooth+ex1;
  delnext=(nsmooth+ex1)*delaz;
  dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
  for(t=0;t<ns-1;t++){
    prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
    passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  delnext=(nsmooth)*delaz;//back to normal grid spacing  

  for(n=nhalf+1;n<nbins;n++){
    k=n*nsmooth+ex1;
    kprev=k-nsmooth;
    dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
    for(t=0;t<nsmooth-1;t++){
      prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
	passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
	psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    }
  }
  kprev=azbins-nsmooth-ex2;
  delnext=(nsmooth+ex2)*delaz;
  k=0;//last bin
  ns=azbins-kprev;
  dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
  for(t=0;t<ns-1;t++){
    prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
    passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  
  
}

void fullk_bc(int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2)
{
  /*interior polar angle*/
  int k=0;
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-2.0*delR*BC*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][azbins-1]) );//for i=0, j=interior, k=0
  passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  for(k=1;k<azbins-1;k++){
    
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-2.0*delR*BC*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][k-1]) );//for i=0, j=interior, k interior
    passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  k=azbins-1;
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-2.0*delR*BC*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][k-1]) );//for i=0, j=interior, k=end
  passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
}
void gridk_bc(int nbins, int nsmooth, int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2)
{
  int k=0;//BC, az=0.
  int ncross=0+nbins/2;
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;

  int bleft, eleft;
  int ex1=0;
  int ex2=0;
  if(azbins%nbins!=0){
    bleft=azbins-(nbins-1)*nsmooth;
    eleft=bleft-nsmooth;
    ex1=eleft/2;
    ex2=eleft-ex1;
  }
  //cout <<"In gridk BC: nbins, "<<nbins<<" nsmooth: "<<nsmooth<<" last bin: "<<bleft<<" leftover: "<<eleft<<" ex1, ex2: "<<ex1<<' '<<ex2<<endl;
  double dpda;
  /*First half of azimuths*/
  int n=0;
  int nnext=0+nsmooth;
  int nprev=azbins-nsmooth-ex2;
  double delnext=nsmooth*delaz;
  int kcross=azbins/2+n*nsmooth;
  //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<endl;
  
  /*Two of the bins can be slightly larger than nsmooth*delaz, affects 4 second derivs*/
  //for unequal spacing
  double h2=delnext;
  double h1=(nsmooth+ex2)*delaz;
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*h1*h2*(h1+h2))*2.0*(h1*prt[i][j][nnext]-(h1+h2)*prt[i][j][k]+h2*prt[i][j][nprev]) );//for i=inter, j=0, k=0
      
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;

  int nhalf=nbins/2;
  double del1=delnext;
  //for(k=nsmooth;k<azbins-nsmooth;k+=nsmooth){
  for(n=1;n<nhalf-1;n++){
    k=n*nsmooth;
    nnext=(n+1)*nsmooth;
    nprev=(n-1)*nsmooth;
    
    //if(n<nbins/2)
    ncross=n+nbins/2;
    kcross=azbins/2+n*nsmooth;
    //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
    //else
    //ncross=n-nbins/2;
    
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][nnext]-2.0*prt[i][j][k]+prt[i][j][nprev]) );//for i=inter, j=0, k interior
    
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  //at nhalf-1, the next bin is possibly larger by ex1 bins
  n=nhalf-1;
  k=n*nsmooth;
  ncross=n+nbins/2;
  kcross=azbins/2+n*nsmooth;
  nnext=(n+1)*nsmooth+ex1;
  nprev=(n-1)*nsmooth;
  //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
  delnext=(nsmooth+ex1)*delaz;
  h2=delnext;
  h1=del1;
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*h1*h2*(h1+h2))*2.0*(h1*prt[i][j][nnext]-(h1+h2)*prt[i][j][k]+h2*prt[i][j][nprev]) );//for i=inter, j=0, k interior
  
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  
  /*last in half*/
  n=nhalf;
  k=nsmooth*nhalf+ex1;
  nnext=(nhalf+1)*nsmooth+ex1;
  nprev=(nhalf-1)*nsmooth;
  ncross=n-nbins/2;
  kcross=0;
  //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
  h1=h2;
  h2=del1;
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*h1*h2*(h1+h2))*2.0*(h1*prt[i][j][nnext]-(h1+h2)*prt[i][j][k]+h2*prt[i][j][nprev]) );//for i=inter, j=0, k interior
    
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  
  
  delnext=(nsmooth)*delaz;//back to normal grid spacing
  /*Now do second half*/
  for(n=nhalf+1;n<nbins-1;n++){
    k=n*nsmooth+ex1;//add in extra nodes added to bind nhalf

    nnext=(n+1)*nsmooth+ex1;
    nprev=(n-1)*nsmooth+ex1;
    ncross=n-nbins/2;
    kcross=ncross*nsmooth;
    //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][nnext]-2.0*prt[i][j][k]+prt[i][j][nprev]) );//for i=inter, j=0, k interior
    
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  n=nbins-1;
  k=azbins-nsmooth-ex2;//size of last bin is nsmooth+ex2
  nprev=k-nsmooth;
  nnext=0;
  ncross=n-nbins/2;
  kcross=ncross*nsmooth;
  //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
  delnext=(nsmooth+ex2)*delaz;
  h2=delnext;
  h1=del1;
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+(prt[i+1][j][k]-BC*2.0*delR*prt[i][j][k]))+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*h1*h2*(h1+h2))*2.0*(h1*prt[i][j][0]-(h1+h2)*prt[i][j][k]+h2*prt[i][j][nprev]) );//for i=inter, j=0, k=end
  
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  
  
  //interpolate all the nodes between
  delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int t;
  int kprev;
  for(n=1;n<nhalf;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
    for(t=0;t<nsmooth-1;t++){
      prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
	passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
	psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    }
  }
  n=nhalf;
  k=n*nsmooth+ex1;
  kprev=(n-1)*nsmooth;
  int ns=nsmooth+ex1;
  delnext=(nsmooth+ex1)*delaz;
  dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
  for(t=0;t<ns-1;t++){
    prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
    passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  delnext=(nsmooth)*delaz;//back to normal grid spacing  

  for(n=nhalf+1;n<nbins;n++){
    k=n*nsmooth+ex1;
    kprev=k-nsmooth;
    dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
    for(t=0;t<nsmooth-1;t++){
      prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
	passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
	psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    }
  }
  kprev=azbins-nsmooth-ex2;
  delnext=(nsmooth+ex2)*delaz;
  k=0;//last bin
  ns=azbins-kprev;
  dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
  for(t=0;t<ns-1;t++){
    prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
    passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  
  
}


void fullk(int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2)
{
  int k=0;
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][azbins-1]) );//for i=inter, j=interior, k=0
  
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;

  for(k=1;k<azbins-1;k++){
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][k-1]) );
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }//done looping over az
  k=azbins-1;
  
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][k-1]) );//for i=inter, j=interior, k=end
  
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;

}
void gridk(int nbins, int nsmooth, int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2)
{
  int k=0;//BC, az=0.
  int ncross=0+nbins/2;
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;

  int bleft, eleft;
  int ex1=0;
  int ex2=0;
  if(azbins%nbins!=0){
    bleft=azbins-(nbins-1)*nsmooth;
    eleft=bleft-nsmooth;
    ex1=eleft/2;
    ex2=eleft-ex1;
  }
  //cout <<"In gridk: nbins, "<<nbins<<" nsmooth: "<<nsmooth<<" last bin: "<<bleft<<" leftover: "<<eleft<<" ex1, ex2: "<<ex1<<' '<<ex2<<endl;
  double dpda;
  /*First half of azimuths*/
  int n=0;
  int nnext=0+nsmooth;
  int nprev=azbins-nsmooth-ex2;
  double delnext=nsmooth*delaz;
  int kcross=azbins/2+n*nsmooth;
  //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<endl;
  
  /*Two of the bins can be slightly larger than nsmooth*delaz, affects 4 second derivs*/
  //for unequal spacing
  double del1=delnext;
  double h2=delnext;
  double h1=(nsmooth+ex2)*delaz;
  if(ex1==0){
    /*For uneven grid spacing, second deriv becomes 2/(h1*h2*(h1+h2))*(h1*p_i+1 -(h1+h2)p_i + h2*p_i-1)*/
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][nnext]-2.0*prt[i][j][k]+prt[i][j][nprev]) );//for i=inter, j=0, k=0
  }else{
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st)*(2.0/(h2*(h1+h2))*prt[i][j][nnext]-2.0/(h1*h2)*prt[i][j][k]+2.0/(h1*(h1+h2))*prt[i][j][nprev]) );//for i=inter, j=0, k=0
  }
      
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;

  int nhalf=nbins/2;

  //for(k=nsmooth;k<azbins-nsmooth;k+=nsmooth){
  for(n=1;n<nhalf-1;n++){
    k=n*nsmooth;
    nnext=(n+1)*nsmooth;
    nprev=(n-1)*nsmooth;
    
    //if(n<nbins/2)
    ncross=n+nbins/2;
    kcross=azbins/2+n*nsmooth;
    //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
    //else
    //ncross=n-nbins/2;
    
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][nnext]-2.0*prt[i][j][k]+prt[i][j][nprev]) );//for i=inter, j=0, k interior
    
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  //at nhalf-1, the next bin is possibly larger by ex1 bins
  n=nhalf-1;
  k=n*nsmooth;
  ncross=n+nbins/2;
  kcross=azbins/2+n*nsmooth;
  nnext=(n+1)*nsmooth+ex1;
  nprev=(n-1)*nsmooth;
  //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
  delnext=(nsmooth+ex1)*delaz;
  h2=delnext;
  h1=del1;
  if(ex1==0){
    /*For uneven grid spacing, second deriv becomes 2/(h1*h2*(h1+h2))*(h1*p_i+1 -(h1+h2)p_i + h2*p_i-1)*/
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][nnext]-2.0*prt[i][j][k]+prt[i][j][nprev]) );//for i=inter, j=0, k=0
  }else{
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st)*(2.0/(h2*(h1+h2))*prt[i][j][nnext]-2.0/(h1*h2)*prt[i][j][k]+2.0/(h1*(h1+h2))*prt[i][j][nprev]) );//for i=inter, j=0, k=0
  }
  
  
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  
  /*last in half*/
  n=nhalf;
  k=nsmooth*nhalf+ex1;
  nnext=(nhalf+1)*nsmooth+ex1;
  nprev=(nhalf-1)*nsmooth;
  ncross=n-nbins/2;
  kcross=0;
  //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
  h1=h2;
  h2=del1;
  if(ex1==0){
    /*For uneven grid spacing, second deriv becomes 2/(h1*h2*(h1+h2))*(h1*p_i+1 -(h1+h2)p_i + h2*p_i-1)*/
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][nnext]-2.0*prt[i][j][k]+prt[i][j][nprev]) );//for i=inter, j=0, k=0
  }else{
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st)*(2.0/(h2*(h1+h2))*prt[i][j][nnext]-2.0/(h1*h2)*prt[i][j][k]+2.0/(h1*(h1+h2))*prt[i][j][nprev]) );//for i=inter, j=0, k=0
  }
  
    
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  
  
  delnext=(nsmooth)*delaz;//back to normal grid spacing
  /*Now do second half*/
  for(n=nhalf+1;n<nbins-1;n++){
    k=n*nsmooth+ex1;//add in extra nodes added to bind nhalf

    nnext=(n+1)*nsmooth+ex1;
    nprev=(n-1)*nsmooth+ex1;
    ncross=n-nbins/2;
    kcross=ncross*nsmooth;
    //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][nnext]-2.0*prt[i][j][k]+prt[i][j][nprev]) );//for i=inter, j=0, k interior
    
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  n=nbins-1;
  k=azbins-nsmooth-ex2;//size of last bin is nsmooth+ex2
  nprev=k-nsmooth;
  nnext=0;
  ncross=n-nbins/2;
  kcross=ncross*nsmooth;
  //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
  delnext=(nsmooth+ex2)*delaz;
  h2=delnext;
  h1=del1;
  if(ex1==0){
    /*For uneven grid spacing, second deriv becomes 2/(h1*h2*(h1+h2))*(h1*p_i+1 -(h1+h2)p_i + h2*p_i-1)*/
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][nnext]-2.0*prt[i][j][k]+prt[i][j][nprev]) );//for i=inter, j=0, k=0
  }else{
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st)*(2.0/(h2*(h1+h2))*prt[i][j][nnext]-2.0/(h1*h2)*prt[i][j][k]+2.0/(h1*(h1+h2))*prt[i][j][nprev]) );//for i=inter, j=0, k=0
  }
  
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  
  
  //interpolate all the nodes between
  delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int t;
  int kprev;
  for(n=1;n<nhalf;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
    for(t=0;t<nsmooth-1;t++){
      prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
	passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
	psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    }
  }
  n=nhalf;
  k=n*nsmooth+ex1;
  kprev=(n-1)*nsmooth;
  
  delnext=(nsmooth+ex1)*delaz;
  int ns=nsmooth+ex1;
  dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
  for(t=0;t<ns-1;t++){
    prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
    passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  delnext=(nsmooth)*delaz;//back to normal grid spacing  

  for(n=nhalf+1;n<nbins;n++){
    k=n*nsmooth+ex1;
    kprev=k-nsmooth;
    dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
    for(t=0;t<nsmooth-1;t++){
      prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
	passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
	psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    }
  }
  kprev=azbins-nsmooth-ex2;
  delnext=(nsmooth+ex2)*delaz;
  k=0;//last bin
  ns=azbins-kprev;
  dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
  for(t=0;t<ns-1;t++){
    prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
    passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  
  
}

void gridk_old(int nsmooth, int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2)
{
  int k=0;
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  int kprev;
  double delgrid;
  double dpda;
  
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][azbins-1]) );//for i=inter, j=interior, k=0
  
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;

  for(k=nsmooth;k<azbins-nsmooth;k+=nsmooth){
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][k-1]) );
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }//done looping over az
  if(nsmooth==1){
    k=azbins-1;
    
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][k-1]) );//for i=inter, j=interior, k=end
    
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }else{
    //skip the last point, it's adjacent to the start node.
    //interpolate all the nodes between
    delgrid=nsmooth*delaz;
    int t;
    for(k=nsmooth;k<azbins-nsmooth;k+=nsmooth){
      kprev=k-nsmooth;
      dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delgrid;
      for(t=0;t<nsmooth-1;t++){
	prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz;
	passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
	psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
      }
      
    }
    //do final set of bins to zero
    kprev=k-nsmooth;
    k=0;
    int ns=azbins-kprev;
    delgrid=ns*delaz;
    dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delgrid;
    for(t=0;t<ns-1;t++){
      prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz;
      passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
      psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    }

  }

}
void fullk_npole(int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2)
{
  int k=0;//BC, az=0.
  int kcross=k+azbins/2;
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j][kcross]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j][kcross])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][azbins-1]) );//for i=inter, j=0, k=0
      
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  
  for(k=1;k<azbins-1;k++){
    if(k<azbins/2)
      kcross=k+azbins/2;
    else
      kcross=k-azbins/2;
    
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j][kcross]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j][kcross])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][k-1]) );//for i=inter, j=0, k interior
    
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  k=azbins-1;//BC, az=2pi-delta.
  kcross=k-azbins/2;
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j][kcross]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j][kcross])+1.0/(st*st*delaz2)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][k-1]) );//for i=inter, j=0, k=end
  
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
}
void gridk_npole(int nbins, int nsmooth, int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2)
{
  int k=0;//BC, az=0.
  int ncross=0+nbins/2;
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;

  int bleft, eleft;
  int ex1=0;
  int ex2=0;
  if(azbins%nbins!=0){
    bleft=azbins-(nbins-1)*nsmooth;
    eleft=bleft-nsmooth;
    ex1=eleft/2;
    ex2=eleft-ex1;
  }
  //cout <<"In gridk Npole: nbins, "<<nbins<<" nsmooth: "<<nsmooth<<" last bin: "<<bleft<<" leftover: "<<eleft<<" ex1, ex2: "<<ex1<<' '<<ex2<<endl;
  double dpda;
  /*First half of azimuths*/
  int n=0;
  int nnext=0+nsmooth;
  int nprev=azbins-nsmooth-ex2;
  double delnext=nsmooth*delaz;
  int kcross=azbins/2+n*nsmooth;
  //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<endl;
  
  /*Two of the bins can be slightly larger than nsmooth*delaz, affects 4 second derivs*/
  //for unequal spacing
  double h2=delnext;
  double h1=(nsmooth+ex2)*delaz;
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j][kcross]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j][kcross])+1.0/(st*st*h1*h2*(h1+h2))*2.0*(h1*prt[i][j][nnext]-(h1+h2)*prt[i][j][k]+h2*prt[i][j][nprev]) );//for i=inter, j=0, k=0
      
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;

  int nhalf=nbins/2;
  double del1=delnext;
  //for(k=nsmooth;k<azbins-nsmooth;k+=nsmooth){
  for(n=1;n<nhalf-1;n++){
    k=n*nsmooth;
    nnext=(n+1)*nsmooth;
    nprev=(n-1)*nsmooth;
    
    //if(n<nbins/2)
    ncross=n+nbins/2;
    kcross=azbins/2+n*nsmooth;
    //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
    //else
    //ncross=n-nbins/2;
    
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j][kcross]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j][kcross])+1.0/(st*st*del1*del1)*(prt[i][j][nnext]-2.0*prt[i][j][k]+prt[i][j][nprev]) );//for i=inter, j=0, k interior
    
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  if(nbins>2){
    //at nhalf-1, the next bin is possibly larger by ex1 bins
    n=nhalf-1;
    k=n*nsmooth;
    ncross=n+nbins/2;
    kcross=azbins/2+n*nsmooth;
    nnext=(n+1)*nsmooth+ex1;
    nprev=(n-1)*nsmooth;
    //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
    delnext=(nsmooth+ex1)*delaz;
    h2=delnext;
    h1=del1;
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j][kcross]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j][kcross])+1.0/(st*st*h1*h2*(h1+h2))*2.0*(h1*prt[i][j][nnext]-(h1+h2)*prt[i][j][k]+h2*prt[i][j][nprev]) );//for i=inter, j=0, k interior
    
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    
    /*last in half*/
    n=nhalf;
    k=nsmooth*nhalf+ex1;
    nnext=(nhalf+1)*nsmooth+ex1;
    nprev=(nhalf-1)*nsmooth;
    ncross=n-nbins/2;
    kcross=0;
    //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
    h1=h2;
    h2=del1;
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j][kcross]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j][kcross])+1.0/(st*st*h1*h2*(h1+h2))*2.0*(h1*prt[i][j][nnext]-(h1+h2)*prt[i][j][k]+h2*prt[i][j][nprev]) );//for i=inter, j=0, k interior
    
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  
  delnext=(nsmooth)*delaz;//back to normal grid spacing
  /*Now do second half*/
  for(n=nhalf+1;n<nbins-1;n++){
    k=n*nsmooth+ex1;//add in extra nodes added to bind nhalf

    nnext=(n+1)*nsmooth+ex1;
    nprev=(n-1)*nsmooth+ex1;
    ncross=n-nbins/2;
    kcross=ncross*nsmooth;
    //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j][kcross]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j][kcross])+1.0/(st*st*del1*del1)*(prt[i][j][nnext]-2.0*prt[i][j][k]+prt[i][j][nprev]) );//for i=inter, j=0, k interior
    
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  n=nbins-1;
  k=azbins-nsmooth-ex2;//size of last bin is nsmooth+ex2
  nprev=k-nsmooth;
  nnext=0;
  ncross=n-nbins/2;
  kcross=ncross*nsmooth;
  //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
  delnext=(nsmooth+ex2)*delaz;
  h2=delnext;
  h1=del1;
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j][kcross]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j][kcross])+1.0/(st*st*h1*h2*(h1+h2))*2.0*(h1*prt[i][j][0]-(h1+h2)*prt[i][j][k]+h2*prt[i][j][nprev]) );//for i=inter, j=0, k=end
  
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  
  
  //interpolate all the nodes between
  delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int t;
  int kprev;
  for(n=1;n<nhalf;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
    for(t=0;t<nsmooth-1;t++){
      prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
	passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
	psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    }
  }
  n=nhalf;
  k=n*nsmooth+ex1;
  kprev=(n-1)*nsmooth;
  delnext=(nsmooth+ex1)*delaz;
  int ns=nsmooth+ex1;
  dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
  for(t=0;t<ns-1;t++){
    prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
    passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  delnext=(nsmooth)*delaz;//back to normal grid spacing  

  for(n=nhalf+1;n<nbins;n++){
    k=n*nsmooth+ex1;
    kprev=k-nsmooth;
    dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
    for(t=0;t<nsmooth-1;t++){
      prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
	passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
	psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    }
  }
  kprev=azbins-nsmooth-ex2;
  delnext=(nsmooth+ex2)*delaz;
  k=0;//last bin
  ns=azbins-kprev;
  dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
  for(t=0;t<ns-1;t++){
    prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
    passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  
  
}

void fullk_spole(int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2)
{
  
  int k=0;//BC, az=0.
  int kcross=k+azbins/2;
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j][kcross]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j][kcross]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][azbins-1]) );//for i=inter, j=end, k=0
  
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  
  for(k=1;k<azbins-1;k++){
    if(k<azbins/2)
      kcross=k+azbins/2;
    else
      kcross=k-azbins/2;
    
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j][kcross]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j][kcross]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][k+1]-2.0*prt[i][j][k]+prt[i][j][k-1]) );//for i=inter, j=end, k interior
    
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  k=azbins-1;//BC, az=2pi-delta.
  kcross=k-azbins/2;
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j][kcross]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j][kcross]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][k-1]) );//for i=inter, j=end, k=end
  
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
}    

void gridk_spole(int nbins, int nsmooth, int i, int j, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ***prt, double ***prt2)
{
  int k=0;//BC, az=0.
  int ncross=0+nbins/2;
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;

  int bleft, eleft;
  int ex1=0;
  int ex2=0;
  if(azbins%nbins!=0){
    bleft=azbins-(nbins-1)*nsmooth;
    eleft=bleft-nsmooth;
    ex1=eleft/2;
    ex2=eleft-ex1;
  }
  //cout <<"In gridk Spole: nbins, "<<nbins<<" nsmooth: "<<nsmooth<<" last bin: "<<bleft<<" leftover: "<<eleft<<" ex1, ex2: "<<ex1<<' '<<ex2<<endl;
  double dpda;
  /*First half of azimuths*/
  int n=0;
  int nnext=0+nsmooth;
  int nprev=azbins-nsmooth-ex2;
  double delnext=nsmooth*delaz;
  int kcross=azbins/2+n*nsmooth;
  //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<endl;
  
  /*Two of the bins can be slightly larger than nsmooth*delaz, affects 4 second derivs*/
  //for unequal spacing
  double h2=delnext;
  double h1=(nsmooth+ex2)*delaz;
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j][kcross]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j][kcross]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*h1*h2*(h1+h2))*2.0*(h1*prt[i][j][nnext]-(h1+h2)*prt[i][j][k]+h2*prt[i][j][nprev]) );//for i=inter, j=0, k=0
      
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;

  int nhalf=nbins/2;
  double del1=delnext;
  //for(k=nsmooth;k<azbins-nsmooth;k+=nsmooth){
  for(n=1;n<nhalf-1;n++){
    k=n*nsmooth;
    nnext=(n+1)*nsmooth;
    nprev=(n-1)*nsmooth;
    
    //if(n<nbins/2)
    ncross=n+nbins/2;
    kcross=azbins/2+n*nsmooth;
    //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
    //else
    //ncross=n-nbins/2;
    
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j][kcross]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j][kcross]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][nnext]-2.0*prt[i][j][k]+prt[i][j][nprev]) );//for i=inter, j=0, k interior
    
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  if(nbins>2){
    //at nhalf-1, the next bin is possibly larger by ex1 bins
    n=nhalf-1;
    k=n*nsmooth;
    ncross=n+nbins/2;
    kcross=azbins/2+n*nsmooth;
    nnext=(n+1)*nsmooth+ex1;
    nprev=(n-1)*nsmooth;
    //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
    delnext=(nsmooth+ex1)*delaz;
    h2=delnext;
    h1=del1;
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j][kcross]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j][kcross]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*h1*h2*(h1+h2))*2.0*(h1*prt[i][j][nnext]-(h1+h2)*prt[i][j][k]+h2*prt[i][j][nprev]) );//for i=inter, j=0, k interior
    
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    
    /*last in half*/
    n=nhalf;
    k=nsmooth*nhalf+ex1;
    nnext=(nhalf+1)*nsmooth+ex1;
    nprev=(nhalf-1)*nsmooth;
    ncross=n-nbins/2;
    kcross=0;
    //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
    h1=h2;
    h2=del1;
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j][kcross]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j][kcross]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*h1*h2*(h1+h2))*2.0*(h1*prt[i][j][nnext]-(h1+h2)*prt[i][j][k]+h2*prt[i][j][nprev]) );//for i=inter, j=0, k interior
    
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  
  delnext=(nsmooth)*delaz;//back to normal grid spacing
  /*Now do second half*/
  for(n=nhalf+1;n<nbins-1;n++){
    k=n*nsmooth+ex1;//add in extra nodes added to bind nhalf

    nnext=(n+1)*nsmooth+ex1;
    nprev=(n-1)*nsmooth+ex1;
    ncross=n-nbins/2;
    kcross=ncross*nsmooth;
    //cout <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
    prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j][kcross]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j][kcross]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][nnext]-2.0*prt[i][j][k]+prt[i][j][nprev]) );//for i=inter, j=0, k interior
    
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  n=nbins-1;
  k=azbins-nsmooth-ex2;//size of last bin is nsmooth+ex2
  nprev=k-nsmooth;
  nnext=0;
  ncross=n-nbins/2;
  kcross=ncross*nsmooth;
  //cout    <<"bin: "<<n<<" k:: "<<k<<" next bin: "<<nnext<<" prev bin: "<<nprev<<" ncross: "<<ncross<<" kcross: "<<kcross<<endl;
  delnext=(nsmooth+ex2)*delaz;
  h2=delnext;
  h1=del1;
  prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j][kcross]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j][kcross]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*h1*h2*(h1+h2))*2.0*(h1*prt[i][j][0]-(h1+h2)*prt[i][j][k]+h2*prt[i][j][nprev]) );//for i=inter, j=0, k=end
  
  psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  
  
  //interpolate all the nodes between
  delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int t;
  int kprev;
  for(n=1;n<nhalf;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
    for(t=0;t<nsmooth-1;t++){
      prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
	passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
	psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    }
  }
  n=nhalf;
  k=n*nsmooth+ex1;
  kprev=(n-1)*nsmooth;
  delnext=(nsmooth+ex1)*delaz;
  int ns=nsmooth+ex1;
  dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
  for(t=0;t<ns-1;t++){
    prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
    passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  delnext=(nsmooth)*delaz;//back to normal grid spacing  

  for(n=nhalf+1;n<nbins;n++){
    k=n*nsmooth+ex1;
    kprev=k-nsmooth;
    dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
    for(t=0;t<nsmooth-1;t++){
      prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
	passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
	psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
    }
  }
  kprev=azbins-nsmooth-ex2;
  delnext=(nsmooth+ex2)*delaz;
  k=0;//last bin
  ns=azbins-kprev;
  dpda=(prt2[i][j][k]-prt2[i][j][kprev])/delnext;
  for(t=0;t<ns-1;t++){
    prt2[i][j][kprev+t+1]=prt2[i][j][kprev]+dpda*delaz*(t+1);
    passoc1+=prt2[0][j][k]*ka*deltat*delpol*delaz*st;
    psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }
  
  
}


void free_prop_sphere(double sigma, int nrbins, int polbins, int azbins, double delR, double delpol, double delaz, double Dtot, double time, double r0, double az0, double pol0, double ***pfree, double &norm)
{
  /*pfree=1/(4piDt)^(3/2)*exp(-(rvec-r0vec)^2/(4Dt))*/
  
  double f1=4.0*M_PI*Dtot*time;
  double cof=pow(f1, -3.0/2.0);
  double x0=r0*cos(az0)*sin(pol0);
  double y0=r0*sin(az0)*sin(pol0);
  double z0=r0*cos(pol0);
  int i, j, k;
  double std=4.0*Dtot*time;
  double dx, dy, dz, drvec;
  double x1, y1, z1;
  double r1, az1, pol1;
   norm=0.0;
  double r2;
  double st, ct;
  for(i=0;i<nrbins;i++){
    r1=i*delR+sigma;
    r2=r1*r1;
    for(j=0;j<polbins;j++){
      pol1=(j+0.5)*delpol;
      st=sin(pol1);
      ct=cos(pol1);
      for(k=0;k<azbins;k++){
	az1=k*delaz;
	x1=r1*cos(az1)*st;
	y1=r1*sin(az1)*st;
	z1=r1*ct;
	dx=x1-x0;
	dy=y1-y0;
	dz=z1-z0;
	drvec=dx*dx+dy*dy+dz*dz;
	pfree[i][j][k]=cof*exp(-1.0/std*drvec);
	norm+=pfree[i][j][k]*r2*delR*st*delpol*delaz;
      }
    }
  }

}

void free_prop_cartesian(double sigma, int xbins, int ybins, int zbins, double delx, double dely, double delz, double Dtot, double time, double r0, double az0, double pol0, double ***pfree)
{
  /*pfree=1/(4piDt)^(3/2)*exp(-(rvec-r0vec)^2/(4Dt))*/
  
  double f1=4.0*M_PI*Dtot*time;
  double cof=pow(f1, -3.0/2.0);
  double x0=r0*cos(az0)*sin(pol0);
  double y0=r0*sin(az0)*sin(pol0);
  double z0=r0*cos(pol0);
  x0=0;
  y0=0;
  z0=0;
  int i, j, k;
  double std=4.0*Dtot*time;
  double dx, dy, dz, drvec;
  double x1, y1, z1;
  double r1, az1, pol1;
  double sum=0.0;
  for(i=0;i<xbins;i++){
    x1=(i-xbins/2)*delx;
    for(j=0;j<ybins;j++){
      y1=(j-ybins/2)*dely;
      for(k=0;k<zbins;k++){
	z1=(k-zbins/2)*delz;
	
	dx=x1-x0;
	dy=y1-y0;
	dz=z1-z0;
	drvec=dx*dx+dy*dy+dz*dz;
	pfree[i][j][k]=cof*exp(-1.0/std*drvec);
	sum+=pfree[i][j][k]*delx*dely*delz;
      }
    }
  }
  cout <<"Total integral in cartesian: "<<sum<<endl;
  ofstream ffile;
  char cfname[300];
  sprintf(cfname, "Cart_sol_%g_%g.dat",r0, time);
  ffile.open(cfname);
  double tstart=time;
  for(i=0;i<xbins;i++){
    for(j=0;j<ybins;j++){
      for(k=0;k<zbins;k++){
	ffile<<pfree[i][j][k]<<'\t';
      }
    }
    ffile<<endl;
  }
  
}

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

capgrid routing propagates PDE in full spherical coordinates, smoothing over 
grid points at the poles where the azimuthal bin size becomes too small. 

The calculation of survival probabilities could be improved by using a better integration
routine, like Simpson's. Currently just using midpoint method. 

New PDE to solve if interested in r2, where r2 is the separation between the rotating vector and 
the partner site.

maggie johnson
*/

#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include "full1.h"
#include "gridk_r2.h"
#include "gf_rot.h"

using namespace std;


void integrate_simple_v2(int tbins, int nrbins, double delR, double deltat, double sigma, double ka, double Dtot, double BC, double **prt, double *passoc);
void integrate_spherical_v2(int tbins, int nrbins, double delR, double deltat, double sigma, double ka, double Dtot, double BC, double ***prt, double ***prt2, double *passoc, int polbins, double delpol, int azbins, double delaz, double r0);

void integrate_spherical_r2_capgrid(int tbins, int nrbins, double delR, double deltat, double sigma, double ka, double Dtot, double BC, double ****prt, double ****prt2, double *passoc, int polbins, double delpol, int azbins, double delaz, double r0, double mindel,  double tstart, double az0, double pol0, double delApol, int polAbins, double d1, double Drot);
void free_prop_sphere(double sigma, int nrbins, int polbins, int azbins, double delR, double delpol, double delaz, double Dtot, double time, double r0, double az0, double pol0, double ***pfree, double &norm);
void free_prop_cartesian(double sigma, int xbins, int ybins, int zbins, double delx, double dely, double delz, double Dtot, double time, double r0, double az0, double pol0, double ***pfree);
void get_initial_config(double ***pfree, double **protfree, double ****prt, int nrbins, double delR, double sigma, double d1, int polbins, double delpol, int polAbins, double delApol, int azbins, double delaz, double rlim);


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
   
  int nrbins;
  double maxR=10;
  double sigma=1.0;
  double delR=atof(argv[1]);//(maxR-sigma)/(1.0*nrbins);
  nrbins=int((maxR-sigma)/delR);
  //double *r=new double[nrbins];
  double r;
  int tbins=10;
  double deltat=atof(argv[2]);//1E-4;//us
  int polbins=50;
  int azbins=51;
  int polAbins=50;
  double delApol=M_PI/(1.0*polAbins);
  double delpol=M_PI/(1.0*polbins);
  /*Az starts at zero, last bin is at pi*/
  double delaz=1.0*M_PI/(1.0*(azbins-1));//Only 1PI for az=phi-phiA!!!!!!!!
  double ****prt=new double***[nrbins+1];
  double ***pfree=new double**[nrbins+1];
  double **protfree=new double*[polbins];
  double ****prt2=new double***[nrbins+1];
  for(i=0;i<nrbins+1;i++){
    prt[i]=new double**[polbins];
    pfree[i]=new double*[polbins];
    prt2[i]=new double**[polbins];
  }
  for(i=0;i<polAbins;i++)
    protfree[i]=new double[(azbins-1)*2];
  
  for(i=0;i<nrbins+1;i++){
    for(j=0;j<polbins;j++){
      prt[i][j]=new double*[polAbins];
      prt2[i][j]=new double*[polAbins];
      
      pfree[i][j]=new double[(azbins-1)*2];
    }
  }
  for(i=0;i<nrbins+1;i++){
    for(j=0;j<polbins;j++){
      for(k=0;k<polAbins;k++){
	prt[i][j][k]=new double[azbins];
	prt2[i][j][k]=new double[azbins];
      }
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
  int fullazbins=(azbins-1)*2;
  
  int Lmax=50;
  double polA0=M_PI/2.0;
  //if azA0=-az0, then w=pi (azA0-az0)
  double azA0=M_PI-az0;
  /*Get transport coefficients from Einstein and Einstein Stokes equations*/
  double d1=1.02;
  double kb=1.3806488E-23;
  double Temp=293;
  double prea=16;
  double nu=0.001;
  Dtot=kb*Temp/(6.0*M_PI*nu*prea)*1E27/1E6;
  
  double Drot=kb*Temp/(8.0*M_PI*nu)*1E27/1E6/(prea*prea*prea);
  cout <<"delaz: "<<delaz<<" delpol: "<<delpol<<" azbins: "<<fullazbins<<endl;
  free_prop_sphere(sigma, nrbins, polbins, fullazbins, delR, delpol, delaz, Dtot,  time,  r0,  az0,pol0, pfree, psurvivef1);  
  
  

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
  
  gf_rot_1time(Lmax,polAbins, fullazbins, cos(polA0), azA0, time, protfree, Drot);
  double psrot=0.0;
  for(j=0;j<polAbins;j++){
    theta=(j+0.5)*delApol;
    st=sin(theta);
    for(k=0;k<fullazbins;k++){
      psrot+=protfree[j][k]*st*delApol*delaz;
    }
  }
  cout <<"Total protfree 1?:" <<psrot<<endl;
  for(j=0;j<polAbins;j++){
    for(k=0;k<fullazbins;k++){
      protfree[j][k]/=psrot;
    }
  }
  
  /*Need to sample initial prt[i][j][l][k] by looping over all nonzero values of the 5 variables,
    and when w=phi-phiA adding to prob via p(phi, phiA)*delphi*delphiA/delw
    with each p(r2), add in p(r1)*dr1/dr2/J, where J=dr1/dr2=r2/(r1+d*anglething)
  */
  double rlim=sigma+d1;
  get_initial_config(pfree,  protfree,  prt, nrbins,  delR,  sigma,  d1,  polbins, delpol, polAbins, delApol, azbins,  delaz,  rlim);
  
  
  psurvive=0;

  /*  for(i=0;i<nrbins+1;i++){
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
  */
  
  // for(i=0;i<nrbins+1;i++){
//     for(j=0;j<polbins;j++){
//       for(k=0;k<polAazbins;k++){
// 	prt[i][j][k]=pfree[i][j][k]/psurvivef1;
//       }
//     }
//   }
  
  
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

  ka=ka/(4.0*M_PI*sigma*sigma);//KAPPA!!
  cout <<" Using kappa: "<<ka<<endl;
  double BC=ka/Dtot;//ka/(4.0*M_PI*Dtot*sigma*sigma);
  
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
  integrate_spherical_r2_capgrid(tbins, nrbins, delR, deltat, sigma, ka, Dtot, BC, prt, prt2, passoc, polbins, delpol, azbins, delaz, r0, mindel,tstart, az0, pol0, delApol, polAbins, d1, Drot);
  
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

void integrate_spherical_r2_capgrid(int tbins, int nrbins, double delR, double deltat, double sigma, double ka, double Dtot, double BC, double ****prt, double ****prt2, double *passoc, int polbins, double delpol, int azbins, double delaz, double r0, double mindel,  double tstart, double az0, double pol0, double delApol, int polAbins, double d1, double Drot)
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
  double ct, cta, sta;
  int l;
  double polA;
  
  for(t=1;t<tbins;t++){
    time=t*deltat+tstart;
    if(t%twrite==0){

      sprintf(fname,"p_r2_spheret_v2numer_r0_%g_dt_%g.dat",r0, time);
      cout <<"new file: "<<fname<<" time: "<<time<<endl;
      outfile.open(fname);
      // sprintf(fname,"pfree_v2numer_r0_%g_dt_%g.dat",r0, time);
//       cout <<"new file: "<<fname<<" time: "<<time<<endl;
//       ffile.open(fname);

    }
    /*First do r=sigma with BC*/
    psurvive=0.0;
    r=sigma;
    r2=r*r;
    i=0;//BC, r=sigma.
    
    j=0;//BC, theta=0+delta
    theta=(j+0.5)*delpol;//need to skip the poles.
    st=sin(theta);
    ct=cos(theta);
    l=0;
    polA=(l+0.5)*delApol;
    sta=sin(polA);
    cta=cos(polA);
    
    //for each grid point, decide which azimuthals to include 
    arc=d1*sta*delaz+abs(r-d1)*st*delaz;
    if(arc<mindel){
      //choose a subset of grid points in k to update, interpolate the rest
      nsmooth=int(ceil(mindel/arc));
      gbins=azbins/nsmooth+1;
      nsmooth=int(azbins/(gbins-1));

      //so every nsmooth steps only update the az
      //which bins to use? just start at zero and loop around.

      //gridk_bc_npole(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      
      gridk_r2_npole1_npole2_bc( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot, nsmooth, gbins);
    }else{
      //BC, j=0, l=0
      fullk_r2_npole1_npole2_bc( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot);
    }
    for(l=1;l<polAbins-1;l++){
      polA=(l+0.5)*delApol;
      sta=sin(polA);
      cta=cos(polA);
      
      //for each grid point, decide which azimuthals to include 
      arc=d1*sta*delaz+abs(r-d1)*st*delaz;
      if(arc<mindel){
	//choose a subset of grid points in k to update, interpolate the rest
	nsmooth=int(ceil(mindel/arc));
	gbins=azbins/nsmooth+1;
	nsmooth=int(azbins/(gbins-1));
     
	//cout <<" at r=: "<<r<<" theta=: "<<theta<<" too small arc: "<<arc<<" nsmooth: "<<nsmooth<<endl;
	//so every nsmooth steps only update the az
	//which bins to use? just start at zero and loop around.
	
	//gridk_bc_npole(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      
	gridk_r2_npole1_bc( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot, nsmooth, gbins);
      }else{
	//BC, j=0, l=middle
	fullk_r2_npole1_bc( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot);
      }
    }
    l=polAbins-1;
    polA=(l+0.5)*delApol;
    sta=sin(polA);
    cta=cos(polA);
    
    //for each grid point, decide which azimuthals to include 
    arc=d1*sta*delaz+abs(r-d1)*st*delaz;
    if(arc<mindel){
      //choose a subset of grid points in k to update, interpolate the rest
      nsmooth=int(ceil(mindel/arc));
      gbins=azbins/nsmooth+1;
      nsmooth=int(azbins/(gbins-1));
     
      //cout <<" at r=: "<<r<<" theta=: "<<theta<<" too small arc: "<<arc<<" nsmooth: "<<nsmooth<<endl;
      //so every nsmooth steps only update the az
      //which bins to use? just start at zero and loop around.
      
      //gridk_bc_npole(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      
      gridk_r2_npole1_spole2_bc( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot, nsmooth, gbins);
    }else{
      //BC, j=0, l=end
      fullk_r2_npole1_spole2_bc( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot);
    }
    
    

    for(j=1;j<polbins-1;j++){
      theta=(j+0.5)*delpol;//need to skip the poles.
      st=sin(theta);
      ct=cos(theta);
      l=0;
      polA=(l+0.5)*delApol;
      sta=sin(polA);
      cta=cos(polA);
      //for each grid point, decide which azimuthals to include 
      arc=d1*sta*delaz+abs(r-d1)*st*delaz;//(r-d1) is smallest size of r1.
      if(arc<mindel){
	//choose a subset of grid points in k to update, interpolate the rest
	nsmooth=int(ceil(mindel/arc));
	gbins=azbins/nsmooth+1;
	nsmooth=int(azbins/(gbins-1));
	
	//cout <<" at r=: "<<r<<" theta=: "<<theta<<" too small arc: "<<arc<<" nsmooth: "<<nsmooth<<endl;
	//so every nsmooth steps only update the az
	//which bins to use? just start at zero and loop around.
	//gridk_bc(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      
	gridk_r2_npole2_bc( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot, nsmooth, gbins);
      }else{
	//BC, j=middle, l=0
	fullk_r2_npole2_bc( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot);
      }
      //now loop over other l values
      for(l=1;l<polAbins-1;l++){
	polA=(l+0.5)*delApol;//need to skip the poles.
	sta=sin(polA);
	cta=cos(polA);
      
	//for each grid point, decide which azimuthals to include 
	arc=d1*sta*delaz+abs(r-d1)*st*delaz;//(r-d1) is smallest size of r1.
	if(arc<mindel){
	  //choose a subset of grid points in k to update, interpolate the rest
	  nsmooth=int(ceil(mindel/arc));
	  gbins=azbins/nsmooth+1;
	  nsmooth=int(azbins/(gbins-1));
	  //cout <<" at r=: "<<r<<" theta=: "<<theta<<" too small arc: "<<arc<<" nsmooth: "<<nsmooth<<endl;
	  //so every nsmooth steps only update the az
	  //which bins to use? just start at zero and loop around.
	  //gridk_bc(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      
	  gridk_r2_bc( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot, nsmooth, gbins);
	}else{
	  //BC, j=middle, l=middle
	  fullk_r2_bc( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot);
	}
      }//end over l
      //now do final l
      l=polAbins-1;
      polA=(l+0.5)*delApol;//need to skip the poles.
      sta=sin(polA);
      cta=cos(polA);
      
      //for each grid point, decide which azimuthals to include 
      arc=d1*sta*delaz+abs(r-d1)*st*delaz;//(r-d1) is smallest size of r1.
      if(arc<mindel){
	//choose a subset of grid points in k to update, interpolate the rest
	nsmooth=int(ceil(mindel/arc));
	gbins=azbins/nsmooth+1;
	nsmooth=int(azbins/(gbins-1));
	//cout <<" at r=: "<<r<<" theta=: "<<theta<<" too small arc: "<<arc<<" nsmooth: "<<nsmooth<<endl;
	//so every nsmooth steps only update the az
	//which bins to use? just start at zero and loop around.
	//gridk_bc(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      
	gridk_r2_spole2_bc( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot, nsmooth, gbins);
      }else{
	//BC, j=middle, l=polAbins
	fullk_r2_spole2_bc( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot);
      }
    }
    j=polbins-1;
    theta=(j+0.5)*delpol;//need to skip the poles.
    st=sin(theta);
    ct=cos(theta);
    //for each grid point, decide which azimuthals to include 
    l=0;
    polA=(l+0.5)*delApol;//need to skip the poles.
    sta=sin(polA);
    cta=cos(polA);
    
    arc=d1*sta*delaz+abs(r-d1)*st*delaz;//(r-d1) is smallest size of r1.
    
    if(arc<mindel){
      //choose a subset of grid points in k to update, interpolate the rest
      nsmooth=int(ceil(mindel/arc));
      gbins=azbins/nsmooth+1;
      nsmooth=int(azbins/(gbins-1));
	//cout <<" at r=: "<<r<<" theta=: "<<theta<<" too small arc: "<<arc<<" nsmooth: "<<nsmooth<<endl;
      //so every nsmooth steps only update the az
      //which bins to use? just start at zero and loop around.
      //gridk_bc_spole(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      
      gridk_r2_spole1_bc( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot, nsmooth, gbins);
    }else{
      fullk_r2_spole1_npole2_bc( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot);
    }
    //now loop over other l's
    for(l=1;l<polAbins-1;l++){
      polA=(l+0.5)*delApol;//need to skip the poles.
      sta=sin(polA);
      cta=cos(polA);
      
      arc=d1*sta*delaz+abs(r-d1)*st*delaz;//(r-d1) is smallest size of r1.
      
      if(arc<mindel){
	//choose a subset of grid points in k to update, interpolate the rest
	nsmooth=int(ceil(mindel/arc));
	gbins=azbins/nsmooth+1;
	nsmooth=int(azbins/(gbins-1));
	//cout <<" at r=: "<<r<<" theta=: "<<theta<<" too small arc: "<<arc<<" nsmooth: "<<nsmooth<<endl;
	//so every nsmooth steps only update the az
	//which bins to use? just start at zero and loop around.
	//gridk_bc_spole(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      
	gridk_r2_spole1_bc( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot, nsmooth, gbins);
      }else{
	fullk_r2_spole1_bc( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot);
      }
    
    }
    l=polAbins-1;
    polA=(l+0.5)*delApol;//need to skip the poles.
    sta=sin(polA);
    cta=cos(polA);
    
    arc=d1*sta*delaz+abs(r-d1)*st*delaz;//(r-d1) is smallest size of r1.
    
    if(arc<mindel){
      //choose a subset of grid points in k to update, interpolate the rest
      nsmooth=int(ceil(mindel/arc));
      gbins=azbins/nsmooth+1;
      nsmooth=int(azbins/(gbins-1));
	      //cout <<" at r=: "<<r<<" theta=: "<<theta<<" too small arc: "<<arc<<" nsmooth: "<<nsmooth<<endl;
      //so every nsmooth steps only update the az
      //which bins to use? just start at zero and loop around.
      //gridk_bc_spole(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      
      gridk_r2_spole1_spole2_bc( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot, nsmooth, gbins);
    }else{
      fullk_r2_spole1_spole2_bc( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot);
    }

    /*End over all angles at the r=sigma boundary*/
    psreal=1.0-passoc1;
    passoc[t]=passoc1;

    
    /*************All r, non BC****************/
    for(i=1;i<nrbins;i++){
      r=i*delR+sigma;
      r2=r*r;
      j=0;//BC, theta=0+delta
      theta=(j+0.5)*delpol;//need to skip the poles.
      st=sin(theta);
      ct=cos(theta);
      l=0;
      polA=(l+0.5)*delApol;
      sta=sin(polA);
      cta=cos(polA);
      
      //for each grid point, decide which azimuthals to include 
      arc=d1*sta*delaz+abs(r-d1)*st*delaz;
      if(arc<mindel){
	//choose a subset of grid points in k to update, interpolate the rest
	nsmooth=int(ceil(mindel/arc));
	gbins=azbins/nsmooth+1;
	nsmooth=int(azbins/(gbins-1));
	
	//so every nsmooth steps only update the az
	//which bins to use? just start at zero and loop around.
	
	//gridk_bc_npole(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      
	gridk_r2_npole1_npole2( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot, nsmooth, gbins);
      }else{
	//BC, j=0, l=0
	fullk_r2_npole1_npole2( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot);
      }
      for(l=1;l<polAbins-1;l++){
	polA=(l+0.5)*delApol;
	sta=sin(polA);
	cta=cos(polA);
	
	//for each grid point, decide which azimuthals to include 
	arc=d1*sta*delaz+abs(r-d1)*st*delaz;
	if(arc<mindel){
	  //choose a subset of grid points in k to update, interpolate the rest
	  nsmooth=int(ceil(mindel/arc));
	  gbins=azbins/nsmooth+1;
	  nsmooth=int(azbins/(gbins-1));
	  
	  //cout <<" at r=: "<<r<<" theta=: "<<theta<<" too small arc: "<<arc<<" nsmooth: "<<nsmooth<<endl;
	  //so every nsmooth steps only update the az
	  //which bins to use? just start at zero and loop around.
	  
	  //gridk_bc_npole(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      
	  gridk_r2_npole1( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot, nsmooth, gbins);
	}else{
	  //BC, j=0, l=middle
	  fullk_r2_npole1( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot);
	}
      }
      l=polAbins-1;
      polA=(l+0.5)*delApol;
      sta=sin(polA);
      cta=cos(polA);
      
      //for each grid point, decide which azimuthals to include 
      arc=d1*sta*delaz+abs(r-d1)*st*delaz;
      if(arc<mindel){
	//choose a subset of grid points in k to update, interpolate the rest
	nsmooth=int(ceil(mindel/arc));
	gbins=azbins/nsmooth+1;
	nsmooth=int(azbins/(gbins-1));
	
	//cout <<" at r=: "<<r<<" theta=: "<<theta<<" too small arc: "<<arc<<" nsmooth: "<<nsmooth<<endl;
	//so every nsmooth steps only update the az
	//which bins to use? just start at zero and loop around.
	
	//gridk_bc_npole(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      
	gridk_r2_npole1_spole2( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot, nsmooth, gbins);
      }else{
	//BC, j=0, l=end
	fullk_r2_npole1_spole2( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot);
      }
      
      
      
      for(j=1;j<polbins-1;j++){
	theta=(j+0.5)*delpol;//need to skip the poles.
	st=sin(theta);
	ct=cos(theta);
	l=0;
	polA=(l+0.5)*delApol;
	sta=sin(polA);
	cta=cos(polA);
	//for each grid point, decide which azimuthals to include 
	arc=d1*sta*delaz+abs(r-d1)*st*delaz;//(r-d1) is smallest size of r1.
	if(arc<mindel){
	  //choose a subset of grid points in k to update, interpolate the rest
	  nsmooth=int(ceil(mindel/arc));
	  gbins=azbins/nsmooth+1;
	  nsmooth=int(azbins/(gbins-1));
	  
	  //cout <<" at r=: "<<r<<" theta=: "<<theta<<" too small arc: "<<arc<<" nsmooth: "<<nsmooth<<endl;
	  //so every nsmooth steps only update the az
	  //which bins to use? just start at zero and loop around.
	  //gridk(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      
	  gridk_r2_npole2( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot, nsmooth, gbins);
	}else{
	  //BC, j=middle, l=0
	  fullk_r2_npole2( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot);
	}
	//now loop over other l values
	for(l=1;l<polAbins-1;l++){
	  polA=(l+0.5)*delApol;//need to skip the poles.
	  sta=sin(polA);
	  cta=cos(polA);
	  
	  //for each grid point, decide which azimuthals to include 
	  arc=d1*sta*delaz+abs(r-d1)*st*delaz;//(r-d1) is smallest size of r1.
	  if(arc<mindel){
	    //choose a subset of grid points in k to update, interpolate the rest
	    nsmooth=int(ceil(mindel/arc));
	    gbins=azbins/nsmooth+1;
	    nsmooth=int(azbins/(gbins-1));
	    //cout <<" at r=: "<<r<<" theta=: "<<theta<<" too small arc: "<<arc<<" nsmooth: "<<nsmooth<<endl;
	    //so every nsmooth steps only update the az
	    //which bins to use? just start at zero and loop around.
	    //gridk(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      
	    gridk_r2( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot, nsmooth, gbins);
	  }else{
	    //BC, j=middle, l=middle
	    fullk_r2( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot);
	  }
	}//end over l
	//now do final l
	l=polAbins-1;
	polA=(l+0.5)*delApol;//need to skip the poles.
	sta=sin(polA);
	cta=cos(polA);
	
	//for each grid point, decide which azimuthals to include 
	arc=d1*sta*delaz+abs(r-d1)*st*delaz;//(r-d1) is smallest size of r1.
	if(arc<mindel){
	  //choose a subset of grid points in k to update, interpolate the rest
	  nsmooth=int(ceil(mindel/arc));
	  gbins=azbins/nsmooth+1;
	  nsmooth=int(azbins/(gbins-1));
	  //cout <<" at r=: "<<r<<" theta=: "<<theta<<" too small arc: "<<arc<<" nsmooth: "<<nsmooth<<endl;
	  //so every nsmooth steps only update the az
	  //which bins to use? just start at zero and loop around.
	  //gridk(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      
	  gridk_r2_spole2( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot, nsmooth, gbins);
	}else{
	  //BC, j=middle, l=polAbins
	  fullk_r2_spole2( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot);
	}
      }
      j=polbins-1;
      theta=(j+0.5)*delpol;//need to skip the poles.
      st=sin(theta);
      ct=cos(theta);
      //for each grid point, decide which azimuthals to include 
      l=0;
      polA=(l+0.5)*delApol;//need to skip the poles.
      sta=sin(polA);
      cta=cos(polA);
      
      arc=d1*sta*delaz+abs(r-d1)*st*delaz;//(r-d1) is smallest size of r1.
      
      if(arc<mindel){
	//choose a subset of grid points in k to update, interpolate the rest
	nsmooth=int(ceil(mindel/arc));
	gbins=azbins/nsmooth+1;
	nsmooth=int(azbins/(gbins-1));
	//cout <<" at r=: "<<r<<" theta=: "<<theta<<" too small arc: "<<arc<<" nsmooth: "<<nsmooth<<endl;
	//so every nsmooth steps only update the az
	//which bins to use? just start at zero and loop around.
	//gridk_bc_spole(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      
	gridk_r2_spole1( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot, nsmooth, gbins);
      }else{
	fullk_r2_spole1_npole2( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot);
      }
      //now loop over other l's
      for(l=1;l<polAbins-1;l++){
	polA=(l+0.5)*delApol;//need to skip the poles.
	sta=sin(polA);
	cta=cos(polA);
	
	arc=d1*sta*delaz+abs(r-d1)*st*delaz;//(r-d1) is smallest size of r1.
	
	if(arc<mindel){
	  //choose a subset of grid points in k to update, interpolate the rest
	  nsmooth=int(ceil(mindel/arc));
	  gbins=azbins/nsmooth+1;
	  nsmooth=int(azbins/(gbins-1));
	  //cout <<" at r=: "<<r<<" theta=: "<<theta<<" too small arc: "<<arc<<" nsmooth: "<<nsmooth<<endl;
	  //so every nsmooth steps only update the az
	  //which bins to use? just start at zero and loop around.
	  //gridk_bc_spole(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      
	  gridk_r2_spole1( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot, nsmooth, gbins);
	}else{
	  fullk_r2_spole1( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot);
	}
	
      }
      l=polAbins-1;
      polA=(l+0.5)*delApol;//need to skip the poles.
      sta=sin(polA);
      cta=cos(polA);
      
      arc=d1*sta*delaz+abs(r-d1)*st*delaz;//(r-d1) is smallest size of r1.
      
      if(arc<mindel){
	//choose a subset of grid points in k to update, interpolate the rest
	nsmooth=int(ceil(mindel/arc));
	gbins=azbins/nsmooth+1;
	nsmooth=int(azbins/(gbins-1));
	//cout <<" at r=: "<<r<<" theta=: "<<theta<<" too small arc: "<<arc<<" nsmooth: "<<nsmooth<<endl;
	//so every nsmooth steps only update the az
	//which bins to use? just start at zero and loop around.
	//gridk_bc_spole(gbins, nsmooth, i,  j,  r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2);      
	gridk_r2_spole1_spole2( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot, nsmooth, gbins);
      }else{
	fullk_r2_spole1_spole2( i,  j,  l, r, r2,  theta, st,  BC,  delR,  deltat, Dtot, dr2, ka,  delpol, delaz, azbins, passoc1, psurvive, prt,  prt2, sta, delApol, ct, cta, d1, Drot);
      }
      
      
      
      
    }//end R
    /*renormalize*/
    cout <<"time: "<<t<<" dt: "<<time<<" normalization: "<<psurvive<<" passoc: "<<passoc[t]<<" 1-passoc: "<<psreal<<endl;
    /*for each time step, integrate over all space to get survival prob*/
    for(i=0;i<nrbins;i++){
      for(j=0;j<polbins;j++){
	for(l=0;l<polAbins;l++){
	  for(k=0;k<azbins;k++){
	    prt[i][j][l][k]=prt2[i][j][l][k]/psurvive*psreal;
	  }
	}
      }
    }
    
    if(t%twrite==0){
      //free_prop_sphere(sigma, nrbins, polbins, azbins, delR, delpol, delaz, Dtot,  time,  r0,  az0,pol0, pfree, norm);  
      //cout <<"free prop normalization at time: "<<time<<" norm: "<<norm<<endl;
      /*for each time step, integrate over all space to get survival prob*/
      for(i=0;i<nrbins;i++){
	for(j=0;j<polbins;j++){
	  for(k=0;k<azbins;k++){
	    outfile<<prt[i][j][k]<<'\t';
	    //ffile<<pfree[i][j][k]/norm<<'\t';
	  }
	}
	outfile<<endl;
	//ffile<<endl;
      }
      outfile.close();
      //ffile.close();
      
    }
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
void get_initial_config(double ***pfree, double **protfree, double ****prt, int nrbins, double delR, double sigma, double d1, int polbins, double delpol, int polAbins, double delApol, int azbins, double delaz, double rlim)
{
  int i, j, k, p, q;
  double jac, b;
  double r, pol, polA, az1, az2;
  double st, sta, ct, cta;
  double omega;
  int lbin, rbin;
  double r2;
  int fullazbins=(azbins-1)*2;
  int l;
  double c, sfact, r1mag, r1sq;
  for(i=0;i<nrbins;i++){
    for(j=0;j<polbins;j++){
      for(k=0;k<polAbins;k++){
	for(l=0;l<azbins;l++){
	  prt[i][j][k][l]=0.0;
	}
      }
    }
  }

  /*loop over initial coordinates to get prt(r, pol, polA, w)*/
  for(i=0;i<nrbins;i++){
    r=i*delR+sigma;
    if(r>rlim){
      cout <<"sum at r1=: "<<r<<endl;
      for(j=0;j<polbins;j++){
	pol=(j+0.5)*delpol;
	st=sin(pol);
	ct=cos(pol);
	for(k=0;k<fullazbins;k++){
	  az1=k*delaz;
	  if(pfree[i][j][k]>0){
	    for(p=0;p<polAbins;p++){
	      polA=(p+0.5)*delApol;
	      st=sin(polA);
	      ct=cos(polA);
	      
	      for(q=0;q<fullazbins;q++){
		az2=q*delaz;
		omega=abs(az1-az2);
		b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
		r2=sqrt(r*r+d1*d1+r*b);//*d1*(st*sta*cos(omega)+ct*cta));
		jac=r2/(r+b/2.0);
		rbin=int((r2-sigma)/delR);
		lbin=int(omega/delaz);
		if(jac==0)cout <<"divide by zero at r1=: "<<r<<" r2: "<<r2<<endl;
		if(isnan(jac))cout <<"jac is NAN: "<<r<<" r2: "<<r2<<endl;
		prt[rbin][j][k][lbin]+=pfree[i][j][k]/jac*protfree[p][q]*delaz;//assume delaz and delomega are the same
	      }
	    }//end sphere
	  }//skip sphere if pfree (translate) is zero.
	}
      }
    }//don't include r1 values that would push r2 into the boundary
  }
  cout <<"Now measure survival prob: "<<endl;
  double psurvive=0.0;

  /*count normalization*/
  for(i=0;i<nrbins;i++){
    r=delR*i+sigma;
    for(j=0;j<polbins;j++){
      pol=(j+0.5)*delpol;
      st=sin(pol);
      ct=cos(pol);
      for(k=0;k<polAbins;k++){
	polA=(k+0.5)*delApol;
	sta=sin(polA);
	cta=cos(polA);
	for(l=0;l<azbins;l++){
	  omega=l*delaz;
	  b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
	  c=-r*r+d1*d1;//r is r2!
	  sfact=sqrt(b*b-4.0*c);
	  r1mag=-b/2.0+sfact/2.0;
	  r1sq=r1mag*r1mag;
	  jac=r/(r1mag+b/2.0);
	  if(jac==0)cout <<"divide by zero at r1=: "<<r1mag<<" r2: "<<r<<endl;
	  if(isnan(jac))cout <<"jac is NAN: "<<r1mag<<" r2: "<<r<<endl;
		
	  psurvive+=r1sq*delR*prt[i][j][k][l]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
	}
      }
    }
  }
  cout <<"Psurvive of initial configuration?: "<<psurvive<<endl;
  
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

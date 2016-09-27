/*


particle-based Reaction diffusion
algorithm with trajectory reweighting.

this program only propagates a single pair of particles,
collects statistics for different starting separations over
many repeats. 

radiation BC

both particles can rotate.
  
  
*/
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "rand_gsl.h"
#include "md_timer.h"
#include "vector_rot_calls.h"

#define MAXIFACE 20
#define MAXPRTNER 20
#define MAXCOMPLEX 20
#define MAXRXN 20
#define MAXOVERLAP 20

using namespace std;

struct MD_Timer totaltime;
struct MD_Timer bimoltime;


class Protein
{
public:
  int ninterface;
  int valiface[MAXIFACE];
  int npropart;
  int propart[MAXPRTNER];
  double Dx;
  double Dy;
  double Dz;
  double radx;
  double rady;
  double radz;
  double Drx;
  double Dry;
  double Drz;
  
  int nint_write;
  int wrlist[MAXIFACE];
};
class Fullmol
{
public:
  double mytime;
  int mybin;
  int mybinind;
  int protype;
  int ninterface;
  int istatus[MAXIFACE];
  int npartner;
  int mycomplex;
  int matrix[MAXIFACE];
  double xcom;
  double ycom;
  double zcom;
  double x[MAXIFACE];
  double y[MAXIFACE];
  double z[MAXIFACE];
  int nbnd;
  int nfree;
  int freelist[MAXIFACE];
  int bndlist[MAXIFACE];
  double massx;
  double massy;
  double massz;
  int partner[MAXIFACE];
  double Dx;
  double Dy;
  double Dz;
  int npropart;
  int propart[MAXPRTNER];
  double radx;
  double rady;
  double radz;
};
class Parms
{
public:
  double dihed0;
  int Nx0;
  int Nthet;
  int thetbins;
  int dihbins;
  int tabbins;
  double x0;
  int Nprotypes;
  int Nifaces;
  double Nit;
  double leglen;
  double eps_scale;
  double dt_scale;
  double rate;
  double Maxtime;
  int restart;
  int statwrite;
  int configwrite;
  int Nrep;
  double deltat;
  double V;
  double X0total;
  int Nspecies;
  int Nrxn;
  double mass;
  double D;
  int nspec_complex;
  double maxsep2;
  int ntotalcomplex;
  double xboxl;
  double yboxl;
  double zboxl;
  int Ntotalmol;
  int Natom;
  int Natomwrite;
};



void read_parms(ifstream &parmfile, Parms &plist);
double GaussV();
double pirr_pfree_ratio_ps(double rcurr, double r0, double tcurr, double Dtot, double bindrad, double alpha, double ps_prev, double rtol);
double survive_irr(double r0, double tcurr, double Dtot, double bindrad, double alpha, double cof);
double pirrev_value(double rcurr, double r0, double tcurr, double Dtot, double bindrad, double alpha);
double pfree_value_norm(double rcurr, double r0, double tcurr, double Dtot, double bindrad,double alpha);


int main(int argc, char *argv[])
{
  
  int i, j, k;
  timeval tim;
  gettimeofday(&tim, 0);
  double t1=tim.tv_sec+tim.tv_usec;
  
  int seed=int(t1);
  
  //seed=1353432282;
  double randmax=pow(2.0, 32);
  cout <<"seed: "<<seed<<" randmax: "<<randmax<<endl;
  srand_gsl(seed);
  double irandmax=1.0/randmax;
  
  ifstream parmfile(argv[1]);
  Parms plist;
  read_parms(parmfile, plist);
  //  write_parms(plist);
  cout.precision(12); 
  initialize_timer(&totaltime);
  initialize_timer(&bimoltime);
  start_timer(&totaltime);
  
  int Nprotypes=plist.Nprotypes;//total distinct protein types, 
  int Nifaces=plist.Nifaces;//this is the number of interfaces

  int Nrxn=plist.Nrxn;
  int Nspecies=plist.Nspecies;//this will include product species



  double cellvol=plist.xboxl*plist.yboxl*plist.zboxl/(1.0*1E9);//sidelength*sidelength*sidelength;
  double V=cellvol; 
  double um_to_L=1E15;
  double avagad=6.022E23;


  int *Ncopy=new int[Nprotypes];

  Ncopy[0]=1;
  Ncopy[1]=1;
  int Ntotalmol=0;
  plist.Natom=0;
  int ntmp;
  for(i=0;i<Nprotypes;i++){

    Ntotalmol+=Ncopy[i];
  }
  cout <<"Ntotal mols: "<<Ntotalmol<<endl;//ASSUMED TO BE 2 MOLECULES BELOW
  plist.Ntotalmol=Ntotalmol;
  Fullmol *bases=new Fullmol[Ntotalmol];//contains information on each protein in the full system
  
  int *numpartners=new int[Nifaces];//this should account for all free interfaces
  int **Speclist=new int*[Nifaces];
  
  for(i=0;i<Nifaces;i++)
    Speclist[i]=new int[MAXPRTNER];
  

  /*A interacts with B*/  
  numpartners[0]=1;
  numpartners[1]=1;
  Speclist[0][0]=1;
  Speclist[1][0]=0;

  char fname[100];
  
  Protein *wholep=new Protein[Nprotypes];
  int *p_home=new int[Nifaces];//this reverses and tells you what protein a given interface belongs to
  int *i_home=new int[Nspecies];//for both free and bound states, what index are you on the protein's list
  double *bindrad=new double[Nrxn];//binding or unbinding radius for each reaction
  int *Ncoup=new int[Nrxn];//list of reactions coupled to this one
  int **mycoupled=new int*[Nrxn];
  for(i=0;i<Nrxn;i++)
    mycoupled[i]=new int[MAXRXN];
  /*The number of reactions is fixed and all the same reactions are possible in each spatial cell*/
  double *kr=new double[Nrxn]; //reaction rate (with dimensions)

  int **Rlist=new int*[Nrxn]; //the identity of the species in the reaction
  int *Npart=new int[Nrxn]; //The number of participant species in a reaction
  int **Del=new int*[Nrxn]; //The coeffiecients of the participants in the reaction
  int *rxntype=new int[Nrxn];//this sets whether the reaction is bimolecular, unimolecular.
  int maxrctant=5;
  for(i=0;i<Nrxn;i++){
    Rlist[i]=new int[maxrctant];
    Del[i]=new int[maxrctant];
  }
  
  /*Each protein has a single site, each protein can have it's own Diffusion constant*/
  
  wholep[0].ninterface=1;
  wholep[1].ninterface=1;
  
  bindrad[0]=1;
  kr[0]=plist.rate;
  cout <<"ACTIVATION RATE: "<<kr[0]<<"  radius: "<<bindrad[0]<<endl;
  double kact=kr[0];

  for(i=0;i<Nprotypes;i++){  
    ntmp=wholep[i].ninterface+1;
    plist.Natom+=Ncopy[i]*ntmp;
  }
  cout <<"N atoms: "<<plist.Natom<<endl;
  cout <<"read reactions "<<endl;

  
  double *savecrds=new double[Ntotalmol*3];//for x, y, z
  
  /*Print out specific reactions*/
  cout <<"Print specific interaction network "<<endl;
  int ncomplex=0;
  for(i=0;i<Nifaces;i++){
    cout <<i<<'\t';
    for(j=0;j<numpartners[i];j++){
      cout <<Speclist[i][j]<<'\t';
      ncomplex++;
    }
    cout <<endl;
  }
  ncomplex/=2;
  plist.nspec_complex=ncomplex;
  

  
  int ind, r1, m;
  int begin, end;
  


  double rnum;
  
  int nfaces=6;//for a cubic volume
  int direction, neighbor;
  
  int rxn;
  
  double curr_time=0;

  int checkpoint=10000000; /*How often to write out full current solution*/
  int stepwrite=100; /*How often to write out species numbers*/
  char fnamemid[100];

  
  
  double h;
  /*Iterate over time steps until you hit a max time*/
  int t, mu;
  int flag;
  double prob;
  double sum;
  double xchg, ychg, zchg;
  int icom;
  int pro_type, mp;
  int whichspecie;
  double rerand;
  double hfact;
  double dx, dy, dz; 
  double box_x=plist.xboxl;
  double box_y=plist.yboxl;
  double box_z=plist.zboxl;//nm
  double xtot, ytot, ztot;
  int nfree, wprot;
  int p, i1, i2;
  int np;
  double r2, r;
  double maxsep2=bindrad[0]*bindrad[0];//plist.maxsep2;
  cout <<"squared distance cutoff: "<<maxsep2<<endl;
  int iind, iind2, ppart;
  int twrite=plist.configwrite;



  
  int it;
  plist.ntotalcomplex=Ntotalmol;

  
  int s1;
  cout <<"Ntotal complexes: "<<plist.ntotalcomplex<<endl;

  int amol,df; 
  double us_to_s=1E-6;
  int statwrite=plist.statwrite;

  double kpi=4.0*M_PI*bindrad[0];
  double leglen=plist.leglen;
  double leglen2=leglen*leglen;
  cout <<"Leg length: "<<leglen<<endl;
  double T=293;//K
  double nu=0.001;//kg/(m*s)
  double scale=3.0;//greater the one to correct for non-spherical
  double x2avg;
  double xavg;
  
  double kb=1.3806488E-23;
  double preT=kb*T/(6.0*M_PI*nu)*1E27/1E6;//1E27 is nm^3 and 1E6 is us
  double crad=preT/plist.D*2.0;//cut D in half because it is Da+Db
  double preR=kb*T/(8.0*M_PI*nu)*1E27/1E6;//1E27 is nm^3 and 1E6 is us
  wholep[1].Drx=preR/(crad*crad*crad);
  wholep[1].Dry=wholep[1].Drx;
  wholep[1].Drz=wholep[1].Drx;

  wholep[1].Dx=plist.D;
  wholep[1].Dy=wholep[1].Dx;
  wholep[1].Dz=wholep[1].Dx;
  cout <<"D: "<<plist.D <<" effective radius: "<<crad<<" D, calc: "<<wholep[1].Dx<<" Drot, calc: "<<wholep[1].Drx<<endl;
  wholep[0].Dx=0;
  wholep[0].Dy=0;
  wholep[0].Dz=0;
  wholep[0].Drx=wholep[1].Drx;
  wholep[0].Dry=wholep[1].Dry;
  wholep[0].Drz=wholep[1].Drz;
  
  double r1x, r1y, r1z;
  double cthet1;

  int ind_thet;
 
  
  
  double Dtot=wholep[0].Dx+wholep[1].Dx;
  double fourpi=4.0*M_PI;
  
  double R2, R1;
  
  double r0, passoc;

  double fact;
  double kdiff;//will be kpi*D
  double aexp;
  double bexp;
  double alpha;
  
  mu=0;
  double Rmax;
  i=0;
  
  r0=bindrad[0];
  int r0bins=1000;
  double delr0;
  int r0ind;
  
  
  double p0_ratio=1.0;
  //int distflag=0;
  double currsep;
  double prevpassoc=0;
  double probvec1;
  //cout <<"Rmax1: "<<Rmax1
  cout <<" Dtot: "<<Dtot<<endl;

  double rnum2;
  double probA=0;
  int flagA=0;
  double pact;
  double tmpx, tmpy, tmpz;
  int p1, p2;
  int c1, c2;
  int ci1, ci2;
  double tremain, tevent;
  int mu_ret;
  int rxn1;
  int go;
  double rate;
  int cancel;

  double maxrad=20;
  int radbins=1000;
  double delrad=maxrad/(1.0*radbins);
  int ind_rad;
  double rad2, rad;

  int nc1, nc2;
  double sep, ratio;
  int flagsep;
  it=1;
  int rep;
  int Nrep=plist.Nrep;
  double realsmall=1E-14;
  double prevnorm=1.0;
  int Nx0=plist.Nx0;
  double dels;
  
  ofstream probfile;
  double space;
  double tval;
  char tname[200];
  double dx0, dy0, dz0;
  double rtol=1E-10;
  int previter;
  
  double currnorm;
  double newprob;
  double prevsep;
  double ps_prev;
  int s;
  double x0;
  double currx, curry, currz;
  double Nhist;
  double *x0vec=new double[Nx0];
  cout <<"Initial separation: "<<plist.x0<<endl;
  x0vec[0]=plist.x0;//for k<inf, this will not go to 1 
  for(i=1;i<Nx0;i++){
    x0vec[i]=plist.x0+0.05*i;
    cout <<"X0, i: "<<i<<" x0: "<<x0vec[i]<<endl;
  }
  
  double kappa=kact/(4.0*M_PI*bindrad[0]*bindrad[0]);
  cout <<"Kact: "<<kact<<" kappa: "<<kappa<<endl;
  double epsilon=plist.eps_scale*Dtot/kappa;
  double tau=epsilon/kappa;
  cout <<"epsscale: "<<plist.eps_scale<<" epsilon: "<<epsilon<<" tau: "<<tau<<endl;
  Rmax=bindrad[0]+epsilon;
  cout <<"Bindrad: "<<bindrad[0]<<" reaction limit: "<<Rmax<<endl;
  cout <<"different time step expansions: "<<epsilon*epsilon/Dtot/100.0<<" other: "<<Rmax*Rmax*0.0001/(2.0*Dtot)<<endl;
  double scaled=plist.dt_scale;
  double deltat_reac=scaled*epsilon*epsilon/(2.0*Dtot);

  double cf=cos(sqrt(4.0*wholep[1].Drx*deltat_reac));
  double Dr1=2.0*leglen2*(1.0-cf);

  Dtot+=2.0*Dr1/(6.0*deltat_reac);//add in rotation for both molecules
  cout <<"add to Dtot from Rotation: "<<2.0*Dr1/(6.0*deltat_reac)<<" original dt: "<<deltat_reac<<" final Dtot: "<<Dtot<<endl;
  /*Now update to new Dtot*/
  epsilon=plist.eps_scale*Dtot/kappa;
  tau=epsilon/kappa;
  Rmax=bindrad[0]+epsilon;
  deltat_reac=scaled*epsilon*epsilon/(2.0*Dtot);
  double deltat=deltat_reac;
  
  kdiff=fourpi*Dtot*bindrad[0];
  fact=1.0+kact/kdiff;
  alpha=fact*sqrt(Dtot)/bindrad[0];
  double cof=kact/(kact+kdiff);
  
  cout <<" deltat: "<<deltat_reac<<" epsilon: "<<epsilon<<" tau: "<<tau<<" Rmax: "<<Rmax<<endl;
  double currtime=0;
  double Maxtime=plist.Maxtime;
  int Nitbin=int(Maxtime/deltat_reac);
  Maxtime=Nitbin*deltat_reac;
  cout <<"number of time bins: "<<Nitbin<<" new max time: "<<Maxtime<<endl;
  double pirrev, pfree;
  //current statwrite is the number of datapoints it will write out.
  if(Nitbin<statwrite)
    statwrite=500;
  else
    statwrite=int(round(Nitbin/statwrite));
  
  cout <<"statwrite: "<<statwrite<<endl;
  
  double theta;
  int Nthet=plist.Nthet;
  double deltheta=2.0/(1.0*(Nthet-1));
  cout <<"Ntheta: "<<Nthet<<" delthet: "<<deltheta<<endl;
  double *costheta=new double[Nthet];
  for(i=0;i<Nthet;i++){
    costheta[i]=1-i*deltheta;
    cout <<"costheta: "<<costheta[i]<<endl;
  }
  
  //int Nwrite=int(ceil(Nitbin/statwrite));
  double *phist=new double[Nitbin];
  for(i=0;i<Nitbin;i++)
    phist[i]=0.0;
  
  
  int cnt;
  
  int n;
  int itmove;
  int itmax;
  double lz, ly, lx;
  double tx, ty, tz;
  double Rrange=2.5*sqrt(6.0*Dtot*Maxtime);
  cout <<"MaxR: "<<Rrange<<endl;
  
  double delR=0.01;
  int Rbins=int(Rrange/delR);
  int rind;
  int thetbins=plist.thetbins;
  int thetbins2=thetbins*thetbins;
  double delthet1=2.0/(1.0*thetbins);
  cout <<"theta bins: "<<thetbins<<" delthet1: "<<delthet1<<endl;
  ofstream thetfile("thetas.out");
  for(i=0;i<thetbins;i++)
    thetfile <<(i+0.5)*delthet1-1<<endl;
  
  int dihbins=plist.dihbins;
  double deldih=M_PI/(1.0*dihbins);
  cout <<"dihedral bins: "<<dihbins<< " delta: "<<deldih<<endl;
  
  double ***prhist=new double**[thetbins2];//
  for(i=0;i<thetbins2;i++)
    prhist[i]=new double*[dihbins];
  for(i=0;i<thetbins2;i++){
    for(k=0;k<dihbins;k++){
      prhist[i][k]=new double[Rbins];
    }
  }
  for(i=0;i<thetbins2;i++){
    for(k=0;k<dihbins;k++){
      for(j=0;j<Rbins;j++)
	prhist[i][k][j]=0.0;
    }
  }
  double *p1hist=new double[Rbins];
  for(i=0;i<Rbins;i++)
    p1hist[i]=0.0;

  //int tabbins=plist.tabbins;//number of bins for angle to orientation
  // double deltab=2.0/(1.0*tabbins);//another indexed in cos(theta).
  // double **prhist2=new double*[tabbins*tabbins];
//   for(i=0;i<tabbins*tabbins;i++)
//     prhist2[i]=new double[Rbins];
//   for(i=0;i<tabbins*tabbins;i++){
//     for(j=0;j<Rbins;j++)
//       prhist2[i][j]=0.0;
//   }
  double *M=new double[9];
  double *v3=new double[3];
  double *v2=new double[3];
  double *v1=new double[3];
  double *v=new double[3];
  ofstream rfile;
  //ofstream tabfile;
  //declare more vars
  double rabx, raby, rabz;
  double rab2, rab;
  double r0x, r0y, r0z;
  double tab0, tab1;
  int tab0ind, tab1ind;
  double cthet2;
  int ind_thet2;
  double *n1=new double[3];
  double *n2=new double[3];
  double *sumr=new double[Rbins];
  double *sumr2=new double[Rbins];
	    
  double cosdih;
  double dih;
  int dihbin;
  //ofstream matchfile;
  int wcnt;
  int ct;
  int bt;
  double lx_0, ly_0, lz_0;
  double ctab_0, ctab1;
  double lR;
  double imx, imy, imz;
  int cts=0;
  int ThetEnd=Nthet;
  double dihed=plist.dihed0;
  double dihlimit;//close to pi, but is more lenient for oriented legs.
  double anglelimit=0.1;
  if(dihed>0){
    cts=1;//skip 1 again.
    ThetEnd=Nthet-1;//also -1 has only a single dihedral.
  }
  for(s=0;s<Nx0;s++){
    x0=x0vec[s];
    delR=(Rrange+x0-bindrad[0])/(1.0*Rbins);
    for(ct=cts;ct<ThetEnd;ct++){
      lz=costheta[ct]*leglen;
      lx=sqrt(leglen2-lz*lz);
      ly=0;
      cout <<"Initial pos for leg1 com: "<<lx<<' '<<ly<<' '<<x0+lz<<"  for leg: "<<0<<' '<<0<<' '<<x0<<endl;
      
      
      for(bt=ct;bt<ThetEnd;bt++){
	
	lz_0=-costheta[bt]*leglen;
	lR=sqrt(leglen2-lz_0*lz_0);
	lx_0=lR*cos(dihed);
	ly_0=lR*sin(dihed);
	cout <<"Initial pos for leg2 com: "<<lx_0<<' '<<ly_0<<' '<<lz_0<<"  for leg: "<<0<<' '<<0<<' '<<0<<endl;
	//vector to rotator centers
	raby=ly-ly_0;
	rabx=lx-lx_0;
	rabz=lz+x0-lz_0;
	rab2=rabx*rabx+raby*raby+rabz*rabz;
	rab=sqrt(rab2);
	//negative is for zero-the center
	ctab_0=-lz_0*rabz-ly_0*raby-lx_0*rabx;
	ctab_0/=(rab*leglen);//costheta, for orientation, theta->0, cos(theta)->1
	ctab1=lz*rabz+lx*rabx+ly*raby;
	ctab1/=(rab*leglen);
	cout <<"Initial leg1 angle to the com divider: "<<ctab1<<" initial leg2 angle to the com divider: "<<ctab_0<<endl;
	for(i=0;i<Rbins;i++){
	  p1hist[i]=0.0;
	  sumr[i]=0.0;
	  sumr2[i]=0.0;
	}
	for(i=0;i<thetbins2;i++){
	  for(k=0;k<dihbins;k++){
	    for(j=0;j<Rbins;j++)
	      prhist[i][k][j]=0.0;
	  }
	}
	// for(i=0;i<tabbins*tabbins;i++){
// 	  for(j=0;j<Rbins;j++)
// 	    prhist2[i][j]=0.0;
// 	}
	
	for(i=0;i<Nitbin;i++)
	  phist[i]=0.0;
	
	cout <<"X0: "<<x0<<" costheta: "<<costheta[ct]<<" costheta2: "<<costheta[bt]<<endl;
      /*Need to sample over different values of the sphere holding
	all the points at leg length. origin of sphere is at z=x0, x=y=0,
	azimuth doesn't matter, choose phi=0, so y=0 always to start.
	then sample evenly over cos(theta),
	x=d1*sin(theta), y=0, z=z0+d1*cos(theta).
	You will perform rotation around this point.
      */
	
	for(rep=0;rep<Nrep;rep++){
	  /*start proteins separated by x0 */
	  //cout <<"Repeat: "<<rep<<endl;
	  
	  plist.ntotalcomplex=2;
	  Ntotalmol=2;
	  bases[0].x[0]=0;
	  bases[0].y[0]=0;
	  bases[0].z[0]=0;
	  bases[0].zcom=lz_0;
	  bases[0].xcom=lx_0;
	  bases[0].ycom=ly_0;
	  //these positions are fixed
	  imx=bases[0].xcom;
	  imy=bases[0].ycom;
	  imz=bases[0].zcom;
	  
	  
	  bases[0].nfree=1;
	  bases[0].nbnd=0;
	  bases[0].npartner=0;
	  
	  
	  bases[1].x[0]=0.0;
	  bases[1].y[0]=0.0;
	  bases[1].z[0]=x0;
	  
	  bases[1].zcom=lz+x0;
	  bases[1].xcom=lx;
	  bases[1].ycom=ly;
	  
	  bases[1].nfree=1;
	  bases[1].nbnd=0;
	  bases[1].npartner=0;
	  
	  /***************************/
	  /*Begin RD simulation*/
	  it=0;
	  currtime=0;
	  //cout <<"rep: "<<rep<<" Maxtime; "<<Maxtime<<endl;
	  while(bases[0].nfree>0 &&currtime<Maxtime){
	    i=0;
	    /*Test bimolecular reactions!*/
	    nfree=bases[i].nfree;
	    j=1;
	    
	    dx=bases[1].x[0]-bases[0].x[0];
	    dy=bases[1].y[0]-bases[0].y[0];
	    dz=bases[1].z[0]-bases[0].z[0];
	    
	    R2=dx*dx+dy*dy+dz*dz;
	    R1=sqrt(R2);
	    prob=2;
	    if(R1<Rmax){
	      itmove=1;
	      /*terminate some trajectories.*/
	      deltat=deltat_reac;
	      /*Now tau is scaled by whether the orientation is correct.*/
	      
	      r0x=bases[0].x[0]-bases[0].xcom;
	      r0y=bases[0].y[0]-bases[0].ycom;
	      r0z=bases[0].z[0]-bases[0].zcom;
	      
	      r1x=bases[1].x[0]-bases[1].xcom;
	      r1y=bases[1].y[0]-bases[1].ycom;
	      r1z=bases[1].z[0]-bases[1].zcom;
	      
	      cthet1=-r1x*dx-r1y*dy-r1z*dz;
	      cthet1/=(R1*leglen);
	      
	      cthet2=r0x*dx+r0y*dy+r0z*dz;
	      cthet2/=(R1*leglen);
	      
	      	      
	      if((1-cthet1)<anglelimit && (1-cthet2)<anglelimit)
		prob=exp(-deltat_reac/tau);
	      else
		prob=2;
	      	      
	      /*might perform this reaction, depending on k_associate*/
	      rnum=1.0*rand_gsl();
	      p1=1;
	      if(prob<rnum){
		//terminate trajectory
		//cout <<"Associate, prob:  "<<rnum<<" rep: "<<rep<<" time: "<<currtime<<endl;
		itmove=0;
		p2=0;
		rxn1=0;
		plist.ntotalcomplex=1;
		bases[p1].nbnd++;
		bases[p2].nbnd++;
		bases[p1].nfree=0;
		bases[p2].nfree=0;
		bases[1].x[0]=bases[0].x[0];
		bases[1].y[0]=bases[0].y[0];
		bases[1].z[0]=bases[0].z[0];
		
		
	      }else {
		/*propagate particles and avoid overlap*/
		
		dx=sqrt(2.0*deltat*wholep[1].Dx)*GaussV();
		dy=sqrt(2.0*deltat*wholep[1].Dy)*GaussV();
		dz=sqrt(2.0*deltat*wholep[1].Dz)*GaussV();
		tx=sqrt(2.0*deltat*wholep[1].Drx)*GaussV();
		ty=sqrt(2.0*deltat*wholep[1].Dry)*GaussV();
		tz=sqrt(2.0*deltat*wholep[1].Drz)*GaussV();
		rotationEuler( tx,  ty,  tz,M);
		v[0]=bases[1].x[0]-bases[1].xcom;
		v[1]=bases[1].y[0]-bases[1].ycom;
		v[2]=bases[1].z[0]-bases[1].zcom;//pivot is the center
		rotate(v, M, v2);//includes the interface that will align
		/*now rotate particles 0*/
		tx=sqrt(2.0*deltat*wholep[0].Drx)*GaussV();
		ty=sqrt(2.0*deltat*wholep[0].Dry)*GaussV();
		tz=sqrt(2.0*deltat*wholep[0].Drz)*GaussV();
		rotationEuler( tx,  ty,  tz,M);
		v1[0]=bases[0].x[0]-bases[0].xcom;
		v1[1]=bases[0].y[0]-bases[0].ycom;
		v1[2]=bases[0].z[0]-bases[0].zcom;//pivot is the center
		rotate(v1, M, v3);//includes the interface that will align
		
		
		currx=(bases[1].xcom+dx+v2[0])-(bases[0].xcom+v3[0]);
		curry=(bases[1].ycom+dy+v2[1])-(bases[0].ycom+v3[1]);
		currz=(bases[1].zcom+dz+v2[2])-(bases[0].zcom+v3[2]);
		
		R2=currx*currx+curry*curry+currz*currz;
		wcnt=0;
		while(R2<bindrad[0]*bindrad[0]){
		  // wcnt++;
// 		  if(wcnt%1000==0){
// 		    cout <<"In while loop stil "<<endl;
// 		    cout <<"starting pos 0: "<<bases[0].xcom<<' '<<bases[0].ycom<<' '<<bases[0].zcom<<endl;
// 		    cout <<"starting pos 0leg: "<<bases[0].x[0]<<' '<<bases[0].y[0]<<' '<<bases[0].z[0]<<endl;
// 		    cout <<"starting pos 1: "<<bases[1].xcom<<' '<<bases[1].ycom<<' '<<bases[1].zcom<<endl;
// 		    cout <<"starting pos 1leg: "<<bases[1].x[0]<<' '<<bases[1].y[0]<<' '<<bases[1].z[0]<<endl;
// 		    currx=(bases[1].x[0])-(bases[0].x[0]);
// 		    curry=(bases[1].y[0])-(bases[0].y[0]);
// 		    currz=(bases[1].z[0])-(bases[0].z[0]);
		    
// 		    R2=currx*currx+curry*curry+currz*currz;
// 		    cout <<"starting sep: "<<sqrt(R2)<<endl;
		    
// 		    cout <<"displacement: "<<dx<<' '<<dy<<' '<<dz<<endl;
// 		    cout <<"rot vec 1: "<<v2[0]<<' '<<v2[1]<<' '<<v2[2]<<endl;
// 		    cout <<"rot vec 2: "<<v3[0]<<' '<<v3[1]<<' '<<v3[2]<<endl;
		    
// 		    currx=(bases[1].xcom+dx+v2[0])-(bases[0].xcom+v3[0]);
// 		    curry=(bases[1].ycom+dy+v2[1])-(bases[0].ycom+v3[1]);
// 		    currz=(bases[1].zcom+dz+v2[2])-(bases[0].zcom+v3[2]);
		    
// 		    R2=currx*currx+curry*curry+currz*currz;
// 		    cout <<"new sep: "<<sqrt(R2)<<endl;
// 		  }
		
		  dx=sqrt(2.0*deltat*wholep[1].Dx)*GaussV();
		  dy=sqrt(2.0*deltat*wholep[1].Dy)*GaussV();
		  dz=sqrt(2.0*deltat*wholep[1].Dz)*GaussV();
		  tx=sqrt(2.0*deltat*wholep[1].Drx)*GaussV();
		  ty=sqrt(2.0*deltat*wholep[1].Dry)*GaussV();
		  tz=sqrt(2.0*deltat*wholep[1].Drz)*GaussV();
		  rotationEuler( tx,  ty,  tz,M);
		  
		  rotate(v, M, v2);//includes the interface that will align
		  
		  /*now rotate particles 0*/
		  tx=sqrt(2.0*deltat*wholep[0].Drx)*GaussV();
		  ty=sqrt(2.0*deltat*wholep[0].Dry)*GaussV();
		  tz=sqrt(2.0*deltat*wholep[0].Drz)*GaussV();
		  rotationEuler( tx,  ty,  tz,M);
		  rotate(v1, M, v3);//includes the interface that will align
		  		  
		  currx=(bases[1].xcom+dx+v2[0])-(bases[0].xcom+v3[0]);
		  curry=(bases[1].ycom+dy+v2[1])-(bases[0].ycom+v3[1]);
		  currz=(bases[1].zcom+dz+v2[2])-(bases[0].zcom+v3[2]);
		  
		  R2=currx*currx+curry*curry+currz*currz;
		  
		}
		/*update particle positions that do not overlap*/
		bases[1].xcom+=dx;
		bases[1].ycom+=dy;
		bases[1].zcom+=dz;
		bases[1].x[0]=bases[1].xcom+v2[0];
		bases[1].y[0]=bases[1].ycom+v2[1];
		bases[1].z[0]=bases[1].zcom+v2[2];
		
		bases[0].x[0]=bases[0].xcom+v3[0];
		bases[0].y[0]=bases[0].ycom+v3[1];
		bases[0].z[0]=bases[0].zcom+v3[2];
		
		
	      }
	    }else{
	      /*free diffusion for both particles 
		outside of 'reaction' zone
	      */
	      //max disp=3*sqrt(6*D*dt), (R1-bindrad)/56/D
	      deltat=scaled*(R1-bindrad[0])*(R1-bindrad[0])/(Dtot*2.0);//+deltat_reac;
	      itmove=int(deltat/deltat_reac);
	      deltat=itmove*deltat_reac;
	      dx=sqrt(2.0*deltat*wholep[1].Dx)*GaussV();
	      dy=sqrt(2.0*deltat*wholep[1].Dy)*GaussV();
	      dz=sqrt(2.0*deltat*wholep[1].Dz)*GaussV();
	      tx=sqrt(2.0*deltat*wholep[1].Drx)*GaussV();
	      ty=sqrt(2.0*deltat*wholep[1].Dry)*GaussV();
	      tz=sqrt(2.0*deltat*wholep[1].Drz)*GaussV();
	      rotationEuler( tx,  ty,  tz,M);
	      v[0]=bases[1].x[0]-bases[1].xcom;
	      v[1]=bases[1].y[0]-bases[1].ycom;
	      v[2]=bases[1].z[0]-bases[1].zcom;//pivot is the center
	      rotate(v, M, v2);//includes the interface that will align
	      
	      bases[1].xcom+=dx;
	      bases[1].ycom+=dy;
	      bases[1].zcom+=dz;
	      bases[1].x[0]=bases[1].xcom+v2[0];
	      bases[1].y[0]=bases[1].ycom+v2[1];
	      bases[1].z[0]=bases[1].zcom+v2[2];
	      /*now rotate particles 0*/
	      tx=sqrt(2.0*deltat*wholep[0].Drx)*GaussV();
	      ty=sqrt(2.0*deltat*wholep[0].Dry)*GaussV();
	      tz=sqrt(2.0*deltat*wholep[0].Drz)*GaussV();
	      rotationEuler( tx,  ty,  tz,M);
	      v1[0]=bases[0].x[0]-bases[0].xcom;
	      v1[1]=bases[0].y[0]-bases[0].ycom;
	      v1[2]=bases[0].z[0]-bases[0].zcom;//pivot is the center
	      rotate(v1, M, v3);//includes the interface that will align
	      
	      bases[0].x[0]=bases[0].xcom+v3[0];
	      bases[0].y[0]=bases[0].ycom+v3[1];
	      bases[0].z[0]=bases[0].zcom+v3[2];
	      
	    }
	    
	    
	    currtime+=deltat;
	    
	    // 	  if(itmove==1){
	    // 	    phist[it+itmove]+=1;
	    // 	  }else{
	    
	    itmax=itmove+1;
	    if((itmax+it)>Nitbin)
	      itmax=Nitbin-it;
	    
	    for(n=1;n<itmax;n++){
	      phist[it+n]+=1;
	    }
	    
	    //cout <<"rep: "<<rep<<" prev it: "<<it<<" itmorve: "<<itmove<<" time: "<<currtime<<" nfree: "<<bases[0].nfree<<" final pos: "<<bases[1].xcom<<' '<<bases[1].ycom<<' '<<bases[1].zcom<<endl;
	    it+=itmove;
	    
	  }
	  if(rep%10000==0)
	    cout <<"finished rep: "<<rep<< " at time: "<<currtime<<endl;
	  /*Keep a histogram of final positions for particles that survived until the end*/
	  if(currtime>Maxtime){
	    dx=bases[1].x[0]-bases[0].x[0];
	    dy=bases[1].y[0]-bases[0].y[0];
	    dz=bases[1].z[0]-bases[0].z[0];
	    
	    R2=dx*dx+dy*dy+dz*dz;
	    R1=sqrt(R2);
	    rind=int((R1-bindrad[0])/delR);
	    rabx=bases[1].xcom-bases[0].xcom;
	    raby=bases[1].ycom-bases[0].ycom;
	    rabz=bases[1].zcom-bases[0].zcom;
	    rab2=rabx*rabx+raby*raby+rabz*rabz;
	    rab=sqrt(rab2);
	    r0x=bases[0].x[0]-bases[0].xcom;
	    r0y=bases[0].y[0]-bases[0].ycom;
	    r0z=bases[0].z[0]-bases[0].zcom;
	    //tab0=(r0x*rabx+r0y*raby+r0z*rabz)/(leglen*rab);
	    //tab0ind=int((tab0+1)/deltab);
	    r1x=bases[1].x[0]-bases[1].xcom;
	    r1y=bases[1].y[0]-bases[1].ycom;
	    r1z=bases[1].z[0]-bases[1].zcom;
	    //tab1=(-r1x*rabx-r1y*raby-r1z*rabz)/(leglen*rab);
	    //tab1ind=int((tab1+1)/deltab);
	    cthet1=-r1x*dx-r1y*dy-r1z*dz;
	    cthet1/=(R1*leglen);
	    ind_thet=int((cthet1+1)/delthet1);
	    
	    
	    cthet2=r0x*dx+r0y*dy+r0z*dz;
	    cthet2/=(R1*leglen);
	    ind_thet2=int((cthet2+1)/delthet1);
	    //calculate dihedral angle.
	    //n1=r1x (x) (-dx)
	    v[0]=-dx;
	    v[1]=-dy;
	    v[2]=-dz;
	    v1[0]=r1x;
	    v1[1]=r1y;
	    v1[2]=r1z;
	    crossproduct(v1, v, n1);
	    v2[0]=-r0x;
	    v2[1]=-r0y;
	    v2[2]=-r0z;
	    crossproduct(v, v2, n2);
	    cosdih=n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
	    dih=acos(cosdih);
	    dihbin=int(dih/deldih);
	    if(isnan(dih)){
	      
	      if(round(cosdih)==-1)dihbin=dihbins-1;
	      else dihbin=0;
	      cout <<"NAN: "<<dih<<" cosdih: "<<cosdih<<" bin: "<<dihbin<<endl;
	    }else if(dihbin==dihbins)
	      dihbin=dihbins-1;
	    
	    if(ind_thet>=thetbins){
	      // 	      cout <<"cos(theta)>1 ?: "<<cthet1<<" ind: "<<ind_thet<<endl;
 	      ind_thet=thetbins-1;
 	    }else if(ind_thet<0){
	      // 	      cout <<"cos(theta)<-1 ?: "<<cthet1<<" ind: "<<ind_thet<<endl;
 	      ind_thet=0;
 	    }
	    if(rind<Rbins){
	      prhist[ind_thet*thetbins+ind_thet2][dihbin][rind]++;
	      //prhist2[tab0ind*tabbins+tab1ind][rind]++;
	      p1hist[rind]++;
	    }
	    
	  }
	}//end over all reps
	cout <<"finished all reps "<<endl;
	/*Now write out final probability histogram*/
	sprintf(tname, "rot_sprob_x0_%g_dh%g_ct%g_c2t_%g_dt%g.dat",x0, dihed, costheta[ct], costheta[bt], deltat_reac);
	cout <<"open file: "<<tname<<endl;
	probfile.open(tname);
	probfile<<0<<' '<<1<<' '<<1<<endl;
	//afile<<x0<<'\t';
	
	for(i=1;i<Nitbin;i++){
	  tval=i*deltat_reac;
	  if(i%statwrite==0){
	    passoc=survive_irr( x0, tval,  Dtot,  bindrad[0],  alpha,  cof);
	    probfile<<tval<<' ' <<phist[i]/(1.0*Nrep)<<' '<<1.0-passoc<<endl;

	  }
	}
	i=Nitbin-1;
	tval=i*deltat_reac;
	passoc=survive_irr( x0, tval,  Dtot,  bindrad[0],  alpha,  cof);
	probfile<<tval<<' ' <<phist[i]/(1.0*Nrep)<<' '<<1.0-passoc<<endl;
	probfile.close();
	
	passoc=survive_irr( x0, Maxtime,  Dtot,  bindrad[0],  alpha,  cof);
	// sprintf(tname, "crt_sprob_x0_%g_dh0_%g_oat%g_obt%g_dt%g.dat",x0, dihed, ctab1, ctab_0, deltat_reac);
// 	tabfile.open(tname);
	  	
	sprintf(tname, "prt_sprob_x0_%g_dh0_%g_cat%g_cbt%g_dt%g.dat",x0, dihed,  costheta[ct], costheta[bt], deltat_reac);
	rfile.open(tname);
	rfile<<0<<'\t'<<phist[Nitbin-1]/(1.0*Nrep)<<'\t'<<1-passoc<<'\t'<<1-passoc<<endl;  
	
	for(i=0;i<Rbins;i++){
	  R1=bindrad[0]+delR*(i+0.5);
	  pirrev=pirrev_value(R1, x0, Maxtime,  Dtot, bindrad[0],  alpha);
	  pfree=pfree_value_norm(R1, x0,  Maxtime,  Dtot, bindrad[0], alpha);
	  
	  rfile<<R1<<'\t'<<pirrev<<'\t'<<pfree*(1-passoc)<<'\t'<<p1hist[i]/(1.0*Nrep*delR*R1*R1*4.0*M_PI)<<'\t';
	  for(k=0;k<dihbins;k++){
	    for(j=0;j<thetbins2;j++){
	      
	      rfile<<prhist[j][k][i]/(1.0*Nrep*delR*R1*R1*delthet1*delthet1*deldih*4.0*M_PI)<<'\t';
	      //sumr[i]+=prhist[j][k][i]/(1.0*Nrep*delR*R1*R1*delthet1*delthet1*deldih);
	      
	    }
	    
	    // if(k==0){
// 	      tabfile<<R1<<'\t'<<pirrev<<'\t';
// 	      for(j=0;j<tabbins*tabbins;j++){
// 		tabfile<<prhist2[j][i]/(1.0*Nrep*delR*R1*R1*deltab*deltab)<<'\t';
// 		sumr2[i]+=prhist2[j][i]/(1.0*Nrep*delR*R1*R1*deltab*deltab);
// 	      }
// 	      tabfile<<endl;
	    //}
	  }//end dihedral bins
	  rfile<<endl;
	  
	  
	}//end over r bins
	rfile.close();
	//tabfile.close();
	//sprintf(tname, "rmatch_x0_%g_dh%g_cat%g_cbt%g.dat",x0, dihed, costheta[ct], costheta[bt]);
	
	//matchfile.open(tname);
	//for(i=0;i<Rbins;i++){
	//R1=bindrad[0]+delR*(i+0.5);
	//matchfile<<R1<<'\t'<<p1hist[i]/(1.0*Nrep*delR*R1*R1*4.0*M_PI)<<'\t'<<sumr[i]<<endl;
	//}
	//matchfile.close();
      }//end second angle

    }//end first angle
  }//end over all x0s

  
  stop_timer(&totaltime);
  cout <<timer_duration(totaltime)<<" total time "<<endl;  
  /*Write out final result*/
  cout <<"End Main, complete run "<<endl;
  
}//end main


double pirr_pfree_ratio_ps(double rcurr, double r0, double tcurr, double Dtot, double bindrad, double alpha, double ps_prev, double rtol)
{
    
  //double pirr=pirrev_value(rcurr, r0, deltat, Dtot, bindrad, alpha_irr);
  //double pfree=pfree_value(rcurr, r0, deltat, Dtot, bindrad, alpha_irr, prevpassoc);
  
  double fDt=4.0*Dtot*tcurr;
  double sq_fDt=sqrt(fDt);

  double f1=1.0/(sqrt(4.0*M_PI*tcurr));
  double f2=1.0/(4.0*M_PI*r0*sqrt(Dtot));

  double sep, dist;
  double sqrt_t=sqrt(tcurr);
  double a2=alpha*alpha;
  double r1, term1, term2, e1, ef1, sum;
  
  r1=rcurr;
  sep=r1+r0-2.0*bindrad;
  dist=r1-r0;
  term1=f1*( exp(-dist*dist/fDt)+exp(-sep*sep/fDt) );
  
  e1=2.0*sep/sq_fDt*sqrt_t*alpha+a2*tcurr;
  ef1=sep/sq_fDt+alpha*sqrt_t;
  term2=alpha*exp(e1)*erfc(ef1);
  
  sum=term1-term2;
  sum*=f2/r1;
  double pirr=sum;

  
  /*Normalization for no diffusion inside binding radius!*/
  double cof=f1*f2;//equal to 1/( 8pir_0 sqrt(piDt))
  double c1=4.0*M_PI*cof;
  double ndist=bindrad-r0;
  double nadist=bindrad+r0;
  double sq_P=sqrt(M_PI);
  term1=-0.5*fDt*exp(-ndist*ndist/fDt)-0.5*sq_fDt*sq_P*r0*erf(-ndist/sq_fDt);
  term2=0.5*fDt*exp(-nadist*nadist/fDt)+0.5*sq_fDt*sq_P*r0*erf(nadist/sq_fDt);
  double pnorm=1.0-c1*(term1+term2);//this is then the normlization from sigma->inf 
  //  cout <<"Normlization: integratl from sigma to infinity: "<<pnorm<<endl;
  
  double adist=r1+r0;
  term1=exp(-dist*dist/fDt)-exp(-adist*adist/fDt);
  double pfree=cof/r1*term1/pnorm;//NORMALIZES TO ONE


  double ratio;
  if(abs(pirr-pfree*ps_prev)<rtol)ratio=1.0;
  else {
    ratio=pirr/(pfree*ps_prev);
    //    cout <<"pirr: "<<pirr<<" pfree*ps: "<<pfree*ps_prev<<" ratio: "<<ratio<<endl;
  }

  return ratio;
}
double pirrev_value(double rcurr, double r0, double tcurr, double Dtot, double bindrad, double alpha)
{
  double fDt=4.0*Dtot*tcurr;
  double sq_fDt=sqrt(fDt);

  double f1=1.0/(sqrt(4.0*M_PI*tcurr));
  double f2=1.0/(4.0*M_PI*r0*sqrt(Dtot));
  int i, j;
  double sep, dist;
  double sqrt_t=sqrt(tcurr);
  double a2=alpha*alpha;
  double r1, term1, term2, e1, ef1, sum;
  
  
  r1=rcurr;
  sep=r1+r0-2.0*bindrad;
  dist=r1-r0;
  
  term1=f1*( exp(-dist*dist/fDt)+exp(-sep*sep/fDt) );
  
  
  e1=2.0*sep/sq_fDt*sqrt_t*alpha+a2*tcurr;
  ef1=sep/sq_fDt+alpha*sqrt_t;
  term2=alpha*exp(e1)*erfc(ef1);
  
  sum=term1-term2;
  //    cout <<"irr 1: "<<term1<<" irr 2: "<<term2<<" sum: "<<sum<<" scaled: "<<sum*f2/r1<<endl;
  sum*=f2/r1;
  double pirrev=sum;
  return pirrev;
  
  
  
}

double pfree_value_norm(double rcurr, double r0, double tcurr, double Dtot, double bindrad,double alpha)
{
  /*Psurvive for the first time step is 1-passoc and should be exact. After
    more than one step the distribution will be off so passoc will be off as well*/
  
  double fDt=4.0*Dtot*tcurr;
  double sq_fDt=sqrt(fDt);

  double f1=1.0/(sqrt(4.0*M_PI*tcurr));
  double f2=1.0/(4.0*M_PI*r0*sqrt(Dtot));
  double cof=f1*f2;//equal to 1/( 8pir_0 sqrt(piDt))
  int i, j;
  double sep, dist;
  double sqrt_t=sqrt(tcurr);
  double a2=alpha*alpha;
  double r1, term1, term2, e1, ef1, sum;
  double adist;
  /*Normalization for no diffusion inside binding radius!*/
  double c1=4.0*M_PI*cof;
  dist=bindrad-r0;
  adist=bindrad+r0;
  term1=-0.5*fDt*exp(-dist*dist/fDt)-0.5*sqrt(4.0*M_PI*Dtot*tcurr)*r0*erf(-dist/sqrt(4.0*Dtot*tcurr));
  term2=0.5*fDt*exp(-adist*adist/fDt)+0.5*sqrt(4.0*M_PI*Dtot*tcurr)*r0*erf(adist/sqrt(4.0*Dtot*tcurr));
  double pnorm=1.0-c1*(term1+term2);//this is then the normlization from sigma->inf 
  
  //  cout <<"Normlization: integratl from sigma to infinity: "<<pnorm<<endl;

  r1=rcurr;
  //    sep=r1+r0-2.0*bindrad;
  dist=r1-r0;
  adist=r1+r0;
  term1=exp(-dist*dist/fDt)-exp(-adist*adist/fDt);
  
  
  double pfree=cof/r1*term1/pnorm;//NORMALIZES TO ONE
  return pfree;

  
    

}
double survive_irr(double r0, double tcurr, double Dtot, double bindrad, double alpha, double cof)
{
  double fDt=4.0*Dtot*tcurr;
  double sq_fDt=sqrt(fDt);

  double f1=cof*bindrad/r0;

  int i, j;
  double sep, dist;
  double sqrt_t=sqrt(tcurr);
  double a2=alpha*alpha;
  double r1, term1, term2, e1, ef1, sum;
  double passoc;
  sep=(r0-bindrad)/sq_fDt;

  e1=2.0*sep*sqrt_t*alpha+a2*tcurr;
  ef1=sep+alpha*sqrt_t;
  term1=erfc(sep);
  term2=exp(e1)*erfc(ef1);
  
  sum=term1-term2;
  sum*=f1;
  passoc=sum;
  return passoc;
  //  cout <<"s_irr: "<<sirr<<" time: "<<tcurr<<" unscaled: "<<term1-term2<<endl;
}

void read_parms(ifstream &parmfile, Parms &plist)
{
  parmfile >>plist.Nprotypes;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Nifaces;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Maxtime;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Nrxn;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Nspecies;
  parmfile.ignore(400,'\n');
  parmfile >>plist.rate;
  parmfile.ignore(400,'\n');
  parmfile >>plist.eps_scale;
  parmfile.ignore(400,'\n');
  parmfile >>plist.dt_scale;
  parmfile.ignore(400,'\n');
  parmfile >>plist.statwrite;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Nrep;
  parmfile.ignore(400,'\n');
  parmfile >>plist.D;
  parmfile.ignore(400,'\n');
  parmfile >>plist.leglen;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Nx0;
  parmfile.ignore(400,'\n');
  parmfile >>plist.x0;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Nthet;
  parmfile.ignore(400,'\n');
  parmfile >>plist.dihed0;
  parmfile.ignore(400,'\n');
  parmfile >>plist.thetbins;
  parmfile.ignore(400,'\n');
  parmfile >>plist.dihbins;
  parmfile.ignore(400,'\n');
  parmfile >>plist.tabbins;
  parmfile.ignore(400,'\n');

}
double GaussV()
{
  /*Box mueller method for gaussian distributed random number from a uniform
    random number generator~ trand()*/
  
  double R=2.0;
  double rnum;
  double V1, V2;

  while(R>=1.0)
    {

      V1=2.0*rand_gsl()-1.0;
      V2=2.0*rand_gsl()-1.0;
      R=V1*V1+V2*V2;
    }
  rnum=V1*sqrt(-2.0*log(R)/R);
  return rnum;

}
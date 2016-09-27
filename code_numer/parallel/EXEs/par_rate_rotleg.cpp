/*

Brownian Dynamics,
H.X. Zhou, JCP 1990, V94, 8794-8800

this program only propagates a single pair of particles,
collects statistics for different starting separations over
many repeats. 


both particles can rotate, they are both spheres but rotating
around a separate center, not there own center. 
------0, so a spherical site at the end of a leg. 

  
  
*/
#include <mpi.h>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "rand_gsl.h"
#include "md_timer.h"
#include "vector_rot_calls.h"

#define MAXIFACE 2
#define MAXPRTNER 2
#define MAXCOMPLEX 2
#define MAXRXN 2
#define MAXOVERLAP 2

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

  double bindrad;
  int Ndih;
  int Nx0;
  int Nthet;
  int thetbins_out;
  int dihbins_out;
  int tabbins;
  double x0;
  int Nprotypes;
  int Nifaces;
  double Nit;
  double leglen;
  double eps_scale;
  double dt_scale;
  double ka;
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
  MPI::Init(argc,argv);
  int rank, nprocs;
  //MPI::Intracomm full_world=MPI::COMM_WORLD;
  rank=MPI::COMM_WORLD.Get_rank();
  nprocs=MPI::COMM_WORLD.Get_size();
  
  int i, j, k;
  timeval tim;
  gettimeofday(&tim, 0);
  double t1=tim.tv_sec+tim.tv_usec;
  
  int seed=int(t1);
  seed = seed + 1000*rank;
  //seed=1353432282;
  double randmax=pow(2.0, 32);
  cout <<"rank: "<<rank<<'\t' <<"seed: "<<seed<<" randmax: "<<randmax<<endl;
  srand_gsl(seed);
  double irandmax=1.0/randmax;
  
  ifstream parmfile(argv[1]);
  Parms plist;
  
  int arrsize=20;
  double *dubparm=new double[arrsize];
  int *intparm=new int[arrsize];
  for(i=0;i<arrsize;i++){
    dubparm[i]=0;
    intparm[i]=0;
  }
  if(rank==0){
    read_parms(parmfile, plist);
    
    intparm[0]=plist.statwrite;
    intparm[1]=plist.Nrep;
    intparm[2]=plist.Ndih;
    intparm[3]=plist.Nthet;
    intparm[4]=plist.thetbins_out;
    intparm[5]=plist.dihbins_out;
    dubparm[0]=plist.Maxtime;
    dubparm[1]=plist.ka;
    dubparm[2]=plist.eps_scale;
    dubparm[3]=plist.dt_scale;
    dubparm[4]=plist.D;
    dubparm[5]=plist.leglen;
    dubparm[6]=plist.bindrad;
    
    
  }
  MPI::COMM_WORLD.Barrier();

  MPI::COMM_WORLD.Bcast(intparm, arrsize, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(dubparm, arrsize, MPI::DOUBLE, 0);

  plist.statwrite=intparm[0];
  plist.Nrep=intparm[1];
  plist.Ndih=intparm[2];
  plist.Nthet=intparm[3];
  plist.thetbins_out=intparm[4];
  plist.dihbins_out=intparm[5];
  plist.Maxtime=dubparm[0];
  plist.ka=dubparm[1];
  plist.eps_scale=dubparm[2];
  plist.dt_scale=dubparm[3];
  plist.D=dubparm[4];
  plist.leglen=dubparm[5];
  plist.bindrad=dubparm[6];
    
    
  cout.precision(12);
  initialize_timer(&totaltime);
  initialize_timer(&bimoltime);
  start_timer(&totaltime);
  
  int Nprotypes=2;
  int Nifaces=2;

  int Nrxn=1;


  int Ntotalmol=2;
  cout <<"rank: "<<rank<<'\t' <<"Ntotal mols: "<<Ntotalmol<<endl;//ASSUMED TO BE 2 MOLECULES BELOW
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

  double bindrad;

  bindrad=plist.bindrad;

  cout <<"rank: "<<rank<<'\t' <<"ACTIVATION RATE: "<<plist.ka<<"  radius: "<<bindrad<<endl;
  double kact=plist.ka;

    
  
  int ind, r1, m;
  double passoc;
  double rnum;
  
  double curr_time=0;

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
  
  double xtot, ytot, ztot;
  int nfree, wprot;
  int p, i1, i2;
  int np;
  double r2, r;

  int iind, iind2, ppart;
  
  int it;
  plist.ntotalcomplex=Ntotalmol;

  
  int s1;
  //  cout <<"rank: "<<rank<<'\t' <<"Ntotal complexes: "<<plist.ntotalcomplex<<endl;

  int amol,df; 
  double us_to_s=1E-6;
  int statwrite=plist.statwrite;

  double kpi=4.0*M_PI*bindrad;
  double leglen=plist.leglen;
  double leglen2=leglen*leglen;
  cout <<"rank: "<<rank<<'\t' <<"Leg length: "<<leglen<<endl;
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
  cout <<"rank: "<<rank<<'\t' <<"D: "<<plist.D <<" effective radius: "<<crad<<" D, calc: "<<wholep[1].Dx<<" Drot, calc: "<<wholep[1].Drx<<endl;
  wholep[0].Dx=0;
  wholep[0].Dy=0;
  wholep[0].Dz=0;
  wholep[0].Drx=wholep[1].Drx;
  wholep[0].Dry=wholep[1].Dry;
  wholep[0].Drz=wholep[1].Drz;
  
  double r1x, r1y, r1z;
  double cthet1;

  int ind_thet;
 
  
  
  double Dtrans=plist.D;//translation diffusion
  double fourpi=4.0*M_PI;
  
  double R2, R1;
  
  double fact;
  double kdiff;//will be kpi*D
  double aexp;
  double bexp;
  double alpha;
  
  mu=0;
  double Rmax;
  i=0;
  
  double rnum2;

  double tmpx, tmpy, tmpz;
  int p1, p2;
  int c1, c2;
  int ci1, ci2;

  int rxn1;

  double rate;

  int ind_rad;
  double rad2, rad;

  int nc1, nc2;
  double sep, ratio;
  int flagsep;
  it=1;
  int rep;
  int Nrep=plist.Nrep;
  int myNrep=int(round(1.0*Nrep/(1.0*nprocs)));
  Nrep=myNrep*nprocs;
  cout <<"rank: "<<rank<<'\t' <<"My reps: "<<myNrep<<" total reps: "<<Nrep<<endl;
  double realsmall=1E-14;
  double prevnorm=1.0;
  int Nx0=2.0;
  double dels;
  
  ofstream probfile;
  double space;
  double tval;
  char tname[200];
    char tname2[200];
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
  cout <<"rank: "<<rank<<'\t' <<"Initial separation: "<<bindrad<<endl;
  x0vec[0]=bindrad;//for k<inf, this will not go to 1 
  
  double kappa=kact/(4.0*M_PI*bindrad*bindrad);
  cout <<"rank: "<<rank<<'\t' <<"Kact: "<<kact<<" kappa: "<<kappa<<endl;
  double epsilon=plist.eps_scale*Dtrans/kappa;
  double tau=epsilon/kappa;
  x0vec[1]=bindrad+epsilon;
  cout <<"rank: "<<rank<<'\t' <<"epsscale: "<<plist.eps_scale<<" epsilon: "<<epsilon<<" tau: "<<tau<<endl;
  Rmax=bindrad+epsilon;
  cout <<"rank: "<<rank<<'\t' <<"Bindrad: "<<bindrad<<" reaction limit: "<<Rmax<<endl;
  //cout <<"rank: "<<rank<<'\t' <<"different time step expansions: "<<epsilon*epsilon/Dtrans/100.0<<" other: "<<Rmax*Rmax*0.0001/(2.0*Dtrans)<<endl;
  double scaled=plist.dt_scale;
  double deltat_reac=scaled*epsilon*epsilon/(2.0*Dtrans);

  double cf=cos(sqrt(4.0*wholep[1].Drx*deltat_reac));
  double Dr1=2.0*leglen2*(1.0-cf);

  double Dtot=Dtrans+2.0*Dr1/(6.0*deltat_reac);//add in rotation for both molecules
  cout <<"rank: "<<rank<<'\t' <<"add to Dtot from Rotation: "<<2.0*Dr1/(6.0*deltat_reac)<<" original dt: "<<deltat_reac<<" final Dtot: "<<Dtot<<endl;
  /*Now update to new Dtot*/
  epsilon=plist.eps_scale*Dtot/kappa;
  tau=epsilon/kappa;
  Rmax=bindrad+epsilon;
  deltat_reac=scaled*epsilon*epsilon/(2.0*Dtot);
  double deltat=deltat_reac;
    
  /*to calculate rate for spheres, where Dtot is Deffective (includes effect of rotation*/
  kdiff=fourpi*Dtot*bindrad;
  fact=1.0+kact/kdiff;
  alpha=fact*sqrt(Dtot)/bindrad;
  double cof=kact/(kact+kdiff);
  
  cout <<"rank: "<<rank<<'\t' <<" deltat: "<<deltat_reac<<" epsilon: "<<epsilon<<" tau: "<<tau<<" Rmax: "<<Rmax<<endl;
  double currtime=0;
  double Maxtime=plist.Maxtime;
  int Nitbin=int(Maxtime/deltat_reac);
  Maxtime=Nitbin*deltat_reac;
  cout <<"rank: "<<rank<<'\t' <<"number of time bins: "<<Nitbin<<" new max time: "<<Maxtime<<endl;
  double pirrev, pfree;
  //current statwrite is the number of datapoints it will write out.
  if(Nitbin<statwrite)
    statwrite=500;
  else
    statwrite=int(round(Nitbin/statwrite));
  

  
  double theta;
  int Nthet=plist.Nthet;
  double deltheta=2.0/(1.0*(Nthet-1));
  cout <<"rank: "<<rank<<'\t' <<"Ntheta: "<<Nthet<<" delthet: "<<deltheta<<endl;
  double *costheta=new double[Nthet];
  for(i=0;i<Nthet;i++){
    costheta[i]=-1+i*deltheta;
    cout <<"rank: "<<rank<<'\t' <<"costheta: "<<costheta[i]<<endl;
  }
  int Ndih=plist.Ndih;
  double deldih=M_PI/(1.0*(Ndih));//dihedral relative to azimuths only range over pi.
  double *dihvec=new double[Ndih];
  for(i=0;i<Ndih;i++){
    dihvec[i]=0+i*deldih;
  }
  int Ndihsim;
  double dihnorm;
  int Nwrite=int(ceil(Nitbin/statwrite));
  cout <<"rank: "<<rank<<'\t' <<"statwrite: "<<statwrite<<" Nwrite bins, for S(t): "<<Nwrite<<endl;
  if((Nitbin-1)==Nwrite*statwrite)
    Nitbin-=1;
  double *phist=new double[Nwrite+1];
  double *phist_sum=new double[Nwrite+1];
  for(i=0;i<Nwrite+1;i++)
    phist[i]=0.0;
  
  
  int cnt;
  
  int n;
  int itmove;
  int itmax;
  double lz, ly, lx;
  double tx, ty, tz;
  /*If we save distribution only need separations out to the reaction zone.*/
  double Rrange=3.0*sqrt(6.0*Dtot*Maxtime);//no +bindrad here because only measure distance beyond sigma.

  
  double delR=0.01;
  int Rbins=int(Rrange/delR);
  if(Rbins>1000){
    Rbins=1000;
    delR=(Rrange)/(1.0*Rbins);
  }
  cout <<"rank: "<<rank<<'\t' <<"MaxR: "<<Rrange<<" Rbins: "<<Rbins<<" delR: "<<delR<<endl;
  int rind;
  double small=1E-9;
  int thetbins_out=plist.thetbins_out;
  int thetbins_out2=thetbins_out*thetbins_out;
  double delthet1=2.0/(1.0*thetbins_out)+small;
  cout <<"rank: "<<rank<<'\t' <<"theta bins: "<<thetbins_out<<" delthet1: "<<delthet1<<endl;
  ofstream thetfile("thetas.out");
  for(i=0;i<thetbins_out;i++)
    thetfile <<(i+0.5)*delthet1-1<<endl;
  
  ofstream ratefile;
  ofstream transfile;
  int dihbins_out=plist.dihbins_out;
  double deldih_out=M_PI/(1.0*dihbins_out)+small;
  cout <<"rank: "<<rank<<'\t' <<"dihedral bins: "<<dihbins_out<< " delta: "<<deldih_out<<endl;
  
  double *prhist=new double[thetbins_out2*dihbins_out*Rbins];
  double *prhist_sum=new double[thetbins_out2*dihbins_out*Rbins];
  // for(i=0;i<thetbins_out2;i++)
//     prhist[i]=new double*[dihbins_out];
//   for(i=0;i<thetbins_out2;i++){
//     for(k=0;k<dihbins_out;k++){
//       prhist[i][k]=new double[Rbins];
//     }
//   }
  for(i=0;i<thetbins_out2;i++){
    for(k=0;k<dihbins_out;k++){
      for(j=0;j<Rbins;j++)
	prhist[i*thetbins_out2+k*dihbins_out+j]=0.0;
    }
  }
  double *p1hist=new double[Rbins];
  double *p1hist_sum=new double[Rbins];
  for(i=0;i<Rbins;i++)
    p1hist[i]=0.0;

  double *M=new double[9];
  double *v3=new double[3];
  double *v2=new double[3];
  double *v1=new double[3];
  double *v=new double[3];
  ofstream rfile;

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
  double *stot=new double[Nwrite+1];
  double *ktime=new double[Nwrite+1];
  double *savg=new double[Nwrite+1];
  int hh;
  double dihed;
  for(i=0;i<Nwrite+1;i++)
      savg[i]=0.0;

  for(s=0;s<Nx0;s++){
    /*For rate calculation, only loop over x0=sigma, and x0=sigma+epsilon. Average survival probs to get rate.*/
    x0=x0vec[s];
    for(i=0;i<Nwrite+1;i++)
      stot[i]=0.0;

    for(ct=0;ct<Nthet;ct++){
      lz=costheta[ct]*leglen;
      lx=sqrt(leglen2-lz*lz);
      ly=0;
      cout <<"rank: "<<rank<<'\t' <<"Initial pos for leg1 com: "<<lx<<' '<<ly<<' '<<x0+lz<<"  for leg: "<<0<<' '<<0<<' '<<x0<<endl;
      
      
      for(bt=0;bt<Nthet;bt++){
	
	if(costheta[ct]==1 ||costheta[bt]==1 || costheta[ct]==-1 || costheta[bt]==-1){
	  Ndihsim=1;
	  dihnorm=M_PI;
	}else{
	  Ndihsim=Ndih;
	  dihnorm=deldih;
	}
	for(hh=0;hh<Ndihsim;hh++){
	  dihed=dihvec[hh];
	  
	  lz_0=-costheta[bt]*leglen;
	  lR=sqrt(leglen2-lz_0*lz_0);
	  lx_0=lR*cos(dihed);
	  ly_0=lR*sin(dihed);
	  cout <<"rank: "<<rank<<'\t' <<"Initial pos for leg2 com: "<<lx_0<<' '<<ly_0<<' '<<lz_0<<"  for leg: "<<0<<' '<<0<<' '<<0<<endl;
	  
	  for(i=0;i<Rbins;i++){
	    p1hist[i]=0.0;
	    p1hist_sum[i]=0.0;
	  }
	  for(i=0;i<thetbins_out2;i++){
	    for(k=0;k<dihbins_out;k++){
	      for(j=0;j<Rbins;j++){
		prhist[i*thetbins_out2+k*dihbins_out+j]=0.0;
		prhist_sum[i*thetbins_out2+k*dihbins_out+j]=0.0;
	      }
	    }
	  }
	  for(i=0;i<Nwrite+1;i++){
	    phist[i]=0.0;
	    phist_sum[i]=0.0;
	  }
	  
	  cout <<"rank: "<<rank<<'\t' <<"X0: "<<x0<<" costheta: "<<costheta[ct]<<" costheta2: "<<costheta[bt]<<endl;
	  /*Need to sample over different values of the sphere holding
	    all the points at leg length. origin of sphere is at z=x0, x=y=0,
	    azimuth doesn't matter, choose phi=0, so y=0 always to start.
	    then sample evenly over cos(theta),
	    x=d1*sin(theta), y=0, z=z0+d1*cos(theta).
	    You will perform rotation around this point.
	  */
	  
	  for(rep=0;rep<myNrep;rep++){
	    /*start proteins separated by x0 */
	    //cout <<"rank: "<<rank<<'\t' <<"Repeat: "<<rep<<endl;
	    
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
	    //cout <<"rank: "<<rank<<'\t' <<"rep: "<<rep<<" Maxtime; "<<Maxtime<<endl;
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
		prob=exp(-deltat_reac/tau);
		
		/*might perform this reaction, depending on k_associate*/
		rnum=1.0*rand_gsl();
		p1=1;
		if(prob<rnum){
		  //terminate trajectory
		  //cout <<"rank: "<<rank<<'\t' <<"Associate, prob:  "<<rnum<<" rep: "<<rep<<" time: "<<currtime<<endl;
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
		  while(R2<bindrad*bindrad){
		    
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
		deltat=scaled*(R1-bindrad)*(R1-bindrad)/(Dtot*2.0);//+deltat_reac;
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
	      
	      itmax=itmove+1;
	      if((itmax+it)>Nitbin)
		itmax=Nitbin-it;
	      
	      for(n=1;n<itmax;n++){
		if((it+n)%statwrite==0){
		  ind=(it+n)/statwrite;
		  phist[ind]+=1;
		}
	      }
	      
	      
	      it+=itmove;
	      
	    }
	    if(rep%100==0)
	      cout <<"rank: "<<rank<<'\t' <<"finished rep: "<<rep<< " at time: "<<currtime<<endl;
	    /*Keep a histogram of final positions for particles that survived until the end*/
	    /*if(currtime>Maxtime){
	      
	      dx=bases[1].x[0]-bases[0].x[0];
	      dy=bases[1].y[0]-bases[0].y[0];
	      dz=bases[1].z[0]-bases[0].z[0];
	      
	      R2=dx*dx+dy*dy+dz*dz;
	      R1=sqrt(R2);
	      rind=int((R1-bindrad)/delR);
	      r0x=bases[0].x[0]-bases[0].xcom;
	      r0y=bases[0].y[0]-bases[0].ycom;
	      r0z=bases[0].z[0]-bases[0].zcom;
	      
	      r1x=bases[1].x[0]-bases[1].xcom;
	      r1y=bases[1].y[0]-bases[1].ycom;
	      r1z=bases[1].z[0]-bases[1].zcom;
	      
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
	      dihbin=int(dih/deldih_out);
	      if(isnan(dih)){
		
		if(round(cosdih)==-1)dihbin=dihbins_out-1;
		else dihbin=0;
		cout <<"rank: "<<rank<<'\t' <<"NAN: "<<dih<<" cosdih: "<<cosdih<<" bin: "<<dihbin<<endl;
	      }else if(dihbin==dihbins_out)
		dihbin=dihbins_out-1;
	      
	      
	      if(ind_thet>=thetbins_out){
		// 	      cout <<"rank: "<<rank<<'\t' <<"cos(theta)>1 ?: "<<cthet1<<" ind: "<<ind_thet<<endl;
		ind_thet=thetbins_out-1;
	      }else if(ind_thet<0){
		// 	      cout <<"rank: "<<rank<<'\t' <<"cos(theta)<-1 ?: "<<cthet1<<" ind: "<<ind_thet<<endl;
		ind_thet=0;
	      }
	      if(rind<Rbins){
		prhist[(ind_thet*thetbins_out+ind_thet2)*thetbins_out2+dihbin*dihbins_out+rind]++;
		//prhist2[tab0ind*tabbins+tab1ind][rind]++;
		p1hist[rind]++;
	      }
	      
	      }*/
	  }//end over all reps
	  /*Sum the energies across processors*/
	  //MPI::COMM_WORLD.Reduce(p1hist, p1hist_sum, Rbins, MPI::DOUBLE, MPI::SUM, 0);
	  MPI::COMM_WORLD.Reduce(phist, phist_sum, Nwrite+1, MPI::DOUBLE, MPI::SUM, 0);
	  //  MPI::COMM_WORLD.Reduce(prhist, prhist_sum, thetbins_out2*dihbins_out*Rbins, MPI::DOUBLE, MPI::SUM, 0);
	  MPI::COMM_WORLD.Barrier();
	  cout <<"rank: "<<rank<<'\t' <<"finished all reps "<<endl;
	  /*Now write out final probability histogram*/
	  //sprintf(tname, "rot_sprob_x0_%g_dh%g_ct%g_c2t_%g_dt%g.dat",x0, dihed, costheta[ct], costheta[bt], deltat_reac);
	  //cout <<"rank: "<<rank<<'\t' <<"open file: "<<tname<<endl;
	  if(rank==0){
	    // probfile.open(tname);
// 	    probfile<<0<<' '<<1<<' '<<1<<endl;
	  
	    
	    double integrate_theta=1.0;
	    if(deltheta*Nthet-2.0>1E-6){
	      /*first and last bins need factor of 0.5 for integration*/
	      if(bt==0)integrate_theta*=0.5;
	      else if(bt==Nthet-1)integrate_theta*=0.5;
	      if(ct==0)integrate_theta*=0.5;
	      else if(ct==Nthet-1)integrate_theta*=0.5;
	      
	    }
	    for(i=statwrite;i<Nitbin;i+=statwrite){
	      tval=i*deltat_reac;
	      //if(i%statwrite==0){
	      ind=i/statwrite;
	      //passoc=survive_irr( x0, tval,  Dtot,  bindrad,  alpha,  cof);
	      //probfile<<tval<<' ' <<phist_sum[ind]/(1.0*Nrep)<<' '<<1.0-passoc<<endl;
	      stot[ind]+=phist_sum[ind]/(1.0*Nrep)*deltheta*deltheta*dihnorm*integrate_theta;//dihnorm is either deldih, which when done Ndih times is pi, or it's pi, done once.
	      
	    }
	   //  i=Nitbin-1;
// 	    tval=i*deltat_reac;
// 	    passoc=survive_irr( x0, tval,  Dtot,  bindrad,  alpha,  cof);
// 	    probfile<<tval<<' ' <<phist_sum[Nwrite]/(1.0*Nrep)<<' '<<1.0-passoc<<endl;//index Nwrite is for the last time point
// 	    probfile.close();
	    
	    //passoc=survive_irr( x0, Maxtime,  Dtot,  bindrad,  alpha,  cof);
	    
	    
	   //  sprintf(tname, "prt_sprob_x0_%g_dh0_%g_cat%g_cbt%g_dt%g.dat",x0, dihed,  costheta[ct], costheta[bt], deltat_reac);
// 	    rfile.open(tname);
// 	    rfile<<0<<'\t'<<phist_sum[Nwrite]/(1.0*Nrep)<<'\t'<<1-passoc<<'\t'<<1-passoc<<endl;  
	    
// 	    for(i=0;i<Rbins;i++){
// 	      R1=bindrad+delR*(i+0.5);
// 	      pirrev=pirrev_value(R1, x0, Maxtime,  Dtot, bindrad,  alpha);
// 	      pfree=pfree_value_norm(R1, x0,  Maxtime,  Dtot, bindrad, alpha);
	      
// 	      rfile<<R1<<'\t'<<pirrev<<'\t'<<pfree*(1-passoc)<<'\t'<<p1hist_sum[i]/(1.0*Nrep*delR*R1*R1*4.0*M_PI)<<'\t';
// 	      for(k=0;k<dihbins_out;k++){
// 		for(j=0;j<thetbins_out2;j++){
		  
// 		  rfile<<prhist_sum[j*thetbins_out2+k*dihbins_out+i]/(1.0*Nrep*delR*R1*R1*delthet1*delthet1*deldih_out*4.0*M_PI)<<'\t';
// 		}
// 	      }//end dihedral bins
// 	      rfile<<endl;
	      
	      
// 	    }//end over r bins
// 	    rfile.close();
	    
	  }//end writing out, rank0
	  MPI::COMM_WORLD.Barrier();

	}//End dihedral angle loop
      }//end second angle

    }//end first angle

    sprintf(tname, "rate_at_x0_%g_dt%g.dat",x0, deltat_reac);
    if(rank==0){
      
      ratefile.open(tname);
      for(i=statwrite;i<Nitbin;i+=statwrite){
	tval=i*deltat_reac;
	ind=i/statwrite;
	/*pi is from dihedral, 
	  2*2 is from integrals over cos(theta). 
	  4pisigma^2*kappa is the integral over r=sigma surface. This vector can rotate fully and particles still align
	*/
	savg[ind]+=stot[ind]*0.5;//Average over both x0
	ktime[ind]=4.0*M_PI*bindrad*bindrad*kappa*stot[ind]/(M_PI*2.0*2.0);//normalizing factor over full angular range noted above
	ratefile<<tval<<'\t'<<ktime[ind]<<endl;
	cout<<" time, rate: "<<tval<<'\t'<<ktime[ind]<<" stot: "<<stot[ind]<<endl;
      }
      ratefile.close();
    }

  }//end over all x0s
  sprintf(tname, "rate_vs_time_calc_dt%g.dat",deltat_reac);
  sprintf(tname2, "rate_for_Dtrans_only_dt%g.dat",deltat_reac);
  if(rank==0){
    
    ratefile.open(tname);
    transfile.open(tname2);
    double kdiff2=fourpi*Dtrans*bindrad;
    double alpha2=(1.0+kact/kdiff2)*sqrt(Dtrans)/bindrad;
    double cof2=kact/(kact+kdiff2);
      
    for(i=statwrite;i<Nitbin;i+=statwrite){
      tval=i*deltat_reac;
      ind=i/statwrite;
      /*pi is from dihedral, 
	2*2 is from integrals over cos(theta). 
	4pisigma^2*kappa is the integral over r=sigma surface. This vector can rotate fully and particles still align
      */
      /*From Zhou, rate is average of calculated survival at x0=sigma and x0=epsilon*/
      ktime[ind]=4.0*M_PI*bindrad*bindrad*kappa*savg[ind]/(M_PI*2.0*2.0);//normalizing factor over full angular range noted above
      
      /*Calc rate for approx (with Deffective)*/
      passoc=survive_irr( bindrad, tval,  Dtot,  bindrad,  alpha,  cof);
      ratefile<<tval<<'\t'<<ktime[ind]<<'\t'<<4.0*M_PI*bindrad*bindrad*kappa*(1.0-passoc)<<endl;
      cout<<" time, rate: "<<tval<<'\t'<<ktime[ind]<<" stot: "<<stot[ind]<<endl;
      
      passoc=survive_irr( bindrad, tval,  Dtrans,  bindrad,  alpha2, cof2 );
      transfile<<tval<<'\t'<<4.0*M_PI*bindrad*bindrad*kappa*(1.0-passoc)<<endl;
    }
    ratefile.close();
  }

  
  stop_timer(&totaltime);
  cout <<"rank: "<<rank<<'\t' <<timer_duration(totaltime)<<" total time "<<endl;  
  /*Write out final result*/
  cout <<"rank: "<<rank<<'\t' <<"End Main, complete run "<<endl;
  MPI::COMM_WORLD.Barrier();
  MPI::Finalize();
  
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

  parmfile >>plist.Maxtime;
  parmfile.ignore(400,'\n');
  parmfile >>plist.ka;
  parmfile.ignore(400,'\n');
  parmfile >>plist.D;
  parmfile.ignore(400,'\n');
  parmfile >>plist.leglen;
  parmfile.ignore(400,'\n');
  parmfile >>plist.bindrad;
  parmfile.ignore(400,'\n');
  parmfile >>plist.eps_scale;
  parmfile.ignore(400,'\n');
  parmfile >>plist.dt_scale;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Nrep;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Nthet;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Ndih;
  parmfile.ignore(400,'\n');
  parmfile >>plist.thetbins_out;
  parmfile.ignore(400,'\n');
  parmfile >>plist.dihbins_out;
  parmfile.ignore(400,'\n');
  parmfile >>plist.statwrite;
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

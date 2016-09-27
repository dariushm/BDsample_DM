/*

this program determines the survival probability
as a function of time averaged over multiple repetitions
for each set of initial separations and orientations.

In addition to translation diffusion, both leg vectors
can rotate, and the orientation of the molecule is also
kept track of. 
Therefull there are full orientational restrictions on the
associaton, legs must orient to one another, and molecules
must be oriented in same direction.

in this version, only loop over initial touching orientations that
have a kappa value of nonzero. 
Either select starting positions from this range of angles, by randomly
sampling cos(theta), cos(theta2), psi, and phi1, or systematically loop over
each one. 
Instead of doing more reps at fewer points, try more points and fewer reps
to better cover the range of starting values. 
If you randomly sample in the range, you end up with S(t|sigma),
which then is multiplied by k(0)=integral over all angles with nonzero kappa.
For systematic stepping, we'd get S(t|r0, theta1, theta2, psi) and then
integrate those against k. 
  
  
*/
#include <mpi.h>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "rand_gsl.h"
#include "md_timer.h"
#include "rd_class_defs.h"
#include "utility_calls.h"
#include "Vrnumer.h"
#include "GF_calls.h"
#include "vector_rot_calls.h"
#include "angle_calc.h"

using namespace std;

struct MD_Timer totaltime;
struct MD_Timer bimoltime;


class survParms
{
public:
  double bindrad;
  int Ndih;
  int Nx0;
  int Nthet;
  int thetbins_out;
  int dihbins_out;
  int psibins_out;
  
  int Npsibins;
  int Nphi;

  double sprob_tol;
  double dihedstart;
  double costheta_range;
  double psi_range;
  double x0;
  int Nprotypes;
  int Nifaces;
  double Nit;
  double leglen1;
  double leglen0;
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



void read_parms(ifstream &parmfile, survParms &plist);
void write_crds(Fullmol *bases, int p1);
void run_sim_at_1angleset(int rank, int nprocs, double psiread, double costheta1, double costheta2, double &norm_psi, double &sprobvalue, int Ndihsim, int Nphi, double x0, double lx, double ly, double lz, double leglen1, double leglen0, survParms &plist, double deltat_reac, double scaled, double bindrad, Fullmol *bases, Protein *wholep, double *dihvec, int myNrep, double Rmax, double tau, double anglelimit, double psi_limit, int Nitbin, int Nwrite, double*phist, double *phist_sum, double *p1hist, double *p1hist_sum,double *prhist, double *prhist_sum, ofstream &rfile, ofstream &probfile, ofstream &pavfile, ofstream &jacfile, ofstream &sonlyfile, double tol, double *dih_intervals);
void loop_all_psi(int rank, int nprocs, double *psivec, double costheta1, double costheta2, double *norm_psi, double *sprobvec, int Ndihsim, int Nphi, double x0, double lx, double ly, double lz, double leglen1, double leglen0, survParms &plist, double deltat_reac, double scaled, double bindrad, Fullmol *bases, Protein *wholep, double *dihvec, int myNrep, double Rmax, double tau, double anglelimit, double psi_limit, int Nitbin, int Nwrite, double*phist, double *phist_sum, double *p1hist, double *p1hist_sum,double *prhist, double *prhist_sum, ofstream &rfile, ofstream &probfile, ofstream &pavfile, ofstream &jacfile, ofstream &sonlyfile, double tol, int ct1ind, int ct2ind, double *dih_intervals);
void loop_over_costheta2_gridpt(int rank, int nprocs, double *psivec, double costheta1, double *costheta, double *norm_psi, double *sprobvec, int Ndihsim, int Nphi, double x0, double lx, double ly, double lz, double leglen1, double leglen0, survParms &plist, double deltat_reac, double scaled, double bindrad, Fullmol *bases, Protein *wholep, double *dihvec, int myNrep, double Rmax, double tau, double anglelimit, double psi_limit, int Nitbin, int Nwrite, double*phist, double *phist_sum, double *p1hist, double *p1hist_sum,double *prhist, double *prhist_sum, ofstream &rfile, ofstream &probfile, ofstream &pavfile, ofstream &jacfile, ofstream &sonlyfile, double tol, int ct1ind, int ctstart, double *dih_intervals );
void loop_over_costheta2_halfgridpt(int rank, int nprocs, double *psivec, double costheta1, double costheta_next, double *norm_psi, double *sprobvec, int Ndihsim, int Nphi, double x0, double lx, double ly, double lz, double leglen1, double leglen0, survParms &plist, double deltat_reac, double scaled, double bindrad, Fullmol *bases, Protein *wholep, double *dihvec, int myNrep, double Rmax, double tau, double anglelimit, double psi_limit, int Nitbin, int Nwrite, double*phist, double *phist_sum, double *p1hist, double *p1hist_sum,double *prhist, double *prhist_sum, ofstream &rfile, ofstream &probfile, ofstream &pavfile, ofstream &jacfile, ofstream &sonlyfile, double tol, int ct1ind, int ct1start, double *dih_intervals);
void loop_over_costheta2_quartergridpt(int rank, int nprocs, double *psivec, double costheta1, double costheta_next, double *norm_psi, double *sprobvec, int Ndihsim, int Nphi, double x0, double lx, double ly, double lz, double leglen1, double leglen0, survParms &plist, double deltat_reac, double scaled, double bindrad, Fullmol *bases, Protein *wholep, double *dihvec, int myNrep, double Rmax, double tau, double anglelimit, double psi_limit, int Nitbin, int Nwrite, double*phist, double *phist_sum, double *p1hist, double *p1hist_sum,double *prhist, double *prhist_sum, ofstream &rfile, ofstream &probfile, ofstream &pavfile, ofstream &jacfile, ofstream &sonlyfile, double tol, int ct1ind, int ct1start, double *dih_intervals);
void loop_over_costheta2_eighthgridpt(int rank, int nprocs, double *psivec, double costheta1, double costheta_next, double *norm_psi, double *sprobvec, int Ndihsim, int Nphi, double x0, double lx, double ly, double lz, double leglen1, double leglen0, survParms &plist, double deltat_reac, double scaled, double bindrad, Fullmol *bases, Protein *wholep, double *dihvec, int myNrep, double Rmax, double tau, double anglelimit, double psi_limit, int Nitbin, int Nwrite, double*phist, double *phist_sum, double *p1hist, double *p1hist_sum,double *prhist, double *prhist_sum, ofstream &rfile, ofstream &probfile, ofstream &pavfile, ofstream &jacfile, ofstream &sonlyfile, double tol, int ct1ind, double *dih_intervals);

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
  
  //seed=1353432282;
  double randmax=pow(2.0, 32);
  cout <<"rank: "<<rank<<"seed: "<<seed<<" randmax: "<<randmax<<endl;
  srand_gsl(seed);
  double irandmax=1.0/randmax;
  
  ifstream parmfile(argv[1]);
  survParms plist;
  int arrsize=20;
  double *dubparm=new double[arrsize];
  int *intparm=new int[arrsize];
  for(i=0;i<arrsize;i++){
    dubparm[i]=0;
    intparm[i]=0;
  }
  if(rank==0){
    read_parms(parmfile, plist);
    intparm[0]=plist.Nprotypes;
    intparm[1]=plist.Nifaces;
    intparm[2]=plist.Nrxn;
    intparm[3]=plist.Nspecies;
    intparm[4]=plist.statwrite;
    intparm[5]=plist.Nrep;
    intparm[6]=plist.Nthet;
    intparm[7]=plist.Npsibins;
    intparm[8]=plist.Nphi;
    intparm[9]=plist.Ndih;
    intparm[10]=plist.thetbins_out;
    intparm[11]=plist.dihbins_out;
    intparm[12]=plist.psibins_out;
    intparm[13]=plist.Nx0;
    
    dubparm[0]=plist.Maxtime;
    dubparm[1]=plist.rate;
    dubparm[2]=plist.eps_scale;
    dubparm[3]=plist.dt_scale;
    dubparm[4]=plist.D;
    dubparm[5]=plist.leglen1;
    dubparm[6]=plist.x0;
    dubparm[7]=plist.costheta_range;
    dubparm[8]=plist.psi_range;
    dubparm[9]=plist.bindrad;
    dubparm[10]=plist.dihedstart;
    dubparm[11]=plist.sprob_tol;
    dubparm[12]=plist.leglen0;
    cout <<"rank 0, bindrad: "<<plist.bindrad<<endl;
    
    
  }
  MPI::COMM_WORLD.Barrier();

  MPI::COMM_WORLD.Bcast(intparm, arrsize, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(dubparm, arrsize, MPI::DOUBLE, 0);

  plist.Nprotypes=intparm[0];
  plist.Nifaces=intparm[1];
  plist.Nrxn=intparm[2];
  plist.Nspecies=intparm[3];
  plist.statwrite=intparm[4];
  plist.Nrep=intparm[5];
  plist.Nthet=intparm[6];
  plist.Npsibins=intparm[7];
  plist.Nphi=intparm[8];
  plist.Ndih=intparm[9];
  plist.thetbins_out=intparm[10];
  plist.dihbins_out=intparm[11];
  plist.psibins_out=intparm[12];
  plist.Nx0=intparm[13];

  plist.Maxtime=dubparm[0];
  plist.rate=dubparm[1];
  plist.eps_scale=dubparm[2];
  plist.dt_scale=dubparm[3];
  plist.D=dubparm[4];
  plist.leglen1=dubparm[5];
  plist.x0=dubparm[6];
  plist.costheta_range=dubparm[7];
  plist.psi_range=dubparm[8];
  plist.bindrad=dubparm[9];
  plist.dihedstart=dubparm[10];
  plist.sprob_tol=dubparm[11];
  plist.leglen0=dubparm[12];
  
  cout.precision(12); 
  initialize_timer(&totaltime);
  initialize_timer(&bimoltime);
  start_timer(&totaltime);
  
  int Nprotypes=plist.Nprotypes;//total distinct protein types, 
  int Nifaces=plist.Nifaces;//this is the number of interfaces

  int Nrxn=plist.Nrxn;
  int Nspecies=plist.Nspecies;//this will include product species


  int *Ncopy=new int[Nprotypes];

  Ncopy[0]=1;
  Ncopy[1]=1;
  int Ntotalmol=0;
  plist.Natom=0;
  int ntmp;
  for(i=0;i<Nprotypes;i++){

    Ntotalmol+=Ncopy[i];
  }
  cout <<"rank: "<<rank<<"Ntotal mols: "<<Ntotalmol<<endl;//ASSUMED TO BE 2 MOLECULES BELOW
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
  Speclist[0][1]=-1;
  Speclist[1][1]=-1;

  char fname[100];
  
  Protein *wholep=new Protein[Nprotypes];
  int *p_home=new int[Nifaces];//this reverses and tells you what protein a given interface belongs to
  int *i_home=new int[Nspecies];//for both free and bound states, what index are you on the protein's list
  //double *bindrad=new double[Nrxn];//binding or unbinding radius for each reaction
  double bindrad;
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
  
  wholep[0].ninterface=2;
  wholep[1].ninterface=2;
  bases[0].ninterface=2;
  bases[1].ninterface=2;
  bindrad=plist.bindrad;
  kr[0]=plist.rate;
  cout <<"rank: "<<rank<<"ACTIVATION RATE: "<<kr[0]<<"  radius: "<<bindrad<<endl;
  double kact=kr[0];

  for(i=0;i<Nprotypes;i++){  
    ntmp=wholep[i].ninterface+1;
    plist.Natom+=Ncopy[i]*ntmp;
  }
  cout <<"rank: "<<rank<<"N atoms: "<<plist.Natom<<endl;
  cout <<"rank: "<<rank<<"read reactions "<<endl;

  
  double *savecrds=new double[Ntotalmol*3];//for x, y, z
  
  /*Print out specific reactions*/
  // cout <<"rank: "<<rank<<"Print specific interaction network "<<endl;
//   int ncomplex=0;
//   for(i=0;i<Nifaces;i++){
//     cout <<"rank: "<<rank<<i<<'\t';
//     for(j=0;j<numpartners[i];j++){
//       cout <<"rank: "<<rank<<Speclist[i][j]<<'\t';
//       ncomplex++;
//     }
//     cout <<"rank: "<<rank<<endl;
//   }
//   ncomplex/=2;
//   plist.nspec_complex=ncomplex;
  

  
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
  double xtot, ytot, ztot;
  int nfree, wprot;
  int p, i1, i2;
  int np;
  double r2, r;
  double maxsep2=bindrad*bindrad;//plist.maxsep2;
  cout <<"rank: "<<rank<<"squared distance cutoff: "<<maxsep2<<endl;
  int iind, iind2, ppart;
  int twrite=plist.configwrite;



  
  int it;
  plist.ntotalcomplex=Ntotalmol;

  
  int s1;
  cout <<"rank: "<<rank<<"Ntotal complexes: "<<plist.ntotalcomplex<<endl;

  int amol,df; 
  double us_to_s=1E-6;
  int statwrite=plist.statwrite;

  double kpi=4.0*M_PI*bindrad;
  double leglen1=plist.leglen1;
  double leglen0=plist.leglen0;
  double leglen1_sq=leglen1*leglen1;
  double leglen0_sq=leglen0*leglen0;
  cout <<"rank: "<<rank<<"Leg length: "<<leglen1<<" leglen0: "<<leglen0<<endl;
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
  cout <<"rank: "<<rank<<"D: "<<plist.D <<" effective radius: "<<crad<<" D, calc: "<<wholep[1].Dx<<" Drot, calc: "<<wholep[1].Drx<<endl;
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
  
  r0=bindrad;
  int r0bins=1000;
  double delr0;
  int r0ind;
  
  
  double p0_ratio=1.0;
  //int distflag=0;
  double currsep;
  double prevpassoc=0;
  double probvec1;
  //cout <<"rank: "<<rank<<"Rmax1: "<<Rmax1
  cout <<"rank: "<<rank<<" Dtot: "<<Dtot<<endl;

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
  int myNrep=int(round(1.0*Nrep/(1.0*nprocs)));
  Nrep=myNrep*nprocs;
  cout <<"rank: "<<rank<<'\t' <<"My reps: "<<myNrep<<" total reps: "<<Nrep<<endl;
  
  double realsmall=1E-14;
  double prevnorm=1.0;
  int Nx0=plist.Nx0;
  
  double dels;
  
  ofstream probfile;
  ofstream ratefile;
  
  ofstream rfile;
  ofstream pavfile;
  
  double space;
  double tval;
  char tname[200];
  char pname[300];
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
  cout <<"rank: "<<rank<<"Initial separation: "<<plist.x0<<" Nx0: "<<Nx0<<endl;
  x0vec[0]=plist.x0;//for k<inf, this will not go to 1 
  
  
  
  double kappa=kact/(4.0*M_PI*bindrad*bindrad);
  cout <<"rank: "<<rank<<"Kact: "<<kact<<" kappa: "<<kappa<<endl;
  double epsilon=plist.eps_scale*Dtot/kappa;
  double tau=epsilon/kappa;
  cout <<"rank: "<<rank<<"epsscale: "<<plist.eps_scale<<" epsilon: "<<epsilon<<" tau: "<<tau<<endl;
  Rmax=bindrad+epsilon;
  cout <<"rank: "<<rank<<"Bindrad: "<<bindrad<<" reaction limit: "<<Rmax<<endl;
  //  cout <<"rank: "<<rank<<"different time step expansions: "<<epsilon*epsilon/Dtot/100.0<<" other: "<<Rmax*Rmax*0.0001/(2.0*Dtot)<<endl;
  double scaled=plist.dt_scale;
  double deltat_reac=scaled*epsilon*epsilon/(2.0*Dtot);

  double cf=cos(sqrt(4.0*wholep[1].Drx*deltat_reac));
  //  double Dr1=2.0*leglen2*(1.0-cf);
  /*Do not add in rotation to Deff for comparison, since it does not change the distance between how
    far they move.
  */
  //Dtot+=2.0*Dr1/(6.0*deltat_reac);//add in rotation for both molecules
  int ee;
  /*Now update to new Dtot*/
  epsilon=plist.eps_scale*Dtot/kappa;
  tau=epsilon/kappa;
  Rmax=bindrad+epsilon;
  deltat_reac=scaled*epsilon*epsilon/(2.0*Dtot);
  double deltat=deltat_reac;
  
  //x0vec[1]=plist.x0+epsilon;
  for(i=1;i<Nx0;i++)
    x0vec[i]=plist.x0+0.05*i;

  cout <<"Values of X0: ";
  for(i=0;i<Nx0;i++)
    cout <<x0vec[i]<<'\t';
  cout <<endl;
  kdiff=fourpi*Dtot*bindrad;
  fact=1.0+kact/kdiff;
  alpha=fact*sqrt(Dtot)/bindrad;
  double cof=kact/(kact+kdiff);
  
  cout <<"rank: "<<rank<<" deltat: "<<deltat_reac<<" epsilon: "<<epsilon<<" tau: "<<tau<<" Rmax: "<<Rmax<<" Dtot: "<<Dtot<<endl;
  double currtime=0;
  double Maxtime=plist.Maxtime;
  int Nitbin=int(Maxtime/deltat_reac);
  Maxtime=Nitbin*deltat_reac;
  plist.Maxtime=Maxtime; 
  cout <<"rank: "<<rank<<"number of time bins: "<<Nitbin<<" new max time: "<<Maxtime<<endl;
  double pirrev, pfree;
  //current statwrite is the number of datapoints it will write out.
  if(Nitbin<statwrite)
    statwrite=500;
  else
    statwrite=int(round(Nitbin/statwrite));
  
  cout <<"rank: "<<rank<<"statwrite: "<<statwrite<<endl;
  plist.statwrite=statwrite;
  /*Establish limits to orientational association.*/
  double anglelimit=plist.costheta_range;
  double psi_limit=plist.psi_range;
  
  cout <<"Angular limits: "<<anglelimit<<" psi_limit: "<<psi_limit<<endl;
  
  double theta;
  int Nthet=plist.Nthet;
  int Nthet_split=(Nthet-1)*8.0+1;
  double deltheta;
  double intervals;
  //deltheta=anglelimit/(1.0*(Nthet));
  intervals=2.0/(1.0*(Nthet-1));//Here we use -1 because we are not integrating.
  deltheta=2.0/(1.0*Nthet);
  cout <<"rank: "<<rank<<"Ntheta: "<<Nthet<<" theta intervals: "<<intervals<<" delthet: "<<deltheta<<endl;
  double *costheta=new double[Nthet];
  for(i=0;i<Nthet;i++){
    costheta[i]=1-i*intervals;
    cout <<"rank: "<<rank<<"costheta: "<<costheta[i]<<endl;
  }
  
  int Nwrite=int(ceil(Nitbin/statwrite));
  cout <<"rank: "<<rank<<"rank: "<<rank<<'\t' <<"statwrite: "<<statwrite<<" Nwrite bins, for S(t): "<<Nwrite<<endl;
  if((Nitbin-1)==Nwrite*statwrite)
    Nitbin-=1;
  double *phist=new double[Nwrite+1];
  double *stot=new double[Nwrite+1];
  double *ktime=new double[Nwrite+1];
  double *phist_sum=new double[Nwrite+1];
  for(i=0;i<Nwrite+1;i++)
    phist[i]=0.0;
    
  
  int cnt;
  
  int n;
  int itmove;
  int itmax;
  double lz, ly, lx;
  double tx, ty, tz;
  
  double Rrange=3.0*sqrt(6.0*Dtot*Maxtime);//no +bindrad here because only measure distance beyond sigma.
  
  double delR=0.01;
  int Rbins=int(Rrange/delR);
  if(Rbins>1000){
    Rbins=1000;
    delR=(Rrange)/(1.0*Rbins);
  }
  cout <<"rank: "<<rank<<"rank: "<<rank<<'\t' <<"MaxR: "<<Rrange<<" Rbins: "<<Rbins<<" delR: "<<delR<<endl;
    
  int rind;
  
  double small=1E-9;
  

  int thetbins_out=plist.thetbins_out;
  int thetbins2_out=thetbins_out*thetbins_out;

  double delthet_out=2.0/(1.0*thetbins_out)+small;
  cout <<"rank: "<<rank<<"theta bins: "<<thetbins_out<<" delthet out: "<<delthet_out<<endl;
  ofstream thetfile("thetas.out");
  for(i=0;i<thetbins_out;i++)
    thetfile <<(i+0.5)*delthet_out-1<<endl;
  
  int dihbins_out=plist.dihbins_out;
  double deldih_out=M_PI/(1.0*dihbins_out)+small;
  cout <<"rank: "<<rank<<"dihedral bins: "<<dihbins_out<< " delta: "<<deldih_out<<endl;
  int psibins_out=plist.psibins_out;
  double delpsi_out=1.0*M_PI/(1.0*psibins_out)+small;
  int psiind;
  int k1, k2, j1, j2;
  double *prhist=new double[thetbins2_out*dihbins_out*psibins_out*Rbins];//
  double *prhist_sum=new double[thetbins2_out*dihbins_out*psibins_out*Rbins];//
  for(i=0;i<thetbins2_out;i++){
    for(k=0;k<dihbins_out;k++){
      for(m=0;m<psibins_out;m++){
	for(j=0;j<Rbins;j++)
	  prhist[i*dihbins_out*psibins_out*Rbins+k*psibins_out*Rbins+m*Rbins+j]=0.0;
      }
    }
  }
  double *p1hist=new double[Rbins];
  double *p1hist_sum=new double[Rbins];
  for(i=0;i<Rbins;i++)
    p1hist[i]=0.0;

  
  double *MA=new double[9];
  double *M=new double[9];
  double *MB=new double[9];
  double *v3=new double[3];
  double *v2=new double[3];
  double *v1=new double[3];
  double *v=new double[3];
  
  double *c0toc1=new double[3];
  double *c1toim=new double[3];
  double *c1toib=new double[3];
  double *c0toim=new double[3];
  double *c0toib=new double[3];
  double *rotv=new double[3];
  double *p0n=new double[3];
  double *p1n=new double[3];
  
  //declare more vars
  double rabx, raby, rabz;
  double rab2, rab;
  
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

  double psimol, psi;
  //ofstream matchfile;
  int wcnt;
  int ct;
  int bt;
  int dd, cc, hh;
  int Ndih=plist.Ndih;
  double dih_interval=M_PI/(1.0*(Ndih-1));
  double *dih_interval_vec=new double[Ndih];
  double deldih=M_PI/(1.0*(Ndih-1));//dihedral relative to azimuths only range over pi.
  double *dihvec=new double[Ndih];
  
  for(i=0;i<Ndih;i++){
    dihvec[i]=plist.dihedstart+i*dih_interval;
    dih_interval_vec[i]=dih_interval;
    
  }
  /*First and last bins are half the size, at the end points*/
  dih_interval_vec[0]=dih_interval/2.0/M_PI;
  dih_interval_vec[Ndih-1]=dih_interval/2.0;
  for(i=0;i<Ndih;i++){
    if(rank==0)
      cout <<"Dih value: "<<dihvec[i]<<" intervals size: "<<dih_interval_vec[i]<<endl;
  }
  int Npsibins=plist.Npsibins; 
  int nsplit=8;
  int Npsi_split=(Npsibins*0.5-1)*1.0*nsplit+1;//each interval can be split into eighths, plus the last point.
  /*for integrating over initial angle values, the del_interval value should multiply each point,
    AND end points should be multiplied by a facto of a half, because their intervals are half the size of interior points
    to treat all points equally, it is better to use delpsi than the actual interval size.  
  */
  double delpsi_interval=1*M_PI/(1.0*(Npsibins*0.5-1));
  double delpsi=2*M_PI/(1.0*(Npsibins));
  //if(abs(M_PI-psi_limit)<1E-10)
  //Npsibins-=1;//don't do pi twice 
  cout <<"Delpsi: "<<delpsi<<" delpsi_interval: "<<delpsi_interval<<endl;  
  double *psivec=new double[Npsibins];
  for(i=0;i<Npsibins;i+=2){
    // psivec[i]=-M_PI+i*delpsi;
//     if(psivec[i]<0)psivec[i]+=2*M_PI;
    psivec[i]=i/2*delpsi_interval;
    psivec[i+1]=2*M_PI-i/2*delpsi_interval;
    if(i==0)
      psivec[i+1]=0;//don't use 2pi
    cout <<"psivec: "<<psivec[i]<<' '<<psivec[i+1]<<endl;
  }
  
  double currnormpsi=0.0;
  ofstream jacfile;
  if(rank==0){
    jacfile.open("psi_jacobians.dat");
  }
  int Nphi=plist.Nphi;
  double phi1;
  /*Phi values always range over 2 pi because this is a free DOF but it needs to be sampled
    becasue it changes the orientation of molecules and could affect the dynamics.  
  */
  
  double delphi1=2.0*M_PI/(1.0*Nphi);//no minus one becasue we don't want to do both 0 and 2*pi.
  cout <<"Nphi: "<<Nphi<<" delphi: "<<delphi1<<endl;
  double n1_x, n1_y, n1_z;
  double psi0;
  double *pn1=new double[3];
  double *par1=new double[3];
  double tol=1E-7;
  double lx_0, ly_0, lz_0;
  double psi1;
  double parcheck;
  double ctab_0, ctab1;
  double lR;
  double dxm, dym, dzm;
  double p1x, p1y, p1z;
  double h0, h1, c1cof, c0cof;
  double n1n2;
  double p1len;
  double nx_0, ny_0, nz_0;
  double nx_1, ny_1, nz_1;
  int cts=0;
  //double norm_psi[Npsibins];
  double hist_add;
  double psi_jac;
  int ThetEnd=Nthet;
  double dihed;
  int molfix, molrot;
  double phi2;
  int bothperp;
  //double dihlimit;//close to pi, but is more lenient for oriented legs.
  ofstream sorgfile;
  ofstream prrfile;
  double *prr=new double[Rbins];
  double totrnorm=0;
  ofstream sonlyfile;
  double totpsinorm;
  int Ndihsim;

  int cit=0;
  int zerotest_bins=40;
  int zeroflag;
  int splitflag;
  double sprob_tol=plist.sprob_tol;//if the survival prob changes by more than this, split the interval.
  int Nsprob=Npsi_split*Nthet_split*Nthet_split;
  double *sprobvec=new double[Nsprob];
  double *norm_psi=new double[Nsprob];
  for(i=0;i<Nsprob;i++)
    sprobvec[i]=-1;
  int ct1ind, ct2ind, psind0; 
  double cost4, cost2, cost6, cost7, cost5, cost3, cost1, cost8, cost0; 
  int sindex0, sindex1, sindex2, sindex3, sindex4, sindex5, sindex6, sindex7, sindex8;
  int ctindex1[8];
  int ctindex2[8];
  int splitflagct2, splitflagct1;
  double psimid, psiquart, change;
  double costheta1, costheta2;
  double costheta_next, costheta_half, costheta_quart, costheta_eighth;
  int ct1start;
  Ndihsim=Ndih;//just loop over all dihedral, even for 1, 1, so they all have same statistics.
  
  int sxloop;
  Nx0=plist.Nx0;
  for(sxloop=0;sxloop<Nx0;sxloop++){
    x0=x0vec[sxloop];
    
    if(rank==0)cout <<"x0: "<<x0<<" s: "<<sxloop<<endl;
    for(i=0;i<Nwrite+1;i++)
      stot[i]=0;
    totrnorm=0;
    for(i=0;i<Rbins;i++)
      prr[i]=0.0;
    
    if(rank==0){
      sprintf(tname, "survival_vs_angles_x0_%g_dt%g.dat",x0, deltat_reac);
      sonlyfile.open(tname);
      
      sprintf(tname, "Order_survival_vs_angles_x0_%g_dt%g.dat",x0, deltat_reac);
      sorgfile.open(tname);
      
      sprintf(pname, "prr_x0_%g_dt%g.dat", x0,  deltat_reac);
      prrfile.open(pname);
      
      sprintf(pname, "pav_x0_%g_dt%g.dat", x0,  deltat_reac);
      pavfile.open(pname);
      pavfile<<0<<'\t'<<0<<'\t'<<0<<'\t';
      for(i=0;i<Rbins;i++){
	R1=bindrad+delR*(i+0.5);
	pavfile<<R1<<'\t';
      }
      pavfile<<endl;
      pavfile<<0<<'\t'<<0<<'\t'<<0<<'\t';
      for(i=0;i<Rbins;i++){
	R1=bindrad+delR*(i+0.5);
	pirrev=pirrev_value(R1, x0, Maxtime,  Dtot, bindrad,  alpha);
	pavfile<<pirrev<<'\t';
      }
      pavfile<<endl;
    }
    /*loop here over dihedral as well*/
    for(ct=0;ct<Nthet;ct++){
      lz=-costheta[ct]*leglen1;//switched sign from COM offset program, now COM-COM distance is farthest
      lx=sqrt(leglen1_sq-lz*lz);
      ly=0;
      ct1ind=ct*8;
      costheta1=costheta[ct];
      splitflagct1=1;
      cout <<"starting pos: "<<lx<<'\t'<<ly<<'\t'<<lz<<endl;
      while(splitflagct1==1){
	/*For the adaptive version, we need to add values of costheta if the previous value and the current value are too different!*/
	loop_over_costheta2_gridpt( rank,  nprocs,psivec, costheta1,costheta, norm_psi,sprobvec,  Ndihsim, Nphi, x0,  lx, ly,  lz, leglen1, leglen0,plist, deltat_reac,  scaled,  bindrad, bases,wholep, dihvec,  myNrep,  Rmax, tau, anglelimit, psi_limit,  Nitbin,  Nwrite, phist, phist_sum,p1hist, p1hist_sum,prhist, prhist_sum, rfile, probfile, pavfile, jacfile,sonlyfile, tol,  ct1ind, ct, dih_interval_vec );
	if(ct==0)splitflagct1=0;
	else{
	  /*split the intervals here*/
	  ct2ind=ct1ind;
	  ctindex1[8]=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+0;//this is your current bin
	  ctindex1[0]=(ct-1)*8*Nthet_split*Npsi_split+(ct-1)*8*Npsi_split+0;//this is previous diag
	  change=sprobvec[ctindex1[8]]-sprobvec[ctindex1[0]];//diagonal entries
	  //double sprob_tmp=0.01;
	  if(abs(change)>sprob_tol){
	    /*split, start with diagonal*/
	    ct1ind=(ct-1)*8+4;//half interval
	    ct2ind=ct1ind;
	    costheta1=(costheta[ct-1]-costheta[ct])*0.5+costheta[ct];
	    lz=-costheta1*leglen1;//switched sign from COM offset program, now COM-COM distance is farthest
	    lx=sqrt(leglen1_sq-lz*lz);
	    ly=0;
	    
	    cost4=costheta1;
	    cout <<"add point 4: "<<cost4<<" change: "<<change<<" ctindex1: "<<ctindex1[8]<<'\t'<<sprobvec[ctindex1[8]]<<" index2; "<<ctindex1[0]<<" sprob: "<<sprobvec[ctindex1[0]]<<endl;
	    cost8=costheta[ct];
	    cost0=costheta[ct-1];//higher number
	    costheta_half=costheta[ct];
	    loop_over_costheta2_halfgridpt( rank,  nprocs, psivec,  costheta1,  costheta_half, norm_psi,sprobvec, Ndihsim, Nphi, x0, lx,  ly,  lz,  leglen1, leglen0, plist,  deltat_reac, scaled, bindrad, bases, wholep,dihvec, myNrep,  Rmax,  tau, anglelimit, psi_limit, Nitbin, Nwrite,phist, phist_sum,p1hist,p1hist_sum,prhist, prhist_sum, rfile,probfile,pavfile, jacfile, sonlyfile, tol,  ct1ind, ct1ind, dih_interval_vec);//first do new diagaonl, to next grid point
	    /*Now do new costheta 1, with rest of grid points.*/
	    loop_over_costheta2_gridpt( rank,  nprocs,psivec, costheta1,costheta, norm_psi,sprobvec,  Ndihsim, Nphi, x0,  lx, ly,  lz, leglen1, leglen0,plist, deltat_reac,  scaled,  bindrad, bases,wholep, dihvec,  myNrep,  Rmax, tau, anglelimit, psi_limit,  Nitbin,  Nwrite, phist, phist_sum,p1hist, p1hist_sum,prhist, prhist_sum, rfile, probfile, pavfile, jacfile,sonlyfile, tol,  ct1ind, ct+1, dih_interval_vec );//ct+1 is because in half grid you do the current value of ct.
	    /*Now test left and right of this interval.*/
	    ctindex1[4]=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+0;//diagonal entry.
	    /*left first*/
	    change=sprobvec[ctindex1[4]]-sprobvec[ctindex1[0]];//diagonal entries
	    if(abs(change)>sprob_tol){
	      /*add point 2.*/
	      ct1ind=(ct-1)*8+2;//half interval
	      ct2ind=ct1ind;
	      costheta1=(cost0-cost4)*0.5+cost4;
	      lz=-costheta1*leglen1;//switched sign from COM offset program, now COM-COM distance is farthest
	      lx=sqrt(leglen1_sq-lz*lz);
	      ly=0;
	      cout <<"add point 2: "<<cost2<<" change: "<<change<<" ctindex1: "<<ctindex1[4]<<'\t'<<sprobvec[ctindex1[4]]<<" index2; "<<ctindex1[0]<<" sprob: "<<sprobvec[ctindex1[0]]<<endl;
	      cost2=costheta1;
	      costheta2=cost2;//start off with diagonal.
	      costheta_quart=cost4;
	      ctindex1[2]=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+0;//diagonal entry.
	      ct1start=ct1ind;
	      loop_over_costheta2_quartergridpt( rank,  nprocs, psivec,  costheta1,  costheta_quart, norm_psi,sprobvec, Ndihsim, Nphi, x0, lx,  ly,  lz,  leglen1, leglen0, plist,  deltat_reac, scaled, bindrad, bases, wholep,dihvec, myNrep,  Rmax,  tau, anglelimit, psi_limit, Nitbin, Nwrite,phist, phist_sum,p1hist,p1hist_sum,prhist, prhist_sum, rfile,probfile,pavfile, jacfile, sonlyfile, tol,  ct1ind, ct1start, dih_interval_vec);//first do new diagaonl, to next grid point
	      ct1start=ct1ind+2;
	      loop_over_costheta2_halfgridpt( rank,  nprocs, psivec,  costheta1,  costheta_half, norm_psi,sprobvec, Ndihsim, Nphi, x0, lx,  ly,  lz,  leglen1, leglen0, plist,  deltat_reac, scaled, bindrad, bases, wholep,dihvec, myNrep,  Rmax,  tau, anglelimit, psi_limit, Nitbin, Nwrite,phist, phist_sum,p1hist,p1hist_sum,prhist, prhist_sum, rfile,probfile,pavfile, jacfile, sonlyfile, tol,  ct1ind, ct1start, dih_interval_vec);//first do new diagaonl, to next grid point
	      /*Now do new costheta 1, with rest of grid points.*/
	      loop_over_costheta2_gridpt( rank,  nprocs,psivec, costheta1,costheta, norm_psi,sprobvec,  Ndihsim, Nphi, x0,  lx, ly,  lz, leglen1, leglen0,plist, deltat_reac,  scaled,  bindrad, bases,wholep, dihvec,  myNrep,  Rmax, tau, anglelimit, psi_limit,  Nitbin,  Nwrite, phist, phist_sum,p1hist, p1hist_sum,prhist, prhist_sum, rfile, probfile, pavfile, jacfile,sonlyfile, tol,  ct1ind, ct+1, dih_interval_vec );//ct+1 is because in half grid you do the current value of ct.
	      /*now split this interval, add points 1 and 3*/
	      change=sprobvec[ctindex1[2]]-sprobvec[ctindex1[0]];//diagonal entries
	      if(abs(change)>sprob_tol){
		/*add point 1*/
		ct1ind=(ct-1)*8+1;
		ct2ind=ct1ind;
		costheta1=(cost0-cost2)*0.5+cost2;
		lz=-costheta1*leglen1;//switched sign from COM offset program, now COM-COM distance is farthest
		lx=sqrt(leglen1_sq-lz*lz);
		ly=0;
		
		cost1=costheta1;
		cout <<"add point 1: "<<cost1<<" change: "<<change<<" ctindex1: "<<ctindex1[2]<<'\t'<<sprobvec[ctindex1[2]]<<" index2; "<<ctindex1[0]<<" sprob: "<<sprobvec[ctindex1[0]]<<endl;
		costheta2=costheta1;//start off diagonal
		costheta_eighth=cost2;
		loop_over_costheta2_eighthgridpt( rank,  nprocs, psivec,  costheta1,  costheta_eighth, norm_psi,sprobvec, Ndihsim, Nphi, x0, lx,  ly,  lz,  leglen1, leglen0, plist,  deltat_reac, scaled, bindrad, bases, wholep,dihvec, myNrep,  Rmax,  tau, anglelimit, psi_limit, Nitbin, Nwrite,phist, phist_sum,p1hist,p1hist_sum,prhist, prhist_sum, rfile,probfile,pavfile, jacfile, sonlyfile, tol,  ct1ind, dih_interval_vec);//first do new diagaonl, to next grid point
		ct1start=ct1ind+1;
		loop_over_costheta2_quartergridpt( rank,  nprocs, psivec,  costheta1,  costheta_quart, norm_psi,sprobvec, Ndihsim, Nphi, x0, lx,  ly,  lz,  leglen1, leglen0, plist,  deltat_reac, scaled, bindrad, bases, wholep,dihvec, myNrep,  Rmax,  tau, anglelimit, psi_limit, Nitbin, Nwrite,phist, phist_sum,p1hist,p1hist_sum,prhist, prhist_sum, rfile,probfile,pavfile, jacfile, sonlyfile, tol,  ct1ind, ct1start, dih_interval_vec);//first do new diagaonl, to next grid point
		ct1start=ct1ind+3;
		loop_over_costheta2_halfgridpt( rank,  nprocs, psivec,  costheta1,  costheta_half, norm_psi,sprobvec, Ndihsim, Nphi, x0, lx,  ly,  lz,  leglen1, leglen0, plist,  deltat_reac, scaled, bindrad, bases, wholep,dihvec, myNrep,  Rmax,  tau, anglelimit, psi_limit, Nitbin, Nwrite,phist, phist_sum,p1hist,p1hist_sum,prhist, prhist_sum, rfile,probfile,pavfile, jacfile, sonlyfile, tol,  ct1ind, ct1start, dih_interval_vec);//first do new diagaonl, to next grid point
		/*Now do new costheta 1, with rest of grid points.*/
		loop_over_costheta2_gridpt( rank,  nprocs,psivec, costheta1,costheta, norm_psi,sprobvec,  Ndihsim, Nphi, x0,  lx, ly,  lz, leglen1, leglen0,plist, deltat_reac,  scaled,  bindrad, bases,wholep, dihvec,  myNrep,  Rmax, tau, anglelimit, psi_limit,  Nitbin,  Nwrite, phist, phist_sum,p1hist, p1hist_sum,prhist, prhist_sum, rfile, probfile, pavfile, jacfile,sonlyfile, tol,  ct1ind, ct+1, dih_interval_vec );//ct+1 is because in half grid you do the current value of ct.
		
	      }
	      /*now left*/
	      change=sprobvec[ctindex1[2]]-sprobvec[ctindex1[4]];//diagonal entries
	      if(abs(change)>sprob_tol){
		/*add point 3*/
		ct1ind=(ct-1)*8+3;
		ct2ind=ct1ind;
		costheta1=(cost2-cost4)*0.5+cost4;
		lz=-costheta1*leglen1;//switched sign from COM offset program, now COM-COM distance is farthest
		lx=sqrt(leglen1_sq-lz*lz);
		ly=0;
	    
		cost3=costheta1;
		cout <<"add point 3: "<<cost3<<" change: "<<change<<" ctindex1: "<<ctindex1[2]<<'\t'<<sprobvec[ctindex1[2]]<<" index2; "<<ctindex1[4]<<" sprob: "<<sprobvec[ctindex1[4]]<<endl;
		costheta2=costheta1;//start off diagonal
		//costheta_eighth=cost2;
		loop_over_costheta2_eighthgridpt( rank,  nprocs, psivec,  costheta1,  costheta_quart, norm_psi,sprobvec, Ndihsim, Nphi, x0, lx,  ly,  lz,  leglen1, leglen0, plist,  deltat_reac, scaled, bindrad, bases, wholep,dihvec, myNrep,  Rmax,  tau, anglelimit, psi_limit, Nitbin, Nwrite,phist, phist_sum,p1hist,p1hist_sum,prhist, prhist_sum, rfile,probfile,pavfile, jacfile, sonlyfile, tol,  ct1ind, dih_interval_vec);//first do new diagaonl, to next grid point
		ct1start=ct1ind+1;
		loop_over_costheta2_halfgridpt( rank,  nprocs, psivec,  costheta1,  costheta_half, norm_psi,sprobvec, Ndihsim, Nphi, x0, lx,  ly,  lz,  leglen1, leglen0, plist,  deltat_reac, scaled, bindrad, bases, wholep,dihvec, myNrep,  Rmax,  tau, anglelimit, psi_limit, Nitbin, Nwrite,phist, phist_sum,p1hist,p1hist_sum,prhist, prhist_sum, rfile,probfile,pavfile, jacfile, sonlyfile, tol,  ct1ind, ct1start, dih_interval_vec);//first do new diagaonl, to next grid point
		/*Now do new costheta 1, with rest of grid points.*/
		loop_over_costheta2_gridpt( rank,  nprocs,psivec, costheta1,costheta, norm_psi,sprobvec,  Ndihsim, Nphi, x0,  lx, ly,  lz, leglen1, leglen0,plist, deltat_reac,  scaled,  bindrad, bases,wholep, dihvec,  myNrep,  Rmax, tau, anglelimit, psi_limit,  Nitbin,  Nwrite, phist, phist_sum,p1hist, p1hist_sum,prhist, prhist_sum, rfile, probfile, pavfile, jacfile,sonlyfile, tol,  ct1ind, ct+1, dih_interval_vec );//ct+1 is because in half grid you do the current value of ct.
		
	      }
	      
	      
	    }
	    /*now right*/
	    change=sprobvec[ctindex1[4]]-sprobvec[ctindex1[8]];//diagonal entries
	    if(abs(change)>sprob_tol){
	      /*add point 6.*/
	      ct1ind=(ct-1)*8+6;//half interval
	      ct2ind=ct1ind;
	      costheta1=(cost4-cost8)*0.5+cost8;
	      lz=-costheta1*leglen1;//switched sign from COM offset program, now COM-COM distance is farthest
	      lx=sqrt(leglen1_sq-lz*lz);
	      ly=0;
	      
	      cost6=costheta1;
	      ctindex1[6]=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+0;//diagonal entry.
	      cout <<"add point 6: "<<cost6<<" change: "<<change<<" ctindex1: "<<ctindex1[8]<<'\t'<<sprobvec[ctindex1[8]]<<" index2; "<<ctindex1[4]<<" sprob: "<<sprobvec[ctindex1[4]]<<endl;
	      costheta2=cost6;//start off with diagonal.
	      costheta_quart=cost6;
	      ct1start=ct1ind;
	      loop_over_costheta2_quartergridpt( rank,  nprocs, psivec,  costheta1,  costheta_half, norm_psi,sprobvec, Ndihsim, Nphi, x0, lx,  ly,  lz,  leglen1, leglen0, plist,  deltat_reac, scaled, bindrad, bases, wholep,dihvec, myNrep,  Rmax,  tau, anglelimit, psi_limit, Nitbin, Nwrite,phist, phist_sum,p1hist,p1hist_sum,prhist, prhist_sum, rfile,probfile,pavfile, jacfile, sonlyfile, tol,  ct1ind, ct1start, dih_interval_vec);//first do new diagaonl, to next grid point
	      //ct1start=ct1ind+2;
	      // loop_over_costheta2_halfgridpt( rank,  nprocs, psivec,  costheta1,  costheta_half, norm_psi,sprobvec, Ndihsim, Nphi, x0, lx,  ly,  lz,  leglen1, leglen0, plist,  deltat_reac, scaled, bindrad, bases, wholep,dihvec, myNrep,  Rmax,  tau, anglelimit, psi_limit, Nitbin, Nwrite,phist, phist_sum,p1hist,p1hist_sum,prhist, prhist_sum, rfile,probfile,pavfile, jacfile, sonlyfile, tol,  ct1ind, ct1start);//first do new diagaonl, to next grid point
	      /*Now do new costheta 1, with rest of grid points.*/
	      loop_over_costheta2_gridpt( rank,  nprocs,psivec, costheta1,costheta, norm_psi,sprobvec,  Ndihsim, Nphi, x0,  lx, ly,  lz, leglen1, leglen0,plist, deltat_reac,  scaled,  bindrad, bases,wholep, dihvec,  myNrep,  Rmax, tau, anglelimit, psi_limit,  Nitbin,  Nwrite, phist, phist_sum,p1hist, p1hist_sum,prhist, prhist_sum, rfile, probfile, pavfile, jacfile,sonlyfile, tol,  ct1ind, ct+1, dih_interval_vec );//ct+1 is because in half grid you do the current value of ct.
	      /*split this, add points 5 and 7*/
	      change=sprobvec[ctindex1[6]]-sprobvec[ctindex1[4]];//diagonal entries
	      if(abs(change)>sprob_tol){
		/*add point 5*/
		ct1ind=(ct-1)*8+5;
		ct2ind=ct1ind;
		costheta1=(cost4-cost6)*0.5+cost6;
		lz=-costheta1*leglen1;//switched sign from COM offset program, now COM-COM distance is farthest
		lx=sqrt(leglen1_sq-lz*lz);
		ly=0;
	    
		cost5=costheta1;
		cout <<"add point 5: "<<cost5<<" change: "<<change<<" ctindex1: "<<ctindex1[6]<<'\t'<<sprobvec[ctindex1[6]]<<" index2; "<<ctindex1[4]<<" sprob: "<<sprobvec[ctindex1[4]]<<endl;
		costheta2=costheta1;//start off diagonal
		costheta_eighth=cost6;
		loop_over_costheta2_eighthgridpt( rank,  nprocs, psivec,  costheta1,  costheta_eighth, norm_psi,sprobvec, Ndihsim, Nphi, x0, lx,  ly,  lz,  leglen1, leglen0, plist,  deltat_reac, scaled, bindrad, bases, wholep,dihvec, myNrep,  Rmax,  tau, anglelimit, psi_limit, Nitbin, Nwrite,phist, phist_sum,p1hist,p1hist_sum,prhist, prhist_sum, rfile,probfile,pavfile, jacfile, sonlyfile, tol,  ct1ind, dih_interval_vec);//first do new diagaonl, to next grid point
		ct1start=ct1ind+1;
		loop_over_costheta2_quartergridpt( rank,  nprocs, psivec,  costheta1,  costheta_half, norm_psi,sprobvec, Ndihsim, Nphi, x0, lx,  ly,  lz,  leglen1, leglen0, plist,  deltat_reac, scaled, bindrad, bases, wholep,dihvec, myNrep,  Rmax,  tau, anglelimit, psi_limit, Nitbin, Nwrite,phist, phist_sum,p1hist,p1hist_sum,prhist, prhist_sum, rfile,probfile,pavfile, jacfile, sonlyfile, tol,  ct1ind, ct1start, dih_interval_vec);//first do new diagaonl, to next grid point
		
		/*Now do new costheta 1, with rest of grid points.*/
		loop_over_costheta2_gridpt( rank,  nprocs,psivec, costheta1,costheta, norm_psi,sprobvec,  Ndihsim, Nphi, x0,  lx, ly,  lz, leglen1, leglen0,plist, deltat_reac,  scaled,  bindrad, bases,wholep, dihvec,  myNrep,  Rmax, tau, anglelimit, psi_limit,  Nitbin,  Nwrite, phist, phist_sum,p1hist, p1hist_sum,prhist, prhist_sum, rfile, probfile, pavfile, jacfile,sonlyfile, tol,  ct1ind, ct+1, dih_interval_vec );//ct+1 is because in half grid you do the current value of ct.
		
	      }
	      change=sprobvec[ctindex1[6]]-sprobvec[ctindex1[8]];//diagonal entries
	      if(abs(change)>sprob_tol){
		/*add point 7*/
		ct1ind=(ct-1)*8+7;
		ct2ind=ct1ind;
		costheta1=(cost6-cost8)*0.5+cost8;
		lz=-costheta1*leglen1;//switched sign from COM offset program, now COM-COM distance is farthest
		lx=sqrt(leglen1_sq-lz*lz);
		ly=0;
	    
		cost7=costheta1;
		cout <<"add point 7: "<<cost7<<" change: "<<change<<" ctindex1: "<<ctindex1[8]<<'\t'<<sprobvec[ctindex1[8]]<<" index2; "<<ctindex1[6]<<" sprob: "<<sprobvec[ctindex1[6]]<<endl;
		costheta2=costheta1;//start off diagonal
		costheta_eighth=cost6;
		loop_over_costheta2_eighthgridpt( rank,  nprocs, psivec,  costheta1,  costheta_half, norm_psi,sprobvec, Ndihsim, Nphi, x0, lx,  ly,  lz,  leglen1, leglen0, plist,  deltat_reac, scaled, bindrad, bases, wholep,dihvec, myNrep,  Rmax,  tau, anglelimit, psi_limit, Nitbin, Nwrite,phist, phist_sum,p1hist,p1hist_sum,prhist, prhist_sum, rfile,probfile,pavfile, jacfile, sonlyfile, tol,  ct1ind, dih_interval_vec);//first do new diagaonl, to next grid point
				
		/*Now do new costheta 1, with rest of grid points.*/
		loop_over_costheta2_gridpt( rank,  nprocs,psivec, costheta1,costheta, norm_psi,sprobvec,  Ndihsim, Nphi, x0,  lx, ly,  lz, leglen1, leglen0,plist, deltat_reac,  scaled,  bindrad, bases,wholep, dihvec,  myNrep,  Rmax, tau, anglelimit, psi_limit,  Nitbin,  Nwrite, phist, phist_sum,p1hist, p1hist_sum,prhist, prhist_sum, rfile, probfile, pavfile, jacfile,sonlyfile, tol,  ct1ind, ct+1, dih_interval_vec );//ct+1 is because in half grid you do the current value of ct.
		
	      }
	      
	      
	    }
	  }
	  
	  splitflagct1=0;
	}
      
      }
      
      //each intergral over cos(theta) gives factor of 2.
    }//end first angle
    //integral over dihedral is over 2 azimuths, so (2*pi)^2, but we average over one full azimuth (2pi), and then only consider difference for other one so up to pi.
    //so dihedral gives norm fact of pi.
    
    
    double dsmall=intervals/(1.0*nsplit);
    double dpsismall=delpsi_interval/(1.0*nsplit);
    if(rank==0){
      /*write out organized survival probs.*/
      sorgfile<<000<<'\t'<<000<<'\t';
      for(k=0;k<Npsi_split;k++)
	sorgfile<<k*dpsismall<<'\t';
      sorgfile<<endl;
      for(i=0;i<Nthet_split;i++){
	for(j=i;j<Nthet_split;j++){
	  sorgfile<<1-i*dsmall<<'\t'<<1-j*dsmall<<'\t';
	  for(k=0;k<Npsi_split;k++){
	    sorgfile<<sprobvec[i*Nthet_split*Npsi_split+j*Npsi_split+k]<<'\t';
	  }
	  sorgfile<<endl;
	}
      }
      
    }
    
    // if(rank==0){
//       for(i=0;i<Rbins;i++){
// 	R1=bindrad+delR*(i+0.5);
// 	prrfile<<R1<<'\t'<<prr[i]/(1.0*totrnorm)<<endl;
//       }
//     }
    MPI::COMM_WORLD.Barrier();
    //  prrfile.close();
    sorgfile.close();
    sonlyfile.close();
    /*Add up survival probs given initial angles, weighted by appropriate Jacobian factor, for angular ocnditions where kappa is non zero.
      divide by full range of angular degrees of freedom possible (for dihedral just use pi) and multiply by 4*pi*sigma^2*kappa
      to get the rate constant.
      this will be the actual rate at sigma, so average it with the rate at epsilon.
    */
    // sprintf(tname, "rate_vs_time_x0_%g_dt%g.dat",x0, deltat_reac);
//     if(rank==0){
      
//       ratefile.open(tname);
//       for(i=statwrite;i<Nitbin;i+=statwrite){
// 	tval=i*deltat_reac;
// 	ind=i/statwrite;
// 	/*pi is from dihedral, although rate is independent of dihedral, survival would depend on this orientation, but 
// 	  only over range of pi, further factor of 2pi and 2 are redundant.
// 	  2*2 is from integrals over cos(theta). 4pi^2 is from integral of psi and phi, each has range of 2pi.
// 	  4pisigma^2*kappa is the integral over r=sigma surface. This vector can rotate fully and particles still align
// 	*/
// 	ktime[ind]=4.0*M_PI*bindrad*bindrad*kappa*stot[ind]/(M_PI*2.0*2.0*4.0*M_PI*M_PI);//normalizing factor over full angular range noted above
// 	ratefile<<tval<<'\t'<<ktime[ind]<<endl;
// 	cout<<" time, rate: "<<tval<<'\t'<<ktime[ind]<<" stot: "<<stot[ind]<<endl;
//       }
//       ratefile.close();
//     }
    cout <<"FINISHED WITH X0: "<<x0<<" curr s: "<<sxloop<<endl;
    MPI::COMM_WORLD.Barrier();
  }//end over all x0s

  
  stop_timer(&totaltime);
  cout <<"rank: "<<rank<<timer_duration(totaltime)<<" total time "<<endl;  
  /*Write out final result*/
  cout <<"rank: "<<rank<<"End Main, complete run "<<endl;
  MPI::COMM_WORLD.Barrier();
  MPI::Finalize();
  
}//end main


void read_parms(ifstream &parmfile, survParms &plist)
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
  parmfile >>plist.leglen1;
  parmfile.ignore(400,'\n');
  parmfile >>plist.leglen0;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Nx0;
  parmfile.ignore(400,'\n');
  parmfile >>plist.x0;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Nthet;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Npsibins;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Nphi;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Ndih;
  parmfile.ignore(400,'\n');
  parmfile >>plist.costheta_range;
  parmfile.ignore(400,'\n');
  parmfile >>plist.psi_range;
  parmfile.ignore(400,'\n');
  parmfile >>plist.bindrad;
  parmfile.ignore(400,'\n');
  parmfile >>plist.dihedstart;
  parmfile.ignore(400,'\n');
  parmfile >>plist.thetbins_out;
  parmfile.ignore(400,'\n');
  parmfile >>plist.dihbins_out;
  parmfile.ignore(400,'\n');
  parmfile >>plist.psibins_out;
  parmfile.ignore(400,'\n');
  parmfile >>plist.sprob_tol;
  parmfile.ignore(400,'\n');

}
void write_crds(Fullmol *bases, int p1)
{
  cout <<bases[p1].xcom<<' '<<bases[p1].ycom<<' '<<bases[p1].zcom<<endl;
  int i;
  for(i=0;i<bases[p1].ninterface;i++)
    cout <<bases[p1].x[i]<<' '<<bases[p1].y[i]<<' '<<bases[p1].z[i]<<endl;

}




void run_sim_at_1angleset(int rank, int nprocs, double psiread, double costheta1, double costheta2, double &norm_psi, double &sprobvalue, int Ndihsim, int Nphi, double x0, double lx, double ly, double lz, double leglen1, double leglen0, survParms &plist, double deltat_reac, double scaled, double bindrad, Fullmol *bases, Protein *wholep, double *dihvec, int myNrep, double Rmax, double tau, double anglelimit, double psi_limit, int Nitbin, int Nwrite, double*phist, double *phist_sum, double *p1hist, double *p1hist_sum,double *prhist, double *prhist_sum, ofstream &rfile, ofstream &probfile, ofstream &pavfile, ofstream &jacfile, ofstream &sonlyfile, double tol, double *dih_intervals)
{
  
  int i, j, k, m;
  int k1, k2, j1, j2;
  double currtime=0;
  double Maxtime=plist.Maxtime;
  double Dtot=plist.D;
  double Rrange=3.0*sqrt(6.0*Dtot*Maxtime);//no +bindrad here because only measure distance beyond sigma.
  double delR=0.01;
  int Rbins=int(Rrange/delR);
  if(Rbins>1000){
    Rbins=1000;
    delR=(Rrange)/(1.0*Rbins);
  }
  int rind;
  double small=1E-9;
  int thetbins_out=plist.thetbins_out;
  int thetbins2_out=thetbins_out*thetbins_out;
  double delthet_out=2.0/(1.0*thetbins_out)+small;
  int dihbins_out=plist.dihbins_out;
  double deldih_out=M_PI/(1.0*dihbins_out)+small;
  int psibins_out=plist.psibins_out;
  double delpsi_out=1.0*M_PI/(1.0*psibins_out)+small;
  int zeroflag;
  int zerotest_bins=40;//to see if it looks like the whole distribution will be zero for the output angles.
  double sprob_final, sprob_finaltot;
  /*initialize arrays that sum over procs, and over phi values*/
  for(i=0;i<Nwrite+1;i++){
    phist_sum[i]=0.0;
    phist[i]=0.0;
  }
  /*initialize arrays that sum over procs, and over phi values*/
  for(i=0;i<Rbins;i++){
    p1hist_sum[i]=0.0;
    p1hist[i]=0.0;
  }
  for(i=0;i<thetbins2_out;i++){
    for(k=0;k<dihbins_out;k++){
      for(m=0;m<psibins_out;m++){
	for(j=0;j<Rbins;j++){
	  prhist_sum[i*dihbins_out*psibins_out*Rbins+k*psibins_out*Rbins+m*Rbins+j]=0.0;
	  prhist[i*dihbins_out*psibins_out*Rbins+k*psibins_out*Rbins+m*Rbins+j]=0.0;
	}
	
      }
    }
  }
  
  double R2, R1;
  double dxm, dym, dzm;
  double dx, dy, dz;
  double r0x, r0y, r0z;
  double r1x, r1y, r1z;
  double n1_x, n1_y, n1_z;
  double currx, curry, currz;
  double leglen0_sq=leglen0*leglen0;
  double leglen1_sq=leglen1*leglen1;
  double deltat;
  double *MA=new double[9];
  double *M=new double[9];
  double *MB=new double[9];
  double *v3=new double[3];
  double *v2=new double[3];
  double *v1=new double[3];
  double *v=new double[3];
  double *n1=new double[3];
  double *n2=new double[3];
  
  
  double *c0toc1=new double[3];
  double *c1toim=new double[3];
  double *c1toib=new double[3];
  double *c0toim=new double[3];
  double *c0toib=new double[3];
  double *rotv=new double[3];
  double *p0n=new double[3];
  double *p1n=new double[3];
  
  double delphi1=2.0*M_PI/(1.0*Nphi);//no minus one becasue we don't want to do both 0 and 2*pi.
  int molfix, molrot;
  double phi1, phi2, dihed, psimol; 
  double lx_0, ly_0, lz_0, lR;
  double currnormpsi;
  double psi_jac, hist_add;
  int rep;
  int Ntotalmol;
  int it, itmove, cit;
  int itmax;
  
  /*for inside BD part*/
  int nfree;
  double prob;
  double cthet1, cthet2, dih, psi_out;
  int ind_thet, ind_thet2;
  int dihbin;
  int psiind;
  double rnum, tval;
  double tx, ty, tz;
  int p1, p2, rxn1;  
  int ind, n;
  int statwrite=plist.statwrite;
  char tname[300];
  int intparm;
  double dubsend;
  /*BEGIN LOOPING OVER BOTH PSI, DIHEDRALS, AND PHI2s */
  int usephi2=0;
  molfix=1;
  if(abs(costheta1)>abs(costheta2))
    molfix=0;//if bt is closer to perp, keep 0 fixed
  
  if(costheta1==0 && costheta2==0){
    usephi2=1;
    molfix=1;
  }
  double dihbinsize;
  double psi0=psiread;
  if(usephi2==1){
    /*in this case only 2 psi values can be sampled.*/
    //  phi2=psivec[dd];//to rotate aroun norm vector regardless.
    if(rank==0){
      sprintf(tname, "prt_sprob_x0_%g_cat%g_cbt%g_PHI%g_dt%g.dat",x0,  costheta1, costheta2, phi2, deltat_reac);
      rfile.open(tname);
      sprintf(tname, "St_x0_%g_ct%g_c2t_%g_PHI%g_dt%g.dat",x0, costheta1, costheta2, phi2, deltat_reac);
      probfile.open(tname);
    }
  }else{
    if(rank==0){
      sprintf(tname, "prt_sprob_x0_%g_cat%g_cbt%g_psi%g_dt%g.dat",x0,  costheta1, costheta2, psi0, deltat_reac);
      rfile.open(tname);
      sprintf(tname, "St_x0_%g_ct%g_c2t_%g_psi%g_dt%g.dat",x0,  costheta1, costheta2, psi0, deltat_reac);
      probfile.open(tname);
    }
  }
  norm_psi=0;
  //totpsinorm=0.0;
  int ee, hh, cc;
  for(ee=0;ee<2;ee++){
    if(ee==1 && psiread!=0 )psi0=2.0*M_PI-psiread;
    
    for(hh=0;hh<Ndihsim;hh++){
      dihed=dihvec[hh];
      dihbinsize=dih_intervals[hh];
      
      lz_0=+costheta2*leglen0;//switched from other program
      
      lR=sqrt(leglen0_sq-lz_0*lz_0);
      lx_0=lR*cos(dihed);
      ly_0=lR*sin(dihed);
      //cout <<"rank: "<<rank<<"Initial pos for leg2 com: "<<lx_0<<' '<<ly_0<<' '<<lz_0<<"  for leg: "<<0<<' '<<0<<' '<<0<<endl;
      //vector to rotator centers
      
      //cout <<"rank: "<<rank<<"Initial leg1 angle to the com divider: "<<ctab1<<" initial leg2 angle to the com divider: "<<ctab_0<<endl;
      
      bases[0].z[0]=lz_0;
      bases[0].x[0]=lx_0;
      bases[0].y[0]=ly_0;
      
      bases[1].z[0]=lz+x0;
      bases[1].x[0]=lx;
      bases[1].y[0]=ly;
      
      bases[0].xcom=0;
      bases[0].ycom=0;
      bases[0].zcom=0;
      
      bases[1].xcom=0.0;
      bases[1].ycom=0.0;
      bases[1].zcom=x0;
      
      if(rank==0){
	cout <<"rank: "<<rank<<" p0 crds, before setting mol orient. "<<endl;
	write_crds(bases, 0);
	cout <<"rank: "<<rank<<" p1 crds, before setting mol orient. "<<endl;
	write_crds(bases, 1);
      
      
	/*Assign orientation of molecules relative to one another (normals, not legs)*/
	
	cout <<"rank: "<<rank<<"target angles: "<<costheta1<<' '<<costheta2<<' '<<dihed<<" target psi: "<<psi0<<endl;
      }
      currnormpsi=0.0;
      sprob_final=0.0;
      if(usephi2==1){
	/*in this case only 2 psi values can be sampled.*/
	phi2=psi0;//to rotate aroun norm vector regardless.
      }
      for(cc=0;cc<Nphi;cc++){
	/*First normal can be selected into the y plane, since the leg is fixed in the x-z plane
	  x and z coordinates are same as COM, length is unit*/
	
	bases[0].z[0]=lz_0;
	bases[0].x[0]=lx_0;
	bases[0].y[0]=ly_0;
	
	bases[1].z[0]=lz+x0;
	bases[1].x[0]=lx;
	bases[1].y[0]=ly;
	
	bases[0].xcom=0;
	bases[0].ycom=0;
	bases[0].zcom=0;
	
	bases[1].xcom=0.0;
	bases[1].ycom=0.0;
	bases[1].zcom=x0;
	
	/*COM to COM vector*/
	dxm=bases[1].xcom-bases[0].xcom;
	dym=bases[1].ycom-bases[0].ycom;
	dzm=bases[1].zcom-bases[0].zcom;
	R2=dxm*dxm+dym*dym+dzm*dzm;
	R1=sqrt(R2);//this should be x0!
	
	c0toc1[0]=dxm/R1;
	c0toc1[1]=dym/R1;
	c0toc1[2]=dzm/R1;
	
	phi1=(cc+0.5)*delphi1;
	
	/*here we might rotate phi around the other molecule first, if perp*/
	if(molfix==1){
	  molrot=0;
	  n1_x=bases[1].xcom+0;
	  n1_y=bases[1].ycom+1;
	  n1_z=bases[1].zcom+0;
	  
	  /*position of the end point of the molecule orientation vector*/
	  bases[1].x[1]=n1_x;
	  bases[1].y[1]=n1_y;
	  bases[1].z[1]=n1_z;
	  
	  /*Vector from COM mol1 to interface for molecule orientation (1)*/
	  c1toim[0]=bases[1].x[1]-bases[1].xcom;
	  c1toim[1]=bases[1].y[1]-bases[1].ycom;
	  c1toim[2]=bases[1].z[1]-bases[1].zcom;
	  /*This vector, c1toim,  should be rotated full 360 to sample all the DOF, around the c1toib vector.*/
	  
	  c1toib[0]=(bases[1].x[0]-bases[1].xcom)/leglen1;
	  c1toib[1]=(bases[1].y[0]-bases[1].ycom)/leglen1;
	  c1toib[2]=(bases[1].z[0]-bases[1].zcom)/leglen1;
	  calc_Rmatrix(c1toib, phi1, M);
	  rotate(c1toim, M, rotv);
	  c1toim[0]=rotv[0];
	  c1toim[1]=rotv[1];
	  c1toim[2]=rotv[2];
	  /*new positions below for molecule orientation end point*/
	  bases[1].x[1]=bases[1].xcom+c1toim[0];
	  bases[1].y[1]=bases[1].ycom+c1toim[1];
	  bases[1].z[1]=bases[1].zcom+c1toim[2];
	  
	  /*now use psi0 to define orientation of bases[0].c0toim, fixed bases[1].c1toim*/
	  if(usephi2==1){
	    /*in this case, only two psis are sampled, rotate phi2 and calc what psi*/
	    psimol=rotate_phi2_mol(phi2, bases, leglen1, leglen0, molfix, molrot, c0toim);
	    cout <<"rotated phi2 instead, calculated psi: "<<psimol<<endl;
	    psi0=psimol;
	  }else
	    psimol=set_norm_to_psi(psi0, bases, leglen1, leglen0, molfix, molrot, c0toim);
	    
	  
	}else{
	  /*rotate around the other mol, first find a starting norm that will be intersection of planes
	    defined by the two norms c0toib and c0toc1.
	  */
	  molfix=0;
	  molrot=1;
	  /*molecule 0 vector to binding interface*/
	  c0toib[0]=(bases[0].x[0]-bases[0].xcom)/leglen0;
	  c0toib[1]=(bases[0].y[0]-bases[0].ycom)/leglen0;
	  c0toib[2]=(bases[0].z[0]-bases[0].zcom)/leglen0;
	  
	  /*COM to COM vector*/
	  dxm=bases[1].xcom-bases[0].xcom;
	  dym=bases[1].ycom-bases[0].ycom;
	  dzm=bases[1].zcom-bases[0].zcom;
	  R2=dxm*dxm+dym*dym+dzm*dzm;
	  R1=sqrt(R2);//this should be x0!
	  
	  
	  c0toc1[0]=dxm/R1;
	  c0toc1[1]=dym/R1;
	  c0toc1[2]=dzm/R1;
	  
	  crossproduct(c0toc1, c0toib, v3);
	  c0toim[0]=v3[0]-bases[0].xcom;
	  c0toim[1]=v3[1]-bases[0].ycom;
	  c0toim[2]=v3[2]-bases[0].zcom;
	  
	  R2=c0toim[0]*c0toim[0]+c0toim[1]*c0toim[1]+c0toim[2]*c0toim[2];
	  R1=sqrt(R2);
	  c0toim[0]/=R1;
	  c0toim[1]/=R1;
	  c0toim[2]/=R1;
	  /*now rotate this molecules vector by phi1*/
	  calc_Rmatrix(c0toib, phi1, M);
	  rotate(c0toim, M, rotv);
	  c0toim[0]=rotv[0];
	  c0toim[1]=rotv[1];
	  c0toim[2]=rotv[2];
	  bases[0].x[1]=bases[0].xcom+c0toim[0];
	  bases[0].y[1]=bases[0].ycom+c0toim[1];
	  bases[0].z[1]=bases[0].zcom+c0toim[2];
	  
	  /*now use psi0 to define orientation of bases[1].c1toim, already fixed bases[0].c0toim*/
	  if(usephi2==1){
	    /*in this case, only two psis are sampled, rotate phi2 and calc what psi*/
	    psi0=rotate_phi2_mol(phi2, bases, leglen1, leglen0, molfix, molrot, c1toim);
	    cout <<"rotated phi2 instead, calculated psi: "<<psi0<<endl;
	  }else
	    psimol=set_norm_to_psi(psi0, bases, leglen1, leglen0, molfix, molrot, c1toim);//just regular psi
	  
	}
	if(ee+hh+cc==0 && rank==0){
	  cout <<"WRITE crds after orient set "<<endl;
	  write_crds(bases, 0);
	  cout <<"WRITE crds after orient set "<<endl;
	  write_crds(bases, 1);
	}
	/*Now that we have the starting orientations of molecule normals, need to calculate what the Jacobian factor is
	  dPsi/dphi_2 , which is a function of Psi, phi_2, phi_1, cthet1, cthet2, dih. and r?.
	*/
	if(usephi2==1){
	  /*rotation will not change psi.
	    for each psi, we'll do 2*pi over phi values, 
	    and we'll do each of the 2 psi values Npsi/2 times, 
	    so Jac(psi)=2*pi*pi/dpsi, or 2*pi*Npsi/2. integral over J(psi)dpsi
	    will still be 4pi^2 (J for all other psis is zero.)
	  */
	  psi_jac=1;
	  hist_add=psi_jac*delphi1/(1.0*myNrep)*dihbinsize;
	  norm_psi+=psi_jac*delphi1*dihbinsize;//this will end up being the jacobian factor for just integrating over psi.
	  currnormpsi+=psi_jac*delphi1*dihbinsize;
	  cout <<"usephi2=1, both perp:hist add: "<<hist_add<<" norm_psi: "<<norm_psi<<endl;
	}else{
	  psi_jac=numer_calc_Jacobian(bases, psi0, leglen1, leglen0, molrot);//this returns |dphi2/dpsi|
	  hist_add=psi_jac*delphi1/(1.0*myNrep)*dihbinsize;
	  norm_psi+=psi_jac*delphi1*dihbinsize;//this will end up being the jacobian factor for just integrating over psi.
	  currnormpsi+=psi_jac*delphi1*dihbinsize;
	  //cout <<" Jacobian, dphi2/dpsi: "<<psi_jac<<" at psi: "<<psi0<<" and phi1: "<<phi1<<endl;
	  //cout <<"rank: "<<rank<<" X0: "<<x0<<" costheta: "<<costheta1<<" costheta2: "<<costheta2<<endl;
	  
	  
	}
	// cout <<"rank: "<<rank<<" p0 crds: "<<endl;
	// 	      write_crds(bases, 0);
	// 	      cout <<"rank: "<<rank<<" p1 crds: "<<endl;
	// 	      write_crds(bases, 1);
	dx=bases[1].xcom-bases[0].xcom;
	dy=bases[1].ycom-bases[0].ycom;
	dz=bases[1].zcom-bases[0].zcom;
	
	R2=dx*dx+dy*dy+dz*dz;
	R1=sqrt(R2);
	
	//cout <<"rank: "<<rank<<"Initial separation: "<<R1<<endl;
	//cout <<"rank: "<<rank<<"target angles: "<<costheta1<<' '<<costheta2<<' '<<dihed<<" target psi: "<<psi0<<endl;
	calc_three_angle( cthet1,  cthet2,  dih,  bases,  dx, dy, dz,  R1, leglen1, leglen0,  v,  v1, n1, n2);
	//cout <<" calculated angles: "<<cthet1<<' '<<cthet2<<' '<<dih<<endl;
	/*calculate molecule dihedral*/
	psimol=calc_psi_2pi(bases, leglen1, leglen0);
	//cout <<"rank: "<<rank<<"target psi: "<<psi0<<" calculated: "<<psimol<<endl;
	if(abs(psimol-psi0)>tol)cout <<"TARGET PSI : "<<psi0<< " MISMATCHED CALC PSI: "<<psimol<<endl;
	
	for(rep=0;rep<myNrep;rep++){
	  /*start proteins separated by x0 */
	  //cout <<"rank: "<<rank<<"Repeat: "<<rep<<endl;
	  
	  plist.ntotalcomplex=2;
	  Ntotalmol=2;
	  bases[0].z[0]=lz_0;
	  bases[0].x[0]=lx_0;
	  bases[0].y[0]=ly_0;
	  
	  bases[1].z[0]=lz+x0;
	  bases[1].x[0]=lx;
	  bases[1].y[0]=ly;
	  
	  bases[0].xcom=0;
	  bases[0].ycom=0;
	  bases[0].zcom=0;
	  
	  bases[1].xcom=0.0;
	  bases[1].ycom=0.0;
	  bases[1].zcom=x0;
	  
	  bases[0].x[1]=bases[0].xcom+c0toim[0];
	  bases[0].y[1]=bases[0].ycom+c0toim[1];
	  bases[0].z[1]=bases[0].zcom+c0toim[2];
	  
	  bases[0].nfree=1;
	  bases[0].nbnd=0;
	  bases[0].npartner=0;
	  
	  
	  bases[1].x[1]=bases[1].xcom+c1toim[0];
	  bases[1].y[1]=bases[1].ycom+c1toim[1];
	  bases[1].z[1]=bases[1].zcom+c1toim[2];
	  
	  bases[1].nfree=1;
	  bases[1].nbnd=0;
	  bases[1].npartner=0;
	  
	  
	  
	  dx=bases[1].x[0]-bases[0].x[0];
	  dy=bases[1].y[0]-bases[0].y[0];
	  dz=bases[1].z[0]-bases[0].z[0];
	  
	  R2=dx*dx+dy*dy+dz*dz;
	  R1=sqrt(R2);
	  
	  // /*Check the leg and the normal are perp*/
	  r0x=bases[0].x[0]-bases[0].xcom;
	  r0y=bases[0].y[0]-bases[0].ycom;
	  r0z=bases[0].z[0]-bases[0].zcom;
	  r1x=bases[0].x[1]-bases[0].xcom;
	  r1y=bases[0].y[1]-bases[0].ycom;
	  r1z=bases[0].z[1]-bases[0].zcom;
	  if(abs(r0x*r1x+r0y*r1y+r0z*r1z)>1E-6){
	    cout<<"LEG AND NORMAL ARE NOT PERP, MOL 0: Rep: "<<rep<<" rank: "<<rank<<endl;
	    cout <<"rank: "<<rank<<" p0 crds: "<<endl;
	    write_crds(bases, 0);
	    cout <<"c0 vector: "<<c0toim[0]<<' '<<c0toim[1]<<' '<<c0toim[2]<<endl;
	  }
	  //    cout <<"rank: "<<rank<<"p0 dot prod, perp 0? : "<<
	  //    cout <<"rank: "<<rank<<"p0 dir length: "<<sqrt(r1x*r1x+r1y*r1y+r1z*r1z)<<endl;
	  r0x=bases[1].x[0]-bases[1].xcom;
	  r0y=bases[1].y[0]-bases[1].ycom;
	  r0z=bases[1].z[0]-bases[1].zcom;
	  r1x=bases[1].x[1]-bases[1].xcom;
	  r1y=bases[1].y[1]-bases[1].ycom;
	  r1z=bases[1].z[1]-bases[1].zcom;
	  if(abs(r0x*r1x+r0y*r1y+r0z*r1z)>1E-6){
	    cout<<"LEG AND NORMAL ARE NOT PERP, MOL 1: Rep: "<<rep<<" rank: "<<rank<<endl;
	    cout <<"rank: "<<rank<<" p1 crds: "<<endl;
	    write_crds(bases, 1);
	    cout <<"c1 vector: "<<c1toim[0]<<' '<<c1toim[1]<<' '<<c1toim[2]<<endl;
	  }
	  
	  // 	    cout <<"rank: "<<rank<<"p1 dot prod, perp 0? : "<<r0x*r1x+r0y*r1y+r0z*r1z<<endl;
	  
	  /***************************/
	  /*Begin RD simulation*/
	  it=0;
	  currtime=0;
	  
	  //cout <<"rank: "<<rank<<"rep: "<<rep<<" Maxtime; "<<Maxtime<<endl;
	  while(bases[0].nfree>0 &&currtime<Maxtime){
	    i=0;
	    /*Test bimolecular reactions!*/
	    nfree=bases[i].nfree;
	    j=1;
	    
	    dx=bases[1].xcom-bases[0].xcom;
	    dy=bases[1].ycom-bases[0].ycom;
	    dz=bases[1].zcom-bases[0].zcom;
	    
	    R2=dx*dx+dy*dy+dz*dz;
	    R1=sqrt(R2);
	    prob=2;
	    if(R1<Rmax){
	      itmove=1;
	      /*terminate some trajectories.*/
	      deltat=deltat_reac;
	      /*Now tau is scaled by whether the orientation is correct.*/
	      calc_three_angle( cthet1,  cthet2,  dih,  bases,  dx, dy, dz,  R1, leglen1, leglen0,  v,  v1, n1, n2);
	      
	      /*calculate molecule dihedral*/
	      psimol=calc_psi_2pi(bases, leglen1, leglen0);
	      /*here we put no limit on chi, it is limited by excluded volume of particles that
		have theta aligned.
	      */
	      if(1-cthet1<anglelimit){
		
		if(1-cthet2<anglelimit){
		  if(psimol>M_PI)psimol=2*M_PI-psimol;//move from zero to pi
		  if(psimol<psi_limit)
		    prob=exp(-deltat_reac/tau);
		  else
		    prob=2;
		}else
		  prob=2;
	      }else
		prob=2;
	      
	      /*might perform this reaction, depending on k_associate*/
	      rnum=1.0*rand_gsl();
	      p1=1;
	      if(prob<rnum){
		//terminate trajectory
		//cout <<"rank: "<<rank<<"Associate, prob:  "<<rnum<<" rep: "<<rep<<" time: "<<currtime<<endl;
		itmove=0;
		p2=0;
		rxn1=0;
		plist.ntotalcomplex=1;
		bases[p1].nbnd++;
		bases[p2].nbnd++;
		bases[p1].nfree=0;
		bases[p2].nfree=0;
		bases[1].xcom=bases[0].xcom;
		bases[1].ycom=bases[0].ycom;
		bases[1].zcom=bases[0].zcom;
		
		      
	      }else {
		/*propagate particles and avoid overlap*/
		
		dx=sqrt(2.0*deltat*wholep[1].Dx)*GaussV();
		dy=sqrt(2.0*deltat*wholep[1].Dy)*GaussV();
		dz=sqrt(2.0*deltat*wholep[1].Dz)*GaussV();
		
		
		currx=(bases[1].xcom+dx)-(bases[0].xcom);
		curry=(bases[1].ycom+dy)-(bases[0].ycom);
		currz=(bases[1].zcom+dz)-(bases[0].zcom);
		
		R2=currx*currx+curry*curry+currz*currz;
		cit=0;
		while(R2<bindrad*bindrad){
		  cit++;
		  dx=sqrt(2.0*deltat*wholep[1].Dx)*GaussV();
		  dy=sqrt(2.0*deltat*wholep[1].Dy)*GaussV();
		  dz=sqrt(2.0*deltat*wholep[1].Dz)*GaussV();
		  
		  currx=(bases[1].xcom+dx)-(bases[0].xcom);
		  curry=(bases[1].ycom+dy)-(bases[0].ycom);
		  currz=(bases[1].zcom+dz)-(bases[0].zcom);
		  
		  R2=currx*currx+curry*curry+currz*currz;
		  if(cit>1000){
		    cout <<" cannot get outside binding radius! current dt: "<<deltat<<" separation: "<<R2<<" crds, p0: "<<'\t' ;
		    write_crds(bases, 0);
		    cout <<"crds, p1: "<<'\t';
		    write_crds(bases, 1);
		    exit(1);
		  }
		}
		/*update particle positions that do not overlap*/
		tx=sqrt(2.0*deltat*wholep[1].Drx)*GaussV();
		ty=sqrt(2.0*deltat*wholep[1].Dry)*GaussV();
		tz=sqrt(2.0*deltat*wholep[1].Drz)*GaussV();
		rotationEuler( tx,  ty,  tz,MA);
		//vector pre-movement of any particles.
		v[0]=bases[1].x[0]-bases[1].xcom;
		v[1]=bases[1].y[0]-bases[1].ycom;
		v[2]=bases[1].z[0]-bases[1].zcom;//pivot is the center
		rotate(v, MA, v2);//includes the interface that will align
		//v2 currently is for i0, binding site
		bases[1].x[0]=bases[1].xcom+dx+v2[0];
		bases[1].y[0]=bases[1].ycom+dy+v2[1];
		bases[1].z[0]=bases[1].zcom+dz+v2[2];
		
		//rotate direction vector, find orientation before motion
		v[0]=bases[1].x[1]-bases[1].xcom;
		v[1]=bases[1].y[1]-bases[1].ycom;
		v[2]=bases[1].z[1]-bases[1].zcom;//pivot is the center
		
		//now v2 is for i1, the molecule orientation site
		rotate(v, MA, v2);
		bases[1].x[1]=bases[1].xcom+dx+v2[0];
		bases[1].y[1]=bases[1].ycom+dy+v2[1];
		bases[1].z[1]=bases[1].zcom+dz+v2[2];
		
		bases[1].xcom+=dx;
		bases[1].ycom+=dy;
		bases[1].zcom+=dz;
		
		/*now rotate particles 0, who's COM does not change*/
		tx=sqrt(2.0*deltat*wholep[0].Drx)*GaussV();
		ty=sqrt(2.0*deltat*wholep[0].Dry)*GaussV();
		tz=sqrt(2.0*deltat*wholep[0].Drz)*GaussV();
		rotationEuler( tx,  ty,  tz,MB);
		v1[0]=bases[0].x[0]-bases[0].xcom;
		v1[1]=bases[0].y[0]-bases[0].ycom;
		v1[2]=bases[0].z[0]-bases[0].zcom;//pivot is the center
		rotate(v1, MB, v3);//includes the interface that will align
		
		bases[0].x[0]=bases[0].xcom+v3[0];
		bases[0].y[0]=bases[0].ycom+v3[1];
		bases[0].z[0]=bases[0].zcom+v3[2];
		//rotate direction vector
		v1[0]=bases[0].x[1]-bases[0].xcom;
		v1[1]=bases[0].y[1]-bases[0].ycom;
		v1[2]=bases[0].z[1]-bases[0].zcom;//pivot is the center
		rotate(v1, MB, v3);//includes the interface that will align
		bases[0].x[1]=bases[0].xcom+v3[0];
		bases[0].y[1]=bases[0].ycom+v3[1];
		bases[0].z[1]=bases[0].zcom+v3[2];
		
		
	      }
	    }else{
	      /*free diffusion for both particles 
		outside of 'reaction' zone
	      */
	      //max disp=3*sqrt(6*D*dt), (R1-bindrad)/56/D
	      deltat=deltat_reac+scaled*(R1-Rmax)*(R1-Rmax)/(Dtot*2.0);//+deltat_reac;
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
	      
	      bases[1].x[0]=bases[1].xcom+dx+v2[0];
	      bases[1].y[0]=bases[1].ycom+dy+v2[1];
	      bases[1].z[0]=bases[1].zcom+dz+v2[2];
	      
	      //store this vector before moving the com.
	      v[0]=bases[1].x[1]-bases[1].xcom;
	      v[1]=bases[1].y[1]-bases[1].ycom;
	      v[2]=bases[1].z[1]-bases[1].zcom;//pivot is the center
	      rotate(v, M, v2);//includes the interface that will align
	      
	      bases[1].x[1]=bases[1].xcom+dx+v2[0];
	      bases[1].y[1]=bases[1].ycom+dy+v2[1];
	      bases[1].z[1]=bases[1].zcom+dz+v2[2];
	      
	      bases[1].xcom+=dx;
	      bases[1].ycom+=dy;
	      bases[1].zcom+=dz;
	      
	      /*now rotate particles 0, whose COM does not change*/
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
	      
	      v1[0]=bases[0].x[1]-bases[0].xcom;
	      v1[1]=bases[0].y[1]-bases[0].ycom;
	      v1[2]=bases[0].z[1]-bases[0].zcom;//pivot is the center
	      rotate(v1, M, v3);//includes the interface that will align
	      bases[0].x[1]=bases[0].xcom+v3[0];
	      bases[0].y[1]=bases[0].ycom+v3[1];
	      bases[0].z[1]=bases[0].zcom+v3[2];
	      
	    }
	    
	    
	    currtime+=deltat;
	    
	    // 	  if(itmove==1){
	    // 	    phist[it+itmove]+=1;
	    // 	  }else{
	    
	    itmax=itmove+1;
	    if((itmax+it)>Nitbin)
	      itmax=Nitbin-it;
	    
	    for(n=1;n<itmax;n++){
	      if((it+n)%statwrite==0){
		ind=(it+n)/statwrite;
		phist[ind]+=hist_add;
	      }
	    }
	    
	    
	    //cout <<"rank: "<<rank<<"rep: "<<rep<<" prev it: "<<it<<" itmorve: "<<itmove<<" time: "<<currtime<<" nfree: "<<bases[0].nfree<<" final pos: "<<bases[1].xcom<<' '<<bases[1].ycom<<' '<<bases[1].zcom<<endl;
	    it+=itmove;
	    
	  }
	  
	  //if(rep+1%10000==0)
	    // cout <<"rank: "<<rank<<"finished rep: "<<rep<< " at time: "<<currtime<<endl;
	  /*Keep a histogram of final positions for particles that survived until the end*/
	  if(currtime>Maxtime){
	    dx=bases[1].xcom-bases[0].xcom;
	    dy=bases[1].ycom-bases[0].ycom;
	    dz=bases[1].zcom-bases[0].zcom;
	    
	    R2=dx*dx+dy*dy+dz*dz;
	    R1=sqrt(R2);
	    rind=int((R1-bindrad)/delR);
	    //		cout <<"final sep: "<<R1<<endl;
	    calc_three_angle( cthet1,  cthet2,  dih,  bases,  dx, dy, dz,  R1, leglen1, leglen0,  v,  v1, n1, n2);
	    psimol=calc_psi_1pi(bases, leglen1, leglen0);
	    ind_thet=int((cthet1+1)/delthet_out);
	    ind_thet2=int((cthet2+1)/delthet_out);
	    dihbin=int(dih/deldih_out);
	    psiind=int(psimol/delpsi_out);
	    
	    if(rind<Rbins){
	      prhist[(ind_thet*thetbins_out+ind_thet2)*dihbins_out*psibins_out*Rbins+dihbin*psibins_out*Rbins+psiind*Rbins+rind]+=hist_add;
	      p1hist[rind]+=hist_add;
	    }
	    sprob_final+=hist_add;
	    
	  }
	}//end over all reps
	
	//cout <<"rank: "<<rank<<"finished all reps "<<" currhist: "<<phist[0]<<' '<<phist[1]<<" norm psi: "<<norm_psi<<endl;
      }//end phi1 values.
      
      MPI::COMM_WORLD.Reduce(&sprob_final, &sprob_finaltot, 1, MPI::DOUBLE, MPI::SUM, 0);
      if(rank==0){
	jacfile <<costheta1<<'\t'<<costheta2<<'\t'<<psi0<<'\t'<<dihed<<'\t'<<currnormpsi<<'\t'<<norm_psi<<'\t'<<sprob_finaltot/(1.0*nprocs*currnormpsi)<<endl;
      }
      
    }//end dihedral
  }//end both psi
  MPI::COMM_WORLD.Reduce(phist, phist_sum, Nwrite+1, MPI::DOUBLE, MPI::SUM, 0);
  MPI::COMM_WORLD.Reduce(p1hist, p1hist_sum, Rbins, MPI::DOUBLE, MPI::SUM, 0);
  MPI::COMM_WORLD.Reduce(prhist, prhist_sum, thetbins2_out*dihbins_out*psibins_out*Rbins, MPI::DOUBLE, MPI::SUM, 0);
  MPI::COMM_WORLD.Barrier();
  
  /*Now write out final probability histogram*/
  
  /*INSTEAD OF AVERAGING OVER NREP, NOW AVERAGE OVER NPROC  */
  //totpsinorm+=delpsi*norm_psi;
  
  
  if(rank==0){
    //cout <<"curr integral over J(psi): "<<totpsinorm<<" psi: "<<psi0<<" dd: "<<dd<<endl;
    //cout <<"rank: "<<rank<<"open file: "<<tname<<endl;
    /*if usephi2=1 (true) then the S(t) would sum over phist_sum/(1*npro*norm_psi)*dpsi/pi*/
    probfile<<0<<' '<<1<<' '<<1<<endl;
    
    for(i=statwrite;i<Nitbin;i+=statwrite){
      tval=i*deltat_reac;
      ind=i/statwrite;
      //passoc=survive_irr( x0, tval,  Dtot,  bindrad,  alpha,  cof);
      probfile<<tval<<' ' <<phist_sum[ind]/(1.0*nprocs*norm_psi)<<endl;//' '<<1.0-passoc<<endl;
      // stot[ind]+=phist_sum[ind]/(1.0*nprocs*norm_psi[dd])*norm_psi[dd]*delpsi*deltheta*deltheta*dihnorm;//dihnorm is either deldih, which when done Ndih times is pi, or it's pi, done once.
      //cout <<"Time: "<<tval<<" S(t): "<<stot[ind]<<endl;
    }
    
    probfile.close();
    
    //passoc=survive_irr( x0, Maxtime,  Dtot,  bindrad,  alpha,  cof);
    psi_out=psi0;
    if(psiread!=0)psi_out=2.0*M_PI-psi0;
    sonlyfile <<costheta1<<'\t'<<costheta2<<'\t'<<psi_out<<'\t'<<phist_sum[Nwrite-1]/(1.0*nprocs*norm_psi)<<'\t'<<norm_psi<<endl;
    sprobvalue=phist_sum[Nwrite-1]/(1.0*nprocs*norm_psi);
    
    rfile<<0<<'\t'<<0<<'\t'<<0<<'\t';
    
    pavfile<<costheta1<<'\t'<<costheta2<<'\t'<<psi_out<<'\t';
    //totrnorm+=norm_psi*delpsi*deltheta*deltheta*dihnorm;
    for(i=0;i<Rbins;i++){
      R1=bindrad+delR*(i+0.5);
      pavfile<<p1hist_sum[i]/(1.0*nprocs*norm_psi*delR*R1*R1*4.0*M_PI)<<'\t';//norm_psi takes care of loop over phi and chi and both psis 
      rfile<<R1<<'\t';
      //prr[i]+=p1hist_sum[i]/(1.0*nprocs*norm_psi*delR*R1*R1*4.0*M_PI)*norm_psi*delpsi*deltheta*deltheta*dihnorm;
    }
    pavfile<<endl;
    rfile<<endl;
    
    for(k1=0;k1<dihbins_out;k1++){
      for(k2=0;k2<psibins_out;k2++){
	for(j1=0;j1<thetbins_out;j1++){
	  for(j2=0;j2<thetbins_out;j2++){
	    //rfile<<(k1+0.5)*deldih_out<<'\t'<<(k2+0.5)*delpsi_out<<'\t'<<(j1+0.5)*delthet_out-1<<'\t'<<(j2+0.5)*delthet_out-1<<'\t';
	    
	    j=j1*thetbins_out+j2;
	    zeroflag=0;
	    for(i=0;i<zerotest_bins;i++){
	      if(prhist_sum[j*dihbins_out*psibins_out*Rbins+k1*psibins_out*Rbins+k2*Rbins+i]!=0)
		zeroflag=1;
	    }
	    if(zeroflag==1){
	      rfile<<(k2+0.5)*delpsi_out<<'\t'<<(j1+0.5)*delthet_out-1<<'\t'<<(j2+0.5)*delthet_out-1<<'\t';
	      for(i=0;i<Rbins;i++){
		R1=bindrad+delR*(i+0.5);
		rfile<<prhist_sum[j*dihbins_out*psibins_out*Rbins+k1*psibins_out*Rbins+k2*Rbins+i]/(1.0*nprocs*norm_psi*delR*R1*R1*delthet_out*delthet_out*deldih_out*delpsi_out*4.0*M_PI)<<'\t';
	      }
	      rfile<<endl;
	    }
	  }
	}
      }//end psi
    }//end dih chi
    
    rfile.close();
    
    //intparm=sprobvalue;
    dubsend=sprobvalue;
  }//end rank=0;
  MPI::COMM_WORLD.Barrier();
  
  MPI::COMM_WORLD.Bcast(&dubsend, 1, MPI::DOUBLE, 0);
  //sprobvalue=intparm;
  sprobvalue=dubsend;
  if(rank==0)
    cout <<"Norm over all psi angles; "<<norm_psi<<" Dihed: "<<dihed<<" ct1: "<<costheta1<<" ct2: "<<costheta2<<" survival prob: "<<sprobvalue<<endl; 
  
  delete[] MA;
  delete[] M;
  delete[] MB;
  delete[] v3;
  delete[] v2;
  delete[] v1;
  delete[] v;
  delete[] n1;
  delete[] n2;
  
  
  delete[] c0toc1;
  delete[] c1toim;
  delete[] c1toib;
  delete[] c0toim;
  delete[] c0toib;
  delete[] rotv;
  delete[] p0n;
  delete[] p1n;
  
}
void loop_all_psi(int rank, int nprocs, double *psivec, double costheta1, double costheta2, double *norm_psi, double *sprobvec, int Ndihsim, int Nphi, double x0, double lx, double ly, double lz, double leglen1, double leglen0, survParms &plist, double deltat_reac, double scaled, double bindrad, Fullmol *bases, Protein *wholep, double *dihvec, int myNrep, double Rmax, double tau, double anglelimit, double psi_limit, int Nitbin, int Nwrite, double*phist, double *phist_sum, double *p1hist, double *p1hist_sum,double *prhist, double *prhist_sum, ofstream &rfile, ofstream &probfile, ofstream &pavfile, ofstream &jacfile, ofstream &sonlyfile, double tol, int ct1ind, int ct2ind, double *dih_intervals )
{

  int nsplit=8;
  int Npsibins=plist.Npsibins;


  int Npsi_split=(Npsibins*0.5-1)*1.0*nsplit+1;//each interval can be split into eighths, plus the last point.
  int Nthet=plist.Nthet;
  int Nthet_split=(Nthet-1)*1.0*nsplit+1;
  
  int dd;
  int splitflag;
  double psi0;
  int psind0; 
  int sindex0, sindex1, sindex2, sindex3, sindex4, sindex5, sindex6, sindex7, sindex8;
  double psimid, psiquart, change;
  double ang0, ang1, ang2, ang3, ang4, ang5, ang6, ang7, ang8;
  double sprob_tol=plist.sprob_tol;
  for(dd=0;dd<Npsibins;dd+=2){
    psi0=psivec[dd];
    splitflag=1;
    while(splitflag==1){
      /*first do your own value of dd.*/
      psind0=(dd-2)/2*nsplit;//previous interval
      sindex8=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+psind0+nsplit;//this is your current bin
      cout <<" angles: "<<costheta1<< " costheta2; "<<costheta2<<" psi0: "<<psi0<<" sindex:" <<sindex8<<endl;
      run_sim_at_1angleset( rank, nprocs,  psi0,  costheta1,  costheta2, norm_psi[sindex8], sprobvec[sindex8],  Ndihsim,  Nphi, x0,  lx,  ly,  lz,  leglen1, leglen0, plist,  deltat_reac, scaled,  bindrad, bases, wholep, dihvec,  myNrep,  Rmax,  tau, anglelimit, psi_limit,  Nitbin, Nwrite,phist,phist_sum, p1hist,p1hist_sum, prhist, prhist_sum, rfile, probfile, pavfile, jacfile, sonlyfile,  tol, dih_intervals);
      
      if(dd==0)splitflag=0;//break to next value of dd.
      else{
	      /*Now decide whether to split again.*/
	      /*split1*/
	sindex0=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+psind0;//this entry will already be filled in from previous step.
	change=sprobvec[sindex8]-sprobvec[sindex0];
	if(abs(change)>sprob_tol){
	  /*check middle interval psi*/
	  ang8=psivec[dd];
	  ang0=psivec[dd-2];
	  psi0=(ang8-ang0)*0.5+ang0;
	  psimid=psi0;
	  ang4=psimid;
	  sindex4=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+psind0+4;
	  cout <<"add middle point 4! "<<psi0<<" sindex: "<<sindex4<<" change: "<<change<<" upper: "<<sprobvec[sindex8]<<" lower: "<<sprobvec[sindex0]<<endl;
	  run_sim_at_1angleset( rank, nprocs,  psi0,  costheta1,  costheta2, norm_psi[sindex4], sprobvec[sindex4],  Ndihsim,  Nphi, x0,  lx,  ly,  lz,  leglen1, leglen0, plist,  deltat_reac, scaled,  bindrad, bases, wholep, dihvec,  myNrep,  Rmax,  tau, anglelimit, psi_limit,  Nitbin, Nwrite,phist,phist_sum, p1hist,p1hist_sum, prhist, prhist_sum, rfile, probfile, pavfile, jacfile, sonlyfile,  tol, dih_intervals);
	  
	  /*split2, left first*/
	  
	  change=sprobvec[sindex4]-sprobvec[sindex0];
	  if(abs(change)>sprob_tol){
	    psi0=(ang4-ang0)*0.5+ang0;
	    psiquart=psi0;
	    ang2=psiquart;
	    sindex2=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+psind0+2;
	    cout <<"add middle point 2! "<<psi0<<" sindex: "<<sindex2<<" change: "<<change<<" upper: "<<sprobvec[sindex4]<<" lower: "<<sprobvec[sindex0]<<endl;
	    run_sim_at_1angleset( rank, nprocs,  psi0,  costheta1,  costheta2, norm_psi[sindex2], sprobvec[sindex2],  Ndihsim,  Nphi, x0,  lx,  ly,  lz,  leglen1, leglen0, plist,  deltat_reac, scaled,  bindrad, bases, wholep, dihvec,  myNrep,  Rmax,  tau, anglelimit, psi_limit,  Nitbin, Nwrite,phist,phist_sum, p1hist,p1hist_sum, prhist, prhist_sum, rfile, probfile, pavfile, jacfile, sonlyfile,  tol, dih_intervals);
	    
	    /*Now split this interval, left first*/
	    change=sprobvec[sindex2]-sprobvec[sindex0];
	    if(abs(change)>sprob_tol){
	      psi0=(ang2-ang0)*0.5+ang0;
	      ang1=psi0;
	      sindex1=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+psind0+1;
	      cout <<"add middle point 1! "<<psi0<<" sindex: "<<sindex1<<" change: "<<change<<" upper: "<<sprobvec[sindex2]<<" lower: "<<sprobvec[sindex0]<<endl;
	      run_sim_at_1angleset( rank, nprocs,  psi0,  costheta1,  costheta2, norm_psi[sindex1], sprobvec[sindex1],  Ndihsim,  Nphi, x0,  lx,  ly,  lz,  leglen1, leglen0, plist,  deltat_reac, scaled,  bindrad, bases, wholep, dihvec,  myNrep,  Rmax,  tau, anglelimit, psi_limit,  Nitbin, Nwrite,phist,phist_sum, p1hist,p1hist_sum, prhist, prhist_sum, rfile, probfile, pavfile, jacfile, sonlyfile,  tol, dih_intervals);
	      
	    }
	    
	    /*Now split this interval, right*/
	    change=sprobvec[sindex4]-sprobvec[sindex2];
	    if(abs(change)>sprob_tol){
	      psi0=(ang4-ang2)*0.5+ang2;
	      ang3=psi0;
	      sindex3=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+psind0+3;
	      cout <<"add middle point 3! "<<psi0<<" sindex: "<<sindex3<<" change: "<<change<<" upper: "<<sprobvec[sindex4]<<" lower: "<<sprobvec[sindex2]<<endl;
	      run_sim_at_1angleset( rank, nprocs,  psi0,  costheta1,  costheta2, norm_psi[sindex3], sprobvec[sindex3],  Ndihsim,  Nphi, x0,  lx,  ly,  lz,  leglen1, leglen0, plist,  deltat_reac, scaled,  bindrad, bases, wholep, dihvec,  myNrep,  Rmax,  tau, anglelimit, psi_limit,  Nitbin, Nwrite,phist,phist_sum, p1hist,p1hist_sum, prhist, prhist_sum, rfile, probfile, pavfile, jacfile, sonlyfile,  tol, dih_intervals);
	      
	    }
	    
	  }//end inside 4-0 interval
	  /*split2, right now*/
	  change=sprobvec[sindex8]-sprobvec[sindex4];
	  if(abs(change)>sprob_tol){
	    psi0=(ang8-ang4)*0.5+ang4;//ang4=psimid
	    ang6=psi0;
	    psiquart=psi0;
	    sindex6=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+psind0+6;
	    cout <<"add middle point 6! "<<psi0<<" sindex: "<<sindex6<<" change: "<<change<<" upper: "<<sprobvec[sindex8]<<" lower: "<<sprobvec[sindex4]<<endl;
	    run_sim_at_1angleset( rank, nprocs,  psi0,  costheta1,  costheta2, norm_psi[sindex6], sprobvec[sindex6],  Ndihsim,  Nphi, x0,  lx,  ly,  lz,  leglen1, leglen0, plist,  deltat_reac, scaled,  bindrad, bases, wholep, dihvec,  myNrep,  Rmax,  tau, anglelimit, psi_limit,  Nitbin, Nwrite,phist,phist_sum, p1hist,p1hist_sum, prhist, prhist_sum, rfile, probfile, pavfile, jacfile, sonlyfile,  tol, dih_intervals);
	    
	    /*Now split this interval, left first*/
	    change=sprobvec[sindex6]-sprobvec[sindex4];
	    if(abs(change)>sprob_tol){
	      psi0=(ang6-ang4)*0.5+ang4;
	      ang5=psi0;
	      sindex5=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+psind0+5;
	      cout <<"add middle point 5! "<<psi0<<" sindex: "<<sindex5<<" change: "<<change<<" upper: "<<sprobvec[sindex6]<<" lower: "<<sprobvec[sindex4]<<endl;
	      run_sim_at_1angleset( rank, nprocs,  psi0,  costheta1,  costheta2, norm_psi[sindex5], sprobvec[sindex5],  Ndihsim,  Nphi, x0,  lx,  ly,  lz,  leglen1, leglen0, plist,  deltat_reac, scaled,  bindrad, bases, wholep, dihvec,  myNrep,  Rmax,  tau, anglelimit, psi_limit,  Nitbin, Nwrite,phist,phist_sum, p1hist,p1hist_sum, prhist, prhist_sum, rfile, probfile, pavfile, jacfile, sonlyfile,  tol, dih_intervals);
	      
	    }
	    
	    /*Now split this interval, right*/
	    change=sprobvec[sindex8]-sprobvec[sindex6];
	    if(abs(change)>sprob_tol){
	      psi0=(ang8-ang6)*0.5+ang6;
	      ang7=psi0;
	      sindex7=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+psind0+7;
	      cout <<"add middle point 7! "<<psi0<<" sindex: "<<sindex7<<" change: "<<change<<" upper: "<<sprobvec[sindex8]<<" lower: "<<sprobvec[sindex6]<<endl;
	      run_sim_at_1angleset( rank, nprocs,  psi0,  costheta1,  costheta2, norm_psi[sindex7], sprobvec[sindex7],  Ndihsim,  Nphi, x0,  lx,  ly,  lz,  leglen1, leglen0, plist,  deltat_reac, scaled,  bindrad, bases, wholep, dihvec,  myNrep,  Rmax,  tau, anglelimit, psi_limit,  Nitbin, Nwrite,phist,phist_sum, p1hist,p1hist_sum, prhist, prhist_sum, rfile, probfile, pavfile, jacfile, sonlyfile,  tol, dih_intervals);
	    }
	    
	  }
	  splitflag=0;//finished checking all internal intervals
	}else
	  splitflag=0;//first interval is already close enough
      }
    }
    
  }//end psi angles	
  
}
void loop_over_costheta2_gridpt(int rank, int nprocs, double *psivec, double costheta1, double *costheta, double *norm_psi, double *sprobvec, int Ndihsim, int Nphi, double x0, double lx, double ly, double lz, double leglen1, double leglen0, survParms &plist, double deltat_reac, double scaled, double bindrad, Fullmol *bases, Protein *wholep, double *dihvec, int myNrep, double Rmax, double tau, double anglelimit, double psi_limit, int Nitbin, int Nwrite, double*phist, double *phist_sum, double *p1hist, double *p1hist_sum,double *prhist, double *prhist_sum, ofstream &rfile, ofstream &probfile, ofstream &pavfile, ofstream &jacfile, ofstream &sonlyfile, double tol, int ct1ind, int ctstart, double *dih_intervals )
{
  int bt;
  int ct2ind;
  int splitflagct2;
  int nsplit=8;
  int Npsibins=plist.Npsibins;

  int Npsi_split=(Npsibins*0.5-1)*1.0*nsplit+1;//each interval can be split into eighths, plus the last point.
  int Nthet=plist.Nthet;
  int Nthet_split=(Nthet-1)*1.0*nsplit+1;
  
  int ctindex2[8];
  double change;
  //double sprob_tmp=0.1;
  double sprob_tol=plist.sprob_tol;
  double cost4, cost2, cost6, cost7, cost5, cost3, cost1, cost8, cost0; 
  double costheta2;
  for(bt=ctstart;bt<Nthet;bt++){
    /*for each set of cos(thetas), loop over all dihedrals.
     */
    splitflagct2=1;
    ct2ind=bt*8;//this is your current index.
    ctindex2[8]=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+0;//for first psi.
    
    costheta2=costheta[bt];
    cout <<" bt: "<<bt<<" costheta2: "<<costheta2<<" ctindex: "<<ctindex2[8]<<" ct1ind: "<<ct1ind<<" ct2ind: "<<ct2ind<<endl;
    cout <<"starting pos: "<<lx<<'\t'<<ly<<'\t'<<lz<<endl;
    while(splitflagct2==1){
      loop_all_psi( rank,  nprocs,  psivec, costheta1,  costheta2, norm_psi,sprobvec,  Ndihsim,  Nphi, x0, lx, ly,  lz,  leglen1, leglen0, plist,  deltat_reac,  scaled,  bindrad, bases, wholep, dihvec, myNrep,  Rmax,  tau,  anglelimit,  psi_limit,  Nitbin, Nwrite, phist, phist_sum, p1hist,p1hist_sum,prhist,prhist_sum, rfile,probfile, pavfile, jacfile, sonlyfile, tol, ct1ind, ct2ind, dih_intervals);
      if(bt==ctstart)  splitflagct2=0;
      else{
	ct2ind=(bt-1)*8;//this is your previous index.
	ctindex2[0]=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+0;//for first psi.
	change=sprobvec[ctindex2[8]]-sprobvec[ctindex2[0]];
	/*consider splitting the interval*/
	if(abs(change)>sprob_tol){
	  ct2ind=(bt-1)*8+4;
	  ctindex2[4]=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+0;//for first psi.
	  cost8=costheta[bt];
	  cost0=costheta[bt-1];//this is the highestnumber though 
	  costheta2=(cost0-cost8)*0.5+cost8;
	  cost4=costheta2;
	  loop_all_psi( rank,  nprocs,  psivec, costheta1,  costheta2, norm_psi,sprobvec,  Ndihsim,  Nphi, x0, lx, ly,  lz,  leglen1, leglen0, plist,  deltat_reac,  scaled,  bindrad, bases, wholep, dihvec, myNrep,  Rmax,  tau,  anglelimit,  psi_limit,  Nitbin, Nwrite, phist, phist_sum, p1hist,p1hist_sum,prhist,prhist_sum, rfile,probfile, pavfile, jacfile, sonlyfile, tol, ct1ind, ct2ind, dih_intervals);
	  /*now can split this interval, left side first*/
	  change=sprobvec[ctindex2[4]]-sprobvec[ctindex2[0]];
	  if(abs(change)>sprob_tol){
	    ct2ind=(bt-1)*8+2;
	    ctindex2[2]=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+0;//for first psi.
	    costheta2=(cost0-cost4)*0.5+cost4;
	    cost2=costheta2;
	    loop_all_psi( rank,  nprocs,  psivec, costheta1,  costheta2, norm_psi,sprobvec,  Ndihsim,  Nphi, x0, lx, ly,  lz,  leglen1, leglen0, plist,  deltat_reac,  scaled,  bindrad, bases, wholep, dihvec, myNrep,  Rmax,  tau,  anglelimit,  psi_limit,  Nitbin, Nwrite, phist, phist_sum, p1hist,p1hist_sum,prhist,prhist_sum, rfile,probfile, pavfile, jacfile, sonlyfile, tol, ct1ind, ct2ind, dih_intervals);
	    /*now split this, left first*/
	    change=sprobvec[ctindex2[2]]-sprobvec[ctindex2[0]];
	    if(abs(change)>sprob_tol){
	      ct2ind=(bt-1)*8+1;
	      ctindex2[1]=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+0;//for first psi.
	      costheta2=(cost0-cost2)*0.5+cost2;
	      cost1=costheta2;
	      loop_all_psi( rank,  nprocs,  psivec, costheta1,  costheta2, norm_psi,sprobvec,  Ndihsim,  Nphi, x0, lx, ly,  lz,  leglen1, leglen0, plist,  deltat_reac,  scaled,  bindrad, bases, wholep, dihvec, myNrep,  Rmax,  tau,  anglelimit,  psi_limit,  Nitbin, Nwrite, phist, phist_sum, p1hist,p1hist_sum,prhist,prhist_sum, rfile,probfile, pavfile, jacfile, sonlyfile, tol, ct1ind, ct2ind, dih_intervals);
	    }
	    /*now split this, right*/
	    change=sprobvec[ctindex2[4]]-sprobvec[ctindex2[2]];
	    if(abs(change)>sprob_tol){
	      ct2ind=(bt-1)*8+3;
	      ctindex2[3]=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+0;//for first psi.
	      costheta2=(cost2-cost4)*0.5+cost4;
	      cost3=costheta2;
	      loop_all_psi( rank,  nprocs,  psivec, costheta1,  costheta2, norm_psi,sprobvec,  Ndihsim,  Nphi, x0, lx, ly,  lz,  leglen1, leglen0, plist,  deltat_reac,  scaled,  bindrad, bases, wholep, dihvec, myNrep,  Rmax,  tau,  anglelimit,  psi_limit,  Nitbin, Nwrite, phist, phist_sum, p1hist,p1hist_sum,prhist,prhist_sum, rfile,probfile, pavfile, jacfile, sonlyfile, tol, ct1ind, ct2ind, dih_intervals);
	    }
	    
	  }
	  /*now can split this interval, right side*/
	  change=sprobvec[ctindex2[8]]-sprobvec[ctindex2[4]];
	  if(abs(change)>sprob_tol){
	    ct2ind=(bt-1)*8+6;
	    ctindex2[6]=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+0;//for first psi.
	    costheta2=(cost4-cost8)*0.5+cost8;
	    cost6=costheta2;
	    loop_all_psi( rank,  nprocs,  psivec, costheta1,  costheta2, norm_psi,sprobvec,  Ndihsim,  Nphi, x0, lx, ly,  lz,  leglen1, leglen0, plist,  deltat_reac,  scaled,  bindrad, bases, wholep, dihvec, myNrep,  Rmax,  tau,  anglelimit,  psi_limit,  Nitbin, Nwrite, phist, phist_sum, p1hist,p1hist_sum,prhist,prhist_sum, rfile,probfile, pavfile, jacfile, sonlyfile, tol, ct1ind, ct2ind, dih_intervals);
	    /*now split this, left first*/
	    change=sprobvec[ctindex2[6]]-sprobvec[ctindex2[4]];
	    if(abs(change)>sprob_tol){
	      ct2ind=(bt-1)*8+5;
	      ctindex2[5]=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+0;//for first psi.
	      costheta2=(cost4-cost6)*0.5+cost6;
	      cost5=costheta2;
	      loop_all_psi( rank,  nprocs,  psivec, costheta1,  costheta2, norm_psi,sprobvec,  Ndihsim,  Nphi, x0, lx, ly,  lz,  leglen1, leglen0, plist,  deltat_reac,  scaled,  bindrad, bases, wholep, dihvec, myNrep,  Rmax,  tau,  anglelimit,  psi_limit,  Nitbin, Nwrite, phist, phist_sum, p1hist,p1hist_sum,prhist,prhist_sum, rfile,probfile, pavfile, jacfile, sonlyfile, tol, ct1ind, ct2ind, dih_intervals);
	    }
	    /*now split this, right*/
	    change=sprobvec[ctindex2[8]]-sprobvec[ctindex2[6]];
	    if(abs(change)>sprob_tol){
	      ct2ind=(bt-1)*8+7;
	      ctindex2[7]=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+0;//for first psi.
	      costheta2=(cost6-cost8)*0.5+cost8;
	      cost7=costheta2;
	      loop_all_psi( rank,  nprocs,  psivec, costheta1,  costheta2, norm_psi,sprobvec,  Ndihsim,  Nphi, x0, lx, ly,  lz,  leglen1, leglen0, plist,  deltat_reac,  scaled,  bindrad, bases, wholep, dihvec, myNrep,  Rmax,  tau,  anglelimit,  psi_limit,  Nitbin, Nwrite, phist, phist_sum, p1hist,p1hist_sum,prhist,prhist_sum, rfile,probfile, pavfile, jacfile, sonlyfile, tol, ct1ind, ct2ind, dih_intervals);
	    }
	    
	  }
	  
	}
	splitflagct2=0;//finished all splits.
      }
      
      
    }//end while loop
  }//end second angle
  
}

void loop_over_costheta2_halfgridpt(int rank, int nprocs, double *psivec, double costheta1, double costheta_next, double *norm_psi, double *sprobvec, int Ndihsim, int Nphi, double x0, double lx, double ly, double lz, double leglen1, double leglen0, survParms &plist, double deltat_reac, double scaled, double bindrad, Fullmol *bases, Protein *wholep, double *dihvec, int myNrep, double Rmax, double tau, double anglelimit, double psi_limit, int Nitbin, int Nwrite, double*phist, double *phist_sum, double *p1hist, double *p1hist_sum,double *prhist, double *prhist_sum, ofstream &rfile, ofstream &probfile, ofstream &pavfile, ofstream &jacfile, ofstream &sonlyfile, double tol, int ct1ind, int ct1start, double *dih_intervals)
{
  int bt;
  int ct2ind;
  int splitflagct2;
  
  int nsplit=8;
  int Npsibins=plist.Npsibins;

  int Npsi_split=(Npsibins*0.5-1)*1.0*nsplit+1;//each interval can be split into eighths, plus the last point.
  int Nthet=plist.Nthet;
  int Nthet_split=(Nthet-1)*1.0*nsplit+1;
  double deltheta=2.0/(1.0*(plist.Nthet-1));
  int ctindex2[8];
  double change;
  //double sprob_tmp=0.01;
  double sprob_tol=plist.sprob_tol;
  double cost4, cost2, cost6, cost7, cost5, cost3, cost1, cost8, cost0; 
  double costheta2;
  
    /*for each set of cos(thetas), loop over all dihedrals.
     */
  int i;
  int i1=0;
  if(ct1start!=ct1ind)i1=1;
  
  ctindex2[0]=ct1ind*Nthet_split*Npsi_split+ct1start*Npsi_split+0;//for first psi.
  for(i=i1;i<2;i++){
    //only 2 bins
    splitflagct2=1;
    ct2ind=ct1start+i*4;
    
    if(i==0)
      costheta2=costheta1;
    else
      costheta2=costheta_next;
    cout <<" half grid:  costheta2: "<<costheta2<<" ctindex: "<<ctindex2[8]<<" ct1ind: "<<ct1ind<<" ct2ind: "<<ct2ind<<endl;
    while(splitflagct2==1){
      loop_all_psi( rank,  nprocs,  psivec, costheta1,  costheta2, norm_psi,sprobvec,  Ndihsim,  Nphi, x0, lx, ly,  lz,  leglen1, leglen0, plist,  deltat_reac,  scaled,  bindrad, bases, wholep, dihvec, myNrep,  Rmax,  tau,  anglelimit,  psi_limit,  Nitbin, Nwrite, phist, phist_sum, p1hist,p1hist_sum,prhist,prhist_sum, rfile,probfile, pavfile, jacfile, sonlyfile, tol, ct1ind, ct2ind,dih_intervals);
      if(i==0)  splitflagct2=0;
      else{
	ctindex2[8]=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+0;//for first psi.
	change=sprobvec[ctindex2[8]]-sprobvec[ctindex2[0]];
	/*consider splitting the interval, but fewer times, left side first.*/

	if(abs(change)>sprob_tol){
	  ct2ind=ct1start+2;
	  ctindex2[4]=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+0;//for first psi.
	  cost8=costheta_next;
	  cost0=costheta_next+deltheta/2.0;//this is the highestnumber though 
	  costheta2=(cost0-cost8)*0.5+cost8;
	  cost4=costheta2;
	  cout <<"add point in costheta2,half : "<<cost4<<" change: "<<change<<" ct1ind: "<<ct1ind<<" ct2ind: "<<ct2ind<<endl;
	  loop_all_psi( rank,  nprocs,  psivec, costheta1,  costheta2, norm_psi,sprobvec,  Ndihsim,  Nphi, x0, lx, ly,  lz,  leglen1, leglen0, plist,  deltat_reac,  scaled,  bindrad, bases, wholep, dihvec, myNrep,  Rmax,  tau,  anglelimit,  psi_limit,  Nitbin, Nwrite, phist, phist_sum, p1hist,p1hist_sum,prhist,prhist_sum, rfile,probfile, pavfile, jacfile, sonlyfile, tol, ct1ind, ct2ind,dih_intervals);
	  /*now can split this interval, left side first*/
	  change=sprobvec[ctindex2[4]]-sprobvec[ctindex2[0]];
	  if(abs(change)>sprob_tol){
	    ct2ind=ct1start+1;
	    ctindex2[2]=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+0;//for first psi.
	    costheta2=(cost0-cost4)*0.5+cost4;
	    cost2=costheta2;
	    cout <<"add point in costheta2,half : "<<cost2<<" change: "<<change<<" ct1ind: "<<ct1ind<<" ct2ind: "<<ct2ind<<endl;
	    loop_all_psi( rank,  nprocs,  psivec, costheta1,  costheta2, norm_psi,sprobvec,  Ndihsim,  Nphi, x0, lx, ly,  lz,  leglen1, leglen0, plist,  deltat_reac,  scaled,  bindrad, bases, wholep, dihvec, myNrep,  Rmax,  tau,  anglelimit,  psi_limit,  Nitbin, Nwrite, phist, phist_sum, p1hist,p1hist_sum,prhist,prhist_sum, rfile,probfile, pavfile, jacfile, sonlyfile, tol, ct1ind, ct2ind,dih_intervals);
	    
	  }
	  /*now can split this interval, right side*/
	  change=sprobvec[ctindex2[8]]-sprobvec[ctindex2[4]];
	  if(abs(change)>sprob_tol){
	    ct2ind=ct1start+3;
	    ctindex2[6]=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+0;//for first psi.
	    costheta2=(cost4-cost8)*0.5+cost8;
	    cost6=costheta2;
	    cout <<"add point in costheta2,half : "<<cost6<<" change: "<<change<<" ct1ind: "<<ct1ind<<" ct2ind: "<<ct2ind<<endl;
	    loop_all_psi( rank,  nprocs,  psivec, costheta1,  costheta2, norm_psi,sprobvec,  Ndihsim,  Nphi, x0, lx, ly,  lz,  leglen1, leglen0, plist,  deltat_reac,  scaled,  bindrad, bases, wholep, dihvec, myNrep,  Rmax,  tau,  anglelimit,  psi_limit,  Nitbin, Nwrite, phist, phist_sum, p1hist,p1hist_sum,prhist,prhist_sum, rfile,probfile, pavfile, jacfile, sonlyfile, tol, ct1ind, ct2ind,dih_intervals);
	    
	  }
	  
	}
	splitflagct2=0;//finished all splits.
      }
      
      
    }//end while loop
  }//end second angle
  
}

void loop_over_costheta2_quartergridpt(int rank, int nprocs, double *psivec, double costheta1, double costheta_next, double *norm_psi, double *sprobvec, int Ndihsim, int Nphi, double x0, double lx, double ly, double lz, double leglen1, double leglen0, survParms &plist, double deltat_reac, double scaled, double bindrad, Fullmol *bases, Protein *wholep, double *dihvec, int myNrep, double Rmax, double tau, double anglelimit, double psi_limit, int Nitbin, int Nwrite, double*phist, double *phist_sum, double *p1hist, double *p1hist_sum,double *prhist, double *prhist_sum, ofstream &rfile, ofstream &probfile, ofstream &pavfile, ofstream &jacfile, ofstream &sonlyfile, double tol, int ct1ind, int ct1start, double *dih_intervals)
{
  int bt;
  int ct2ind;
  int splitflagct2;
  
  int nsplit=8;
  int Npsibins=plist.Npsibins;

  int Npsi_split=(Npsibins*0.5-1)*1.0*nsplit+1;//each interval can be split into eighths, plus the last point.
  int Nthet=plist.Nthet;
  int Nthet_split=(Nthet-1)*1.0*nsplit+1;
  double deltheta=2.0/(1.0*(plist.Nthet-1));
  
  int ctindex2[8];
  double change;
  //double sprob_tmp=0.01;
  double sprob_tol=plist.sprob_tol;
  double cost4, cost2, cost6, cost7, cost5, cost3, cost1, cost8, cost0; 
  double costheta2;
  
    /*for each set of cos(thetas), loop over all dihedrals.
     */
  int i;
  int i1=0;
  if(ct1ind!=ct1start)i1=1;
  ctindex2[0]=ct1ind*Nthet_split*Npsi_split+ct1start*Npsi_split+0;//for first psi.
  for(i=i1;i<2;i++){
    //only 2 bins
    splitflagct2=1;
    ct2ind=ct1start+i*2;

    if(i==0)
      costheta2=costheta1;
    else
      costheta2=costheta_next;
    cout <<" quartergrid costheta2: "<<costheta2<<" ctindex: "<<ctindex2[8]<<" ct1ind: "<<ct1ind<<" ct2ind: "<<ct2ind<<endl;
    while(splitflagct2==1){
      loop_all_psi( rank,  nprocs,  psivec, costheta1,  costheta2, norm_psi,sprobvec,  Ndihsim,  Nphi, x0, lx, ly,  lz,  leglen1, leglen0, plist,  deltat_reac,  scaled,  bindrad, bases, wholep, dihvec, myNrep,  Rmax,  tau,  anglelimit,  psi_limit,  Nitbin, Nwrite, phist, phist_sum, p1hist,p1hist_sum,prhist,prhist_sum, rfile,probfile, pavfile, jacfile, sonlyfile, tol, ct1ind, ct2ind,dih_intervals);
      if(i==0)  splitflagct2=0;
      else{
	ctindex2[8]=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+0;//for first psi.
	change=sprobvec[ctindex2[8]]-sprobvec[ctindex2[0]];
	/*consider splitting the interval, but fewer times, left side first.*/

	if(abs(change)>sprob_tol){
	  
	  ct2ind=ct1start+1;
	  ctindex2[4]=ct1ind*Nthet_split*Npsi_split+ct2ind*Npsi_split+0;//for first psi.
	  cost8=costheta_next;
	  cost0=costheta_next+0.25*deltheta;//this is the highestnumber though 
	  costheta2=(cost0-cost8)*0.5+cost8;
	  cost4=costheta2;
	  cout <<"add point in costheta2, : "<<cost4<<" change: "<<change<<" ct1ind: "<<ct1ind<<" ct2ind: "<<ct2ind<<endl;
	  loop_all_psi( rank,  nprocs,  psivec, costheta1,  costheta2, norm_psi,sprobvec,  Ndihsim,  Nphi, x0, lx, ly,  lz,  leglen1, leglen0, plist,  deltat_reac,  scaled,  bindrad, bases, wholep, dihvec, myNrep,  Rmax,  tau,  anglelimit,  psi_limit,  Nitbin, Nwrite, phist, phist_sum, p1hist,p1hist_sum,prhist,prhist_sum, rfile,probfile, pavfile, jacfile, sonlyfile, tol, ct1ind, ct2ind,dih_intervals);
	  	    
	  }
	  
	
	splitflagct2=0;//finished all splits.
      }
      
      
    }//end while loop
  }//end second angle
  
}
void loop_over_costheta2_eighthgridpt(int rank, int nprocs, double *psivec, double costheta1, double costheta_next, double *norm_psi, double *sprobvec, int Ndihsim, int Nphi, double x0, double lx, double ly, double lz, double leglen1, double leglen0, survParms &plist, double deltat_reac, double scaled, double bindrad, Fullmol *bases, Protein *wholep, double *dihvec, int myNrep, double Rmax, double tau, double anglelimit, double psi_limit, int Nitbin, int Nwrite, double*phist, double *phist_sum, double *p1hist, double *p1hist_sum,double *prhist, double *prhist_sum, ofstream &rfile, ofstream &probfile, ofstream &pavfile, ofstream &jacfile, ofstream &sonlyfile, double tol, int ct1ind, double *dih_intervals)
{
  int bt;
  int ct2ind;
  int splitflagct2;
  
  int nsplit=8;
  int Npsibins=plist.Npsibins;

  int Npsi_split=(Npsibins*0.5-1)*1.0*nsplit+1;//each interval can be split into eighths, plus the last point.
  int Nthet=plist.Nthet;
  int Nthet_split=(Nthet-1)*1.0*nsplit+1;
  
  int ctindex2[8];
  double change;
  //double sprob_tmp=0.01;
  double sprob_tol=plist.sprob_tol;
  double cost4, cost2, cost6, cost7, cost5, cost3, cost1, cost8, cost0; 
  double costheta2;
  
    /*for each set of cos(thetas), loop over all dihedrals.
     */
  int i;
  ctindex2[0]=ct1ind*Nthet_split*Npsi_split+ct1ind*Npsi_split+0;//for first psi.
  for(i=0;i<2;i++){
    //only 2 bins
    splitflagct2=1;
    ct2ind=ct1ind+i*1;

    if(i==0)
      costheta2=costheta1;
    else
      costheta2=costheta_next;
    cout <<" eighthgrid costheta2: "<<costheta2<<" ct1ind: "<<ct1ind<<" ct2ind: "<<ct2ind<<endl;
    
    loop_all_psi( rank,  nprocs,  psivec, costheta1,  costheta2, norm_psi,sprobvec,  Ndihsim,  Nphi, x0, lx, ly,  lz,  leglen1, leglen0, plist,  deltat_reac,  scaled,  bindrad, bases, wholep, dihvec, myNrep,  Rmax,  tau,  anglelimit,  psi_limit,  Nitbin, Nwrite, phist, phist_sum, p1hist,p1hist_sum,prhist,prhist_sum, rfile,probfile, pavfile, jacfile, sonlyfile, tol, ct1ind, ct2ind,dih_intervals);
            
    
  }//end second angle
  
}

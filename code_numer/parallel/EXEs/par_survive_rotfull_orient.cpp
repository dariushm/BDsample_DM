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
  double dihed0;
  int Nx0;
  int Nthet;
  int thetbins;
  int dihbins;
  int Npsibins;
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
void write_crds(Fullmol *bases, int p1);
double pirr_pfree_ratio_ps(double rcurr, double r0, double tcurr, double Dtot, double bindrad, double alpha, double ps_prev, double rtol);
double survive_irr(double r0, double tcurr, double Dtot, double bindrad, double alpha, double cof);
double pirrev_value(double rcurr, double r0, double tcurr, double Dtot, double bindrad, double alpha);
double pfree_value_norm(double rcurr, double r0, double tcurr, double Dtot, double bindrad,double alpha);
double calc_psi(Fullmol *bases, double *v, double *v1, double *n1, double *n2);
void calc_three_angle(double &cthet1, double &cthet2, double &dih, Fullmol *bases, double dx, double dy, double dz, double R1, double leglen, double *v, double *v1, double *n1, double *n2);


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
  Parms plist;
  int arrsize=12;
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
    intparm[6]=plist.Nx0;
    intparm[7]=plist.Nthet;
    intparm[8]=plist.thetbins;
    intparm[9]=plist.dihbins;
    intparm[10]=plist.Npsibins;
    dubparm[0]=plist.Maxtime;
    dubparm[1]=plist.rate;
    dubparm[2]=plist.eps_scale;
    dubparm[3]=plist.dt_scale;
    dubparm[4]=plist.D;
    dubparm[5]=plist.leglen;
    dubparm[6]=plist.x0;
    dubparm[7]=plist.dihed0;
    
    
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
  plist.Nx0=intparm[6];
  plist.Nthet=intparm[7];
  plist.thetbins=intparm[8];
  plist.dihbins=intparm[9];
  plist.Npsibins=intparm[9];
  plist.Maxtime=dubparm[0];
  plist.rate=dubparm[1];
  plist.eps_scale=dubparm[2];
  plist.dt_scale=dubparm[3];
  plist.D=dubparm[4];
  plist.leglen=dubparm[5];
  plist.x0=dubparm[6];
  plist.dihed0=dubparm[7];

  
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
  
  wholep[0].ninterface=2;
  wholep[1].ninterface=2;
  bases[0].ninterface=2;
  bases[1].ninterface=2;
  bindrad[0]=1;
  kr[0]=plist.rate;
  cout <<"rank: "<<rank<<"ACTIVATION RATE: "<<kr[0]<<"  radius: "<<bindrad[0]<<endl;
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
  double box_x=plist.xboxl;
  double box_y=plist.yboxl;
  double box_z=plist.zboxl;//nm
  double xtot, ytot, ztot;
  int nfree, wprot;
  int p, i1, i2;
  int np;
  double r2, r;
  double maxsep2=bindrad[0]*bindrad[0];//plist.maxsep2;
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

  double kpi=4.0*M_PI*bindrad[0];
  double leglen=plist.leglen;
  double leglen2=leglen*leglen;
  cout <<"rank: "<<rank<<"Leg length: "<<leglen<<endl;
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
  
  r0=bindrad[0];
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
  cout <<"rank: "<<rank<<"Initial separation: "<<plist.x0<<endl;
  x0vec[0]=plist.x0;//for k<inf, this will not go to 1 
  for(i=1;i<Nx0;i++){
    x0vec[i]=plist.x0+0.05*i;
    cout <<"rank: "<<rank<<"X0, i: "<<i<<" x0: "<<x0vec[i]<<endl;
  }
  
  double kappa=kact/(4.0*M_PI*bindrad[0]*bindrad[0]);
  cout <<"rank: "<<rank<<"Kact: "<<kact<<" kappa: "<<kappa<<endl;
  double epsilon=plist.eps_scale*Dtot/kappa;
  double tau=epsilon/kappa;
  cout <<"rank: "<<rank<<"epsscale: "<<plist.eps_scale<<" epsilon: "<<epsilon<<" tau: "<<tau<<endl;
  Rmax=bindrad[0]+epsilon;
  cout <<"rank: "<<rank<<"Bindrad: "<<bindrad[0]<<" reaction limit: "<<Rmax<<endl;
  //  cout <<"rank: "<<rank<<"different time step expansions: "<<epsilon*epsilon/Dtot/100.0<<" other: "<<Rmax*Rmax*0.0001/(2.0*Dtot)<<endl;
  double scaled=plist.dt_scale;
  double deltat_reac=scaled*epsilon*epsilon/(2.0*Dtot);

  double cf=cos(sqrt(4.0*wholep[1].Drx*deltat_reac));
  double Dr1=2.0*leglen2*(1.0-cf);

  Dtot+=2.0*Dr1/(6.0*deltat_reac);//add in rotation for both molecules
  cout <<"rank: "<<rank<<"add to Dtot from Rotation: "<<2.0*Dr1/(6.0*deltat_reac)<<" original dt: "<<deltat_reac<<" final Dtot: "<<Dtot<<endl;
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
  
  cout <<"rank: "<<rank<<" deltat: "<<deltat_reac<<" epsilon: "<<epsilon<<" tau: "<<tau<<" Rmax: "<<Rmax<<endl;
  double currtime=0;
  double Maxtime=plist.Maxtime;
  int Nitbin=int(Maxtime/deltat_reac);
  Maxtime=Nitbin*deltat_reac;
  cout <<"rank: "<<rank<<"number of time bins: "<<Nitbin<<" new max time: "<<Maxtime<<endl;
  double pirrev, pfree;
  //current statwrite is the number of datapoints it will write out.
  if(Nitbin<statwrite)
    statwrite=500;
  else
    statwrite=int(round(Nitbin/statwrite));
  
  cout <<"rank: "<<rank<<"statwrite: "<<statwrite<<endl;
  
  double theta;
  int Nthet=plist.Nthet;
  double deltheta=2.0/(1.0*(Nthet-1));
  cout <<"rank: "<<rank<<"Ntheta: "<<Nthet<<" delthet: "<<deltheta<<endl;
  double *costheta=new double[Nthet];
  for(i=0;i<Nthet;i++){
    costheta[i]=1-i*deltheta;
    cout <<"rank: "<<rank<<"costheta: "<<costheta[i]<<endl;
  }
  
  //int Nwrite=int(ceil(Nitbin/statwrite));
  int Nwrite=int(ceil(Nitbin/statwrite));
  cout <<"rank: "<<rank<<"rank: "<<rank<<'\t' <<"statwrite: "<<statwrite<<" Nwrite bins, for S(t): "<<Nwrite<<endl;
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
  
  double Rrange=3.0*sqrt(6.0*Dtot*Maxtime);//no +bindrad here because only measure distance beyond sigma.
  
  double delR=0.01;
  int Rbins=int(Rrange/delR);
  if(Rbins>1000){
    Rbins=1000;
    delR=(Rrange)/(1.0*Rbins);
  }
  cout <<"rank: "<<rank<<"rank: "<<rank<<'\t' <<"MaxR: "<<Rrange<<" Rbins: "<<Rbins<<" delR: "<<delR<<endl;
    
  int rind;
  int thetbins=plist.thetbins;
  int thetbins2=thetbins*thetbins;
  double small=1E-9;
  double delthet1=2.0/(1.0*thetbins)+small;
  cout <<"rank: "<<rank<<"theta bins: "<<thetbins<<" delthet1: "<<delthet1<<endl;
  ofstream thetfile("thetas.out");
  for(i=0;i<thetbins;i++)
    thetfile <<(i+0.5)*delthet1-1<<endl;
  
  int dihbins=plist.dihbins;
  double deldih=M_PI/(1.0*dihbins)+small;
  cout <<"rank: "<<rank<<"dihedral bins: "<<dihbins<< " delta: "<<deldih<<endl;
  
  double *prhist=new double[thetbins2*dihbins*Rbins];//
  double *prhist_sum=new double[thetbins2*dihbins*Rbins];//
  // for(i=0;i<thetbins2;i++)
//     prhist[i]=new double*[dihbins];
//   for(i=0;i<thetbins2;i++){
//     for(k=0;k<dihbins;k++){
//       prhist[i][k]=new double[Rbins];
//     }
//   }
  for(i=0;i<thetbins2;i++){
    for(k=0;k<dihbins;k++){
      for(j=0;j<Rbins;j++)
	prhist[i*thetbins2+k*dihbins+j]=0.0;
    }
  }
  double *p1hist=new double[Rbins];
  double *p1hist_sum=new double[Rbins];
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
  double *MA=new double[9];
  double *M=new double[9];
  double *MB=new double[9];
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

  double psimol, psi;
  double psi_limit;
  psi_limit=M_PI/3.0;//very accomodating on molecule orientation.
  //ofstream matchfile;
  int wcnt;
  int ct;
  int bt;
  int dd;
  int Npsibins=plist.Npsibins; 
  double delpsi=M_PI/(1.0*(Npsibins-1));
  double n1_x, n1_y, n1_z;
  double psi0;
  double *pn1=new double[3];
  double *par1=new double[3];
  double tol=1E-10;
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
  int ThetEnd=Nthet;
  double dihed=plist.dihed0;
  double dihlimit;//close to pi, but is more lenient for oriented legs.
  double anglelimit=0.1;
  ofstream sonlyfile;

  if(dihed>0){
    cts=1;//skip 1 again.
    ThetEnd=Nthet-1;//also -1 has only a single dihedral.
  }
  for(s=0;s<Nx0;s++){
    x0=x0vec[s];
    delR=(Rrange+x0-bindrad[0])/(1.0*Rbins);
    sprintf(tname, "survival_vs_angles_x0_%g_dh%g_dt%g.dat",x0, dihed, deltat_reac);
    sonlyfile.open(tname);
    for(ct=cts;ct<ThetEnd;ct++){
      lz=costheta[ct]*leglen;
      lx=sqrt(leglen2-lz*lz);
      ly=0;
      cout <<"rank: "<<rank<<"Initial pos for leg1 com: "<<lx<<' '<<ly<<' '<<x0+lz<<"  for leg: "<<0<<' '<<0<<' '<<x0<<endl;
      
      
      for(bt=ct;bt<ThetEnd;bt++){
	
	lz_0=-costheta[bt]*leglen;
	lR=sqrt(leglen2-lz_0*lz_0);
	lx_0=lR*cos(dihed);
	ly_0=lR*sin(dihed);
	cout <<"rank: "<<rank<<"Initial pos for leg2 com: "<<lx_0<<' '<<ly_0<<' '<<lz_0<<"  for leg: "<<0<<' '<<0<<' '<<0<<endl;
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
	cout <<"rank: "<<rank<<"Initial leg1 angle to the com divider: "<<ctab1<<" initial leg2 angle to the com divider: "<<ctab_0<<endl;
	for(dd=0;dd<Npsibins;dd++){
	  /*Assign orientation of molecules relative to one another (normals, not legs)*/
	  psi0=dd*delpsi;
	  /*Determine locations of the normal*/
	  /*n1 will be into the y plane, since the leg is fixed in the x-z plane
	   x and z coordinates are same as COM, length is unit*/
	  n1_x=lx;
	  n1_y=1.0;
	  n1_z=lz+x0;
	  bases[1].x[1]=n1_x;
	  bases[1].y[1]=n1_y;
	  bases[1].z[1]=n1_z;
	  v1[0]=0;//subtract off r0 of lx, 0, lz+x0
	  v1[1]=n1_y;
	  v1[2]=0;
	  /*now n0, vector is intersection of two planes, normal to leg, and 
	   normal to psi angle vector normal*/
	  bases[0].zcom=lz_0;
	  bases[0].xcom=lx_0;
	  bases[0].ycom=ly_0;
	  
	  bases[1].zcom=lz+x0;
	  bases[1].xcom=lx;
	  bases[1].ycom=ly;
	  
	  bases[0].x[0]=0;
	  bases[0].y[0]=0;
	  bases[0].z[0]=0;
	  bases[1].x[0]=0.0;
	  bases[1].y[0]=0.0;
	  bases[1].z[0]=x0;
	    
	  //COM to COM vector
	  dxm=bases[1].xcom-bases[0].xcom;
	  dym=bases[1].ycom-bases[0].ycom;
	  dzm=bases[1].zcom-bases[0].zcom;
	  R2=dxm*dxm+dym*dym+dzm*dzm;
	  R1=sqrt(R2);
	  //rotate norm around this vector by psi.
	  v[0]=dxm/R1;
	  v[1]=dym/R1;
	  v[2]=dzm/R1;
	  crossproduct(v, v1, par1);//par1 is normal to un-rotated
	  calc_Rmatrix(v, psi0, M);
	  rotate(v1, M, n1);//now n1 points in the direction of the normal for leg0 plane, fix to be perp to leg.
	  cout <<"rank: "<<rank<<"direction of rotated by psi mol angle; "<<n1[0]<<' '<<n1[1]<<' '<<n1[2]<<endl;
	  //then use that vector to set target plane via normal
	  crossproduct(v, n1, pn1);//pn1 is normal to plane need to intersect.
	  cout <<"rank: "<<rank<<"direction of normal to plane; "<<pn1[0]<<' '<<pn1[1]<<' '<<pn1[2]<<endl;
	  //move central point for that plane to bases[0].com
	  h1=pn1[0]*bases[0].xcom+pn1[1]*bases[0].ycom+pn1[2]*bases[0].zcom;

	  //plane 2 normal is leg_0 vector
	  r0x=bases[0].x[0]-bases[0].xcom;
	  r0y=bases[0].y[0]-bases[0].ycom;
	  r0z=bases[0].z[0]-bases[0].zcom;
	  //plane 2 point is bases[0].com
	  n2[0]=r0x/leglen;
	  n2[1]=r0y/leglen;
	  n2[2]=r0z/leglen;
	  parcheck=par1[0]*n2[0]+par1[1]*n2[1]+par1[2]*n2[2];
	  //n2 dot r0
	  h0=n2[0]*bases[0].xcom+n2[1]*bases[0].ycom+n2[2]*bases[0].zcom;
	  n1n2=pn1[0]*n2[0]+pn1[1]*n2[1]+pn1[2]*n2[2];
	  if(parcheck==1 || parcheck==-1){
	    //they are parallel
	    cout <<"rank: "<<rank<<"Parallel "<<endl;
	    //define second normal by rotating initial orient by psi around leg.
	    calc_Rmatrix(n2, psi0, M);
	    v1[0]=0;//subtract off r0 of lx, 0, lz+x0
	    v1[1]=n1_y;
	    v1[2]=0;
	    rotate(v1, M, n1);//now n1 points in the direction of the normal for leg0 plane, fix to be perp to leg.
	    p1x=n1[0];
	    p1y=n1[1];
	    p1z=n1[2];
	    R2=p1x*p1x+p1y*p1y+p1z*p1z;
	    p1len=sqrt(R2);
	    bases[0].x[1]=bases[0].xcom+p1x/p1len;
	    bases[0].y[1]=bases[0].ycom+p1y/p1len;
	    bases[0].z[1]=bases[0].zcom+p1z/p1len;
	    psimol=calc_psi(bases, v, v1, n1, n2);
	    
	  }else{
	    c1cof=(h1-h0*n1n2)/(1-n1n2*n1n2);
	    c0cof=(h0-h1*n1n2)/(1-n1n2*n1n2);
	    crossproduct(pn1, n2, v3);
	    p1x=c1cof*pn1[0]+c0cof*n2[0]+v3[0]-bases[0].xcom;
	    p1y=c1cof*pn1[1]+c0cof*n2[1]+v3[1]-bases[0].ycom;
	    p1z=c1cof*pn1[2]+c0cof*n2[2]+v3[2]-bases[0].zcom;
	    
	    R2=p1x*p1x+p1y*p1y+p1z*p1z;
	    p1len=sqrt(R2);
	    bases[0].x[1]=bases[0].xcom+p1x/p1len;
	    bases[0].y[1]=bases[0].ycom+p1y/p1len;
	    bases[0].z[1]=bases[0].zcom+p1z/p1len;
	    psimol=calc_psi(bases, v, v1, n1, n2);
	    if(abs(psi0-psimol)>tol){
	      psi1=psimol;
	      p1x*=-1;
	      p1y*=-1;
	      p1z*=-1;
	      bases[0].x[1]=bases[0].xcom+p1x/p1len;
	      bases[0].y[1]=bases[0].ycom+p1y/p1len;
	      bases[0].z[1]=bases[0].zcom+p1z/p1len;
	      psimol=calc_psi(bases, v, v1, n1, n2);
	      cout <<"rank: "<<rank<<"switch direction! original "<<psi1<<" new psi: "<<psimol<<" target "<<psi0<<endl;
	    }
	  }
	  for(i=0;i<Rbins;i++){
	    p1hist[i]=0.0;
	    p1hist_sum[i]=0.0;
	  }
	  
	  for(i=0;i<thetbins2;i++){
	    for(k=0;k<dihbins;k++){
	      for(j=0;j<Rbins;j++){
		prhist[i*thetbins2+k*dihbins+j]=0.0;
		prhist_sum[i*thetbins2+k*dihbins+j]=0.0;
	      }
	    }
	  }
	
	  for(i=0;i<Nwrite+1;i++){
	    phist[i]=0.0;
	    phist_sum[i]=0.0;
	  }
	  
	  cout <<"rank: "<<rank<<"X0: "<<x0<<" costheta: "<<costheta[ct]<<" costheta2: "<<costheta[bt]<<endl;
	  
	  cout <<"rank: "<<rank<<"p0 crds: "<<endl;
	  write_crds(bases, 0);
	  cout <<"rank: "<<rank<<"p1 crds: "<<endl;
	  write_crds(bases, 1);
	  dx=bases[1].x[0]-bases[0].x[0];
	  dy=bases[1].y[0]-bases[0].y[0];
	  dz=bases[1].z[0]-bases[0].z[0];
	  
	  R2=dx*dx+dy*dy+dz*dz;
	  R1=sqrt(R2);
	  
	  cout <<"rank: "<<rank<<"Initial separation: "<<R1<<endl;
	  cout <<"rank: "<<rank<<"target angles: "<<costheta[ct]<<' '<<costheta[bt]<<' '<<dihed<<endl;
	  calc_three_angle( cthet1,  cthet2,  dih,  bases,  dx, dy, dz,  R1, leglen,  v,  v1, n1, n2);
	  cout <<"rank: "<<rank<<" calculated angles: "<<cthet1<<' '<<cthet2<<' '<<dih<<endl;
	  /*calculate molecule dihedral*/
	  psimol=calc_psi(bases, v, v1, n1, n2);
	  cout <<"rank: "<<rank<<"target psi: "<<psi0<<" calculated: "<<psimol<<endl;
	  
	  /*Need to sample over different values of the sphere holding
	    all the points at leg length. origin of sphere is at z=x0, x=y=0,
	    azimuth doesn't matter, choose phi=0, so y=0 always to start.
	    then sample evenly over cos(theta),
	    x=d1*sin(theta), y=0, z=z0+d1*cos(theta).
	    You will perform rotation around this point.
	  */
	  
	  for(rep=0;rep<myNrep;rep++){
	    /*start proteins separated by x0 */
	    //cout <<"rank: "<<rank<<"Repeat: "<<rep<<endl;
	    
	    plist.ntotalcomplex=2;
	    Ntotalmol=2;
	    bases[0].x[0]=0;
	    bases[0].y[0]=0;
	    bases[0].z[0]=0;
	    /*p0_com does not move at all throughout sim.*/
	    bases[0].zcom=lz_0;
	    bases[0].xcom=lx_0;
	    bases[0].ycom=ly_0;
	    bases[0].x[1]=bases[0].xcom+p1x/p1len;
	    bases[0].y[1]=bases[0].ycom+p1y/p1len;
	    bases[0].z[1]=bases[0].zcom+p1z/p1len;
	    
	    bases[0].nfree=1;
	    bases[0].nbnd=0;
	    bases[0].npartner=0;
	    
	    
	    bases[1].x[0]=0.0;
	    bases[1].y[0]=0.0;
	    bases[1].z[0]=x0;
	    
	    bases[1].zcom=lz+x0;
	    bases[1].xcom=lx;
	    bases[1].ycom=ly;
	    
	    bases[1].x[1]=n1_x;
	    bases[1].y[1]=n1_y;
	    bases[1].z[1]=n1_z;
	    
	    bases[1].nfree=1;
	    bases[1].nbnd=0;
	    bases[1].npartner=0;
	    
	    
	    
// 	    dx=bases[1].x[0]-bases[0].x[0];
// 	    dy=bases[1].y[0]-bases[0].y[0];
// 	    dz=bases[1].z[0]-bases[0].z[0];
	    
// 	    R2=dx*dx+dy*dy+dz*dz;
// 	    R1=sqrt(R2);
	    
	    // /*Check the leg and the normal are perp*/
// 	    r0x=bases[0].x[0]-bases[0].xcom;
// 	    r0y=bases[0].y[0]-bases[0].ycom;
// 	    r0z=bases[0].z[0]-bases[0].zcom;
// 	    r1x=bases[0].x[1]-bases[0].xcom;
// 	    r1y=bases[0].y[1]-bases[0].ycom;
// 	    r1z=bases[0].z[1]-bases[0].zcom;
// 	    cout <<"rank: "<<rank<<"p0 dot prod, perp 0? : "<<r0x*r1x+r0y*r1y+r0z*r1z<<endl;
// 	    cout <<"rank: "<<rank<<"p0 dir length: "<<sqrt(r1x*r1x+r1y*r1y+r1z*r1z)<<endl;
// 	    r0x=bases[1].x[0]-bases[1].xcom;
// 	    r0y=bases[1].y[0]-bases[1].ycom;
// 	    r0z=bases[1].z[0]-bases[1].zcom;
// 	    r1x=bases[1].x[1]-bases[1].xcom;
// 	    r1y=bases[1].y[1]-bases[1].ycom;
// 	    r1z=bases[1].z[1]-bases[1].zcom;
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
		calc_three_angle( cthet1,  cthet2,  dih,  bases,  dx, dy, dz,  R1, leglen,  v,  v1, n1, n2);
		
		/*calculate molecule dihedral*/
		psimol=calc_psi(bases, v, v1, n1, n2);
				
		if(abs(cthet1-cthet2)<anglelimit){
		  dihlimit=abs(cthet1)*abs(cthet2)*M_PI;
		  if(M_PI-dih<dihlimit){
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
		  rotationEuler( tx,  ty,  tz,MA);
		  v[0]=bases[1].x[0]-bases[1].xcom;
		  v[1]=bases[1].y[0]-bases[1].ycom;
		  v[2]=bases[1].z[0]-bases[1].zcom;//pivot is the center
		  rotate(v, MA, v2);//includes the interface that will align
		  /*now rotate particles 0*/
		  tx=sqrt(2.0*deltat*wholep[0].Drx)*GaussV();
		  ty=sqrt(2.0*deltat*wholep[0].Dry)*GaussV();
		  tz=sqrt(2.0*deltat*wholep[0].Drz)*GaussV();
		  rotationEuler( tx,  ty,  tz,MB);
		  v1[0]=bases[0].x[0]-bases[0].xcom;
		  v1[1]=bases[0].y[0]-bases[0].ycom;
		  v1[2]=bases[0].z[0]-bases[0].zcom;//pivot is the center
		  rotate(v1, MB, v3);//includes the interface that will align
		  
		  
		  currx=(bases[1].xcom+dx+v2[0])-(bases[0].xcom+v3[0]);
		  curry=(bases[1].ycom+dy+v2[1])-(bases[0].ycom+v3[1]);
		  currz=(bases[1].zcom+dz+v2[2])-(bases[0].zcom+v3[2]);
		  
		  R2=currx*currx+curry*curry+currz*currz;
		  //  wcnt=0;
		  while(R2<bindrad[0]*bindrad[0]){
		    // wcnt++;
		    // 		  if(wcnt%1000==0){
		    // 		    cout <<"rank: "<<rank<<"In while loop stil "<<endl;
		    // 		    cout <<"rank: "<<rank<<"starting pos 0: "<<bases[0].xcom<<' '<<bases[0].ycom<<' '<<bases[0].zcom<<endl;
		    // 		    cout <<"rank: "<<rank<<"starting pos 0leg: "<<bases[0].x[0]<<' '<<bases[0].y[0]<<' '<<bases[0].z[0]<<endl;
		    // 		    cout <<"rank: "<<rank<<"starting pos 1: "<<bases[1].xcom<<' '<<bases[1].ycom<<' '<<bases[1].zcom<<endl;
		    // 		    cout <<"rank: "<<rank<<"starting pos 1leg: "<<bases[1].x[0]<<' '<<bases[1].y[0]<<' '<<bases[1].z[0]<<endl;
		    // 		    currx=(bases[1].x[0])-(bases[0].x[0]);
		    // 		    curry=(bases[1].y[0])-(bases[0].y[0]);
		    // 		    currz=(bases[1].z[0])-(bases[0].z[0]);
		    
		    // 		    R2=currx*currx+curry*curry+currz*currz;
		    // 		    cout <<"rank: "<<rank<<"starting sep: "<<sqrt(R2)<<endl;
		    
		    // 		    cout <<"rank: "<<rank<<"displacement: "<<dx<<' '<<dy<<' '<<dz<<endl;
		    // 		    cout <<"rank: "<<rank<<"rot vec 1: "<<v2[0]<<' '<<v2[1]<<' '<<v2[2]<<endl;
		    // 		    cout <<"rank: "<<rank<<"rot vec 2: "<<v3[0]<<' '<<v3[1]<<' '<<v3[2]<<endl;
		    
		    // 		    currx=(bases[1].xcom+dx+v2[0])-(bases[0].xcom+v3[0]);
		    // 		    curry=(bases[1].ycom+dy+v2[1])-(bases[0].ycom+v3[1]);
		    // 		    currz=(bases[1].zcom+dz+v2[2])-(bases[0].zcom+v3[2]);
		    
		    // 		    R2=currx*currx+curry*curry+currz*currz;
		    // 		    cout <<"rank: "<<rank<<"new sep: "<<sqrt(R2)<<endl;
		    // 		  }
		    
		    dx=sqrt(2.0*deltat*wholep[1].Dx)*GaussV();
		    dy=sqrt(2.0*deltat*wholep[1].Dy)*GaussV();
		    dz=sqrt(2.0*deltat*wholep[1].Dz)*GaussV();
		    tx=sqrt(2.0*deltat*wholep[1].Drx)*GaussV();
		    ty=sqrt(2.0*deltat*wholep[1].Dry)*GaussV();
		    tz=sqrt(2.0*deltat*wholep[1].Drz)*GaussV();
		    rotationEuler( tx,  ty,  tz,MA);
		    
		    rotate(v, MA, v2);//includes the interface that will align
		    
		    /*now rotate particles 0*/
		    tx=sqrt(2.0*deltat*wholep[0].Drx)*GaussV();
		    ty=sqrt(2.0*deltat*wholep[0].Dry)*GaussV();
		    tz=sqrt(2.0*deltat*wholep[0].Drz)*GaussV();
		    rotationEuler( tx,  ty,  tz,MB);
		    rotate(v1, MB, v3);//includes the interface that will align
		    
		    currx=(bases[1].xcom+dx+v2[0])-(bases[0].xcom+v3[0]);
		    curry=(bases[1].ycom+dy+v2[1])-(bases[0].ycom+v3[1]);
		    currz=(bases[1].zcom+dz+v2[2])-(bases[0].zcom+v3[2]);
		    
		    R2=currx*currx+curry*curry+currz*currz;
		    
		  }
		  /*update particle positions that do not overlap*/
		  //rotate direction vector, find orientation before motion
		  v[0]=bases[1].x[1]-bases[1].xcom;
		  v[1]=bases[1].y[1]-bases[1].ycom;
		  v[2]=bases[1].z[1]-bases[1].zcom;//pivot is the center
		  
		  bases[1].xcom+=dx;
		  bases[1].ycom+=dy;
		  bases[1].zcom+=dz;
		  //v2 currently is for i0
		  bases[1].x[0]=bases[1].xcom+v2[0];
		  bases[1].y[0]=bases[1].ycom+v2[1];
		  bases[1].z[0]=bases[1].zcom+v2[2];
		  //now v2 is for i1
		  rotate(v, MA, v2);//includes the interface that will align
		  bases[1].x[1]=bases[1].xcom+v2[0];
		  bases[1].y[1]=bases[1].ycom+v2[1];
		  bases[1].z[1]=bases[1].zcom+v2[2];
		  

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
		v[0]=bases[1].x[1]-bases[1].xcom;
		v[1]=bases[1].y[1]-bases[1].ycom;
		v[2]=bases[1].z[1]-bases[1].zcom;//pivot is the center
		
		bases[1].xcom+=dx;
		bases[1].ycom+=dy;
		bases[1].zcom+=dz;
		//v2 is for i0
		bases[1].x[0]=bases[1].xcom+v2[0];
		bases[1].y[0]=bases[1].ycom+v2[1];
		bases[1].z[0]=bases[1].zcom+v2[2];
		//now v2 is for i1
		rotate(v, M, v2);//includes the interface that will align
		bases[1].x[1]=bases[1].xcom+v2[0];
		bases[1].y[1]=bases[1].ycom+v2[1];
		bases[1].z[1]=bases[1].zcom+v2[2];
		
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
		  phist[ind]+=1;
		}
	      }
	    
	      
	      //cout <<"rank: "<<rank<<"rep: "<<rep<<" prev it: "<<it<<" itmorve: "<<itmove<<" time: "<<currtime<<" nfree: "<<bases[0].nfree<<" final pos: "<<bases[1].xcom<<' '<<bases[1].ycom<<' '<<bases[1].zcom<<endl;
	      it+=itmove;
	      
	    }
	    if(rep%10000==0)
	      cout <<"rank: "<<rank<<"finished rep: "<<rep<< " at time: "<<currtime<<endl;
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
		cout <<"rank: "<<rank<<"NAN: "<<dih<<" cosdih: "<<cosdih<<" bin: "<<dihbin<<endl;
	      }else if(dihbin==dihbins)
		dihbin=dihbins-1;
	    
	      //psimol=calc_psi(bases, v, v1, n1, n2);
	      //psibin=int(psimol/del1psi);
	      
	      if(ind_thet>=thetbins){
		// 	      cout <<"rank: "<<rank<<"cos(theta)>1 ?: "<<cthet1<<" ind: "<<ind_thet<<endl;
		ind_thet=thetbins-1;
	      }else if(ind_thet<0){
		// 	      cout <<"rank: "<<rank<<"cos(theta)<-1 ?: "<<cthet1<<" ind: "<<ind_thet<<endl;
		ind_thet=0;
	      }
	      if(rind<Rbins){
		prhist[(ind_thet*thetbins+ind_thet2)*thetbins2+dihbin*dihbins+rind]++;
		//prhist2[tab0ind*tabbins+tab1ind][rind]++;
		p1hist[rind]++;
	      }
	      
	    }
	  }//end over all reps
	  MPI::COMM_WORLD.Reduce(p1hist, p1hist_sum, Rbins, MPI::DOUBLE, MPI::SUM, 0);
	  MPI::COMM_WORLD.Reduce(phist, phist_sum, Nwrite+1, MPI::DOUBLE, MPI::SUM, 0);
	  MPI::COMM_WORLD.Reduce(prhist, prhist_sum, thetbins2*dihbins*Rbins, MPI::DOUBLE, MPI::SUM, 0);
	  MPI::COMM_WORLD.Barrier();
	  cout <<"rank: "<<rank<<"rank: "<<rank<<'\t' <<"finished all reps "<<endl;
	  
	  cout <<"rank: "<<rank<<"finished all reps "<<endl;
	  
	  /*Now write out final probability histogram*/
	  sprintf(tname, "rotprob_x0_%g_dh%g_ct%g_c2t_%g_psi%g_dt%g.dat",x0, dihed, costheta[ct], costheta[bt], psi0, deltat_reac);
	  if(rank==0){
	    cout <<"rank: "<<rank<<"open file: "<<tname<<endl;
	    probfile.open(tname);
	    probfile<<0<<' '<<1<<' '<<1<<endl;
	    //afile<<x0<<'\t';
	    for(i=statwrite;i<Nitbin;i+=statwrite){
	      tval=i*deltat_reac;
	      ind=i/statwrite;
	      passoc=survive_irr( x0, tval,  Dtot,  bindrad[0],  alpha,  cof);
	      probfile<<tval<<' ' <<phist_sum[ind]/(1.0*Nrep)<<' '<<1.0-passoc<<endl;
	      
	    }
	    i=Nitbin-1;
	    tval=i*deltat_reac;
	    passoc=survive_irr( x0, tval,  Dtot,  bindrad[0],  alpha,  cof);
	    probfile<<tval<<' ' <<phist_sum[Nwrite]/(1.0*Nrep)<<' '<<1.0-passoc<<endl;
	    probfile.close();
	    
	    passoc=survive_irr( x0, Maxtime,  Dtot,  bindrad[0],  alpha,  cof);
	    // sprintf(tname, "crt_sprob_x0_%g_dh0_%g_oat%g_obt%g_dt%g.dat",x0, dihed, ctab1, ctab_0, deltat_reac);
	    // 	tabfile.open(tname);
	    
	    sprintf(tname, "prt_sprob_x0_%g_dh0_%g_cat%g_cbt%g_psi%g_dt%g.dat",x0, dihed,  costheta[ct], costheta[bt], psi0, deltat_reac);
	    rfile.open(tname);
	    rfile<<0<<'\t'<<phist_sum[Nwrite]/(1.0*Nrep)<<'\t'<<1-passoc<<'\t'<<1-passoc<<endl;  
	    sonlyfile <<costheta[ct]<<'\t'<<costheta[bt]<<'\t'<<psi0<<'\t'<<phist_sum[Nwrite]/(1.0*Nrep)<<endl;
	    for(i=0;i<Rbins;i++){
	      R1=bindrad[0]+delR*(i+0.5);
	      pirrev=pirrev_value(R1, x0, Maxtime,  Dtot, bindrad[0],  alpha);
	      pfree=pfree_value_norm(R1, x0,  Maxtime,  Dtot, bindrad[0], alpha);
	      
	      rfile<<R1<<'\t'<<pirrev<<'\t'<<pfree*(1-passoc)<<'\t'<<p1hist_sum[i]/(1.0*Nrep*delR*R1*R1*4.0*M_PI)<<'\t';
	      for(k=0;k<dihbins;k++){
		for(j=0;j<thetbins2;j++){
		  
		  rfile<<prhist_sum[j*thetbins2+k*dihbins+i]/(1.0*Nrep*delR*R1*R1*delthet1*delthet1*deldih*4.0*M_PI)<<'\t';
		}
	      }//end dihedral bins
	      rfile<<endl;
	      
	      
	    }//end over r bins
	    rfile.close();
	  }//end rank=0;
	  MPI::COMM_WORLD.Barrier();
	
	}//end psi angles	
      }//end second angle
      
    }//end first angle
    sonlyfile.close();
  }//end over all x0s

  
  stop_timer(&totaltime);
  cout <<"rank: "<<rank<<timer_duration(totaltime)<<" total time "<<endl;  
  /*Write out final result*/
  cout <<"rank: "<<rank<<"End Main, complete run "<<endl;
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
  parmfile >>plist.Npsibins;
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
double calc_psi(Fullmol *bases, double *v, double *v1, double *n1, double *n2)
{
  /*the vectors that are read in are just place holders, they are overwritten
    by molecule geometry.
  */
  double dxm, dym, dzm;
  double r0x, r0y, r0z;
  double r1x, r1y, r1z;
  
  /*calculate molecule dihedral*/
  dxm=bases[1].xcom-bases[0].xcom;
  dym=bases[1].ycom-bases[0].ycom;
  dzm=bases[1].zcom-bases[0].zcom;
  
  r0x=bases[0].x[1]-bases[0].xcom;
  r0y=bases[0].y[1]-bases[0].ycom;
  r0z=bases[0].z[1]-bases[0].zcom;
  
  r1x=bases[1].x[1]-bases[1].xcom;
  r1y=bases[1].y[1]-bases[1].ycom;
  r1z=bases[1].z[1]-bases[1].zcom;
  
  
  //calculate dihedral angle.
  //n1=(1-0)com (x) (r1)
  v[0]=dxm;
  v[1]=dym;
  v[2]=dzm;
  v1[0]=r1x;
  v1[1]=r1y;
  v1[2]=r1z;
  double dp=v[0]*v1[0]+v[1]*v1[1]+v[2]*v1[2];
  double psimol;
  if(dp==1 || dp==-1){
    psimol=0;//by convention, if paralell set angle to zero.
    cout <<"parallel orient and com connector, set psi zero " <<endl;
  }else{
    crossproduct(v, v1, n1);
    //n2=(-r0) (x) (1-0)com
    v1[0]=-r0x;
    v1[1]=-r0y;
    v1[2]=-r0z;
    dp=v[0]*v1[0]+v[1]*v1[1]+v[2]*v1[2];
    if(dp==1 || dp==-1){
      psimol=0;
      cout <<"parallel orient and com connector, set psi zero " <<endl;
    }else{
      crossproduct(v1, v, n2);
      double cosdih=n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
      psimol=acos(cosdih);//if this is close to zero, molecules are oriented same
      if(isnan(psimol)){
	if(round(cosdih)==-1)psimol=M_PI;
	else psimol=0;
	cout <<"NAN: "<<psimol<<" cosdih: "<<cosdih<<endl;
      }  
      
    }
  }
  return psimol;
}
void calc_three_angle(double &cthet1, double &cthet2, double &dih, Fullmol *bases, double dx, double dy, double dz, double R1, double leglen, double *v, double *v1, double *n1, double *n2)
{
  
  double r0x, r0y, r0z;
  double r1x, r1y, r1z;
  
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
  
  //calculate dihedral angle.
  //n1=r1x (x) (0i-1i)
  v[0]=-dx;
  v[1]=-dy;
  v[2]=-dz;
  v1[0]=r1x;
  v1[1]=r1y;
  v1[2]=r1z;
  crossproduct(v1, v, n1);
  //n2= (0i-1i) (x) (-r0)
  v1[0]=-r0x;
  v1[1]=-r0y;
  v1[2]=-r0z;
  crossproduct(v, v1, n2);
  double cosdih=n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
  dih=acos(cosdih);
  if(cthet1==1 || cthet1==-1)
    dih=M_PI;
  else if(cthet2==1 || cthet2==-1)
    dih=M_PI;
 
 
}
void write_crds(Fullmol *bases, int p1)
{
  cout <<bases[p1].xcom<<' '<<bases[p1].ycom<<' '<<bases[p1].zcom<<endl;
  int i;
  for(i=0;i<bases[p1].ninterface;i++)
    cout <<bases[p1].x[i]<<' '<<bases[p1].y[i]<<' '<<bases[p1].z[i]<<endl;

}

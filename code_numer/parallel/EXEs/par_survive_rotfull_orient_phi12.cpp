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
#include "rd_class_defs.h"
#include "utility_calls.h"
#include "GF_calls.h"
#include "vector_rot_calls.h"
#include "angle_calc.h"

using namespace std;

struct MD_Timer totaltime;
struct MD_Timer bimoltime;


class survParms
{
public:
  double dihed0;
  int Nx0;
  int Nthet;
  int thetbins_out;
  int dihbins_out;
  int Npsibins;
  int Nphi;
  int psibins_out;

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



void read_parms(ifstream &parmfile, survParms &plist);
void write_crds(Fullmol *bases, int p1);


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
    intparm[6]=plist.Nx0;
    intparm[7]=plist.Nthet;
    intparm[8]=plist.thetbins_out;
    intparm[9]=plist.dihbins_out;
    intparm[10]=plist.Npsibins;
    intparm[11]=plist.psibins_out;
    intparm[12]=plist.Nphi;
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
  plist.thetbins_out=intparm[8];
  plist.dihbins_out=intparm[9];
  plist.Npsibins=intparm[10];
  plist.psibins_out=intparm[11];
  plist.Nphi=intparm[12];
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
  double deltheta;

  deltheta=2.0/(1.0*(Nthet-1));
  cout <<"rank: "<<rank<<"Ntheta: "<<Nthet<<" delthet: "<<deltheta<<endl;
  double *costheta=new double[Nthet];
  for(i=0;i<Nthet;i++){
    costheta[i]=1-i*deltheta;
    cout <<"rank: "<<rank<<"costheta: "<<costheta[i]<<endl;
  }
  costheta[0]=0.9;
  statwrite=1;
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
  int thetbins_out=plist.thetbins_out;
  int thetbins2_out=thetbins_out*thetbins_out;
  double small=1E-9;
  double delthet_out=2.0/(1.0*thetbins_out)+small;
  cout <<"rank: "<<rank<<"theta bins: "<<thetbins_out<<" delthet out: "<<delthet_out<<endl;
  ofstream thetfile("thetas.out");
  for(i=0;i<thetbins_out;i++)
    thetfile <<(i+0.5)*delthet_out-1<<endl;
  
  int dihbins_out=plist.dihbins_out;
  double deldih_out=M_PI/(1.0*dihbins_out)+small;
  cout <<"rank: "<<rank<<"dihedral bins: "<<dihbins_out<< " delta: "<<deldih_out<<endl;
  int psibins_out=plist.psibins_out;
  double delpsi_out=2.0*M_PI/(1.0*psibins_out)+small;
  int psiind;
  int k1, k2, j1, j2;
  double *prhist=new double[thetbins2_out*dihbins_out*psibins_out*Rbins];//
  double *prhist_sum=new double[thetbins2_out*dihbins_out*psibins_out*Rbins];//
  // for(i=0;i<thetbins2;i++)
//     prhist[i]=new double*[dihbins];
//   for(i=0;i<thetbins2;i++){
//     for(k=0;k<dihbins;k++){
//       prhist[i][k]=new double[Rbins];
//     }
//   }
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
  
  double *c0toc1=new double[3];
  double *c1toim=new double[3];
  double *c1toib=new double[3];
  double *c0toim=new double[3];
  double *c0toib=new double[3];
  double *rotv=new double[3];
  double *p0n=new double[3];
  double *p1n=new double[3];
  ofstream rfile;
  ofstream pavfile;

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
  int dd, cc;
  int Npsibins=plist.Npsibins; 
  double delpsi=2.0*M_PI/(1.0*(Npsibins));
  int Nphi=plist.Nphi;
  int Nphi2=plist.Npsibins;
  double phi1;
  double delphi1=2.0*M_PI/(1.0*Nphi);
  cout <<"Nphi: "<<Nphi<<" delphi: "<<delphi1<<endl;
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
  double norm_psi[Npsibins];
  double hist_add;
  double psi_jac;
  int ThetEnd=Nthet;
  double dihed=plist.dihed0;
  double dihlimit;//close to pi, but is more lenient for oriented legs.
  double anglelimit=0.6;
  ofstream sonlyfile;
  double delang;
  int ee;
  double phi2;
  // if(dihed>0){
//     cts=1;//skip 1 again.
//     ThetEnd=Nthet-1;//also -1 has only a single dihedral.
//   }
  for(s=0;s<Nx0;s++){
    x0=x0vec[s];
    delR=(Rrange+x0-bindrad[0])/(1.0*Rbins);
    sprintf(tname, "survival_vs_angles_x0_%g_dh%g_dt%g.dat",x0, dihed, deltat_reac);
    sonlyfile.open(tname);
    sprintf(pname, "pav_x0_%g_dh%g_dt%g.dat", x0, dihed, deltat);
    pavfile.open(pname);
    pavfile<<0<<'\t'<<0<<'\t'<<0<<'\t';
    for(i=0;i<Rbins;i++){
      R1=bindrad[0]+delR*(i+0.5);
      pavfile<<R1<<'\t';
    }
    pavfile<<endl;
    pavfile<<0<<'\t'<<0<<'\t'<<0<<'\t';
    for(i=0;i<Rbins;i++){
      R1=bindrad[0]+delR*(i+0.5);
      pirrev=pirrev_value(R1, x0, Maxtime,  Dtot, bindrad[0],  alpha);
      pavfile<<pirrev<<'\t';
    }
    pavfile<<endl;
    
    for(ct=cts;ct<ThetEnd;ct++){
      lz=costheta[ct]*leglen;
      lx=sqrt(leglen2-lz*lz);
      ly=0;
      //cout <<"rank: "<<rank<<"Initial pos for leg1 com: "<<lx<<' '<<ly<<' '<<x0+lz<<"  for leg: "<<0<<' '<<0<<' '<<x0<<endl;
      
      
      for(bt=ct;bt<ThetEnd;bt++){
	
	//lz_0=-costheta[bt]*leglen;
	lz_0=-0.8*leglen;
	lR=sqrt(leglen2-lz_0*lz_0);
	lx_0=lR*cos(dihed);
	ly_0=lR*sin(dihed);
	//cout <<"rank: "<<rank<<"Initial pos for leg2 com: "<<lx_0<<' '<<ly_0<<' '<<lz_0<<"  for leg: "<<0<<' '<<0<<' '<<0<<endl;
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
	//cout <<"rank: "<<rank<<"Initial leg1 angle to the com divider: "<<ctab1<<" initial leg2 angle to the com divider: "<<ctab_0<<endl;
	
	//for(dd=0;dd<Npsibins;dd++){
	  
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
	
	
	/*Assign orientation of molecules relative to one another (normals, not legs)*/
	for(ee=0;ee<Nphi2;ee++){
	  psi0=0;
	  for(i=0;i<Nwrite+1;i++){
	  phist_sum[i]=0.0;
	  phist[i]=0.0;
	}
	
	  norm_psi[ee]=0;
	  
	  /*First normal can be selected into the y plane, since the leg is fixed in the x-z plane
	    x and z coordinates are same as COM, length is unit*/
	  
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
	  
	  phi1=0;//cc*delphi1;
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
	  
	  /*COM to COM vector*/
	  dxm=bases[1].xcom-bases[0].xcom;
	  dym=bases[1].ycom-bases[0].ycom;
	  dzm=bases[1].zcom-bases[0].zcom;
	  R2=dxm*dxm+dym*dym+dzm*dzm;
	  R1=sqrt(R2);
	  
	  
	  c0toc1[0]=dxm/R1;
	  c0toc1[1]=dym/R1;
	  c0toc1[2]=dzm/R1;
	  
	  c1toib[0]=(bases[1].x[0]-bases[1].xcom)/leglen;
	  c1toib[1]=(bases[1].y[0]-bases[1].ycom)/leglen;
	  c1toib[2]=(bases[1].z[0]-bases[1].zcom)/leglen;
	  
	  /*rotate norm around COM-COM vector by psi to get plane for mol orientation vector of molecule 0 */
	  calc_Rmatrix(c0toc1, psi0, M);
	  rotate(c1toim, M, rotv);//now rotv points in the direction of the normal for leg0 plane, fix to be perp to leg.
	  
	  //cout <<"rank: "<<rank<<"direction of rotated by psi mol angle; "<<rotv[0]<<' '<<rotv[1]<<' '<<rotv[2]<<" original vector: "<<c1toim[0]<<' '<<c1toim[1]<<' '<<c1toim[2]<<endl;
	  
	  //then use that vector to set target plane via normal
	  crossproduct(c0toc1, rotv, p0n);//p0n is normal to plane need to intersect.
	  
	  //cout <<"rank: "<<rank<<"direction of normal to plane for molecule 0 orientation; "<<p0n[0]<<' '<<p0n[1]<<' '<<p0n[2]<<endl;
	  /*set central point for that plane to bases[0].com*/
	  h1=p0n[0]*bases[0].xcom+p0n[1]*bases[0].ycom+p0n[2]*bases[0].zcom;
	  
	  /*molecule 0 vector to binding interface*/
	  c0toib[0]=(bases[0].x[0]-bases[0].xcom)/leglen;
	  c0toib[1]=(bases[0].y[0]-bases[0].ycom)/leglen;
	  c0toib[2]=(bases[0].z[0]-bases[0].zcom)/leglen;
	  
	  /*c0toib dot to COM position vector*/
	  h0=c0toib[0]*bases[0].xcom+c0toib[1]*bases[0].ycom+c0toib[2]*bases[0].zcom;
	  n1n2=p0n[0]*c0toib[0]+p0n[1]*c0toib[1]+p0n[2]*c0toib[2];
	  
	  /*test to see if the norm we will looking for is in a plane, rather than a vector
	    at the intersection of 2 planes. 
	  */
	  crossproduct(c0toc1, c1toim, p1n);//p1n is normal to plane of un-rotated mol 1 orientation vector
	  parcheck=p1n[0]*c0toib[0]+p1n[1]*c0toib[1]+p1n[2]*c0toib[2];
	  if(parcheck==1 || parcheck==-1){
	    //they are parallel
	    cout <<"rank: "<<rank<<"Parallel plane and interface vector, mol 0. Psi is independent of rotation, but still need to sample the different leg angles."<<endl;
	    //define second normal by rotating initial orient by psi around leg.
	    calc_Rmatrix(c0toib, psi0, M);
	    
	    rotate(c1toim, M, rotv);//now n1 points in the direction of the normal for leg0 plane, fix to be perp to leg.
	    R2=rotv[0]*rotv[0]+rotv[1]*rotv[1]+rotv[2]*rotv[2];
	    R1=sqrt(R2);
	    c0toim[0]=rotv[0]/R1;
	    c0toim[1]=rotv[1]/R1;
	    c0toim[2]=rotv[2]/R1;
	    
	    bases[0].x[1]=bases[0].xcom+c0toim[0];
	    bases[0].y[1]=bases[0].ycom+c0toim[1];
	    bases[0].z[1]=bases[0].zcom+c0toim[2];
	    psimol=calc_psi_2pi(bases);
	    
	  }else{
	    c1cof=(h1-h0*n1n2)/(1-n1n2*n1n2);
	    c0cof=(h0-h1*n1n2)/(1-n1n2*n1n2);
	    crossproduct(p0n, c0toib, v3);
	    c0toim[0]=c1cof*p0n[0]+c0cof*c0toib[0]+v3[0]-bases[0].xcom;
	    c0toim[1]=c1cof*p0n[1]+c0cof*c0toib[1]+v3[1]-bases[0].ycom;
	    c0toim[2]=c1cof*p0n[2]+c0cof*c0toib[2]+v3[2]-bases[0].zcom;
	    /*this vector is through a line of points intersecting the two planes,
	      we choose an abritrary point on the line so the dihedral angle
	      could potentially end up rotated by pi, hence the tolerance
	      criterion.
	    */
	    R2=c0toim[0]*c0toim[0]+c0toim[1]*c0toim[1]+c0toim[2]*c0toim[2];
	    R1=sqrt(R2);
	    c0toim[0]/=R1;
	    c0toim[1]/=R1;
	    c0toim[2]/=R1;
	    
	    bases[0].x[1]=bases[0].xcom+c0toim[0];
	    bases[0].y[1]=bases[0].ycom+c0toim[1];
	    bases[0].z[1]=bases[0].zcom+c0toim[2];
	    psimol=calc_psi_2pi(bases);
	    if(abs(psi0-psimol)>tol){
	      psi1=psimol;
	      c0toim[0]*=-1;//point vector in opposite direction from bases[0].com
	      c0toim[1]*=-1;
	      c0toim[2]*=-1;
	      bases[0].x[1]=bases[0].xcom+c0toim[0];
	      bases[0].y[1]=bases[0].ycom+c0toim[1];
	      bases[0].z[1]=bases[0].zcom+c0toim[2];
	      
	      psimol=calc_psi_2pi(bases);
	      cout <<"rank: "<<rank<<"switch direction! original "<<psi1<<" new psi: "<<psimol<<" target "<<psi0<<endl;
	    }
	  }
	  /*Currently the im vectors are set at an angle of psi0=0, and c1im is along y axis.
	    first rotate c0im by phi2, which mimics changes in dihedral. Then loop over changes to phi1.
	  */
	  
	  //c0toib should be correct from above

	  
	  phi2=ee*delpsi;
	  calc_Rmatrix(c0toib, phi2, M);//as you return to each new phi2, you reinitialize legs to y axis and psi0.
	  cout <<"current phi2: "<<phi2<<endl;

	  rotate(c0toim, M, rotv);//now n1 points in the direction of the normal for leg0 plane, fix to be perp to leg.
	  R2=rotv[0]*rotv[0]+rotv[1]*rotv[1]+rotv[2]*rotv[2];
	  R1=sqrt(R2);
	  c0toim[0]=rotv[0]/R1;
	  c0toim[1]=rotv[1]/R1;
	  c0toim[2]=rotv[2]/R1;
	  
	  
	  for(cc=0;cc<Nphi;cc++){
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
	    
	    phi1=cc*delphi1;
	    /*This vector, c1toim,  should be rotated full 360 to sample all the DOF, around the c1toib vector.*/
	    if(cc>0)delang=delphi1;
	    else delang=0;
	    calc_Rmatrix(c1toib, delang, M);
	    rotate(c1toim, M, rotv);
	    c1toim[0]=rotv[0];
	    c1toim[1]=rotv[1];
	    c1toim[2]=rotv[2];
	    /*new positions below for molecule orientation end point*/
	    bases[1].x[1]=bases[1].xcom+c1toim[0];
	    bases[1].y[1]=bases[1].ycom+c1toim[1];
	    bases[1].z[1]=bases[1].zcom+c1toim[2];
	
	    psi_jac=1;//numer_calc_Jacobian(bases, psi0, leglen);//this returns |dphi2/dpsi|
	    hist_add=psi_jac*delphi1/(1.0*myNrep);
	    norm_psi[ee]+=psi_jac*delphi1;//this will end up being the jacobian factor for just integrating over psi.
	    cout <<" Jacobian, dphi2/dpsi: "<<psi_jac<<" at psi: "<<psi0<<" and phi1: "<<phi1<<endl;
	    cout <<"rank: "<<rank<<" X0: "<<x0<<" costheta: "<<costheta[ct]<<" costheta2: "<<costheta[bt]<<endl;
	    
	    cout <<"rank: "<<rank<<" p0 crds: "<<endl;
	    write_crds(bases, 0);
	    cout <<"rank: "<<rank<<" p1 crds: "<<endl;
	    write_crds(bases, 1);
	    dx=bases[1].x[0]-bases[0].x[0];
	    dy=bases[1].y[0]-bases[0].y[0];
	    dz=bases[1].z[0]-bases[0].z[0];
	    
	    R2=dx*dx+dy*dy+dz*dz;
	    R1=sqrt(R2);
	    
	    //cout <<"rank: "<<rank<<"Initial separation: "<<R1<<endl;
	    cout <<"rank: "<<rank<<"target angles: "<<costheta[ct]<<' '<<costheta[bt]<<' '<<dihed<<'\t';
	    calc_three_angle( cthet1,  cthet2,  dih,  bases,  dx, dy, dz,  R1, leglen,  v,  v1, n1, n2);
	    cout <<" calculated angles: "<<cthet1<<' '<<cthet2<<' '<<dih<<endl;
	    /*calculate molecule dihedral*/
	    psimol=calc_psi_2pi(bases);
	    cout <<"rank: "<<rank<<"target psi: "<<psi0<<" calculated: "<<psimol<<endl;
	    
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
	      bases[0].x[1]=bases[0].xcom+c0toim[0];
	      bases[0].y[1]=bases[0].ycom+c0toim[1];
	      bases[0].z[1]=bases[0].zcom+c0toim[2];
	      
	      bases[0].nfree=1;
	      bases[0].nbnd=0;
	      bases[0].npartner=0;
	      
	      
	      bases[1].x[0]=0.0;
	      bases[1].y[0]=0.0;
	      bases[1].z[0]=x0;
	      
	      bases[1].zcom=lz+x0;
	      bases[1].xcom=lx;
	      bases[1].ycom=ly;
	      
	      bases[1].x[1]=bases[1].xcom+c1toim[0];
	      bases[1].y[1]=bases[1].ycom+c1toim[1];
	      bases[1].z[1]=bases[1].zcom+c1toim[2];
	      
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
		  psimol=calc_psi_2pi(bases);
		  
		  if(abs(cthet1-cthet2)<anglelimit){
		    dihlimit=abs(cthet1)*abs(cthet2)*M_PI;
		    if(M_PI-dih<dihlimit){
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
		    
		    while(R2<bindrad[0]*bindrad[0]){
		      
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
		    //v2 currently is for i0, binding site
		    bases[1].x[0]=bases[1].xcom+v2[0];
		    bases[1].y[0]=bases[1].ycom+v2[1];
		    bases[1].z[0]=bases[1].zcom+v2[2];
		    //now v2 is for i1, the molecule orientation site
		    rotate(v, MA, v2);
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
		  //store this vector before moving the com.
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
		    phist[ind]+=hist_add;
		  }
		}
		
		
		//cout <<"rank: "<<rank<<"rep: "<<rep<<" prev it: "<<it<<" itmorve: "<<itmove<<" time: "<<currtime<<" nfree: "<<bases[0].nfree<<" final pos: "<<bases[1].xcom<<' '<<bases[1].ycom<<' '<<bases[1].zcom<<endl;
		it+=itmove;
		
	      }
	      if(rep+1%10000==0)
		cout <<"rank: "<<rank<<"finished rep: "<<rep<< " at time: "<<currtime<<endl;
	      /*Keep a histogram of final positions for particles that survived until the end*/
	      if(currtime>Maxtime){
		dx=bases[1].x[0]-bases[0].x[0];
		dy=bases[1].y[0]-bases[0].y[0];
		dz=bases[1].z[0]-bases[0].z[0];
		
		R2=dx*dx+dy*dy+dz*dz;
		R1=sqrt(R2);
		rind=int((R1-bindrad[0])/delR);
		//		cout <<"final sep: "<<R1<<endl;
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
		ind_thet=int((cthet1+1)/delthet_out);
		
		
		cthet2=r0x*dx+r0y*dy+r0z*dz;
		cthet2/=(R1*leglen);
		ind_thet2=int((cthet2+1)/delthet_out);
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
		  //cout <<"rank: "<<rank<<"NAN: "<<dih<<" cosdih: "<<cosdih<<" bin: "<<dihbin<<endl;
		}else if(dihbin==dihbins_out)
		  dihbin=dihbins_out-1;
		
		psimol=calc_psi_2pi(bases);
		psiind=int(psimol/delpsi_out);
		
		if(ind_thet>=thetbins_out){
		  // 	      cout <<"rank: "<<rank<<"cos(theta)>1 ?: "<<cthet1<<" ind: "<<ind_thet<<endl;
		  ind_thet=thetbins_out-1;
		}else if(ind_thet<0){
		  // 	      cout <<"rank: "<<rank<<"cos(theta)<-1 ?: "<<cthet1<<" ind: "<<ind_thet<<endl;
		  ind_thet=0;
		}
		if(rind<Rbins){
		  prhist[(ind_thet*thetbins_out+ind_thet2)*dihbins_out*psibins_out*Rbins+dihbin*psibins_out*Rbins+psiind*Rbins+rind]+=hist_add;
		  //prhist2[tab0ind*tabbins+tab1ind][rind]++;
		  p1hist[rind]+=hist_add;
		}
		
	      }
	    }//end over all reps
	    
	    cout <<"rank: "<<rank<<"finished all reps "<<" currhist: "<<phist[0]<<' '<<phist[1]<<endl;
	  }//end phi1 values.
	  MPI::COMM_WORLD.Reduce(p1hist, p1hist_sum, Rbins, MPI::DOUBLE, MPI::SUM, 0);
	  MPI::COMM_WORLD.Reduce(phist, phist_sum, Nwrite+1, MPI::DOUBLE, MPI::SUM, 0);
	  MPI::COMM_WORLD.Reduce(prhist, prhist_sum, thetbins2_out*dihbins_out*psibins_out*Rbins, MPI::DOUBLE, MPI::SUM, 0);
	  MPI::COMM_WORLD.Barrier();
	    
	  /*Now write out final probability histogram*/
	  sprintf(tname, "rotprob_x0_%g_dh%g_ct%g_c2t_%g_psi%g_dt%g.dat",x0, dihed, costheta[ct], costheta[bt], psi0, deltat_reac);
	  /*INSTEAD OF AVERAGING OVER NREP, NOW AVERAGE OVER NPROC  */
	  if(rank==0){
	    //cout <<"rank: "<<rank<<"open file: "<<tname<<endl;
	    probfile.open(tname);
	    probfile<<0<<' '<<1<<' '<<1<<endl;
	    //afile<<x0<<'\t';
	    for(i=statwrite;i<Nitbin;i+=statwrite){
	      tval=i*deltat_reac;
	      ind=i/statwrite;
	      passoc=survive_irr( x0, tval,  Dtot,  bindrad[0],  alpha,  cof);
	      probfile<<tval<<' ' <<phist_sum[ind]/(1.0*nprocs*norm_psi[ee])<<' '<<1.0-passoc<<endl;
	      
	    }
	    i=Nitbin-1;
	    tval=i*deltat_reac;
	    passoc=survive_irr( x0, tval,  Dtot,  bindrad[0],  alpha,  cof);
	    //probfile<<tval<<' ' <<phist_sum[Nwrite]/(1.0*nprocs*norm_psi[dd])<<' '<<1.0-passoc<<endl;
	    probfile.close();
	    
	    passoc=survive_irr( x0, Maxtime,  Dtot,  bindrad[0],  alpha,  cof);
	    // sprintf(tname, "crt_sprob_x0_%g_dh0_%g_oat%g_obt%g_dt%g.dat",x0, dihed, ctab1, ctab_0, deltat_reac);
	    // 	tabfile.open(tname);
	    /*Create a pavfile?*/
	    rfile<<0<<'\t'<<0<<'\t'<<0<<'\t'<<0<<'\t';
	    pavfile<<costheta[ct]<<'\t'<<costheta[bt]<<'\t'<<psi0<<'\t';
	    for(i=0;i<Rbins;i++){
	      R1=bindrad[0]+delR*(i+0.5);
	      pavfile<<p1hist_sum[i]/(1.0*nprocs*norm_psi[ee]*delR*R1*R1*4.0*M_PI)<<'\t';
	      rfile<<R1<<'\t';
	    }
	    pavfile<<endl;
	    rfile<<endl;
	

	    sprintf(tname, "prt_sprob_x0_%g_dh0_%g_cat%g_cbt%g_psi%g_dt%g.dat",x0, dihed,  costheta[ct], costheta[bt], psi0, deltat_reac);
	    rfile.open(tname);
	    //rfile<<0<<'\t'<<phist_sum[Nwrite]/(1.0*Nrep)<<'\t'<<1-passoc<<'\t'<<1-passoc<<endl;  
	    sonlyfile <<costheta[ct]<<'\t'<<costheta[bt]<<'\t'<<phi2<<'\t'<<phist_sum[Nwrite-1]/(1.0*nprocs*norm_psi[dd])<<'\t'<<norm_psi[dd]<<endl;
	    for(k1=0;k1<dihbins_out;k1++){
	      for(k2=0;k2<psibins_out;k2++){
		for(j1=0;j1<thetbins_out;j1++){
		  for(j2=0;j2<thetbins_out;j2++){
		    rfile<<(k1+0.5)*deldih_out<<'\t'<<(k2+0.5)*delpsi_out<<'\t'<<(j1+0.5)*delthet_out-1<<'\t'<<(j2+0.5)*delthet_out-1<<'\t';
		    j=j1*thetbins_out+j2;
		    
		    for(i=0;i<Rbins;i++){
		      R1=bindrad[0]+delR*(i+0.5);
		      rfile<<prhist_sum[j*dihbins_out*psibins_out*Rbins+k1*psibins_out*Rbins+k2*Rbins+i]/(1.0*nprocs*norm_psi[dd]*delR*R1*R1*delthet_out*delthet_out*deldih_out*delpsi_out*4.0*M_PI)<<'\t';
		    }
		    rfile<<endl;
		  }
		}
	      }//end psi
	    }//end dih chi

	    rfile.close();
	  }//end rank=0;
	  MPI::COMM_WORLD.Barrier();
	
	}//end phi2 values 
      }//end second angle
      
    }//end first angle
    sonlyfile.close();
    pavfile.close();
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
  parmfile >>plist.thetbins_out;
  parmfile.ignore(400,'\n');
  parmfile >>plist.dihbins_out;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Npsibins;
  parmfile.ignore(400,'\n');
  parmfile >>plist.psibins_out;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Nphi;
  parmfile.ignore(400,'\n');
  
}
void write_crds(Fullmol *bases, int p1)
{
  cout <<bases[p1].xcom<<' '<<bases[p1].ycom<<' '<<bases[p1].zcom<<endl;
  int i;
  for(i=0;i<bases[p1].ninterface;i++)
    cout <<bases[p1].x[i]<<' '<<bases[p1].y[i]<<' '<<bases[p1].z[i]<<endl;

}

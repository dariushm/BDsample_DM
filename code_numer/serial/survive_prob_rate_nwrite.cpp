/*

this program determines the survival probability
as a function of time averaged over multiple repetitions
for each set of initial separations and orientations.

This version survive_prob.cpp
only performs translational diffusion, should reproduce
exactly solvable solution.
  
*/
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "rand_gsl.h"
#include "md_timer.h"
#include "Faddeeva.hh"

#define MAXIFACE 10
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
  double x0;
  int Nx0;
  int Nprotypes;
  int Nifaces;
  double Nit;
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
double survive_irrF(double r0, double tcurr, double Dtot, double bindrad, double alpha, double cof);

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
  /*Keep center protein fixed, diffuse other protein with Dtot*/
  wholep[0].Dx=0;
  wholep[1].Dx=plist.D;
  wholep[0].Dy=0;
  wholep[1].Dy=plist.D;
  wholep[0].Dz=0;
  wholep[1].Dz=plist.D;

  bindrad[0]=1;
  kr[0]=plist.rate;
  cout <<"ACTIVATION RATE: "<<kr[0]<<"  radius: "<<bindrad[0]<<endl;
  

  for(i=0;i<Nprotypes;i++){  
    ntmp=wholep[i].ninterface+1;
    plist.Natom+=Ncopy[i]*ntmp;
  }
  cout <<"N atoms: "<<plist.Natom<<endl;
  cout <<"read reactions "<<endl;

  
  //double *savecrds=new double[Ntotalmol*3];//for x, y, z
  
  /*Print out specific reactions*/
  // cout <<"Print specific interaction network "<<endl;
//   int ncomplex=0;
//   for(i=0;i<Nifaces;i++){
//     cout <<i<<'\t';
//     for(j=0;j<numpartners[i];j++){
//       cout <<Speclist[i][j]<<'\t';
//       ncomplex++;
//     }
//     cout <<endl;
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

  int iind, iind2, ppart;
  
  int it;
  plist.ntotalcomplex=Ntotalmol;

  
  int s1;
  cout <<"Ntotal complexes: "<<plist.ntotalcomplex<<endl;

  int amol,df; 
  double us_to_s=1E-6;
  int statwrite=plist.statwrite;

  double kpi=4.0*M_PI*bindrad[0];
  double Dtot=wholep[0].Dx+wholep[1].Dx;
  double fourpi=4.0*M_PI;
  
  double R2, R1;
  
  double r0, passoc;
  double kact;
  double fact;
  double kdiff;//will be kpi*D
  double aexp;
  double bexp;
  double alpha;
  
  mu=0;
  kdiff=fourpi*Dtot*bindrad[mu];
  kact=kr[mu];
  fact=1.0+kact/kdiff;
  alpha=fact*sqrt(Dtot)/bindrad[mu];
  
  double cof=kr[0]/(kr[0]+kdiff);

  double Rmax;
  i=0;
    
  double currsep;
  double prevpassoc=0;
  double probvec1;

  cout <<" Dtot: "<<Dtot<<endl;

  double rnum2;

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
  double realsmall=1E-14;
  it=1;
  int rep;

  /*Number of times to repeat initial conditions*/
  int Nrep=plist.Nrep;
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
  
  /*Determine initial starting separations
    for rate calcluation, only do two starting separations, sigma and sigma+epsilon
  */
  Nx0=2;
  double *x0vec=new double[Nx0];
  x0vec[0]=bindrad[0];
  // for(i=1;i<Nx0;i++){
//     x0vec[i]=plist.x0+0.05*i;
//     cout <<"X0, i: "<<i<<" x0: "<<x0vec[i]<<endl;
//   }

  double kappa=kact/(4.0*M_PI*bindrad[0]*bindrad[0]);
  double epsilon=plist.eps_scale*Dtot/kappa;
  double tau=epsilon/kappa;
  x0vec[1]=bindrad[0]+epsilon;
  Rmax=bindrad[0]+epsilon;
  cout <<"Bindrad: "<<bindrad[0]<<" reaction limit: "<<Rmax<<endl;
  cout <<"time step expansions, using 3sqrt(6Dt): "<<epsilon*epsilon/Dtot/100.0<<endl;
  double scaled=plist.dt_scale;
  double deltat_reac=scaled*epsilon*epsilon/(2.0*Dtot);
  double deltat;
  cout <<" deltat: "<<deltat_reac<<" epsilon: "<<epsilon<<" tau: "<<tau<<endl;
  double currtime=0;
  double Maxtime=plist.Maxtime;
  int Nitbin=int(Maxtime/deltat_reac);
  Maxtime=Nitbin*deltat_reac;
  cout <<"number of time bins: "<<Nitbin<<" new max time: "<<Maxtime<<endl;
  double pirrev, pfree;
  
  if(Nitbin<statwrite)
    statwrite=500;
  else
    statwrite=int(round(Nitbin/statwrite));
  
  int Nwrite=int(ceil(Nitbin/statwrite));
  cout <<"statwrite: "<<statwrite<<" Nwrite bins, for S(t): "<<Nwrite<<endl;
  if((Nitbin-1)==Nwrite*statwrite)
    Nitbin-=1;
  double **phist=new double*[Nwrite+1];
  for(i=0;i<Nwrite+1;i++)
    phist[i]=new double[Nx0];
  
  
  for(i=0;i<Nwrite+1;i++){
    for(j=0;j<Nx0;j++)
      phist[i][j]=0.0;
  }
  int n;
  int itmove;
  int itmax;

  int Nacc=0;
  
  int cnt;
  double Rrange=3.0*sqrt(6.0*Dtot*Maxtime);
  cout <<"MaxR: "<<Rrange<<endl;
  
  double delR=0.01;
  int Rbins=int(Rrange/delR);
  if(Rbins>1000){
    Rbins=1000;
    delR=(Rrange)/(1.0*Rbins);
  }
  cout <<"MaxR: "<<Rrange<<" Rbins: "<<Rbins<<" delR: "<<delR<<endl;
  
  int rind;
  double *prhist=new double[Rbins];
  ofstream rfile;
  
  /*BEGIN SAMPLING OVER INITIAL PARTICLE SEPARATIONS*/  
  for(s=0;s<Nx0;s++){
    /*determine initial separation*/
    x0=x0vec[s];

    for(i=0;i<Rbins;i++)
      prhist[i]=0.0;

    for(i=0;i<Nwrite+1;i++)
      phist[i][s]=0.0;
    
    cout <<"X0: "<<x0<<endl;
    
    for(rep=0;rep<Nrep;rep++){
      /*start proteins separated by x0 */
      //cout <<"Repeat: "<<rep<<endl;

      plist.ntotalcomplex=2;
      Ntotalmol=2;
      bases[0].xcom=0;
      bases[0].ycom=0;
      bases[0].zcom=0;
      
      bases[0].nfree=1;
      bases[0].nbnd=0;
      bases[0].npartner=0;
      
      
      bases[1].xcom=0.0;
      bases[1].ycom=0.0;
      bases[1].zcom=x0;
      bases[1].nfree=1;
      bases[1].nbnd=0;
      bases[1].npartner=0;
      
      /***************************/
      /*Begin RD simulation*/
      it=0;
      currtime=0;
      //cout <<"rep: "<<rep<<" Maxtime; "<<Maxtime<<endl;
      while(bases[0].nfree>0 &&currtime<Maxtime)
	{
	  
	  
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
	    /*terminate some trajectories, using method of H.-X. Zhou.*/
	    deltat=deltat_reac;
	    prob=exp(-deltat_reac/tau);

	    rnum=1.0*rand_gsl();
	    p1=1;
	    if(prob<rnum){
	      //terminate trajectory
	      
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
	      /*propagate both particles and avoid overlap*/
	      /*Only translational diffusion*/
	      
	      dx=sqrt(2.0*deltat*wholep[1].Dx)*GaussV();
	      dy=sqrt(2.0*deltat*wholep[1].Dy)*GaussV();
	      dz=sqrt(2.0*deltat*wholep[1].Dz)*GaussV();
	      
	      currx=(bases[1].xcom+dx);//-(bases[0].xcom);
	      curry=(bases[1].ycom+dy);//-(bases[0].ycom);
	      currz=(bases[1].zcom+dz);//-(bases[0].zcom);
	      
	      R2=currx*currx+curry*curry+currz*currz;
	      
	      while(R2<bindrad[0]*bindrad[0]){
		
		dx=sqrt(2.0*deltat*wholep[1].Dx)*GaussV();
		dy=sqrt(2.0*deltat*wholep[1].Dy)*GaussV();
		dz=sqrt(2.0*deltat*wholep[1].Dz)*GaussV();
		
		currx=(bases[1].xcom+dx);//-(bases[0].xcom);
		curry=(bases[1].ycom+dy);//-(bases[0].ycom);
		currz=(bases[1].zcom+dz);//-(bases[0].zcom);
		
		R2=currx*currx+curry*curry+currz*currz;
		
	      }
	      /*update particle positions that do not overlap*/
	      bases[1].xcom+=dx;
	      bases[1].ycom+=dy;
	      bases[1].zcom+=dz;
	      	      
	      
	    }
	  }else{
	    /*free diffusion for both particles 
	      outside of 'reaction' zone
	    */
	    /*take larger steps at larger separations*/
	    deltat=scaled*(R1-bindrad[0])*(R1-bindrad[0])/(Dtot*2.0);//+deltat_reac;
	    itmove=int(deltat/deltat_reac);
	    deltat=itmove*deltat_reac;
	    dx=sqrt(2.0*deltat*wholep[1].Dx)*GaussV();
	    dy=sqrt(2.0*deltat*wholep[1].Dy)*GaussV();
	    dz=sqrt(2.0*deltat*wholep[1].Dz)*GaussV();
	    bases[1].xcom+=dx;
	    bases[1].ycom+=dy;
	    bases[1].zcom+=dz;
	    	    
	  }
	  
	  
	  currtime+=deltat;
	  
	  itmax=itmove+1;
	  if((itmax+it)>Nitbin)
	    itmax=Nitbin-it;
	  
	  for(n=1;n<itmax;n++){
	    if((it+n)%statwrite==0){
	      ind=(it+n)/statwrite;
	      phist[ind][s]+=1;
	    }
	  }
	  
	  it+=itmove;
	  
	}
      if(rep%10000==0)
	cout <<"finished rep: "<<rep<< " at time: "<<currtime<<endl;
      if(currtime>Maxtime){
	dx=bases[1].xcom-bases[0].xcom;
	dy=bases[1].ycom-bases[0].ycom;
	dz=bases[1].zcom-bases[0].zcom;
	
	R2=dx*dx+dy*dy+dz*dz;
	R1=sqrt(R2);
	rind=int((R1-bindrad[0])/delR);
	if(rind<Rbins)
	  prhist[rind]++;
	//else
	//cout <<"Beyond Rmax: "<<R1<<endl;
      }
    }//end over all reps
    cout <<"finished all reps "<<endl;
    /*Now write out final probability histogram*/
    sprintf(tname, "survive_prob_x0_%g_deltat%g.dat",x0, deltat_reac);
    cout <<"open file: "<<tname<<endl;
    probfile.open(tname);
    probfile<<0<<' '<<1<<' '<<1<<endl;
    //afile<<x0<<'\t';
    cnt=0;
    for(i=statwrite;i<Nitbin;i+=statwrite){
      tval=i*deltat_reac;
      
      ind=i/statwrite;
      passoc=survive_irr( x0, tval,  Dtot,  bindrad[0],  alpha,  cof);
      probfile<<tval<<' ' <<phist[ind][s]/(1.0*Nrep)<<' '<<1.0-passoc<<endl;

    }
    i=Nitbin-1;
    tval=i*deltat_reac;
    passoc=survive_irrF( x0, tval,  Dtot,  bindrad[0],  alpha,  cof);
    probfile<<tval<<' ' <<phist[Nwrite][s]/(1.0*Nrep)<<' '<<1.0-passoc<<endl;

    probfile.close();
    sprintf(tname, "prt_survive_prob_x0_%g_deltat%g.dat",x0, deltat_reac);
    rfile.open(tname);
    passoc=survive_irrF( x0, Maxtime,  Dtot,  bindrad[0],  alpha,  cof);
    rfile<<0<<'\t'<<phist[Nwrite][s]/(1.0*Nrep)<<'\t'<<1-passoc<<'\t'<<1-passoc<<endl;
    for(i=0;i<Rbins;i++){
      R1=bindrad[0]+delR*(i+0.5);
      pirrev=pirrev_value(R1, x0, Maxtime,  Dtot, bindrad[0],  alpha);
      pfree=pfree_value_norm(R1, x0,  Maxtime,  Dtot, bindrad[0], alpha);
      
      rfile<<R1<<'\t'<<prhist[i]/(1.0*Nrep*delR*R1*R1*4.0*M_PI)<<'\t'<<pirrev<<'\t'<<pfree*(1-passoc)<<endl;
    }
    rfile.close();
  }//end over all x0s
  x0=bindrad[0];
  double kt0=4.0*M_PI*bindrad[0]*bindrad[0]*kappa;
  double ikon=1.0/kt0+1.0/(4.0*M_PI*Dtot*bindrad[0]);
  double kon=1.0/ikon;
  sprintf(tname, "rate_vs_time_eps%g_dt%g_kon_%g.dat",epsilon, deltat_reac, kon);
  cout <<"open file: "<<tname<<endl;
  ofstream ratefile;
  ratefile.open(tname);
  
  cout <<"k(0): "<<kt0<<endl;
  ratefile<<0<<' '<<kt0<<' '<<kt0<<endl;
  for(i=statwrite;i<Nitbin;i+=statwrite){
    tval=i*deltat_reac;
	
    ind=i/statwrite;
	
    passoc=survive_irrF( x0, tval,  Dtot,  bindrad[0],  alpha,  cof);
    ratefile<<tval<<' ' <<kt0*0.5*(phist[ind][0]/(1.0*Nrep)+phist[ind][1]/(1.0*Nrep))<<' '<<kt0*(1.0-passoc)<<endl;
  }
  i=Nitbin-1;
  tval=i*deltat_reac;
  passoc=survive_irrF( x0, tval,  Dtot,  bindrad[0],  alpha,  cof);
  ratefile<<tval<<' ' <<kt0*0.5*(phist[Nwrite][0]/(1.0*Nrep)+phist[Nwrite][1]/(1.0*Nrep))<<' '<<kt0*(1.0-passoc)<<endl;
  
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
double survive_irrF(double r0, double tcurr, double Dtot, double bindrad, double alpha, double cof)
{
  double fDt=4.0*Dtot*tcurr;
  double sq_fDt=sqrt(fDt);

  double f1=cof*bindrad/r0;

  int i, j;
  double sep, dist;
  double sqrt_t=sqrt(tcurr);
  double a2=alpha*alpha;
  double r1, term1, term2, e1, ef1, sum;
  double onemsirr;
  sep=(r0-bindrad)/sq_fDt;//a

  e1=2.0*sep*sqrt_t*alpha+a2*tcurr;
  ef1=sep+alpha*sqrt_t; //a+b
  term1=erfc(sep);
  double ep1=exp(e1);
  if(isinf(ep1)){
    std::complex<double> z;
    real(z)=0.0;
    imag(z)=ef1;
    //cout <<"Complex number: "<<z<<endl;
    std::complex<double> value;
    double relerr=0;
    value=Faddeeva::w(z, relerr);
    double ea2=exp(-sep*sep);
    term2=ea2*real(value);
    
    
  }else{
    term2=exp(e1)*erfc(ef1);
  }
  sum=term1-term2;
  sum*=f1;
  onemsirr=sum;
  return onemsirr;//1-sirr=passoc
  //  cout <<"s_irr: "<<sirr<<" time: "<<tcurr<<" unscaled: "<<term1-term2<<endl;
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
  parmfile >>plist.x0;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Nx0;
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

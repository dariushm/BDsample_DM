/*


read in multiple, average  
  
*/
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>


using namespace std;



//void read_parms(ifstream &parmfile, Parms &plist);
double GaussV();
double pirr_pfree_ratio_ps(double rcurr, double r0, double tcurr, double Dtot, double bindrad, double alpha, double ps_prev, double rtol);
double survive_irr(double r0, double tcurr, double Dtot, double bindrad, double alpha, double cof);
double pirrev_value(double rcurr, double r0, double tcurr, double Dtot, double bindrad, double alpha);
double pfree_value_norm(double rcurr, double r0, double tcurr, double Dtot, double bindrad,double alpha);


int main(int argc, char *argv[])
{
  
  int i, j, k;
  

  
  char fname[300];
  ifstream probfile;
  int ig;
  int t=0;
  
  /*get the length from the first one*/
  i=0;
  double dummy;
  double sp_ex;
  int Nthet=atoi(argv[1]);
  double deltat=atof(argv[2]);
  double x0=atof(argv[3]);
  int thetbins1=atoi(argv[4]);
  int thetbins=thetbins1*thetbins1;
  cout <<"theta bins2: "<<thetbins<<endl;
  double del1thet=2.0/(1.0*thetbins1);
  double dih0=0;
  double dih1;
  int dih1bins=atoi(argv[5]);
  double deldih=M_PI/(1.0*dih1bins);
  dih1=0.5*deldih;
  double deltheta=2.0/(1.0*(Nthet-1));
  double *costheta=new double[Nthet];
  for(i=0;i<Nthet;i++)
    costheta[i]=1-i*deltheta;
  
  int jk;
  double cat0=1;
  double cbt0=1;
  int dir1=1;
  int dih0bins=atoi(argv[6]);
  char **dirlist=new char*[dih0bins];
  for(i=0;i<dih0bins;i++)
    dirlist[i]=new char[200];
  
  dirlist[0]="dih_0";
  dirlist[1]="dih_qt1";
  dirlist[2]="dih_half";
  dirlist[3]="dih_qt2";
  dirlist[4]="dih_pi";
  

  t=0;

  ofstream pfile;
  int maxbins=2000;
  double *pex=new double[maxbins];
  double *pfree=new double[maxbins];
  double *pravg=new double[maxbins];
  double *r=new double[maxbins];
  
  double psim;
  double psimthet;
  ifstream infile;
  double fact, thet1, thet2;
  char fnameout[300];
  int fnum=1;
  int Nruns=10;
  sprintf(fname, "%s/r%d/prt_sprob_x0_%g_dh0_%g_cat%g_cbt%g_dt%g.dat", dirlist[0], fnum, x0, dih0,  cat0, cbt0, deltat);
  probfile.open(fname);
  if(probfile.is_open())
    cout <<"able to open file ! "<<endl;
  else{
    while(!probfile.is_open() &&fnum<Nruns+1){
      cout <<"no file exist, try next num !"<<endl;
      fnum++;
      sprintf(fname, "%s/r%d/prt_sprob_x0_%g_dh0_%g_cat%g_cbt%g_dt%g.dat", dirlist[0], fnum, x0, dih0,  cat0, cbt0, deltat);
      probfile.open(fname);
    }
    if(fnum<Nruns+1)
      cout <<" opened file: "<<fname<<endl;
    else{
      exit(1);
      cerr<<" exiting, no file found "<<endl;
    }
  }
  double sps;
  probfile >>ig>>dummy>>dummy>>dummy;

  t=0;
  while(!probfile.eof()){
    probfile >>dummy>>dummy>>dummy>>dummy;
    for(i=0;i<thetbins*dih1bins;i++)
      probfile>>dummy;
    t++;
  }
  int rbins=t-1;
  cout <<"Rbins: "<<rbins<<endl;
  double **pirr_sim=new double*[thetbins];
  double **pirr_sim2=new double*[thetbins];
  for(i=0;i<thetbins;i++){
    pirr_sim[i]=new double[rbins];
    pirr_sim2[i]=new double[rbins];
  }
  double *dih0v=new double[dih0bins];
  for(i=0;i<dih0bins;i++){
    dih0v[i]=i*M_PI/(1.0*(dih0bins-1));
  }
  int m;
  int ct;
  int t1;
  int bt;
  int ct1;
  int Nacruns;
  int Nend;
  int c0=0;
  for(ct=c0;ct<Nthet;ct++){
    for(bt=ct;bt<Nthet;bt++){
      
      /*array to average over dihedral starts and finishes, for given starting angles*/
      for(i=0;i<thetbins;i++){
	for(t=0;t<rbins;t++){
	  pirr_sim[i][t]=0.0;
	  //pirr_sim2[i][t]=0.0;
	}
      }
  

      for(m=0;m<dih0bins;m++){

	sprintf(fnameout,"prt_trot2avg_x0_%g_dh0_%g_cat%g_cbt%g_dt%g.dat", x0, dih0v[m], costheta[ct], costheta[bt], deltat);
	infile.open(fnameout);
	
	
	cout <<"costheta: "<<costheta[ct]<<endl;
	
	
	if(infile.is_open()){
	  
	  infile >>ig>>ig>>ig;
	  
	  for(t=0;t<rbins;t++){
	    infile >>r[t] ;
	  }
	  for(jk=0;jk<dih1bins;jk++){
	    for(i=0;i<thetbins;i++){
	      
	      infile>>dih1>>thet1>>thet2;
	      for(t=0;t<rbins;t++){
		infile >>psim;
		//for each thetas and r, average over dihedral angles.
		fact=1;
		if(ct==0 || ct==Nthet-1)fact=fact*5;
		//if(ct!=bt)fact=fact*2;
		pirr_sim[i][t]+=psim*fact;//so each one can be divided by Ndih*dih1bins
		//pirr_sim2[k*thetbins+i][t]+=psim*psim;
	      }
	    }
	  }
	  infile.close();
	}else{
	  cout <<"COULD NOT OPEN: "<<fnameout<<endl;
	}
      }//end looping over dihedral bins.
      sprintf(fnameout,"dav_rt_x0_%g_cat%g_cbt%g_dt%g.dat", x0,  costheta[ct], costheta[bt], deltat);
      pfile.open(fnameout);
	
      pfile <<0<<'\t'<<0<<'\t';
      for(t=0;t<rbins;t++)
	pfile <<r[t]<<'\t';
      pfile<<endl;
	  
      for(i=0;i<thetbins1;i++){
	for(j=0;j<thetbins1;j++){
	  t1=i*thetbins1+j;
	  pfile <<(i+0.5)*del1thet-1<<'\t'<<(j+0.5)*del1thet-1<<'\t';
	  for(t=0;t<rbins;t++)
	    pfile <<pirr_sim[t1][t]/(1.0*dih0bins*dih1bins)<<'\t';
	  pfile<<endl;
	      
	}
      }
	
      pfile.close();
	
    }//end cbt
    
  }//end cat
  
  
  
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


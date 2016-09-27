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
  
  int Nruns=atoi(argv[1]);
  char fname[300];
  ifstream probfile;
  int ig;
  int t=0;
  
  /*get the length from the first one*/
  i=0;
  double dummy;
  double sp_ex;

  double deltat=atof(argv[2]);
  double x0=atof(argv[3]);
  
  int dir1=1;
  sprintf(fname, "r%d/survive_prob_x0_%g_deltat%g.dat", dir1, x0, deltat);
  probfile.open(fname);
  
  double zero, one;
  probfile >>zero>>one>>one;
  cout <<"open file: "<<fname<<" first line: "<<zero<<' '<<one<<' '<<one<<endl;
  t=0;

  ofstream pfile;
  int maxbins=35000;
  double *ps_ex=new double[maxbins];
  double *ps=new double[maxbins];
  double *time=new double[maxbins];
  
  double psim, pex;;
  while(!probfile.eof()){
    probfile >>time[t]>>psim>>pex;
    
    t++;
  }
  cout <<"tbins: "<<t-1<<endl;
  probfile.close();
  char fnameout[300];
  sprintf(fnameout,"Surv_trans_avg_x0_%g_deltat_%g.dat", x0, deltat);
  pfile.open(fnameout);
  
  int tbins=t-1;
  
  int ct;
  double *pirr_sim=new double[tbins];
  double *pirr_sim2=new double[tbins];
  double sp_sim, sp_avg;
  double sp_avg2;
  double pmore;
  ifstream infile;

  char sname[300];

  for(t=0;t<tbins;t++){
    pirr_sim[t]=0.0;
    pirr_sim2[t]=0.0;
  }

  for(i=1;i<Nruns+1;i++){
    sprintf(fname, "r%d/survive_prob_x0_%g_deltat%g.dat", dir1+i-1, x0, deltat);
    cout <<"file: "<<fname<<endl;
    infile.open(fname);
    
    for(t=0;t<tbins;t++){
      infile >>time[t] >>psim>>ps_ex[t];
      pirr_sim[t]+=psim;
      pirr_sim2[t]+=psim*psim;
      
    }
    infile.close();
  }//end runs
  double var;
  for(t=0;t<tbins;t++){
    pfile <<time[t]<<'\t'<<pirr_sim[t]/(1.0*Nruns)<<'\t';
    
    var=pirr_sim2[t]/(1.0*Nruns)-pirr_sim[t]/(1.0*Nruns)*pirr_sim[t]/(1.0*Nruns);
    pfile <<sqrt(var)<<endl;
  }
  


  
  
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


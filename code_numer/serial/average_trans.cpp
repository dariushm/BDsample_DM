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
  
  
  sprintf(fname, "r%d/prt_survive_prob_x0_%g_deltat%g.dat", 1, x0, deltat);
  probfile.open(fname);
  
  probfile >>ig>>dummy>>sp_ex>>dummy;
  cout <<"open file: "<<fname<<" first line: "<<ig<<' '<<sp_ex<<' '<<dummy<<' '<<dummy<<endl;
  t=0;
  cout <<"ASSUMES THAT THE SIMULATED SURVIVAL IS THE SECOND COLUMN, THEN EXACT, THEN FREE, VERIFY FROM PROGRAM! "<<endl
  ofstream pfile;
  int maxbins=1000;
  double *pex=new double[maxbins];
  double *pfree=new double[maxbins];
  double *r=new double[maxbins];
  
  double psim;
  while(!probfile.eof()){
    probfile >>r[t]>>psim>>pex[t]>>pfree[t];
    
    t++;
  }
  cout <<"rbins: "<<t-1<<endl;
  probfile.close();
  char fnameout[300];
  sprintf(fnameout,"prt_sp_trans_avg_x0_%g_deltat_%g.dat", x0, deltat);
  pfile.open(fnameout);
  
  int rbins=t-1;
  for(i=0;i<rbins;i++)
    pfile <<r[i]<<'\t';
  pfile <<endl;
  for(i=0;i<rbins;i++)
    pfile <<pex[i]<<'\t';
  pfile <<endl;
  for(i=0;i<rbins;i++)
    pfile <<pfree[i]<<'\t';
  pfile <<endl;
  
  int ct;
  double *pirr_sim=new double[rbins];
  double *pirr_sim2=new double[rbins];
  double sp_sim, sp_avg;
  double sp_avg2;
  double pmore;
  ifstream infile;
  ofstream sfile;
  char sname[300];
  sprintf(sname,"P_surviveTrans_x0_%g_dt%g.dat",x0,deltat);
  sfile.open(sname);
  sfile <<sp_ex<<'\t'<<"000"<<endl;

  for(t=0;t<rbins;t++){
    pirr_sim[t]=0.0;
    pirr_sim2[t]=0.0;
  }
  sp_avg=0.0;
  sp_avg2=0.0;
  for(i=1;i<Nruns+1;i++){
    sprintf(fname, "r%d/prt_survive_prob_x0_%g_deltat%g.dat", i, x0, deltat);
    cout <<"file: "<<fname<<endl;
    infile.open(fname);
    
    infile >>ig>>sp_sim>>sp_ex>>sp_ex;
    sp_avg+=sp_sim;
    sp_avg2+=sp_sim*sp_sim;
    for(t=0;t<rbins;t++){
      infile >>r[t] >>psim>>pmore>>pmore;
      pirr_sim[t]+=psim;
      pirr_sim2[t]+=psim*psim;
      
    }
    infile.close();
  }//end runs
  sp_avg/=(1.0*Nruns);
  sfile<<sp_avg<<'\t'<<sp_avg2/(1.0*Nruns)-sp_avg*sp_avg<<endl;
  for(t=0;t<rbins;t++)
    pfile <<pirr_sim[t]/(1.0*Nruns)<<'\t';
  pfile<<endl;
  for(t=0;t<rbins;t++)
    pfile <<pirr_sim2[t]/(1.0*Nruns)-pirr_sim[t]/(1.0*Nruns)*pirr_sim[t]/(1.0*Nruns)<<'\t';
  pfile<<endl;
  


  
  
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

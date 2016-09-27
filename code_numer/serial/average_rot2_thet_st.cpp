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
  int Nthet=atoi(argv[2]);
  double deltat=atof(argv[3]);
  double x0=atof(argv[4]);
  double dih0=0;
  double dih1;

  double deltheta=2.0/(1.0*(Nthet-1));
  double *costheta=new double[Nthet];
  for(i=0;i<Nthet;i++)
    costheta[i]=1-i*deltheta;
  
  int jk;
  double cat0=1;
  double cbt0=1;
  int dir1=1;
  int dih0bins=atoi(argv[5]);
  char **dirlist=new char*[dih0bins];
  for(i=0;i<dih0bins;i++)
    dirlist[i]=new char[200];
  
  dirlist[0]="dih_0";
  dirlist[1]="dih_qt1";
  dirlist[2]="dih_half";
  dirlist[3]="dih_qt2";
  dirlist[4]="dih_pi";
  
  int fnum=1;
  sprintf(fname, "%s/r%d/rot_sprob_x0_%g_dh%g_ct%g_c2t_%g_dt%g.dat", dirlist[0], fnum, x0, dih0,  cat0, cbt0, deltat);
  probfile.open(fname);
  if(probfile.is_open())
    cout <<"able to open file ! "<<endl;
  else{
    while(!probfile.is_open() &&fnum<Nruns+1){
      cout <<"no file exist, try next num !"<<endl;
      fnum++;
      sprintf(fname, "%s/r%d/rot_sprob_x0_%g_dh%g_ct%g_c2t_%g_dt%g.dat", dirlist[0], fnum, x0, dih0,  cat0, cbt0, deltat);
      probfile.open(fname);
    }
    if(fnum<Nruns+1)
      cout <<" opened file: "<<fname<<endl;
    else{
      exit(1);
      cerr<<" exiting, no file found "<<endl;
    }
  }
  double one;
  double tig;
  probfile >>tig>>one>>one;
  cout <<"open file: "<<fname<<" first line: "<<ig<<' '<<one<<' '<<one<<endl;
  t=1;//read in first time point, 0

  ofstream pfile;
  int maxbins=35000;
  double *time=new double[maxbins];

  double *pstrans=new double[maxbins];
  double psim, ps_ex;
  double psimthet;
  while(!probfile.eof()){
    probfile >>tig>>psim>>pstrans[t];

    t++;
  }
  cout <<"tbins: "<<t-1<<endl;
  probfile.close();
  char fnameout[300];

  double psend1, psend;
  int tbins=t-1;

  int ct;
  ifstream infile;
  ofstream sfile;
  char sname[300];
  int n;
  
  
  ofstream ravfile;
  sprintf(sname,"Psurv_vs_angle_x0_%g_dt%g.dat",x0,deltat);
  ravfile.open(sname);
  ofstream rav2file;
  sprintf(sname,"Var_Psurv_vs_angle_x0_%g_dt%g.dat",x0,deltat);
  rav2file.open(sname);
  
  double ***psav=new double**[dih0bins];
  double ***psav2=new double**[dih0bins];
  int m;
  for(m=0;m<dih0bins;m++){
    psav[m]=new double*[Nthet*Nthet];
    psav2[m]=new double*[Nthet*Nthet];
  }
  for(m=0;m<dih0bins;m++){
    for(k=0;k<Nthet*Nthet;k++){
      psav[m][k]=new double[tbins];
      psav2[m][k]=new double[tbins];
    }
  }
  for(m=0;m<dih0bins;m++){
    for(k=0;k<Nthet*Nthet;k++){
      for(i=0;i<tbins;i++){
	psav[m][k][i]=0.0;
	psav2[m][k][i]=0.0;
      }
    }
  }


  double *dih0v=new double[dih0bins];
  for(i=0;i<dih0bins;i++){
    dih0v[i]=i*M_PI/(1.0*(dih0bins-1));
  }

  int t1;
  int bt;
  int ct1;
  int Nacruns;
  int Nend;
  for(m=0;m<dih0bins;m++){
    ct1=0;
    Nend=Nthet;
    if(m>0){
      ct1=1;
      Nend=Nthet-1;
    }
    for(ct=ct1;ct<Nend;ct++){
      for(bt=ct;bt<Nend;bt++){
	//sprintf(fnameout,"Surv_rot2avg_x0_%g_dh0_%g_cat%g_cbt%g_dt%g.dat", x0, dih0v[m], costheta[ct], costheta[bt], deltat);
	//pfile.open(fnameout);

	// for(i=0;i<tbins;i++)
// 	  pfile <<time[i]<<'\t';
// 	pfile <<endl;
// 	for(t=0;t<tbins;t++){
// 	  pirr_sim[t]=0.0;
// 	  pirr_sim2[t]=0.0;
	  
// 	}
	
	Nacruns=Nruns;
	for(n=1;n<Nruns+1;n++){
	  sprintf(fname, "%s/r%d/rot_sprob_x0_%g_dh%g_ct%g_c2t_%g_dt%g.dat", dirlist[m], n+dir1-1, x0, dih0v[m],  costheta[ct], costheta[bt], deltat);
	  
	  infile.open(fname);
	  if(infile.is_open()){
	    // cout <<"file open: "<<fname<<endl;
	    	    
	    for(t=0;t<tbins-2;t++){
	      infile >>time[t] >>psim>>pstrans[t];
	      psav[m][ct*Nthet+bt][t]+=psim;
	      psav2[m][ct*Nthet+bt][t]+=psim*psim;
		
	    }
	    cout <<"curr t: "<<t<<endl;
	    infile >>time[t] >>psim>>pstrans[t];
	    psav[m][ct*Nthet+bt][t]+=psim*0.5;
	    psav2[m][ct*Nthet+bt][t]+=psim*psim*0.5*0.5;
	    t=t+1;
	    infile >>time[t] >>psim>>pstrans[t];
	    psav[m][ct*Nthet+bt][t]+=psim*0.5;
	    psav2[m][ct*Nthet+bt][t]+=psim*psim*0.5*0.5;
	    
	    infile.close();
	  }else{
	    cout <<"COULD NOT OPEN: "<<fname<<endl;
	    Nacruns--;
	  }
	}//end runs at this theta
	cout <<"costheta1: "<<costheta[ct]<<" costheta2: "<<costheta[bt]<<'\t';
	cout <<"total runs: "<<Nacruns<<" final fname: "<<fname<<endl;
	
		
	for(t=0;t<tbins;t++){
	  psav[m][ct*Nthet+bt][t]/=(1.0*Nacruns);
	  psav2[m][ct*Nthet+bt][t]/=(1.0*Nacruns);
	}
	//pfile.close();
	
      }//end cbt
    
    }//end cat
  
  }//end dih0
  
  int ind, mbin;
  double var, std;

  
  //write out to ravfile
  for(i=0;i<tbins;i++){
    ravfile<<time[i]<<'\t'<<pstrans[i]<<'\t';
    rav2file<<time[i]<<'\t'<<pstrans[i]<<'\t';
    /*For each dihedral angle, loop over all costheta1, costheta2.*/
    
    m=0;//dihedral
    for(ct=0;ct<Nthet;ct++){
      for(bt=0;bt<Nthet;bt++){
	ind=ct*Nthet+bt;
	if(bt<ct)
	  ind=bt*Nthet+ct;
	
	ravfile<<psav[m][ind][i]<<'\t';//only one dihedral angle 
	var=psav2[m][ind][i]-psav[m][ind][i]*psav[m][ind][i];
	rav2file<<sqrt(var)<<'\t';//only one dihedral angle 
      }
    }
    for(m=1;m<dih0bins;m++){
      for(ct=0;ct<Nthet;ct++){
	for(bt=0;bt<Nthet;bt++){
	  if(bt==0 || bt==Nthet-1)
	    mbin=0;
	  else if(ct==0 || ct==Nthet-1)
	    mbin=0;
	  else
	    mbin=m;
	      
	  ind=ct*Nthet+bt;
	  if(bt<ct)
	    ind=bt*Nthet+ct;
	  
	  ravfile<<psav[mbin][ind][i]<<'\t';//only one dihedral angle 
	  var=psav2[mbin][ind][i]-psav[mbin][ind][i]*psav[mbin][ind][i];
	  rav2file<<sqrt(var)<<'\t';//only one dihedral angle w
	}
      }
    }
    
    
    ravfile<<endl;
    rav2file<<endl;
  }
  //now write out survival probs.

  
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


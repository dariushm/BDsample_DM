#include "gridk_r2.h"



void gridk_r2(int i, int j, int l, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ****prt, double ****prt2, double sta, double delApol, double ct, double cta, double d1, double Drot, int nsmooth, int gbins)
{
  /*Cannot let r2 be equal to d, otherwise r1 will go to zero*/
  int k=0;
  
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  double delApol2=delApol*delApol;
  
  double c=-r*r+d1*d1;//r is r2!
  double omega=k*delaz;
  double b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
  double sfact=sqrt(b*b-4.0*c);
  double r1mag=-b/2.0+sfact/2.0;
  double r1sq=r1mag*r1mag;
  double derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
  double r2cube=r*r*r;
  double derv2=1.0/r-sfact*sfact/(4.0*r2cube);
  int knext=k+nsmooth;
  int kprev=k+nsmooth;
  double del1=nsmooth*delaz;
  /*az now goes from 0 to pi, so use reflecting BC.*/
  prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][knext]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0
  
  double jac=r/(r1mag+b/2.0);
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  int n;
  for(n=1;n<gbins-1;n++){
    k=n*nsmooth;
    if(n==gbins-2)
      knext=azbins-1;//
    else
      knext=nsmooth*(n+1);
    kprev=nsmooth*(n-1);
    omega=k*delaz;
    b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
    sfact=sqrt(b*b-4.0*c);
    r1mag=-b/2.0+sfact/2.0;
    r1sq=r1mag*r1mag;
    derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
    derv2=1.0/r-sfact*sfact/(4.0*r2cube);
    jac=r/(r1mag+b/2.0);
    prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

    //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][knext]-2.0*prt[i][j][k]+prt[i][j][k-1]) );
    psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
    //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }//done looping over az

  k=azbins-1;
  kprev=(gbins-2)*nsmooth;
  double del2=(k-kprev)*delaz;
  knext=kprev;
  omega=k*delaz;
   b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
   sfact=sqrt(b*b-4.0*c);
   r1mag=-b/2.0+sfact/2.0;
   r1sq=r1mag*r1mag;
   derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
   derv2=1.0/r-sfact*sfact/(4.0*r2cube);
   jac=r/(r1mag+b/2.0);
   prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del2)+ Drot*deltat/(sta*sta*del1*del2))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

   //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][k-1]) );//for i=inter, j=interior, k=end
   psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  
  /*Now interpolate*/
  double delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int ns=nsmooth;
  int t;
  double dpda;
  for(n=1;n<gbins;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    if(n==gbins-1){
      k=azbins-1;
      kprev=(n-1)*nsmooth;
      ns=k-kprev;
      delnext=ns*delaz;
    }
	       
    dpda=(prt2[i][j][l][k]-prt2[i][j][l][kprev])/delnext;
    for(t=0;t<ns-1;t++){
      omega=(kprev+t+1)*delaz;
      b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
      sfact=sqrt(b*b-4.0*c);
      r1mag=-b/2.0+sfact/2.0;
      jac=r/(r1mag+b/2.0);
    
      prt2[i][j][l][kprev+t+1]=prt2[i][j][l][kprev]+dpda*delaz*(t+1);
      psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
      
    }
  }


}

void gridk_r2_npole2(int i, int j, int l, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ****prt, double ****prt2, double sta, double delApol, double ct, double cta, double d1, double Drot, int nsmooth, int gbins)
{
  /*Cannot let r2 be equal to d, otherwise r1 will go to zero*/
  int k=0;
  
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  double delApol2=delApol*delApol;
  
  double c=-r*r+d1*d1;//r is r2!
  double omega=k*delaz;
  double b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
  double sfact=sqrt(b*b-4.0*c);
  double r1mag=-b/2.0+sfact/2.0;
  double r1sq=r1mag*r1mag;
  double derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
  double r2cube=r*r*r;
  double derv2=1.0/r-sfact*sfact/(4.0*r2cube);
  int kcross=azbins-1-k;
  int knext=k+nsmooth;
  int kprev=k+nsmooth;
  double del1=nsmooth*delaz;
  
  /*az now goes from 0 to pi, so use reflecting BC.*/
  prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l][kcross]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) ) ;//for i=inter, j=interior, k=0
  
  double jac=r/(r1mag+b/2.0);
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  int n;
  for(n=1;n<gbins-1;n++){
    k=n*nsmooth;
    if(n==gbins-2)
      knext=azbins-1;//
    else
      knext=nsmooth*(n+1);
    kprev=nsmooth*(n-1);
    
    omega=k*delaz;
    kcross=(gbins-1-n)*nsmooth;
    b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
    sfact=sqrt(b*b-4.0*c);
    r1mag=-b/2.0+sfact/2.0;
    r1sq=r1mag*r1mag;
    derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
    derv2=1.0/r-sfact*sfact/(4.0*r2cube);
    jac=r/(r1mag+b/2.0);
    prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l][kcross]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) ) ;//for i=inter, j=interior, k=0

    //prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][l][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][knext]-2.0*prt[i][j][l][k]+prt[i][j][k-1]) );
    psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
    //psurvive+=r2*delR*prt2[i][j][l][k]*delpol*delaz*st;
  }//done looping over az
   k=azbins-1;
   
   kcross=0;
   kprev=(gbins-2)*nsmooth;
   double del2=(k-kprev)*delaz;
   knext=kprev;
  
   omega=k*delaz;
   b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
   sfact=sqrt(b*b-4.0*c);
   r1mag=-b/2.0+sfact/2.0;
   r1sq=r1mag*r1mag;
   derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
   derv2=1.0/r-sfact*sfact/(4.0*r2cube);
   jac=r/(r1mag+b/2.0);
   prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del2)+ Drot*deltat/(sta*sta*del2*del1))*(prt[i][j][l][kprev]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l][kcross]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) ) ;//for i=inter, j=interior, k=0

   //prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][kprev]) );//for i=inter, j=interior, k=end
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  /*Now interpolate*/
  double delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int ns=nsmooth;
int t;
  double dpda;
  
  for(n=1;n<gbins;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    if(n==gbins-1){
      k=azbins-1;
      kprev=(n-1)*nsmooth;
      ns=k-kprev;
      delnext=ns*delaz;
    }
	       
    dpda=(prt2[i][j][l][k]-prt2[i][j][l][kprev])/delnext;
    for(t=0;t<ns-1;t++){
      omega=(kprev+t+1)*delaz;
      b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
      sfact=sqrt(b*b-4.0*c);
      r1mag=-b/2.0+sfact/2.0;
      jac=r/(r1mag+b/2.0);
    
      prt2[i][j][l][kprev+t+1]=prt2[i][j][l][kprev]+dpda*delaz*(t+1);
      psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
      
    }
  }

}

void gridk_r2_spole2(int i, int j, int l, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ****prt, double ****prt2, double sta, double delApol, double ct, double cta, double d1, double Drot, int nsmooth, int gbins)
{
  /*Cannot let r2 be equal to d, otherwise r1 will go to zero*/
  int k=0;
  
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  double delApol2=delApol*delApol;
  
  double c=-r*r+d1*d1;//r is r2!
  double omega=k*delaz;
  double b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
  double sfact=sqrt(b*b-4.0*c);
  double r1mag=-b/2.0+sfact/2.0;
  double r1sq=r1mag*r1mag;
  double derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
  double r2cube=r*r*r;
  double derv2=1.0/r-sfact*sfact/(4.0*r2cube);
  int kcross=azbins-1-k;
  int knext=k+nsmooth;
  int kprev=k+nsmooth;
  double del1=nsmooth*delaz;
  
  /*az now goes from 0 to pi, so use reflecting BC.*/
  prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][knext]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l][kcross]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0
  
  double jac=r/(r1mag+b/2.0);
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?


  int n;
  for(n=1;n<gbins-1;n++){
    k=n*nsmooth;
    if(n==gbins-2)
      knext=azbins-1;//
    else
      knext=nsmooth*(n+1);
    kprev=nsmooth*(n-1);
    
    omega=k*delaz;
    kcross=(gbins-1-n)*nsmooth;
  
    b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
    sfact=sqrt(b*b-4.0*c);
    r1mag=-b/2.0+sfact/2.0;
    r1sq=r1mag*r1mag;
    derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
    derv2=1.0/r-sfact*sfact/(4.0*r2cube);
    jac=r/(r1mag+b/2.0);
    prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l][kcross]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

    //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][knext]-2.0*prt[i][j][k]+prt[i][j][kprev]) );
    psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
    //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }//done looping over az
   k=azbins-1;
   
   kprev=(gbins-2)*nsmooth;
   double del2=(k-kprev)*delaz;
   knext=kprev;

   kcross=azbins-1-k;
   omega=k*delaz;
   b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
   sfact=sqrt(b*b-4.0*c);
   r1mag=-b/2.0+sfact/2.0;
   r1sq=r1mag*r1mag;
   derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
   derv2=1.0/r-sfact*sfact/(4.0*r2cube);
   jac=r/(r1mag+b/2.0);
   prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del2)+ Drot*deltat/(sta*sta*del1*del2))*(prt[i][j][l][kprev]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l][kcross]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

   //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][kprev]) );//for i=inter, j=interior, k=end
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  /*Now interpolate*/
  double delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int ns=nsmooth;
int t;
  double dpda;
  
  for(n=1;n<gbins;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    if(n==gbins-1){
      k=azbins-1;
      kprev=(n-1)*nsmooth;
      ns=k-kprev;
      delnext=ns*delaz;
    }
	       
    dpda=(prt2[i][j][l][k]-prt2[i][j][l][kprev])/delnext;
    for(t=0;t<ns-1;t++){
      omega=(kprev+t+1)*delaz;
      b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
      sfact=sqrt(b*b-4.0*c);
      r1mag=-b/2.0+sfact/2.0;
      jac=r/(r1mag+b/2.0);
    
      prt2[i][j][l][kprev+t+1]=prt2[i][j][l][kprev]+dpda*delaz*(t+1);
      psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
      
    }
  }

}
void gridk_r2_npole1(int i, int j, int l, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ****prt, double ****prt2, double sta, double delApol, double ct, double cta, double d1, double Drot, int nsmooth, int gbins)
{
  /*Cannot let r2 be equal to d, otherwise r1 will go to zero*/
  /*At npole, the previous bin in theta, is the same value of theta, but omega'=pi-omega.*/
  int k=0;
  
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  double delApol2=delApol*delApol;
  
  double c=-r*r+d1*d1;//r is r2!
  double omega=k*delaz;
  double b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
  double sfact=sqrt(b*b-4.0*c);
  double r1mag=-b/2.0+sfact/2.0;
  double r1sq=r1mag*r1mag;
  double derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
  double r2cube=r*r*r;
  double derv2=1.0/r-sfact*sfact/(4.0*r2cube);
  int kcross=azbins-k-1;
  int knext=k+nsmooth;
  int kprev=k+nsmooth;
  double del1=nsmooth*delaz;
  
  /*az now goes from 0 to pi, so use reflecting BC.*/
  prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j][l][kcross]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][knext]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0
  
  double jac=r/(r1mag+b/2.0);
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?


  int n;
  for(n=1;n<gbins-1;n++){
    k=n*nsmooth;
    if(n==gbins-2)
      knext=azbins-1;//
    else
      knext=nsmooth*(n+1);
    kprev=nsmooth*(n-1);
    
    omega=k*delaz;
    kcross=(gbins-1-n)*nsmooth;
  
    b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
    sfact=sqrt(b*b-4.0*c);
    r1mag=-b/2.0+sfact/2.0;
    r1sq=r1mag*r1mag;
    derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
    derv2=1.0/r-sfact*sfact/(4.0*r2cube);
    jac=r/(r1mag+b/2.0);
    prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j][l][kcross]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

    //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][knext]-2.0*prt[i][j][k]+prt[i][j][kprev]) );
    psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
    //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }//done looping over az
   k=azbins-1;
   kcross=0;
   kprev=(gbins-2)*nsmooth;
   double del2=(k-kprev)*delaz;
   knext=kprev;
  
   omega=k*delaz;
   b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
   sfact=sqrt(b*b-4.0*c);
   r1mag=-b/2.0+sfact/2.0;
   r1sq=r1mag*r1mag;
   derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
   derv2=1.0/r-sfact*sfact/(4.0*r2cube);
   jac=r/(r1mag+b/2.0);
   prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j][l][kcross]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del2)+ Drot*deltat/(sta*sta*del1*del2))*(prt[i][j][l][kprev]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

   //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*delaz2)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][kprev]) );//for i=inter, j=interior, k=end
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  /*Now interpolate*/
  double delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int ns=nsmooth;
int t;
  double dpda;
  
  for(n=1;n<gbins;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    if(n==gbins-1){
      k=azbins-1;
      kprev=(n-1)*nsmooth;
      ns=k-kprev;
      delnext=ns*delaz;
    }
	       
    dpda=(prt2[i][j][l][k]-prt2[i][j][l][kprev])/delnext;
    for(t=0;t<ns-1;t++){
      omega=(kprev+t+1)*delaz;
      b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
      sfact=sqrt(b*b-4.0*c);
      r1mag=-b/2.0+sfact/2.0;
      jac=r/(r1mag+b/2.0);
    
      prt2[i][j][l][kprev+t+1]=prt2[i][j][l][kprev]+dpda*delaz*(t+1);
      psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
      
    }
  }

}
void gridk_r2_npole1_npole2(int i, int j, int l, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ****prt, double ****prt2, double sta, double delApol, double ct, double cta, double d1, double Drot, int nsmooth, int gbins)
{
  /*Cannot let r2 be equal to d, otherwise r1 will go to zero*/
  /*At npole, the previous bin in theta, is the same value of theta, but omega'=pi-omega.*/
  int k=0;
  
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  double delApol2=delApol*delApol;
  
  double c=-r*r+d1*d1;//r is r2!
  double omega=k*delaz;
  double b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
  double sfact=sqrt(b*b-4.0*c);
  double r1mag=-b/2.0+sfact/2.0;
  double r1sq=r1mag*r1mag;
  double derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
  double r2cube=r*r*r;
  double derv2=1.0/r-sfact*sfact/(4.0*r2cube);
  int kcross=azbins-k-1;
  int knext=k+nsmooth;
  int kprev=k+nsmooth;
  double del1=nsmooth*delaz;
  
  /*az now goes from 0 to pi, so use reflecting BC.*/
  prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j][l][kcross]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][knext]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l][kcross]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) ) ;//for i=inter, j=interior, k=0
  
  double jac=r/(r1mag+b/2.0);
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?


  int n;
  for(n=1;n<gbins-1;n++){
    k=n*nsmooth;
    if(n==gbins-2)
      knext=azbins-1;//
    else
      knext=nsmooth*(n+1);
    kprev=nsmooth*(n-1);
    
    omega=k*delaz;
    kcross=(gbins-1-n)*nsmooth;
  
    b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
    sfact=sqrt(b*b-4.0*c);
    r1mag=-b/2.0+sfact/2.0;
    r1sq=r1mag*r1mag;
    derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
    derv2=1.0/r-sfact*sfact/(4.0*r2cube);
    jac=r/(r1mag+b/2.0);
    prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j][l][kcross]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l][kcross]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) ) ;//for i=inter, j=interior, k=0

    //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][knext]-2.0*prt[i][j][k]+prt[i][j][kprev]) );
    psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
    //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }//done looping over az
   k=azbins-1;
   kcross=azbins-k-1;
   kprev=(gbins-2)*nsmooth;
   double del2=(k-kprev)*delaz;
   knext=kprev;
  
   omega=k*delaz;
   b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
   sfact=sqrt(b*b-4.0*c);
   r1mag=-b/2.0+sfact/2.0;
   r1sq=r1mag*r1mag;
   derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
   derv2=1.0/r-sfact*sfact/(4.0*r2cube);
   jac=r/(r1mag+b/2.0);
   prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j][l][kcross]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del2)+ Drot*deltat/(sta*sta*del1*del2))*(prt[i][j][l][kprev]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l][kcross]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) ) ;//for i=inter, j=interior, k=0

   //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][kprev]) );//for i=inter, j=interior, k=end
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  /*Now interpolate*/
  double delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int ns=nsmooth;
int t;
  double dpda;
  
  for(n=1;n<gbins;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    if(n==gbins-1){
      k=azbins-1;
      kprev=(n-1)*nsmooth;
      ns=k-kprev;
      delnext=ns*delaz;
    }
	       
    dpda=(prt2[i][j][l][k]-prt2[i][j][l][kprev])/delnext;
    for(t=0;t<ns-1;t++){
      omega=(kprev+t+1)*delaz;
      b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
      sfact=sqrt(b*b-4.0*c);
      r1mag=-b/2.0+sfact/2.0;
      jac=r/(r1mag+b/2.0);
    
      prt2[i][j][l][kprev+t+1]=prt2[i][j][l][kprev]+dpda*delaz*(t+1);
      psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
      
    }
  }

}
void gridk_r2_npole1_spole2(int i, int j, int l, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ****prt, double ****prt2, double sta, double delApol, double ct, double cta, double d1, double Drot, int nsmooth, int gbins)
{
  /*Cannot let r2 be equal to d, otherwise r1 will go to zero*/
  /*At npole, the previous bin in theta, is the same value of theta, but omega'=pi-omega.*/
  int k=0;
  
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  double delApol2=delApol*delApol;
  
  double c=-r*r+d1*d1;//r is r2!
  double omega=k*delaz;
  double b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
  double sfact=sqrt(b*b-4.0*c);
  double r1mag=-b/2.0+sfact/2.0;
  double r1sq=r1mag*r1mag;
  double derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
  double r2cube=r*r*r;
  double derv2=1.0/r-sfact*sfact/(4.0*r2cube);
  int kcross=azbins-k-1;
  int knext=k+nsmooth;
  int kprev=k+nsmooth;
  double del1=nsmooth*delaz;
  
  /*az now goes from 0 to pi, so use reflecting BC.*/
  prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j][l][kcross]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][knext]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l][kcross]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0
  
  double jac=r/(r1mag+b/2.0);
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?


  int n;
  for(n=1;n<gbins-1;n++){
    k=n*nsmooth;
    if(n==gbins-2)
      knext=azbins-1;//
    else
      knext=nsmooth*(n+1);
    kprev=nsmooth*(n-1);
    
    omega=k*delaz;
    kcross=(gbins-1-n)*nsmooth;
   
    b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
    sfact=sqrt(b*b-4.0*c);
    r1mag=-b/2.0+sfact/2.0;
    r1sq=r1mag*r1mag;
    derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
    derv2=1.0/r-sfact*sfact/(4.0*r2cube);
    jac=r/(r1mag+b/2.0);
    prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j][l][kcross]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l][kcross]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

    //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][knext]-2.0*prt[i][j][k]+prt[i][j][kprev]) );
    psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
    //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }//done looping over az
   k=azbins-1;
   kcross=0;
   kprev=(gbins-2)*nsmooth;
   double del2=(k-kprev)*delaz;
   knext=kprev;
  
   omega=k*delaz;
   b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
   sfact=sqrt(b*b-4.0*c);
   r1mag=-b/2.0+sfact/2.0;
   r1sq=r1mag*r1mag;
   derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
   derv2=1.0/r-sfact*sfact/(4.0*r2cube);
   jac=r/(r1mag+b/2.0);
   prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j][l][kcross]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del2)+ Drot*deltat/(sta*sta*del1*del2))*(prt[i][j][l][kprev]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l][kcross]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

   //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][kprev]) );//for i=inter, j=interior, k=end
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  /*Now interpolate*/
  double delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int ns=nsmooth;
int t;
  double dpda;
  
  for(n=1;n<gbins;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    if(n==gbins-1){
      k=azbins-1;
      kprev=(n-1)*nsmooth;
      ns=k-kprev;
      delnext=ns*delaz;
    }
	       
    dpda=(prt2[i][j][l][k]-prt2[i][j][l][kprev])/delnext;
    for(t=0;t<ns-1;t++){
      omega=(kprev+t+1)*delaz;
      b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
      sfact=sqrt(b*b-4.0*c);
      r1mag=-b/2.0+sfact/2.0;
      jac=r/(r1mag+b/2.0);
    
      prt2[i][j][l][kprev+t+1]=prt2[i][j][l][kprev]+dpda*delaz*(t+1);
      psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
      
    }
  }

}
void gridk_r2_spole1(int i, int j, int l, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ****prt, double ****prt2, double sta, double delApol, double ct, double cta, double d1, double Drot, int nsmooth, int gbins)
{
  /*Cannot let r2 be equal to d, otherwise r1 will go to zero*/
  /*At npole, the previous bin in theta, is the same value of theta, but omega'=pi-omega.*/
  int k=0;
  
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  double delApol2=delApol*delApol;
  
  double c=-r*r+d1*d1;//r is r2!
  double omega=k*delaz;
  double b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
  double sfact=sqrt(b*b-4.0*c);
  double r1mag=-b/2.0+sfact/2.0;
  double r1sq=r1mag*r1mag;
  double derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
  double r2cube=r*r*r;
  double derv2=1.0/r-sfact*sfact/(4.0*r2cube);
  int kcross=azbins-k-1;
  int knext=k+nsmooth;
  int kprev=k+nsmooth;
  double del1=nsmooth*delaz;
  
  /*az now goes from 0 to pi, so use reflecting BC.*/
  prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j][l][kcross]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][knext]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0
  
  double jac=r/(r1mag+b/2.0);
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?


  int n;
  for(n=1;n<gbins-1;n++){
    k=n*nsmooth;
    if(n==gbins-2)
      knext=azbins-1;//
    else
      knext=nsmooth*(n+1);
    kprev=nsmooth*(n-1);
    
    omega=k*delaz;
    kcross=(gbins-1-n)*nsmooth;
  
    b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
    sfact=sqrt(b*b-4.0*c);
    r1mag=-b/2.0+sfact/2.0;
    r1sq=r1mag*r1mag;
    derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
    derv2=1.0/r-sfact*sfact/(4.0*r2cube);
    jac=r/(r1mag+b/2.0);
    prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j][l][kcross]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

    //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][knext]-2.0*prt[i][j][k]+prt[i][j][kprev]) );
    psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
    //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }//done looping over az
   k=azbins-1;
   kcross=0;
   kprev=(gbins-2)*nsmooth;
   double del2=(k-kprev)*delaz;
   knext=kprev;
  
   omega=k*delaz;
   b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
   sfact=sqrt(b*b-4.0*c);
   r1mag=-b/2.0+sfact/2.0;
   r1sq=r1mag*r1mag;
   derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
   derv2=1.0/r-sfact*sfact/(4.0*r2cube);
   jac=r/(r1mag+b/2.0);
   prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j][l][kcross]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del2)+ Drot*deltat/(sta*sta*del1*del2))*(prt[i][j][l][kprev]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

   //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][kprev]) );//for i=inter, j=interior, k=end
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  /*Now interpolate*/
  double delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int ns=nsmooth;
int t;
  double dpda;
  
  for(n=1;n<gbins;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    if(n==gbins-1){
      k=azbins-1;
      kprev=(n-1)*nsmooth;
      ns=k-kprev;
      delnext=ns*delaz;
    }
	       
    dpda=(prt2[i][j][l][k]-prt2[i][j][l][kprev])/delnext;
    for(t=0;t<ns-1;t++){
      omega=(kprev+t+1)*delaz;
      b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
      sfact=sqrt(b*b-4.0*c);
      r1mag=-b/2.0+sfact/2.0;
      jac=r/(r1mag+b/2.0);
    
      prt2[i][j][l][kprev+t+1]=prt2[i][j][l][kprev]+dpda*delaz*(t+1);
      psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
      
    }
  }

}
void gridk_r2_spole1_npole2(int i, int j, int l, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ****prt, double ****prt2, double sta, double delApol, double ct, double cta, double d1, double Drot, int nsmooth, int gbins)
{
  /*Cannot let r2 be equal to d, otherwise r1 will go to zero*/
  /*At npole, the previous bin in theta, is the same value of theta, but omega'=pi-omega.*/
  int k=0;
  
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  double delApol2=delApol*delApol;
  
  double c=-r*r+d1*d1;//r is r2!
  double omega=k*delaz;
  double b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
  double sfact=sqrt(b*b-4.0*c);
  double r1mag=-b/2.0+sfact/2.0;
  double r1sq=r1mag*r1mag;
  double derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
  double r2cube=r*r*r;
  double derv2=1.0/r-sfact*sfact/(4.0*r2cube);
  int kcross=azbins-k-1;
  int knext=k+nsmooth;
  int kprev=k+nsmooth;
  double del1=nsmooth*delaz;
  
  /*az now goes from 0 to pi, so use reflecting BC.*/
  prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j][l][kcross]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][knext]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l][kcross]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) ) ;//for i=inter, j=interior, k=0
  
  double jac=r/(r1mag+b/2.0);
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?


  int n;
  for(n=1;n<gbins-1;n++){
    k=n*nsmooth;
    if(n==gbins-2)
      knext=azbins-1;//
    else
      knext=nsmooth*(n+1);
    kprev=nsmooth*(n-1);
    
    omega=k*delaz;
    kcross=(gbins-1-n)*nsmooth;
  
    b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
    sfact=sqrt(b*b-4.0*c);
    r1mag=-b/2.0+sfact/2.0;
    r1sq=r1mag*r1mag;
    derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
    derv2=1.0/r-sfact*sfact/(4.0*r2cube);
    jac=r/(r1mag+b/2.0);
    prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j][l][kcross]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l][kcross]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) ) ;//for i=inter, j=interior, k=0

    //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][knext]-2.0*prt[i][j][k]+prt[i][j][kprev]) );
    psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
    //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }//done looping over az
   k=azbins-1;
   kcross=azbins-k-1;
   kprev=(gbins-2)*nsmooth;
   double del2=(k-kprev)*delaz;
   knext=kprev;
  
   omega=k*delaz;
   b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
   sfact=sqrt(b*b-4.0*c);
   r1mag=-b/2.0+sfact/2.0;
   r1sq=r1mag*r1mag;
   derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
   derv2=1.0/r-sfact*sfact/(4.0*r2cube);
   jac=r/(r1mag+b/2.0);
   prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j][l][kcross]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del2)+ Drot*deltat/(sta*sta*del1*del2))*(prt[i][j][l][kprev]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l][kcross]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) ) ;//for i=inter, j=interior, k=0

   //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][kprev]) );//for i=inter, j=interior, k=end
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  /*Now interpolate*/
  double delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int ns=nsmooth;
int t;
  double dpda;
  
  for(n=1;n<gbins;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    if(n==gbins-1){
      k=azbins-1;
      kprev=(n-1)*nsmooth;
      ns=k-kprev;
      delnext=ns*delaz;
    }
	       
    dpda=(prt2[i][j][l][k]-prt2[i][j][l][kprev])/delnext;
    for(t=0;t<ns-1;t++){
      omega=(kprev+t+1)*delaz;
      b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
      sfact=sqrt(b*b-4.0*c);
      r1mag=-b/2.0+sfact/2.0;
      jac=r/(r1mag+b/2.0);
    
      prt2[i][j][l][kprev+t+1]=prt2[i][j][l][kprev]+dpda*delaz*(t+1);
      psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
      
    }
  }

}
void gridk_r2_spole1_spole2(int i, int j, int l, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ****prt, double ****prt2, double sta, double delApol, double ct, double cta, double d1, double Drot, int nsmooth, int gbins)
{
  /*Cannot let r2 be equal to d, otherwise r1 will go to zero*/
  /*At npole, the previous bin in theta, is the same value of theta, but omega'=pi-omega.*/
  int k=0;
  
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  double delApol2=delApol*delApol;
  
  double c=-r*r+d1*d1;//r is r2!
  double omega=k*delaz;
  double b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
  double sfact=sqrt(b*b-4.0*c);
  double r1mag=-b/2.0+sfact/2.0;
  double r1sq=r1mag*r1mag;
  double derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
  double r2cube=r*r*r;
  double derv2=1.0/r-sfact*sfact/(4.0*r2cube);
  int kcross=azbins-k-1;
  int knext=k+nsmooth;
  int kprev=k+nsmooth;
  double del1=nsmooth*delaz;
  
  /*az now goes from 0 to pi, so use reflecting BC.*/
  prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j][l][kcross]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][knext]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l][kcross]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0
  
  double jac=r/(r1mag+b/2.0);
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?


  int n;
  for(n=1;n<gbins-1;n++){
    k=n*nsmooth;
    if(n==gbins-2)
      knext=azbins-1;//
    else
      knext=nsmooth*(n+1);
    kprev=nsmooth*(n-1);
    
    omega=k*delaz;
    kcross=(gbins-1-n)*nsmooth;
  
      b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
    sfact=sqrt(b*b-4.0*c);
    r1mag=-b/2.0+sfact/2.0;
    r1sq=r1mag*r1mag;
    derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
    derv2=1.0/r-sfact*sfact/(4.0*r2cube);
    jac=r/(r1mag+b/2.0);
    prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j][l][kcross]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l][kcross]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

    //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][knext]-2.0*prt[i][j][k]+prt[i][j][kprev]) );
    psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
    //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }//done looping over az
   k=azbins-1;
   kcross=0;
   kprev=(gbins-2)*nsmooth;
   double del2=(k-kprev)*delaz;
   knext=kprev;
  
   omega=k*delaz;
   b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
   sfact=sqrt(b*b-4.0*c);
   r1mag=-b/2.0+sfact/2.0;
   r1sq=r1mag*r1mag;
   derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
   derv2=1.0/r-sfact*sfact/(4.0*r2cube);
   jac=r/(r1mag+b/2.0);
   prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-prt[i-1][j][l][k])+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+prt[i-1][j][l][k])+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j][l][kcross]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del2)+ Drot*deltat/(sta*sta*del1*del2))*(prt[i][j][l][kprev]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l][kcross]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

   //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][kprev]) );//for i=inter, j=interior, k=end
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  /*Now interpolate*/
  double delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int ns=nsmooth;
int t;
  double dpda;
  
  for(n=1;n<gbins;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    if(n==gbins-1){
      k=azbins-1;
      kprev=(n-1)*nsmooth;
      ns=k-kprev;
      delnext=ns*delaz;
    }
	       
    dpda=(prt2[i][j][l][k]-prt2[i][j][l][kprev])/delnext;
    for(t=0;t<ns-1;t++){
      omega=(kprev+t+1)*delaz;
      b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
      sfact=sqrt(b*b-4.0*c);
      r1mag=-b/2.0+sfact/2.0;
      jac=r/(r1mag+b/2.0);
    
      prt2[i][j][l][kprev+t+1]=prt2[i][j][l][kprev]+dpda*delaz*(t+1);
      psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
      
    }
  }

}

void gridk_r2_bc(int i, int j, int l, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ****prt, double ****prt2, double sta, double delApol, double ct, double cta, double d1, double Drot, int nsmooth, int gbins)
{
  /*Cannot let r2 be equal to d, otherwise r1 will go to zero*/
  int k=0;
  
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  double delApol2=delApol*delApol;
  
  double c=-r*r+d1*d1;//r is r2!
  double omega=k*delaz;
  double b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
  double sfact=sqrt(b*b-4.0*c);
  double r1mag=-b/2.0+sfact/2.0;
  double r1sq=r1mag*r1mag;
  double derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
  double r2cube=r*r*r;
  double derv2=1.0/r-sfact*sfact/(4.0*r2cube);
  int knext=k+nsmooth;
  int kprev=k+nsmooth;
  double del1=nsmooth*delaz;
  
  /*az now goes from 0 to pi, so use reflecting BC.*/
  prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][knext]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0
  
  double jac=r/(r1mag+b/2.0);
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;

  int n;
  for(n=1;n<gbins-1;n++){
    k=n*nsmooth;
    if(n==gbins-2)
      knext=azbins-1;//
    else
      knext=nsmooth*(n+1);
    kprev=nsmooth*(n-1);
    
    omega=k*delaz;

  
    b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
    sfact=sqrt(b*b-4.0*c);
    r1mag=-b/2.0+sfact/2.0;
    r1sq=r1mag*r1mag;
    derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
    derv2=1.0/r-sfact*sfact/(4.0*r2cube);
    jac=r/(r1mag+b/2.0);
    prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

    //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][knext]-2.0*prt[i][j][k]+prt[i][j][kprev]) );
    psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
    passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
    //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }//done looping over az
   k=azbins-1;
   omega=k*delaz;
   kprev=(gbins-2)*nsmooth;
   double del2=(k-kprev)*delaz;
   knext=kprev;
  
   b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
   sfact=sqrt(b*b-4.0*c);
   r1mag=-b/2.0+sfact/2.0;
   r1sq=r1mag*r1mag;
   derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
   derv2=1.0/r-sfact*sfact/(4.0*r2cube);
   jac=r/(r1mag+b/2.0);
   prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del2))*(prt[i][j][l][kprev]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

   //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][kprev]) );//for i=inter, j=interior, k=end
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
  //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  /*Now interpolate*/
  double delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int ns=nsmooth;
int t;
  double dpda;
  
  for(n=1;n<gbins;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    if(n==gbins-1){
      k=azbins-1;
      kprev=(n-1)*nsmooth;
      ns=k-kprev;
      delnext=ns*delaz;
    }
	       
    dpda=(prt2[i][j][l][k]-prt2[i][j][l][kprev])/delnext;
    for(t=0;t<ns-1;t++){
      omega=(kprev+t+1)*delaz;
      b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
      sfact=sqrt(b*b-4.0*c);
      r1mag=-b/2.0+sfact/2.0;
      jac=r/(r1mag+b/2.0);
    
      prt2[i][j][l][kprev+t+1]=prt2[i][j][l][kprev]+dpda*delaz*(t+1);
      psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?

      passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
      

    }
  }

}

/*BC*/

void gridk_r2_npole2_bc(int i, int j, int l, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ****prt, double ****prt2, double sta, double delApol, double ct, double cta, double d1, double Drot, int nsmooth, int gbins)
{
  /*Cannot let r2 be equal to d, otherwise r1 will go to zero*/
  int k=0;
  
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  double delApol2=delApol*delApol;
  
  double c=-r*r+d1*d1;//r is r2!
  double omega=k*delaz;
  double b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
  double sfact=sqrt(b*b-4.0*c);
  double r1mag=-b/2.0+sfact/2.0;
  double r1sq=r1mag*r1mag;
  double derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
  double r2cube=r*r*r;
  double derv2=1.0/r-sfact*sfact/(4.0*r2cube);
  int kcross=azbins-1-k;
  int knext=k+nsmooth;
  int kprev=k+nsmooth;
  double del1=nsmooth*delaz;
  
  /*az now goes from 0 to pi, so use reflecting BC.*/
  prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][knext]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l][kcross]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) ) ;//for i=inter, j=interior, k=0
  
  double jac=r/(r1mag+b/2.0);
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;

  int n;
  for(n=1;n<gbins-1;n++){
    k=n*nsmooth;
    if(n==gbins-2)
      knext=azbins-1;//
    else
      knext=nsmooth*(n+1);
    kprev=nsmooth*(n-1);
    
    omega=k*delaz;
    kcross=(gbins-1-n)*nsmooth;
  
  
    b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
    sfact=sqrt(b*b-4.0*c);
    r1mag=-b/2.0+sfact/2.0;
    r1sq=r1mag*r1mag;
    derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
    derv2=1.0/r-sfact*sfact/(4.0*r2cube);
    jac=r/(r1mag+b/2.0);
    prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l][kcross]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) ) ;//for i=inter, j=interior, k=0

    //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][knext]-2.0*prt[i][j][k]+prt[i][j][kprev]) );
    psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
    passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
    //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }//done looping over az
   k=azbins-1;
   kcross=0;
   kprev=(gbins-2)*nsmooth;
   double del2=(k-kprev)*delaz;
   knext=kprev;
  
   omega=k*delaz;
   b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
   sfact=sqrt(b*b-4.0*c);
   r1mag=-b/2.0+sfact/2.0;
   r1sq=r1mag*r1mag;
   derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
   derv2=1.0/r-sfact*sfact/(4.0*r2cube);
   jac=r/(r1mag+b/2.0);
   prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del2)+ Drot*deltat/(sta*sta*del1*del2))*(prt[i][j][l][kprev]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l][kcross]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) ) ;//for i=inter, j=interior, k=0

   //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][kprev]) );//for i=inter, j=interior, k=end
   
   psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
   passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
   //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  /*Now interpolate*/
  double delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int ns=nsmooth;
int t;
  double dpda;
  
  for(n=1;n<gbins;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    if(n==gbins-1){
      k=azbins-1;
      kprev=(n-1)*nsmooth;
      ns=k-kprev;
      delnext=ns*delaz;
    }
	       
    dpda=(prt2[i][j][l][k]-prt2[i][j][l][kprev])/delnext;
    for(t=0;t<ns-1;t++){
      omega=(kprev+t+1)*delaz;
      b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
      sfact=sqrt(b*b-4.0*c);
      r1mag=-b/2.0+sfact/2.0;
      jac=r/(r1mag+b/2.0);
    
      prt2[i][j][l][kprev+t+1]=prt2[i][j][l][kprev]+dpda*delaz*(t+1);
      psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?

      passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;

    }
  }

}

void gridk_r2_spole2_bc(int i, int j, int l, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ****prt, double ****prt2, double sta, double delApol, double ct, double cta, double d1, double Drot, int nsmooth, int gbins)
{
  /*Cannot let r2 be equal to d, otherwise r1 will go to zero*/
  int k=0;
  
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  double delApol2=delApol*delApol;
  
  double c=-r*r+d1*d1;//r is r2!
  double omega=k*delaz;
  double b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
  double sfact=sqrt(b*b-4.0*c);
  double r1mag=-b/2.0+sfact/2.0;
  double r1sq=r1mag*r1mag;
  double derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
  double r2cube=r*r*r;
  double derv2=1.0/r-sfact*sfact/(4.0*r2cube);
  int kcross=azbins-1-k;
  int knext=k+nsmooth;
  int kprev=k+nsmooth;
  double del1=nsmooth*delaz;
  
  /*az now goes from 0 to pi, so use reflecting BC.*/
  prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][knext]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l][kcross]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0
  
  double jac=r/(r1mag+b/2.0);
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
  
  int n;
  for(n=1;n<gbins-1;n++){
    k=n*nsmooth;
    if(n==gbins-2)
      knext=azbins-1;//
    else
      knext=nsmooth*(n+1);
    kprev=nsmooth*(n-1);
    
    omega=k*delaz;
    kcross=(gbins-1-n)*nsmooth;
  

    b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
    sfact=sqrt(b*b-4.0*c);
    r1mag=-b/2.0+sfact/2.0;
    r1sq=r1mag*r1mag;
    derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
    derv2=1.0/r-sfact*sfact/(4.0*r2cube);
    jac=r/(r1mag+b/2.0);
    prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l][kcross]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

    //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][knext]-2.0*prt[i][j][k]+prt[i][j][kprev]) );
    psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
    passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
    //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }//done looping over az
   k=azbins-1;
   kcross=azbins-1-k;
   kprev=(gbins-2)*nsmooth;
   double del2=(k-kprev)*delaz;
   knext=kprev;
  
   omega=k*delaz;
   b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
   sfact=sqrt(b*b-4.0*c);
   r1mag=-b/2.0+sfact/2.0;
   r1sq=r1mag*r1mag;
   derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
   derv2=1.0/r-sfact*sfact/(4.0*r2cube);
   jac=r/(r1mag+b/2.0);
   prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del2)+ Drot*deltat/(sta*sta*del1*del2))*(prt[i][j][l][kprev]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l][kcross]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

   //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][kprev]) );//for i=inter, j=interior, k=end
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
  //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  /*Now interpolate*/
  double delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int ns=nsmooth;
int t;
  double dpda;
  
  for(n=1;n<gbins;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    if(n==gbins-1){
      k=azbins-1;
      kprev=(n-1)*nsmooth;
      ns=k-kprev;
      delnext=ns*delaz;
    }
	       
    dpda=(prt2[i][j][l][k]-prt2[i][j][l][kprev])/delnext;
    for(t=0;t<ns-1;t++){
      omega=(kprev+t+1)*delaz;
      b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
      sfact=sqrt(b*b-4.0*c);
      r1mag=-b/2.0+sfact/2.0;
      jac=r/(r1mag+b/2.0);
    
      prt2[i][j][l][kprev+t+1]=prt2[i][j][l][kprev]+dpda*delaz*(t+1);
      psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?

      passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;

    }
  }

}
void gridk_r2_npole1_bc(int i, int j, int l, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ****prt, double ****prt2, double sta, double delApol, double ct, double cta, double d1, double Drot, int nsmooth, int gbins)
{
  /*Cannot let r2 be equal to d, otherwise r1 will go to zero*/
  /*At npole, the previous bin in theta, is the same value of theta, but omega'=pi-omega.*/
  int k=0;
  
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  double delApol2=delApol*delApol;
  
  double c=-r*r+d1*d1;//r is r2!
  double omega=k*delaz;
  double b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
  double sfact=sqrt(b*b-4.0*c);
  double r1mag=-b/2.0+sfact/2.0;
  double r1sq=r1mag*r1mag;
  double derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
  double r2cube=r*r*r;
  double derv2=1.0/r-sfact*sfact/(4.0*r2cube);
  int kcross=azbins-k-1;
  int knext=k+nsmooth;
  int kprev=k+nsmooth;
  double del1=nsmooth*delaz;
  
  /*az now goes from 0 to pi, so use reflecting BC.*/
  prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j][l][kcross]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][knext]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0
  
  double jac=r/(r1mag+b/2.0);
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;

  int n;
  for(n=1;n<gbins-1;n++){
    k=n*nsmooth;
    if(n==gbins-2)
      knext=azbins-1;//
    else
      knext=nsmooth*(n+1);
    kprev=nsmooth*(n-1);
    
    omega=k*delaz;
    kcross=(gbins-1-n)*nsmooth;


    b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
    sfact=sqrt(b*b-4.0*c);
    r1mag=-b/2.0+sfact/2.0;
    r1sq=r1mag*r1mag;
    derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
    derv2=1.0/r-sfact*sfact/(4.0*r2cube);
    jac=r/(r1mag+b/2.0);
    prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j][l][kcross]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

    //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][knext]-2.0*prt[i][j][k]+prt[i][j][kprev]) );
    psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
    passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
    //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }//done looping over az
   k=azbins-1;
   kcross=0;
   kprev=(gbins-2)*nsmooth;
   double del2=(k-kprev)*delaz;
   knext=kprev;
  
   omega=k*delaz;
   b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
   sfact=sqrt(b*b-4.0*c);
   r1mag=-b/2.0+sfact/2.0;
   r1sq=r1mag*r1mag;
   derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
   derv2=1.0/r-sfact*sfact/(4.0*r2cube);
   jac=r/(r1mag+b/2.0);
   prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j][l][kcross]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del2)+ Drot*deltat/(sta*sta*del1*del2))*(prt[i][j][l][kprev]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

   //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][kprev]) );//for i=inter, j=interior, k=end
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
  //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  /*Now interpolate*/
  double delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int ns=nsmooth;
int t;
  double dpda;
  
  for(n=1;n<gbins;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    if(n==gbins-1){
      k=azbins-1;
      kprev=(n-1)*nsmooth;
      ns=k-kprev;
      delnext=ns*delaz;
    }
	       
    dpda=(prt2[i][j][l][k]-prt2[i][j][l][kprev])/delnext;
    for(t=0;t<ns-1;t++){
      omega=(kprev+t+1)*delaz;
      b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
      sfact=sqrt(b*b-4.0*c);
      r1mag=-b/2.0+sfact/2.0;
      jac=r/(r1mag+b/2.0);
    
      prt2[i][j][l][kprev+t+1]=prt2[i][j][l][kprev]+dpda*delaz*(t+1);
      psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?

      passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
      

    }
  }

}
void gridk_r2_npole1_npole2_bc(int i, int j, int l, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ****prt, double ****prt2, double sta, double delApol, double ct, double cta, double d1, double Drot, int nsmooth, int gbins)
{
  /*Cannot let r2 be equal to d, otherwise r1 will go to zero*/
  /*At npole, the previous bin in theta, is the same value of theta, but omega'=pi-omega.*/
  int k=0;
  
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  double delApol2=delApol*delApol;
  
  double c=-r*r+d1*d1;//r is r2!
  double omega=k*delaz;
  double b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
  double sfact=sqrt(b*b-4.0*c);
  double r1mag=-b/2.0+sfact/2.0;
  double r1sq=r1mag*r1mag;
  double derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
  double r2cube=r*r*r;
  double derv2=1.0/r-sfact*sfact/(4.0*r2cube);
  int kcross=azbins-k-1;
  int knext=k+nsmooth;
  int kprev=k+nsmooth;
  double del1=nsmooth*delaz;
  
  /*az now goes from 0 to pi, so use reflecting BC.*/
  prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j][l][kcross]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][knext]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l][kcross]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) ) ;//for i=inter, j=interior, k=0
  
  double jac=r/(r1mag+b/2.0);
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;

  int n;
  for(n=1;n<gbins-1;n++){
    k=n*nsmooth;
    if(n==gbins-2)
      knext=azbins-1;//
    else
      knext=nsmooth*(n+1);
    kprev=nsmooth*(n-1);
    
    omega=k*delaz;
    kcross=(gbins-1-n)*nsmooth;

    b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
    sfact=sqrt(b*b-4.0*c);
    r1mag=-b/2.0+sfact/2.0;
    r1sq=r1mag*r1mag;
    derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
    derv2=1.0/r-sfact*sfact/(4.0*r2cube);
    jac=r/(r1mag+b/2.0);
    prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j][l][kcross]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l][kcross]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) ) ;//for i=inter, j=interior, k=0

    //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][knext]-2.0*prt[i][j][k]+prt[i][j][kprev]) );
    psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
    passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
    //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }//done looping over az
   k=azbins-1;
   kcross=azbins-k-1;
   kprev=(gbins-2)*nsmooth;
   double del2=(k-kprev)*delaz;
   knext=kprev;
  
   omega=k*delaz;
   b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
   sfact=sqrt(b*b-4.0*c);
   r1mag=-b/2.0+sfact/2.0;
   r1sq=r1mag*r1mag;
   derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
   derv2=1.0/r-sfact*sfact/(4.0*r2cube);
   jac=r/(r1mag+b/2.0);
   prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j][l][kcross]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del2)+ Drot*deltat/(sta*sta*del1*del2))*(prt[i][j][l][kprev]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l][kcross]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) ) ;//for i=inter, j=interior, k=0

   //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][kprev]) );//for i=inter, j=interior, k=end
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
  //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  /*Now interpolate*/
  double delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int ns=nsmooth;
int t;
  double dpda;
  
  for(n=1;n<gbins;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    if(n==gbins-1){
      k=azbins-1;
      kprev=(n-1)*nsmooth;
      ns=k-kprev;
      delnext=ns*delaz;
    }
	       
    dpda=(prt2[i][j][l][k]-prt2[i][j][l][kprev])/delnext;
    for(t=0;t<ns-1;t++){
      omega=(kprev+t+1)*delaz;
      b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
      sfact=sqrt(b*b-4.0*c);
      r1mag=-b/2.0+sfact/2.0;
      jac=r/(r1mag+b/2.0);
    
      prt2[i][j][l][kprev+t+1]=prt2[i][j][l][kprev]+dpda*delaz*(t+1);
      psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?

      passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;

    }
  }

}
void gridk_r2_npole1_spole2_bc(int i, int j, int l, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ****prt, double ****prt2, double sta, double delApol, double ct, double cta, double d1, double Drot, int nsmooth, int gbins)
{
  /*Cannot let r2 be equal to d, otherwise r1 will go to zero*/
  /*At npole, the previous bin in theta, is the same value of theta, but omega'=pi-omega.*/
  int k=0;
  
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  double delApol2=delApol*delApol;
  
  double c=-r*r+d1*d1;//r is r2!
  double omega=k*delaz;
  double b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
  double sfact=sqrt(b*b-4.0*c);
  double r1mag=-b/2.0+sfact/2.0;
  double r1sq=r1mag*r1mag;
  double derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
  double r2cube=r*r*r;
  double derv2=1.0/r-sfact*sfact/(4.0*r2cube);
  int kcross=azbins-k-1;
  int knext=k+nsmooth;
  int kprev=k+nsmooth;
  double del1=nsmooth*delaz;
  
  /*az now goes from 0 to pi, so use reflecting BC.*/
  prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j][l][kcross]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][knext]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l][kcross]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0
  
  double jac=r/(r1mag+b/2.0);
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;

  int n;
  for(n=1;n<gbins-1;n++){
    k=n*nsmooth;
    if(n==gbins-2)
      knext=azbins-1;//
    else
      knext=nsmooth*(n+1);
    kprev=nsmooth*(n-1);
    
    omega=k*delaz;
    kcross=(gbins-1-n)*nsmooth;
  
      b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
    sfact=sqrt(b*b-4.0*c);
    r1mag=-b/2.0+sfact/2.0;
    r1sq=r1mag*r1mag;
    derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
    derv2=1.0/r-sfact*sfact/(4.0*r2cube);
    jac=r/(r1mag+b/2.0);
    prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j][l][kcross]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l][kcross]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

    //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][knext]-2.0*prt[i][j][k]+prt[i][j][kprev]) );
    psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
    passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
    //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }//done looping over az
   k=azbins-1;
   kcross=0;
   kprev=(gbins-2)*nsmooth;
   double del2=(k-kprev)*delaz;
   knext=kprev;
  
   omega=k*delaz;
   b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
   sfact=sqrt(b*b-4.0*c);
   r1mag=-b/2.0+sfact/2.0;
   r1sq=r1mag*r1mag;
   derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
   derv2=1.0/r-sfact*sfact/(4.0*r2cube);
   jac=r/(r1mag+b/2.0);
   prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j+1][l][k]-prt[i][j][l][kcross]) + 1.0/delpol2*(prt[i][j+1][l][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del2)+ Drot*deltat/(sta*sta*del1*del2))*(prt[i][j][l][kprev]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l][kcross]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

   //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][kprev]) );//for i=inter, j=interior, k=end
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
  //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  /*Now interpolate*/
  double delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int ns=nsmooth;
int t;
  double dpda;
  
  for(n=1;n<gbins;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    if(n==gbins-1){
      k=azbins-1;
      kprev=(n-1)*nsmooth;
      ns=k-kprev;
      delnext=ns*delaz;
    }
	       
    dpda=(prt2[i][j][l][k]-prt2[i][j][l][kprev])/delnext;
    for(t=0;t<ns-1;t++){
      omega=(kprev+t+1)*delaz;
      b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
      sfact=sqrt(b*b-4.0*c);
      r1mag=-b/2.0+sfact/2.0;
      jac=r/(r1mag+b/2.0);
    
      prt2[i][j][l][kprev+t+1]=prt2[i][j][l][kprev]+dpda*delaz*(t+1);
      psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?

      passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
      

    }
  }

}
void gridk_r2_spole1_bc(int i, int j, int l, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ****prt, double ****prt2, double sta, double delApol, double ct, double cta, double d1, double Drot, int nsmooth, int gbins)
{
  /*Cannot let r2 be equal to d, otherwise r1 will go to zero*/
  /*At npole, the previous bin in theta, is the same value of theta, but omega'=pi-omega.*/
  int k=0;
  
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  double delApol2=delApol*delApol;
  
  double c=-r*r+d1*d1;//r is r2!
  double omega=k*delaz;
  double b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
  double sfact=sqrt(b*b-4.0*c);
  double r1mag=-b/2.0+sfact/2.0;
  double r1sq=r1mag*r1mag;
  double derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
  double r2cube=r*r*r;
  double derv2=1.0/r-sfact*sfact/(4.0*r2cube);
  int kcross=azbins-k-1;
  int knext=k+nsmooth;
  int kprev=k+nsmooth;
  double del1=nsmooth*delaz;
  
  /*az now goes from 0 to pi, so use reflecting BC.*/
  prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j][l][kcross]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][knext]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0
  
  double jac=r/(r1mag+b/2.0);
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;

  int n;
  for(n=1;n<gbins-1;n++){
    k=n*nsmooth;
    if(n==gbins-2)
      knext=azbins-1;//
    else
      knext=nsmooth*(n+1);
    kprev=nsmooth*(n-1);
    
    omega=k*delaz;
    kcross=(gbins-1-n)*nsmooth;
  

    b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
    sfact=sqrt(b*b-4.0*c);
    r1mag=-b/2.0+sfact/2.0;
    r1sq=r1mag*r1mag;
    derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
    derv2=1.0/r-sfact*sfact/(4.0*r2cube);
    jac=r/(r1mag+b/2.0);
    prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j][l][kcross]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

    //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][knext]-2.0*prt[i][j][k]+prt[i][j][kprev]) );
    psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
    passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
    //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }//done looping over az
   k=azbins-1;
   kcross=0;
   kprev=(gbins-2)*nsmooth;
   double del2=(k-kprev)*delaz;
   knext=kprev;
  
   omega=k*delaz;
   b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
   sfact=sqrt(b*b-4.0*c);
   r1mag=-b/2.0+sfact/2.0;
   r1sq=r1mag*r1mag;
   derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
   derv2=1.0/r-sfact*sfact/(4.0*r2cube);
   jac=r/(r1mag+b/2.0);
   prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j][l][kcross]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del2)+ Drot*deltat/(sta*sta*del1*del2))*(prt[i][j][l][kprev]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

   //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][kprev]) );//for i=inter, j=interior, k=end
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
//psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  /*Now interpolate*/
  double delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int ns=nsmooth;
int t;
  double dpda;
  
  for(n=1;n<gbins;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    if(n==gbins-1){
      k=azbins-1;
      kprev=(n-1)*nsmooth;
      ns=k-kprev;
      delnext=ns*delaz;
    }
	       
    dpda=(prt2[i][j][l][k]-prt2[i][j][l][kprev])/delnext;
    for(t=0;t<ns-1;t++){
      omega=(kprev+t+1)*delaz;
      b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
      sfact=sqrt(b*b-4.0*c);
      r1mag=-b/2.0+sfact/2.0;
      jac=r/(r1mag+b/2.0);
    
      prt2[i][j][l][kprev+t+1]=prt2[i][j][l][kprev]+dpda*delaz*(t+1);
      psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?

      passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
      

    }
  }

}
void gridk_r2_spole1_npole2_bc(int i, int j, int l, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ****prt, double ****prt2, double sta, double delApol, double ct, double cta, double d1, double Drot, int nsmooth, int gbins)
{
  /*Cannot let r2 be equal to d, otherwise r1 will go to zero*/
  /*At npole, the previous bin in theta, is the same value of theta, but omega'=pi-omega.*/
  int k=0;
  
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  double delApol2=delApol*delApol;
  
  double c=-r*r+d1*d1;//r is r2!
  double omega=k*delaz;
  double b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
  double sfact=sqrt(b*b-4.0*c);
  double r1mag=-b/2.0+sfact/2.0;
  double r1sq=r1mag*r1mag;
  double derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
  double r2cube=r*r*r;
  double derv2=1.0/r-sfact*sfact/(4.0*r2cube);
  int kcross=azbins-k-1;
  int knext=k+nsmooth;
  int kprev=k+nsmooth;
  double del1=nsmooth*delaz;
  
  /*az now goes from 0 to pi, so use reflecting BC.*/
  prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j][l][kcross]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][knext]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l][kcross]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) ) ;//for i=inter, j=interior, k=0
  
  double jac=r/(r1mag+b/2.0);
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
  
  int n;
  for(n=1;n<gbins-1;n++){
    k=n*nsmooth;
    if(n==gbins-2)
      knext=azbins-1;//
    else
      knext=nsmooth*(n+1);
    kprev=nsmooth*(n-1);
    
    omega=k*delaz;
    kcross=(gbins-1-n)*nsmooth;
  
    b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
    sfact=sqrt(b*b-4.0*c);
    r1mag=-b/2.0+sfact/2.0;
    r1sq=r1mag*r1mag;
    derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
    derv2=1.0/r-sfact*sfact/(4.0*r2cube);
    jac=r/(r1mag+b/2.0);
    prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j][l][kcross]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l][kcross]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) ) ;//for i=inter, j=interior, k=0

    //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][knext]-2.0*prt[i][j][k]+prt[i][j][kprev]) );
    psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
    passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
    //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }//done looping over az
   k=azbins-1;
   kcross=azbins-k-1;
   kprev=(gbins-2)*nsmooth;
   double del2=(k-kprev)*delaz;
   knext=kprev;
  
   omega=k*delaz;
   b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
   sfact=sqrt(b*b-4.0*c);
   r1mag=-b/2.0+sfact/2.0;
   r1sq=r1mag*r1mag;
   derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
   derv2=1.0/r-sfact*sfact/(4.0*r2cube);
   jac=r/(r1mag+b/2.0);
   prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j][l][kcross]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del2)+ Drot*deltat/(sta*sta*del1*del2))*(prt[i][j][l][kprev]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l+1][k]-prt[i][j][l][kcross]) +1.0/delApol2*(prt[i][j][l+1][k]-2.0*prt[i][j][l][k]+prt[i][j][l][kcross]) ) ;//for i=inter, j=interior, k=0

   //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][kprev]) );//for i=inter, j=interior, k=end
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
//psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  /*Now interpolate*/
  double delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int ns=nsmooth;
int t;
  double dpda;
  
  for(n=1;n<gbins;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    if(n==gbins-1){
      k=azbins-1;
      kprev=(n-1)*nsmooth;
      ns=k-kprev;
      delnext=ns*delaz;
    }
	       
    dpda=(prt2[i][j][l][k]-prt2[i][j][l][kprev])/delnext;
    for(t=0;t<ns-1;t++){
      omega=(kprev+t+1)*delaz;
      b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
      sfact=sqrt(b*b-4.0*c);
      r1mag=-b/2.0+sfact/2.0;
      jac=r/(r1mag+b/2.0);
    
      prt2[i][j][l][kprev+t+1]=prt2[i][j][l][kprev]+dpda*delaz*(t+1);
      psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?

      passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
      

    }
  }

}
void gridk_r2_spole1_spole2_bc(int i, int j, int l, double r, double r2, double theta, double st, double BC, double delR, double deltat, double Dtot, double dr2, double ka, double delpol, double delaz, int azbins, double &passoc1, double &psurvive, double ****prt, double ****prt2, double sta, double delApol, double ct, double cta, double d1, double Drot, int nsmooth, int gbins)
{
  /*Cannot let r2 be equal to d, otherwise r1 will go to zero*/
  /*At npole, the previous bin in theta, is the same value of theta, but omega'=pi-omega.*/
  int k=0;
  
  double delaz2=delaz*delaz;
  double delpol2=delpol*delpol;
  double delApol2=delApol*delApol;
  
  double c=-r*r+d1*d1;//r is r2!
  double omega=k*delaz;
  double b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
  double sfact=sqrt(b*b-4.0*c);
  double r1mag=-b/2.0+sfact/2.0;
  double r1sq=r1mag*r1mag;
  double derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
  double r2cube=r*r*r;
  double derv2=1.0/r-sfact*sfact/(4.0*r2cube);
  int kcross=azbins-k-1;
  int knext=k+nsmooth;
  int kprev=k+nsmooth;
  double del1=nsmooth*delaz;
  
  /*az now goes from 0 to pi, so use reflecting BC.*/
  prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j][l][kcross]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][knext]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l][kcross]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0
  
  double jac=r/(r1mag+b/2.0);
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;

  int n;
  for(n=1;n<gbins-1;n++){
    k=n*nsmooth;
    if(n==gbins-2)
      knext=azbins-1;//
    else
      knext=nsmooth*(n+1);
    kprev=nsmooth*(n-1);
    
    omega=k*delaz;
    kcross=(gbins-1-n)*nsmooth;

    b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
    sfact=sqrt(b*b-4.0*c);
    r1mag=-b/2.0+sfact/2.0;
    r1sq=r1mag*r1mag;
    derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
    derv2=1.0/r-sfact*sfact/(4.0*r2cube);
    jac=r/(r1mag+b/2.0);
    prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j][l][kcross]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del1)+ Drot*deltat/(sta*sta*del1*del1))*(prt[i][j][l][knext]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l][kcross]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

    //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1)*(prt[i][j][knext] -2.0*prt[i][j][k]+prt[i][j][kprev]) );
    psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
    passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
    //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  }//done looping over az
   k=azbins-1;
   kcross=0;
   kprev=(gbins-2)*nsmooth;
   double del2=(k-kprev)*delaz;
   knext=kprev;
  
   omega=k*delaz;
   b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
   sfact=sqrt(b*b-4.0*c);
   r1mag=-b/2.0+sfact/2.0;
   r1sq=r1mag*r1mag;
   derv1=1.0/(2.0*r)*sfact;//(2.0*r1mag+b);=sfact  dr2/dr1
   derv2=1.0/r-sfact*sfact/(4.0*r2cube);
   jac=r/(r1mag+b/2.0);
   prt2[i][j][l][k]=prt[i][j][l][k]+Dtot*deltat/delR*( sfact/(r*r1mag) + derv2 )*(prt[i+1][j][l][k]-(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/dr2*derv1*derv1*(prt[i+1][j][l][k]-2.0*prt[i][j][l][k]+(prt[i+1][j][l][k]-BC*2.0*delR*prt[i][j][l][k]))+Dtot*deltat/r1sq*( ct/st/delpol*(prt[i][j][l][kcross]-prt[i][j-1][l][k]) + 1.0/delpol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j-1][l][k]) )+ (Dtot*deltat/r1sq*1.0/(st*st*del1*del2)+ Drot*deltat/(sta*sta*del1*del2))*(prt[i][j][l][kprev]-2.0*prt[i][j][l][k]+prt[i][j][l][kprev]) + Drot*deltat*( cta/(sta*delApol)*(prt[i][j][l][kcross]-prt[i][j][l-1][k]) +1.0/delApol2*(prt[i][j][l][kcross]-2.0*prt[i][j][l][k]+prt[i][j][l-1][k]) ) ;//for i=inter, j=interior, k=0

   //prt2[i][j][k]=prt[i][j][k]+Dtot*deltat/(r*delR)*(prt[i+1][j][k]-prt[i-1][j][k])+Dtot*deltat/dr2*(prt[i+1][j][k]-2.0*prt[i][j][k]+prt[i-1][j][k])+Dtot*deltat/r2*(cos(theta)/st/delpol*(prt[i][j+1][k]-prt[i][j-1][k]) + 1.0/delpol2*(prt[i][j+1][k]-2.0*prt[i][j][k]+prt[i][j-1][k])+1.0/(st*st*del1*del1 )*(prt[i][j][0]-2.0*prt[i][j][k]+prt[i][j][kprev]) );//for i=inter, j=interior, k=end
  psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?
  passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;
  //psurvive+=r2*delR*prt2[i][j][k]*delpol*delaz*st;
  /*Now interpolate*/
  double delnext=(nsmooth)*delaz;//back to normal grid spacing  
  int ns=nsmooth;
int t;
  double dpda;
  
  for(n=1;n<gbins;n++){
    k=n*nsmooth;
    kprev=k-nsmooth;
    if(n==gbins-1){
      k=azbins-1;
      kprev=(n-1)*nsmooth;
      ns=k-kprev;
      delnext=ns*delaz;
    }
	       
    dpda=(prt2[i][j][l][k]-prt2[i][j][l][kprev])/delnext;
    for(t=0;t<ns-1;t++){
      omega=(kprev+t+1)*delaz;
      b=2.0*d1*(st*sta*cos(omega)+ct*cta);//2.0*d1*(sin(pol)*sin(Apol)*cos(omega)+cos(pol)*cos(Apol));
      sfact=sqrt(b*b-4.0*c);
      r1mag=-b/2.0+sfact/2.0;
      jac=r/(r1mag+b/2.0);
    
      prt2[i][j][l][kprev+t+1]=prt2[i][j][l][kprev]+dpda*delaz*(t+1);
      psurvive+=r1sq*delR*prt2[i][j][l][k]*delpol*delaz*st*delApol*sta*jac;//Is extra 2pi from other azimuth needed? and extra 2 from cutting az down to pi?

      passoc1+=prt2[0][j][l][k]*ka*deltat*delpol*delaz*st*delApol*sta;


    }
  }

}


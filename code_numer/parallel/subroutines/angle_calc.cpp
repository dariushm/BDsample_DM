#include "rd_class_defs.h"
#include "angle_calc.h"
#include "vector_rot_calls.h"

using namespace std;

double calc_psi(Fullmol *bases)
{
  /*the vectors that are read in are just place holders, they are overwritten
    by molecule geometry.
  */
  double dxm, dym, dzm;
  
  
  /*calculate molecule dihedral*/
  dxm=bases[1].xcom-bases[0].xcom;
  dym=bases[1].ycom-bases[0].ycom;
  dzm=bases[1].zcom-bases[0].zcom;
  double R2=dxm*dxm+dym*dym+dzm*dzm;
  double R1=sqrt(R2);
  
  
  
  //calculate dihedral angle.
  //n1=(1-0)com (x) (r1)
  double *c0toc1=new double[3];
  double *c1toim=new double[3];
  double *c0toim=new double[3];
  double *n1=new double[3];
  double *n2=new double[3];
  c0toc1[0]=dxm/R1;
  c0toc1[1]=dym/R1;
  c0toc1[2]=dzm/R1;
  
  c0toim[0]=-(bases[0].x[1]-bases[0].xcom);//should be length 1
  c0toim[1]=-(bases[0].y[1]-bases[0].ycom);
  c0toim[2]=-(bases[0].z[1]-bases[0].zcom);
  
  c1toim[0]=bases[1].x[1]-bases[1].xcom;//should be length 1
  c1toim[1]=bases[1].y[1]-bases[1].ycom;
  c1toim[2]=bases[1].z[1]-bases[1].zcom;
  
    
  double dp=c0toc1[0]*c1toim[0]+c0toc1[1]*c1toim[1]+c0toc1[2]*c1toim[2];
  double psimol;
  if(dp==1 || dp==-1){
    psimol=0;//by convention, if paralell set angle to zero.
    cout <<"parallel orient and com connector, set psi zero " <<endl;
  }else{
    crossproduct(c0toc1, c1toim, n1);
    //n2=(-r0) (x) (1-0)com
    dp=c0toc1[0]*c0toim[0]+c0toc1[1]*c0toim[1]+c0toc1[2]*c0toim[2];

    if(dp==1 || dp==-1){
      psimol=0;
      cout <<"parallel orient and com connector, set psi zero " <<endl;
    }else{
      crossproduct(c0toim, c0toc1, n2);
      double cospsi=n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
      psimol=acos(cospsi);//if this is close to zero, molecules are oriented same
      if(isnan(psimol)){
	if(round(cospsi)==-1)psimol=M_PI;
	else psimol=0;
	cout <<"NAN: "<<psimol<<" cospsi: "<<cospsi<<endl;
      }  
      
    }
  }
  delete[] c0toc1;
  delete[] c1toim;
  delete[] c0toim;
  delete[] n1;
  delete[] n2;
  return psimol;
}
double calc_psi_2pi(Fullmol *bases, double leglen, int it)
{
  /*the vectors that are read in are just place holders, they are overwritten
    by molecule geometry.
  */
  double dxm, dym, dzm;

  
  /*calculate molecule dihedral*/
  dxm=bases[1].xcom-bases[0].xcom;
  dym=bases[1].ycom-bases[0].ycom;
  dzm=bases[1].zcom-bases[0].zcom;
  double R2=dxm*dxm+dym*dym+dzm*dzm;
  double R1=sqrt(R2);

  double *c0toc1=new double[3];
  double *c1toim=new double[3];
  double *c0toim=new double[3];
  double *n1=new double[3];
  double *n2=new double[3];
  
  c0toc1[0]=dxm/R1;
  c0toc1[1]=dym/R1;
  c0toc1[2]=dzm/R1;
  
  c0toim[0]=-(bases[0].x[1]-bases[0].xcom);//should be length 1
  c0toim[1]=-(bases[0].y[1]-bases[0].ycom);
  c0toim[2]=-(bases[0].z[1]-bases[0].zcom);
  
  c1toim[0]=bases[1].x[1]-bases[1].xcom;//should be length 1
  c1toim[1]=bases[1].y[1]-bases[1].ycom;
  c1toim[2]=bases[1].z[1]-bases[1].zcom;
  //calculate dihedral angle.
  //n1=(1-0)com (x) (r1)
  double *M=new double[9];
  double *rotv=new double[3];
  double tol=1E-7;
  double dp1=c0toc1[0]*c1toim[0]+c0toc1[1]*c1toim[1]+c0toc1[2]*c1toim[2];
  double dp0=c0toc1[0]*c0toim[0]+c0toc1[1]*c0toim[1]+c0toc1[2]*c0toim[2];
  double psimol;
  double angle=0.01;
  double *c1toib=new double[3];
  double *c0toib=new double[3];
  
  if(1-abs(dp1)<tol){
    /*This occurs when the binding interface is perp to the c0toc1, then
      one of the im vectors will be parallel to c0toc1.
      Can define the psi value based on rotations about, since psi is independent
      of rotations about this binding interface.
      Will need to choose sign of psi, which switches by pi.
    */
    /*need to calc c1toib*/
    //cout<<" ctoim parallel to comvec, mol1: "<<dp1<<" iter: "<<it<<'\t';
    c1toib[0]=(bases[1].x[0]-bases[1].xcom)/leglen;
    c1toib[1]=(bases[1].y[0]-bases[1].ycom)/leglen;
    c1toib[2]=(bases[1].z[0]-bases[1].zcom)/leglen;
     // cout <<"com crds: "<<bases[1].xcom<<' '<<bases[1].ycom<<' '<<bases[1].zcom<<endl;
//      cout <<"ib crds: "<<bases[1].x[0]<<' '<<bases[1].y[0]<<' '<<bases[1].z[0]<<endl;
//      cout <<"im crds: "<<bases[1].x[1]<<' '<<bases[1].y[1]<<' '<<bases[1].z[1]<<endl;
    calc_Rmatrix(c1toib, angle, M);
    rotate(c1toim, M, rotv);
    /*calculate psi using vector rotv.*/
    crossproduct(c0toc1, rotv, n1);
    if(1-abs(dp0)<tol){
      /*need to calc c0toib*/
      // cout<<" ctoim parallel to comvec, mol0: "<<dp0<<" iter: "<<it<<'\t';
//        cout <<"com crds: "<<bases[0].xcom<<' '<<bases[0].ycom<<' '<<bases[0].zcom<<endl;
//        cout <<"ib crds: "<<bases[0].x[0]<<' '<<bases[0].y[0]<<' '<<bases[0].z[0]<<endl;
//        cout <<"im crds: "<<bases[0].x[1]<<' '<<bases[0].y[1]<<' '<<bases[0].z[1]<<endl;
    
      c0toib[0]=(bases[0].x[0]-bases[0].xcom)/leglen;
      c0toib[1]=(bases[0].y[0]-bases[0].ycom)/leglen;
      c0toib[2]=(bases[0].z[0]-bases[0].zcom)/leglen;
    
      calc_Rmatrix(c0toib, angle, M);
      rotate(c0toim, M, rotv);
      /*calculate psi using vector rotv.*/
      crossproduct(rotv, c0toc1, n2);
    }else{
      crossproduct(c0toim,  c0toc1, n2);//returns normalized vectors
    }
      
    
  }else{
    crossproduct(c0toc1, c1toim, n1);
    //n2=(-r0) (x) (1-0)com
    
    if(1-abs(dp0)<tol){
      // cout<<" ctoim parallel to comvec, mol0: "<<dp0<<" iter: "<<it<<'\t';
//      cout <<"com crds: "<<bases[0].xcom<<' '<<bases[0].ycom<<' '<<bases[0].zcom<<endl;
//      cout <<"ib crds: "<<bases[0].x[0]<<' '<<bases[0].y[0]<<' '<<bases[0].z[0]<<endl;
//      cout <<"im crds: "<<bases[0].x[1]<<' '<<bases[0].y[1]<<' '<<bases[0].z[1]<<endl;
     
      c0toib[0]=(bases[0].x[0]-bases[0].xcom)/leglen;
      c0toib[1]=(bases[0].y[0]-bases[0].ycom)/leglen;
      c0toib[2]=(bases[0].z[0]-bases[0].zcom)/leglen;
    
      calc_Rmatrix(c0toib, angle, M);
      rotate(c0toim, M, rotv);
      /*calculate psi using vector rotv.*/
      crossproduct(rotv, c0toc1, n2);
    }else{
      crossproduct(c0toim,  c0toc1, n2);//returns normalized vectors
    }
  }
  double cospsi=n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
  psimol=acos(cospsi);//if this is close to zero, molecules are oriented same
  if(isnan(psimol)){
    if(round(cospsi)==-1)psimol=M_PI;
    else psimol=0;
    //cout <<"NAN: "<<psimol<<" cospsi: "<<cospsi<<endl;
  }  
  
  /*Now we want to differentiate between +/-psi
   */
  calc_Rmatrix(c0toc1, psimol, M);
  rotate(n1, M, rotv);//if rotating the same direction as initial rotation aligns, then +.  
    
  double dp=rotv[0]*n2[0]+rotv[1]*n2[1]+rotv[2]*n2[2];
  tol=1E-6;
  if(abs(dp-1)>tol){
    //cout <<"minus.orig psi: "<<psimol<<'\t';
    psimol=2*M_PI-psimol;
    //cout <<" angle: "<<acos(dp)<<" new psi: "<<psimol<<endl;
  }
  delete[] c0toc1;
  delete[] c1toim;
  delete[] c0toim;
  delete[] c1toib;
  delete[] c0toib;
  delete[] n1;
  delete[] n2;
  delete[] rotv;
  delete[] M;

  return psimol;
}
double calc_psi_2pi(Fullmol *bases, double leglen1, double leglen0)
{
  /*the vectors that are read in are just place holders, they are overwritten
    by molecule geometry.
  */
  double dxm, dym, dzm;

  
  /*calculate molecule dihedral*/
  dxm=bases[1].xcom-bases[0].xcom;
  dym=bases[1].ycom-bases[0].ycom;
  dzm=bases[1].zcom-bases[0].zcom;
  double R2=dxm*dxm+dym*dym+dzm*dzm;
  double R1=sqrt(R2);

  double *c0toc1=new double[3];
  double *c1toim=new double[3];
  double *c0toim=new double[3];
  double *n1=new double[3];
  double *n2=new double[3];
  
  c0toc1[0]=dxm/R1;
  c0toc1[1]=dym/R1;
  c0toc1[2]=dzm/R1;
  
  c0toim[0]=-(bases[0].x[1]-bases[0].xcom);//should be length 1
  c0toim[1]=-(bases[0].y[1]-bases[0].ycom);
  c0toim[2]=-(bases[0].z[1]-bases[0].zcom);
  
  c1toim[0]=bases[1].x[1]-bases[1].xcom;//should be length 1
  c1toim[1]=bases[1].y[1]-bases[1].ycom;
  c1toim[2]=bases[1].z[1]-bases[1].zcom;
  //calculate dihedral angle.
  //n1=(1-0)com (x) (r1)
  double *M=new double[9];
  double *rotv=new double[3];
  double tol=1E-7;
  double dp1=c0toc1[0]*c1toim[0]+c0toc1[1]*c1toim[1]+c0toc1[2]*c1toim[2];
  double dp0=c0toc1[0]*c0toim[0]+c0toc1[1]*c0toim[1]+c0toc1[2]*c0toim[2];
  double psimol;
  double angle=0.01;
  double *c1toib=new double[3];
  double *c0toib=new double[3];
  
  if(1-abs(dp1)<tol){
    /*This occurs when the binding interface is perp to the c0toc1, then
      one of the im vectors will be parallel to c0toc1.
      Can define the psi value based on rotations about, since psi is independent
      of rotations about this binding interface.
      Will need to choose sign of psi, which switches by pi.
    */
    /*need to calc c1toib*/
    //    cout<<" ctoim parallel to comvec, mol1: "<<dp1<<endl;
    c1toib[0]=(bases[1].x[0]-bases[1].xcom)/leglen1;
    c1toib[1]=(bases[1].y[0]-bases[1].ycom)/leglen1;
    c1toib[2]=(bases[1].z[0]-bases[1].zcom)/leglen1;
    calc_Rmatrix(c1toib, angle, M);
    rotate(c1toim, M, rotv);
    /*calculate psi using vector rotv.*/
    crossproduct(c0toc1, rotv, n1);
    if(1-abs(dp0)<tol){
      /*need to calc c0toib*/
      //cout<<" ctoim parallel to comvec, mol0: "<<dp0<<endl;
      
      c0toib[0]=(bases[0].x[0]-bases[0].xcom)/leglen0;
      c0toib[1]=(bases[0].y[0]-bases[0].ycom)/leglen0;
      c0toib[2]=(bases[0].z[0]-bases[0].zcom)/leglen0;
    
      calc_Rmatrix(c0toib, angle, M);
      rotate(c0toim, M, rotv);
      /*calculate psi using vector rotv.*/
      crossproduct(rotv, c0toc1, n2);
    }else{
      crossproduct(c0toim,  c0toc1, n2);//returns normalized vectors
    }
      
    
  }else{
    crossproduct(c0toc1, c1toim, n1);
    //n2=(-r0) (x) (1-0)com
    
    if(1-abs(dp0)<tol){
      //cout<<" ctoim parallel to comvec, mol0: "<<dp0<<endl;
      c0toib[0]=(bases[0].x[0]-bases[0].xcom)/leglen0;
      c0toib[1]=(bases[0].y[0]-bases[0].ycom)/leglen0;
      c0toib[2]=(bases[0].z[0]-bases[0].zcom)/leglen0;
    
      calc_Rmatrix(c0toib, angle, M);
      rotate(c0toim, M, rotv);
      /*calculate psi using vector rotv.*/
      crossproduct(rotv, c0toc1, n2);
    }else{
      crossproduct(c0toim,  c0toc1, n2);//returns normalized vectors
    }
  }
  double cospsi=n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
  psimol=acos(cospsi);//if this is close to zero, molecules are oriented same
  if(isnan(psimol)){
    if(round(cospsi)==-1)psimol=M_PI;
    else psimol=0;
    //cout <<"NAN: "<<psimol<<" cospsi: "<<cospsi<<endl;
  }  
  
  /*Now we want to differentiate between +/-psi
   */
  calc_Rmatrix(c0toc1, psimol, M);
  rotate(n1, M, rotv);//if rotating the same direction as initial rotation aligns, then +.  
    
  double dp=rotv[0]*n2[0]+rotv[1]*n2[1]+rotv[2]*n2[2];
  tol=1E-6;
  if(abs(dp-1)>tol){
    //cout <<"minus.orig psi: "<<psimol<<'\t';
    psimol=2*M_PI-psimol;
    //cout <<" angle: "<<acos(dp)<<" new psi: "<<psimol<<endl;
  }
  delete[] c0toc1;
  delete[] c1toim;
  delete[] c0toim;
  delete[] c1toib;
  delete[] c0toib;
  delete[] n1;
  delete[] n2;
  delete[] rotv;
  delete[] M;

  return psimol;
}
double calc_psi_1pi(Fullmol *bases, double leglen1, double leglen0)
{
  /*the vectors that are read in are just place holders, they are overwritten
    by molecule geometry.
  */
  double dxm, dym, dzm;

  
  /*calculate molecule dihedral*/
  dxm=bases[1].xcom-bases[0].xcom;
  dym=bases[1].ycom-bases[0].ycom;
  dzm=bases[1].zcom-bases[0].zcom;
  double R2=dxm*dxm+dym*dym+dzm*dzm;
  double R1=sqrt(R2);

  double *c0toc1=new double[3];
  double *c1toim=new double[3];
  double *c0toim=new double[3];
  double *n1=new double[3];
  double *n2=new double[3];
  
  c0toc1[0]=dxm/R1;
  c0toc1[1]=dym/R1;
  c0toc1[2]=dzm/R1;
  
  c0toim[0]=-(bases[0].x[1]-bases[0].xcom);//should be length 1
  c0toim[1]=-(bases[0].y[1]-bases[0].ycom);
  c0toim[2]=-(bases[0].z[1]-bases[0].zcom);
  
  c1toim[0]=bases[1].x[1]-bases[1].xcom;//should be length 1
  c1toim[1]=bases[1].y[1]-bases[1].ycom;
  c1toim[2]=bases[1].z[1]-bases[1].zcom;
  //calculate dihedral angle.
  //n1=(1-0)com (x) (r1)
  double *M=new double[9];
  double *rotv=new double[3];
  double tol=1E-7;
  double dp1=c0toc1[0]*c1toim[0]+c0toc1[1]*c1toim[1]+c0toc1[2]*c1toim[2];
  double dp0=c0toc1[0]*c0toim[0]+c0toc1[1]*c0toim[1]+c0toc1[2]*c0toim[2];
  double psimol;
  double angle=0.01;
  double *c1toib=new double[3];
  double *c0toib=new double[3];
  
  if(1-abs(dp1)<tol){
    /*This occurs when the binding interface is perp to the c0toc1, then
      one of the im vectors will be parallel to c0toc1.
      Can define the psi value based on rotations about, since psi is independent
      of rotations about this binding interface.
      Will need to choose sign of psi, which switches by pi.
    */
    /*need to calc c1toib*/
    //    cout<<" ctoim parallel to comvec, mol1: "<<dp1<<endl;
    c1toib[0]=(bases[1].x[0]-bases[1].xcom)/leglen1;
    c1toib[1]=(bases[1].y[0]-bases[1].ycom)/leglen1;
    c1toib[2]=(bases[1].z[0]-bases[1].zcom)/leglen1;
    calc_Rmatrix(c1toib, angle, M);
    rotate(c1toim, M, rotv);
    /*calculate psi using vector rotv.*/
    crossproduct(c0toc1, rotv, n1);
    if(1-abs(dp0)<tol){
      /*need to calc c0toib*/
      //cout<<" ctoim parallel to comvec, mol0: "<<dp0<<endl;
      
      c0toib[0]=(bases[0].x[0]-bases[0].xcom)/leglen0;
      c0toib[1]=(bases[0].y[0]-bases[0].ycom)/leglen0;
      c0toib[2]=(bases[0].z[0]-bases[0].zcom)/leglen0;
    
      calc_Rmatrix(c0toib, angle, M);
      rotate(c0toim, M, rotv);
      /*calculate psi using vector rotv.*/
      crossproduct(rotv, c0toc1, n2);
    }else{
      crossproduct(c0toim,  c0toc1, n2);//returns normalized vectors
    }
      
    
  }else{
    crossproduct(c0toc1, c1toim, n1);
    //n2=(-r0) (x) (1-0)com
    
    if(1-abs(dp0)<tol){
      //cout<<" ctoim parallel to comvec, mol0: "<<dp0<<endl;
      c0toib[0]=(bases[0].x[0]-bases[0].xcom)/leglen0;
      c0toib[1]=(bases[0].y[0]-bases[0].ycom)/leglen0;
      c0toib[2]=(bases[0].z[0]-bases[0].zcom)/leglen0;
    
      calc_Rmatrix(c0toib, angle, M);
      rotate(c0toim, M, rotv);
      /*calculate psi using vector rotv.*/
      crossproduct(rotv, c0toc1, n2);
    }else{
      crossproduct(c0toim,  c0toc1, n2);//returns normalized vectors
    }
  }
  double cospsi=n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
  psimol=acos(cospsi);//if this is close to zero, molecules are oriented same
  if(isnan(psimol)){
    if(round(cospsi)==-1)psimol=M_PI;
    else psimol=0;
    //cout <<"NAN: "<<psimol<<" cospsi: "<<cospsi<<endl;
  }  
  
  delete[] c0toc1;
  delete[] c1toim;
  delete[] c0toim;
  delete[] c1toib;
  delete[] c0toib;
  delete[] n1;
  delete[] n2;
  delete[] rotv;
  delete[] M;

  return psimol;
}
void calc_three_angle(double &cthet1, double &cthet2, double &dih, Fullmol *bases, double dx, double dy, double dz, double R1, double leglen1, double leglen0, double *v, double *v1, double *n1, double *n2)
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
  cthet1/=(R1*leglen1);
  
  cthet2=r0x*dx+r0y*dy+r0z*dz;
  cthet2/=(R1*leglen0);
  
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
  if(isnan(dih)){
    /*this can happen if there is numerical noise such that cos(dih) is 1 or -1 but a little beyond*/
    if(round(cosdih)==-1)dih=M_PI;
    else dih=0;
  }
  /*in the cases below one vector is oriented along separating vector, such that dih is not defined.
    choose to set it to PI, such that they are assumed oriented.
  */
  if(cthet1==1 || cthet1==-1)
    dih=M_PI;
  else if(cthet2==1 || cthet2==-1)
    dih=M_PI;
 
 
}
double numer_calc_Jacobian(Fullmol *bases, double psi0,  double leglen1, double leglen0, int molrot)
{

  /*Given the current position of the molecule orientation vector, 
    rotate molecule 0 around it's interface vector, to change phi2.
    measure the corersponding change in psi (the dihedral), for the current
    value of psi and phi1.  (dpsi/dphi2 as a functino of psi, and phi1).
  */
  double *c0toim=new double[3];
  double *c0toib=new double[3];
  double *n1=new double[3];
  double *n2=new double[3];
  double *M=new double[9];
  c0toim[0]=bases[molrot].x[1]-bases[molrot].xcom;
  c0toim[1]=bases[molrot].y[1]-bases[molrot].ycom;
  c0toim[2]=bases[molrot].z[1]-bases[molrot].zcom;
  double leglenr;
  if(molrot==0)leglenr=leglen0;
  else leglenr=leglen1;
  c0toib[0]=(bases[molrot].x[0]-bases[molrot].xcom)/leglenr;
  c0toib[1]=(bases[molrot].y[0]-bases[molrot].ycom)/leglenr;
  c0toib[2]=(bases[molrot].z[0]-bases[molrot].zcom)/leglenr;
  //cout <<"ROTATED DIHEDRAL AGAIN, previous: "<<bases[1].x[1]<<' '<<bases[1].y[1]<<' '<<bases[1].z[1]<<" COM: "<<bases[1].xcom<<' '<<bases[1].ycom<<' '<<bases[1].zcom<<endl;
  double dphi2=1E-3;//small displacement
  calc_Rmatrix(c0toib, dphi2, M);
  rotate(c0toim, M, n1);
  
  calc_Rmatrix(c0toib, -dphi2, M);
  rotate(c0toim, M, n2);
  
  bases[molrot].x[1]=n1[0]+bases[molrot].xcom;
  bases[molrot].y[1]=n1[1]+bases[molrot].ycom;
  bases[molrot].z[1]=n1[2]+bases[molrot].zcom;
  //=abs(psimol_plus-psi0)/dphi2;
  //=abs(psimol_minus-psi0)/dphi2;
  
  double psimol_plus=calc_psi_2pi(bases, leglen1, leglen0);
  double dfup=abs(psi0-psimol_plus);
  if(dfup>M_PI){
    //cout <<"psimol+dphi rotated: "<<psimol_plus<<" corrected: "<<2*M_PI+psimol_plus<<endl;
    //    psimol_plus=2*M_PI+psimol_plus;
    dfup=2*M_PI-dfup;
  }
  bases[molrot].x[1]=n2[0]+bases[molrot].xcom;
  bases[molrot].y[1]=n2[1]+bases[molrot].ycom;
  bases[molrot].z[1]=n2[2]+bases[molrot].zcom;
  double psimol_minus=calc_psi_2pi(bases, leglen1, leglen0);
  double dfdown=abs(psimol_minus-psi0);
  if(dfdown>M_PI){
    //cout <<"psimol-dphi rotated: "<<psimol_minus<<" corrected: "<<psimol_minus-2*M_PI<<endl;
    //psimol_minus=psimol_minus-2*M_PI;//negative number.
    dfdown=2*M_PI-dfdown;
  }
  double tol=dphi2*100;
  
  if(dfup/dphi2-dfdown/dphi2>tol)cout <<"CAREFUL, dphi/dpsi has different values approaching up/down "<<dfup/dphi2<<'\t'<<dfdown/dphi2<<endl; 
  double derv=(dfup+dfdown)/(2*dphi2);
  /*return orientation back to start*/
  bases[molrot].x[1]=c0toim[0]+bases[molrot].xcom;
  bases[molrot].y[1]=c0toim[1]+bases[molrot].ycom;
  bases[molrot].z[1]=c0toim[2]+bases[molrot].zcom;
  /*the factor we need is |dphi2/dpsi|, or derivative inverse. */
  return 1.0/derv;
}


double set_norm_to_psi(double psi0, Fullmol *bases, double leglen1, double leglen0, int molfix, int molrot, double *csettoim)
{
  //COM to COM vector
  double dxm=bases[molfix].xcom-bases[molrot].xcom;
  double dym=bases[molfix].ycom-bases[molrot].ycom;
  double dzm=bases[molfix].zcom-bases[molrot].zcom;
  double R2=dxm*dxm+dym*dym+dzm*dzm;
  double R1=sqrt(R2);
  double tol=1E-7;//to account for variations around exact matching psi
  //rotate norm around this vector by psi.
  double *comv=new double[3];
  double *cfixtoim=new double[3];
  comv[0]=dxm/R1;
  comv[1]=dym/R1;
  comv[2]=dzm/R1;
  cfixtoim[0]=bases[molfix].x[1]-bases[molfix].xcom;//should be length 1.
  cfixtoim[1]=bases[molfix].y[1]-bases[molfix].ycom;
  cfixtoim[2]=bases[molfix].z[1]-bases[molfix].zcom;
  double *csettoib=new double[3];
  double *cfixtoib=new double[3];
  double *M=new double[9];
  double *plane_n1=new double[3];
  double *rotv=new double[3];
  double *v3=new double[3];
  double angle=0.01;
  double dp1=comv[0]*cfixtoim[0]+comv[1]*cfixtoim[1]+comv[2]*cfixtoim[2];
  double leglenf, leglenr;
  if(molfix==0){
    leglenf=leglen0;
    leglenr=leglen1;
  }else{
    leglenf=leglen1;
    leglenr=leglen0;
  }
  
  if(1-abs(dp1)<tol){
    //cout <<"mol orient points along COM vec! "<<endl;
    /*normal only points either towards ib or away from it */
    
    cfixtoib[0]=(bases[molfix].x[0]-bases[molfix].xcom)/leglenf;
    cfixtoib[1]=(bases[molfix].y[0]-bases[molfix].ycom)/leglenf;
    cfixtoib[2]=(bases[molfix].z[0]-bases[molfix].zcom)/leglenf;
  
    calc_Rmatrix(cfixtoib, angle, M);
    rotate(cfixtoim, M, v3);//temporary, to determine a normal;
    //instead of using cfix to im, which won't rotate around a parallel vec, use an offset.
    calc_Rmatrix(comv, psi0, M);
    rotate(v3, M, rotv);//now n1 points in the direction of the normal for leg0 plane, fix to be perp to leg.
    
    //then use that vector to set target plane via normal
    crossproduct(comv, rotv, plane_n1);//plane_n1 is normal to plane need to intersect.
    //cout <<"direction of normal to mol orient plane: "<<plane_n1[0]<<'\t'<<plane_n1[1]<<'\t'<<plane_n1[2]<<endl;
  
  }else{
    calc_Rmatrix(comv, psi0, M);
    rotate(cfixtoim, M, rotv);//now n1 points in the direction of the normal for leg0 plane, fix to be perp to leg.
    
    //then use that vector to set target plane via normal
    crossproduct(comv, rotv, plane_n1);//plane_n1 is normal to plane need to intersect.
    //cout <<"direction of normal to mol orient plane: "<<plane_n1[0]<<'\t'<<plane_n1[1]<<'\t'<<plane_n1[2]<<endl;
  }
  //  cout <<"direction of normal to plane; "<<plane_n1[0]<<' '<<plane_n1[1]<<' '<<plane_n1[2]<<endl;
  //move central point for that plane to bases[0].com
  double h1=plane_n1[0]*bases[molrot].xcom+plane_n1[1]*bases[molrot].ycom+plane_n1[2]*bases[molrot].zcom;
  
  //plane 2 normal is leg_0 vector
  csettoib[0]=(bases[molrot].x[0]-bases[molrot].xcom)/leglenr;
  csettoib[1]=(bases[molrot].y[0]-bases[molrot].ycom)/leglenr;
  csettoib[2]=(bases[molrot].z[0]-bases[molrot].zcom)/leglenr;

  //csettoib dot r0
  double h0=csettoib[0]*bases[molrot].xcom+csettoib[1]*bases[molrot].ycom+csettoib[2]*bases[molrot].zcom;
  double n1n2=plane_n1[0]*csettoib[0]+plane_n1[1]*csettoib[1]+plane_n1[2]*csettoib[2];
  double p1x, p1y, p1z;
  double p1len, psi1, psimol;
  /*parallel check should be done prior to setting dihedral. If both c0toib vectors
    are perpendicular to c0toc1, then they represent normals to the mol orientations for all rotations,
    meaning psi only takes two values, separated by pi.
  */
  
  double c1cof=(h1-h0*n1n2)/(1-n1n2*n1n2);
  double c0cof=(h0-h1*n1n2)/(1-n1n2*n1n2);
  crossproduct(plane_n1, csettoib, v3);//v3 is vector in both planes (normal to both normals.)
  csettoim[0]=c1cof*plane_n1[0]+c0cof*csettoib[0]+v3[0]-bases[molrot].xcom;
  csettoim[1]=c1cof*plane_n1[1]+c0cof*csettoib[1]+v3[1]-bases[molrot].ycom;
  csettoim[2]=c1cof*plane_n1[2]+c0cof*csettoib[2]+v3[2]-bases[molrot].zcom;
  /*this vector is through a line of points intersecting the two planes,
    we choose an abritrary point on the line so the dihedral angle
    could potentially end up rotated by pi, hence the tolerance
    criterion.
  */
  R2=csettoim[0]*csettoim[0]+csettoim[1]*csettoim[1]+csettoim[2]*csettoim[2];
  R1=sqrt(R2);
  csettoim[0]/=R1;//we'll return this vector
  csettoim[1]/=R1;
  csettoim[2]/=R1;
  
  bases[molrot].x[1]=bases[molrot].xcom+csettoim[0];
  bases[molrot].y[1]=bases[molrot].ycom+csettoim[1];
  bases[molrot].z[1]=bases[molrot].zcom+csettoim[2];
  psimol=calc_psi_2pi(bases, leglen1, leglen0);
  
  if(abs(psi0-psimol)>tol){
    psi1=psimol;
    csettoim[0]*=-1;
    csettoim[1]*=-1;
    csettoim[2]*=-1;
    bases[molrot].x[1]=bases[molrot].xcom+csettoim[0];
    bases[molrot].y[1]=bases[molrot].ycom+csettoim[1];
    bases[molrot].z[1]=bases[molrot].zcom+csettoim[2];
    
    psimol=calc_psi_2pi(bases, leglen1, leglen0);
    //cout <<"switch direction! original "<<psi1<<" new psi: "<<psimol<<" target "<<psi0<<endl;
  }
  

  delete[] plane_n1;
  delete[] rotv;
  delete[] comv;

  delete[] v3;
  delete[] csettoib;
  delete[] M;
  delete[] cfixtoim;
  return psimol;

}

double rotate_phi2_mol(double psi0, Fullmol *bases, double leglen1, double leglen0, int molfix, int molrot, double *csettoim)
{
  double psimol;
  double *csettoib=new double[3];
  double *rotv=new double[3];
  double *comv=new double[3];
  double *M=new double[9];
  double leglenr;
  if(molrot==0)leglenr=leglen0;
  else leglenr=leglen1;
  
  
  csettoib[0]=(bases[molrot].x[0]-bases[molrot].xcom)/leglenr;
  csettoib[1]=(bases[molrot].y[0]-bases[molrot].ycom)/leglenr;
  csettoib[2]=(bases[molrot].z[0]-bases[molrot].zcom)/leglenr;
  double dxm=bases[molfix].xcom-bases[molrot].xcom;
  double dym=bases[molfix].ycom-bases[molrot].ycom;
  double dzm=bases[molfix].zcom-bases[molrot].zcom;
  double R2=dxm*dxm+dym*dym+dzm*dzm;
  double R1=sqrt(R2);
  comv[0]=dxm/R1;
  comv[1]=dym/R1;
  comv[2]=dzm/R1;
  cout <<"In rotate phi2, phi2: "<<psi0<<" molfix: "<<molfix<<" fixed im pos: "<<bases[molfix].x[1]<<'\t'<<bases[molfix].y[1]<<'\t'<<bases[molfix].z[1]<<'\t';
  calc_Rmatrix(csettoib, psi0, M);
  
  rotate(comv, M, rotv);//comv is perpendicular to cset toib, so it represents one mol orientation direction
  R2=rotv[0]*rotv[0]+rotv[1]*rotv[1]+rotv[2]*rotv[2];
  R1=sqrt(R2);
  csettoim[0]=rotv[0]/R1;
  csettoim[1]=rotv[1]/R1;
  csettoim[2]=rotv[2]/R1;
  
  bases[molrot].x[1]=bases[molrot].xcom+csettoim[0];
  bases[molrot].y[1]=bases[molrot].ycom+csettoim[1];
  bases[molrot].z[1]=bases[molrot].zcom+csettoim[2];
  cout <<" rotated im pos: "<<bases[molrot].x[1]<<'\t'<<bases[molrot].y[1]<<'\t'<<bases[molrot].z[1]<<'\t';
  psimol=calc_psi_2pi(bases, leglen1, leglen0);
  
  
  
   delete[] csettoib;
   delete[] rotv;
   delete[] comv;
   delete[] M;
  
  return psimol;

}

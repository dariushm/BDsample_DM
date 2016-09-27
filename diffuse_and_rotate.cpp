/*
generate coordinates for clathrin molecules, then propagate them
undergoing rotational and translational diffusion, without any boundaries,
without any reactions. So it is just moving them to measure the total squared
displacement with time, as a test of Dtrans and Drot sampling. 
2013
*/

#include <fstream>
#include <iostream>
#include <cmath>
#include "rand_gsl.h"
#include "vector_rot_calls.h"
#include "utility_calls.h"

using namespace std;

int main(int argc, char *argv[])
{


  int i, j, k;
  int jnext, t2, t;
  
  timeval tim;
  gettimeofday(&tim, 0);
  double t1=tim.tv_sec+tim.tv_usec;
  int seed=int(t1);
  
  //seed=1364043256;
  cout <<"seed: "<<seed<<endl;
  srand_gsl(seed);
  

  int N; //number of particles
  int Tsteps; //number of time steps you want to measure Dt over
  
  int Ncycles; //Should be number of configurations
  int Nsteps; //Number of timesteps between config saves
  double dx, dy, dz;
  int timezero;

  double deltat;
  double delT;
  double leglen;
  double thetaleg;
  double Dtrans;
  double Drot;
  /*
    Examples:
    
    N=100 molecules
    Nsteps=1000 steps
    deltat=1 us
    leglen=15 nm
    Dtrans=10 nm^2/us
    Drot=0.02 rad^2/us

  */
  
  N=atoi(argv[1]);//total molecules
  Nsteps=atoi(argv[2]);//total steps
  deltat=atof(argv[3]);//units us
  leglen=atof(argv[4]);//units nm
  Dtrans=atof(argv[5]);//units nm^2/us
  Drot=atof(argv[6]);//units rad^2/us
  
  thetaleg=120.0*M_PI/180.0;//three legs, 120 degrees apart
  
  double **xcm=new double*[Nsteps+1];
  double **ycm=new double*[Nsteps+1];
  double **zcm=new double*[Nsteps+1];
  double **xleg1=new double*[Nsteps+1];
  double **yleg1=new double*[Nsteps+1];
  double **zleg1=new double*[Nsteps+1];
  double **xleg2=new double*[Nsteps+1];
  double **yleg2=new double*[Nsteps+1];
  double **zleg2=new double*[Nsteps+1];
  double **xleg3=new double*[Nsteps+1];
  double **yleg3=new double*[Nsteps+1];
  double **zleg3=new double*[Nsteps+1];
  
  /*Each protein is a rigid molecule, so geometry of legs will not change.*/
  double *leg1=new double[3];
  double *leg2=new double[3];
  double *leg3=new double[3];
  double *Mz=new double[9];
  double *M=new double[9];
  
  leg1[0]=0.0;
  leg1[1]=leglen;//orientation of first leg is just along y axis.
  leg1[2]=0.0;
  /*Now to create the other legs, rotate knees by theta around the z-axis*/
  double *zaxis=new double[3];
  zaxis[0]=0.0;
  zaxis[1]=0.0;
  zaxis[2]=1.0;

  calc_Rmatrix(zaxis, thetaleg,Mz);
  cout <<" leg1 vec: "<<leg1[0]<<' '<<leg1[1]<<' '<<leg1[2]<<endl;
  rotate(leg1, Mz, leg2);
  cout <<" leg1 vec: "<<leg2[0]<<' '<<leg2[1]<<' '<<leg2[2]<<endl;
  rotate(leg2, Mz, leg3);
  cout <<" leg3 vec: "<<leg3[0]<<' '<<leg3[1]<<' '<<leg3[2]<<endl;
  

  for(i=0;i<Nsteps+1;i++){
    xcm[i]=new double[N];
    ycm[i]=new double[N];
    zcm[i]=new double[N];
    xleg1[i]=new double[N];
    yleg1[i]=new double[N];
    zleg1[i]=new double[N];
    xleg2[i]=new double[N];
    yleg2[i]=new double[N];
    zleg2[i]=new double[N];
    xleg3[i]=new double[N];
    yleg3[i]=new double[N];
    zleg3[i]=new double[N];
    
  }
  cout <<"set...\n";
  double ty, tx, tz;
  double *v=new double[3];
  double *v2=new double[3];
  /*Generate initial coordinates*/
  double boxl=100.0;//this isn't important since particles will not be resstricted to motion inside any boundaries
  for(i=0;i<N;i++){
    xcm[0][i]=boxl*rand_gsl();
    ycm[0][i]=boxl*rand_gsl();
    zcm[0][i]=boxl*rand_gsl();
    
    xleg1[0][i]=xcm[0][i]+leg1[0];
    yleg1[0][i]=ycm[0][i]+leg1[1];
    zleg1[0][i]=zcm[0][i]+leg1[2];
    
    xleg2[0][i]=xcm[0][i]+leg2[0];
    yleg2[0][i]=ycm[0][i]+leg2[1];
    zleg2[0][i]=zcm[0][i]+leg2[2];
    
    xleg3[0][i]=xcm[0][i]+leg3[0];
    yleg3[0][i]=ycm[0][i]+leg3[1];
    zleg3[0][i]=zcm[0][i]+leg3[2];
  }

  double time;
  /*Propagate clathrins via translational and rotational diffusion*/
  for(t=1;t<Nsteps+1;t++){
    time=t*deltat;
    
    for(i=0;i<N;i++){
      /*Rotate and translate each leg, then the COM will be translated only*/
      dx=sqrt(2.0*deltat*Dtrans)*GaussV();
      dy=sqrt(2.0*deltat*Dtrans)*GaussV();
      dz=sqrt(2.0*deltat*Dtrans)*GaussV();
      
      tx=sqrt(2.0*deltat*Drot)*GaussV();
      ty=sqrt(2.0*deltat*Drot)*GaussV();
      tz=sqrt(2.0*deltat*Drot)*GaussV();
      
      rotationEuler(tx, ty, tz, M);
      
      /*New coordinates are rotation around previous COM, plus translation*/
      /*leg1*/
      v[0]=xleg1[t-1][i]-xcm[t-1][i];
      v[1]=yleg1[t-1][i]-ycm[t-1][i];
      v[2]=zleg1[t-1][i]-zcm[t-1][i];
      
      rotate(v, M, v2);
      
      xleg1[t][i]=xcm[t-1][i]+v2[0]+dx;      
      yleg1[t][i]=ycm[t-1][i]+v2[1]+dy;      
      zleg1[t][i]=zcm[t-1][i]+v2[2]+dz;      
      /*leg2*/
      v[0]=xleg2[t-1][i]-xcm[t-1][i];
      v[1]=yleg2[t-1][i]-ycm[t-1][i];
      v[2]=zleg2[t-1][i]-zcm[t-1][i];
      
      rotate(v, M, v2);
      
      xleg2[t][i]=xcm[t-1][i]+v2[0]+dx;      
      yleg2[t][i]=ycm[t-1][i]+v2[1]+dy;      
      zleg2[t][i]=zcm[t-1][i]+v2[2]+dz;      
      /*leg3*/
      v[0]=xleg3[t-1][i]-xcm[t-1][i];
      v[1]=yleg3[t-1][i]-ycm[t-1][i];
      v[2]=zleg3[t-1][i]-zcm[t-1][i];
      
      rotate(v, M, v2);
      
      xleg3[t][i]=xcm[t-1][i]+v2[0]+dx;      
      yleg3[t][i]=ycm[t-1][i]+v2[1]+dy;      
      zleg3[t][i]=zcm[t-1][i]+v2[2]+dz;      
      
      xcm[t][i]=xcm[t-1][i]+dx;
      ycm[t][i]=ycm[t-1][i]+dy;
      zcm[t][i]=zcm[t-1][i]+dz;
    }//loop over all coordinates
    
  }//loop over all time steps
    
  
  
  /*
    Now measure Root Mean Square displacement, averaged over all coordinates and
    over all time separations.
    Initialize rmsd arrays to zero, many of the values will remain
    zero if there are no configurations a givenn deltat apart.
  */ 
  double xcm0[N];
  double ycm0[N];
  double zcm0[N];
  
  
  int maxintdt=Nsteps/10;//maximum time separation between configs to measure RMSD
   
  double *rmsd=new double[maxintdt];
  int *ntimes=new int[maxintdt];
  

  
  /*Initialize rmsd arrays to zero before measuring for each Chunk of configs*/
      
  for(t=0;t<maxintdt;t++)
    {
      rmsd[t]=0;
      ntimes[t]=0; //to keep track of how many times you measured this t
    }
  cout << "initialized rmsd array. \n";
  
  
  
  int dint;
  for(t=0;t<Nsteps+1;t++){
    
      time=t*deltat;
      cout <<"time: "<<time<<endl;
      
      /*Set reference time zero*/
      
      for(i=0;i<N;i++)
	{
	  xcm0[i]=xleg1[t][i];
	  ycm0[i]=yleg1[t][i];
	  zcm0[i]=zleg1[t][i];
	}
      
      /*Now loop over possible time new*/
      
      for(t2=t+1;t2<Nsteps+1;t2++)
	{
	  
	  /*now we are at time2*/
	  
	  dint=t2-t;//change in time measured in steps.
	  
	  /*If this deltat is within range of interest, 
	    calculate rmsd over the N particles and keep
	    count of calculations at this deltat
	  */
	  if(dint<maxintdt)
	    {
	      for(i=0;i<N;i++)
		{
		  dx=xleg1[t2][i]-xcm0[i];
		  dy=yleg1[t2][i]-ycm0[i];
		  dz=zleg1[t2][i]-zcm0[i];
		  rmsd[dint]+=dx*dx+dy*dy+dz*dz;
		}
	      ntimes[dint]++;
	    }
	}//END LOOP OVER SECOND TIME STEP, NOW CHOOSE NEW TIMEZERO
      
      
    }//End choosing different time zeros
  
  
  /*Output results*/
  char file2[80];
  
  sprintf(file2, "rmsd_Dtran%g_Drot%g.out",Dtrans, Drot );
  ofstream outfile(file2);
  ofstream nfile("times.msd.out");
  for(t=1;t<maxintdt;t++)
    {
      delT=t*deltat;
	  
      outfile << delT<<' '<<rmsd[t]/(N*ntimes[t])<<endl;
      nfile<<delT<<' '<<ntimes[t]<<endl;
    }
  outfile.close();
  nfile.close();
  
  
  cout << "done. \n";
  
} //End main

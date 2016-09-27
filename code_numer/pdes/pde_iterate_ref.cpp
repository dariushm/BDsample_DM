/*


calculate the solution to the Smoluchowski equation
that is, the probability distribution as it evolves in space
and time
p(x,t) should approach peq=exp(-beta*v) at long times
also calculates the flux as appropriate derivative of p(x,t)
jx=-Dxx*(beta*dVdx*p+dpdt) (if Dxy=0)

uses reflecting boundary conditions
also uses symmetrized form of equation (p is scaled by peq^1/2)

oct 2010
maggie johnson
*/

#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>



using namespace std;
class Parms
{
public:
  int M;
  int Ntraj;
  int Nit;
  int Tstep;
  double deltatau;
  double kap1;
  double lambda;

};
void indexx(unsigned long n, double arr[], unsigned long indx[]);
double calc_potential(double x, double y);
double calc_bistable_potential(double x, double y, double beta);
double calc_bistable_xgradient(double x, double y, double beta);
double calc_bistable_ygradient(double x, double y, double beta);  
void write_potential(ofstream &outfile, ofstream &xfile, ofstream &yfile );
void write_coords(double *xs,int M, ofstream &outfile, ofstream &yfile);
void write_vector(double *xavg,double *xvec, int M, ofstream &outfile, ofstream &yfile);
void write_bistable_potential(ofstream &outfile, ofstream &xfile, ofstream &yfile, double beta);
void write_string(double *phi, int M, ofstream &outfile);
void calc_gradient(double *xs, double *fs, int M);
void calc_bistable_gradient(double *xs, double *fs, int M, double beta);

void integrate(double *xs, double *xstar, double *fs, int M, double beta, double deltat);
void voronoi(double *xs, double *xstar, double *phi, int M, double *xavg);
void voronoi_end(double *xs, double *phi, int M);
void diffusion1(double x, double y, double &Dxx, double &Dyy,double beta);
void diffusion_cons(double *xs,  double *Dxx, double *Dyy, double *dDx, double *dDy, double beta, int M);
void diffusionx(double x, double y, double &Dxx, double &dDx,double beta);
void diffusiony(double x, double y, double &Dyy, double &dDy,double beta);
void laplaceV(double x, double y, double beta, double &d2vx, double &d2vy);
double GaussV();
void read_parms(ifstream &parmfile, Parms &plist);
void write_parms(Parms &plist);
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define K 7
#define NSTACK 50

int main(int argc, char *argv[])
{
  int i,j;

  srand(static_cast<unsigned>(time(0)));

  
  /*calculate the flux*/
  //in the simpler version, calculate the flux without the eigenvector
  //peq(x,y)*D(x,y)*dphi/dx(x,y);
  double beta=1;
  double x, y;
  double v;
  double Dxx, Dyy;
  int M=131;
  int N=131;
  /*remember that for the lapack calls, the indexing of the matrices in array form is such
   * that the first N values are the first column, so A(i*N+j) is column i,
   * row j
   */
  double x1=-2;
  double x2=2;
  double xrange=x2-x1;
  double y1=-3.1;
  double y2=3.1;
  double yrange=y2-y1;
  double hx=xrange/((M-1)*1.0);
  double hy=yrange/((N-1)*1.0);
  double hx2=hx*hx;
  double hy2=hy*hy;
  double dDx, dDy;

  int size=M*N;//N*M is the size of each row/column
  double *p=new double[size];
  double *pnew=new double[size];
  double d2V, d2vx, d2vy;
  int m, n;
  int col, row;
  double dvx, dvy;
  ofstream pfile;
  char fname[100];
  /**********************************/
  /*Construct the operator matrix
   *operator gamma= -div(D*peq*grad(peq^-1)) [div is divergence]  
   *in symmetric version use the operator gamma_twidle=peq^-1/2*gamma*peq^1/2
   *calculate x and y components of the operator separately, no cross terms
   */
  /*create initial conditions*/
  //the boundary is where m=0 and where n=0
  for(n=0;n<N;n++){
    p[0*M+n]=0;
    p[(M-1)*M+n]=0;
  }
  for(m=0;m<M;m++){
    p[m*N+0]=0;
    p[m*N+N-1]=0;
  }
  /*calc peq*/
  double a=1;
  double del=5;
  double pi=acos(-1.0);
  double alpha=140*pi/180;
  
  double ct=1.0/tan(alpha/2);
  double normy=sqrt(pi*a*a/(beta*2*del*ct*ct)); //a=1
  double normx=0.8340588041;
  double norm=normx*normy;
  double *peq=new double[M*N];
  double *sqrtpeq=new double[M*N];
  for(m=0;m<M;m++){
    for(n=0;n<N;n++){
      x=(m)*hx+x1;
      y=(n)*hy+y1;
      v=calc_bistable_potential(x, y, beta);
      peq[m*N+n]=exp(-beta*v)/norm;
      sqrtpeq[m*N+n]=sqrt(peq[m*N+n]);
    }
  }
  ofstream eqfile("Peq.dat");

  for(m=0;m<M;m++){
    for(n=0;n<N;n++){
      eqfile<<peq[m*N+n]<<'\t';
    }
    eqfile<<endl;
  }

  /*these are the initial conditions for p(x,t=0)*/
  
  double ex=-1;
  double ey=0;
  double ax=0.1*0.1;
  double ay=0.25*0.25;
  double ell, xex, yey;
  double sx=ax;
  double sy=ay;
  double tot=0;
  
  double fact=1.0/(sqrt(4*pi*pi)*sx*sy);
  for(m=0;m<M;m++){
    for(n=0;n<N;n++){
      x=(m)*hx+x1;
      y=(n)*hy+y1;
      //put all the probability within the ellipse defined by (x-ex)^2/ax+(y-ey)^2/ay=1
      xex=x-ex;
      yey=y-ey;
      p[m*N+n]=0;
      ell=xex*xex/ax+yey*yey/ay;
      if(ell<=1){
	//within this ellipse
	p[m*N+n]=fact*exp(-xex*xex/(2*sx*sx)-yey*yey/(2*sy*sy));
	tot+=p[m*N+n]*hx*hy;
	p[m*N+n]/=sqrtpeq[m*N+n];
      }
    }
  }
  /*normalize p0*/
  for(m=0;m<M;m++){
    for(n=0;n<N;n++){
      p[m*N+n]/=(1.0*tot);
    }
  }
  double jx,jy; 
  double temp, temp2, temp0;
  ofstream pyfile;
  ofstream pxfile;
  char pname[100]; 

  //write out the answer
  sprintf(fname,"prob.initial");
  pfile.open(fname);
  for(m=0;m<M;m++){
    for(n=0;n<N;n++){
      pfile<<p[m*N+n]*sqrtpeq[m*N+n]<<'\t';
    }
    pfile<<endl;
  }
  pfile.close();
  
  /*by using the symmetric version we are actually calculating the evolution
   *of ptwidle = peq^(-1/2)*P(x,t)
   */
  int t;
  int Nit=1000000;
  double dt=1E-5;
  for(t=1;t<Nit;t++){
    //do boundary, enforce reflecting boundary conditions, x[-1]=x[+1]
    for(m=1;m<M-1;m++){
      //when y==0
      n=0;
      x=(m)*hx+x1;
      y=(n)*hy+y1;
      dvx=calc_bistable_xgradient(x,y,beta);
      diffusionx(x, y, Dxx, dDx, beta);
      diffusiony(x, y, Dyy, dDy, beta);
      dvy=calc_bistable_ygradient(x,y,beta);
      laplaceV(x,y,beta, d2vx, d2vy);
      pnew[m*N+n]=p[m*N+n];
      pnew[m*N+n]+=dt/hx2*Dxx*(p[(m+1)*N+n]+p[(m-1)*N+n]-2*p[m*N+n]);//first add in second derivative points at the old values
      pnew[m*N+n]+=dt/hy2*Dyy*(p[m*N+n+1]+p[m*N+n+1]-2*p[m*N+n]);//first add in second derivative points at the old values
      pnew[m*N+n]+=dt/(2*hx)*dDx*(p[(m+1)*N+n]-p[(m-1)*N+n]);//dpdx
      pnew[m*N+n]+=dt/(2*hy)*dDy*(p[m*N+n+1]-p[m*N+n+1]);//dpdy
      pnew[m*N+n]+=dt*(beta/2.0*(dvx*dDx+dvy*dDy+Dxx*d2vx+Dyy*d2vy)-beta*beta/4.0*(Dxx*dvx*dvx+Dyy*dvy*dvy))*p[m*N+n];
      
      n=N-1;
      x=(m)*hx+x1;
      y=(n)*hy+y1;
      dvx=calc_bistable_xgradient(x,y,beta);
      diffusionx(x, y, Dxx, dDx, beta);
      diffusiony(x, y, Dyy, dDy, beta);
      dvy=calc_bistable_ygradient(x,y,beta);
      laplaceV(x,y,beta, d2vx, d2vy);
      pnew[m*N+n]=p[m*N+n];
      pnew[m*N+n]+=dt/hx2*Dxx*(p[(m+1)*N+n]+p[(m-1)*N+n]-2*p[m*N+n]);//first add in second derivative points at the old values
      pnew[m*N+n]+=dt/hy2*Dyy*(p[m*N+n-1]+p[m*N+n-1]-2*p[m*N+n]);//first add in second derivative points at the old values
      pnew[m*N+n]+=dt/(2*hx)*dDx*(p[(m+1)*N+n]-p[(m-1)*N+n]);//dpdx
      pnew[m*N+n]+=dt/(2*hy)*dDy*(p[m*N+n-1]-p[m*N+n-1]);//dpdy
      pnew[m*N+n]+=dt*(beta/2*(dvx*dDx+dvy*dDy+Dxx*d2vx+Dyy*d2vy)-beta*beta/4*(Dxx*dvx*dvx+Dyy*dvy*dvy))*p[m*N+n];
      
    }
    for(n=1;n<N-1;n++){
      m=0;
      x=(m)*hx+x1;
      y=(n)*hy+y1;
      dvx=calc_bistable_xgradient(x,y,beta);
      diffusionx(x, y, Dxx, dDx, beta);
      diffusiony(x, y, Dyy, dDy, beta);
      dvy=calc_bistable_ygradient(x,y,beta);
      laplaceV(x,y,beta, d2vx, d2vy);
      pnew[m*N+n]=p[m*N+n];
      pnew[m*N+n]+=dt/hx2*Dxx*(p[(m+1)*N+n]+p[(m+1)*N+n]-2*p[m*N+n]);//first add in second derivative points at the old values
      pnew[m*N+n]+=dt/hy2*Dyy*(p[m*N+n+1]+p[m*N+n-1]-2*p[m*N+n]);//first add in second derivative points at the old values
      pnew[m*N+n]+=dt/(2*hx)*dDx*(p[(m+1)*N+n]-p[(m+1)*N+n]);//dpdx
      pnew[m*N+n]+=dt/(2*hy)*dDy*(p[m*N+n+1]-p[m*N+n-1]);//dpdy
      pnew[m*N+n]+=dt*(beta/2*(dvx*dDx+dvy*dDy+Dxx*d2vx+Dyy*d2vy)-beta*beta/4*(Dxx*dvx*dvx+Dyy*dvy*dvy))*p[m*N+n];
      m=M-1;
      x=(m)*hx+x1;
      y=(n)*hy+y1;
      dvx=calc_bistable_xgradient(x,y,beta);
      diffusionx(x, y, Dxx, dDx, beta);
      diffusiony(x, y, Dyy, dDy, beta);
      dvy=calc_bistable_ygradient(x,y,beta);
      laplaceV(x,y,beta, d2vx, d2vy);
      pnew[m*N+n]=p[m*N+n];
      pnew[m*N+n]+=dt/hx2*Dxx*(p[(m-1)*N+n]+p[(m-1)*N+n]-2*p[m*N+n]);//first add in second derivative points at the old values
      pnew[m*N+n]+=dt/hy2*Dyy*(p[m*N+n+1]+p[m*N+n-1]-2*p[m*N+n]);//first add in second derivative points at the old values
      pnew[m*N+n]+=dt/(2*hx)*dDx*(p[(m-1)*N+n]-p[(m-1)*N+n]);//dpdx
      pnew[m*N+n]+=dt/(2*hy)*dDy*(p[m*N+n+1]-p[m*N+n-1]);//dpdy
      pnew[m*N+n]+=dt*(beta/2*(dvx*dDx+dvy*dDy+Dxx*d2vx+Dyy*d2vy)-beta*beta/4*(Dxx*dvx*dvx+Dyy*dvy*dvy))*p[m*N+n];
    }
    //now do the corners
    //m=0, n=0
    m=0;
    n=0;
    x=(m)*hx+x1;
    y=(n)*hy+y1;
    dvx=calc_bistable_xgradient(x,y,beta);
    diffusionx(x, y, Dxx, dDx, beta);
    diffusiony(x, y, Dyy, dDy, beta);
    dvy=calc_bistable_ygradient(x,y,beta);
    laplaceV(x,y,beta, d2vx, d2vy);
    pnew[m*N+n]=p[m*N+n];
    pnew[m*N+n]+=dt/hx2*Dxx*(p[(m+1)*N+n]+p[(m+1)*N+n]-2*p[m*N+n]);//first add in second derivative points at the old values
    pnew[m*N+n]+=dt/hy2*Dyy*(p[m*N+n+1]+p[m*N+n+1]-2*p[m*N+n]);//first add in second derivative points at the old values
    pnew[m*N+n]+=dt/(2*hx)*dDx*(p[(m+1)*N+n]-p[(m+1)*N+n]);//dpdx
    pnew[m*N+n]+=dt/(2*hy)*dDy*(p[m*N+n+1]-p[m*N+n+1]);//dpdy
    pnew[m*N+n]+=dt*(beta/2*(dvx*dDx+dvy*dDy+Dxx*d2vx+Dyy*d2vy)-beta*beta/4*(Dxx*dvx*dvx+Dyy*dvy*dvy))*p[m*N+n];
    
    //m=0, n=N-1
    m=0;
    n=N-1;
    x=(m)*hx+x1;
    y=(n)*hy+y1;
    dvx=calc_bistable_xgradient(x,y,beta);
    diffusionx(x, y, Dxx, dDx, beta);
    diffusiony(x, y, Dyy, dDy, beta);
    dvy=calc_bistable_ygradient(x,y,beta);
    laplaceV(x,y,beta, d2vx, d2vy);
    pnew[m*N+n]=p[m*N+n];
    pnew[m*N+n]+=dt/hx2*Dxx*(p[(m+1)*N+n]+p[(m+1)*N+n]-2*p[m*N+n]);//first add in second derivative points at the old values
    pnew[m*N+n]+=dt/hy2*Dyy*(p[m*N+n-1]+p[m*N+n-1]-2*p[m*N+n]);//first add in second derivative points at the old values
    pnew[m*N+n]+=dt/(2*hx)*dDx*(p[(m+1)*N+n]-p[(m+1)*N+n]);//dpdx
    pnew[m*N+n]+=dt/(2*hy)*dDy*(p[m*N+n-1]-p[m*N+n-1]);//dpdy
    pnew[m*N+n]+=dt*(beta/2*(dvx*dDx+dvy*dDy+Dxx*d2vx+Dyy*d2vy)-beta*beta/4*(Dxx*dvx*dvx+Dyy*dvy*dvy))*p[m*N+n];
    
    //m=M-1, n=0
    m=M-1;
    n=0;
    x=(m)*hx+x1;
    y=(n)*hy+y1;
    dvx=calc_bistable_xgradient(x,y,beta);
    diffusionx(x, y, Dxx, dDx, beta);
    diffusiony(x, y, Dyy, dDy, beta);
    dvy=calc_bistable_ygradient(x,y,beta);
    laplaceV(x,y,beta, d2vx, d2vy);
    pnew[m*N+n]=p[m*N+n];
    pnew[m*N+n]+=dt/hx2*Dxx*(p[(m-1)*N+n]+p[(m-1)*N+n]-2*p[m*N+n]);//first add in second derivative points at the old values
    pnew[m*N+n]+=dt/hy2*Dyy*(p[m*N+n+1]+p[m*N+n+1]-2*p[m*N+n]);//first add in second derivative points at the old values
    pnew[m*N+n]+=dt/(2*hx)*dDx*(p[(m-1)*N+n]-p[(m-1)*N+n]);//dpdx
    pnew[m*N+n]+=dt/(2*hy)*dDy*(p[m*N+n+1]-p[m*N+n+1]);//dpdy
    pnew[m*N+n]+=dt*(beta/2*(dvx*dDx+dvy*dDy+Dxx*d2vx+Dyy*d2vy)-beta*beta/4*(Dxx*dvx*dvx+Dyy*dvy*dvy))*p[m*N+n];
    
    //m=M-1, n=N-1
    m=M-1;
    n=N-1;
    x=(m)*hx+x1;
    y=(n)*hy+y1;
    dvx=calc_bistable_xgradient(x,y,beta);
    diffusionx(x, y, Dxx, dDx, beta);
    diffusiony(x, y, Dyy, dDy, beta);
    dvy=calc_bistable_ygradient(x,y,beta);
    laplaceV(x,y,beta, d2vx, d2vy);
    pnew[m*N+n]=p[m*N+n];
    pnew[m*N+n]+=dt/hx2*Dxx*(p[(m-1)*N+n]+p[(m-1)*N+n]-2*p[m*N+n]);//first add in second derivative points at the old values
    pnew[m*N+n]+=dt/hy2*Dyy*(p[m*N+n-1]+p[m*N+n-1]-2*p[m*N+n]);//first add in second derivative points at the old values
    pnew[m*N+n]+=dt/(2*hx)*dDx*(p[(m-1)*N+n]-p[(m-1)*N+n]);//dpdx
    pnew[m*N+n]+=dt/(2*hy)*dDy*(p[m*N+n-1]-p[m*N+n-1]);//dpdy
    pnew[m*N+n]+=dt*(beta/2*(dvx*dDx+dvy*dDy+Dxx*d2vx+Dyy*d2vy)-beta*beta/4*(Dxx*dvx*dvx+Dyy*dvy*dvy))*p[m*N+n];
      
    
    /*ALL INTERIOR MESH POINTS*/
    
    for(m=1;m<M-1;m++){
      
      for(n=1;n<N-1;n++){
	/*calculate new time points based on previous time points*/
	x=(m)*hx+x1;
	y=(n)*hy+y1;
	dvx=calc_bistable_xgradient(x,y,beta);
	diffusionx(x, y, Dxx, dDx, beta);
	diffusiony(x, y, Dyy, dDy, beta);
	dvy=calc_bistable_ygradient(x,y,beta);
	laplaceV(x,y,beta, d2vx, d2vy);
	pnew[m*N+n]=p[m*N+n];
	pnew[m*N+n]+=dt/hx2*Dxx*(p[(m+1)*N+n]+p[(m-1)*N+n]-2*p[m*N+n]);//first add in second derivative points at the old values
	pnew[m*N+n]+=dt/hy2*Dyy*(p[m*N+n+1]+p[m*N+n-1]-2*p[m*N+n]);//first add in second derivative points at the old values
	pnew[m*N+n]+=dt/(2*hx)*dDx*(p[(m+1)*N+n]-p[(m-1)*N+n]);//dpdx
	pnew[m*N+n]+=dt/(2*hy)*dDy*(p[m*N+n+1]-p[m*N+n-1]);//dpdy
	pnew[m*N+n]+=dt*(beta/2*(dvx*dDx+dvy*dDy+Dxx*d2vx+Dyy*d2vy)-beta*beta/4*(Dxx*dvx*dvx+Dyy*dvy*dvy))*p[m*N+n];
      }
    }//end looping over x and y
    for(m=1;m<M-1;m++){
      for(n=1;n<N-1;n++){
	p[m*N+n]=pnew[m*N+n];
      }
    }
    if(t%10000==0){
      //write out the answer
      sprintf(fname,"prob.%d",t);
      pfile.open(fname);
      for(m=0;m<M;m++){
	for(n=0;n<N;n++){
	  pfile<<p[m*N+n]*sqrtpeq[m*N+n]<<'\t';
	}
	pfile<<endl;
      }
      pfile.close();
      
      sprintf(fname,"xflux.%d",t);
      pxfile.open(fname);
      sprintf(fname,"yflux.%d",t);
      pyfile.open(fname);
      for(m=1;m<M-1;m++){
	for(n=1;n<N-1;n++){
	  x=(m)*hx+x1;
	  y=(n)*hy+y1;
	  dvx=calc_bistable_xgradient(x,y,beta);
	  diffusionx(x, y, Dxx, dDx, beta);
	  diffusiony(x, y, Dyy, dDy, beta);
	  dvy=calc_bistable_ygradient(x,y,beta);
	  
	  temp0=p[m*N+n]*sqrtpeq[m*N+n];
	  temp=p[(m+1)*N+n]*sqrtpeq[(m+1)*N+n];
	  temp2=p[(m-1)*N+n]*sqrtpeq[(m-1)*N+n];
	  jx=-Dxx*(temp-temp2)/(2*hx)-Dxx*beta*dvx*temp0;
	  pxfile<<jx<<'\t';
	  
	  temp=p[m*N+n+1]*sqrtpeq[m*N+n+1];
	  temp2=p[m*N+n+1]*sqrtpeq[m*N+n-1];
	  jy=-Dyy*(temp-temp2)/(2*hy)-Dyy*beta*dvy*temp0;
	  pyfile<<jy<<'\t';
	  
	}
	pxfile<<endl;
	pyfile<<endl;
      }
      pxfile.close();
      pyfile.close();
    }
  }

  //write out the answer
  sprintf(fname,"prob.final");
  pfile.open(fname);
  for(m=0;m<M;m++){
    for(n=0;n<N;n++){
      pfile<<p[m*N+n]*sqrtpeq[m*N+n]<<'\t';
    }
    pfile<<endl;
  }
  pfile.close();
  
  sprintf(fname,"xflux.final");
  pxfile.open(fname);
  sprintf(fname,"yflux.final");
  pyfile.open(fname);
  for(m=1;m<M-1;m++){
    for(n=1;n<N-1;n++){
      x=(m)*hx+x1;
      y=(n)*hy+y1;
      dvx=calc_bistable_xgradient(x,y,beta);
      diffusionx(x, y, Dxx, dDx, beta);
      diffusiony(x, y, Dyy, dDy, beta);
      dvy=calc_bistable_ygradient(x,y,beta);
      
      temp0=p[m*N+n]*sqrtpeq[m*N+n];
      temp=p[(m+1)*N+n]*sqrtpeq[(m+1)*N+n];
      temp2=p[(m-1)*N+n]*sqrtpeq[(m-1)*N+n];
      jx=-Dxx*(temp-temp2)/(2*hx)-Dxx*beta*dvx*temp0;
      pxfile<<jx<<'\t';
      
      temp=p[m*N+n+1]*sqrtpeq[m*N+n+1];
      temp2=p[m*N+n+1]*sqrtpeq[m*N+n-1];
      jy=-Dyy*(temp-temp2)/(2*hy)-Dyy*beta*dvy*temp0;
      pyfile<<jy<<'\t';
      
    }
    pxfile<<endl;
    pyfile<<endl;
  }
  pxfile.close();
  pyfile.close();
  
       
}/*End main*/


double calc_potential(double x, double y)
{
  /*Muller potential*/
  double xmin1=x-1;
  double ymin1=y-1;
  double y2=y*y;
  double x2=x*x;
  double a=(y-0.5);
  double b=(x+0.5);
  double c=(y-1.5);
  double d=(x+1);
  double pi=acos(-1.0);
  double vx=-200*exp(-xmin1*xmin1-10*y2)-100*exp(-x2-10*a*a)-170*exp(-6.5*b*b+11*b*c-6.5*c*c)+15*exp(0.7*d*d+0.6*d*ymin1+0.7*ymin1*ymin1);
  vx+=9*sin(2*pi*5*x)*sin(2*pi*5*y);
  return vx;

}
double calc_bistable_potential(double x, double y, double beta)
{

  double a=1;
  double del=5;
  double pi=acos(-1.0);
  double alpha=140*pi/180;
  
  double xa=x/a;
  double x21=(xa*xa-1);
  double yc=(y/a)*1.0/tan(alpha/2);
  double vx=(del*x21*x21+2*del*yc*yc)/beta;
  
  return vx;

}
void calc_gradient(double *xs, double *fs, int M)
{
  
  double xmin1;
  double ymin1;
  double y2;
  double x2;
  double a;
  double b;
  double c;
  double d;
  double pi=acos(-1.0);
  double vx;
  double sinx;
  double siny;
  int i;
  double x, y;
  for(i=0;i<M;i++)
    {
      x=xs[0*M+i];
      y=xs[1*M+i];
      xmin1=x-1;
      ymin1=y-1;
      y2=y*y;
      x2=x*x;
      a=(y-0.5);
      b=(x+0.5);
      c=(y-1.5);
      d=(x+1);
      sinx=sin(2*pi*5*x);
      siny=sin(2*pi*5*y);
      //  vx=-200*exp(-xmin1*xmin1-10*y2)-100*exp(-x2-10*a*a)-170*exp(-6.5*b*b+11*b*c-6.5*c*c)+15*exp(0.7*d*d+0.6*d*ymin1+0.7*ymin1*ymin1);
      //vx+=9*sinx*siny;
      fs[0*M+i]=-200*exp(-xmin1*xmin1-10*y2)*(-2*xmin1)-100*exp(-x2-10*a*a)*(-2*x)-170*exp(-6.5*b*b+11*b*c-6.5*c*c)*(-13*b+11*c)+15*exp(0.7*d*d+0.6*d*ymin1+0.7*ymin1*ymin1)*(1.4*d+0.6*ymin1);//gradx
      fs[0*M+i]+=9*2*pi*5*cos(2*pi*5*x)*siny;
      fs[1*M+i]=-200*exp(-xmin1*xmin1-10*y2)*(-20*y)-100*exp(-x2-10*a*a)*(-20*a)-170*exp(-6.5*b*b+11*b*c-6.5*c*c)*(11*b-13*c)+15*exp(0.7*d*d+0.6*d*ymin1+0.7*ymin1*ymin1)*(0.6*d+1.4*ymin1);//gradY
      fs[1*M+i]+=9*2*pi*5*cos(2*pi*5*y)*sinx;
      fs[0*M+i]*=-1;
      fs[1*M+i]*=-1;
      
    }
  
  

}
void calc_bistable_gradient(double *xs, double *fs, int M, double beta)
{
  
  double a=1;
  double del=5;
  double pi=acos(-1.0);
  double alpha=140*pi/180;
  
  double xa;
  double x21;
  double ct;
  double ya;
  int i;
  double x, y;
  for(i=0;i<M;i++)
    {
      x=xs[0*M+i];
      y=xs[1*M+i];
      xa=x/a;
      ya=y/a;
      x21=(xa*xa-1);
      ct=1.0/tan(alpha/2);
      
      //take negative gradient for force
      fs[0*M+i]=-4*del*x21*x/(a*a)/beta;
      fs[1*M+i]=-4*del*ya*ct*ct/a/beta;
      
    }
  
  

}
double calc_bistable_xgradient(double x, double y, double beta) 
{
  
  double a=1;
  double del=5;
  double pi=acos(-1.0);
  double alpha=140*pi/180;
  
  double xa;
  double x21;
  double ct;
  double ya;
  int i;

  xa=x/a;
  ya=y/a;
  x21=(xa*xa-1);
  ct=1.0/tan(alpha/2);
  
  double dvx=4*del*x21*x/(a*a)/beta;
  //fs[1*M+i]=-4*del*ya*ct*ct/a/beta;
  return dvx;
}
double calc_bistable_ygradient(double x, double y, double beta) 
{
  
  double a=1;
  double del=5;
  double pi=acos(-1.0);
  double alpha=140*pi/180;
  
  double xa;
  double x21;
  double ct;
  double ya;
  int i;

  xa=x/a;
  ya=y/a;
  x21=(xa*xa-1);
  ct=1.0/tan(alpha/2);
  
  //dvx=4*del*x21*x/(a*a)/beta;
  double dvy=4*del*ya*ct*ct/a/beta;
  return dvy;
  
}
void laplaceV(double x, double y, double beta, double &d2vx, double &d2vy) 
{
  
  double a=1;
  double del=5;
  double pi=acos(-1.0);
  double alpha=140*pi/180;
  
  double xa;
  double x21;
  double ct;
  double ya;
  int i;

  xa=x/a;
  //  ya=y/a;
  //  x21=(xa*xa-1);
  ct=1.0/tan(alpha/2);
  double a2=a*a;
  d2vx=4*del/a2/beta*(3*xa*xa-1);
  d2vy=4*del*ct*ct/a2/beta;
  
}

void integrate(double *xs, double *xstar, double *fs, int M, double beta, double deltat)
{
  int i;
  for(i=0;i<M;i++)
    {
      xstar[0*M+i]=xs[0*M+i]-deltat*fs[0*M+i]+sqrt(2.0/beta*deltat)*GaussV();
      xstar[1*M+i]=xs[1*M+i]-deltat*fs[1*M+i]+sqrt(2.0/beta*deltat)*GaussV();
    }
}
void voronoi(double *xs, double *xstar, double *phi, int M, double *xavg)
{
  double xtmp, ytmp;
  int i, j;
  double dx, dy;
  double d2, g2;
  for(i=0;i<M;i++)
    {
      dx=xstar[0*M+i]-phi[0*M+i];
      dy=xstar[1*M+i]-phi[1*M+i];
      d2=dx*dx+dy*dy; //distance between you and your string
      j=0;
      xtmp=xstar[0*M+i];
      ytmp=xstar[1*M+i];//if d2 is minimum distance, then set xs equal to xstar

      while(j<M)
	{
	  if(j!=i)
	    {
	      dx=xstar[0*M+i]-phi[0*M+j];
	      dy=xstar[1*M+i]-phi[1*M+j];
	      g2=dx*dx+dy*dy; //distance between you and other string elements
	      
	      if(d2>g2)
		{
		  //we are not in the voronoi region, so reject this move
		  j=M;//exit from loop over other string elements
		  //xs does not change! 
		  xtmp=xs[0*M+i];
		  ytmp=xs[1*M+i];
		}
	      j++;
	    }
	  else 
	    j++;
	}
      //update position of xs, to either xstar, or staying at previous xs 
      xs[0*M+i]=xtmp;
      xs[1*M+i]=ytmp;
      /*According to paper, should only update averages to previous coordinates*/
      //update running average of position of xs

      //      xavg[0*M+i]+=xs[0*M+i];
      //xavg[1*M+i]+=xs[1*M+i];
      
    }
	  
	      



} 
void voronoi_end(double *xs, double *phi, int M)
{
  double xtmp, ytmp;
  int i, j;
  double dx, dy;
  double d2, g2;
  for(i=0;i<M;i++)
    {
      dx=xs[0*M+i]-phi[0*M+i];
      dy=xs[1*M+i]-phi[1*M+i];
      d2=dx*dx+dy*dy; //distance between you and your string
      j=0;
      xtmp=xs[0*M+i];
      ytmp=xs[1*M+i];//if d2 is minimum distance, then set xs equal to xstar
      while(j<M)
	{
	  if(j!=i)
	    {
	      dx=xs[0*M+i]-phi[0*M+j];
	      dy=xs[1*M+i]-phi[1*M+j];
	      g2=dx*dx+dy*dy; //distance between you and other string elements
	      if(d2>g2)
		{
		  //we are not in the voronoi region, so reject this move
		  j=M;//exit from loop over other string elements
		  //xs does not change! 
		  xtmp=phi[0*M+i];
		  ytmp=phi[1*M+i];
		}
	      j++;
	    }
	  else 
	    j++;
	}
      //update position of xs, to either xstar, or staying at previous xs 
      xs[0*M+i]=xtmp;
      xs[1*M+i]=ytmp;
    }

} 

double GaussV()
{
  /*Box mueller method for gaussian distributed random number from a uniform
    random number generator~ rand()*/
  
  double R=2.0;
  double rnum;
  double V1, V2;

  while(R>=1.0)
    {

      V1=2.0*rand()/RAND_MAX-1.0;
      V2=2.0*rand()/RAND_MAX-1.0;
      R=V1*V1+V2*V2;
    }
  rnum=V1*sqrt(-2.0*log(R)/R);
  return rnum;

}
void write_potential(ofstream &outfile, ofstream &xfile, ofstream &yfile )
{
  int i, j;
  double xmin=-1.5;
  double xmax=1;
  double ymin=-0.2;
  double ymax=2;
  int N=200;
  double delx=(xmax-xmin)/(N*1.0);
  double dely=(ymax-ymin)/(N*1.0);


  double x, y, vx;
  for(i=0;i<N;i++)
    {
      for(j=0;j<N;j++)
	{
	  x=xmin+j*delx;
	  y=ymin+i*dely;
	  //vx=calc_bistable_potential(x, y);
	  vx=calc_potential(x, y);
	  if(vx>100)vx=100;
	  outfile <<vx<<'\t';
	}
      outfile<<endl;
    }
  for(i=0;i<N;i++)
    {
      xfile<<xmin+i*delx<<endl;
      yfile<<ymin+i*dely<<endl;
    }
  //clear out upper corner of matrix where values are so large
  
}
void write_bistable_potential(ofstream &outfile, ofstream &xfile, ofstream &yfile, double beta)
{
  int i, j;
  double xmin=-2;
  double xmax=2;
  double ymin=-2;
  double ymax=2;
  int N=200;
  double delx=(xmax-xmin)/(N*1.0);
  double dely=(ymax-ymin)/(N*1.0);


  double x, y, vx;
  for(i=0;i<N;i++)
    {
      for(j=0;j<N;j++)
	{
	  x=xmin+j*delx;
	  y=ymin+i*dely;
	  vx=calc_bistable_potential(x, y, beta);
	  //vx=calc_potential(x, y);
	  //if(vx>100)vx=100;
	  outfile <<vx<<'\t';
	}
      outfile<<endl;
    }
  for(i=0;i<N;i++)
    {
      xfile<<xmin+i*delx<<endl;
      yfile<<ymin+i*dely<<endl;
    }
  //clear out upper corner of matrix where values are so large
  
}
void write_string(double *phi, int M, ofstream &outfile)
{
  int i;
  for(i=0;i<M;i++)
    {
      outfile<<phi[0*M+i]<<' '<<phi[1*M+i]<<endl;
    }

}

void diffusion1(double x, double y, double &Dxx, double &Dyy,double beta)
{
  //two by two, Dxx, Dxy, Dyx, Dyy
  int i, j;
  double kT=1.0/beta;
  double norm;
  double Eoff=0;
  double Exx, Eyy;
  double Einf=1;
  double delE=8;
  double sig=0.2;
  double sig2=sig*sig;
  double e1;
  
  norm=x*x+y*y;
  e1=exp(-norm/(2*sig2));
  Exx=Einf+delE*e1;
  Eyy=Einf+delE*e1;
  Dxx=kT/Exx;
  Dyy=kT/Eyy;
  //Dxy=0.0;
  //dDx[i]=kT*2*x*delE/(2*sig2)*e1/(Exx*Exx);
  //dDy[i]=kT*2*y*delE/(2*sig2)*e1/(Eyy*Eyy);
    
}
void diffusionx(double x, double y, double &Dxx, double &dDx,double beta)
{
  //two by two, Dxx, Dxy, Dyx, Dyy
  int i, j;
  double kT=1.0/beta;
  double norm;
  double Eoff=0;
  double Exx, Eyy;
  double Einf=1;
  double delE=8;
  double sig=0.2;
  double sig2=sig*sig;
  double e1;
  
  norm=x*x+y*y;
  e1=exp(-norm/(2*sig2));
  Exx=Einf+delE*e1;
  //Eyy=Einf+delE*e1;
  Dxx=kT/Exx;
  //Dyy=kT/Eyy;
  // Dxy=0.0;
  dDx=kT*2*x*delE/(2*sig2)*e1/(Exx*Exx);
  //dDy[i]=kT*2*y*delE/(2*sig2)*e1/(Eyy*Eyy);
    
}
void diffusiony(double x, double y, double &Dyy, double &dDy,double beta)
{
  //two by two, Dxx, Dxy, Dyx, Dyy
  int i, j;
  double kT=1.0/beta;
  double norm;
  
  double Eoff=0;
  double Exx, Eyy;
  double Einf=1;
  double delE=8;
  double sig=0.2;
  double sig2=sig*sig;
  double e1;
  
  norm=x*x+y*y;
  e1=exp(-norm/(2*sig2));
  //Exx=Einf+delE*e1;
  Eyy=Einf+delE*e1;
  //Dxx=kT/Exx;
  Dyy=kT/Eyy;
  // Dxy=0.0;
  //dDx=kT*2*x*delE/(2*sig2)*e1/(Exx*Exx);
  dDy=kT*2*y*delE/(2*sig2)*e1/(Eyy*Eyy);
    
}

void diffusion_cons(double *xs, double *Dxx, double *Dyy, double *dDx, double *dDy, double beta, int M)
{
  //two by two, Dxx, Dxy, Dyx, Dyy
  int i, j;
  double kT=1.0/beta;
  double norm;
  double x, y;
  double Eoff=0;
  double Exx, Eyy;
  double Einf=1;
  double delE=8;
  double sig=0.2;
  double sig2=sig*sig;
  double e1;
  for(i=0;i<M;i++)
    {
      x=xs[0*M+i];
      y=xs[1*M+i];
      norm=x*x+y*y;
      //      e1=exp(-norm/(2*sig2));
      Exx=1;
      Eyy=1;
      Dxx[i]=kT/Exx;
      Dyy[i]=kT/Eyy;
      
      dDx[i]=0;
      dDy[i]=0;
    }
}
void read_parms(ifstream &parmfile, Parms &plist)
{
  parmfile >>plist.M;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Ntraj;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Tstep;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Nit;
  parmfile.ignore(400,'\n');
  parmfile >>plist.deltatau;
  parmfile.ignore(400,'\n');
  parmfile >>plist.kap1;
  parmfile.ignore(400,'\n');
  parmfile >>plist.lambda;
  parmfile.ignore(400,'\n');

}
void write_parms(Parms &plist)
{
  cout <<plist.M<<endl;
  cout <<plist.Ntraj<<endl;
  cout <<plist.Tstep<<endl;
  cout <<plist.Nit<<endl;
  cout <<plist.deltatau<<endl;
  cout <<plist.kap1<<endl;
  cout <<plist.lambda<<endl;

}


void indexx(unsigned long n, double arr[], unsigned long indx[])
{
  unsigned long i, indxt, ir=n,itemp, j, k, l=1;
  int jstack=0;
  double a;
  int *istack=new int[NSTACK];

  for(j=1;j<=n;j++)indx[j]=j;
  for(;;){
    if(ir-l<K){
      for(j=l+1;j<=ir;j++){
	indxt=indx[j];
	a=arr[indxt];
	for(i=j-1;i>=l;i--){
	  if(arr[indx[i]]<=a)break;
	  indx[i+1]=indx[i];
	}
	indx[i+1]=indxt;
      }
      if(jstack==0)break;
      ir=istack[jstack--];
      l=istack[jstack--];
    }else {
      k=(l+ir)>>1;
      SWAP(indx[k],indx[l+1]);
      if(arr[indx[l]]>arr[indx[ir]]){
	SWAP(indx[l],indx[ir]);
      }
      if(arr[indx[l+1]]>arr[indx[ir]]){
	SWAP(indx[l+1],indx[ir]);
      }
      if(arr[indx[l]]>arr[indx[l+1]]){
	SWAP(indx[l],indx[l+1]);
      }
      i=l+1;
      j=ir;
      indxt=indx[l+1];
      a=arr[indxt];
      for(;;){
	do i++; while(arr[indx[i]]<a);
	do j--; while(arr[indx[j]]>a);
	if(j<i)break;
	SWAP(indx[i],indx[j]);
      }
      indx[l+1]=indx[j];
      indx[j]=indxt;
      jstack+=2;
      if(jstack>NSTACK)printf("NSTACK too small in indexx");
      if(ir-i+1 >= j-l){
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
  delete[] istack;
}


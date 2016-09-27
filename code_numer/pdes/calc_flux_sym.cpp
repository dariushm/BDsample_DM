/*


program to calculate the flux from the smoluchowski equation,
by first finding the eigenvectors of the Smol. time evolution 
operator via finite difference on a grid.
Then use J=peq*dE/dx*D to find the flux field

It would be faster instead of using the LAPACK routines to use
a sparse matrix solver. The matrix is block tridiagonal but not 
quite symmetric because of the boundary conditions and therefore
also cannot be solved using the band symmetric solvers.
ARPACK has algorithms for sparse matrices which might
be much faster for diagonalizing this matrix

Uses symmetrized form of Smoluchowski
operator gamma= -div(D*peq*grad(peq^-1)) [div is divergence]  
in symmetric version use the operator gamma_twidle=peq^-1/2*gamma*peq^1/2
Means that eigenvectors changed by factor peq^1/2 
calculate x and y components of the operator separately, no cross terms
  
oct 2010
maggie johnson
*/

#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include "/usr/local/intel/Compiler/mkl_current/include/mkl.h"


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
  /*perform finite difference to evaluate the 1st eigenvector of the smoluchowski time evolution operator*/
  /*L*psi=lambda*psi, where psi is phi*peq, and phi is the eigenvector
   *the eigenvector is dependent on x and y, so it will be a column vector of length N*M, 
   *where N is length of x grid points, and M is y grid points
   */
  int M=131;
  int N=131;
  /*remember that for the lapack calls, the indexing of the matrices in array form is such
   * that the first N values are the first column, so A(i*N+j) is column i,
   * row j
   */
  double x1=-2;
  double x2=2;
  double xrange=x2-x1;
  double y1=-2.7;
  double y2=2.7;
  double yrange=y2-y1;
  double hx=xrange/((M-1)*1.0);
  double hy=yrange/((N-1)*1.0);
  double hx2=hx*hx;
  double hy2=hy*hy;
  double dDx, dDy;

  int size=M*N;//N*M is the size of each row/column
  double d2V, d2vx, d2vy;
  double *A=new double[size*size];//but A is sparse!
  for(i=0;i<size*size;i++)
    A[i]=0.0;
  int m, n;
  int col, row;
  double dvx, dvy;
  /**********************************/
  /*Construct the operator matrix
   *operator gamma= -div(D*peq*grad(peq^-1)) [div is divergence]  
   *in symmetric version use the operator gamma_twidle=peq^-1/2*gamma*peq^1/2
   *calculate x and y components of the operator separately, no cross terms
   */
      
  for(m=1;m<M-1;m++){
    /*this m loops over the M diagonal block
     * then n loops within this N*N block along the diagonal 
     */
    for(n=0;n<N;n++){
      row=m*N+n;//start of this block=m*N
      col=m*N+n;
      x=(m)*hx+x1;
      y=(n)*hy+y1;
      dvx=calc_bistable_xgradient(x,y,beta);
      diffusionx(x, y, Dxx, dDx, beta);
      
      /*first do second derivatives d2p/dx2=(p(x+h)+p(x-h)-2p(x))/h^2*/
      A[col*size+row]+=2*Dxx/hx2;//diagonal
      A[(col+N)*size+row]-=1*Dxx/hx2; //+h
      A[(col-N)*size+row]-=1*Dxx/hx2; //-h 
      
      /*now do first derivative in x: dp/dx=(p(x+h)-p(x-h))/2h */
      A[(col+N)*size+row]-=1*(dDx)/(2*hx);//+h
      A[(col-N)*size+row]+=1*(dDx)/(2*hx);//-h
      
    }
  }
  
  /*Now do the first block, m=0, incorporating the reflecgin BC, and do the 
   *final block, m=M-1, also with reflecting BC
   */
    
  for(n=0;n<N;n++){
    /*FIRST BLOCK m=0*/
    row=0*N+n;
    col=0*N+n;
    x=(0)*hx+x1;
    y=(n)*hy+y1;
    diffusionx(x, y, Dxx, dDx, beta);
    dvx=calc_bistable_xgradient(x,y,beta);
    
    /*second derv in x*/
    A[col*size+row]+=2*Dxx/hx2;
    A[(col+N)*size+row]-=1*Dxx/hx2; //+h 
      
    /*add in reflecting boundary condition for 2nd derv*/
    A[(col+N)*size+row]-=1*Dxx/hx2;//the -h effects the +h term
    
    /*Now first derv*/
    A[(col+N)*size+row]-=1*(dDx)/(2*hx); //+h
    
    /*Now add in reflecting BC for 1st derv*/
    A[(col+N)*size+row]+=1*(dDx)/(2*hx); //-h effects +h
        
    /*FINAL BLOCK m=M-1*/
    row=(M-1)*N+n;
    col=(M-1)*N+n;
    x=((M-1))*hx+x1;
    y=(n)*hy+y1;
    diffusionx(x, y, Dxx, dDx, beta);
    dvx=calc_bistable_xgradient(x,y,beta);
    
    /*second derv in x*/
    A[col*size+row]+=2*Dxx/hx2;
    A[(col-N)*size+row]-=1*Dxx/hx2; //-h 
      
    /*add in reflecting boundary condition for 2nd derv*/
    A[(col-N)*size+row]-=1*Dxx/hx2;//the +h effects the -h term
    
    /*Now first derv*/
    A[(col-N)*size+row]+=1*(dDx)/(2*hx);
    
    /*Now add in reflecting BC for 1st derv*/
    A[(col-N)*size+row]-=1*(dDx)/(2*hx);
    
  }

  /*Y dimension*/
  for(m=0;m<M;m++){
    /*this m loops over the M diagonal block
     * then n loops within this N*N block along the diagonal 
     */
    for(n=1;n<N-1;n++){
      row=m*N+n;
      col=m*N+n;
      x=(m)*hx+x1;
      y=(n)*hy+y1;
      diffusiony(x, y, Dyy, dDy, beta);
      dvy=calc_bistable_ygradient(x,y,beta);

      /*now do the second derivative in y d2p/dy2=(p(y+h)+p(y-h)-2p(y))/h^2*/
      A[col*size+row]+=2*Dyy/hy2;//diag
      A[(col+1)*size+row]-=1*Dyy/hy2;//+h
      A[(col-1)*size+row]-=1*Dyy/hy2;//-h
      
      /*now first derivative in y*/
      A[(col+1)*size+row]-=1*(dDy)/(2*hy);//+h
      A[(col-1)*size+row]+=1*(dDy)/(2*hy);//-h
    }
  }
  
  /*last N element of each M block*/
  for(m=0;m<M;m++){
    /*FIRST COLUMN n=0 */
    row=m*N+0;
    col=m*N+0;
    x=(m)*hx+x1;
    y=(0)*hy+y1;
    diffusiony(x, y, Dyy, dDy, beta);
    dvy=calc_bistable_ygradient(x,y,beta);
    
    /*second derv in y*/
    A[col*size+row]+=2*Dyy/hy2;
    A[(col+1)*size+row]-=1*Dyy/hy2;//+h
    
    /*reflecting BC*/
    A[(col+1)*size+row]-=1*Dyy/hy2;//-h effects +h
    
    /*Now first derv*/
    A[(col+1)*size+row]-=1*(dDy)/(2*hy); //+h
    
    /*Now add in reflecting BC for 1st derv*/
    A[(col+1)*size+row]+=1*(dDy)/(2*hy); //-h effects +h
    
    /*FINAL COLUMN n=N-1*/
    row=m*N+N-1;
    col=m*N+N-1;
    x=(m)*hx+x1;
    y=((N-1))*hy+y1;
    diffusiony(x, y, Dyy, dDy, beta);
    dvy=calc_bistable_ygradient(x,y,beta);
    
    /*second derv in y*/
    A[col*size+row]+=2*Dyy/hy2;
    A[(col-1)*size+row]-=1*Dyy/hy2;//-h
    
    /*reflecting BC*/
    A[(col-1)*size+row]-=1*Dyy/hy2;//+h effects -h
    
    /*Now first derv*/
    A[(col-1)*size+row]+=1*(dDy)/(2*hy); //-h
    
    /*Now add in reflecting BC for 1st derv*/
    A[(col-1)*size+row]-=1*(dDy)/(2*hy); //+h effects -h
    
  }
  
  /*Now position part of operator*/
  for(m=0;m<M;m++){
    for(n=0;n<N;n++){
      row=m*N+n;
      col=m*N+n;
      x=(m)*hx+x1;
      y=(n)*hy+y1;
      diffusiony(x, y, Dyy, dDy, beta);
      diffusionx(x, y, Dxx, dDx, beta);
      dvy=calc_bistable_ygradient(x,y,beta);
      dvx=calc_bistable_xgradient(x,y,beta);
      laplaceV(x,y,beta, d2vx, d2vy);
      /*coefficient includes laplacian of potential*/
      /*function value itself*/
      A[col*size+row]-=(beta/2*(dvx*dDx+dvy*dDy+Dxx*d2vx+Dyy*d2vy)-beta*beta/4*(Dxx*dvx*dvx+Dyy*dvy*dvy));
    }
  }
  ofstream afile("Amatrix.out");
  cout <<"finished with assigning the A matrix"<<endl;
  for(i=0;i<size;i++){
    for(j=0;j<size;j++){
      afile<<A[j*size+i]<<'\t';
    }
    afile<<endl;
  }
  /*Enforce boundary conditions, for dirichlet boundary conditions,
   * zero at the boundray, the matrix doesnt' change. for reflectin
   *boundary conditions, need to add terms to all p(0,y), p(x, 0)
   * p(M,y), p(x,N) values. such that dPdx=0, dpdy=0;
   */
  
  /*Now the A matrix is filled out, calculate the eigenvectors
   *although we only need the first one
   */
  /*Test whether A*peq=0, as it should*/
  /*calc peq*/
  double a=1;
  double del=5;
  double pi=acos(-1.0);
  double alpha=140*pi/180;
  
  double ct=cot(alpha/2);
  double normy=sqrt(pi*a*a/(beta*2*del*ct*ct)); //a=1 analytic integral over dy 
  double normx=0.8340588041; //numeric integral over dx
  double norm=normx*normy;
  double *peq=new double[M*N];
  double *sqrpeq=new double[M*N];
  for(m=0;m<M;m++){
    for(n=0;n<N;n++){
      x=(m)*hx+x1;
      y=(n)*hy+y1;
      v=calc_bistable_potential(x, y, beta);
      peq[m*N+n]=exp(-beta*v)/norm;
      sqrpeq[m*N+n]=sqrt(peq[m*N+n]);
    }
  }
  ofstream eqfile("Peq.dat");
  ofstream veqfile("Peqvec.dat");
  for(m=0;m<M;m++){
    for(n=0;n<N;n++){
      eqfile<<peq[m*N+n]<<'\t';
      veqfile<<peq[m*N+n]<<endl;
    }
    eqfile<<endl;
  }
  //dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
  char trans='N';//transpose is T
  double ap=1;
  double *yvec=new double[size];
  for(m=0;m<size;m++)
    yvec[m]=0.0;

  int inc=1;
  /*In the symmetric version instead of Apeq=0, A is symmetrized by peq^1/2, so 
   *the zero eigenvalue is from A(twdle)*peq^1/2=0; A(twdle)=peq^-1/2*A*peq^1/2 
   */
  //dgemv(&trans, &size, &size, &ap, A, &size, sqrpeq,&inc, 0, yvec, &inc);
  /*now write out the solution to A*peq, it should be 0!*/
//   ofstream zfile("zeroeig.dat");
//   for(m=0;m<M;m++){
//     for(n=0;n<N;n++){
//       zfile<<yvec[m*N+n]<<'\t';
//     }
//     zfile<<endl;
//  }
  /*Calculate the eigenvalues and vectors for the transition matrix*/
  char jobvl='N';//left eigenvalues are computed for 'V'
  char jobvr='V';//right eigenvalues are computed for 'V'
  int info, lwork;
  lwork=4*size;
  int one=1;
  double *work=new double[lwork];
  double *eigvalue=new double[size];
  double *eigimag=new double[size];
  double *leftvec=new double[size];
  double *rightvec=new double[size*size];
//dgeev_(&jobvl, &jobvr, (__CLPK_integer*)&size, Amat, (__CLPK_integer*)&size, eigvalue, eigimag,leftvec, (__CLPK_integer*)&size, rightvec, (__CLPK_integer*)&size, work, (__CLPK_integer*)&lwork, (__CLPK_integer*)&info);
  dgeev(&jobvl, &jobvr, (MKL_INT*)&size, A, (MKL_INT*)&size, eigvalue, eigimag,leftvec, (MKL_INT*)&one, rightvec, (MKL_INT*)&size, work, (MKL_INT*)&lwork,(MKL_INT*)&info);
  
  /*Now write out the eigenvector, first nonzero value
   *and divide it by peq.
   */
  double *ev=new double[size+1];
  long unsigned int *index=new long unsigned int[size+1];
  
  for(i=0;i<size;i++)
    ev[i+1]=eigvalue[i];

  indexx(size, ev, index);
  for(i=1;i<size+1;i++)
    cout<< "index: "<<index[i]-1<<" eigvalue: "<<ev[index[i]]<<endl;
  
  int large=index[size-1]-1;//largest eigevalue that isn't 1 
  int second=index[size-2]-1;
  printf("Large index, %d second:  %d\n",large,second);
  /*write out the first 50 eigenvectors*/
  int nv=10;
  char fname[80];
  ofstream vfile;
  int start, id;
  double fact2=0;
  double fact;
  double phi;
  for(i=1;i<nv+1;i++){
    
    sprintf(fname, "sym.%d.%lf.dat",i-1,eigvalue[index[i]-1]);
    vfile.open(fname);
    start=index[i]-1;
    /*first normalize the eigenvectors*/
    fact2=0;
    for(m=0;m<M;m++){
      for(n=0;n<N;n++){
	id=m*N+n;
	phi=rightvec[start*size+id]/sqrpeq[id];
	fact2+=phi*peq[id]*phi*hx*hy;
      }
    }
    fact=sqrt(fact2);
    for(m=0;m<M;m++){
      for(n=0;n<N;n++){
	id=m*N+n;
	vfile<<rightvec[start*size+id]/sqrpeq[id]/fact<<'\t';
      }
      vfile<<endl;
    }
    vfile.close();
  }
      
  //for(i=0;i<size;i++)
  //  cout <<eigvalue[i]<<' '<<eigimag[i]<<endl;
  
  /*******************************************************/
  /*now calculate flux*/
//     double xs=-1.5;
//   double xf=1.5;
//   int Nx=200;
//   double dx=(xf-xs)/(1.0*Nx);
//   double ys=-1;
//   double yf=1;
//   int Ny=200;
//   double dy=(yf-ys)/(1.0*Ny);

  // for(i=0;i<Nx;i++){
//     for(j=0;j<Ny;j++){
//       x=xs+i*dx;
//       y=ys+j*dy;
//       v=calc_bistable_potential(x, y, beta);
//       diffusion1(x, y, Dxx, Dyy, beta);
//       Jx=exp(-beta*v)*Dxx;
//       Jy=exp(-beta*v)*Dyy;
//     }
//   }

  
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
  double yc=(y/a)*cot(alpha/2);
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
      ct=cot(alpha/2);
      
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
  ct=cot(alpha/2);
  
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
  ct=cot(alpha/2);
  
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
  ct=cot(alpha/2);
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


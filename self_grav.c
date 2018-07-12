#include "paul.h"

//const int L = 64;
//const int M = 1;
const int N_pp = 16;
const int n = 12;

struct ppoint{  //Point in Poisson Grid

  double x[3];
  double density;
  double potential;
};

struct poisson{

//Poisson Grid
  double L;       //Number of grid points in r direction
  double M;       //Number of grid points in z direction for cylindrical coords
  double N_pp;    //Number of phi grid points at each annulus--can become an array
  //double *r_pj;   //r coords for poisson grid
  //double *z_pk;   //z coords for poisson grid

  struct ppoint **thePoints; //array of points on poisson grid

//FFT Things
  //int n_trunc;  //FFT truncation length/term

//Poisson Equation
  //double **potential;  //Potential at each phi in a given annulus
  //double **density;    //Density at each phi point in an annulus

};

void solve_tri_diag( double *a, double *b, double *c, double *d, double *x, int N ){
  //solves linear systems of the form  a_i*x_i-1 + b_i*x_i + c_i*x_i+1 = d_i

  c[0] = c[0]/b[0];
  d[0] = d[0]/b[0];
//Forward Sweep
  int i;
  for( i=1; i<N; i++ ){
    double m = 1/( b[i]-a[i]*c[i-1] );
    c[i] *= m;
    d[i]  = ( d[i] - a[i]*d[i-1] )*m;
  }
//Back Substitution
  x[N-1] = d[N-1];
  for( i=N-2; i>=0; i-- )
    x[i] = d[i] - c[i]*x[i+1];

  return;
}

void setupPoission( struct domain *theDomain ){

  struct poisson *thePoisson = theDomain->thePoisson;
  thePoisson->L = 64; //theDomain->Nr;
  thePoisson->M = 1;  //theDomain->Nz;

  thePoisson->N_pp = N_pp; //input this from .par file?
  //thePoisson->n_trunc = n;
}

void initializePoisson( struct poisson *thePoisson ){
  //TODO: Put all of this into setupPoisson and then it'll get initialized
  //      in interpolation from disco grid

  int L = thePoisson->L;
  int M = thePoisson->M;
  //N_pp = thePoisson->N_pp;
  //thePoisson->density   = (double **) malloc( L*M*sizeof(double *) );
  //thePoisson->potential = (double **) malloc( L*M*sizeof(double *) );
  thePoisson->thePoints = (struct ppoint **) malloc( L*M*sizeof( struct ppoint *) );


  int jk;
  for( jk=0; jk<L*M; jk++ ){
    //thePoisson->density[ jk ] = (double *) malloc( N_pp*sizeof(double) );
    //thePoisson->potential[jk] = (double *) malloc( N_pp*sizeof(double) );
    thePoisson->thePoints[jk] = (struct ppoint *) malloc( N_pp*sizeof(struct ppoint) );
  }
}

void solve_sg_ode( struct poisson *thePoisson, *fft_buffer ){
  //TODO: Need to include boundary conditions
  //TODO: Need to make this 3D: e.g. include z-direction
  //TODO: ***** fix indexing to account for r2r transform ******

  int L = thePoisson->L;
  int M = thePoisson->M;
  int N_pp = thePoisson->N_pp;

  //Coefficient for solving tri-diag
  double dr = 2./(2*L+1);
  double *a = (double *) malloc( L*sizeof(double) );
  double *b = (double *) malloc( L*sizeof(double) );
  double *c = (double *) malloc( L*sizeof(double) );
  double *F_n = (double *) malloc( L*sizeof(double) );

  double **Phi = (double **) malloc( N_pp*sizeof(double) ); //N_pp pointers

  //Calculate coeffs for U_{i-1} and U_{i+1}--depend only on r_i and dr
  int i;
  double dr2_inv = 1/dr/dr;
  for( i=1; i<L+1; i++ ){
    r_i = dr*(i - 0.5);
    a[i-1] = dr2_inv - 0.5/del_r/r_i;
    c[i-1] = dr2_inv + 0.5/del_r/r_i;
  }

  //Calculate coeffs for U_i and for each n, get array F_n[i] where i is at r_i
  //Then solve tri-diag for each n across all i (annuli)
  int n;
  for( n=0; n<N_pp; n++ ){
    for( i=1; i<L+1; i++ ){
      int idx = n +(i-1)*N_pp;
      r_i = dr*(i - 0.5);
      b[i-1] = -2*dr2_inv - n*n/r_i/r_i;
      F_n[i-1] = fft_buffer[idx];
    }
    //Have phi be pointer to starting point in fft_buffer?
    //Then, solve_tri_diag--should run through only relevant elements of fft_buff?
    Phi[n] = (double *) malloc( L*sizeof(double) );  //an n for every annulus
    solve_tri_diag( a, b, c, F_n, Phi[n], L );
  }

  //Refill real part of fft_buffer with phi coeffs
  swap_fft_buffer( fft_buffer, Phi, L, M, N_pp );

}

freePoisson( struct poisson *thePoisson ){

  int L = thePoisson->L;
  int M = thePoisson->M;
  int N_pp = thePoisson->N_pp;

  int jk;
  for( jk=0; jk<L*M; jk++ )
    free( thePoisson->thePoints[jk] );

  free( thePoisson->thePoints );
}

void self_grav( struct domain *theDomain ){

  struct poisson *thePoisson = theDomain->thePoisson;
  int L = thePoisson->L;
  int M = thePoisson->M;
  int N_pp = thePoisson->N_pp;

  //Interp to Poisson Grid

  //FFT density
  // - build buffer, fft, and then build back?
  double *fft_buffer = (double *) malloc( L*M*N_pp*sizeof(int) );
  build_fft_buffer( thePoisson, fft_buffer );
  density_fft( thePoisson, fft_buffer)

  //Solve Tri-Diag systems
  // - calculate coeffs, then solve
  solve_sg_ode( thePoisson, fft_buffer );

  //Inverse FFT -> potential
  potential_fft( thePoisson, fft_buffer );
  unpack_fft_buffer( thePoisson, fft_buffer );

  //Calc Fore <-> Interpolate back

  free( fft_buffer );
  freePoisson( thePoisson );
  free( thePoisson );
}

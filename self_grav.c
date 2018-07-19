#include "paul.h"


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

void swap_fft_buffer( double *, double **, int, int, int );

void solve_sg_ode( struct poisson *thePoisson, double *fft_buffer ){
  //TODO: Need to include boundary conditions
  //TODO: Need to make this 3D: e.g. include z-direction
  //TODO: ***** fix indexing to account for r2r transform? ******

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
    double r_i = dr*(i - 0.5);
    a[i-1] = dr2_inv - 0.5/dr/r_i;
    c[i-1] = dr2_inv + 0.5/dr/r_i;
  }

  //Calculate coeffs for U_i and for each n, get array F_n[i] where i is at r_i
  //Then solve tri-diag for each n across all i (annuli)
  int n;
  for( n=0; n<N_pp; n++ ){
    for( i=1; i<L+1; i++ ){
      int idx = n +(i-1)*N_pp;
      double r_i = dr*(i - 0.5);
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

void calc_sg_phi_force( struct domain *theDomain ){

  struct cell **theCells = theDomain->theCells;
  int Nr = theDomain->Nr;
  int Nz = theDomain->Nz;
  int * Np = theDomain->Np;

  int i,j,k;
  for( j=0 ; j<Nr ; ++j ){
     for( k=0 ; k<Nz ; ++k ){
        int jk = j+Nr*k;
        for( i=0 ; i<Np[jk] ; ++i ){
           int im = i-1;
           if( im == -1 )
            im = Np[jk]-1;
           int ip = (i+1)%Np[jk];
           struct cell * c  = &(theCells[jk][i]);
           struct cell * cL = &(theCells[jk][im]);
           struct cell * cR = &(theCells[jk][ip]);
           //double psiC = c->Psi;
           double psiL = cL->Psi;
           double psiR = cR->Psi;
           double dpL = cL->dphi;
           double dpC = c->dphi;
           double dpR = cR->dphi;
           c->gradPsi[1] = ( psiR - psiL )/( 0.5*dpL + dpC + 0.5*dpR );
        }
     }
  }

}

double get_signed_dp( double, double );
double get_dA( double *, double *, int );

void calc_sg_trans_force( struct domain *theDomain, struct face *theFaces, int Nf, int dim ){

  struct cell ** theCells = theDomain->theCells;
  int Nr  = theDomain->Nr;
  int Nz  = theDomain->Nz;
  int *Np = theDomain->Np;
  double *r_jph = theDomain->r_jph;
  double *z_kph = theDomain->z_kph;

  int DIR;
  if( dim==1 )
    DIR = 0;
  if( dim==2 )
    DIR = 2;

  //Add weighted slopes
  int n;
  for( n=0 ; n<Nf ; ++n ){
     struct face * f  = &( theFaces[n] );
     double phi = f->cm[1];
     struct cell * cL = f->L;
     struct cell * cR = f->R;
     double dxL = f->dxL;
     double dxR = f->dxR;
     double phiL = cL->piph - .5*cL->dphi;
     double phiR = cR->piph - .5*cR->dphi;
     double dpL = get_signed_dp(phi,phiL);
     double dpR = get_signed_dp(phiR,phi);
     double psiL = cL->Psi + dpL*cL->gradPsi[1];
     double psiR = cR->Psi - dpR*cR->gradPsi[1];

     double dA = f->dA;
     double S  = ( psiR - psiL )/( dxR + dxL );
     cL->gradPsi[DIR] += S*dA;
     cR->gradPsi[DIR] += S*dA;
  }

  //Divide by total weight
  int i,j,k;
  for( j=0 ; j<Nr ; ++j ){
     for( k=0 ; k<Nz ; ++k ){
        int jk = j+Nr*k;
        for( i=0 ; i<Np[jk] ; ++i ){
           struct cell * c = &(theCells[jk][i]);
           double phip = c->piph;
           double phim = phip - c->dphi;
           double xp[3] = {r_jph[j  ],phip,z_kph[k  ]};
           double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
           if( dim==1 )
              xm[0] = r_jph[j];
           else xm[2] = z_kph[k];
           double dAp = get_dA(xp,xm,dim);
           if( dim==1 ){
              xp[0] = r_jph[j-1];
              xm[0] = r_jph[j-1];
           }else{
              xp[2] = z_kph[k-1];
              xm[2] = z_kph[k-1];
           }
           double dAm = get_dA(xp,xm,dim);
           double dAtot = dAp+dAm;
           if( (dim==1 && j==0   ) || (dim==2 && k==0   ) )
              dAtot = dAp;
           if( (dim==1 && j==Nr-1) || (dim==2 && k==Nz-1) )
              dAtot = dAm;
            c->gradPsi[DIR] /= dAtot;
        }
     }
  }

}


void self_grav_src( double *gradPsi, double *prim, double *cons, double *xp, double *xm, double dVdt ){

  double rp = xp[0];
  double rm = xm[0];
  double rho = prim[RHO];
  double vr  = prim[URR];
  double vz  = prim[UZZ];
  double omega = prim[UPP];

  double r = 0.5*(rp+rm);
  double vp = r*omega;

  double Fr = gradPsi[0];
  double Fp = gradPsi[1];
  double Fz = gradPsi[2];

  cons[SRR] += rho*Fr*dVdt;
  cons[SZZ] += rho*Fz*dVdt;
  cons[LLL] += rho*Fp*r*dVdt;
  cons[TAU] += rho*( Fr*vr + Fz*vz + Fp*vp )*dVdt;
}

void freePoisson( struct poisson *thePoisson ){

  int L = thePoisson->L;
  int M = thePoisson->M;
  //int N_pp = thePoisson->N_pp;

  int jk;
  for( jk=0; jk<L*M; jk++ )
    free( thePoisson->thePoints[jk] );

  free( thePoisson->thePoints );
}

void initializePoisson( struct domain * );
void build_fft_buffer( struct poisson *, double * );
void density_fft( struct poisson *, double * );
void potential_fft( struct poisson *, double * );
void unpack_fft_buffer( struct poisson *, double * );
void finalizePoisson( struct domain * );

void self_grav( struct domain *theDomain ){

  struct poisson *thePoisson = theDomain->thePoisson;
  int L = thePoisson->L;
  int M = thePoisson->M;
  int N_pp = thePoisson->N_pp;

  //Interp to Poisson Grid
  initializePoisson( theDomain );

  //FFT density
  // - build buffer, fft, and then build back?
  double *fft_buffer = (double *) malloc( L*M*N_pp*sizeof(int) );
  build_fft_buffer( thePoisson, fft_buffer );
  density_fft( thePoisson, fft_buffer);

  //Solve Tri-Diag systems
  // - calculate coeffs, then solve
  solve_sg_ode( thePoisson, fft_buffer );

  //Inverse FFT -> potential
  potential_fft( thePoisson, fft_buffer );
  unpack_fft_buffer( thePoisson, fft_buffer );

  //Interpolate back
  finalizePoisson( theDomain );

  //Calculate gravitational force
  int Nf;
  struct face *theFaces;
  calc_sg_phi_force( theDomain );
  Nf = theDomain->fIndex_r[theDomain->N_ftracks_r];
  theFaces = theDomain->theFaces_1;
  calc_sg_trans_force( theDomain, theFaces, Nf, 1 );
  Nf = theDomain->fIndex_z[theDomain->N_ftracks_z];
  theFaces = theDomain->theFaces_2;
  calc_sg_trans_force( theDomain, theFaces, Nf, 2 );

  free( fft_buffer );
  freePoisson( thePoisson );
  free( thePoisson );
}

/*
void get_cell( struct cell *annulus, double Np, double phi, double phi_max, struct cell *c ){

  int i;
  for( i=0; i<Np; i++ ){
    struct cell *temp = &(annulus[i]);
    if( check_phi(phi, temp->piph, temp->dphi, phi_max) ){
      c = temp;
      return;
    }
  }
  c = NULL;
  printf(" Error: phi search for poisson grid didn't work; will SegFault");

}
*/

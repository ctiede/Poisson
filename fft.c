#include "paul.h"

void build_fft_buffer( struct poisson *thePoisson, double *fft_buffer ){

  int L = thePoisson->L;
  int M = thePoisson->M;
  int N_pp = thePoisson->N_pp;
  struct ppoint **thePoints = thePoints->thePoints;

  int i,j,k;
  for( j=0; j<L; j++ ){
    for( k=0; k<M; k++ ){
      int jk = j + L*k;
      for( i=0; i<N_pp; i++ ){
        struct ppoint *pt = &(thePoints[jk][i]);
        int idx = jk*N_pp + i;
        fft_buffer[idx] = 4*pi*pt->density;        //TODO: Double check indexing
      }
    }
  }
/*
  int i,j;
  for( j=0; j<L; j++ ){
    for( i=0; i<N_pp; i++ ){
      struct ppoint *pt = &(thePoints[j][i]);
      int idx = j + i;
      fft_buffer[idx] = 4*pi*pt->density;
    }
  }
*/

}

void unpack_fft_buffer( struct poisson *thePoisson, double *fft_buffer ){

  int L = thePoisson->L;
  int M = thePoisson->M;
  int N_pp = thePoisson->N_pp;
  struct ppoint **thePoints = thePoints->thePoints;

  int i,j,k;
  for( j=0; j<L; j++ ){
    for( k=0; k<M; k++ ){
      int jk = j + L*k;
      for( i=0; i<N_pp; i++ ){
        struct ppoint *pt = &(thePoints[jk][i]);
        int idx = jk*N_pp + i;
        pt->potential = fft_buffer[idx];
      }
    }
  }
}



void swap_fft_buffer( double *fft_buffer, double **Phi, int L, int M, int N_pp ){
  //rn only 2D
  //TODO: eventually change to acoomodate 3D algorithm
  //TODO: double check indexing

  int i,n;
  for( i=0; i<L; i++ ){
    for( n=0; n<N_pp; n++){
      idx = i*N_pp + n;
      fft_buffer[idx] = Phi[n][i];
    }
  }

}


void density_fft( struct poisson *thePoisson, double *fft_buffer ){

  int L = thePoisson->L;
  int M = thePoisson->M;
  int N_pp = thePoisson->N_pp;

  //fft_buffer = (double *) malloc( Nr*Nz*N_pp*sizeof(int) );
  //struct ppoint *thePoints = thePoisson->thePoints;
  //build_fft_buffer( fft_buffer, thePoints, L, M, N_pp );

  dim = 1;
	const int n[1] = { N_pp };
	const fftw_r2r_kind kind[1] = { FFTW_R2HC }; //Real to Half-Complex
	fftw_plan my_plan = fftw_plan_many_r2r( dim, n, L*M, fft_buffer, n, 1, N_pp,
                                          fft_buffer, n, 1, N_pp, kind, FFTW_ESTIMATE );
	//fftw_plan_many_r2r( dim, size, # of transforms (annuli), input, NULL or n, stride=1, # of elements
 	//                    output, NULL or n, stride=1, # of elements, kind, flags)
  fftw_execute( my_plan );
	fftw_destroy_plan( my_plan );
	fftw_cleanup();

}

void potential_fft( struct poisson *thePoisson, double *fft_buffer ){

  int L = thePoisson->L;
  int M = thePoisson->M;
  int N_pp = thePoisson->N_pp;

  dim = 1;
  const int n[1] = { N_pp };
  const fftw_r2r_kind kind[1] = { FFTW_HC2R };  //Half-Complex to Real
  fftw_plan my_plan = fftw_plan_many_r2r( dim, n, L*M, fft_buffer, n, 1, N_pp,
                                          fft_buffer, n, 1, N_pp, kind, FFTW_ESTIMATE );
	//fftw_plan_many_r2r( dim, size, # of transforms (annuli), input, NULL or n, stride=1, # of elements
 	//                    output, NULL or n, stride=1, # of elements, kind, flags)
  fftw_execute( my_plan );
	fftw_destroy_plan( my_plan );
	fftw_cleanup();

  //Need to divide by number of elements for inverse fft b/c fftw doesn't
  int idx;
  for( idx=0; i<L*M*N_pp; i++ )
    fft_buffer /= N_pp;

}

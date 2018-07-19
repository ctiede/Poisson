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

  struct ppoint **thePoints; //array of points on poisson grid

};

void setupPoission( struct domain *theDomain ){

  struct poisson *thePoisson = theDomain->thePoisson;
  thePoisson->L = theDomain->Nr; //64;
  thePoisson->M = 1;  //theDomain->Nz;
  thePoisson->N_pp = N_pp; //input this from .par file?

  thPoisson->thePoints = (struct ppoint **) malloc( L*M*sizeof(sturct ppoint *) );
  int jk;
  for( jk=0; jk<L*M; jk++ ){
    //thePoisson->density[ jk ] = (double *) malloc( N_pp*sizeof(double) );
    //thePoisson->potential[jk] = (double *) malloc( N_pp*sizeof(double) );
    thePoisson->thePoints[jk] = (struct ppoint *) malloc( N_pp*sizeof(struct ppoint) );
  }
}

/*
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
*/
int check_phi( double, double, double, double );

void interpolate( struct ppoint *points, double N_pp, struct cell *cells, double Np, double phi_max, int dir ){

  int i = 0;
  double x_run = 0.0;
  double offset = 0.0;
  double dpp = phi_max/N_pp;

  if( dir==0 ){  //interpolating density from disco grid to poisson grid
    int p;       //p for poisson grid
    for( p=0; p<N_pp; p++ ){         //Loop through all points in annulus of poisson grid
      struct ppoint *pt = points+p;  // &(points[p]);
      pt->density = 0.0;
      double x_p = (p+1)*dpp;
      //run through disco cells that overlap with poisson "cell"
      while( x_run < x_p ){
        if( i > Np )
          printf("Error: overran disco array; will SegFault");

        struct cell *c = cells+i;  // &(cells[i]);
        double crho = c->prim[RHO];
        double dx = c->dphi - offset;
        if( x_run + dx > x_p ){
          dx  = x_p - x_run; //should be last iteration
          offset += dx;
        }
        else{
          offset = 0.0;
          i++;
        }
        pt->density += (dx/dpp)*crho;
        c->Psi = 0.0;  //reset potential for later
        x_run += dx;
      }
    }
  }
  else if( dir==1 ){  //interpolating potential from poisson grid to disco grid
    int p;
    for( p=0; p<N_pp; p++ ){
      struct ppoint *pt = points+p;  // &(points[p]);
      double Psi = pt->potential;
      double x_p = (p+1)*dpp;
      while( x_run < x_p ){
        if( i > Np )
          printf("Error: overran disco array; will SegFault");

        struct cell *c = cells+i;  // &(cells[i]);
        double dx = c->dphi - offset;
        if( x_run + dx > x_p ){
          dx  = x_p - x_run;
          offset += dx;
        }
        else{
          offset = 0.0;
          i++;
        }
        c->Psi += (dx/c->dphi)*Psi;
        x_run += dx;
      }
    }
  }
  else
    printf("\nError: illegal direction for interpolation\n");


/*
  else{     //interpolating from poisson grid back to disco grid
    int d;  //d for disco grid
    struct ppoint *pt = points;
    for( d=0; d<Np; d++ ){
      struct cell *c = cells+d; //&(cells[d])
      double dx = c->dphi;
      flag = 1;
      c->Psi = 0.0;
      double xp = (i+1)*dpp;
      if( x_run + dx > xp ){
        dx1 = xp - x_run;
        dx2 =
      }
      else
        c->Psi = pt->potential

      while( flag ){
        if( i > N_pp )
          //TODO: figure out what to do here
        struct ppoint *pt = points+i; // &(points[i]);
        double xp = (i+1)*dpp;
        if( x_run + dx >= xp ){
          dx = xp - x_run;
          i++;
          flag = 1;
        }
        else{
          offset = 0.0;
          flag = 0;
          //c->Psi = pt->potential;
        }
        c->Psi += (dx/c->dphi)*pt->potential;
        x_run += dx;
      }
    }
  }
*/
}

void initializePoisson( struct domain *theDomain ){
  //Sets values of density for poisson grid points by interpolating from
  //the Disco grid

  //TODO: don't need to keep the points ->can just build buffer right away
  //
  //TODO: Double check sign on dx.

  struct poisson *thePoisson = theDomain->thePoisson;
  struct ppoint  **thePoints = thePoisson->thePoints;
  int L = thePoisson->L;
  int M = thePoisson->M;
  int N_pp = thePoisson->N_pp;

  struct cell **theCells = theDomain->theCells;
  int *Np = theDomain->Np;
  double phi_max = theDomain->phi_max;

  int j, k;
  for( j=0; j<L; j++ ){
    for( k=0; k<M; k++ ){
      int jk = j +k*L;
      interpolate( thePoints[jk], N_pp, theCells[jk], Np[jk], phi_max, 0 );
    }
  }

}

void finalizePoisson( struct domain *theDomain ){

  struct poisson *thePoisson = theDomain->thePoisson;
  struct ppoint **thePoints = thePoisson->thePoints;
  int L = thePoisson->L;
  int M = thePoisson->M;
  int N_pp = thePoisson->N_pp;

  struct cell **theCells = theDomain->theCells;
  int *Np = theDomain->Np;
  double phi_max = theDomain->phi_max;

  int j, k;
  for( j=0; j<L; j++ ){
    for( k=0; k<M; k++ ){
      int jk = j +k*L;
      interpolate( thePoints[jk], N_pp, theCells[jk], Np[jk], phi_max, 1 );
    }
  }
}

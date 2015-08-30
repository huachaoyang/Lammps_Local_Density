/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Tanmoy Sanyal, M.Scott Shell, UC Santa Barbara
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_localdensity.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "domain.h" // HERE

using namespace LAMMPS_NS;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

PairLOCALDENSITY::PairLOCALDENSITY(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  one_coeff = 1; 	
  single_enable = 1;
  DEBUG = 1; // turn this off if display of parsed and splined arrays is not required	

  // read from file
  nLD = 0;
  nrho = 0;
  rho_min = NULL;
  rho_max = NULL;
  a = NULL;
  b = NULL;
  c0 = NULL;
  c2 = NULL;
  c4 = NULL;
  c6 = NULL;
  uppercut = NULL;
  lowercut = NULL;
  uppercutsq = NULL;
  lowercutsq = NULL;
  frho = NULL;
  rho = NULL;
  
  // splined arrays
  frho_spline = NULL;
  
  // per-atom arrays
  nmax = 0;
  fp = NULL;
  localrho = NULL;  

  // set comm size needed by this Pair

  comm_forward = 1;
  comm_reverse = 1;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairLOCALDENSITY::~PairLOCALDENSITY()
{

  memory->destroy(localrho);
  memory->destroy(fp);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }

  memory->destroy(frho_spline);
  
  memory->destroy(rho_min);  
  memory->destroy(rho_max);
  memory->destroy(delta_rho);	
  memory->destroy(c0);
  memory->destroy(c2);
  memory->destroy(c4);
  memory->destroy(c6);
  memory->destroy(uppercut);
  memory->destroy(lowercut);
  memory->destroy(uppercutsq);
  memory->destroy(lowercutsq);
  memory->destroy(frho);
  memory->destroy(rho);

  delete [] a;
  delete [] b;
  
  
}

/* ---------------------------------------------------------------------- */

void PairLOCALDENSITY::compute(int eflag, int vflag)
{
  
  int i,j,ii,jj,m,k,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double rsqinv, phi, uLD, dphi, evdwl,fpair;
  double p, *coeff;
  int *ilist,*jlist,*numneigh,**firstneigh;

  phi = uLD = evdwl = fpair = rsqinv = 0.0;		

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = eflag_global = eflag_atom = 0;

  // grow local-density and fp arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(localrho);
    memory->destroy(fp);
    nmax = atom->nmax; 
    memory->create(localrho, nLD, nmax, "pairLD:localrho");
    memory->create(fp, nLD, nmax, "pairLD:fp");
  }

  double **x = atom->x; 
  double **f = atom->f;
  int *type = atom->type; 
  int nlocal = atom->nlocal; 
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero out density and fp

  if (newton_pair) {
    m = nlocal + atom->nghost;
    for (k = 0; k < nLD; k++) { 
        for (i = 0; i < m; i++)	{ 		
            localrho[k][i] = 0.0;
            fp[k][i] = 0.0;
        }
    }	
  } 
  else {
    for (k = 0; k < nLD; k++){
        for (i = 0; i < nlocal; i++) {
            localrho[k][i] = 0.0;
            fp[k][i] = 0.0;
        }
    }
   }


  // localrho = local density at each atom
  // loop over neighbors of my atoms and types of local-densities
  
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
   
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];	

      // calculate distance-squared between i,j atom-types
      
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;	
      
      // calculating local densities based on central and neighbor filters

      for (k = 0; k < nLD; k++) {
        if (rsq < lowercutsq[k]) {
             phi = 1.0;
        }
          else if (rsq > uppercutsq[k]) {
             phi = 0.0;
        }
          else {
             phi = c0[k] + rsq * (c2[k] + rsq * (c4[k] + c6[k]*rsq));
        }
        localrho[k][i] += (phi * b[k][jtype]); 
        //checking for both i,j is necessary since a half neighbor list is processed.
        if (newton_pair || j<nlocal) {
            localrho[k][j] += (phi * b[k][itype]); 
        }
      }
    }	       
  }

  // communicate and sum densities
  if (newton_pair) comm->reverse_comm_pair(this);

  // uLD = embedding energy of each atom due to each LD potential type 
  // fp = derivative of embedding energy at each atom for each LD potential type

  for (ii = 0; ii < inum; ii++) {	
    i = ilist[ii];
    itype = type[i];
    uLD = 0.0;	

    for (k = 0; k < nLD; k++) {

        // skip over this loop if the LD potential is not for itype
        if (!(a[k][itype])) continue; 
            
        // linear extrapolation at rho_min and rho_max
            
        if (localrho[k][i] <= rho_min[k]) {
            coeff = frho_spline[k][0];
            fp[k][i] = coeff[2];
            uLD +=  a[k][itype] * ( coeff[6] + fp[k][i]*(localrho[k][i] - rho_min[k]) );
        }
        else if (localrho[k][i] >= rho_max[k]) {
            coeff = frho_spline[k][nrho-2];
            fp[k][i] = coeff[0] + coeff[1]  + coeff[2];
            uLD +=  a[k][itype] * ( (coeff[3] + coeff[4] + coeff[5] + coeff[6]) + fp[k][i]*(localrho[k][i] - rho_max[k]) );
        }
        else {
            p = (localrho[k][i] - rho_min[k]) / delta_rho[k];
            m = static_cast<int> (p);
            m = MAX(0, MIN(m, nrho-2));
            p -= m;
            p = MIN(p, 1.0);
            coeff = frho_spline[k][m];
            fp[k][i] = (coeff[0]*p + coeff[1])*p + coeff[2];
            uLD +=  a[k][itype] * (((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6]);
        }
    }

    if (eflag) {
        if (eflag_global) eng_vdwl += uLD;
        if (eflag_atom) eatom[i] += uLD;	
    }
 }

  // communicate derivatives of embedding function and localdensity

  comm->forward_comm_pair(this);
  
  // compute forces on each atom
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];	
      j &= NEIGHMASK;
      jtype = type[j];

      // calculate square of distance between i,j atoms

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
                
      // calculate force between two atoms  
      fpair = 0.0;
      if (rsq < cutforcesq) {   // global cutoff check
        rsqinv = 1.0/rsq;
        for (k = 0; k < nLD; k++) {
            if (rsq >= lowercutsq[k] && rsq < uppercutsq[k]) {
               dphi = rsq * (2.0*c2[k] + rsq * (4.0*c4[k] + 6.0*c6[k]*rsq));
               fpair += -(a[k][itype]*b[k][jtype]*fp[k][i] + a[k][jtype]*b[k][itype]*fp[k][j]) * dphi; 
            }
        }	
        fpair *= rsqinv; 
        
        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
            f[j][0] -= delx*fpair;
            f[j][1] -= dely*fpair;
            f[j][2] -= delz*fpair;
        }
      
      if (eflag) evdwl = 0.0; // eng_vdwl has already been completely built, so no need to add anything here
      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }

    }
  }
  
  if (vflag_fdotr) virial_fdotr_compute();
}



/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLOCALDENSITY::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
 
  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLOCALDENSITY::settings(int narg, char **arg)
{
  if (narg > 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for all type pairs
   read LD file
------------------------------------------------------------------------- */

void PairLOCALDENSITY::coeff(int narg, char **arg)
{
  int i, j;
  if (!allocated) allocate();

  if (narg != 3) error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // parse LD file

  parse_file(arg[2]);


 // clear setflag since coeff() called once with I,J = * *

  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++)
      setflag[i][j] = 0;

  // set setflag for all i,j type pairs

  int count = 0;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
        setflag[i][j] = 1;
        count++;
      }
    }
  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLOCALDENSITY::init_style()
{
  // spline rho and frho arrays
  // request half neighbor list

  array2spline();
  if (DEBUG)
    display();

  // half neighbor request
  neighbor->request(this);

  /* //full neighbor list (needed ?)
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1; */
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLOCALDENSITY::init_one(int i, int j)
{
  // single global cutoff = max of all uppercuts read in from LD file

  cutmax = 0.0;
  for (int k = 0; k < nLD; k++)
    cutmax = MAX(cutmax,uppercut[k]);
    
  cutforcesq = cutmax*cutmax;

  return cutmax;
}


/*--------------------------------------------------------------------------
  pair_write functionality for this pair style that gives just a snap-shot 
  of the LD potential without doing an actual MD run
 ---------------------------------------------------------------------------*/

double PairLOCALDENSITY::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{
    int m, k, index;
    double rsqinv, p, uLD;
    double *coeff, **LD;
    double dFdrho, phi, dphi;

    uLD = dFdrho = dphi = 0.0;

    memory->create(LD, nLD, 3, "pairLD:LD");
    for (k = 0; k < nLD; k++) {
        LD[k][1] = 0.0; // itype:- 1
        LD[k][2] = 0.0; // jtype:- 2
        }	

    rsqinv = 1.0/rsq;
    for (k = 0; k < nLD; k++) {
        if (rsq < lowercutsq[k]) {
             phi = 1.0;
        }
        else if (rsq > uppercutsq[k]) {
             phi = 0.0;
        }
        else {
             phi = c0[k] + rsq * (c2[k] + rsq * (c4[k] + c6[k]*rsq));
        }
        LD[k][1] += (phi * b[k][jtype]);
        LD[k][2] += (phi * b[k][itype]);           
    }

    for (k = 0; k < nLD; k++) {
        if (a[k][itype]) index = 1;
        if (a[k][jtype]) index = 2;
        
        if (LD[k][index] <= rho_min[k]) {
            coeff = frho_spline[k][0];
            dFdrho = coeff[2];
            uLD +=  a[k][itype] * ( coeff[6] + dFdrho*(LD[k][index] - rho_min[k]) );
        }
        else if (LD[k][index] >= rho_max[k]) {
            coeff = frho_spline[k][nrho-1];
            dFdrho = coeff[0] + coeff[1]  + coeff[2];
            uLD +=  a[k][itype] * ( (coeff[3] + coeff[4] + coeff[5] + coeff[6]) + dFdrho*(LD[k][index] - rho_max[k]) );
        }
        else {
            p = (LD[k][index] - rho_min[k]) / delta_rho[k];
            m = static_cast<int> (p);
            m = MAX(0, MIN(m, nrho-2));
            p -= m;
            p = MIN(p, 1.0);
            coeff = frho_spline[k][m];
            dFdrho = (coeff[0]*p + coeff[1])*p + coeff[2];
            uLD +=  a[k][itype] * (((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6]);
        }

        if (rsq < lowercutsq[k]) {
           dphi = 0.0;
        }
        else if (rsq > uppercutsq[k]) {
           dphi = 0.0;
        }
        else {
           dphi = rsq * (2.0*c2[k] + rsq * (4.0*c4[k] + 6.0*c6[k]*rsq));
        }
        fforce +=  -(a[k][itype]*b[k][jtype]*dFdrho + a[k][jtype]*b[k][itype]*dFdrho) * dphi *rsqinv;
    }
    memory->destroy(LD);

    return uLD;
}

/*--------------------------------------------------------------------
   Spline the array frho read in from the file to create
   frho_spline
---------------------------------------------------------------------- */

void PairLOCALDENSITY::array2spline() {


  memory->destroy(frho_spline);
  memory->create(frho_spline, nLD, nrho, 7, "pairLD:frho_spline");

  for (int k = 0; k < nLD; k++)
    interpolate_clamped(nrho, delta_rho[k], frho[k], frho_spline[k]);

}


/* ---------------------------------------------------------------------- 
  Cubic Spline interpolation sub-routines
 ------------------------------------------------------------------------*/

void PairLOCALDENSITY::interpolate_clamped(int n, double delta, double *f, double **spline) {
/*   TAKEN FROM OTHER CODE:
     PURPOSE:
          determine the coefficients for the clamped
          cubic spline for a given set of data


     CALLING SEQUENCE:
          cubic_clamped ( n, x, f, b, c, d, fpa, fpb );


     INPUTS:
          n		number of interpolating points
          x		array containing interpolating points
          f		array containing function values to
			be interpolated;  f[i] is the function
			value corresponding to x[i]
          b		array of size at least n; contents will
			be overwritten
          c		array of size at least n; contents will
			be overwritten
          d		array of size at least n; contents will
			be overwritten
          fpa		derivative of function at x=a
          fpb		derivative of function at x=b


     OUTPUTS:
          b		coefficients of linear terms in cubic 
			spline
	  c		coefficients of quadratic terms in
			cubic spline
	  d		coefficients of cubic terms in cubic
			spline

     REMARK:
          remember that the constant terms in the cubic spline
          are given by the function values being interpolated;
          i.e., the contents of the f array are the constant
          terms

          to evaluate the cubic spline, use the routine
          'spline_eval'
*/
                     
     double *dl, *dd, *du;
     double *b, *c, *d;
     double fpa, fpb;

     int i;
     
     b = new double [n];
     c = new double [n];
     d = new double [n];
     dl = new double [n];
     dd = new double [n];
     du = new double [n];

     // initialize values
     for ( i = 0; i<n; i++) {
         b[i] = c[i] = d[i] = 0.;
         dl[i] = dd[i] = du[i] = 0.;
     }

     // set slopes at beginning and end
     fpa = 0.;
     fpb = 0.;
     
     for ( i = 0; i < n-1; i++ ) {
         dl[i] = du[i] = delta;
     }
     
     dd[0] = 2.0 * delta;
     dd[n-1] = 2.0 * delta;
     c[0] = ( 3.0 / delta ) * ( f[1] - f[0] ) - 3.0 * fpa;
     c[n-1] = 3.0 * fpb - ( 3.0 / delta ) * ( f[n-1] - f[n-2] );
     for ( i = 0; i < n-2; i++ ) {
         dd[i+1] = 4.0 * delta;
         c[i+1] = ( 3.0 / delta ) * ( f[i+2] - f[i+1] ) -
                  ( 3.0 / delta ) * ( f[i+1] - f[i] );
     }
     
     // tridiagonal solver
     for ( i = 0; i < n-1; i++ ) {
         du[i] /= dd[i];
         dd[i+1] -= dl[i]*du[i];
     }
     
     c[0] /= dd[0];
     for ( i = 1; i < n; i++ )
         c[i] = ( c[i] - dl[i-1] * c[i-1] ) / dd[i];
         
     for ( i = n-2; i >= 0; i-- )
         c[i] -= c[i+1] * du[i];
     
     for ( i = 0; i < n-1; i++ ) {
         d[i] = ( c[i+1] - c[i] ) / ( 3.0 * delta );
         b[i] = ( f[i+1] - f[i] ) / delta - delta * ( c[i+1] + 2.0*c[i] ) / 3.0;
     }

     // normalize
     for ( i = 0; i < n-1; i++ ) {
         b[i] = b[i] * delta ;
         c[i] = c[i] * delta*delta ;
         d[i] = d[i] * delta*delta*delta;
     }

     //copy to coefficient matrix
     for ( i = 0; i < n; i++) {
         spline[i][3] = d[i];
         spline[i][4] = c[i];
         spline[i][5] = b[i];
         spline[i][6] = f[i];
         spline[i][2] = spline[i][5]/delta;
         spline[i][1] = 2.0*spline[i][4]/delta;
         spline[i][0] = 3.0*spline[i][3]/delta;
     }
     
     delete [] b;
     delete [] c;
     delete [] d;
     delete [] du;
     delete [] dd;
     delete [] dl;
}


void PairLOCALDENSITY::interpolate_natural(int n, double delta, double *f, double **spline) {
  for (int m = 0; m < n; m++) 
      spline[m][6] = f[m];

  spline[0][5] = spline[1][6] - spline[0][6];
  spline[1][5] = 0.5 * (spline[2][6] - spline[0][6]);
  spline[n-2][5] = 0.5 * (spline[n-1][6] - spline[n-3][6]);
  spline[n-1][5] = spline[n-1][6] - spline[n-2][6];

  for (int m = 2; m < n-2; m++)
    spline[m][5] = ((spline[m-2][6] - spline[m+2][6]) +
                    8.0*(spline[m+1][6] - spline[m-1][6])) / 12.0;

  for (int m = 0; m < n-1; m++) {
    spline[m][4] = 3.0*(spline[m+1][6]-spline[m][6]) -
      2.0*spline[m][5] - spline[m+1][5];
    spline[m][3] = spline[m][5] + spline[m+1][5] -
      2.0*(spline[m+1][6] - spline[m][6]);
  }

  spline[n-1][4] = 0.0;
  spline[n-1][3] = 0.0;

  for (int m = 0; m < n; m++) {
    spline[m][2] = spline[m][5]/delta;
    spline[m][1] = 2.0*spline[m][4]/delta;
    spline[m][0] = 3.0*spline[m][3]/delta;
  }
}



/* ----------------------------------------------------------------------
   read potential values from a single Local Density file
------------------------------------------------------------------------- */

void PairLOCALDENSITY::parse_file(char *filename) {
    
  int k, n;
  int me = comm->me;
  FILE *fptr;
  char line[MAXLINE];
  double ratio, lc2, uc2, denom;


  if (me == 0) {
    fptr = fopen(filename, "r");
    if (fptr == NULL) {
      char str[128];
      sprintf(str,"Cannot open Local Density potential file %s",filename);
      error->one(FLERR,str);
    }
  }


 double *ftmp; // temprary variable to extract the complete 2D frho array from file
   
 // broadcast number of LD potentials and number of (rho,frho) pairs
 if (me == 0) {
    
    // first 2 comment lines ignored	
    fgets(line,MAXLINE,fptr);
    fgets(line,MAXLINE,fptr);
    
    // extract number of potentials and number of (frho, rho) points
    fgets(line,MAXLINE,fptr);	
    sscanf(line, "%d %d", &nLD, &nrho);
    fgets(line,MAXLINE,fptr);
  }

  MPI_Bcast(&nLD,1,MPI_INT,0,world);
  MPI_Bcast(&nrho,1,MPI_INT,0,world);
  
  // setting up all arrays to be read from files and broadcasted
  memory->create(uppercut, nLD, "pairLD:uppercut");
  memory->create(lowercut, nLD, "pairLD:lowercut");
  memory->create(uppercutsq, nLD, "pairLD:uppercutsq");
  memory->create(lowercutsq, nLD, "pairLD:lowercutsq");
  memory->create(c0, nLD, "pairLD:c0");
  memory->create(c2, nLD, "pairLD:c2");
  memory->create(c4, nLD, "pairLD:c4");
  memory->create(c6, nLD, "pairLD:c6");
  memory->create(rho_min,  nLD, "pairLD:rho_min"); 
  memory->create(rho_max,  nLD, "pairLD:rho_max");
  memory->create(delta_rho, nLD,"pairLD:delta_rho");
  memory->create(ftmp, nrho*nLD, "pairLD:ftmp");
  
  // setting up central and neighbor atom filters		
  memory->create(a, nLD, atom->ntypes+1 , "pairLD:a");
  memory->create(b, nLD, atom->ntypes+1, "pairLD:b"); 	
  if (me == 0) {
    for (n = 1; n <= atom->ntypes; n++){
        for (k = 0; k < nLD; k++) {
            a[k][n] = 0;
            b[k][n] = 0;
        }
    }
  }	
  
 // read file block by block
  
  if (me == 0) {
    for (k = 0; k < nLD; k++) {
    
        // parse upper and lower cut values	
        if (fgets(line,MAXLINE,fptr)==NULL) break;
        sscanf(line, "%lf %lf", &lowercut[k], &uppercut[k]);
    
        // parse and broadcast central atom filter
        fgets(line, MAXLINE, fptr);
        char *tmp = strtok(line, " /t/n/r/f");
        while (tmp != NULL) {
            a[k][atoi(tmp)] = 1;
            tmp = strtok(NULL, " /t/n/r/f");
        }
        
        // parse neighbor atom filter
        fgets(line, MAXLINE, fptr);
        tmp = strtok(line, " /t/n/r/f");
        while (tmp != NULL) {			
            b[k][atoi(tmp)] = 1;
            tmp = strtok(NULL, " /t/n/r/f");
        }
    
        // parse min, max and delta rho values
        fgets(line, MAXLINE, fptr);
        sscanf(line, "%lf %lf %lf", &rho_min[k], &rho_max[k], &delta_rho[k]);
        // recompute delta_rho from scratch for precision
        delta_rho[k] = (rho_max[k] - rho_min[k]) / (nrho - 1);
        
        // parse tabulated frho values from each line into temporary array
        for (n = 0; n < nrho; n++) {  
            fgets(line,MAXLINE,fptr);
            sscanf(line, "%lf", &ftmp[k*nrho + n]);
        }
        
        // ignore blank line at the end of every block
        fgets(line,MAXLINE,fptr);

        // set coefficients for local density indicator function
        uc2 = uppercut[k] * uppercut[k];
        uppercutsq[k] = uc2;
        lc2 = lowercut[k] * lowercut[k];
        lowercutsq[k] = lc2;
        ratio = lc2/uc2;
        denom = 1.0 - ratio;
        denom = denom*denom*denom;
        c0[k] = (1 - 3.0 * ratio) / denom;
        c2[k] = (6.0 * ratio) / (uc2 * denom);
        c4[k] = -(3.0 + 3.0*ratio) / (uc2*uc2 * denom);
        c6[k] = 2.0 / (uc2*uc2*uc2 * denom);
      }
  }

  // Broadcast all parsed arrays	
  MPI_Bcast(&lowercut[0], nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&uppercut[0], nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&lowercutsq[0], nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&uppercutsq[0], nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&c0[0], nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&c2[0], nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&c4[0], nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&c6[0], nLD, MPI_DOUBLE, 0, world);
  for (k = 0; k < nLD; k++) {
      MPI_Bcast(&a[k][1], atom->ntypes, MPI_INT, 0, world);
      MPI_Bcast(&b[k][1], atom->ntypes, MPI_INT, 0, world);
  }
  MPI_Bcast(&rho_min[0],  nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&rho_max[0],  nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&delta_rho[0], nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&ftmp[0], nLD*nrho, MPI_DOUBLE, 0, world);

  if (me == 0) fclose(fptr);

  // set up rho and frho arrays
  memory->create(rho, nLD, nrho, "pairLD:rho");
  memory->create(frho, nLD, nrho, "pairLD:frho"); 
  
  for (k = 0; k < nLD; k++) {
    for (n = 0; n < nrho; n++) {
        rho[k][n] = rho_min[k] + n*delta_rho[k];
        frho[k][n] = ftmp[k*nrho + n];
    }
 }

  // delete temporary array
  memory->destroy(ftmp);

}
 

/* ----------------------------------------------------------------------
   communication routines
------------------------------------------------------------------------- */


int PairLOCALDENSITY::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc) {
  int i,j,k;
  int m; 	

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i]; 
    for (k = 0; k < nLD; k++) {
      buf[m++] = fp[k][j];  
    }		
  }
  
  return nLD;
}

/* ---------------------------------------------------------------------- */

void PairLOCALDENSITY::unpack_comm(int n, int first, double *buf) {

  int i,k,m,last;
  
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    for (k = 0; k < nLD; k++) {
      fp[k][i] = buf[m++];
    }
 }		
}

/* ---------------------------------------------------------------------- */

int PairLOCALDENSITY::pack_reverse_comm(int n, int first, double *buf) {

  int i,k,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    for (k = 0; k < nLD; k++) {
      buf[m++] = localrho[k][i];
    }
  }
  return nLD;
}

/* ---------------------------------------------------------------------- */

void PairLOCALDENSITY::unpack_reverse_comm(int n, int *list, double *buf) {

  int i,j,k;
  int m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    for (k = 0; k < nLD; k++) {
      localrho[k][j] += buf[m++];
    }	
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairLOCALDENSITY::memory_usage()
{
  double bytes = maxeatom * sizeof(double);
  bytes += maxvatom*6 * sizeof(double);
  bytes += 2 * (nmax*nLD) * sizeof(double);
  return bytes;
}


/* ----------------------------------------------------------------------
   displaying parsed and splined data (for debugging)
------------------------------------------------------------------------- */

void PairLOCALDENSITY::display(){
 int i,j,k, m;
 int me = comm->me;
 double p, LD, uLD, duLD, dduLD, delta, *coeff;

 if (me == 0) {
    FILE *parselog = fopen("parselog.txt", "w");
    if (parselog == NULL) {
        printf("Error opening file!\n");
            exit(1);
    } 
    fprintf(parselog, "#LD-FILE PARSED, FRHO ARRAY SPLINED\n");
    fprintf(parselog, "#DATA FLOW TILL NOW (FOR DEBUGGING PURPOSES)\n");
    fprintf(parselog, "\n\n -------------------------------------------\n\n");
    fprintf(parselog, "N_LD = %d\n", nLD);
    fprintf(parselog, "N_RHO = %d\n", nrho);
 
    fprintf(parselog, "UPPERCUT: ");	
    for (k = 0; k < nLD; k++)
        fprintf(parselog, "%lf\t", uppercut[k]);
        fprintf(parselog, "\nLOWERCUT: ");	
    for (k = 0; k < nLD; k++)
        fprintf(parselog, "%lf\t", lowercut[k]);

    fprintf(parselog, "\nCENTRAL ATOM FILTER\n");
    for (k = 0; k < nLD; k++){
        for (i = 1; i <= atom->ntypes; i++)
            fprintf(parselog, "%d\t", a[k][i]);
        fprintf(parselog, "\n");
    }

    fprintf(parselog, "\nNEIGHBOR ATOM FILTER\n");
    for (k = 0; k < nLD; k++){
        for (i = 1; i <= atom->ntypes; i++)
            fprintf(parselog, "%d\t", b[k][i]);
        fprintf(parselog, "\n");
    }

    fprintf(parselog, "\nRHO_MIN\tRHO_MAX\tDELTA_RHO\n");
    for (k = 0; k < nLD; k++) 
        fprintf(parselog, "%lf\t%lf\t%lf\n", rho_min[k], rho_max[k], delta_rho[k]);

    fprintf(parselog, "\n(RHO, FRHO) AS READ FROM FILE\n");
    for (k = 0; k < nLD; k++) {
        fprintf(parselog, ">>> LD POTENTIAL %d\n", k);
        for (i = 0; i < nrho; i++)
            fprintf(parselog, "%lf\t%lf\n", rho[k][i], frho[k][i]);
    }

    fprintf(parselog, "\nFRHO SPLINE COEFFICIENTS\n");
    for (k = 0; k < nLD; k++) {
        fprintf(parselog, ">>> LD POTENTIAL %d\n", k);
        for (i = 0; i < nrho; i++){
            for (j = 0; j <= 6; j++)
                fprintf(parselog, "%lf\t", frho_spline[k][i][j]);
            fprintf(parselog, "\n");
        }
    }	

    fprintf(parselog, "\nFRHO FROM SPLINE\n");
    for (k = 0; k < nLD; k++) {
        fprintf(parselog, ">>> LD POTENTIAL %d\n", k);
        j = 3*nrho;
        delta = (rho_max[k] - rho_min[k]) / j;
        for (i = 0; i <= j; i++){
            LD = rho_min[k] + i * delta ;
            p = (LD - rho_min[k]) / delta_rho[k];
            m = static_cast<int> (p);
            m = MAX(0, MIN(m, nrho-2));
            p -= m;
            p = MIN(p, 1.0);
            coeff = frho_spline[k][m];
            uLD = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
            duLD = (coeff[0]*p + coeff[1])*p + coeff[2];
            dduLD = 2. * coeff[0] * p + coeff[1];
            fprintf(parselog, "%12.5e\t%12.5e\t%12.5e\t%12.5e\n", LD, uLD, duLD, dduLD);
        }
    }	

    fprintf(parselog, "\nFIRST PART OF FRHO FROM SPLINE\n");
    for (k = 0; k < nLD; k++) {
        fprintf(parselog, ">>> LD POTENTIAL %d\n", k);
        j = 20*nrho;
        delta = (rho_max[k] - rho_min[k]) / j;
        for (i = 0; i <= 60; i++){
            LD = rho_min[k] + i * delta ;
            p = (LD - rho_min[k]) / delta_rho[k];
            m = static_cast<int> (p);
            m = MAX(0, MIN(m, nrho-2));
            p -= m;
            p = MIN(p, 1.0);
            coeff = frho_spline[k][m];
            uLD = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
            duLD = (coeff[0]*p + coeff[1])*p + coeff[2];
            dduLD = 2. * coeff[0] * p + coeff[1];
            fprintf(parselog, "%12.5e\t%12.5e\t%12.5e\t%12.5e\n", LD, uLD, duLD, dduLD);
        }
    }

    fclose(parselog);

 }

}
 



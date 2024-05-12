/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Vikram Thapar, thapar.09@gmail.com
------------------------------------------------------------------------- */

#include "pair_tbsw.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairTBSW::PairTBSW(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}


PairTBSW::~PairTBSW()
{

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(biga);
    memory->destroy(bigb);
    memory->destroy(powerp);
    memory->destroy(powerq);
    memory->destroy(littlea);
    memory->destroy(c1);
    memory->destroy(c2);
    memory->destroy(c3);
    memory->destroy(c4);
    memory->destroy(c5);
    memory->destroy(c6);
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairTBSW::allocate()
{
  allocated = 1;
  int n = atom->ntypes + 1;

  memory->create(setflag, n, n, "pair:setflag");
  for (int i = 1; i < n; i++)
    for (int j = i; j < n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n, n, "pair:cutsq");

  memory->create(cut, n, n, "pair:cut");
  memory->create(epsilon, n, n, "pair:epsilon");
  memory->create(sigma, n, n, "pair:sigma");
  memory->create(biga, n, n, "pair:biga");
  memory->create(bigb, n, n, "pair:bigb");
  memory->create(powerp, n, n, "pair:powerp");
  memory->create(powerq, n, n, "pair:powerq");
  memory->create(littlea, n, n, "pair:littlea");
  memory->create(c1, n, n, "pair:c1");
  memory->create(c2, n, n, "pair:c2");
  memory->create(c3, n, n, "pair:c3");
  memory->create(c4, n, n, "pair:c4");
  memory->create(c5, n, n, "pair:c5");
  memory->create(c6, n, n, "pair:c6");

}


/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairTBSW::settings(int narg, char **arg)
{
  if (narg != 0) error->all(FLERR, "Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairTBSW::coeff(int narg, char **arg)
{
  if (narg != 9) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double epsilon_one = utils::numeric(FLERR, arg[2], false, lmp);
  double sigma_one = utils::numeric(FLERR, arg[3], false, lmp);
  double biga_one = utils::numeric(FLERR, arg[4], false, lmp);
  double bigb_one = utils::numeric(FLERR, arg[5], false, lmp);
  double powerp_one = utils::numeric(FLERR, arg[6], false, lmp);
  double powerq_one = utils::numeric(FLERR, arg[7], false, lmp);
  double littlea_one = utils::numeric(FLERR, arg[8], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      biga[i][j] = biga_one;
      bigb[i][j] = bigb_one;
      powerp[i][j] = powerp_one;
      powerq[i][j] = powerq_one;
      littlea[i][j] = littlea_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairTBSW::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

  cut[i][j] = sigma[i][j]*littlea[i][j];
  c1[i][j] = biga[i][j]*epsilon[i][j]*
    powerp[i][j]*bigb[i][j] *
    pow(sigma[i][j],powerp[i][j]);
  c2[i][j] = biga[i][j]*epsilon[i][j]*powerq[i][j] *
    pow(sigma[i][j],powerq[i][j]);
  c3[i][j] = biga[i][j]*epsilon[i][j]*bigb[i][j] *
    pow(sigma[i][j],powerp[i][j]+1.0);
  c4[i][j] = biga[i][j]*epsilon[i][j] *
    pow(sigma[i][j],powerq[i][j]+1.0);
  c5[i][j] = biga[i][j]*epsilon[i][j]*bigb[i][j] *
    pow(sigma[i][j],powerp[i][j]);
  c6[i][j] = biga[i][j]*epsilon[i][j] *
    pow(sigma[i][j],powerq[i][j]);

  epsilon[j][i] = epsilon[i][j];
  sigma[j][i] = sigma[i][j];
  biga[j][i] = biga[i][j];
  bigb[j][i] = bigb[i][j];
  powerp[j][i] = powerp[i][j];
  powerq[j][i] = powerq[i][j];
  littlea[j][i] = littlea[i][j];
  cut[j][i] = cut[i][j];
  c1[j][i] = c1[i][j];
  c2[j][i] = c2[i][j];
  c3[j][i] = c3[i][j];
  c4[j][i] = c4[i][j];
  c5[j][i] = c5[i][j];
  c6[j][i] = c6[i][j];

  //printf("(epsilon,sigma,biga,bigb,powerp,powerq,littlea) - %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",epsilon[i][j],sigma[i][j],biga[i][j],bigb[i][j],powerp[i][j],powerq[i][j],littlea[i][j]);
  
  //printf("(c1,c2,c3,c4,c5,c6) - %.17g %.17g %.17g %.17g %.17g %.17g\n",c1[i][j],c2[i][j],c3[i][j],c4[i][j],c5[i][j],c6[i][j]);


  return cut[i][j];

}

/* ---------------------------------------------------------------------- */
void PairTBSW::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double r,rinvsq,rp,rq,rainv,rainvsq,expsrainv;
  double rsq, factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

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
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
	r = sqrt(rsq);
	rinvsq = 1.0/rsq;
	rp = pow(r,-powerp[itype][jtype]);
	rq = pow(r,-powerq[itype][jtype]);
	rainv = 1.0 / (r - cut[itype][jtype]);
	rainvsq = rainv*rainv*r;
	expsrainv = exp(sigma[itype][jtype] * rainv);
	fpair = (c1[itype][jtype]*rp - c2[itype][jtype]*rq +
            (c3[itype][jtype]*rp -c4[itype][jtype]*rq) * rainvsq) * expsrainv * rinvsq;
	fpair *= factor_lj;
	
        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;
	//printf("(fpair- %.15g %d %d\n",fpair,i,j);
	//printf("(fx,fy,fz,i,j,fpair - %.15g %.15g %.15g %.15g %d %d\n",f[i][0],f[i][1],f[i][2],fpair,i,j);
	//printf("(fx,fy,fz,i,j - %.15g %.15g %.15g %d %d\n",f[i][0],f[i][1],f[i][2],i,j);
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }

	if (eflag) evdwl = factor_lj*((c5[itype][jtype]*rp - c6[itype][jtype]*rq) * expsrainv);
        if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
    //printf("(fx,fy,fz,i,j - %.15g %.15g %.15g %d\n",f[i][0],f[i][1],f[i][2],i);
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

double PairTBSW::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                         double /*factor_coul*/, double factor_lj, double &fforce)
{
  double r,rinvsq,rp,rq,rainv,rainvsq,expsrainv,eng;

  r = sqrt(rsq);
  rinvsq = 1.0/rsq;
  rp = pow(r,-powerp[itype][jtype]);
  rq = pow(r,-powerq[itype][jtype]);
  rainv = 1.0 / (r - cut[itype][jtype]);
  rainvsq = rainv*rainv*r;
  expsrainv = exp(sigma[itype][jtype] * rainv);
  fforce = (c1[itype][jtype]*rp - c2[itype][jtype]*rq +
            (c3[itype][jtype]*rp -c4[itype][jtype]*rq) * rainvsq) * expsrainv * rinvsq;
  fforce *= factor_lj;

  eng = factor_lj*((c5[itype][jtype]*rp - c6[itype][jtype]*rq) * expsrainv);
  return eng;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairTBSW::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j], sizeof(double), 1, fp);
        fwrite(&sigma[i][j], sizeof(double), 1, fp);
        fwrite(&biga[i][j], sizeof(double), 1, fp);
        fwrite(&bigb[i][j], sizeof(double), 1, fp);
        fwrite(&powerp[i][j], sizeof(double), 1, fp);
        fwrite(&powerq[i][j], sizeof(double), 1, fp);
        fwrite(&littlea[i][j], sizeof(double), 1, fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairTBSW::write_restart_settings(FILE *fp)
{
}


/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairTBSW::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR, &epsilon[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &sigma[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &biga[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &bigb[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &powerp[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &powerq[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &littlea[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&epsilon[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&sigma[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&biga[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&bigb[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&powerp[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&powerq[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&littlea[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}


/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairTBSW::read_restart_settings(FILE *fp)
{
}


/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairTBSW::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++) 
    fprintf(fp, "%d %g %g %g %g %g %g %g\n", i, epsilon[i][i], sigma[i][i], biga[i][i], bigb[i][i],
	    powerp[i][i],powerq[i][i],littlea[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairTBSW::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g %g %g %g %g\n", i, j, epsilon[i][j], sigma[i][j], biga[i][j], 
	      bigb[i][j],powerp[i][j],powerq[i][j],littlea[i][j]);
}

/* ---------------------------------------------------------------------- */

void *PairTBSW::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "epsilon") == 0) return (void *) epsilon;
  if (strcmp(str, "sigma") == 0) return (void *) sigma;
  return nullptr;
}

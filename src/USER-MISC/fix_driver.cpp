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
   Contributing author:  Axel Kohlmeyer (ICTP)
------------------------------------------------------------------------- */

#include "fix_driver.h"
#include "atom.h"
#include "domain.h"
#include "comm.h"
#include "update.h"
#include "force.h"
#include "error.h"
#include "group.h"
#include "memory.h"

#include <stdlib.h>
#include <string.h>

using namespace LAMMPS_NS;
using namespace FixConst;

/****************************************************************************/


/***************************************************************
 * create class and parse arguments in LAMMPS script. Syntax:
 * fix ID group-ID qmmm [couple <group-ID>]
 ***************************************************************/
FixDriver::FixDriver(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{

  if (narg > 3)
    error->all(FLERR,"Illegal fix driver command");

  // allocate arrays
  memory->create(add_force,3*atom->natoms,"driver:add_force");
  for (int i=0; i< 3*atom->natoms; i++) {
    add_force[i] = 0.0;
  }

  if (screen) fprintf(screen,"add_force atoms: %i\n",atom->natoms);
  if (logfile) fprintf(logfile,"add_force atoms: %i\n",atom->natoms);

}

/*********************************
 * Clean up on deleting the fix. *
 *********************************/
FixDriver::~FixDriver()
{

}

/* ---------------------------------------------------------------------- */
int FixDriver::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDriver::exchange_forces()
{
  if (screen) fprintf(screen,"In exchange forces\n");
  if (logfile) fprintf(logfile,"In exchange forces\n");

  double **f = atom->f;
  const int * const mask  = atom->mask;
  const int nlocal = atom->nlocal;

  if (screen) fprintf(screen,"Forces before exchange\n");
  if (logfile) fprintf(logfile,"Forces before exchange\n");
  for (int i=0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      if (screen) fprintf(screen,"f: %i %f %f %f\n",i+1,f[i][0],f[i][1],f[i][2]);
      if (logfile) fprintf(logfile,"f: %i %f %f %f\n",i+1,f[i][0],f[i][1],f[i][2]);
    }
  }

  // add the forces from the driver
  if (screen) fprintf(screen,"Total forces\n");
  if (logfile) fprintf(logfile,"Total forces\n");
  for (int i=0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      f[i][0] += add_force[3*(atom->tag[i]-1)+0];
      f[i][1] += add_force[3*(atom->tag[i]-1)+1];
      f[i][2] += add_force[3*(atom->tag[i]-1)+2];

      if (screen) fprintf(screen,"f: %i %f %f %f\n",i+1,f[i][0],f[i][1],f[i][2]);
      if (logfile) fprintf(logfile,"f: %i %f %f %f\n",i+1,f[i][0],f[i][1],f[i][2]);
    }
  }

}

/* ---------------------------------------------------------------------- */

void FixDriver::init()
{

  /*
  // allocate arrays
  memory->create(add_force,3*atom->natoms,"driver:add_force");
  for (int i=0; i< 3*atom->natoms; i++) {
    add_force[i] = 0.0;
  }

  if (screen) fprintf(screen,"add_force atoms: %i\n",atom->natoms);
  if (logfile) fprintf(logfile,"add_force atoms: %i\n",atom->natoms);
  */

  return;

}

/* ---------------------------------------------------------------------- */

void FixDriver::setup(int)
{
  exchange_forces();
}

/* ---------------------------------------------------------------------- */

void FixDriver::post_force(int vflag)
{
  exchange_forces();
}


// Local Variables:
// mode: c++
// compile-command: "make -j4 openmpi"
// c-basic-offset: 2
// fill-column: 76
// indent-tabs-mode: nil
// End:

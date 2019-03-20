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
   Contributing author: Taylor Barnes (MolSSI)
   MolSSI Driver Interface (MDI) support for LAMMPS
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
#include "modify.h"
#include "compute.h"

#include <stdlib.h>
#include <string.h>

using namespace LAMMPS_NS;
using namespace FixConst;

/****************************************************************************/


/***************************************************************
 * create class and parse arguments in LAMMPS script. Syntax:
 * fix ID group-ID driver [couple <group-ID>]
 ***************************************************************/
FixMDI::FixMDI(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  id_pe(NULL), pe(NULL)
{

  if (narg > 3)
    error->all(FLERR,"Illegal fix mdi command");

  // allocate arrays
  memory->create(add_force,3*atom->natoms,"mdi:add_force");
  for (int i=0; i< 3*atom->natoms; i++) {
    add_force[i] = 0.0;
  }

  // create a new compute pe style
  // id = fix-ID + pe, compute group = all

  int n = strlen(id) + 4;
  id_pe = new char[n];
  strcpy(id_pe,id);
  strcat(id_pe,"_pe");

  char **newarg = new char*[3];
  newarg[0] = id_pe;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pe";
  modify->add_compute(3,newarg);
  delete [] newarg;

}

/*********************************
 * Clean up on deleting the fix. *
 *********************************/
FixMDI::~FixMDI()
{
  modify->delete_compute(id_pe);
  delete [] id_pe;
}

/* ---------------------------------------------------------------------- */
int FixMDI::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMDI::exchange_forces()
{
  double **f = atom->f;
  const int * const mask  = atom->mask;
  const int nlocal = atom->nlocal;

  // add the forces from the driver
  for (int i=0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      f[i][0] += add_force[3*(atom->tag[i]-1)+0];
      f[i][1] += add_force[3*(atom->tag[i]-1)+1];
      f[i][2] += add_force[3*(atom->tag[i]-1)+2];
    }
  }

}

/* ---------------------------------------------------------------------- */

void FixMDI::init()
{

  int icompute = modify->find_compute(id_pe);
  if (icompute < 0)
    error->all(FLERR,"Potential energy ID for fix mdi does not exist");
  pe = modify->compute[icompute];

  return;

}

/* ---------------------------------------------------------------------- */

void FixMDI::setup(int)
{
  exchange_forces();

  //compute the potential energy
  potential_energy = pe->compute_scalar();

  // trigger potential energy computation on next timestep
  pe->addstep(update->ntimestep+1);
}

/* ---------------------------------------------------------------------- */

void FixMDI::post_force(int vflag)
{
  // calculate the energy
  potential_energy = pe->compute_scalar();

  exchange_forces();

  // trigger potential energy computation on next timestep
  pe->addstep(update->ntimestep+1);
}

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

#include <mpi.h>
#include <string.h>
#include "driver.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "min.h"
#include "minimize.h"
#include "modify.h"
#include "output.h"
#include "timer.h"
#include "error.h"
#include "comm.h"
#include "fix_driver.h"
#include "mdi.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

CommandMDI::CommandMDI(LAMMPS *lmp) : Pointers(lmp) {
  return;
}

CommandMDI::~CommandMDI() {
  return;
}

/* ---------------------------------------------------------------------- */

void CommandMDI::command(int narg, char **arg)
{

  // identify the driver fix
  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"mdi") == 0) {
      mdi_fix = static_cast<FixMDI*>(modify->fix[i]);
    }
  }

  /* format for MDI command:
   * mdi
   */
  if (narg > 0) error->all(FLERR,"Illegal MDI command");

  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use MDI command without atom IDs");

  if (atom->tag_consecutive() == 0)
    error->all(FLERR,"MDI command requires consecutive atom IDs");

  // begin engine_mode
  mdi_fix->engine_mode(0);

  return;

}

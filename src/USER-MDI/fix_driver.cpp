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

#include "irregular.h"
#include "min.h"
#include "minimize.h"
#include "modify.h"
#include "output.h"
#include "timer.h"
#include "verlet.h"
#include "mdi.h"

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

  master = (comm->me==0) ? 1 : 0;

  // create instance of the Irregular class
  irregular = new Irregular(lmp);

  most_recent_init = 0;
  exit_flag = false;
  local_exit_flag = false;
  target_node = 0;
  target_command = new char[MDI_COMMAND_LENGTH+1];
  command = new char[MDI_COMMAND_LENGTH+1];

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

  // accept a communicator to the driver
  int ierr;
  if (master) {
    driver_socket = MDI_Accept_Communicator();
    if (driver_socket <= 0)
      error->all(FLERR,"Unable to connect to driver");
  } else driver_socket=0;

}

/*********************************
 * Clean up on deleting the fix. *
 *********************************/
FixMDI::~FixMDI()
{
  modify->delete_compute(id_pe);
  delete irregular;
  delete [] id_pe;
  delete [] target_command;
  delete [] command;
}

/* ---------------------------------------------------------------------- */
int FixMDI::setmask()
{
  int mask = 0;

  // MD masks
  mask |= POST_INTEGRATE;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;

  // Minimizer masks
  mask |= MIN_PRE_FORCE;
  mask |= MIN_POST_FORCE;

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
  //compute the potential energy
  potential_energy = pe->compute_scalar();

  // trigger potential energy computation on next timestep
  pe->addstep(update->ntimestep+1);

  if ( most_recent_init == 1 ) { // md
    // @PRE-FORCES
    engine_mode(2);
  }
}

/* ---------------------------------------------------------------------- */

void FixMDI::post_integrate()
{
  engine_mode(1);
}

/* ---------------------------------------------------------------------- */

void FixMDI::post_force(int vflag)
{
  // calculate the energy
  potential_energy = pe->compute_scalar();

  // @PRE-FORCES
  engine_mode(2);

  // trigger potential energy computation on next timestep
  pe->addstep(update->ntimestep+1);
}


/* ---------------------------------------------------------------------- */

void FixMDI::min_pre_force(int vflag)
{
  // @COORDS
  engine_mode(1);
}

/* ---------------------------------------------------------------------- */

void FixMDI::min_post_force(int vflag)
{
  // calculate the energy
  potential_energy = pe->compute_scalar();

  // @FORCES
  engine_mode(3);

  // trigger potential energy computation on next timestep
  pe->addstep(update->ntimestep+1);
}

/* ---------------------------------------------------------------------- */

void FixMDI::end_of_step()
{
  if ( most_recent_init == 1 ) { // md
    // @FORCES
    engine_mode(3);
  }
  else if ( most_recent_init == 2 ) { // optg
    // @FORCES
    engine_mode(3);
  }
}

/* ---------------------------------------------------------------------- */

char *FixMDI::engine_mode(int node)
{
  if (screen)
    fprintf(screen,"MDI ENGINE MODE: %i\n",node);
  if (logfile)
    fprintf(logfile,"MDI ENGINE MODE: %i\n",node);

  // flag to indicate whether the engine should continue listening for commands at this node
  current_node = node;
  if ( target_node != 0 and target_node != current_node ) {
    local_exit_flag = true;
  }

  /* ----------------------------------------------------------------- */
  // Answer commands from the driver
  /* ----------------------------------------------------------------- */

  while (not exit_flag and not local_exit_flag) {

    if (master) { 
      // read the next command from the driver
      ierr = MDI_Recv_Command(command, driver_socket);
      if (ierr != 0)
        error->all(FLERR,"Unable to receive command from driver");
      command[MDI_COMMAND_LENGTH]=0;
    }
    // broadcast the command to the other tasks
    MPI_Bcast(command,MDI_COMMAND_LENGTH,MPI_CHAR,0,world);

    if (screen)
      fprintf(screen,"MDI command: %s\n",command);
    if (logfile)
      fprintf(logfile,"MDI command: %s:\n",command);

    if (strcmp(command,"STATUS      ") == 0 ) {
      // send the calculation status to the driver
      if (master) {
	ierr = MDI_Send_Command("READY", driver_socket);
        if (ierr != 0)
          error->all(FLERR,"Unable to return status to driver");
      }
    }
    else if (strcmp(command,">NATOMS") == 0 ) {
      // receive the number of atoms from the driver
      if (master) {
        ierr = MDI_Recv((char*) &atom->natoms, 1, MDI_INT, driver_socket);
        if (ierr != 0)
          error->all(FLERR,"Unable to receive number of atoms from driver");
      }
      MPI_Bcast(&atom->natoms,1,MPI_INT,0,world);
    }
    else if (strcmp(command,"<NATOMS") == 0 ) {
      // send the number of atoms to the driver
      if (master) {
        ierr = MDI_Send((char*) &atom->natoms, 1, MDI_INT, driver_socket);
        if (ierr != 0)
          error->all(FLERR,"Unable to send number of atoms to driver");
      }
    }
    else if (strcmp(command,"<NTYPES") == 0 ) {
      // send the number of atom types to the driver
      if (master) {
        ierr = MDI_Send((char*) &atom->ntypes, 1, MDI_INT, driver_socket);
        if (ierr != 0)
          error->all(FLERR,"Unable to send number of atom types to driver");
      }
    }
    else if (strcmp(command,"<TYPES") == 0 ) {
      // send the atom types
      send_types(error);
    }
    else if (strcmp(command,"<MASSES") == 0 ) {
      // send the atom types
      send_masses(error);
    }
    else if (strcmp(command,"<CELL") == 0 ) {
      // send the cell dimensions to the driver
      send_cell(error);
    }
    else if (strcmp(command,">COORDS") == 0 ) {
      // receive the coordinate information
      receive_coordinates(error);
    }
    else if (strcmp(command,"<COORDS") == 0 ) {
      // send the coordinate information
      send_coordinates(error);
    }
    else if (strcmp(command,"<CHARGES") == 0 ) {
      // send the charges
      send_charges(error);
    }
    else if (strcmp(command,"<ENERGY") == 0 ) {
      // send the potential energy to the driver
      send_energy(error);
    }
    else if (strcmp(command,"<FORCES") == 0 ) {
      // send the forces to the driver
      send_forces(error);
    }
    else if (strcmp(command,">FORCES") == 0 ) {
      // receive the forces from the driver
      receive_forces(error);
    }
    else if (strcmp(command,"+FORCES") == 0 ) {
      // receive additional forces from the driver
      // these are added prior to SHAKE or other post-processing
      add_forces(error);
    }
    else if (strcmp(command,"MD_INIT") == 0 ) {
      // initialize a new MD simulation
      //md_init(error);
      most_recent_init = 1;
      local_exit_flag = true;
    }
    else if (strcmp(command,"OPTG_INIT") == 0 ) {
      // initialize a new geometry optimization
      optg_init(error);
    }
    else if (strcmp(command,"@") == 0 ) {
      target_node = 0;
      local_exit_flag = true;
    }
    else if (strcmp(command,"<@") == 0 ) {
      if (master) {
	if ( current_node == 1 ) {
	  ierr = MDI_Send("@COORDS", MDI_NAME_LENGTH, MDI_CHAR, driver_socket);
	}
	else if ( current_node == 2 ) {
	  ierr = MDI_Send("@PRE-FORCES", MDI_NAME_LENGTH, MDI_CHAR, driver_socket);
	}
	else if (current_node == 3 ) {
	  ierr = MDI_Send("@FORCES", MDI_NAME_LENGTH, MDI_CHAR, driver_socket);
	}
	if (ierr != 0)
	  error->all(FLERR,"Unable to send node to driver");
      }
    }
    else if (strcmp(command,"@COORDS") == 0 ) {
      target_node = 1;
      local_exit_flag = true;
    }
    else if (strcmp(command,"@PRE-FORCES") == 0 ) {
      target_node = 2;
      local_exit_flag = true;
    }
    else if (strcmp(command,"@FORCES") == 0 ) {
      target_node = 3;
      local_exit_flag = true;
    }
    else if (strcmp(command,"MD_EXIT") == 0 ) {
      most_recent_init = 0;

      // proceed to the @FORCES node, which corresponds to the original engine_mode call
      target_node = 3;
      local_exit_flag = true;
    }
    else if (strcmp(command,"OPTG_EXIT") == 0 ) {
      most_recent_init = 0;
    }
    else if (strcmp(command,"EXIT") == 0 ) {
      // exit the driver code
      exit_flag = true;

      // if doing a geometry optimization, set the maximum number of evaluations to 0
      if ( most_recent_init == 2 ) {
	update->max_eval = 0;
      }
    }
    else {
      // the command is not supported
      error->all(FLERR,"Unknown command from driver");
    }

    // check if the target node is something other than the current node
    if ( target_node != 0 and target_node != current_node ) {
      local_exit_flag = true;
    }

  }

  // a local exit has completed, so turn off the local exit flag
  local_exit_flag = false;

  return command;

}


void FixMDI::receive_coordinates(Error* error)
{
  double posconv;
  posconv=force->angstrom/MDI_ANGSTROM_TO_BOHR;

  // create a buffer to hold the coordinates
  double *buffer;
  buffer = new double[3*atom->natoms];

  if (master) {
    ierr = MDI_Recv((char*) buffer, 3*atom->natoms, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"Unable to receive coordinates from driver");
  }
  MPI_Bcast(buffer,3*atom->natoms,MPI_DOUBLE,0,world);

  // pick local atoms from the buffer
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    x[i][0]=buffer[3*(atom->tag[i]-1)+0]*posconv;
    x[i][1]=buffer[3*(atom->tag[i]-1)+1]*posconv;
    x[i][2]=buffer[3*(atom->tag[i]-1)+2]*posconv;
  }

  // ensure atoms are in current box & update box via shrink-wrap
  // has to be be done before invoking Irregular::migrate_atoms() 
  //   since it requires atoms be inside simulation box
  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // move atoms to new processors via irregular()
  // only needed if migrate_check() says an atom moves to far
  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  if (irregular->migrate_check()) irregular->migrate_atoms();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  delete [] buffer;
}


void FixMDI::send_coordinates(Error* error)
{
  double posconv;
  posconv=force->angstrom/MDI_ANGSTROM_TO_BOHR;

  double *coords;
  double *coords_reduced;

  coords = new double[3*atom->natoms];
  coords_reduced = new double[3*atom->natoms];

  // pick local atoms from the buffer
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    coords[3*(atom->tag[i]-1)+0] = x[i][0]/posconv;
    coords[3*(atom->tag[i]-1)+1] = x[i][1]/posconv;
    coords[3*(atom->tag[i]-1)+2] = x[i][2]/posconv;
  }

  MPI_Reduce(coords, coords_reduced, 3*atom->natoms, MPI_DOUBLE, MPI_SUM, 0, world);

  if (master) {
    ierr = MDI_Send((char*) coords_reduced, 3*atom->natoms, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"Unable to send coordinates to driver");
  }

  delete [] coords;
  delete [] coords_reduced;
}


void FixMDI::send_charges(Error* error)
{
  double *charges;
  double *charges_reduced;

  charges = new double[atom->natoms];
  charges_reduced = new double[atom->natoms];

  // pick local atoms from the buffer
  double *charge = atom->q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    charges[atom->tag[i]-1] = charge[i];
  }

  MPI_Reduce(charges, charges_reduced, atom->natoms, MPI_DOUBLE, MPI_SUM, 0, world);

  if (master) { 
    ierr = MDI_Send((char*) charges_reduced, atom->natoms, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"Unable to send charges to driver");
  }

  delete [] charges;
  delete [] charges_reduced;
}


void FixMDI::send_energy(Error* error)
{

  double pe;
  double *send_pe = &pe;

  pe = potential_energy;

  // convert the energy to atomic units
  pe *= MDI_KELVIN_TO_HARTREE/force->boltz;

  if (master) {
    ierr = MDI_Send((char*) send_pe, 1, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"Unable to send potential energy to driver");
  }
}


void FixMDI::send_types(Error* error)
{
  int * const type = atom->type;

  if (master) { 
    ierr = MDI_Send((char*) type, atom->natoms, MDI_INT, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"Unable to send atom types to driver");
  }
}


void FixMDI::send_masses(Error* error)
{
  double * const mass = atom->mass;
  int * const type = atom->type;

  if (master) { 
    double *mass_by_atom = new double[atom->natoms];
    for (int iatom=0; iatom < atom->natoms; iatom++) {
      mass_by_atom[iatom] = mass[ type[iatom] ];
    }
    ierr = MDI_Send((char*) mass_by_atom, atom->natoms, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"Unable to send atom masses to driver");
    delete [] mass_by_atom;
  }

}


void FixMDI::send_forces(Error* error)
{
  double potconv, posconv, forceconv;
  potconv=MDI_KELVIN_TO_HARTREE/force->boltz;
  posconv=force->angstrom/MDI_ANGSTROM_TO_BOHR;
  forceconv=potconv*posconv;

  double *forces;
  double *forces_reduced;
  double *x_buf;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  forces = new double[3*atom->natoms];
  forces_reduced = new double[3*atom->natoms];
  x_buf = new double[3*atom->natoms];

  // if not at a node, calculate the forces
  if ( current_node == 0 ) {
    // certain fixes, such as shake, move the coordinates
    // to ensure that the coordinates do not change, store a copy
    double **x = atom->x;
    for (int i = 0; i < nlocal; i++) {
      x_buf[3*i+0] = x[i][0];
      x_buf[3*i+1] = x[i][1];
      x_buf[3*i+2] = x[i][2];
    }

    // calculate the forces
    update->whichflag = 1; // 1 for dynamics
    update->nsteps = 1;
    lmp->init();
    update->integrate->setup_minimal(1);
  }

  // pick local atoms from the buffer
  double **f = atom->f;
  for (int i = 0; i < nlocal; i++) {
    forces[3*(atom->tag[i]-1)+0] = f[i][0]*forceconv;
    forces[3*(atom->tag[i]-1)+1] = f[i][1]*forceconv;
    forces[3*(atom->tag[i]-1)+2] = f[i][2]*forceconv;
  }

  // reduce the forces onto rank 0
  MPI_Reduce(forces, forces_reduced, 3*atom->natoms, MPI_DOUBLE, MPI_SUM, 0, world);

  // send the forces through MDI
  if (master) {
    ierr = MDI_Send((char*) forces_reduced, 3*atom->natoms, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"Unable to send atom forces to driver");
  }

  if ( current_node == 0 ) {
    // restore the original set of coordinates
    double **x_new = atom->x;
    for (int i = 0; i < nlocal; i++) {
      x_new[i][0] = x_buf[3*i+0];
      x_new[i][1] = x_buf[3*i+1];
      x_new[i][2] = x_buf[3*i+2];
    }
  }

  delete [] forces;
  delete [] forces_reduced;
  delete [] x_buf;

}


void FixMDI::receive_forces(Error* error)
{
  double potconv, posconv, forceconv;
  potconv=MDI_KELVIN_TO_HARTREE/force->boltz;
  posconv=force->angstrom/MDI_ANGSTROM_TO_BOHR;
  forceconv=potconv*posconv;

  double *forces;
  forces = new double[3*atom->natoms];

  if (master) {
    ierr = MDI_Recv((char*) forces, 3*atom->natoms, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"Unable to receive atom forces to driver");
  }
  MPI_Bcast(forces,3*atom->natoms,MPI_DOUBLE,0,world);

  // pick local atoms from the buffer
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    f[i][0] = forces[3*(atom->tag[i]-1)+0]/forceconv;
    f[i][1] = forces[3*(atom->tag[i]-1)+1]/forceconv;
    f[i][2] = forces[3*(atom->tag[i]-1)+2]/forceconv;
  }

  delete [] forces;
}


void FixMDI::add_forces(Error* error)
{
  double potconv, posconv, forceconv;
  potconv=MDI_KELVIN_TO_HARTREE/force->boltz;
  posconv=force->angstrom * MDI_Conversion_Factor("angstrom","bohr");
  forceconv=potconv*posconv;

  double *forces;
  forces = new double[3*atom->natoms];

  if (master) {
    ierr = MDI_Recv((char*) forces, 3*atom->natoms, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"Unable to receive atom forces to driver");
  }
  MPI_Bcast(forces,3*atom->natoms,MPI_DOUBLE,0,world);

  // pick local atoms from the buffer
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    f[i][0] += forces[3*(atom->tag[i]-1)+0]/forceconv;
    f[i][1] += forces[3*(atom->tag[i]-1)+1]/forceconv;
    f[i][2] += forces[3*(atom->tag[i]-1)+2]/forceconv;
  }

  delete [] forces;
}


void FixMDI::send_cell(Error* error)
{
  double celldata[12];

  celldata[0] = domain->boxhi[0] - domain->boxlo[0];
  celldata[1] = 0.0;
  celldata[2] = 0.0;
  celldata[3] = domain->xy;
  celldata[4] = domain->boxhi[1] - domain->boxlo[1];
  celldata[5] = 0.0;
  celldata[6] = domain->xz;
  celldata[7] = domain->yz;
  celldata[8] = domain->boxhi[2] - domain->boxlo[2];
  celldata[9 ] = domain->boxlo[0];
  celldata[10] = domain->boxlo[1];
  celldata[11] = domain->boxlo[2];

  // convert the units to bohr
  double unit_conv = force->angstrom * MDI_Conversion_Factor("angstrom","bohr");
  for (int icell=0; icell < 12; icell++) {
    celldata[icell] *= unit_conv;
  }

  if (master) { 
    ierr = MDI_Send((char*) celldata, 12, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"Unable to send cell dimensions to driver");
  }
}


void FixMDI::optg_init(Error* error)
{
  if ( most_recent_init != 0 ) {
      error->all(FLERR,"Atomic propagation method already initialized");
  }

  // create instance of Minimizer class
  minimizer = new Minimize(lmp);

  // initialize the minimizer in a way that ensures optimization will continue until the driver exits
  int narg = 4;
  char* arg[] = {"1.0e-100","1.0e-100","1000000000","1000000000"};
  update->etol = force->numeric(FLERR,arg[0]);
  update->ftol = force->numeric(FLERR,arg[1]);
  update->nsteps = force->inumeric(FLERR,arg[2]);
  update->max_eval = force->inumeric(FLERR,arg[3]);

  update->whichflag = 2; // 2 for minimization
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + update->nsteps;
  lmp->init();
  update->minimize->setup();

  current_node = -1; // after OPTG_INIT
  most_recent_init = 2;

  update->minimize->iterate(1000000000);
}

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
#include "modify.h"
#include "output.h"
#include "timer.h"
#include "error.h"
#include "comm.h"
#include "irregular.h"
#include "verlet.h"
#include "fix_driver.h"
extern "C" {
#include "mdi.h"
}

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Driver::Driver(LAMMPS *lmp) : Pointers(lmp) {
  md_initialized = false;

  // create instance of Irregular class
  irregular = new Irregular(lmp);
}

Driver::~Driver() {
  delete irregular;
}

/* ---------------------------------------------------------------------- */

void Driver::command(int narg, char **arg)
{
  /* format for driver command:
   * driver hostname port mdi_name [unix]
   */
  if (narg < 3) error->all(FLERR,"Illegal driver command");

  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use driver command without atom IDs");

  if (atom->tag_consecutive() == 0)
    error->all(FLERR,"Driver command requires consecutive atom IDs");

  // obtain host information from the command arguments
  host = strdup(arg[0]);
  port = force->inumeric(FLERR,arg[1]);
  const char* mdi_name = strdup(arg[2]);
  inet   = ((narg > 3) && (strcmp(arg[3],"unix") == 0) ) ? 0 : 1;

  master = (comm->me==0) ? 1 : 0;

  // open the socket
  int ierr;
  if (master) {
    driver_socket = MDI_Open(inet, port, host);
    if (driver_socket <= 0)
      error->all(FLERR,"Unable to connect to driver");
  } else driver_socket=0;

  /* ----------------------------------------------------------------- */
  // Answer commands from the driver
  /* ----------------------------------------------------------------- */
  char command[MDI_COMMAND_LENGTH+1];

  bool exit_flag = false;

  while (not exit_flag) {

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
    else if (strcmp(command,"<NAME       ") == 0 ) {
      // send the calculation name to the driver
      if (master) {
	ierr = MDI_Send_Command(mdi_name, driver_socket);
        if (ierr != 0)
          error->all(FLERR,"Unable to send name to driver");
      }
    }
    else if (strcmp(command,">NATOMS     ") == 0 ) {
      // receive the number of atoms from the driver
      if (master) {
        ierr = MDI_Recv((char*) &atom->natoms, 1, MDI_INT, driver_socket);
        if (ierr != 0)
          error->all(FLERR,"Unable to receive number of atoms from driver");
      }
      MPI_Bcast(&atom->natoms,1,MPI_INT,0,world);
    }
    else if (strcmp(command,"<NATOMS     ") == 0 ) {
      // send the number of atoms to the driver
      if (master) {
        ierr = MDI_Send((char*) &atom->natoms, 1, MDI_INT, driver_socket);
        if (ierr != 0)
          error->all(FLERR,"Unable to send number of atoms to driver");
      }
    }
    else if (strcmp(command,"<NTYPES     ") == 0 ) {
      // send the number of atom types to the driver
      if (master) {
        ierr = MDI_Send((char*) &atom->ntypes, 1, MDI_INT, driver_socket);
        if (ierr != 0)
          error->all(FLERR,"Unable to send number of atom types to driver");
      }
    }
    else if (strcmp(command,"<TYPES      ") == 0 ) {
      // send the atom types
      send_types(error);
    }
    else if (strcmp(command,"<MASS       ") == 0 ) {
      // send the atom types
      send_masses(error);
    }
    else if (strcmp(command,"<CELL       ") == 0 ) {
      // send the cell dimensions to the driver
      send_cell(error);
    }
    else if (strcmp(command,">COORDS     ") == 0 ) {
      // receive the coordinate information
      receive_coordinates(error);
    }
    else if (strcmp(command,"<COORDS     ") == 0 ) {
      // send the coordinate information
      send_coordinates(error);
    }
    else if (strcmp(command,"<CHARGES    ") == 0 ) {
      // send the charges
      send_charges(error);
    }
    else if (strcmp(command,"<ENERGY     ") == 0 ) {
      // send the potential energy to the driver
      send_energy(error);
    }
    else if (strcmp(command,"<FORCES     ") == 0 ) {
      // send the forces to the driver
      send_forces(error);
    }
    else if (strcmp(command,">FORCES     ") == 0 ) {
      // receive the forces from the driver
      receive_forces(error);
    }
    else if (strcmp(command,"+PRE-FORCES     ") == 0 ) {
      // receive additional forces from the driver
      // these are added prior to SHAKE or other post-processing
      add_forces(error);
    }
    else if (strcmp(command,"MD_INIT     ") == 0 ) {
      // initialize a new MD simulation
      md_init(error);
    }
    else if (strcmp(command,"TIMESTEP    ") == 0 ) {
      // perform an MD timestep
      timestep(error);
    }
    else if (strcmp(command,"EXIT        ") == 0 ) {
      // exit the driver code
      exit_flag = true;
    }
    else {
      // the command is not supported
      error->all(FLERR,"Unknown command from driver");
    }

  }

}


void Driver::receive_coordinates(Error* error)
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


void Driver::send_coordinates(Error* error)
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


void Driver::send_charges(Error* error)
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


void Driver::send_energy(Error* error)
{

  double pe;
  double *potential_energy = &pe;

  // be certain that the MD simulation has been initialized
  if ( not md_initialized ) {
    md_init(error);
  }

  // identify the driver fix
  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"driver") == 0) {
      FixDriver *fixd = static_cast<FixDriver*>(modify->fix[i]);
      pe = fixd->potential_energy;
    }
  }

  // convert the energy to atomic units
  pe *= MDI_KELVIN_TO_HARTREE/force->boltz;

  if (master) {
    ierr = MDI_Send((char*) potential_energy, 1, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"Unable to send potential energy to driver");
  }
}


void Driver::send_types(Error* error)
{
  int * const type = atom->type;

  if (master) { 
    ierr = MDI_Send((char*) type, atom->natoms, MDI_INT, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"Unable to send atom types to driver");
  }
}


void Driver::send_masses(Error* error)
{
  double * const mass = atom->mass;

  if (master) { 
    ierr = MDI_Send((char*) mass, atom->ntypes+1, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"Unable to send atom masses to driver");
  }
}


void Driver::send_forces(Error* error)
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

  //certain fixes, such as shake, move the coordinates
  //to ensure that the coordinates do not change, store a copy
  double **x = atom->x;
  for (int i = 0; i < nlocal; i++) {
    x_buf[3*i+0] = x[i][0];
    x_buf[3*i+1] = x[i][1];
    x_buf[3*i+2] = x[i][2];
  }


  // calculate the forces
  update->whichflag = 1; // 0 for forces
  update->nsteps = 1;
  lmp->init();
  update->integrate->setup_minimal(1);

  // pick local atoms from the buffer
  double **f = atom->f;
  for (int i = 0; i < nlocal; i++) {
    forces[3*(atom->tag[i]-1)+0] = f[i][0]*forceconv;
    forces[3*(atom->tag[i]-1)+1] = f[i][1]*forceconv;
    forces[3*(atom->tag[i]-1)+2] = f[i][2]*forceconv;
  }

  MPI_Reduce(forces, forces_reduced, 3*atom->natoms, MPI_DOUBLE, MPI_SUM, 0, world);

  if (master) {
    ierr = MDI_Send((char*) forces_reduced, 3*atom->natoms, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"Unable to send atom forces to driver");
  }

  //restore the original set of coordinates
  double **x_new = atom->x;
  for (int i = 0; i < nlocal; i++) {
    x_new[i][0] = x_buf[3*i+0];
    x_new[i][1] = x_buf[3*i+1];
    x_new[i][2] = x_buf[3*i+2];
  }

  delete [] forces;
  delete [] forces_reduced;
  delete [] x_buf;

}


void Driver::receive_forces(Error* error)
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


void Driver::add_forces(Error* error)
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
      error->all(FLERR,"Unable to receive atom +forces to driver");
  }
  MPI_Bcast(forces,3*atom->natoms,MPI_DOUBLE,0,world);
  for (int i = 0; i < 3*atom->natoms; i++) {
    forces[i] /= forceconv;
  }

  //identify the driver fix
  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"driver") == 0) {
      FixDriver *fixd = static_cast<FixDriver*>(modify->fix[i]);
      for (int j = 0; j < 3*atom->natoms; j++) {
	fixd->add_force[j] = forces[j];
      }
    }
  }

  delete [] forces;
}


void Driver::send_cell(Error* error)
{
  double celldata[9];

  celldata[0] = domain->boxlo[0];
  celldata[1] = domain->boxlo[1];
  celldata[2] = domain->boxlo[2];
  celldata[3] = domain->boxhi[0];
  celldata[4] = domain->boxhi[1];
  celldata[5] = domain->boxhi[2];
  celldata[6] = domain->xy;
  celldata[7] = domain->xz;
  celldata[8] = domain->yz;

  if (master) { 
    ierr = MDI_Recv((char*) celldata, 9, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"Unable to receive cell dimensions from driver");
  }
}


void Driver::md_init(Error* error)
{
  // calculate the forces
  update->whichflag = 1; // 1 for dynamics
  timer->init_timeout();
  update->nsteps = 1;
  update->ntimestep = 0;
  update->firststep = update->ntimestep;
  update->laststep = update->ntimestep + update->nsteps;
  update->beginstep = update->firststep;
  update->endstep = update->laststep;
  lmp->init();
  update->integrate->setup(1);

  md_initialized = true;
}


void Driver::timestep(Error* error)
{
  // calculate the forces
  update->whichflag = 1; // 1 for dynamics
  timer->init_timeout();
  update->nsteps += 1;
  update->laststep += 1;
  update->endstep = update->laststep;
  output->next = update->ntimestep + 1;

  update->integrate->run(1);

}

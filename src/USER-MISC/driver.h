/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(driver,Driver)

#else

#ifndef LMP_DRIVER_H
#define LMP_DRIVER_H

#include "pointers.h"

namespace LAMMPS_NS {

class Driver : protected Pointers {
 public:
  Driver(class LAMMPS *);
  void command(int, char **);
  //double *add_force; // stores forces added using +FORCE command

 protected:
  void send_types(Error *);
  void send_masses(Error *);
  void read_coordinates(Error *);
  void send_coordinates(Error *);
  void send_charges(Error *);
  void write_forces(Error *);
  void add_forces(Error *);
  void receive_forces(Error *);
  void send_cell(Error *);
  void md_init(Error *);
  void timestep(Error *);
  char *host; int port; int inet, master;
  int driver_socket;
  
  int nat;

private:
  class Irregular *irregular;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Run command before simulation box is defined

The run command cannot be used before a read_data, read_restart, or
create_box command.

E: Invalid run command N value

The number of timesteps must fit in a 32-bit integer.  If you want to
run for more steps than this, perform multiple shorter runs.

E: Invalid run command upto value

Self-explanatory.

E: Invalid run command start/stop value

Self-explanatory.

E: Run command start value is after start of run

Self-explanatory.

E: Run command stop value is before end of run

Self-explanatory.

E: Too many timesteps

The cumulative timesteps must fit in a 64-bit integer.

*/

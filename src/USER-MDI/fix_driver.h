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


#ifdef FIX_CLASS

FixStyle(mdi,FixMDI)

#else

#ifndef LMP_FIX_MDI_H
#define LMP_FIX_MDI_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMDI : public Fix {
 public:
  FixMDI(class LAMMPS *, int, char **);
  ~FixMDI();
  int setmask();
  void init();

  void engine_mode(int node);

  // receive and update forces
  void setup(int);
  void post_integrate();
  void post_force(int);
  void end_of_step();

  double *add_force; // stores forces added using +FORCE command
  double potential_energy; // stores potential energy

 protected:
  void exchange_forces();

 private:
  int master, ierr;
  int driver_socket;
  int most_recent_init; // which MDI init command was most recently received?
                        // 0 - none
                        // 1 - MD
                        // 2 - OPTG
  bool exit_flag;
  bool local_exit_flag;
  int current_node;
  int target_node;      // is the code supposed to advance to a particular node?
                        // 0 - none
                        // 1 - @COORDS (before pre-force calculation)
                        // 2 - @PREFORCES (before final force calculation)
                        // 3 - @FORCES (before time integration)
                        // 4 - after MD_INIT command

  // command to be executed at the target node
  char *target_command;

  char *id_pe;
  class Irregular *irregular;
  class Minimize *minimizer;
  class Compute *pe;
  void send_types(Error *);
  void send_masses(Error *);
  void receive_coordinates(Error *);
  void send_coordinates(Error *);
  void send_charges(Error *);
  void send_energy(Error *);
  void send_forces(Error *);
  void add_forces(Error *);
  void receive_forces(Error *);
  void send_cell(Error *);
  void md_init(Error *);
  void md_setup(Error *);
  void timestep(Error *);
  void optg_init(Error *);

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.

E: Potential energy ID for fix mdi does not exist

Self-explanatory.

*/

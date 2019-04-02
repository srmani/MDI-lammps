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

CommandStyle(mdi,CommandMDI)

#else

#ifndef LMP_COMMAND_MDI_H
#define LMP_COMMAND_MDI_H

#include "pointers.h"

namespace LAMMPS_NS {

class CommandMDI : protected Pointers {
 public:
  CommandMDI(class LAMMPS *);
  virtual ~CommandMDI();
  void command(int, char **);

 protected:
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
  void timestep(Error *);
  void optg_init(Error *);
  int inet, master, ierr;
  int driver_socket;
  
  int nat;
  bool md_initialized; // has the MD simulation been initialized yet?

private:
  class Irregular *irregular;
  class Minimize *minimizer;
  class FixMDI *mdi_fix;
  int most_recent_init; // which MDI init command was most recently received?
                        // 0 - none
                        // 1 - MD
                        // 2 - OPTG

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use MDI command without atom IDs

Self-explanatory.

E: MDI command requires consecutive atom IDs

Self-explanatory.

E: Unable to connect to driver

Self-explanatory.

E: Unable to ... driver

Self-explanatory.

E: Unknown command from driver

The driver sent a command that is not supported by the LAMMPS 
interface.  In some cases this might be because a nonsensical 
command was sent (i.e. "SCF").  In other cases, the LAMMPS 
interface might benefit from being expanded.

*/

/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(tbsw,PairTBSW);
// clang-format on
#else

#ifndef LMP_PAIR_TBSW_H
#define LMP_PAIR_TBSW_H

#include "pair.h"

namespace LAMMPS_NS {

class PairTBSW : public Pair {
 public:
  PairTBSW(class LAMMPS *);
  virtual ~PairTBSW();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

 protected:
  // double cut_global;
  double **cut;
  double **epsilon, **sigma, **biga, **bigb, **powerp, **powerq, **littlea;
  double **c1, **c2, **c3, **c4, **c5, **c6;
  
  virtual void allocate();

};

}   // namespace LAMMPS_NS

#endif
#endif

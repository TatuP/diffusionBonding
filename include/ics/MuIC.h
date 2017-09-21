/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef MUIC_H
#define MUIC_H

// MOOSE Includes
#include "InitialCondition.h"

// Forward Declarations
class MuIC;

template <>
InputParameters validParams<MuIC>();

/**
 * MuIC just returns a constant value.
 */
class MuIC : public InitialCondition
{
public:
  /**
   * Constructor: Same as the rest of the MOOSE Objects
   */
  MuIC(const InputParameters & parameters);

  /**
   * The value of the variable at a point.
   *
   * This must be overriden by derived classes.
   */
  virtual Real value(const Point & p) override;

private:
  const Real _energy_scale; 
  const Real _c_init_1;
  const Real _c_init_2;
  const Real _c_init_3;
  const Real _c_eq_1;
  const Real _c_eq_2;
  const Real _c_eq_3;
  const Real _G_c_1;
  const Real _G_c_2;
  const Real _G_c_3;
  const Real _G_cc_1;
  const Real _G_cc_2;
  const Real _G_cc_3;
  const VariableValue & _phi_1;
  const VariableValue & _phi_2;
  const bool _initialize_step;
  const Real _step_position; 
};

#endif // EXAMPLEIC_H

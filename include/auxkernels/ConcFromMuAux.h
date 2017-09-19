/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CONCFROMMUAUX_H
#define CONCFROMMUAUX_H

#include "AuxKernel.h"

// Forward Declarations
class ConcFromMuAux;

template <>
InputParameters validParams<ConcFromMuAux>();

/**
 * Visualize the location of grain boundaries in a polycrystalline simulation.
 */
class ConcFromMuAux : public AuxKernel
{
public:
  ConcFromMuAux(const InputParameters & parameters);

protected:
  virtual Real computeValue();
private: 
  const Real _c_eq_1; 
  const Real _c_eq_2; 
  const Real _c_eq_3; 
  const Real _G_c_1; 
  const Real _G_c_2; 
  const Real _G_c_3; 
  const Real _inv_G_cc_1; 
  const Real _inv_G_cc_2; 
  const Real _inv_G_cc_3; 
  const VariableValue & _phi_1;
  const VariableValue & _phi_2;
  const VariableValue & _mu;
};

#endif // CONCFROMMUAUX_H 

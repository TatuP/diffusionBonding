/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef DIFFUSIONCOEFFICIENTSMATERIAL_H
#define DIFFUSIONCOEFFICIENTSMATERIAL_H

#include "Material.h"

// Forward Declarations
class DiffusionCoefficientsMaterial;

template <>
InputParameters validParams<DiffusionCoefficientsMaterial>();

/**
 * Grain boundary energy parameters for isotropic uniform grain boundary energies
 */
class DiffusionCoefficientsMaterial : public Material//GBEvolutionBaseUDriven
{
public:
  DiffusionCoefficientsMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

private:
  //MaterialProperty<Real> & _tau_anisotropy; 
  //MaterialProperty<Real> & _dtau_anisotropy_dop; 

  // earlier all of these were const
  const Real _diffusion_coefficient_1;
  const Real _diffusion_coefficient_2;
  const Real _diffusion_coefficient_3;

  const VariableValue & _phi_1;
  const VariableValue & _phi_2;
  const VariableValue & _phi_3;

  MaterialProperty<Real> & _diffusion_coefficient; 

//  const VariableValue & _op;
//  const VariableGradient & _grad_op;

//  Real _GBEnergy;
};

#endif // GBEVOLUTIONUDRIVEN_H

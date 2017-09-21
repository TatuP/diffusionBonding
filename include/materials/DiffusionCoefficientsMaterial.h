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
  const Real _R;
  const Real _length_scale;
  const Real _time_scale;
  const Real _non_dimensionalizer; 
  const Real _D_1_pre;
  const Real _D_2_pre;
  const Real _D_3_pre;
  const Real _Q_1;
  const Real _Q_2;
  const Real _Q_3;
  const VariableValue & _phi_1;
  const VariableValue & _phi_2;
  Function & _TemperatureRampFunction; 
  MaterialProperty<Real> & _diffusion_coefficient; 

  const Real _c_eq_1; 
  const Real _c_eq_2; 
  const Real _c_eq_3; 

  const Real _energy_scale;
  const Real _G_c_1; 
  const Real _G_c_2; 
  const Real _G_c_3; 
  const Real _G_cc_1; 
  const Real _G_cc_2; 
  const Real _G_cc_3; 
  const Real _inv_G_cc_1; 
  const Real _inv_G_cc_2; 
  const Real _inv_G_cc_3; 
  const Real _G_ccc_1; 
  const Real _G_ccc_2; 
  const Real _G_ccc_3; 
  const VariableValue & _mu;
  MaterialProperty<Real> & _concFromMu; 
  MaterialProperty<Real> & _mobility; 



//  const VariableValue & _op;
//  const VariableGradient & _grad_op;

//  Real _GBEnergy;
};

#endif // GBEVOLUTIONUDRIVEN_H

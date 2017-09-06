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

#ifndef COEFDIFFUSION_H
#define COEFDIFFUSION_H

#include "Kernel.h"

class CoefDiffusion;

template <>
InputParameters validParams<CoefDiffusion>();

/**
 * This kernel implements the Laplacian operator:
 * $\nabla u \cdot \nabla \phi_i$
 */
class CoefDiffusion : public Kernel
{
public:
  CoefDiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
private:
  const MaterialProperty<Real> & _diffusion_coefficient_dimless; 
  /*
  const Real _kinetic_coefficient;
  const Real _interface_width;
  const Real _a1;
  const Real _a2;
  const Real _diffusion_coefficient;
  const Real _capillary_length;
  const Real _coupling_constant;
  const Real _tau_0;
  const Real _diffusion_coefficient_dimless;
  */
};

#endif /* COEFDIFFUSION_H */

/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "DiffusionCoefficientsMaterial.h"
#include "MathUtils.h"
#include "MooseMesh.h"

template <>
InputParameters
validParams<DiffusionCoefficientsMaterial>()
{
  //InputParameters params = validParams<GBEvolutionBaseUDriven>();
  InputParameters params = validParams<Material>();
  params.addRequiredParam<Real>("diffusion_coefficient_1","Diffusion coefficient of phase 1 in m^2/s");
  params.addRequiredParam<Real>("diffusion_coefficient_2","Diffusion coefficient of phase 2 in m^2/s");
  params.addRequiredParam<Real>("diffusion_coefficient_3","Diffusion coefficient of phase 3 in m^2/s");
  params.addRequiredCoupledVar("phi_1", "Phase field of species 1");
  params.addRequiredCoupledVar("phi_2", "Phase field of species 2");
  params.addRequiredCoupledVar("phi_3", "Phase field of species 3");
  return params;
}

DiffusionCoefficientsMaterial::DiffusionCoefficientsMaterial(const InputParameters & parameters)
  : Material(parameters), //GBEvolutionBaseUDriven(parameters), //, _GBEnergy(getParam<Real>("GBenergy"))
//    _tau_anisotropy(declareProperty<Real>("tau_anisotropy")),
    _diffusion_coefficient_1(getParam<Real>("diffusion_coefficient_1")),
    _diffusion_coefficient_2(getParam<Real>("diffusion_coefficient_2")),
    _diffusion_coefficient_3(getParam<Real>("diffusion_coefficient_3")),
	 _phi_1(coupledValue("phi_1")),
	 _phi_2(coupledValue("phi_2")),
	 _phi_3(coupledValue("phi_3")),
    _diffusion_coefficient(declareProperty<Real>("diffusion_coefficient"))
{
}

// Kinetic anisotropy should be implemented here. 
// It should depend on constant scalar material parameters such as W, beta, a1, a2, anisotropy strength epsilon_k
// It should also depend on order parameter unit normal gradient
void
DiffusionCoefficientsMaterial::computeQpProperties()
{
   _diffusion_coefficient[_qp] = _phi_1[_qp]*_diffusion_coefficient_1
                                +_phi_2[_qp]*_diffusion_coefficient_2
                                +_phi_3[_qp]*_diffusion_coefficient_3;
}

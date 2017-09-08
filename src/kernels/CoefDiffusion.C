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

#include "CoefDiffusion.h"

template <>
InputParameters
validParams<CoefDiffusion>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The Laplacian operator ($-\\nabla \\cdot \\nabla u$), with the weak "
                             "form of $(\\nabla \\phi_i, \\nabla u_h)$.");
  params.addRequiredParam<Real>("length_scale","Length scale");
  params.addRequiredParam<Real>("time_scale",  "Time scale");
  //params.addRequiredParam<std::vector<Real> >("diffusion_coefficient", "Value of the constant diffusion coefficient");
//  params.addRequiredParam<Real>("kinetic_coefficient",   "Kinetic coefficient beta");
//  params.addRequiredParam<Real>("interface_width","Numerical solid-liquid interface, usually 10-100 nm");
//  params.addRequiredParam<Real>("a1","Numerical coefficient dependant on PF interpolation function");
//  params.addRequiredParam<Real>("a2","Numerical coefficient dependant on PF interpolation function");
//  params.addRequiredParam<Real>("diffusion_coefficient","Diffusion coefficient in the solid and the liquid in m^2/s, usually around 1e-9 m^2/s");
//  params.addRequiredParam<Real>("capillary_length","capillary length");
  //params.set<Real>("diffusion_coefficient") = 1.0;
  return params;
}

CoefDiffusion::CoefDiffusion(const InputParameters & parameters) 
	: Kernel(parameters), 
    _length_scale(getParam<Real>("length_scale")),
    _time_scale(getParam<Real>("time_scale")),
    _non_dimensionalizer( _time_scale/( _length_scale*_length_scale ) ),
    _diffusion_coefficient( getMaterialProperty<Real>("diffusion_coefficient") )
	/*
    _kinetic_coefficient(getParam<Real>("kinetic_coefficient")),
    _interface_width(getParam<Real>("interface_width")),
    _a1(getParam<Real>("a1")),
    _a2(getParam<Real>("a2")),
    _diffusion_coefficient(getParam<Real>("diffusion_coefficient")),
    _capillary_length(getParam<Real>("capillary_length")),
    _coupling_constant( _a1*_interface_width/_capillary_length),
    _tau_0( _interface_width*_coupling_constant* (_kinetic_coefficient/_a1 + _a2*_interface_width/_diffusion_coefficient ) ),
    _diffusion_coefficient_dimless( _diffusion_coefficient*_tau_0/(_interface_width*_interface_width) ) 
	 */
{
    std::cout << "time_scale / length_scale^2 = " << _time_scale/(_length_scale*_length_scale) << std::endl; 
}

Real
CoefDiffusion::computeQpResidual()
{
  return _non_dimensionalizer*_diffusion_coefficient[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];
}

Real
CoefDiffusion::computeQpJacobian()
{
  return _non_dimensionalizer*_diffusion_coefficient[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}

/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ConcFromMuAux.h"

template <>
InputParameters
validParams<ConcFromMuAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Calculate location of grain boundaries in a polycrystalline sample");
  params.addRequiredParam<Real>("c_eq_1",  "equilibrium concentration for phase 1");
  params.addRequiredParam<Real>("c_eq_2",  "equilibrium concentration for phase 2");
  params.addRequiredParam<Real>("c_eq_3",  "equilibrium concentration for phase 2");
  params.addRequiredParam<Real>("G_c_1",  "dG/dc for phase 1");
  params.addRequiredParam<Real>("G_c_2",  "dG/dc for phase 2");
  params.addRequiredParam<Real>("G_c_3",  "dG/dc for phase 2");
  params.addRequiredParam<Real>("G_cc_1",  "d^2G/dc^2 for phase 1");
  params.addRequiredParam<Real>("G_cc_2",  "d^2G/dc^2 for phase 2");
  params.addRequiredParam<Real>("G_cc_3",  "d^2G/dc^2 for phase 2");
  params.addRequiredCoupledVar("phi_1", "Phase field of phase 1");
  params.addRequiredCoupledVar("phi_2", "Phase field of phase 2");
  params.addRequiredCoupledVar("mu", "Diffusion potential of chromium");
  return params;
}

ConcFromMuAux::ConcFromMuAux(const InputParameters & parameters)
  : AuxKernel(parameters), 
	_c_eq_1( getParam<Real>("c_eq_1") ),
	_c_eq_2( getParam<Real>("c_eq_2") ),
	_c_eq_3( getParam<Real>("c_eq_3") ),
	_G_c_1( getParam<Real>("G_c_1") ),
	_G_c_2( getParam<Real>("G_c_2") ),
	_G_c_3( getParam<Real>("G_c_3") ),
	_inv_G_cc_1( 1./( getParam<Real>("inv_G_cc_1") ) ),
	_inv_G_cc_2( 1./( getParam<Real>("inv_G_cc_2") ) ),
	_inv_G_cc_3( 1./( getParam<Real>("inv_G_cc_3") ) ),
	_phi_1( coupledValue("phi_1") ),
	_phi_2( coupledValue("phi_2") ),
	_mu( coupledValue("mu") )
{
}

Real
ConcFromMuAux::computeValue()
{
  // G = 0.5*G_cc*(c-c_eq)^2 + G_c*(c-c_eq) + G_eq
  // mu = dG/dc = G_cc*(c-c_eq) + G_c
  // <=> c_phase(mu) = (mu - G_c)/G_cc + c_eq for each phase concentration
  // 
  Real _phi_3 = 1 - _phi_1[_qp] - _phi_2[_qp];
  return   _phi_1[_qp]*( _inv_G_cc_1*( _mu[_qp] - _G_c_1 ) + _c_eq_1 )
         + _phi_2[_qp]*( _inv_G_cc_2*( _mu[_qp] - _G_c_2 ) + _c_eq_2 )
         + _phi_3     *( _inv_G_cc_3*( _mu[_qp] - _G_c_3 ) + _c_eq_3 );
}

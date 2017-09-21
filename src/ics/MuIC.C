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

#include "MuIC.h"

template <>
InputParameters
validParams<MuIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<Real>("energy_scale", "Energy scale in J/mol");
  params.addRequiredParam<Real>("c_init_1", "The value of the initial condition");
  params.addRequiredParam<Real>("c_init_2", "The value of the initial condition");
  params.addRequiredParam<Real>("c_init_3", "The value of the initial condition");
  params.addRequiredParam<Real>("c_eq_1", "The value of the initial condition");
  params.addRequiredParam<Real>("c_eq_2", "The value of the initial condition");
  params.addRequiredParam<Real>("c_eq_3", "The value of the initial condition");
  params.addRequiredParam<Real>("G_c_1", "The value of the initial condition");
  params.addRequiredParam<Real>("G_c_2", "The value of the initial condition");
  params.addRequiredParam<Real>("G_c_3", "The value of the initial condition");
  params.addRequiredParam<Real>("G_cc_1", "The value of the initial condition");
  params.addRequiredParam<Real>("G_cc_2", "The value of the initial condition");
  params.addRequiredParam<Real>("G_cc_3", "The value of the initial condition");
  params.addRequiredCoupledVar("phi_1", "The value of the initial condition");
  params.addRequiredCoupledVar("phi_2", "The value of the initial condition");
  params.addParam<bool>("initialize_step", false, "The value of the initial condition");
  params.addParam<Real>("step_position", 30000, "The value of the initial condition");
  return params;
}

MuIC::MuIC(const InputParameters & parameters)
  : InitialCondition(parameters), 
	 _energy_scale( getParam<Real>("energy_scale") ),
	 _c_init_1( getParam<Real>("c_init_1") ),
	 _c_init_2( getParam<Real>("c_init_2") ),
	 _c_init_3( getParam<Real>("c_init_3") ),
	 _c_eq_1( getParam<Real>("c_eq_1") ),
	 _c_eq_2( getParam<Real>("c_eq_2") ),
	 _c_eq_3( getParam<Real>("c_eq_3") ),
	 _G_c_1( getParam<Real>("G_c_1")/_energy_scale ),
	 _G_c_2( getParam<Real>("G_c_2")/_energy_scale ),
	 _G_c_3( getParam<Real>("G_c_3")/_energy_scale ),
	 _G_cc_1( getParam<Real>("G_cc_1")/_energy_scale ),
	 _G_cc_2( getParam<Real>("G_cc_2")/_energy_scale ),
	 _G_cc_3( getParam<Real>("G_cc_3")/_energy_scale ),
	 _phi_1( coupledValue("phi_1")),
	 _phi_2( coupledValue("phi_2")),
	 _initialize_step( getParam<bool>("initialize_step") ),
	 _step_position( getParam<Real>("step_position") )
{
}

Real
MuIC::value(const Point & p)
{
  /**
   * _value * x
   * The Point class is defined in libMesh.  The spatial
   * coordinates x,y,z can be accessed individually using
   * the parenthesis operator and a numeric index from 0..2
   */
  //return 2. * _coefficient * std::abs(p(0));
	// mu = G_cc*( c_init - c_eq ) + G_c
  Real mu_1 = _G_cc_1*( _c_init_1 - _c_eq_1 ) + _G_c_1; 
  Real mu_2 = _G_cc_2*( _c_init_2 - _c_eq_2 ) + _G_c_2; 
  Real mu_3 = _G_cc_3*( _c_init_3 - _c_eq_3 ) + _G_c_3; 
  if( _initialize_step )
  {
	  Real position_scaling = 0.5*( std::tanh( (p(0) - _step_position)/(0.1*_step_position) ) + 1 );
	  mu_1 = position_scaling*mu_1; 
	  mu_2 = position_scaling*mu_2; 
	  mu_3 = position_scaling*mu_3; 
  }

  Real _phi_3 = 1 - _phi_1[_qp] - _phi_2[_qp];
  return  _phi_1[_qp]*mu_1  //_G_c_1 
	     + _phi_2[_qp]*mu_2  //_G_c_2
	     + _phi_3     *mu_3; //_G_c_3;

}

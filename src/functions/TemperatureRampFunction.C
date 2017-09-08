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

#include "TemperatureRampFunction.h"

template <>
InputParameters
validParams<TemperatureRampFunction>()
{
  InputParameters params = validParams<Function>();
  params.addRequiredParam<Real>("T_init_final", "Temperature initially and finally ");
  params.addRequiredParam<Real>("T_high", "Maximum temperature ");
  params.addRequiredParam<Real>("time_increase_decrease", "Time it takes to increase and decrease the temperature");
  params.addRequiredParam<Real>("time_hold_high", "Temperature held at the maximum");
  return params;
}

TemperatureRampFunction::TemperatureRampFunction(const InputParameters & parameters)
  : Function(parameters), 
	_T_init_final(getParam<Real>("T_init_final")),
	_T_high(getParam<Real>("T_high")),
	_time_increase_decrease(getParam<Real>("time_increase_decrease")),
	_time_hold_high(getParam<Real>("time_hold_high")),
   _T_slope( (_T_high - _T_init_final)/_time_increase_decrease)
{
}

Real
TemperatureRampFunction::value(Real t, const Point & /*p*/)
{
   if (t < _time_increase_decrease)
		return _T_init_final + _T_slope*t; 
	else if (t < _time_increase_decrease + _time_hold_high)
		return _T_high;
	else 
	{
		Real temp_decrease = _T_high - _T_slope*(t - _time_increase_decrease - _time_hold_high);
		if ( temp_decrease > _T_init_final )
			return temp_decrease;
		else
			return _T_init_final; // do not go below the original temperature 
	}
}

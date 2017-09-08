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

#ifndef TEMPERATURERAMPFUNCTION_H
#define TEMPERATURERAMPFUNCTION_H

#include "Function.h"

class TemperatureRampFunction;

template <>
InputParameters validParams<TemperatureRampFunction>();

class TemperatureRampFunction : public Function
{
public:
  TemperatureRampFunction(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p) override;

protected:
  Real _T_init_final;
  Real _T_high;
  Real _time_increase_decrease;
  Real _time_hold_high;
  Real _T_slope; 
};

#endif // EXAMPLEFUNCTION_H

/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "DiffusionCoefficientsMaterial.h"
#include "Function.h"
#include "MathUtils.h"
#include "MooseMesh.h"

template <>
InputParameters
validParams<DiffusionCoefficientsMaterial>()
{
  //InputParameters params = validParams<GBEvolutionBaseUDriven>();
  InputParameters params = validParams<Material>();
  params.addRequiredParam<Real>("R","Molar gas constant");
  params.addRequiredParam<Real>("length_scale","Length scale");
  params.addRequiredParam<Real>("time_scale",  "Time scale");
  params.addRequiredParam<Real>("D_1_pre","Arrhenius diffusion precoefficient of phase 1 in m^2/s");
  params.addRequiredParam<Real>("D_2_pre","Arrhenius diffusion precoefficient of phase 2 in m^2/s");
  params.addRequiredParam<Real>("D_3_pre","Arrhenius diffusion precoefficient of phase 3 in m^2/s");
  params.addRequiredParam<Real>("Q_1","Activation energy for phase 1");
  params.addRequiredParam<Real>("Q_2","Activation energy for phase 2");
  params.addRequiredParam<Real>("Q_3","Activation energy for phase 3");
  params.addRequiredParam<FunctionName>("TemperatureRampFunction","Temperature ramp function name");
  params.addRequiredCoupledVar("phi_1", "Phase field of species 1");
  params.addRequiredCoupledVar("phi_2", "Phase field of species 2");
  return params;
}

DiffusionCoefficientsMaterial::DiffusionCoefficientsMaterial(const InputParameters & parameters)
  : Material(parameters), //GBEvolutionBaseUDriven(parameters), //, _GBEnergy(getParam<Real>("GBenergy"))
//    _tau_anisotropy(declareProperty<Real>("tau_anisotropy")),
	 _R(getParam<Real>("R")),
	 _length_scale(getParam<Real>("length_scale")),
	 _time_scale(getParam<Real>("time_scale")),
	 _non_dimensionalizer( _time_scale/( _length_scale*_length_scale ) ),
    _D_1_pre( getParam<Real>("D_1_pre")),
    _D_2_pre( getParam<Real>("D_2_pre")),
    _D_3_pre( getParam<Real>("D_3_pre")),
    _Q_1(getParam<Real>("Q_1")),
    _Q_2(getParam<Real>("Q_2")),
    _Q_3(getParam<Real>("Q_3")),
	 _phi_1(coupledValue("phi_1")),
	 _phi_2(coupledValue("phi_2")),
	 _TemperatureRampFunction( getFunction("TemperatureRampFunction") ),
    _diffusion_coefficient(declareProperty<Real>("diffusion_coefficient"))
{
	std::cout << "_non_dimensionalizer = " << _non_dimensionalizer << std::endl; 

	int current_node = 0;

	Real temperature = _TemperatureRampFunction.value( _t, current_node ); //, *_current_node ); 
	Real inv_RT = 1./(_R*temperature); 
   std::cout << "initial D_1 = " << _D_1_pre*std::exp( -_Q_1*inv_RT ) << std::endl; 
   std::cout << "initial D_2 = " << _D_2_pre*std::exp( -_Q_2*inv_RT ) << std::endl; 

	temperature = _TemperatureRampFunction.value( _t, current_node ); //, *_current_node ); 
	inv_RT = 1./(_R*temperature); 
   std::cout << "initial D_1 = " << _D_1_pre*std::exp( -_Q_1*inv_RT ) << std::endl; 
   std::cout << "initial D_2 = " << _D_2_pre*std::exp( -_Q_2*inv_RT ) << std::endl; 
}

// Kinetic anisotropy should be implemented here. 
// It should depend on constant scalar material parameters such as W, beta, a1, a2, anisotropy strength epsilon_k
// It should also depend on order parameter unit normal gradient
void
DiffusionCoefficientsMaterial::computeQpProperties()
{
	int current_node = 0;
	const Real temperature = _TemperatureRampFunction.value( _t, current_node ); //, *_current_node ); 
	const Real inv_RT = 1./(_R*temperature); 
   _diffusion_coefficient[_qp] = (
				                                    _phi_1[_qp] *_D_1_pre*std::exp( -_Q_1*inv_RT )
                                +               _phi_2[_qp] *_D_2_pre*std::exp( -_Q_2*inv_RT )
                                +(1-_phi_1[_qp]-_phi_2[_qp])*_D_3_pre*std::exp( -_Q_3*inv_RT ) );
}

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


  params.addRequiredParam<Real>("c_eq_1",  "equilibrium concentration for phase 1");
  params.addRequiredParam<Real>("c_eq_2",  "equilibrium concentration for phase 2");
  params.addRequiredParam<Real>("c_eq_3",  "equilibrium concentration for phase 3");
  params.addRequiredParam<Real>("energy_scale",  "Energy scale in J/mol");
  params.addRequiredParam<Real>("G_c_1",  "dG/dc for phase 1");
  params.addRequiredParam<Real>("G_c_2",  "dG/dc for phase 2");
  params.addRequiredParam<Real>("G_c_3",  "dG/dc for phase 3");
  params.addRequiredParam<Real>("G_cc_1",  "d^2G/dc^2 for phase 1");
  params.addRequiredParam<Real>("G_cc_2",  "d^2G/dc^2 for phase 2");
  params.addRequiredParam<Real>("G_cc_3",  "d^2G/dc^2 for phase 3");
  params.addRequiredParam<Real>("G_ccc_1",  "d^3G/dc^3 for phase 1");
  params.addRequiredParam<Real>("G_ccc_2",  "d^3G/dc^3 for phase 2");
  params.addRequiredParam<Real>("G_ccc_3",  "d^3G/dc^3 for phase 3");
  params.addRequiredCoupledVar("mu", "Diffusion potential of chromium");
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
    _diffusion_coefficient(declareProperty<Real>("diffusion_coefficient")),

    _c_eq_1( getParam<Real>("c_eq_1") ),
    _c_eq_2( getParam<Real>("c_eq_2") ),
    _c_eq_3( getParam<Real>("c_eq_3") ),

	 _energy_scale( getParam<Real>("energy_scale") ),
    _G_c_1( getParam<Real>("G_c_1")/_energy_scale ),
    _G_c_2( getParam<Real>("G_c_2")/_energy_scale ),
    _G_c_3( getParam<Real>("G_c_3")/_energy_scale ),
    _G_cc_1( getParam<Real>("G_cc_1")/_energy_scale ),
    _G_cc_2( getParam<Real>("G_cc_2")/_energy_scale ),
    _G_cc_3( getParam<Real>("G_cc_3")/_energy_scale ),
    _inv_G_cc_1( 1./( getParam<Real>("G_cc_1")/_energy_scale ) ),
    _inv_G_cc_2( 1./( getParam<Real>("G_cc_2")/_energy_scale ) ),
    _inv_G_cc_3( 1./( getParam<Real>("G_cc_3")/_energy_scale ) ),
    _G_ccc_1( getParam<Real>("G_ccc_1")/_energy_scale ),
    _G_ccc_2( getParam<Real>("G_ccc_2")/_energy_scale ),
    _G_ccc_3( getParam<Real>("G_ccc_3")/_energy_scale ),
    _mu( coupledValue("mu") ),
    _concFromMu(declareProperty<Real>("concFromMu")),
    _mobility(declareProperty<Real>("mobility"))

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
   Real _phi_3 = 1 - _phi_1[_qp] - _phi_2[_qp];
   _diffusion_coefficient[_qp] = (  _phi_1[_qp]*_D_1_pre*std::exp( -_Q_1*inv_RT )
                                  + _phi_2[_qp]*_D_2_pre*std::exp( -_Q_2*inv_RT )
                                  + _phi_3     *_D_3_pre*std::exp( -_Q_3*inv_RT ) );

// G    = 1/6*G_ccc*(c-c_eq)^3 + 0.5*G_cc*(c-c_eq)^2 + G_c*(c-c_eq) + G_eq
// mu   = 0.5*G_ccc*(c-c_eq)^2 + G_cc*(c-c_eq) + G_c
// curv =     G_ccc*(c-c_eq) + G_cc
// 0.5*G_ccc*(c-c_eq)^2 + G_cc*(c-c_eq) + G_c - mu = 0
// 0.5*G_ccc*(c^2 - 2*c*c_eq + c_eq^2) + G_cc*c - G_cc*c_eq + G_c - mu = 0
// 0.5*G_ccc*c^2 + (-0.5*G_ccc*2*c_eq + G_cc)*c + 0.5*G_ccc*c_eq^2 + G_cc*c_eq + G_c - mu = 0
// 0.5*G_ccc*c^2 + (-G_ccc*c_eq + G_cc)*c + 0.5*G_ccc*c_eq^2 - G_cc*c_eq + G_c - mu = 0
//         A*c^2 +                    B*c                                    +    C = 0

// Forget the positive branch, because A = 0.5*G_ccc seems to be always negative 
// c = ( -B +/- sqrt( B^2 - 4*A*C ) )/( 2*A )

// f(C) = sqrt(B^2 - 4*A*C)
// f(C) ~ f(A=0) + f'(C=0)*C
// f(C=0) = |B|
// f'(C) = 0.5*1./sqrt( B^2 - 4*A*C) * 4*A
// f'(0) = 2*A/|B|
// f(C) ~ (-B + 2*A
// c ~ -0.5*B/A + 
// if A << 1 --> sqrt( B^2 - 4*A*C ) ~ 
	Real A, B, C, Bsq_4AC, conc_1, conc_2, conc_3, inv_curv_1, inv_curv_2, inv_curv_3;

	if ( _G_ccc_1*_G_ccc_1 > 1e-5 )
	{
		A = 0.5*_G_ccc_1;
		B = -_G_ccc_1*_c_eq_1 + _G_cc_1;
		C = 0.5*_G_ccc_1*_c_eq_1*_c_eq_1 - _G_cc_1*_c_eq_1 +  _G_c_1 - _mu[_qp]; 
		Bsq_4AC = B*B - 4*A*C; 
		if (Bsq_4AC < 0)
			std::cout << "Bsq_4AC_1 = " << Bsq_4AC << " < 0!!!" << std::endl;

		conc_1 = 1./(2*A)*( -B + std::sqrt( Bsq_4AC ) ); 
		inv_curv_1 = 1./(_G_ccc_1*( _mu[_qp] - _G_c_1 ) + _G_cc_1); 
	}
	else 
	{
		conc_1 = _inv_G_cc_1*( _mu[_qp] - _G_c_1 ) + _c_eq_1;
		inv_curv_1 = _inv_G_cc_1; 
	}

	if ( _G_ccc_2*_G_ccc_2 > 1e-5 )
	{
		A = 0.5*_G_ccc_2;
		B =-_G_ccc_2*_c_eq_2 + _G_cc_2;
		C = 0.5*_G_ccc_2*_c_eq_2*_c_eq_2 -_G_cc_2*_c_eq_2 +  _G_c_2 - _mu[_qp]; 
		Bsq_4AC = B*B - 4*A*C; 
		if (Bsq_4AC < 0)
			std::cout << "Bsq_4AC_2 = " << Bsq_4AC << " < 0!!!" << std::endl;

		conc_2 = 1./(2*A)*( -B + std::sqrt( Bsq_4AC ) ); 
		inv_curv_2 = 1./(_G_ccc_2*( _mu[_qp] - _G_c_2 ) + _G_cc_2); 
	}
	else 
	{
		conc_2 = _inv_G_cc_2*( _mu[_qp] - _G_c_2 ) + _c_eq_2;
		inv_curv_2 = _inv_G_cc_2; 
	}


	if ( _G_ccc_3*_G_ccc_3 > 1e-5 )
	{
		A = 0.5*_G_ccc_3;
		B =-_G_ccc_3*_c_eq_3 + _G_cc_3;
		C = 0.5*_G_ccc_3*_c_eq_3*_c_eq_3 - _G_cc_3*_c_eq_3 +  _G_c_3 - _mu[_qp]; 
		Bsq_4AC = B*B - 4*A*C; 
		if (Bsq_4AC < 0)
			std::cout << "Bsq_4AC_2 = " << Bsq_4AC << " < 0!!!" << std::endl;

		conc_3 = 1./(2*A)*( -B + std::sqrt( Bsq_4AC ) ); 
		inv_curv_3 = 1./(_G_ccc_3*( _mu[_qp] - _G_c_3 ) + _G_cc_3); 
	}
	else 
	{
		conc_3 = _inv_G_cc_3*( _mu[_qp] - _G_c_3 ) + _c_eq_3;
		inv_curv_3 = _inv_G_cc_3; 
	}

   _concFromMu[_qp] =  _phi_1[_qp]*conc_1  // ( _inv_G_cc_1*( _mu[_qp] - _G_c_1 ) + _c_eq_1 )
                     + _phi_2[_qp]*conc_2  // ( _inv_G_cc_2*( _mu[_qp] - _G_c_2 ) + _c_eq_2 )
                     + _phi_3     *conc_3; // ( _inv_G_cc_3*( _mu[_qp] - _G_c_3 ) + _c_eq_3 );


   _mobility[_qp] = (  _phi_1[_qp]*_D_1_pre*std::exp( -_Q_1*inv_RT )*inv_curv_1   ///G_cc_1
                     + _phi_2[_qp]*_D_2_pre*std::exp( -_Q_2*inv_RT )*inv_curv_2   //inv_G_cc_2
                     + _phi_3     *_D_3_pre*std::exp( -_Q_3*inv_RT )*inv_curv_3); //G_cc_3 );
}

#include "diffusionBondingApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

#include "CoefDiffusion.h"
#include "DiffusionCoefficientsMaterial.h"

#include "TemperatureRampFunction.h"
#include "MuDiffusion.h"
//#include "ConcFromMuAux.h"

template <>
InputParameters
validParams<diffusionBondingApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

diffusionBondingApp::diffusionBondingApp(InputParameters parameters) : MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  diffusionBondingApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  diffusionBondingApp::associateSyntax(_syntax, _action_factory);
}

diffusionBondingApp::~diffusionBondingApp() {}

// External entry point for dynamic application loading
extern "C" void
diffusionBondingApp__registerApps()
{
  diffusionBondingApp::registerApps();
}
void
diffusionBondingApp::registerApps()
{
  registerApp(diffusionBondingApp);
}

// External entry point for dynamic object registration
extern "C" void
diffusionBondingApp__registerObjects(Factory & factory)
{
  diffusionBondingApp::registerObjects(factory);
}
void
diffusionBondingApp::registerObjects(Factory & factory)
{
	registerKernel(CoefDiffusion); 
	registerKernel(MuDiffusion); 
//	registerAux(ConcFromMuAux); 
	registerMaterial(DiffusionCoefficientsMaterial); 
	registerFunction(TemperatureRampFunction); 
}

// External entry point for dynamic syntax association
extern "C" void
diffusionBondingApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  diffusionBondingApp::associateSyntax(syntax, action_factory);
}
void
diffusionBondingApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}

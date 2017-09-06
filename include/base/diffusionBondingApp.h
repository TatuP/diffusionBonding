#ifndef DIFFUSION_BONDINGAPP_H
#define DIFFUSION_BONDINGAPP_H

#include "MooseApp.h"

class diffusionBondingApp;

template <>
InputParameters validParams<diffusionBondingApp>();

class diffusionBondingApp : public MooseApp
{
public:
  diffusionBondingApp(InputParameters parameters);
  virtual ~diffusionBondingApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* DIFFUSION_BONDINGAPP_H */

# Mesh adaptivity and time step adaptivity are used
# Diffusion coefficient for 3 phases in Arrhenius form
# D_i(T) = D_0,i * exp( -Q_i/(R*T) )

[GlobalParams]
  R = 8.3144598
  length_scale = 10e-9 # 0.5 microns is the unit length
  time_scale = 60 # 60 seconds

  T_init_final = 800 # Room temperature
  T_high = 1223 # 273 + 960 K 
  time_increase_decrease = 180 # min = 3*60 min = 3 h
  time_hold_high = 240 # min = 4 h
# Used to non-dimensionalize diffusion coefficient
# D = D_dimless * length_scale^2 / time_scale 
# D_dimless = D * time_scale / length_scale^2
# Solid diffusion ~ 1e-15 m^2/s
# D_dimless ~ 1e-15 * 60 * 1e16 ~ 600 

  energy_scale = 1e4 # Chemical potential and G curvature given in units of energy_scale 

# steel
  
  c_init_1 = 0.0084726

  c_eq_1 = 0.0084726 # steel equilibrium mole fraction
  G_c_1  = -34749.7639282 # J/mol
  G_cc_1 = 2594777.04155 # steel Gibbs free energy curvature
#  G_cc_1_extrap = 1.8327e-05 # Approaching a miscibility gap
#  c_3_extrap = 0.013763
  G_ccc_1 = 0 # -490468970.5 

# copper

  c_init_2 = 0.00976 # corresponding to 0.8 wt%

#  c_eq_2 = 0.0038
#  G_c_2  = 30541.6447 # J/mol
#  G_cc_2         =  2696551.0 # J/mol

# Equilibrium with 0.8 wt%
#  c_eq_2 = 0.00976 # corresponding to 0.8 wt%
#  G_c_2 = 50114.627803 # 
#  G_cc_2 = 2482503.10 # J/mol

# Equilibrium with 0.1 wt%
  c_eq_2 = 0.0012218 # mol frac
  G_c_2  = 35335.6754 # J/mol
  G_cc_2 = 8390770.76 # J/mol
  
# Equilibrium with 0.3 wt%
#  c_eq_2 = 0.00366375 # mol frac
#  G_c_2  = 46581.7361 # J/mol
#  G_cc_2 = 2798165.82 # J/mol

  G_ccc_2 =  0

  c_init_3 = 0.1
  c_eq_3 = 0.0138
  G_c_3  = 30541.6447 # J/mol
  G_cc_3 = 4696551.0 # J/mol
  G_ccc_3 = 0 # -490468970.5 
#  G_cc_3_extrap = 4696551.0 # J/mol
#  c_3_extrap = 0

# G    = 1/6*G_ccc*(c-c_eq)^3 + 0.5*G_cc*(c-c_eq)^2 + G_c*(c-c_eq) + G_eq
# mu   = 0.5*G_ccc*(c-c_eq)^2 + G_cc*(c-c_eq) + G_c
# curv =     G_ccc*(c-c_eq) + G_cc 
# curv( c=c_extrap ) = G_cc_extrap 
# G_ccc*( c_extrap-c_eq ) + G_cc = G_cc_extrap 
# G_ccc = ( G_cc_extrap - G_cc )/( c_extrap - c_eq )

[]

[Mesh]
  # Mesh block.  Meshes can be read in or automatically generated
  type = GeneratedMesh
  dim = 1 # Problem dimension
  nx = 2000 # Number of elements in the x-direction
  xmin = 0    # minimum x-coordinate of the mesh
  xmax = 60000 # maximum x-coordinate of the mesh. In units of length_scale 
               # 6e4*10e-9 = 6e5*10e-9 = 600e-6 m = 0.6 mm
  elem_type = EDGE  # Type of elements used in the mesh
  uniform_refine = 4 # Initial uniformly conducted mesh adaption. 
                     # uniform_refine = max_h_level means full initial refinement
[]
# Initial grid spacing in W units: 110/15
# if 4 levels of adaption are allowed, 
# the maximum refinement spacing is 110/(2^4*15) = 110/(16*15) = 110/240  

[ICs]
  [./steelIC]
    type = SmoothCircleIC
    variable = steel
    int_width = 100
    x1 = 0
    y1 = 0
    radius = 30000
    outvalue = 0 #0
    invalue  = 1
  [../]

  [./copperIC]
    type = SmoothCircleIC
    variable = copper
    int_width = 100
    x1 = 0
    y1 = 0
    radius = 30000
    outvalue = 1 #0
    invalue  = 0
  [../]

#  [./barrierIC]
#    type = ConstantIC
#    variable = barrier
#    value = 0
#  [../]

  [./muIC]
    type = MuIC
    variable = chemPot 
    phi_1 = steel
    phi_2 = copper
  [../]
[]

[Variables]
  [./chemPot]
  [../]
[]

[AuxVariables]
  [./steel]
  [../]
  [./copper]
  [../]
[]


[BCs]
  [./ZeroFlux_left]
    type = NeumannBC
    boundary = left
    variable = chemPot
    value = 0
  [../]

  [./ZeroFlux_right]
    type = NeumannBC
    boundary = right
    variable = chemPot
    value = 0
  [../]

[]

[Kernels]
  # Kernel block, where the kernels defining the residual equations are set up.
  [./mu_dot]
    type = TimeDerivative
    variable = chemPot
  [../]

  [./mu_diffusion]
    type = MuDiffusion
    variable = chemPot
  [../] 

[]

[Functions]
  [./TemperatureRamp]
    type = TemperatureRampFunction    
  [../]
[]

[Materials]

#  [./diffusion_coefficient_parser]
#    type = ParsedMaterial
#    f_name = diffusion_coefficient
#    args = 'steel copper barrier'
#    constant_names       = 'D_copper D_barrier length_scale time_scale'
#    constant_expressions = '1e-15 	
#									 10e-15 	
#									 0.1e-15
#									 10e-9
#                            60
#									' 
#    function =steel*D_steel+copper*D_copper+(1-steel-copper)*D_barrier
#    outputs = exodus
#  [../]

  [./diffusion_coefficients_material]
     type = DiffusionCoefficientsMaterial
     phi_1 = steel     
     phi_2 = copper
     D_1_pre = 3.59e-6
     Q_1 = 179e3  # J/mol
     D_2_pre = 0.337e-4
     Q_2 = 195e3
     D_3_pre = 3.59e-7
     Q_3 = 179e3
     TemperatureRampFunction = TemperatureRamp
     mu = chemPot
     outputs = exodus
  [../]

[]

[Postprocessors]
  # Scalar postprocessors
  [./dt]
    # Outputs the current time step
    type = TimestepSize
  [../]
  [./dofs]
    # Outputs the current time step
    type = NumDOFs
  [../]
[]

[Executioner]
  type = Transient # Type of executioner, here it is transient with an adaptive time step
#  scheme = bdf2 # Type of time integration (2nd order backward euler), 
  scheme = explicit-euler
                # defaults to 1st order backward euler

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK' # preconditioned jacobian free newton krylov
#  solve_type = 'JFNK' # not preconditioned. Much worse convergence

  # Uses newton iteration to solve the problem.
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -mat_mffd_type'
  petsc_options_value = 'hypre boomeramg 101 ds'

  l_max_its = 30 # Max number of linear iterations
  l_tol = 0 # 1e-8 # Relative tolerance for linear solves
  nl_max_its = 40 # Max number of nonlinear iterations
  nl_rel_tol = 1e-9 # 1e-10 # Absolute tolerance for nonlienar solves
 
  start_time = 0.0
  end_time = 600 # minutes = 10 h

  [./TimeStepper]
    type = IterationAdaptiveDT
#    type = IterationAdaptiveDT_maxDT # added a maximum DT
    dt = 0.01 # Initial time step.  In this simulation it changes.
    optimal_iterations = 6 # Time step will adapt 
                           # to maintain this number of nonlinear iterations
    growth_factor = 1.2
    cutback_factor = 0.5
#    max_dt = 150.0
  [../]

[]

[Adaptivity]
  max_h_level = 6 # Max number of refinements used, starting from initial mesh 
                  # (before uniform refinement)
                  # Maximum refinement: initial mesh spacing/
  marker = errorfrac
  steps = 2
  [./Indicators]
    [./error]
      type = GradientJumpIndicator
      variable = chemPot # convected
    [../]
  [../]
  [./Markers]
    [./errorfrac]
      type = ErrorToleranceMarker
      indicator = error
      refine = 0.00002
      coarsen = 0.00001
    [../]
  [../] # Markers
[] # Adaptivity

[Outputs]
  interval = 2
  outputs = all
  exodus = true # Exodus file will be outputted
  csv = true
  [./console]
    type = Console
    max_rows = 20 # Will print the 20 most recent postprocessor values to the screen
  [../]
[]

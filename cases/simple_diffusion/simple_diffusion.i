# Mesh adaptivity and time step adaptivity are used
# Diffusion coefficient for 3 phases in Arrhenius form
# D_i(T) = D_0,i * exp( -Q_i/(R*T) )

[GlobalParams]
  R = 8.3144598
  length_scale = 0.1e-6 # 0.1 microns
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
[]

[Mesh]
  # Mesh block.  Meshes can be read in or automatically generated
  type = GeneratedMesh
  dim = 1 # Problem dimension
  nx = 100 # Number of elements in the x-direction
  xmin = 0    # minimum x-coordinate of the mesh
  xmax = 1000 # maximum x-coordinate of the mesh. In units of length_scale 
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
    int_width = 5
    x1 = 0
    y1 = 0
    radius = 500
    outvalue = 0 #0
    invalue  = 1
  [../]

  [./copperIC]
    type = SmoothCircleIC
    variable = copper
    int_width = 5
    x1 = 0
    y1 = 0
    radius = 500
    outvalue = 1 #0
    invalue  = 0
  [../]

#  [./barrierIC]
#    type = ConstantIC
#    variable = barrier
#    value = 0
#  [../]

  [./chromiumIC]
    type = SmoothCircleIC
    variable = chromium 
    int_width = 1
    x1 = 250
    y1 = 0
    radius = 250
    outvalue = 0 #0
    invalue  = 0.5
  [../]
[]

[Variables]
  [./chromium]
  [../]
[]

[AuxVariables]
  [./steel]
  [../]
  [./copper]
  [../]
[]


[BCs]
  # Boundary Condition block
   [
  [./Periodic]
    [./top_bottom]
      auto_direction = 'x' # Makes problem periodic in the x and y directions
    [../]
  [../]
#  [./chromium_left]
#      type = DirichletBC
#      variable = chromium
#      boundary = left
#      value = 0.5
#  [../]
#
#  [./chromium_right]
#      type = DirichletBC
#      variable = chromium
#      boundary = right
#      value = 0
#  [../]

[]

[Kernels]
  # Kernel block, where the kernels defining the residual equations are set up.
  [./chromium_dot]
    type = TimeDerivative
    variable = chromium
  [../]

# dU/dt = diffusion_coefficient*grad^2(U) + d(op_sum)/dt
  [./chromium_diffusion]
    type = CoefDiffusion # Custom diffusion laplacian kernel, with a predefined coefficient. 
    variable = chromium
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
  l_tol = 1e-4 # Relative tolerance for linear solves
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
    growth_factor = 1.05
    cutback_factor = 0.8
#    max_dt = 150.0
  [../]

[]

[Adaptivity]
  max_h_level = 5 # Max number of refinements used, starting from initial mesh 
                  # (before uniform refinement)
                  # Maximum refinement: initial mesh spacing/
  marker = errorfrac
  steps = 2
  [./Indicators]
    [./error]
      type = GradientJumpIndicator
      variable = chromium # convected
    [../]
  [../]
  [./Markers]
    [./errorfrac]
      type = ErrorToleranceMarker
      indicator = error
      refine = 0.0002
      coarsen = 0.0001
    [../]
  [../] # Markers
[] # Adaptivity

[Outputs]
  interval = 5
  outputs = all
  exodus = true # Exodus file will be outputted
  csv = true
  [./console]
    type = Console
    max_rows = 20 # Will print the 20 most recent postprocessor values to the screen
  [../]
[]

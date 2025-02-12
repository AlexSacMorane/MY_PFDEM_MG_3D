[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 99
  ny = 49
  nz = 99
  xmin = -220.00000000000003
  xmax = 220.00000000000003
  ymin = -110.00000000000001
  ymax = 110.00000000000001
  zmin = -220.00000000000003
  zmax = 220.00000000000003
  elem_type = HEX8
[]

[Variables]
	[./eta1]
		order = FIRST
		family = LAGRANGE
		[./InitialCondition]
			type = FunctionIC
			function = eta1_txt
		[../]
	[../]
	[./eta2]
		order = FIRST
		family = LAGRANGE
		[./InitialCondition]
			type = FunctionIC
			function = eta2_txt
		[../]
	[../]
	[./eta3]
		order = FIRST
		family = LAGRANGE
		[./InitialCondition]
			type = FunctionIC
			function = eta3_txt
		[../]
	[../]
	[./eta4]
		order = FIRST
		family = LAGRANGE
		[./InitialCondition]
			type = FunctionIC
			function = eta4_txt
		[../]
	[../]
  [./c]
    [./InitialCondition]
      type = FunctionIC
      function = c_txt
    [../]
  [../]
[]

[Kernels]
	# Order parameter eta1
	[./deta1dt]
		type = TimeDerivative
		variable = eta1
	[../]
	[./ACBulk1]
		type = AllenCahn
		variable = eta1
		mob_name = L
		f_name = F_total
	[../]
	[./ACInterface1]
		type = ACInterface
		variable = eta1
		mob_name = L
		kappa_name = 'kappa_eta'
	[../]
	# Order parameter eta2
	[./deta2dt]
		type = TimeDerivative
		variable = eta2
	[../]
	[./ACBulk2]
		type = AllenCahn
		variable = eta2
		mob_name = L
		f_name = F_total
	[../]
	[./ACInterface2]
		type = ACInterface
		variable = eta2
		mob_name = L
		kappa_name = 'kappa_eta'
	[../]
	# Order parameter eta3
	[./deta3dt]
		type = TimeDerivative
		variable = eta3
	[../]
	[./ACBulk3]
		type = AllenCahn
		variable = eta3
		mob_name = L
		f_name = F_total
	[../]
	[./ACInterface3]
		type = ACInterface
		variable = eta3
		mob_name = L
		kappa_name = 'kappa_eta'
	[../]
	# Order parameter eta4
	[./deta4dt]
		type = TimeDerivative
		variable = eta4
	[../]
	[./ACBulk4]
		type = AllenCahn
		variable = eta4
		mob_name = L
		f_name = F_total
	[../]
	[./ACInterface4]
		type = ACInterface
		variable = eta4
		mob_name = L
		kappa_name = 'kappa_eta'
	[../]
  # Order parameter c
  [./dcdt]
    type = TimeDerivative
    variable = c
  [../]
	[./eta1_c]
		type = CoefCoupledTimeDerivative
		v = 'eta1'
		variable = c
		coef = 62.266500622665
	[../]
	[./eta2_c]
		type = CoefCoupledTimeDerivative
		v = 'eta2'
		variable = c
		coef = 62.266500622665
	[../]
	[./eta3_c]
		type = CoefCoupledTimeDerivative
		v = 'eta3'
		variable = c
		coef = 62.266500622665
	[../]
	[./eta4_c]
		type = CoefCoupledTimeDerivative
		v = 'eta4'
		variable = c
		coef = 62.266500622665
	[../]
  [./c_diffusion]
    type = ACInterface
    kappa_name = kc
    mob_name = L_c
    variable = c
  [../]
[]

[Materials]
  [./consts]
    type = GenericConstantMaterial
    prop_names  = 'L kappa_eta L_c'
    prop_values ='100.0 72.61225443079469 1'
  [../]
  [./kcmap]
    type = GenericFunctionMaterial
    block = 0
    prop_names = kc
    prop_values = kc_txt
    outputs = exodus
  [../]
  [./asmap]
    type = GenericFunctionMaterial
    block = 0
    prop_names = as
    prop_values = as_txt
    outputs = exodus
  [../]
  [./free_energy_etai]
    type = DerivativeParsedMaterial
    block = 0
    property_name = F
    coupled_variables = 'eta1 eta2 eta3 eta4'
    constant_names = 'h'
    constant_expressions = '1'
    expression = 'h*(eta1^2*(1-eta1)^2)+h*(eta2^2*(1-eta2)^2)+h*(eta3^2*(1-eta3)^2)+h*(eta4^2*(1-eta4)^2)'
    enable_jit = true
    derivative_order = 2
    # outputs = exodus
  [../]
  [./Ed]
  type = DerivativeParsedMaterial
    block = 0
    property_name = Ed
    coupled_variables = 'eta1 eta2 eta3 eta4 c'
    material_property_names = 'as'
    constant_names = 'c_eq k_diss k_prec'
    constant_expressions ='1.0 0.001121186440677966 0.001121186440677966'
    expression = 'if(c<c_eq*as,k_diss,k_prec)*as*(1-c/(c_eq*as))*(3*eta1^2-2*eta1^3+3*eta2^2-2*eta2^3+3*eta3^2-2*eta3^3+3*eta4^2-2*eta4^3)'
    enable_jit = true
    derivative_order = 2
    # outputs = exodus
  [../]
  [./free_energy_and_ed]
    type = DerivativeParsedMaterial
    block = 0
    property_name = F_total
    coupled_variables = 'eta1 eta2 eta3 eta4 c'
    material_property_names = 'F(eta1,eta2,eta3,eta4) Ed(eta1,eta2,eta3,eta4,c)'
    expression = 'F+Ed'
    enable_jit = true
    derivative_order = 2
    # outputs = exodus
  [../]
[]

[Functions]
	[eta1_txt]
		type = PiecewiseMultilinear
		data_file = data/eta_1.txt
	[]
	[eta2_txt]
		type = PiecewiseMultilinear
		data_file = data/eta_2.txt
	[]
	[eta3_txt]
		type = PiecewiseMultilinear
		data_file = data/eta_3.txt
	[]
	[eta4_txt]
		type = PiecewiseMultilinear
		data_file = data/eta_4.txt
	[]
  [c_txt]
    type = PiecewiseMultilinear
    data_file = data/c.txt
  []
	[as_txt]
		type = PiecewiseMultilinear
		data_file = data/as.txt
	[]
	[kc_txt]
		type = PiecewiseMultilinear
		data_file = data/kc.txt
	[]
[]

[Preconditioning]
  # This preconditioner makes sure the Jacobian Matrix is fully populated. Our
  # kernels compute all Jacobian matrix entries.
  # This allows us to use the Newton solver below.
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  scheme = 'bdf2'

  # Automatic differentiation provides a _full_ Jacobian in this example
  # so we can safely use NEWTON for a fast solve
  solve_type = 'NEWTON'

  l_max_its = 20
  l_tol =  0.001
  l_abs_tol =  0.001
  nl_max_its = 10
  nl_rel_tol =  0.001
  nl_abs_tol =  0.001

  start_time = 0.0
  end_time   = 0.009999999999999998

  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 0.009999999999999998
  [../]
[]

[Outputs]
  execute_on = 'initial timestep_end'
  exodus = true
  [./other]
    type = VTK
    execute_on = 'initial final'
  [../]
  [console]
    type = Console
    execute_on = 'nonlinear'
    all_variable_norms = true
  []
[]

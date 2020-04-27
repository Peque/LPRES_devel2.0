--------------------------------------------------------------------------------
-- EcosimPro Simulation Source Code
-- FILE NAME: LPRES_CompHeat.el
-- DESCRIPTION: Defines components in which there are heat exchanges
-- NOTES:
-- AUTHOR: Pablo Sierra Heras
-- CREATION DATE: 19/02/2015
--------------------------------------------------------------------------------

-- Libraries
USE MATH VERSION "3.1"


-- Components
--------------------------------------------------------------------------------
-- Component that represents a cooling jacket
--------------------------------------------------------------------------------
COMPONENT CoolingJacket (ENUM Type_Cooling Type = Darcy)

  "Cooling jacket"

   PORTS
      IN Fluid l							"Inlet fluid port"
      OUT Fluid g							"Outlet fluid port"
		IN Heat h							"Inlet heat exchange port"
		
	DATA
		REAL L = 1.				UNITS u_m 			"Length of the channels" 
		REAL a = 0.002			UNITS u_m 			"Average width of the channels" 
		REAL b = 0.004			UNITS u_m 			"Average height of the channels" 
		INTEGER N = 100		UNITS no_units		"Number of channels"
		REAL rug = 5.e-5		UNITS u_m			"Absolute rugosity"
		REAL k_w = 370			UNITS	u_W_mK		"Thermal conductivity of the nozzle wall material"
		REAL t_w = 0.003		UNITS u_m 			"Nozzle wall thickness" 
		REAL dp = 1.e5			UNITS u_Pa			"Imposed total pressure drop if Type=Known_dp"
		
	DECLS
		REAL Q							UNITS u_W 				"Heat flux"
		REAL T_w_hot					UNITS u_K				"Wall temperature of the combustion gases side"				
		REAL T_w_cold					UNITS u_K				"Wall temperature of the cooling liquid side"	
		REAL A_wet_cooling			UNITS u_m2				"Cooling wet area"
		REAL A_wet_nozzle				UNITS u_m2				"Nozzle wet area"
		REAL h_l							UNITS u_W_m2K			"Cooling liquid heat transfer coefficient"
		REAL Pr							UNITS no_units			"Cooling liquid Prandtl number"
		REAL Nu							UNITS no_units			"Cooling liquid Nusselt number"
		REAL Re							UNITS no_units			"Cooling liquid Reynolds number"
		REAL v							UNITS u_m_s				"Cooling liquid speed"
		REAL f							UNITS no_units			"Darcy friction factor"
		REAL D_eq						UNITS u_m				"Circular equivalent diameter of a rectangular duct for equal friction and flow capacity"
		REAL D_hy						UNITS u_m				"Hydraulic diameter"
		
	INIT PRIORITY 100
		--Init_fluid(Vaporisation_mix(l.fluid), g.fluid)
	
		T_w_cold = 500.
		
	DISCRETE
		ASSERT (Liquid == State(l.fluid)) FATAL "ONLY LIQUIDS CAN ENTER TO THE COOLING JACKET!"

	CONTINUOUS
		--Init_fluid(Vaporisation_mix(l.fluid), g.fluid)
		Vaporisation_mix(l.fluid, g.fluid)
		g.W = l.W	
		h.Q = Q
		h.T = T_w_hot
		h.A = A_wet_nozzle
		
		A_wet_cooling = N * 2 * (a + b) * L
		
		Q = h_l * (T_w_cold - l.Tt) * A_wet_cooling
		Q = k_w / t_w * (T_w_hot - T_w_cold) * A_wet_nozzle
		Q = l.W * cp(g.fluid) * (g.Tt - l.Tt)
		
		h_l = Nu * cond(l.fluid) / D_hy
		Re = rho(l.fluid) * v * D_hy / visc(l.fluid)
		Pr = visc(l.fluid) * cp(l.fluid) / cond(l.fluid)
		-- Sieder-Tate correlation for turbulent flow
		Nu = 0.027 * Re**0.8 * Pr**0.33	
		
		v = l.W / (a * b * rho(l.fluid)) / N
		
		IF (Type == Known_dp) INSERT
			g.pt = l.pt - dp
		ELSE
			g.pt = l.pt - f * L / D_eq * 0.5 *  rho(l.fluid) * v**2
		END IF
		f = hdc_fric(D_eq, rug, Re)
		
		D_hy = 2. * a * b / (a + b)
		-- Huebscher (1948)
		D_eq = 1.3 * (a * b)**0.625 / (a + b)**0.25

END COMPONENT


--------------------------------------------------------------------------------
-- Component that represents a cooling jacket
--------------------------------------------------------------------------------
COMPONENT CoolingJacket_New (ENUM Type_CJ Type_HE = Liquid2Liquid, ENUM Type_Cooling Type_dp = Darcy)

  "Cooling jacket"

   PORTS
      IN Fluid f_in								"Inlet fluid port"
      OUT Fluid f_out							"Outlet fluid port"
		IN Heat h
		
	DATA
		REAL L = 1.				UNITS u_m 			"Length of the channels" 
		REAL a = 0.002			UNITS u_m 			"Average width of the channels" 
		REAL b = 0.004			UNITS u_m 			"Average height of the channels" 
		INTEGER N = 100		UNITS no_units		"Number of channels"
		REAL rug = 5.e-5		UNITS u_m			"Absolute rugosity"
		REAL k_w = 370			UNITS	u_W_mK		"Thermal conductivity of the nozzle wall material"
		REAL t_w = 0.003		UNITS u_m 			"Nozzle wall thickness" 
		REAL dp = 1.e5			UNITS u_Pa			"Imposed total pressure drop if Type=Known_dp"
		REAL visc_P = 0.001		UNITS u_Pas			"Change in combustion products viscosity due to change of state"
		REAL cond_P = 0.001		UNITS u_W_mK		"Change in combustion products thermal conductivity due to change of state"
		REAL rho_P =1000			UNITS u_kg_m3		"Change in fluid density if applied"
		
	DECLS
		REAL Q									UNITS u_W 				"Heat flux"
		REAL T_w_hot							UNITS u_K				"Wall temperature of the combustion gases side"				
		REAL T_w_cold							UNITS u_K				"Wall temperature of the cooling liquid side"	
		REAL A_wet_cooling			UNITS u_m2				"Cooling wet area"
		REAL A_wet_nozzle				UNITS u_m2				"Nozzle wet area"
		HIDDEN REAL fluid_out[ChemName]	UNITS no_units			"Outlet fluid"		
		REAL h_l									UNITS u_W_m2K			"Cooling liquid heat transfer coefficient"
		REAL T_av								UNITS u_K				"Average temperature"
		REAL P_av								UNITS u_Pa				"Average pressure"
		REAL Pr									UNITS no_units			"Cooling liquid Prandtl number"
		REAL Nu									UNITS no_units			"Cooling liquid Nusselt number"
		REAL Re									UNITS no_units			"Cooling liquid Reynolds number"
		REAL v									UNITS u_m_s				"Cooling liquid speed"
		REAL f									UNITS no_units			"Darcy friction factor"
		REAL D_eq								UNITS u_m				"Circular equivalent diameter of a rectangular duct for equal friction and flow capacity"
		REAL D_hy								UNITS u_m				"Hydraulic diameter"
		
		HIDDEN REAL Cp_v[setofSize(LV_extended)+2]			UNITS u_J_kgK			"Vector for vaporisation or condensation calculations. Cp of changing substances"
		HIDDEN REAL Lv_v[setofSize(LV_extended)+2]			UNITS u_J_kg			"Vector for vaporisation or condensation calculations. Heat of vaporisation of changing substances"
		HIDDEN REAL T_boil_v[setofSize(LV_extended)+2]			UNITS u_K			"Vector for vaporisation or condensation calculations. Boiling temperature of changing substances"
		HIDDEN REAL mass_frac_v[setofSize(LV_extended)+2]			UNITS u_kg_s			"Vector for vaporisation or condensation calculations. Mass fraction of changing substances in composition vector"
		
	INIT PRIORITY 100
		--Init_fluid(Vaporisation_mix(l.fluid), g.fluid)
	
		T_w_cold = 500.
		
	DISCRETE
		ASSERT (Liquid == State(f_in.fluid)) FATAL "ONLY LIQUIDS CAN ENTER TO THE COOLING JACKET!"

	CONTINUOUS
		T_av = (f_in.Tt + f_out.Tt)/2
		P_av = (f_in.pt + f_out.pt)/2
		Obtain_Vap_Prop(f_in.fluid, Cp_v, Lv_v, T_boil_v, mass_frac_v)
		f_in.W = f_out.W
	
		IF (Type_HE == Liquid2Liquid) INSERT
		
			fluid_out = f_out.fluid
			f_in.fluid = f_out.fluid
			h.Q = Q
			h.T = T_w_hot
			h.A = A_wet_nozzle
			
			Q = h_l * (T_w_cold - T_av) * A_wet_cooling
			Q = k_w / t_w * (T_w_hot - T_w_cold) * A_wet_nozzle
			Q = f_in.W * cp(f_in.fluid) * (f_out.Tt - f_in.Tt)
			
			h_l = Nu * cond(f_in.fluid) / D_hy
			
			Re = rho(f_in.fluid) * v * D_hy / visc(f_in.fluid)
			Pr = visc(f_in.fluid) * cp(f_in.fluid) / cond(f_in.fluid)
			-- Sieder-Tate correlation for turbulent flow
			Nu = 0.027 * Re**0.8 * Pr**0.33	
			
			v = f_in.W / (a * b * rho(f_in.fluid)) / N
			
			IF (Type_dp == Known_dp) INSERT
				f_out.pt = f_in.pt - dp
			ELSE
				f_out.pt = f_in.pt - f * L / D_eq * 0.5 *  rho(f_in.fluid) * v**2
			END IF
			f = hdc_fric(D_eq, rug, Re)
			
			D_hy = 2. * a * b / (a + b)
			-- Huebscher (1948)
			D_eq = 1.3 * (a * b)**0.625 / (a + b)**0.25
			
		ELSEIF (Type_HE == Liquid2Gas) INSERT
			
			EXPAND (i IN ChemName EXCEPT LV_extended) fluid_out[i] = f_in.fluid[i]
			EXPAND (i IN LV) fluid_out[i] = 0
			fluid_out[O2] = f_in.fluid[LOX]
			fluid_out[H2] = f_in.fluid[LH2]
			fluid_out[CH4] = f_in.fluid[LCH4]
			fluid_out[MMH_vapour] = f_in.fluid[MMH]
			
			
			EXPAND (i IN ChemName EXCEPT Comb_prod_visc, Comb_prod_cond, Comb_prod_rho) f_out.fluid[i] = fluid_out[i]
			f_out.fluid[Comb_prod_visc] = visc_P
			f_out.fluid[Comb_prod_cond] = cond_P
			f_out.fluid[Comb_prod_rho] = rho_P
			
			h.Q = Q
			h.T = T_w_hot
			h.A = A_wet_nozzle
			
			Q = f_in.W * SUM (i IN 1,(setofSize(LV_extended)+2)/2; mass_frac_v[2*i-1] * Cp_v[2*i-1] * (T_boil_v[2*i-1]-f_in.Tt) \
			+ mass_frac_v[2*i-1] * Cp_v[2*i] * (f_out.Tt-T_boil_v[2*i-1]) + mass_frac_v[2*i-1] * Lv_v[2*i-1]) 
			Q = h_l * (T_w_cold - T_av) * A_wet_cooling
			Q = k_w / t_w * (T_w_hot - T_w_cold) * A_wet_nozzle
			
			h_l = Nu * cond(f_in.fluid) / D_hy
			
			Re = ((f_out.pt/(R(f_out.fluid)*f_out.Tt)) * v * D_hy / visc(f_out.fluid))+(rho(f_in.fluid) * v * D_hy / visc(f_in.fluid))
			Pr = (visc(f_out.fluid) * cp(f_out.fluid) / cond(f_out.fluid)) + (visc(f_in.fluid) * cp(f_in.fluid) / cond(f_in.fluid))
			-- Sieder-Tate correlation for turbulent flow
			Nu = 0.027 * Re**0.8 * Pr**0.33	
			
			v = f_in.W / (a * b * rho(f_in.fluid)) / N
			
			f_out.pt = f_in.pt - dp
			
			f = 0
			
			D_hy = 2. * a * b / (a + b)
			-- Huebscher (1948)
			D_eq = 1.3 * (a * b)**0.625 / (a + b)**0.25
		END IF

END COMPONENT


--------------------------------------------------------------------------------------------
-- Component that represents a liquid heat exchanger without vaporisation. Liquid to liquid
--------------------------------------------------------------------------------------------
COMPONENT Heat_exchanger (ENUM Type_Heat_Ex Type_HE = Liquid2Liquid, ENUM Type_Cooling Type_dp = Darcy)

	PORTS
      IN Fluid f_in								"Inlet fluid port"
      OUT Fluid f_out							"Outlet fluid port"
		IN Heat h			
   		
	DATA
		REAL L = 1.					UNITS u_m 			"Length of the channels" 
		REAL a = 0.002				UNITS u_m 			"Average width of the channels" 
		REAL b = 0.004				UNITS u_m 			"Average height of the channels" 
		INTEGER N = 100			UNITS no_units		"Number of channels"
		REAL rug = 5.e-5			UNITS u_m			"Absolute rugosity"
		REAL dp = 1.e5				UNITS u_Pa			"Imposed total pressure drop if Type=Known_dp"
		REAL visc_P = 0.001		UNITS u_Pas			"Change in combustion products viscosity due to change of state"
		REAL cond_P = 0.001		UNITS u_W_mK		"Change in combustion products thermal conductivity due to change of state"
		REAL rho_P =1000			UNITS u_kg_m3		"Change in fluid density if applied"
		
		
	DECLS
		REAL Q									UNITS u_W 				"Heat flux"
		REAL T_w									UNITS u_K				"Wall temperature of the cooling liquid side"	
		REAL A_wet_cooling					UNITS u_m2				"Cooling wet area"
		HIDDEN REAL fluid_out[ChemName]	UNITS no_units			"Outlet fluid"		
		REAL h_l									UNITS u_W_m2K			"Cooling liquid heat transfer coefficient"
		REAL T_av								UNITS u_K				"Average temperature"
		REAL P_av								UNITS u_Pa				"Average pressure"
		REAL Pr									UNITS no_units			"Cooling liquid Prandtl number"
		REAL Nu									UNITS no_units			"Cooling liquid Nusselt number"
		REAL Re									UNITS no_units			"Cooling liquid Reynolds number"
		REAL v									UNITS u_m_s				"Cooling liquid speed"
		REAL f									UNITS no_units			"Darcy friction factor"
		REAL D_eq								UNITS u_m				"Circular equivalent diameter of a rectangular duct for equal friction and flow capacity"
		REAL D_hy								UNITS u_m				"Hydraulic diameter"
		
		HIDDEN REAL Cp_v[setofSize(LV_extended)+2]			UNITS u_J_kgK			"Vector for vaporisation or condensation calculations. Cp of changing substances"
		HIDDEN REAL Lv_v[setofSize(LV_extended)+2]			UNITS u_J_kg			"Vector for vaporisation or condensation calculations. Heat of vaporisation of changing substances"
		HIDDEN REAL T_boil_v[setofSize(LV_extended)+2]			UNITS u_K			"Vector for vaporisation or condensation calculations. Boiling temperature of changing substances"
		HIDDEN REAL mass_frac_v[setofSize(LV_extended)+2]			UNITS u_kg_s			"Vector for vaporisation or condensation calculations. Mass fraction of changing substances in composition vector"
		
	INIT PRIORITY 100
		
		ASSERT((Type_HE == Liquid2Liquid AND State(f_in.fluid) == Liquid) OR \
		((Type_HE == Gas2Gas AND State(f_in.fluid) == Gas))) WARNING "FLUID CHANGED ITS STATE!"
		
	DISCRETE
		ASSERT((Type_HE == Liquid2Liquid AND State(f_in.fluid) == Liquid) OR \
		(Type_HE == Gas2Gas AND State(f_in.fluid) == Gas) OR \
		(Type_HE == Liquid2Gas AND State(f_in.fluid) == Liquid) OR \
		(Type_HE == Gas2Liquid AND State(f_in.fluid) == Gas)) KILLPOINT "CHECK HEAT EXCHANGER TYPE! SELECTED TYPE DOES NOT MATCH ACTUAL FLUID STATE THROUGH HEAT EXCHANGER"	

	CONTINUOUS
		T_av = (f_in.Tt + f_out.Tt)/2
		P_av = (f_in.pt + f_out.pt)/2
		Obtain_Vap_Prop(f_in.fluid, Cp_v, Lv_v, T_boil_v, mass_frac_v)
		f_in.W = f_out.W

		
		IF (Type_HE == Liquid2Liquid) INSERT
		
			fluid_out = f_out.fluid
			f_in.fluid = f_out.fluid
			h.Q = Q
			h.T = T_w
			h.A = A_wet_cooling
			
			Q = f_in.W * cp(f_in.fluid) * (f_out.Tt - f_in.Tt)
			Q = h_l * (T_w - T_av) * A_wet_cooling
			
			h_l = Nu * cond(f_in.fluid) / D_hy
			
			Re = rho(f_in.fluid) * v * D_hy / visc(f_in.fluid)
			Pr = visc(f_in.fluid) * cp(f_in.fluid) / cond(f_in.fluid)
			-- Sieder-Tate correlation for turbulent flow
			Nu = 0.027 * Re**0.8 * Pr**0.33	
			
			v = f_in.W / (a * b * rho(f_in.fluid)) / N
			
			IF (Type_dp == Known_dp) INSERT
				f_out.pt = f_in.pt - dp
			ELSE
				f_out.pt = f_in.pt - f * L / D_eq * 0.5 *  rho(f_in.fluid) * v**2
			END IF
			f = hdc_fric(D_eq, rug, Re)
			
			D_hy = 2. * a * b / (a + b)
			-- Huebscher (1948)
			D_eq = 1.3 * (a * b)**0.625 / (a + b)**0.25
			
		ELSEIF (Type_HE == Gas2Gas) INSERT
		
			fluid_out = f_out.fluid
			f_in.fluid = f_out.fluid
			h.Q = Q
			h.T = T_w
			h.A = A_wet_cooling
			
			Q = f_in.W * cp(f_in.fluid) * (f_out.Tt - f_in.Tt)
			Q = h_l * (T_w - T_av) * A_wet_cooling
			
			h_l = Nu * cond(f_in.fluid) / D_hy
			
			Re = (P_av/(R(f_in.fluid)*T_av)) * v * D_hy / visc(f_in.fluid)
			Pr = visc(f_in.fluid) * cp(f_in.fluid) / cond(f_in.fluid)
			-- Sieder-Tate correlation for turbulent flow
			Nu = 0.027 * Re**0.8 * Pr**0.33	
			
			v = f_in.W / (a * b * rho(f_in.fluid)) / N
			
			f_out.pt = f_in.pt - dp
			
			f = hdc_fric(D_eq, rug, Re)
			
			D_hy = 2. * a * b / (a + b)
			-- Huebscher (1948)
			D_eq = 1.3 * (a * b)**0.625 / (a + b)**0.25
			
		ELSEIF (Type_HE == Gas2Liquid) INSERT
			
			EXPAND (i IN ChemName EXCEPT LV_extended) fluid_out[i] = f_in.fluid[i]
			EXPAND (i IN LV_gases) fluid_out[i] = 0
			fluid_out[LOX] = f_in.fluid[O2]
			fluid_out[LH2] = f_in.fluid[H2]
			fluid_out[LCH4] = f_in.fluid[CH4]
			fluid_out[MMH] = f_in.fluid[MMH_vapour]
			
			EXPAND (i IN ChemName EXCEPT Comb_prod_visc, Comb_prod_cond, Comb_prod_rho) f_out.fluid[i] = fluid_out[i]
			f_out.fluid[Comb_prod_visc] = visc_P
			f_out.fluid[Comb_prod_cond] = cond_P
			f_out.fluid[Comb_prod_rho] = rho_P
			h.Q = Q
			h.T = T_w
			h.A = A_wet_cooling
			
			
			Q = f_in.W * SUM (i IN 1,(setofSize(LV_extended)+2)/2; mass_frac_v[2*i] * Cp_v[2*i-1] * (f_out.Tt-T_boil_v[2*i]) \
			+ mass_frac_v[2*i] * Cp_v[2*i] * (T_boil_v[2*i]-f_in.Tt) - mass_frac_v[2*i] * Lv_v[2*i]) 
			Q = h_l * (T_w - T_av) * A_wet_cooling
			
			h_l = Nu * cond(f_in.fluid) / D_hy
			
			Re = ((f_in.pt/(R(f_in.fluid)*f_in.Tt)) * v * D_hy / visc(f_in.fluid))+(rho(f_out.fluid) * v * D_hy / visc(f_out.fluid))
			Pr = (visc(f_in.fluid) * cp(f_in.fluid) / cond(f_in.fluid)) + (visc(f_out.fluid) * cp(f_out.fluid) / cond(f_out.fluid))
			-- Sieder-Tate correlation for turbulent flow
			Nu = 0.027 * Re**0.8 * Pr**0.33	
			
			v = f_in.W / (a * b * rho(f_in.fluid)) / N
			
			f_out.pt = f_in.pt - dp
			
			f = 0
			
			D_hy = 2. * a * b / (a + b)
			-- Huebscher (1948)
			D_eq = 1.3 * (a * b)**0.625 / (a + b)**0.25
		
		ELSEIF (Type_HE == Liquid2Gas) INSERT
			
			EXPAND (i IN ChemName EXCEPT LV_extended) fluid_out[i] = f_in.fluid[i]
			EXPAND (i IN LV) fluid_out[i] = 0
			fluid_out[O2] = f_in.fluid[LOX]
			fluid_out[H2] = f_in.fluid[LH2]
			fluid_out[CH4] = f_in.fluid[LCH4]
			fluid_out[MMH_vapour] = f_in.fluid[MMH]
			
			
			EXPAND (i IN ChemName EXCEPT Comb_prod_visc, Comb_prod_cond, Comb_prod_rho) f_out.fluid[i] = fluid_out[i]
			f_out.fluid[Comb_prod_visc] = visc_P
			f_out.fluid[Comb_prod_cond] = cond_P
			f_out.fluid[Comb_prod_rho] = rho_P
			
			h.Q = Q
			h.T = T_w
			h.A = A_wet_cooling
			
			Q = f_in.W * SUM (i IN 1,(setofSize(LV_extended)+2)/2; mass_frac_v[2*i-1] * Cp_v[2*i-1] * (T_boil_v[2*i-1]-f_in.Tt) \
			+ mass_frac_v[2*i-1] * Cp_v[2*i] * (f_out.Tt-T_boil_v[2*i-1]) + mass_frac_v[2*i-1] * Lv_v[2*i-1]) 
			Q = h_l * (T_w - T_av) * A_wet_cooling
			
			h_l = Nu * cond(f_in.fluid) / D_hy
			
			Re = ((f_out.pt/(R(f_out.fluid)*f_out.Tt)) * v * D_hy / visc(f_out.fluid))+(rho(f_in.fluid) * v * D_hy / visc(f_in.fluid))
			Pr = (visc(f_out.fluid) * cp(f_out.fluid) / cond(f_out.fluid)) + (visc(f_in.fluid) * cp(f_in.fluid) / cond(f_in.fluid))
			-- Sieder-Tate correlation for turbulent flow
			Nu = 0.027 * Re**0.8 * Pr**0.33	
			
			v = f_in.W / (a * b * rho(f_in.fluid)) / N
			
			f_out.pt = f_in.pt - dp
			
			f = 0
			
			D_hy = 2. * a * b / (a + b)
			-- Huebscher (1948)
			D_eq = 1.3 * (a * b)**0.625 / (a + b)**0.25
		END IF
END COMPONENT

COMPONENT Contact_Wall 

	PORTS
		OUT Heat h_in
		OUT Heat h_out
		
	DATA
		REAL k_w = 0.1 							UNITS u_W_mK 	"Wall material thermal conductivity"
		REAL t_w = 0.1								UNITS u_m		"Wall thickness"
		REAL A_eff_1 = 1							UNITS no_units	"Ratio between effective area 1 and conductive area"
		REAL A_eff_2 = 1							UNITS no_units	"Ratio between effective area 2 and conductive area"
	DECLS
		REAL A_wet
	   REAL Q
	CONTINUOUS
	
		h_in.Q = -h_out.Q
		h_in.A =  A_eff_1 * A_wet
		h_out.A =  A_eff_2 * A_wet
		Q = k_w / t_w * (h_out.T - h_in.T) * A_wet
		h_in.Q = Q
		
END COMPONENT

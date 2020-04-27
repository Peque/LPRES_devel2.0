--------------------------------------------------------------------------------
-- EcosimPro Simulation Source Code
-- FILE NAME: LPRES_CompGasGen.el
-- DESCRIPTION: Defines combustion chamber type components with their inyectors
-- NOTES:
-- AUTHOR: Pablo Sierra Heras
-- CREATION DATE: 05/12/2014
--------------------------------------------------------------------------------

-- Libraries
USE MATH VERSION "3.1"


-- Components
--------------------------------------------------------------------------------
-- Component that represents an injector
--------------------------------------------------------------------------------


COMPONENT Injector (ENUM ChemState Type = Gas)

   "Injector"
	
	PORTS
      IN Fluid f_in						"Inlet fluid port"
      OUT FluidInj f_out				"Outlet fluid injection port"
	
	DATA
		REAL C_D = 0.5			UNITS no_units    	"Discharge coefficient (used only in liquid state)"
		REAL A = 0.05		 	UNITS u_m2				"Output area"
	
	DECLS
		REAL A_ch				UNITS u_m2			"Choked area (calculated only for gases)"
		
		REAL PR = 10			UNITS no_units		"Pressure ratio"
		REAL PR_ch				UNITS no_units		"Choked pressure ratio (calculated only for gases)"
			
		REAL M_out				UNITS no_units		"Outlet Mach number (calculated only for gases)"		
		
		REAL p_out_ch			UNITS u_Pa			"Choked outlet pressure (calculated only for gases)"		
		
		REAL v_ideal			UNITS u_m_s			"Ideal outlet speed"
		REAL Re					UNITS no_units		"Outlet Reynolds number (calculated only for liquids)"
		
	INIT
		
		
	DISCRETE
		ASSERT(A > 0) WARNING "NEGATIVE INJECTOR AREA"
		ASSERT((Type == Gas AND State(f_in.fluid) == Gas) OR ((Type == Liquid AND State(f_in.fluid) == Liquid))) WARNING "CHECK INJECTOR TYPE! SELECTED TYPE DOES NOT MATCH ACTUAL FLUID STATE THROUGH INJECTOR"	
		
		
	CONTINUOUS
		f_in.fluid = f_out.fluid
		f_in.W = f_out.W 		
		
		PR = f_in.pt / f_out.p_c
						
						
			IF (Type == Liquid) INSERT 
			
			PR_ch = 0
			p_out_ch = 0
			A_ch = 0
			M_out = 0
			f_out.p = f_out.p_c
			f_out.T = f_in.Tt
			f_in.W = A * C_D * sqrt(2. * rho(f_in.fluid) * (f_in.pt - f_out.p))	
			v_ideal = f_in.W / (rho(f_in.fluid) * A * C_D)
			Re = rho(f_in.fluid) * v_ideal * C_D * sqrt(4. * A / PI) / visc(f_in.fluid) 
			
			END IF
			
			
			IF (Type == Gas) INSERT 
			
			PR_ch = ((gamma(f_in.fluid) + 1.) / 2.)**(gamma(f_in.fluid) / (gamma(f_in.fluid) - 1.))
			p_out_ch = f_in.pt / PR_ch	
			A_ch = sqrt(f_in.Tt * R(f_in.fluid)) * f_in.W / (FGAMMA(f_in.fluid) * f_in.pt)
			M_out = min(sqrt(2. * (PR**((gamma(f_in.fluid) - 1.) / gamma(f_in.fluid)) - 1.) / (gamma(f_in.fluid) - 1.)),1)
			f_out.p = max(f_out.p_c, p_out_ch)
			f_out.T = f_in.Tt / (1. + (gamma(f_in.fluid) - 1.) / 2. * M_out**2)
			f_in.W = A * FGAMMA(f_in.fluid) * f_in.pt / sqrt(f_in.Tt * R(f_in.fluid)) / (((2. + (gamma(f_in.fluid) - 1.) * M_out**2) / (gamma(f_in.fluid) + 1.))**((gamma(f_in.fluid) + 1.) / (2. * (gamma(f_in.fluid) - 1.))) / M_out)
			v_ideal = M_out * sqrt(gamma(f_in.fluid) * R(f_in.fluid) * f_out.T)
			Re = 0 
			
			END IF
						
			
		
END COMPONENT


--------------------------------------------------------------------------------
-- Component that represents a combustion chamber
--------------------------------------------------------------------------------
COMPONENT CombCha (ENUM YesNo Cooled = No)

	"Combustion chamber"

	PORTS
      IN FluidInj f_O                			"Inlet fluid injection port"
      IN FluidInj f_F              				"Inlet fluid injection port"
		OUT GasNozzle g								"Outlet gas port through a nozzle"
		OUT Info (n = 1) i							"Outlet information port"
		OUT Heat h										"Outlet heat exchange port"
		
	DATA
		REAL eta_d = 0.9						UNITS no_units		"Design combustion efficiency"
		REAL OF_st = 8.						UNITS no_units		"Stoichiometric mixture ratio"
		REAL Q_comb	= 2000000.				UNITS u_J_kg		"Heat of combustion per oxidant mass flow unit"
		REAL cp_P = 4182.						UNITS u_J_kgK		"Specific heat at constant pressure of the products using a stoichiometric mixture. Gas state"
		REAL M_P	= 32.							UNITS u_g_mol		"Molar mass of the products using a stoichiometric mixture"
		REAL c_P	= 4182						UNITS u_J_kgK		"Specific heat of the products using a stoichiometric mixture. Liquid state"	
		REAL Lv_P = 800.						UNITS u_J_kg		"Heat of vaporisation of the products using a stoichiometric mixture"	
		REAL visc_P = 0.001					UNITS u_Pas			"Viscosity of combusion products"
		REAL cond_P = 0.001					UNITS u_W_mK		"Thermal conductivity of combustion products"
		REAL T_boil_P = 200					UNITS u_K			"Boiling temperature of combustion products"
		REAL AR = 10.							UNITS no_units		"Area ratio"	
		REAL A_th = 0.05						UNITS u_m2			"Throat area"
		REAL AR_r = 10. / 2.					UNITS no_units		"Area at the characteristic section of heat exchange divided by the throat area"			
		REAL A_wet = 1.						UNITS u_m2			"Nozzle wet area of the cooled zone"
		ENUM ConDiv Zone = Divergent								"Convergent if the characteristic section of heat exchange is placed in the convergent zone of the nozzle, Divergent if it is placed in the divergent zone"
				
	DECLS
		REAL T_in								UNITS u_K			"Average inlet temperature"
		REAL T_c									UNITS u_K			"Combustion temperature"
	
		REAL eta									UNITS no_units		"Combustion efficiency"
		REAL OF									UNITS no_units		"Mixture ratio"
		REAL phi	= 1							UNITS no_units		"Equivalence ratio"
		REAL W_F_st								UNITS u_kg_s		"Fuel mass flow that provides a stoichiometric mixture when combined with the oxidant mass flow"
		ALG REAL W_F							UNITS u_kg_s		"Fuel mass flow"
		REAL W_O									UNITS u_kg_s		"Oxidant mass flow"
		REAL W_IF								UNITS u_kg_s		"Inert mass flow through the fuel port"
		REAL W_IO								UNITS u_kg_s		"Inert mass flow through the oxidant port"
		REAL p_c									UNITS u_Pa			"Combustion pressure"
		
		HIDDEN REAL fluid_O[ChemName]		UNITS no_units		"Oxidant fluid"		
		HIDDEN REAL fluid_F[ChemName]		UNITS no_units		"Fuel fluid"
		HIDDEN REAL fluid_P[ChemName]		UNITS no_units		"Combustion products fluid"
		HIDDEN REAL fluid_P_vap[ChemName]		UNITS no_units		"Vaporised combustion products fluid"
		HIDDEN REAL fluid_P_decomp [Decomposition]	UNITS no_units "Properties of the gases in vaporised combustion products fluid"
		
		REAL cp_R								UNITS u_J_kgK		"Reactant specific heat at constant pressure"
		
		BOOLEAN Combustion											"TRUE if there is combustion and FALSE if either oxidant or fuel is lacking"	
		REAL Q_comb_effective				UNITS u_J_kg		"Effective heat of combustion per oxidant mass flow unit"
				
		REAL A_out								UNITS u_m2			"Output area"
				
		ALG REAL p_out_ch = 100.			UNITS u_Pa			"Choked outlet pressure"	
		
		REAL c_star								UNITS u_m_s			"Characteristic exhaust velocity"	
				
		REAL A_r									UNITS u_m2			"Area at the characteristic section of heat exchange"
		REAL M_r									UNITS no_units		"Mach number at the characteristic section of heat exchange"
		REAL Pr_r								UNITS no_units		"Prandtl number at the characteristic section of heat exchange"	
		REAL visc_r								UNITS u_Pas			"Dynamic viscosity at the characteristic section of heat exchange"
		REAL T_aw								UNITS u_K			"Adiabatic wall temperature of combustion gases"
		REAL h_g									UNITS u_W_m2K		"Combustion gases heat transfer coefficient"

	INIT PRIORITY 100
		ASSERT (AR_r >= 1) KILLPOINT "AR_r CAN NOT BE LOWER THAN 1!"
		
		Init_fluid(Comb_prod, g.fluid)		
		
		Combustion = TRUE
		
		IF (Zone == Convergent) THEN
			M_r = 0.001
		ELSE
			M_r = 2.8012363276683829
		END IF
		
	DISCRETE
		WHEN (SUM (i IN LiquidsGases ; fluid_O[i] * fluid_F[i]) != 0) THEN
			Combustion = FALSE
		END WHEN
		
		WHEN (SUM (i IN LiquidsGases ; fluid_O[i] * fluid_F[i]) == 0) THEN
			Combustion = TRUE
		END WHEN
		
	CONTINUOUS	
		g.A_out = A_out
		AR = A_out / A_th		
		
		eta = eta_d		
		
		g.pt = f_O.p_c
		g.pt = f_F.p_c
		g.pt = p_c
		
		f_O.W + f_F.W = g.W
		W_O = f_O.W * (1 - f_O.fluid[Comb_prod])
		W_F = f_F.W * (1 - f_F.fluid[Comb_prod])
		W_IO = f_O.W - W_O
		W_IF = f_F.W - W_F
		OF = W_O / W_F
		phi = OF_st / OF
		W_F_st = W_F / phi
		
		Q_comb_effective = 	ZONE (Combustion)		Q_comb
									OTHERS					0
				
		eta * Q_comb_effective = ((1. + OF) / min(OF,OF_st)) * (cp(fluid_P_vap) * (T_c - T_ref) - cp_R * (T_in - T_ref))
		(1 + phi / OF_st) * cp_R * T_in = cp(fluid_O) * f_O.T + phi / OF_st * cp(fluid_F) * f_F.T		
		
		(W_O + W_F)*fluid_P_decomp[CombProd]/ g.W * fluid_P_decomp[CombProd_cp] * (g.Tt - T_c)                         \
		+ ((W_IO / g.W) * f_O.fluid[Comb_prod_cp]) * (g.Tt - f_O.T)       \
		+ ((W_IF / g.W) * f_F.fluid[Comb_prod_cp]) * (g.Tt - f_F.T) = 0
		
		EXPAND (i IN LiquidsGases) fluid_O[i] = f_O.fluid[i] / (1 - f_O.fluid[Comb_prod])
		EXPAND (i IN Comb_prod_prop) fluid_O[i] = 0
		
		
		EXPAND (i IN LiquidsGases) fluid_F[i] = f_F.fluid[i] / (1 - f_F.fluid[Comb_prod])
		EXPAND (i IN Comb_prod_prop) fluid_F[i] = 0
				
		cp_R = (W_O * cp(fluid_O) + W_F * cp(fluid_F)) / (W_O + W_F)
				
		EXPAND (i IN LiquidsGases) fluid_P[i] = 	IF (Combustion)	(fluid_O[i] * max(1 - phi,0) * W_O + fluid_F[i] * max(phi - 1,0) * W_F_st) / (W_O + W_F)
																ELSE					(fluid_O[i] * W_O + fluid_F[i] * W_F) / (W_O + W_F)
		fluid_P[Comb_prod] = 	IF (Combustion)	(fluid_O[Comb_prod] * max(1 - phi,0) * W_O + fluid_F[Comb_prod] * max(phi - 1,0) * W_F_st + (1 - max(1 - phi,0)) * (W_O + W_F_st)) / (W_O + W_F)
										ELSE					0
		fluid_P[Comb_prod_M] = 	IF (Combustion)	M_P
										ELSE					0
		fluid_P[Comb_prod_cp] = IF (Combustion)	cp_P
										ELSE					0
		fluid_P[Comb_prod_c] = 	IF (Combustion)	c_P
										ELSE					0
		fluid_P[Comb_prod_cp_g] = 	IF (Combustion)	cp_P
											ELSE	            0
		fluid_P[Comb_prod_Lv] = IF (Combustion)   Lv_P
										ELSE 					0
		fluid_P[Comb_prod_visc] = IF (Combustion)   visc_P
										ELSE 					0
		fluid_P[Comb_prod_cond] = IF (Combustion)   cond_P
										ELSE 					0
		fluid_P[Comb_prod_T_boil] = IF (Combustion)   T_boil_P
										ELSE 					0
		fluid_P[Comb_prod_rho] = 0				
		
		Vaporisation_mix(fluid_P, fluid_P_vap)
		Decomposition_fun(fluid_P_vap, fluid_P_decomp)
		
		EXPAND (i IN LiquidsGases) g.fluid[i] = fluid_P_vap[i] * (W_O + W_F) / g.W
		g.fluid[Comb_prod] = (fluid_P_vap[Comb_prod] * (W_O + W_F) + (W_IO + W_IF)) / g.W		
		EXPAND (i IN Comb_prod_prop EXCEPT Comb_prod) g.fluid[i] = fluid_P_vap[i]


		
		g.W = g.pt * A_th * FGAMMA(g.fluid) / sqrt(R(g.fluid) * g.Tt) --R(g.fluid)
		
		A_out / A_th = FGAMMA(g.fluid) / ((p_out_ch / g.pt)**(1. / gamma(g.fluid)) * sqrt(2 * gamma(g.fluid) * (1 - (p_out_ch / g.pt)**((gamma(g.fluid) - 1.) / gamma(g.fluid))) / (gamma(g.fluid) - 1)))
		
		c_star = sqrt(R(g.fluid) * g.Tt) / FGAMMA(g.fluid)
		i.Data[1] = c_star	
		g.c_star = c_star
		g.A_th = A_th
		
		
		-- Cooling of the nozzle
		
		T_aw = g.Tt * ((1. + Pr_r**0.33 * (gamma(g.fluid) - 1.) * M_r**2 / 2.) / (1. + (gamma(g.fluid) - 1.) * M_r**2 / 2.))
		Pr_r = 4. * gamma(g.fluid) / (9. * gamma(g.fluid) - 5.)
		visc_r = 1.184e-7 * M(g.fluid)**0.5 * T_aw**0.6		
		
		IF (Cooled == No) INSERT
			h.A = 0
			A_r = A_th
			M_r = 1.
			
			h_g = 0
			h.T = T_aw
			h.Q = 0
			
		ELSE		
			h.A = A_wet
			A_r = A_th * AR_r
			AR_r = 1 / M_r * FGAMMA(g.fluid) / sqrt(gamma(g.fluid)) * (1. + (gamma(g.fluid) - 1.) * M_r**2 / 2.)**((gamma(g.fluid) + 1.) / (2. * (gamma(g.fluid) - 1.)))
		
			-- Bartz correlation
			h_g = 0.026 / sqrt(4 * A_th / PI)**0.2 * (visc_r**0.2 * cp(g.fluid) / Pr_r**0.6) * (g.pt / c_star)**0.8 * (A_th / A_r)**0.9
			h.Q = h_g * (T_aw - h.T) * h.A
			
		END IF		
		
END COMPONENT

COMPONENT MCombCha (ENUM Reaction_mono Reaction = Any, ENUM YesNo Cooled = No)

	"Combustion chamber"

	PORTS
      IN FluidInj f	                					"Inlet fluid injection port"
		OUT GasNozzle g										"Outlet gas port through a nozzle"
		OUT Info (n = 1) i									"Outlet information port"
		OUT Heat h												"Outlet heat exchange port"
		
	DATA
		REAL eta_d = 0.9								UNITS no_units		"Design combustion efficiency"
		
		REAL AR = 10.									UNITS no_units		"Area ratio"	
		REAL A_th = 0.05								UNITS u_m2			"Throat area"
		REAL Cp_data = 1								UNITS u_J_kgK		"Average Specific heat of the products"
		REAL M_data = 14.1							UNITS u_g_mol		"Average molar mass of the products"
		REAL c_P_data	= 4182						UNITS u_J_kgK		"Specific heat of the products using a stoichiometric mixture. Liquid state"	
		REAL Lv_P_data = 800.						UNITS u_J_kg		"Molar mass of the products using a stoichiometric mixture"	
		
		REAL Q_data = 2725840.657					UNITS u_J_kg		"Heat per mass flow unit"
		ENUM ConDiv Zone = Divergent								"Convergent if the characteristic section of heat exchange is placed in the convergent zone of the nozzle, Divergent if it is placed in the divergent zone"
		REAL AR_r = 10. / 2.							UNITS no_units		"Area at the characteristic section of heat exchange divided by the throat area"			
		REAL A_wet = 1.								UNITS u_m2			"Nozzle wet area of the cooled zone"
		
				
	DECLS
		REAL T_in										UNITS u_K			"Average inlet temperature"
		REAL T_c	= 4000								UNITS u_K			"Combustion temperature"
		REAL W = 100									UNITS u_kg_s		"Mass flow through nozzle"
		DISCR REAL Mono_M[Mono_LiquidsGases] = {32.04516, 17.031, 28.01348, 2.01594, 18, 31.9988}		\
															UNITS u_g_mol		"Molar mass of each monopropellant chemical"
	 	HIDDEN REAL fluid_P[ChemName]						UNITS no_units		"Combustion products fluid"
		HIDDEN DISCR REAL fluid_N2H4[ChemName]			UNITS no_units		"Combustion products in hidrazine reaction"
		HIDDEN DISCR REAL fluid_H2O2[ChemName]			UNITS no_units		"Combustion products in hydrogen peroxide reaction"
		REAL M_P											UNITS u_g_mol		"Molar mass of the products"				
		REAL Cp_P										UNITS u_J_kgK		"Reactant specific heat at constant pressure"
		REAL c_P											UNITS u_J_kgK		"Specific heat of the products using a stoichiometric mixture. Liquid state"	
		REAL Lv_P 										UNITS u_J_kg		"Molar mass of the products using a stoichiometric mixture"	
		REAL visc_P 									UNITS u_Pas			"Viscosity of combusion products"
		REAL cond_P 									UNITS u_W_mK		"Thermal conductivity of combustion products"
		REAL T_boil_P									UNITS u_K			"Boiling temperature of combustion products"
		REAL p_c	= 5000000.							UNITS u_Pa			"Design combustion pressure"
		REAL Q_comb_effective						UNITS u_J_kg		"Effective combustion heat"
		REAL eta											UNITS no_units		"Combustion efficiency"
		REAL A_out										UNITS u_m2			"Output area"
				
		ALG REAL p_out_ch = 100.					UNITS u_Pa			"Choked outlet pressure"	
		
		REAL c_star										UNITS u_m_s			"Characteristic exhaust velocity"	
		
		REAL A_r											UNITS u_m2			"Area at the characteristic section of heat exchange"
		REAL M_r											UNITS no_units		"Mach number at the characteristic section of heat exchange"
		REAL Pr_r										UNITS no_units		"Prandtl number at the characteristic section of heat exchange"	
		REAL visc_r										UNITS u_Pas			"Dynamic viscosity at the characteristic section of heat exchange"
		REAL T_aw										UNITS u_K			"Adiabatic wall temperature of combustion gases"
		REAL h_g											UNITS u_W_m2K		"Combustion gases heat transfer coefficient"
		
	INIT PRIORITY 100
	
	ASSERT (AR_r >= 1) KILLPOINT "AR_r CAN NOT BE LOWER THAN 1!"
		
		--Init_fluid(Comb_prod, g.fluid)		
		
	
		
		IF (Zone == Convergent) THEN
			M_r = 0.001
		ELSE
			M_r = 100
		END IF
		
		
		
	DISCRETE
	
		EXPAND (i IN ChemName) 
		WHEN (i != N2 OR i != NH3 OR i != H2 ) THEN
			fluid_N2H4[i] = 0
		END WHEN
		EXPAND (i IN ChemName)  
		WHEN (i == NH3) THEN
			fluid_N2H4[i] =(Mono_M[NH3]*4/5)/(Mono_M[NH3]*4/5 + Mono_M[N2]*3/5 + Mono_M[H2]*4/5)
		END WHEN
		EXPAND (i IN ChemName) 
		WHEN (i == N2) THEN
			fluid_N2H4[i] =(Mono_M[N2]*3/5)/(Mono_M[NH3]*4/5 + Mono_M[N2]*3/5 + Mono_M[H2]*4/5)
		END WHEN
		EXPAND (i IN ChemName) 
		WHEN (i == H2) THEN
			fluid_N2H4[i] =(Mono_M[H2]*4/5)/(Mono_M[NH3]*4/5 + Mono_M[N2]*3/5 + Mono_M[H2]*4/5)
		END WHEN
		
		EXPAND (i IN ChemName) 
		WHEN (i != O2 OR i != H2O_vapour) THEN
			fluid_H2O2[i] = 0
		END WHEN
		EXPAND (i IN ChemName)  
		WHEN (i == H2O_vapour) THEN
			fluid_H2O2[i] =(Mono_M[H2O_vapour])/(Mono_M[H2O_vapour]+ Mono_M[O2]*1/2)
		END WHEN
		EXPAND (i IN ChemName) 
		WHEN (i == O2) THEN
			fluid_H2O2[i] =(Mono_M[O2]*1/2)/(Mono_M[H2O_vapour]+ Mono_M[O2]*1/2)
		END WHEN
		
	CONTINUOUS	
		
		g.A_out = A_out
		AR = A_out / A_th	
		g.W = f.W
		T_in = f.T
		eta = eta_d		
		
		EXPAND (i IN LiquidsGases) fluid_P[i] = 0
		fluid_P[Comb_prod] = 1
		fluid_P[Comb_prod_M] = M_P
		fluid_P[Comb_prod_cp] = Cp_P	
		fluid_P[Comb_prod_cp_g] = Cp_P	
		fluid_P[Comb_prod_c] = c_P
		fluid_P[Comb_prod_Lv] = Lv_P
		fluid_P[Comb_prod_visc] = visc_P
		fluid_P[Comb_prod_cond] = cond_P
		fluid_P[Comb_prod_T_boil] = T_boil_P
		fluid_P[Comb_prod_rho] = 0
		fluid_P = g.fluid
		
		IF (Reaction == Hidrazine) INSERT
			M_P = M(fluid_N2H4)
			Cp_P = cp(fluid_N2H4)
			c_P = 0
			Lv_P = 0
			visc_P = 0
			cond_P = 0
			T_boil_P = 0
			Q_comb_effective = 2725840.657
		END IF
		
		IF (Reaction == H_Peroxide) INSERT
			M_P = M(fluid_H2O2)
			Cp_P = cp(fluid_H2O2)
			c_P = 0
			Lv_P = 0
			visc_P = 0
			cond_P= 0
			T_boil_P = 0
			Q_comb_effective = 1581669.102
		END IF
		
		IF (Reaction == Any) INSERT
			M_P = M_data
			Cp_P = Cp_data
			c_P = c_P_data
			Lv_P = Lv_P_data
			visc_P = 0
			cond_P = 0
			T_boil_P = 0
			Q_comb_effective = Q_data
		END IF 
			
		g.pt = p_c
		
		eta * Q_comb_effective =  Cp_P * (T_c - T_ref) - cp(f.fluid) * (T_in - T_ref)
		T_c = g.Tt
		g.pt = f.p_c
		
		g.W = g.pt * A_th * FGAMMA(g.fluid) / sqrt(R(g.fluid) * g.Tt)
		
		A_out / A_th = FGAMMA(g.fluid) / ((p_out_ch / g.pt)**(1. / gamma(g.fluid)) * sqrt(2 * gamma(g.fluid) * (1 - (p_out_ch / g.pt)**((gamma(g.fluid) - 1.) / gamma(g.fluid))) / (gamma(g.fluid) - 1)))
		
		c_star = sqrt(R(g.fluid) * g.Tt) / FGAMMA(g.fluid)
		i.Data[1] = c_star
		g.c_star = c_star
		g.A_th = A_th
		
		-- Cooling of the nozzle
		
		T_aw = g.Tt * ((1. + Pr_r**0.33 * (gamma(g.fluid) - 1.) * M_r**2 / 2.) / (1. + (gamma(g.fluid) - 1.) * M_r**2 / 2.))
		Pr_r = 4. * gamma(g.fluid) / (9. * gamma(g.fluid) - 5.)
		visc_r = 1.184e-7 * M(g.fluid)**0.5 * T_aw**0.6		
		
		IF (Cooled == No) INSERT
			h.A = 0
			A_r = A_th
			M_r = 1.
			
			h_g = 0
			h.T = T_aw
			h.Q = 0
			
		ELSE		
			h.A = A_wet
			A_r = A_th * AR_r
			AR_r = 1 / M_r * FGAMMA(g.fluid) / sqrt(gamma(g.fluid)) * (1. + (gamma(g.fluid) - 1.) * M_r**2 / 2.)**((gamma(g.fluid) + 1.) / (2. * (gamma(g.fluid) - 1.)))
		
			-- Bartz correlation
			h_g = 0.026 / sqrt(4 * A_th / PI)**0.2 * (visc_r**0.2 * cp(g.fluid) / Pr_r**0.6) * (g.pt / c_star)**0.8 * (A_th / A_r)**0.9
			h.Q = h_g * (T_aw - h.T) * h.A
			
		END IF		
		
		
		
END COMPONENT
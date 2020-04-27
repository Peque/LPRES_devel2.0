--------------------------------------------------------------------------------
-- EcosimPro Simulation Source Code
-- FILE NAME: LPRES_CompGasGen.el
-- DESCRIPTION: Defines gas generator type components
-- NOTES:
-- AUTHOR: Pablo Sierra Heras
-- CREATION DATE: 05/12/2014
--------------------------------------------------------------------------------

-- Libraries
USE MATH VERSION "3.1"


-- Components
--------------------------------------------------------------------------------
-- Component that represents a gas generator
--------------------------------------------------------------------------------


COMPONENT GasGen 

	"Gas generator"

	PORTS
      IN FluidInj f_O                				"Inlet fluid port"
      IN FluidInj f_F           	   				"Inlet fluid port"
      OUT Fluid g          						"Outlet fluid port"
		OUT Info (n = 1) i							"Outlet information port"
		
	DATA
		REAL TPL	= 1							UNITS no_units		"Total pressure loss"
		REAL eta_d = 0.9						UNITS no_units		"Design combustion efficiency"
		REAL OF_st = 8.						UNITS no_units		"Stoichiometric mixture ratio"
		REAL Q_comb	= 2000000.				UNITS u_J_kg		"Heat of combustion per oxidant mass flow unit"
		REAL cp_P = 4182.						UNITS u_J_kgK		"Specific heat at constant pressure of the products using a stoichiometric mixture. Gas state"
		REAL M_P	= 32.							UNITS u_g_mol		"Molar mass of the products using a stoichiometric mixture"
		REAL c_P	= 4182.						UNITS u_J_kgK		"Specific heat of the products using a stoichiometric mixture. Liquid state"	
		REAL Lv_P = 800.						UNITS u_J_kg		"Molar mass of the products using a stoichiometric mixture"	
		REAL visc_P = 0.001					UNITS u_Pas			"Viscosity of combusion products"
		REAL cond_P = 0.001					UNITS u_W_mK		"Thermal conductivity of combustion products"
		REAL T_boil_P = 200					UNITS u_K			"Boiling temperature of combustion products"
		
	DECLS
		
		REAL T_in								UNITS u_K			"Average inlet temperature"
		REAL T_c	= 4000						UNITS u_K			"Combustion temperature"
		REAL p_c									UNITS u_Pa			"Chamber pressure"
		
		
		REAL eta									UNITS no_units		"Combustion efficiency"
		REAL phi									UNITS no_units		"Equivalence ratio"
		REAL W_F_st								UNITS u_kg_s		"Fuel mass flow that provides a stoichiometric mixture when combined with the oxidant mass flow"
		REAL OF									UNITS no_units		"Mixture ratio"
		ALG REAL W_F = 0.5					UNITS u_kg_s		"Fuel mass flow"
		REAL W_O									UNITS u_kg_s		"Oxidant mass flow"
		REAL W_IF								UNITS u_kg_s		"Inert mass flow through the fuel port"
		REAL W_IO								UNITS u_kg_s		"Inert mass flow through the oxidant port"
		
		HIDDEN REAL fluid_O[ChemName]		UNITS no_units		"Oxidant fluid"		
		HIDDEN REAL fluid_F[ChemName]		UNITS no_units		"Fuel fluid"
		HIDDEN REAL fluid_P[ChemName]		UNITS no_units		"Combustion products fluid"
		HIDDEN REAL fluid_P_vap[ChemName]		UNITS no_units		"Vaporised combustion products fluid"
		
		REAL cp_R								UNITS u_J_kgK		"Reactant specific heat at constant pressure"
		
		BOOLEAN Combustion											"TRUE if there is combustion and FALSE if either oxidant or fuel is lacking"	
		REAL Q_comb_effective				UNITS u_J_kg		"Effective heat of combustion per oxidant mass flow unit"
		
		REAL c_star								UNITS u_m_s			"Characteristic exhaust velocity"	
		
	INIT PRIORITY 100
		Init_fluid(Comb_prod, g.fluid)
		Vaporisation_mix(fluid_P, fluid_P_vap)
		Combustion = TRUE
		
	DISCRETE
		WHEN (SUM (i IN LiquidsGases ; fluid_O[i] * fluid_F[i]) != 0) THEN
			Combustion = FALSE
		END WHEN
		
		WHEN (SUM (i IN LiquidsGases ; fluid_O[i] * fluid_F[i]) == 0) THEN
			Combustion = TRUE
		END WHEN
		
	CONTINUOUS		
		eta = eta_d
		
		p_c = f_O.p_c
		p_c = f_F.p_c
		g.pt = p_c * TPL
		
		f_O.W + f_F.W = g.W
		W_O = f_O.W * (1 - f_O.fluid[Comb_prod])
		W_F = f_F.W * (1 - f_F.fluid[Comb_prod])
		W_IO = f_O.W - W_O
		W_IF = f_F.W - W_F
		OF = W_O / W_F
		phi = OF_st / OF
		W_F_st = W_F / phi

		
		Q_comb_effective = 	IF (Combustion)		Q_comb
									ELSE						0
								
		eta * Q_comb_effective = ((1. + OF) / min(OF,OF_st)) * (cp(fluid_P_vap) * (T_c - T_ref) - cp_R * (T_in - T_ref))
		(1 + phi / OF_st) * cp_R * T_in = cp(fluid_O) * f_O.T + phi / OF_st * cp(fluid_F) * f_F.T		
		
		(W_O + W_F)*fluid_P_vap[Comb_prod]/ g.W * fluid_P_vap[Comb_prod_cp] * (g.Tt - T_c)                         \
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
		fluid_P[Comb_prod_T_boil] = IF (Combustion) T_boil_P
										ELSE 					0
		fluid_P[Comb_prod_rho] = 0							
								
	
		
		--Vaporisation_mix(fluid_P, fluid_P_vap)
		EXPAND (i IN ChemName EXCEPT LV_extended) fluid_P_vap[i] = fluid_P[i]
		EXPAND (i IN LV) fluid_P_vap[i] = 0
		fluid_P_vap[O2] = fluid_P[LOX] + fluid_P[O2]
		fluid_P_vap[H2] = fluid_P[LH2] + fluid_P[H2]
		fluid_P_vap[CH4] = fluid_P[LCH4] + fluid_P[CH4]
		fluid_P_vap[MMH_vapour] = fluid_P[MMH] + fluid_P[MMH_vapour]
		
		
		
		EXPAND (i IN LiquidsGases) g.fluid[i] = fluid_P_vap[i] * (W_O + W_F) / g.W
		g.fluid[Comb_prod] = (fluid_P_vap[Comb_prod] * (W_O + W_F) + (W_IO + W_IF)) / g.W		
		EXPAND (i IN Comb_prod_prop EXCEPT Comb_prod) g.fluid[i] = fluid_P_vap[i]
		
		c_star = sqrt(R(g.fluid) * g.Tt) / FGAMMA(g.fluid)
		i.Data[1] = c_star		

END COMPONENT
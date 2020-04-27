--------------------------------------------------------------------------------
-- EcosimPro Simulation Source Code
-- FILE NAME: LPRES_Common.el
-- DESCRIPTION: Defines basic components
-- NOTES:
-- AUTHOR: Pablo Sierra Heras
-- CREATION DATE: 05/12/2014
--------------------------------------------------------------------------------

-- Libraries
USE MATH VERSION "3.1"


-- Global constants
CONST REAL p_std = 101325.						UNITS u_Pa			"Standard pressure"
CONST REAL T_std = 288.15 						UNITS u_K			"Standard temperature"

CONST REAL g_0 = 9.80665    					UNITS u_m_s2		"Gravity acceleration at the Earth surface"

CONST REAL R_u = 8314.							UNITS u_J_kmolK	"Universal Gas Constant"

CONST REAL T_ref = 298.15						UNITS u_K			"Reference temperature"
CONST REAL math_zero = 1e-6					UNITS no_units    "Zero to avoid calculation errors"


-- Global variables
BOUND REAL Altitude = 0.						UNITS u_m			"Geometric altitude"


-- Global enumerations
ENUM ChemName = {LOX, NTO, H2O2, NH3, HNO3, LF2, RP_1, LCH4, LH2, N2H4, UDMH, MMH, JP_10, Kerox, Oil, H2O, IPA, Air, Ar, CH4, CO, CO2, H2, He, N2, O2, MMH_vapour, H2O_vapour, Comb_prod, Comb_prod_M, Comb_prod_cp, Comb_prod_c, Comb_prod_cp_g, Comb_prod_Lv, Comb_prod_visc, Comb_prod_cond, Comb_prod_T_boil, Comb_prod_rho}		\
								"Names of available chemicals"
SET_OF(ChemName) Comb_prod_prop = {Comb_prod, Comb_prod_M, Comb_prod_cp, Comb_prod_c, Comb_prod_cp_g, Comb_prod_Lv, Comb_prod_visc, Comb_prod_cond, Comb_prod_T_boil, Comb_prod_rho} "Properties of combustion products"
SET_OF(ChemName) Liquids = {LOX, NTO, H2O2, HNO3, LF2, RP_1, LCH4, LH2, N2H4, UDMH, MMH, JP_10, Kerox, Oil, H2O, IPA}						"Names of available liquids"
SET_OF(ChemName) Gases = {Air, Ar, CH4, CO, CO2, H2, He, N2, O2, MMH_vapour, NH3, H2O_vapour}																						"Names of available gases"
SET_OF(ChemName) LiquidsGases = {LOX, NTO, H2O2, HNO3, LF2, RP_1, LCH4, LH2, N2H4, UDMH, MMH, JP_10, Kerox, Oil, H2O, IPA, Air, Ar, CH4, CO, CO2, H2, He, N2, O2, MMH_vapour, NH3, H2O_vapour}			\
								"Names of available chemicals except Comb_prod"
SET_OF(ChemName) Mono_LiquidsGases = {N2H4, NH3, N2, H2, H2O_vapour, O2}	"Liquids and gases considered in Monopropellants"							
SET_OF(ChemName) LV = {LOX, LCH4, LH2, MMH}						"Liquids that can be vaporised"
SET_OF(ChemName) LV_gases={O2, CH4, H2, MMH_vapour}
SET_OF(ChemName) LV_extended = {LOX, O2, LCH4, CH4, LH2, H2, MMH, MMH_vapour}			"Liquids that can be vaporised and their respective gases"

SET_OF(ChemName) Fuels = {LH2,H2}				"List of fuels"
SET_OF(ChemName) Oxidants  ={LOX, O2}			"List of oxidants"

ENUM ChemState = {Liquid, Gas}							"State of available chemicals"
ENUM Decomposition ={CombProd, CombProd_cp, RFuel, RFuel_cp, ROxidant, ROxidant_cp} "Parts of combustion decomposition"


ENUM ConDiv = {Convergent, Divergent}
ENUM YesNo = {Yes, No}
ENUM Type_Inlet = {All, Unknown_W}
ENUM Type_Area = {Area, No_Area}
ENUM Type_Heat_Ex={Liquid2Liquid, Gas2Gas, Liquid2Gas, Gas2Liquid}
SET_OF(Type_Heat_Ex) Type_CJ = {Liquid2Liquid, Liquid2Gas}
ENUM Type_Tank = {Pressure_fed, Blowdown}
ENUM Reaction_mono = {H_Peroxide, Hidrazine, Any}

ENUM Type_All = {Design, Known_pi, Known_W, Off_design, Known_pt_out, Known_dp, Known_dpr, Darcy, Gas_valve, Liq_valve, Known_Qp, Known_k, Compressible, Incompressible_l, Incompressible_g}
SET_OF(Type_All) OnOffDesign = {Design, Off_design}
SET_OF(Type_All) Type_Turbines = {Known_pi, Known_W, Off_design}
SET_OF(Type_All) Type_Turbine_liq = {Known_dp, Known_W, Off_design}
SET_OF(Type_All) Type_Regulator2 = {Known_pt_out, Known_dp, Known_dpr, Gas_valve, Liq_valve}
SET_OF(Type_All) Type_Regulator = {Design, Known_pt_out, Known_dp, Known_dpr}
SET_OF(Type_All) Type_Cooling = {Darcy, Known_dp}
SET_OF(Type_All) Type_Deposit = {Known_Qp, Known_k}
SET_OF(Type_All) Type_Junction = {Compressible, Incompressible_g, Incompressible_l}
SET_OF(Type_All) Type_Pipe = {Darcy,Known_dp, Known_dpr}
ENUM AngCoef = {Angles, Coefficients}

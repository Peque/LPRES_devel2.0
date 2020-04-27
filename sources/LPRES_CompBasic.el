--------------------------------------------------------------------------------
-- EcosimPro Simulation Source Code
-- FILE NAME: LPRES_CompBasic.el
-- DESCRIPTION: Defines basic components
-- NOTES:
-- AUTHOR: Pablo Sierra Heras
-- CREATION DATE: 05/12/2014
--------------------------------------------------------------------------------

-- Libraries
USE MATH VERSION "3.1"


-- Abstract components
--------------------------------------------------------------------------------
-- Abstract component for definition of components with a fluid inlet and a
-- fluid outlet
--------------------------------------------------------------------------------
ABSTRACT COMPONENT FluidInFluidOut

  "Abstract component for definition of components with a fluid inlet and a fluid outlet"

   PORTS
      IN Fluid f_in						"Inlet fluid port"
      OUT Fluid f_out					"Outlet fluid port"

   CONTINUOUS	
		f_in.W = f_out.W
		f_in.fluid = f_out.fluid

END COMPONENT

--------------------------------------------------------------------------------
-- Abstract component for definition of gas turbomachinery components
--------------------------------------------------------------------------------
ABSTRACT COMPONENT GasTurbo IS_A FluidInFluidOut

   "Abstract component for definition of gas turbomachinery components"
	
	PORTS      
      IN Mechanical m				"Inlet mechanical port"

   DECLS
      REAL Power			UNITS u_W			"Mechanical power"	

   CONTINUOUS
	   Power = f_in.W * cp(f_in.fluid) * (f_in.Tt - f_out.Tt)
		
		m.Power = Power

END COMPONENT


-- Components
--------------------------------------------------------------------------------
-- Component that simulates a liquid discharge to the atmosphere
--------------------------------------------------------------------------------
COMPONENT Ambient 

	"Liquid discharge to the atmosphere"
	
	PORTS
		IN Fluid l									"Inlet fluid port"	
		
	DATA
		REAL A = 0.01						UNITS u_m2		"Discharge area"
		
	DECLS		
		REAL p_amb 							UNITS u_Pa		"Ambient pressure"
		
	DISCRETE
		ASSERT (Liquid == State(l.fluid)) FATAL "ONLY LIQUIDS CAN GO THROUGH THE AMBIENT!"
						
	CONTINUOUS
		p_amb = ISA_Pressure(Altitude)
		
		l.pt = p_amb + (l.W / A)**2 / (2. * rho(l.fluid))	
				
END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a turbine with a choked inlet
--------------------------------------------------------------------------------
COMPONENT Turbine_ch IS_A GasTurbo

   "Turbine with a choked inlet"
		
	DATA
		REAL eta_d = 0.5			UNITS no_units		"Design efficiency"	
		REAL A_in = 0.001			UNITS u_m2			"Input area"
		
		
		
	DECLS
		REAL eta				UNITS no_units		"Efficiency"
		REAL alpha			UNITS no_units		"Total temperature ratio"
		REAL rpm				UNITS "rpm"			"Design rotational speed [rpm]"
		REAL pi				UNITS no_units		"Design expansion ratio"
		REAL tau				UNITS u_J_kg			"Specific work"
		REAL W				UNITS u_kg_s		"Design mass flow"	
		REAL MF_param				UNITS "Kg/s*K^2/Pa"	"Mass flow parameter"
		
	DISCRETE
		ASSERT (Gas == State(f_in.fluid)) FATAL "ONLY GASES CAN GO THROUGH THE TURBINE!"
	
	CONTINUOUS		
		alpha = 1. - eta * (1. - (f_out.pt / f_in.pt)**((gamma(f_in.fluid) - 1.) / gamma(f_in.fluid)))		
		alpha = f_out.Tt / f_in.Tt
		
		FGAMMA(f_in.fluid) = f_in.W * sqrt(f_in.Tt * R(f_in.fluid)) / (f_in.pt * A_in) --cos(alpha_2 *180*PI) removed
		
		eta = eta_d		
		
		m.N = rpm * 2. * PI / 60.
		pi = f_in.pt / f_out.pt		

		tau = Power/W

		f_in.W = W
		
		MF_param =  W *sqrt(f_in.Tt)/f_in.pt
		
END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a turbine (the inlet may be choked or not)
--------------------------------------------------------------------------------
COMPONENT Turbine IS_A GasTurbo (ENUM AngCoef Type_AC = Coefficients)

   "Turbine (the inlet may be choked or not)"
		
	DATA
		REAL eta_d = 0.8					UNITS no_units		"Design efficiency"
		REAL phi_d = 0.05					UNITS no_units		"Design flow coefficient if Type_AC=Coefficients"
		REAL psi_d = 0.7					UNITS no_units		"Design loading coefficient if Type_AC=Coefficients"
		REAL alpha_2 = 45.				UNITS u_deg			"Flow angle in section 2"
		REAL alpha_4r = -30.				UNITS u_deg			"Relative flow angle in section 4 if Type_AC=Angles"
		REAL A_in = 0.005					UNITS u_m2			"Input area"
		REAL r_m	= 0.01					UNITS u_m			"Average radius"
		
	DECLS
		REAL eta						UNITS no_units		"Efficiency"
		REAL psi						UNITS no_units		"Loading coefficient"
		REAL phi						UNITS no_units		"Flow coefficient"
		
		ALG REAL U	= 500			UNITS u_m_s			"Blade speed"
		REAL alpha					UNITS no_units		"Total temperature ratio"	
		REAL tau						UNITS u_J_kg		"Specific work"
		REAL W						UNITS u_kg_s		"Mass flow"	
		
		REAL M_2	= 0.9				UNITS no_units		"Mach number in section 2"
		REAL V_z2					UNITS u_m_s			"Axial speed in section 2"
		REAL V_2						UNITS u_m_s			"Speed in section 2"
		REAL rpm						UNITS "rpm"			"Blade speed in rpm"
		REAL pi 						UNITS no_units		"Design expansion ratio"
		REAL MF_param				UNITS "Kg/s*K^0.5/Pa"	"Mass flow parameter"
		
		
	DISCRETE
		ASSERT (Gas == State(f_in.fluid)) FATAL "ONLY GASES CAN GO THROUGH THE TURBINE!"
	
	CONTINUOUS		
		alpha = 1. - eta * (1. - (f_out.pt / f_in.pt)**((gamma(f_in.fluid) - 1.) / gamma(f_in.fluid)))		
		alpha = f_out.Tt / f_in.Tt	

		f_in.W * sqrt(f_in.Tt * R(f_in.fluid)) / (f_in.pt * A_in) = 	sqrt(gamma(f_in.fluid))* min(M_2, 1.) * (1. + (gamma(f_in.fluid) - 1.) / 2. *  min(M_2, 1.)**2)**(- (gamma(f_in.fluid) + 1.) / 2. / (gamma(f_in.fluid) - 1.))
																						
		U = m.N * r_m
		psi = tau / U**2
		phi = V_z2 / U
		
		V_z2 = V_2 * cos(alpha_2 * PI / 180.)
		V_2 = M_2 * sqrt(gamma(f_in.fluid) * R(f_in.fluid) * f_in.Tt / (1. + (gamma(f_in.fluid) - 1.) / 2. * M_2**2))
		IF (Type_AC == Coefficients) INSERT
		psi = (1. + psi_d) / phi_d * phi - 1.
		ELSEIF (Type_AC == Angles) INSERT
			psi = phi * (tan(alpha_2 * PI / 180.) - tan(alpha_4r * PI / 180.)) - 1.	
		END IF			
		eta = eta_d
		
		tau = Power / f_in.W
		
		pi = f_in.pt / f_out.pt
		m.N = rpm * 2. * PI / 60.
		f_in.W = W
		
		MF_param =  W *sqrt(f_in.Tt)/f_in.pt
		
		
END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a compressor
--------------------------------------------------------------------------------
COMPONENT Compressor IS_A GasTurbo (ENUM AngCoef Type_AC = Coefficients)

   "Compressor"
			
	DATA
		REAL eta_d = 0.8			UNITS no_units		"Design efficiency"
		REAL phi_d = 0.05			UNITS no_units		"Design flow coefficient"
		REAL psi_d = 0.7			UNITS no_units		"Design loading coefficient if Type_AC=Coefficients"
		REAL alpha_2r = -30.		UNITS u_deg			"Relative current angle in section 2 if Type_AC=Angles"
		REAL r_m = 0.3				UNITS u_m			"Average radius"
		REAL A_in = 0.35			UNITS u_m2			"Input area"
		
		
	DECLS	
		
		REAL eta						UNITS no_units		"Efficiency"
		REAL psi						UNITS no_units		"Loading coefficient"
		REAL phi						UNITS no_units		"Flow coefficient"
		REAL pi						UNITS no_units		"Design compression ratio"
		
		ALG REAL U=500				UNITS u_m_s			"Blade speed"
		REAL tau						UNITS u_J_kg		"Specific work"
		
		REAL rho_in					UNITS u_kg_m3		"Inlet density"
		ALG REAL M_in = 0.001	UNITS no_units		"Inlet Mach number"
		REAL MF_param				UNITS "Kg/s*K^2/Pa"	"Mass flow parameter"
		REAL W						UNITS u_kg_s      "Mass flow"
		
	DISCRETE
		ASSERT (Gas == State(f_in.fluid)) FATAL "ONLY GASES CAN GO THROUGH THE COMPRESSOR!"		
		ASSERT (M_in <= 1.) KILLPOINT "THE COMPRESSOR INLET CANNOT BE CHOKED!"
		
	CONTINUOUS
		(f_out.Tt / f_in.Tt)  = ((f_out.pt / f_in.pt)**((gamma(f_in.fluid) - 1.) / gamma(f_in.fluid)) - 1.) / eta + 1.
		
		U = m.N * r_m		
		psi = tau / U**2
		phi = f_in.W /(A_in * U * rho_in)
		
		f_in.W * sqrt(R(f_in.fluid) * f_in.Tt) / A_in / f_in.pt = sqrt(gamma(f_in.fluid)) * M_in * (1. + (gamma(f_in.fluid) - 1.) / 2. * M_in**2)**(- (gamma(f_in.fluid) + 1.) / 2. / (gamma(f_in.fluid) - 1.))
		((f_in.pt / R(f_in.fluid) / f_in.Tt) / rho_in)**(gamma(f_in.fluid) - 1.) = 1. + (gamma(f_in.fluid) - 1.) / 2. * M_in**2		
		
		IF (Type_AC == Coefficients) INSERT
			psi = 1. - (1. - psi_d) / phi_d * phi
		ELSE
			psi = phi * tan(alpha_2r * PI / 180.) + 1.
		END IF			
		eta = eta_d
				
		tau = - Power / f_in.W
		
		pi = f_out.pt / f_in.pt
		
		f_in.W = W
		MF_param =  W *sqrt(f_in.Tt)/f_in.pt
				
END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a liquid turbine
--------------------------------------------------------------------------------
COMPONENT Turbine_liq IS_A FluidInFluidOut --(ENUM Type_Turbine_liq Type = Known_dp)

	"Liquid turbine"
	
	PORTS      
      IN Mechanical m				"Inlet mechanical port"
		
	DATA
		REAL eta_d = 0.8			UNITS no_units		"Design efficiency"
		REAL phi_d = 0.05			UNITS no_units		"Design flow coefficient"
		REAL psi_d = 0.7			UNITS no_units		"Design loading coefficient"
		REAL A_in = 0.01			UNITS u_m2			"Input area"
		REAL r_m	= 0.05			UNITS u_m			"Average radius"
		
	DECLS		
		
		REAL eta						UNITS no_units		"Efficiency"
		REAL psi						UNITS no_units		"Loading coefficient"
		REAL phi						UNITS no_units		"Flow coefficient"	
		REAL dp 						UNITS u_Pa			"Design pressure decrease"
		REAL W						UNITS u_kg_s		"Design mass flow"	
		REAL rpm 					UNITS "rpm"			"Design rotational speed [rpm]"
		
		ALG REAL U					UNITS u_m_s			"Blade speed"	
		REAL tau						UNITS u_J_kg		"Specific work"
      REAL Power					UNITS u_W			"Mechanical power"			


	DISCRETE
		--ASSERT (Liquid == State(f_in.fluid)) FATAL "ONLY LIQUIDS CAN GO THROUGH THE LIQUID TURBINE!"
		
	CONTINUOUS
		tau = (f_in.pt - f_out.pt) * eta / rho(f_in.fluid)
		--f_out.pt = f_in.pt - tau * rho(f_in.fluid) / eta
		
		U = m.N * r_m	
		psi = tau / U**2
		phi = f_in.W /(A_in * U * rho(f_in.fluid))
			
		psi = (1. + psi_d) / phi_d * phi - 1.
		eta = eta_d
		
		tau = Power / f_in.W		
				
		m.N = rpm * 2. * PI / 60.
		dp = f_in.pt - f_out.pt
		f_in.W = W
		
	   dp / rho(f_in.fluid) * (1. - eta) = cp(f_in.fluid) * (f_out.Tt - f_in.Tt)
		
		m.Power = Power
				
END COMPONENT	

--------------------------------------------------------------------------------
-- Component that represents a pump
--------------------------------------------------------------------------------
COMPONENT Pump IS_A FluidInFluidOut

	"Pump"
	
	PORTS      
      IN Mechanical m				"Inlet mechanical port"
		
	DATA
		REAL eta_d = 0.8			UNITS no_units		"Design efficiency"
		REAL phi_d = 0.05			UNITS no_units		"Design flow coefficient"
		REAL psi_d = 0.7			UNITS no_units		"Design loading coefficient"
		REAL A_in = 0.01			UNITS u_m2			"Input area"
		REAL r_m	= 0.05			UNITS u_m			"Average radius"
		
	DECLS		
		
		REAL eta						UNITS no_units		"Efficiency"
		REAL psi						UNITS no_units		"Loading coefficient"
		REAL phi						UNITS no_units		"Flow coefficient"
		REAL dp 						UNITS u_Pa			"Pressure increase"
		
		ALG REAL U = 500			UNITS u_m_s			"Blade speed"	
		REAL H						UNITS u_m			"Head"
		REAL tau						UNITS u_J_kg		"Specific work"
      REAL Power					UNITS u_W			"Mechanical power"
		
	DISCRETE
		ASSERT (Liquid == State(f_in.fluid)) WARNING "ONLY LIQUIDS CAN GO THROUGH THE PUMP!"
		
	CONTINUOUS
		H = (f_out.pt - f_in.pt) / (g_0 * rho(f_in.fluid))
		tau = (f_out.pt - f_in.pt) / (eta * rho(f_in.fluid))
		
		U = m.N * r_m		
		psi = tau / U**2
		phi = f_in.W /(A_in * U * rho(f_in.fluid))
			
		psi = 1. - (1. - psi_d) / phi_d * phi
		eta = eta_d
		
		tau = - Power / f_in.W	
		dp = f_out.pt - f_in.pt
		
	   (f_out.pt - f_in.pt) / rho(f_in.fluid) * (1. / eta - 1.) = cp(f_in.fluid) * (f_out.Tt - f_in.Tt)
		
		m.Power = Power
		
				
END COMPONENT	

--------------------------------------------------------------------------------
-- Component that represents a liquid pipe with pressure drop
--------------------------------------------------------------------------------
COMPONENT Pipe IS_A FluidInFluidOut (ENUM Type_Pipe Type = Darcy)

   "Liquid pipe with pressure drop"
	
	DATA	
		REAL L = 1.				UNITS u_m 			"Length"	
		REAL D = 0.1			UNITS u_m 			"Diameter"	
		REAL K = 5				UNITS no_units		"Additional pressure losses"
		REAL rug	= 1.5e-6		UNITS u_m			"Absolute rugosity"
		REAL dp = 10000		UNITS u_Pa			"Imposed pressure drop"
		REAL dpr = 0.05			UNITS no_units		"Imposed pressure drop ratio"
		
	DECLS
		REAL v					UNITS u_m_s				"Liquid speed"
		REAL f					UNITS no_units			"Darcy friction factor"		
		REAL Re					UNITS no_units			"Reynolds number"
		
	DISCRETE
		--ASSERT ((Liquid == State(f_in.fluid) AND Type == Darcy)OR(Type == Known_dp)OR(Type == Known_dpr)) FATAL "ONLY LIQUIDS CAN GO THROUGH THE PIPE!"
	CONTINUOUS
	
		IF (Type == Darcy) INSERT
		f_out.pt = f_in.pt - (f * L / D + K) * 0.5 * rho(f_in.fluid) * v**2
		f = hdc_fric(D, rug, Re)
		
		v = f_in.W / (PI * D**2 / 4 * rho(f_in.fluid))
		Re = rho(f_in.fluid) * v * D / visc(f_in.fluid)
		END IF
		
		IF (Type == Known_dp) INSERT
			f_out.pt = f_in.pt - dp
			v=0
			f=0
			Re=0
		END IF
		
		IF (Type == Known_dpr) INSERT
			f_out.pt = f_in.pt * (1 - dpr)
			v=0
			f=0
			Re=0
		END IF	
		
		f_out.Tt = f_in.Tt
		
END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a flow splitter with pressure drop
--------------------------------------------------------------------------------
COMPONENT SplitFrac

   "Flow splitter with pressure drop"

	PORTS
      IN Fluid f_in                 				"Inlet fluid port"
      OUT Fluid f_out               				"Outlet fluid port"
      OUT Fluid f_b           						"Branch fluid port"
		
	DATA	
		REAL TPL = 0.9			UNITS no_units			"Total pressure loss"

   CONTINUOUS	
		f_in.W = f_out.W + f_b.W
		
		f_in.Tt = f_out.Tt
		f_in.Tt = f_b.Tt
		
		f_in.pt = f_out.pt / TPL
		f_in.pt = f_b.pt / TPL 
		
		f_in.fluid = f_out.fluid
		f_in.fluid = f_b.fluid
	
END COMPONENT


--------------------------------------------------------------------------------
-- Component that represents a pressure regulator
--------------------------------------------------------------------------------
COMPONENT Regulator IS_A FluidInFluidOut

   "Pressure regulator"
		
	DATA
		REAL dp_min = 1500. 			UNITS u_Pa			"Minumum total pressure drop"
		REAL pt_out = 1200000.		UNITS u_Pa			"Outlet total pressure"
		
   CONTINUOUS	
		f_out.Tt = f_in.Tt
		f_out.pt = min(pt_out, f_in.pt - dp_min)
		
	
END COMPONENT

--------------------------------------------------------------------------------
-- Component that specifies the conditions of a fluid inlet in the system
--------------------------------------------------------------------------------
COMPONENT Inlet 

	"Conditions of a fluid inlet in the system"
	
	PORTS
		OUT Fluid f						"Outlet fluid port"
	
	DATA
		ENUM LiquidsGases fluid = LOX											"Working fluid name"
  		
	DECLS
		BOUND REAL Tt 									UNITS u_K				"Total temperature" 
		BOUND REAL pt 									UNITS u_Pa				"Total pressure"
		REAL W 											UNITS u_kg_s			"Mass flow"	
		

	INIT PRIORITY 100
		Init_fluid(fluid, f.fluid)

	CONTINUOUS
		Init_fluid(fluid, f.fluid)
		f.pt = pt
		f.Tt = Tt
		f.W = W

END COMPONENT

COMPONENT Outlet 

	"Conditions of a fluid inlet in the system"
	
	PORTS
		IN Fluid f						"Outlet fluid port"
	
	DATA
	DECLS
		BOUND REAL Tt 									UNITS u_K				"Total temperature" 
		BOUND REAL pt 									UNITS u_Pa				"Total pressure"
		REAL W 											UNITS u_kg_s			"Mass flow"	

	INIT PRIORITY 100

	CONTINUOUS
		f.pt = pt
		f.Tt = Tt
		f.W = W
	

END COMPONENT
--------------------------------------------------------------------------------
-- Component that represents a liquid tank
--------------------------------------------------------------------------------
COMPONENT Tank (ENUM Type_Area Type = Area)

	"Liquid tank"
	
	PORTS
		IN Fluid g						"Inlet fluid port"
		OUT Fluid l						"Outlet fluid port"
	
	DATA		
		ENUM Liquids fluid_l	= LOX									"Working liquid name"
		REAL T_d = 288.15					UNITS u_K				"Tank temperature"		
		REAL Ag = 0.001					UNITS u_m2				"Gas input area"
		
		
	DECLS		
		REAL p_g								UNITS u_Pa			"Pressurisation gas pressure"
		REAL rho_g							UNITS u_kg_m3		"Pressurisation gas density"
		ALG REAL M_g = 0.1				UNITS no_units		"Pressurisation gas Mach number"
		REAL p_d								UNITS u_Pa			"Design tank pressure"
		
		
	INIT PRIORITY 100
		Init_fluid(fluid_l, l.fluid)
		
	DISCRETE
		ASSERT (Gas == State(g.fluid)) FATAL "ONLY GASES CAN PRESSURISE THE TANK!"
		ASSERT (M_g <= 1.) KILLPOINT "THE PRESSURISATION GAS INLET CANNOT BE CHOKED!"
			
	CONTINUOUS
		Init_fluid(fluid_l, l.fluid)
		l.Tt = T_d
		
		l.pt = p_g
		
		l.W / rho(l.fluid) = g.W / rho_g
		IF (Type == Area) INSERT
			g.W * sqrt(R(g.fluid) * g.Tt) / Ag / g.pt = sqrt(gamma(g.fluid)) * M_g * (1. + (gamma(g.fluid) - 1.) / 2. * M_g**2)**(- (gamma(g.fluid) + 1.) / 2. / (gamma(g.fluid) - 1.))
			((g.pt / R(g.fluid) / g.Tt) / rho_g)**(gamma(g.fluid) - 1.) = 1. + (gamma(g.fluid) - 1.) / 2. * M_g**2
			(g.pt / p_g)**((gamma(g.fluid) - 1.) / gamma(g.fluid)) = 1. + (gamma(g.fluid) - 1.) / 2. * M_g**2	
		ELSEIF (Type == No_Area) INSERT 		
			rho_g = g.pt / ( R(g.fluid) * g.Tt )
			p_g = g.pt
			M_g = 0
		END IF
		
		p_g = p_d
END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a liquid tank that is pressurised by the atmosphere
--------------------------------------------------------------------------------
COMPONENT TankOpen

	"Liquid tank that is pressurised by the atmosphere"
	
	PORTS
		OUT Fluid l				"Outlet fluid port"
	
	DATA
		ENUM Liquids fluid = LOX									"Working liquid name"
		REAL T_d = 288.15					UNITS u_K				"Tank temperature"
		
	DECLS		
		REAL p_d						UNITS u_Pa		"Tank pressure"
		
	INIT PRIORITY 100
		Init_fluid(fluid, l.fluid)
			
	CONTINUOUS
		Init_fluid(fluid, l.fluid)
		l.Tt = T_d
		
		l.pt = p_d
		
		p_d = ISA_Pressure(Altitude)

END COMPONENT


--------------------------------------------------------------------------------
-- Component that represents a mechanical shaft without acceleration
--------------------------------------------------------------------------------
COMPONENT Shaft

   "Mechanical shaft without acceleration"

   PORTS
      OUT Mechanical m_1				"Outlet mechanical port"
      OUT Mechanical m_2				"Outlet mechanical port"
		
	DATA
		REAL eta	= 1.	UNITS no_units		"Efficiency"
	
	DECLS
		REAL rpm			UNITS "rpm"			"Rotational speed [rpm]"

   CONTINUOUS
      m_1.N = m_2.N		
		m_1.N = rpm * 2. * PI / 60.
		
		m_1.Power = ZONE (m_2.Power > 0)	- m_2.Power * eta
						OTHERS					- m_2.Power / eta

END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a mechanical shaft with a gearbox and without 
-- acceleration
--------------------------------------------------------------------------------
COMPONENT Gearbox

   "Mechanical shaft with a gearbox and without acceleration"

   PORTS
      OUT Mechanical m_in				"Outlet mechanical port"
      OUT Mechanical m_out				"Outlet mechanical port"
		
	DATA
		REAL GR = 1.	UNITS no_units		"Gear ratio"	
		REAL eta	= 1.	UNITS no_units		"Efficiency"

   CONTINUOUS
      m_out.N = GR * m_in.N
		
      m_out.Power = 	ZONE (m_in.Power > 0) 	- m_in.Power * eta
							OTHERS					 	- m_in.Power / eta

END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a convergent nozzle with the ambient already 
-- connected
--------------------------------------------------------------------------------
COMPONENT Nozzle 

	"Convergent nozzle with the ambient already connected"
	
	PORTS
		IN Fluid g						"Inlet fluid port"
		OUT Info (n = 1) i			"Outlet information port"
		
	DATA
		REAL A = 0.02			UNITS u_m2			"Discharge area"
		
	DECLS
		REAL p_amb 				UNITS u_Pa			"Ambient pressure"
		REAL A_sl				UNITS u_m2			"Sonic lock area"		
		
		REAL PR					UNITS no_units		"Pressure ratio"
		REAL PR_sl				UNITS no_units		"Sonic lock pressure ratio"		
		REAL p_out				UNITS u_Pa			"Outlet pressure"	
		REAL p_out_ch			UNITS u_Pa			"Choked outlet pressure"
			
		REAL M_out				UNITS no_units		"Outlet Mach number"	
		REAL T_out			   UNITS u_K			"Outlet temperature"
		REAL v_out				UNITS u_m_s			"Outlet speed"
		REAL Thrust				UNITS u_N			"Thrust"
		
	INIT
		g.pt = 10000000
		
	DISCRETE
		ASSERT (Gas == State(g.fluid)) FATAL "ONLY GASES CAN GO THROUGH THE NOZZLE!"	
						
	CONTINUOUS
		p_amb = ISA_Pressure(Altitude)
		
		PR = g.pt / p_amb
		PR_sl = ((gamma(g.fluid) + 1.) / 2.)**(gamma(g.fluid) / (gamma(g.fluid) - 1.))
		
		p_out_ch = g.pt / PR_sl
		
		g.W = A_sl * FGAMMA(g.fluid) * g.pt / sqrt(g.Tt * R(g.fluid))	
		
		A_sl = ZONE(M_out > 0) A * M_out / ((2. + (gamma(g.fluid) - 1.) * M_out**2) / (gamma(g.fluid) + 1.))**((gamma(g.fluid) + 1.) / (2. * (gamma(g.fluid) - 1.)))
				 OTHERS          A * (-M_out) / ((2. + (gamma(g.fluid) - 1.) * M_out**2) / (gamma(g.fluid) + 1.))**((gamma(g.fluid) + 1.) / (2. * (gamma(g.fluid) - 1.)))
		
		M_out =	ZONE (PR < PR_sl)		sqrt(2. * (PR**((gamma(g.fluid) - 1.) / gamma(g.fluid)) - 1.) / (gamma(g.fluid) - 1.))
					OTHERS					1.
					
		p_out = max(p_amb, p_out_ch)
							
		
		T_out = g.Tt / (1. + (gamma(g.fluid) - 1.) / 2. * M_out**2)
		v_out = M_out * sqrt(gamma(g.fluid) * R(g.fluid) * T_out)
		
		Thrust = g.W * v_out + A * (p_out - p_amb)
		
		i.Data[1] = Thrust
		
END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a choked convergent-divergent nozzle
--------------------------------------------------------------------------------
COMPONENT NozzleConDiv

   "Choked Convergent-divergent nozzle"
	
	PORTS
		IN Fluid g_in						"Inlet fluid port"
		OUT GasNozzle g_out				"Outlet gas port through a nozzle"
		
	DATA
		REAL AR = 10.					UNITS no_units		"Area ratio"
		REAL A_th = 0.05				UNITS u_m2			"Throat area"
		
		
	DECLS
		REAL A_out						UNITS u_m2			"Output area"	
		REAL W							UNITS u_kg_s		"Mass flow"
				
		ALG REAL p_out_ch = 100.	UNITS u_Pa			"Choked outlet pressure"
		
	DISCRETE
		ASSERT (Gas == State(g_in.fluid)) FATAL "ONLY GASES CAN GO THROUGH THE CHOKED CONVERGENT-DIVERGENT NOZZLE!"	
		
	CONTINUOUS
		g_in.fluid = g_out.fluid
		g_in.W = g_out.W
		
		g_out.A_out = A_out 
		AR = A_out / A_th
		
		FGAMMA(g_in.fluid) = g_in.W * sqrt(g_in.Tt * R(g_in.fluid)) / (g_in.pt * A_th) 
		g_in.W = W
		
		AR = FGAMMA(g_out.fluid) / ((p_out_ch / g_in.pt)**(1. / gamma(g_out.fluid)) * sqrt(2 * gamma(g_out.fluid) * (1 - (p_out_ch / g_in.pt)**((gamma(g_out.fluid) - 1.) / gamma(g_out.fluid))) / (gamma(g_out.fluid) - 1)))
		
		g_out.pt = g_in.pt 
		g_out.Tt = g_in.Tt
	
END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a nozzle extension
--------------------------------------------------------------------------------
COMPONENT NozzleExt (ENUM YesNo Cooled = No)
	
	"Nozzle extension"
	
	PORTS
		IN GasNozzle g_in				"Inlet gas port through a nozzle"
		OUT GasNozzle g_out			"Outlet gas port through a nozzle"
		OUT Heat h						"Outlet heat exchange port"
			
	DATA
		REAL AR = 2.									UNITS no_units		"Extension area ratio"
		
		REAL AR_r = 10. / 2.							UNITS no_units		"Area at the characteristic section of heat exchange divided by the throat area"			
		REAL A_wet = 1.								UNITS u_m2			"Nozzle wet area of the cooled zone"
		
	DECLS
		REAL A_in										UNITS u_m2			"Extension input area"
		REAL A_out 										UNITS u_m2			"Extension output area"
		
		REAL A_r											UNITS u_m2			"Area at the characteristic section of heat exchange"
		REAL M_r											UNITS no_units		"Mach number at the characteristic section of heat exchange"
		REAL Pr_r										UNITS no_units		"Prandtl number at the characteristic section of heat exchange"	
		REAL visc_r										UNITS u_Pas			"Dynamic viscosity at the characteristic section of heat exchange"
		REAL T_aw										UNITS u_K			"Adiabatic wall temperature of combustion gases"
		REAL h_g											UNITS u_W_m2K		"Combustion gases heat transfer coefficient"
		
	DISCRETE
		ASSERT (Gas == State(g_in.fluid)) FATAL "ONLY GASES CAN GO THROUGH THE NOZZLE EXTENSION!"
		
	CONTINUOUS
		g_in.fluid = g_out.fluid
		g_in.W = g_out.W

		g_in.A_out = A_in
		g_out.A_out = A_out
		AR = A_out / A_in
		
		g_out.pt = g_in.pt 
		g_out.Tt = g_in.Tt 
		g_out.c_star = g_in.c_star
		g_out.A_th = g_in.A_th
		
		-- Cooling of the nozzle
		
		T_aw = g_in.Tt * ((1. + Pr_r**0.33 * (gamma(g_in.fluid) - 1.) * M_r**2 / 2.) / (1. + (gamma(g_in.fluid) - 1.) * M_r**2 / 2.))
		Pr_r = 4. * gamma(g_in.fluid) / (9. * gamma(g_in.fluid) - 5.)
		visc_r = 1.184e-7 * M(g_in.fluid)**0.5 * T_aw**0.6		
		
		IF (Cooled == No) INSERT
			h.A = 0
			A_r = g_in.A_th
			M_r = 1.
			
			h_g = 0
			h.T = T_aw
			h.Q = 0
			
		ELSE		
			h.A = A_wet
			A_r = g_in.A_th * AR_r
			AR_r = 1 / M_r * FGAMMA(g_in.fluid) / sqrt(gamma(g_in.fluid)) * (1. + (gamma(g_in.fluid) - 1.) * M_r**2 / 2.)**((gamma(g_in.fluid) + 1.) / (2. * (gamma(g_in.fluid) - 1.)))
		
			-- Bartz correlation
			h_g = 0.026 / sqrt(4 * g_in.A_th / PI)**0.2 * (visc_r**0.2 * cp(g_in.fluid) / Pr_r**0.6) * (g_in.pt / g_in.c_star)**0.8 * (g_in.A_th / A_r)**0.9
			h.Q = h_g * (T_aw - h.T) * h.A
			
		END IF		
		
END COMPONENT


-- Measuring components
--------------------------------------------------------------------------------
-- Component that represents a monitor to measure the thrust
--------------------------------------------------------------------------------
COMPONENT ThrustMonitor

	"Monitor to measure the thrust"
	
	PORTS
		IN GasNozzle g						"Inlet gas port through a nozzle"
		OUT Info (n = 1) i				"Outlet information port"
		
	DECLS
		REAL Thrust						UNITS u_N			"Thrust"
		
		REAL p_amb 						UNITS u_Pa			"Ambient pressure"
		
		REAL T_out						UNITS u_K			"Outlet temperature"			
		REAL p_out						UNITS u_Pa			"Outlet pressure"	
		REAL v_out						UNITS u_m_s			"Outlet speed"
		REAL A_out						UNITS u_m2			"Output area"
		ALG REAL M_out = 100.		UNITS no_units		"Gas Mach number"
		
	DISCRETE
		ASSERT (Gas == State(g.fluid)) FATAL "ONLY GASES CAN GO THROUGH THE THRUST MONITOR!"
		ASSERT (M_out >= 1.) KILLPOINT "THE CONVERGENT-DIVERGENT NOZZLE MUST BE CHOKED!"
		ASSERT (p_out / p_amb > (1.88 * M_out - 1)**(-0.64)) WARNING "ACCORDING TO SCHMUCKER CRITERION, THE NOZZLE HAS A REGION OF DETACHMENT!"
						
	CONTINUOUS
		g.A_out = A_out
		p_amb = ISA_Pressure(Altitude)
		
		g.W * sqrt(R(g.fluid) * g.Tt) / A_out / g.pt = sqrt(gamma(g.fluid)) * M_out * (1. + (gamma(g.fluid) - 1.) / 2. * M_out**2)**(- (gamma(g.fluid) + 1.) / 2. / (gamma(g.fluid) - 1.))
		g.Tt / T_out = 1. + (gamma(g.fluid) - 1.) / 2. * M_out**2
		(g.pt / p_out)**((gamma(g.fluid) - 1.) / gamma(g.fluid)) = 1. + (gamma(g.fluid) - 1.) / 2. * M_out**2
		
		v_out = M_out * sqrt(gamma(g.fluid) * R(g.fluid) * T_out)	
		
		Thrust = g.W * v_out + A_out * (p_out - p_amb)
		
		i.Data[1] = Thrust

END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a mass flow meter
--------------------------------------------------------------------------------
COMPONENT FlowMeter IS_A FluidInFluidOut

	"Mass flow meter"
	
	PORTS
		OUT Info (n = 1) i							"Outlet information port"
						
	CONTINUOUS
		f_in.pt = f_out.pt
		f_in.Tt = f_out.Tt
		
		i.Data[1] = abs(f_in.W)

END COMPONENT

--------------------------------------------------------------------------------
-- Component that performs the calculations with the measurements done by other
-- components
--------------------------------------------------------------------------------
COMPONENT ControlPanel

	"Calculations with the measurements done by other components"
	
	PORTS
		IN Info (n = 1) i_W				"Inlet information port"
		IN Info (n = 1) i_Thrust		"Inlet information port"
		IN Info (n = 1) i_Comb			"Inlet information port"
		
	DECLS
		REAL Thrust				UNITS u_N			"Thrust"
		REAL W_tot 				UNITS u_kg_s		"Total mass flow"
		REAL c_star				UNITS u_m_s			"Characteristic exhaust velocity"
		
		REAL Isp					UNITS u_m_s			"Specific impulse [m/s]"
		REAL Isp_0				UNITS u_s			"Specific impulse [s]"
		REAL C_E					UNITS no_units		"Thrust coefficient"
						
	CONTINUOUS		
		i_Thrust.Data[1] = Thrust
		i_W.Data[1] = W_tot
		i_Comb.Data[1] = c_star
		
		Isp = Thrust / W_tot
		Isp = C_E * c_star
		Isp = g_0 * Isp_0

END COMPONENT





--------------------------------------------------------------------------------
-- Component that represents a gas deposit
--------------------------------------------------------------------------------


COMPONENT Deposit (ENUM Type_Deposit Type = Known_k)

	"Constant Volume Gas Deposit"
	
	PORTS 
		OUT Fluid g 									"Outlet information port"
		
	DATA
		ENUM Gases fluid = He						"Working gas name"
		REAL m_0	= 35					UNITS u_kg			"Initial gas mass in the deposit"
		REAL p_0	= 50e5				UNITS u_Pa			"Gas initial pressure"
		REAL T_0	= 273.15				UNITS u_K			"Gas initial temperature"
		REAL Qp_data = 0				UNITS u_W			"Heat per second data"	
		REAL k_data = 1.				UNITS no_units		"Polytropic exponent"
		ENUM t_interp intMet= AKIMA						"Interpolation method"
		
	DECLS
		REAL m							UNITS u_kg			"Gas mass in the deposit"
		REAL pt							UNITS u_Pa 			"Gas pressure"
		REAL Tt							UNITS u_K			"Gas temperature"
		DISCR REAL Vol_d				UNITS u_m3			"Deposit Volume"
		DISCR REAL Rgas				UNITS u_J_kgK		"Gas R constant"
		DISCR REAL Cp					UNITS u_J_kgK		"Specific heat constant pressure"
		DISCR REAL Cv					UNITS u_J_kgK		"Specific heat constant volume"
		REAL Q 							UNITS u_J			"Total heat added to the gas"
		DISCR REAL Q_0 = 0			UNITS u_J	 		"Initial heat"
		REAL Qp 							UNITS u_W			"Heat per second added to the gas"
		REAL k		 					UNITS no_units		"Polytropic exponent"
		REAL flux						UNITS u_J			"Total flux through outlet area"
		DISCR REAL flux_0 = 0 		UNITS u_J			"Initial flux through outlet area"
		REAL U							UNITS u_J			"Gas internal energy"
		DISCR REAL U_0					UNITS u_J			"Initial Gas internal energy"
		DISCR REAL Qp_num, xsaved[4], timePrev[4] = 0
		
	INIT PRIORITY 100
		Init_fluid(fluid, g.fluid)
		Rgas = R(g.fluid)
		Cp = cp(g.fluid)
		Cv = cv(g.fluid)
		Vol_d = m_0 * Rgas * T_0 / p_0
		m = m_0
		g.pt = p_0
		U_0 = Cv * m_0 * T_0
		Q=Q_0
		flux = flux_0
		
		xsaved[1]= 0
      xsaved[2]= Q

		
	DISCRETE
		
		ASSERT( calculateDerivAkima(Q, Qp_num, TIME, xsaved, timePrev, intMet) ) NOTE ""

		
		ASSERT( m > 1e-3 ) KILLPOINT "The deposit is empty"
	
		
		
	CONTINUOUS
		Init_fluid(fluid, g.fluid)	
		
		IF (Type == Known_Qp) INSERT
		Qp_data = Qp
		Q' =  Qp
		END IF 
		
		IF (Type == Known_k) INSERT
		k_data = k
		Qp = Qp_num
		END IF
		
		m' = -g.W
		U= Cv * Tt * m
		U- U_0 = Q - flux
		flux' = g.W * Cp * Tt
		g.pt = pt
		g.Tt = Tt
		pt * Vol_d = m * Rgas * Tt
		(m / m_0)**k = pt / p_0

END COMPONENT


--------------------------------------------------------------------------------
-- Component that represents a gas or liquid valve 
--------------------------------------------------------------------------------


COMPONENT Valve IS_A FluidInFluidOut (ENUM ChemState Type = Gas)
	DATA
		

	DECLS 
		ALG REAL M_th						UNITS no_units						"Mach number in gas valve throat"
		REAL rho_f							UNITS	u_kg_m3						"Fluid density"
		BOUND REAL A =0.001						UNITS u_m2							"Area in gas or liquid valve"
		
		REAL dp 								UNITS	u_Pa							"Pressure drop"
		
	INIT
		
		
		
	DISCRETE
		ASSERT(M_th <=1 ) KILLPOINT "THE VALVE IS CHOKED!"
		ASSERT((Type == Gas AND State(f_in.fluid) == Gas) OR ((Type == Liquid AND State(f_in.fluid) == Liquid))) KILLPOINT "CHECK VALVE TYPE! SELECTED TYPE DOES NOT MATCH ACTUAL FLUID STATE THROUGH VALVE"
	CONTINUOUS 
	f_in.Tt = f_out.Tt
	
	IF (Type == Gas) INSERT
		f_in.W * sqrt(R(f_in.fluid) * f_in.Tt) / A / f_in.pt = sqrt(gamma(f_in.fluid)) * M_th * (1. + (gamma(f_in.fluid) - 1.) / 2. * M_th**2)**(- (gamma(f_in.fluid) + 1.) / 2. / (gamma(f_in.fluid) - 1.))
		((f_in.pt / R(f_in.fluid) / f_in.Tt) / rho_f)**(gamma(f_in.fluid) - 1.) = 1. + (gamma(f_in.fluid) - 1.) / 2. * M_th**2
		(f_in.pt / f_out.pt)**((gamma(f_in.fluid) - 1.) / gamma(f_in.fluid)) = 1. + (gamma(f_in.fluid) - 1.) / 2. * M_th**2	
	END IF
		
	IF (Type == Liquid) INSERT
		M_th = 0
		rho_f = rho(f_in.fluid)
		f_out.pt = f_in.pt - 1/2 * f_in.W**2 / ( rho_f * A**2 ) 
	END IF
		
	dp = -f_out.pt + f_in.pt
END COMPONENT

	
COMPONENT New_Junction (ENUM Type_Junction Type = Incompressible_l)

   "Junction with pressure drop"

	PORTS
      IN Fluid f_in1             		    				"Inlet fluid port"
		IN Fluid f_in2                						"Inlet fluid port"
      OUT Fluid f_out               						"Outlet fluid port"
		
	DATA	
		REAL A_ratio1 = 0.2			UNITS no_units			"Area ratio A_in1 / A_out"
		REAL TPL = 0.98				UNITS no_units			"Total pressure load"
		REAL PR = 1 					UNITS no_units 		"Ratio between inlets total pressures"
		REAL M_A1 = 0.3 				UNITS no_units			"Known mach inlet 1"
		
	DECLS
		REAL M_in1						UNITS no_units			"Mach number in inlet 1"
		ALG REAL M_in2	= 0.3			UNITS no_units			"Mach number in inlet 2"
		ALG REAL M_out = 0.3			UNITS no_units			"Mach number in outlet"
		REAL A_in1						UNITS u_m2				"Area inlet 1"
		REAL A_in2						UNITS u_m2				"Area inlet 2"
		REAL A_out 						UNITS u_m2				"Area outlet 3"
		REAL p_in1						UNITS u_Pa
		REAL p_in2						UNITS u_Pa
		REAL p_out						UNITS u_Pa
		REAL T_in1						UNITS u_K
		REAL T_in2						UNITS u_K
		REAL T_out						UNITS u_K
		REAL V_in1						UNITS u_m_s
		REAL V_in2						UNITS u_m_s
		REAL V_out						UNITS u_m_s
		REAL A_ratio2					UNITS no_units
		
		REAL pt_out 					UNITS no_units
		
   CONTINUOUS	
		f_in1.W + f_in2.W = f_out.W
		
		IF (Type == Compressible ) INSERT
			A_ratio2 = 0
			M_A1 = M_in1
			A_in1 + A_in2 = A_out
		
			V_in1 = M_in1 * sqrt(gamma(f_in1.fluid)* T_in1 * R(f_in1.fluid))
			V_in2 = M_in2 * sqrt(gamma(f_in2.fluid)* T_in2 * R(f_in2.fluid))
			V_out = M_out * sqrt(gamma(f_out.fluid)* T_out * R(f_out.fluid))
			
			f_in1.W * sqrt(R(f_in1.fluid) * f_in1.Tt) / A_in1 / f_in1.pt = sqrt(gamma(f_in1.fluid)) * M_in1 * (1. + (gamma(f_in1.fluid) - 1.) / 2. * M_in1**2)**(- (gamma(f_in1.fluid) + 1.) / 2. / (gamma(f_in1.fluid) - 1.))
			f_in2.W * sqrt(R(f_in2.fluid) * f_in2.Tt) / A_in2 / f_in2.pt = sqrt(gamma(f_in2.fluid)) * M_in2 * (1. + (gamma(f_in2.fluid) - 1.) / 2. * M_in2**2)**(- (gamma(f_in2.fluid) + 1.) / 2. / (gamma(f_in2.fluid) - 1.))
			f_out.W * sqrt(R(f_out.fluid) * f_out.Tt) / A_out / pt_out = sqrt(gamma(f_out.fluid)) * M_out * (1. + (gamma(f_out.fluid) - 1.) / 2. * M_out**2)**(- (gamma(f_out.fluid) + 1.) / 2. / (gamma(f_out.fluid) - 1.))
			
			f_out.W * V_out + p_out * A_out = f_in1.W * V_in1 + p_in1 * A_in1 + f_in2.W * V_in2 + p_in2 * A_in2
			
			p_in1 = p_in2
			
			
			f_in1.pt/p_in1 =  (1 + (gamma(f_in1.fluid) - 1)/2 * M_in1**2 ) ** ( gamma(f_in1.fluid)/(gamma(f_in1.fluid) - 1) )
			f_in2.pt/p_in2 =  (1 + (gamma(f_in2.fluid) - 1)/2 * M_in2**2 ) ** ( gamma(f_in2.fluid)/(gamma(f_in2.fluid) - 1) )
			pt_out/p_out =  (1 + (gamma(f_out.fluid) - 1)/2 * M_out**2 ) ** ( gamma(f_out.fluid)/(gamma(f_out.fluid) - 1) )
			
			f_in1.Tt/T_in1 = ( f_in1.pt / p_in1 ) ** ( (gamma(f_in1.fluid) - 1) / gamma(f_in1.fluid) )
			f_in2.Tt/T_in2 = ( f_in2.pt / p_in2 ) ** ( (gamma(f_in2.fluid) - 1) / gamma(f_in2.fluid) )
			f_out.Tt/T_out = ( pt_out / p_out ) ** ( (gamma(f_out.fluid) - 1) / gamma(f_out.fluid) )
			
			f_in1.W * cp(f_in1.fluid) * f_in1.Tt + f_in2.W * cp(f_in1.fluid) * f_in2.Tt = f_out.W * cp(f_in1.fluid) * f_out.Tt
			pt_out = f_out.pt
			
		END IF
		
		IF (Type == Incompressible_l) INSERT
		
			M_in1 = 0
			M_in2 = 0
			M_out = 0
			p_in1 = 0
			p_in2 = 0
			p_out = 0
			T_in1 = 0
			T_in2 = 0
			T_out = 0
			V_in1 = 0
			V_in2 = 0
			V_out = 0
			
			A_ratio1 = A_in1 / A_out
			A_ratio2 = A_in2 / A_out
			PR =  f_in1.pt / f_in2.pt
			(f_in1.W / f_out.W) * cp(f_in1.fluid) * (f_out.Tt - f_in1.Tt) + (f_in2.W / f_out.W) * cp(f_in2.fluid) * (f_out.Tt - f_in2.Tt) = 0
			f_out.pt / TPL = pt_out 
			f_in1.pt * A_ratio1 + f_in2.pt * A_ratio2 = pt_out
			A_ratio2 = 1 - A_ratio1
			
			
		
		END IF
		
		EXPAND (i IN LiquidsGases) f_out.fluid[i] = (f_in1.fluid[i]*f_in1.W + f_in2.fluid[i]*f_in2.W)/f_out.W
		f_out.fluid[Comb_prod] = (f_in1.fluid[Comb_prod]*f_in1.W + f_in2.fluid[Comb_prod]*f_in2.W)/f_out.W
		f_out.fluid[Comb_prod_M] = f_in1.fluid[Comb_prod_M]
		f_out.fluid[Comb_prod_cp] = f_in1.fluid[Comb_prod_cp]
	
END COMPONENT	


COMPONENT Junction (ENUM Type_Junction Type = Compressible)

   "Junction with pressure drop"

	PORTS
      IN Fluid f_in1             		    				"Inlet fluid port"
		IN Fluid f_in2                						"Inlet fluid port"
      OUT Fluid f_out               						"Outlet fluid port"
		
	DATA	
		REAL TPL = 0.98					UNITS no_units			"Total pressure load"
		REAL A_in1 = 0.001				UNITS u_m2				"Area inlet 1"
		REAL A_in2 = 0.001				UNITS u_m2				"Area inlet 2"
		
	DECLS
		REAL M_in1						UNITS no_units			RANGE 0, 1"Mach number in inlet 1"
		ALG REAL M_in2	= 0.3			UNITS no_units			RANGE 0, 1"Mach number in inlet 2"
		ALG REAL M_out = 0.3			UNITS no_units			RANGE 0, 1 "Mach number in outlet"
		REAL A_out 						UNITS u_m2				"Area outlet 3"
		REAL p_in1						UNITS u_Pa				"Static pressure in inlet 1"
		REAL p_in2						UNITS u_Pa				"Static pressure in inlet 2"
		REAL p_out						UNITS u_Pa				"Static pressure in outlet"
		REAL T_in1						UNITS u_K				"Static temperature in inlet 1"
		REAL T_in2						UNITS u_K				"Static temperature in inlet 2"
		REAL T_out						UNITS u_K				"Static temperature in outlet"
		REAL V_in1						UNITS u_m_s				"Injection speed in inlet 1"
		REAL V_in2						UNITS u_m_s				"Injection speed in inlet 2"
		REAL V_out						UNITS u_m_s				"Outlet speed"
		
		REAL pt_out 					UNITS no_units       "Total outlet pressure"
		
		
		
	INIT
		f_out.fluid[Comb_prod]=0
		f_in1.fluid[Comb_prod]=0
		f_in2.fluid[Comb_prod]=0
		f_out.fluid[Comb_prod_M]=0
		f_in1.fluid[Comb_prod_M]=0
		f_in2.fluid[Comb_prod_M]=0
		f_out.fluid[Comb_prod_cp]=0
		f_in1.fluid[Comb_prod_cp]=0
		f_in2.fluid[Comb_prod_cp]=0
		f_out.fluid[Comb_prod_c]=0
		f_in1.fluid[Comb_prod_c]=0
		f_in2.fluid[Comb_prod_c]=0
		f_out.fluid[Comb_prod_cp_g]=0
		f_in1.fluid[Comb_prod_cp_g]=0
		f_in2.fluid[Comb_prod_cp_g]=0
		f_out.fluid[Comb_prod_Lv]=0
		f_in1.fluid[Comb_prod_Lv]=0
		f_in2.fluid[Comb_prod_Lv]=0
		f_out.fluid[Comb_prod_cond]=0
		f_in1.fluid[Comb_prod_cond]=0
		f_in2.fluid[Comb_prod_cond]=0
		f_out.fluid[Comb_prod_visc]=0
		f_in1.fluid[Comb_prod_visc]=0
		f_in2.fluid[Comb_prod_visc]=0
		f_out.fluid[Comb_prod_T_boil]=0
		f_in1.fluid[Comb_prod_T_boil]=0
		f_in2.fluid[Comb_prod_T_boil]=0
		f_out.fluid[Comb_prod_rho]=0
		f_in1.fluid[Comb_prod_rho]=0
		f_in2.fluid[Comb_prod_rho]=0
		
		
	DISCRETE

		ASSERT(M_out < 1) WARNING "MACH NUMBER IN OUTLET SECTION CONVERGED IN SUPERSONIC SOLUTION"
		ASSERT(M_in1 < 1) WARNING "MACH NUMBER IN INLET 1 SECTION CONVERGED IN SUPERSONIC SOLUTION"
		ASSERT(M_in2 < 1) WARNING "MACH NUMBER IN INLET 2 SECTION CONVERGED IN SUPERSONIC SOLUTION"
		
		ASSERT( (State(f_in1.fluid) == Gas AND State(f_in2.fluid) ==  Gas) OR (State(f_in1.fluid) == Liquid AND State(f_in2.fluid) ==  Liquid) ) WARNING "ONLY FLUIDS IN THE SAME STATE CAN BE MIXED" 

		
   CONTINUOUS	
		f_in1.W + f_in2.W = f_out.W
		A_in1 + A_in2 = A_out
		
		EXPAND (i IN LiquidsGases) f_out.fluid[i] = (f_in1.fluid[i]*f_in1.W + f_in2.fluid[i]*f_in2.W)/f_out.W
		f_out.fluid[Comb_prod] = (f_in1.fluid[Comb_prod]*f_in1.W + f_in2.fluid[Comb_prod]*f_in2.W)/f_out.W
		
		f_out.fluid[Comb_prod_M] = IF (f_in1.fluid[Comb_prod]<1E-4 AND f_in2.fluid[Comb_prod]<1E-4) 0
											ELSE (f_in1.fluid[Comb_prod]*f_in1.W + f_in2.fluid[Comb_prod]*f_in2.W)/(f_in1.fluid[Comb_prod]*f_in1.W/f_in1.fluid[Comb_prod_M] + f_in2.fluid[Comb_prod]*f_in2.W/f_in2.fluid[Comb_prod_M])
		f_out.fluid[Comb_prod_cp] = IF (f_in1.fluid[Comb_prod]<1E-4 AND f_in2.fluid[Comb_prod]<1E-4) 0
									       ELSE (f_in1.fluid[Comb_prod]*f_in1.W*f_in1.fluid[Comb_prod_cp] + f_in2.fluid[Comb_prod]*f_in2.W*f_in2.fluid[Comb_prod_cp])/(f_in1.fluid[Comb_prod]*f_in1.W + f_in2.fluid[Comb_prod]*f_in2.W)
		f_out.fluid[Comb_prod_c] = IF (f_in1.fluid[Comb_prod]<1E-4 AND f_in2.fluid[Comb_prod]<1E-4) 0
									       ELSE (f_in1.fluid[Comb_prod]*f_in1.W*f_in1.fluid[Comb_prod_c] + f_in2.fluid[Comb_prod]*f_in2.W*f_in2.fluid[Comb_prod_c])/(f_in1.fluid[Comb_prod]*f_in1.W + f_in2.fluid[Comb_prod]*f_in2.W)
		f_out.fluid[Comb_prod_cp_g] = IF (f_in1.fluid[Comb_prod]<1E-4 AND f_in2.fluid[Comb_prod]<1E-4) 0
									       ELSE (f_in1.fluid[Comb_prod]*f_in1.W*f_in1.fluid[Comb_prod_cp_g] + f_in2.fluid[Comb_prod]*f_in2.W*f_in2.fluid[Comb_prod_cp_g])/(f_in1.fluid[Comb_prod]*f_in1.W + f_in2.fluid[Comb_prod]*f_in2.W)										 								 
		f_out.fluid[Comb_prod_Lv] = IF (f_in1.fluid[Comb_prod]<1E-4 AND f_in2.fluid[Comb_prod]<1E-4) 0
									       ELSE (f_in1.fluid[Comb_prod]*f_in1.W*f_in1.fluid[Comb_prod_Lv] + f_in2.fluid[Comb_prod]*f_in2.W*f_in2.fluid[Comb_prod_Lv])/(f_in1.fluid[Comb_prod]*f_in1.W + f_in2.fluid[Comb_prod]*f_in2.W)	
		f_out.fluid[Comb_prod_visc] = IF (f_in1.fluid[Comb_prod]<1E-4 AND f_in2.fluid[Comb_prod]<1E-4) 0
									       ELSE (f_in1.fluid[Comb_prod]*f_in1.W*f_in1.fluid[Comb_prod_visc] + f_in2.fluid[Comb_prod]*f_in2.W*f_in2.fluid[Comb_prod_visc])/(f_in1.fluid[Comb_prod]*f_in1.W + f_in2.fluid[Comb_prod]*f_in2.W)			
		f_out.fluid[Comb_prod_cond] = IF (f_in1.fluid[Comb_prod]<1E-4 AND f_in2.fluid[Comb_prod]<1E-4) 0
									       ELSE (f_in1.fluid[Comb_prod]*f_in1.W*f_in1.fluid[Comb_prod_cond] + f_in2.fluid[Comb_prod]*f_in2.W*f_in2.fluid[Comb_prod_cond])/(f_in1.fluid[Comb_prod]*f_in1.W + f_in2.fluid[Comb_prod]*f_in2.W)
		f_out.fluid[Comb_prod_T_boil] = IF (f_in1.fluid[Comb_prod]<1E-4 AND f_in2.fluid[Comb_prod]<1E-4) 0
									       ELSE (f_in1.fluid[Comb_prod]*f_in1.W*f_in1.fluid[Comb_prod_T_boil] + f_in2.fluid[Comb_prod]*f_in2.W*f_in2.fluid[Comb_prod_T_boil])/(f_in1.fluid[Comb_prod]*f_in1.W + f_in2.fluid[Comb_prod]*f_in2.W)									 
		f_out.fluid[Comb_prod_rho] = IF (f_in1.fluid[Comb_prod]<1E-4 AND f_in2.fluid[Comb_prod]<1E-4) 0
									       ELSE (f_in1.fluid[Comb_prod]*f_in1.W*f_in1.fluid[Comb_prod_rho] + f_in2.fluid[Comb_prod]*f_in2.W*f_in2.fluid[Comb_prod_rho])/(f_in1.fluid[Comb_prod]*f_in1.W + f_in2.fluid[Comb_prod]*f_in2.W)																	
											
		IF (Type == Compressible ) INSERT
		
			V_in1 = M_in1 * sqrt(gamma(f_in1.fluid)* T_in1 * R(f_in1.fluid))
			V_in2 = M_in2 * sqrt(gamma(f_in2.fluid)* T_in2 * R(f_in2.fluid))
			V_out = M_out * sqrt(gamma(f_out.fluid)* T_out * R(f_out.fluid))
			
			f_in1.W * sqrt(R(f_in1.fluid) * f_in1.Tt) / A_in1 / f_in1.pt = sqrt(gamma(f_in1.fluid)) * M_in1 * (1. + (gamma(f_in1.fluid) - 1.) / 2. * M_in1**2)**(- (gamma(f_in1.fluid) + 1.) / 2. / (gamma(f_in1.fluid) - 1.))
			f_in2.W * sqrt(R(f_in2.fluid) * f_in2.Tt) / A_in2 / f_in2.pt = sqrt(gamma(f_in2.fluid)) * M_in2 * (1. + (gamma(f_in2.fluid) - 1.) / 2. * M_in2**2)**(- (gamma(f_in2.fluid) + 1.) / 2. / (gamma(f_in2.fluid) - 1.))
			f_out.W * sqrt(R(f_out.fluid) * f_out.Tt) / A_out / pt_out = sqrt(gamma(f_out.fluid)) * M_out * (1. + (gamma(f_out.fluid) - 1.) / 2. * M_out**2)**(- (gamma(f_out.fluid) + 1.) / 2. / (gamma(f_out.fluid) - 1.))
			
			f_out.W * V_out + p_out * A_out = f_in1.W * V_in1 + p_in1 * A_in1 + f_in2.W * V_in2 + p_in2 * A_in2
			
			p_in1 = p_in2
			
			
			f_in1.pt/p_in1 =  (1 + (gamma(f_in1.fluid) - 1)/2 * M_in1**2 ) ** ( gamma(f_in1.fluid)/(gamma(f_in1.fluid) - 1) )
			f_in2.pt/p_in2 =  (1 + (gamma(f_in2.fluid) - 1)/2 * M_in2**2 ) ** ( gamma(f_in2.fluid)/(gamma(f_in2.fluid) - 1) )
			pt_out/p_out =  (1 + (gamma(f_out.fluid) - 1)/2 * M_out**2 ) ** ( gamma(f_out.fluid)/(gamma(f_out.fluid) - 1) )
			
			f_in1.Tt/T_in1 = ( f_in1.pt / p_in1 ) ** ( (gamma(f_in1.fluid) - 1) / gamma(f_in1.fluid) )
			f_in2.Tt/T_in2 = ( f_in2.pt / p_in2 ) ** ( (gamma(f_in2.fluid) - 1) / gamma(f_in2.fluid) )
			f_out.Tt/T_out = ( pt_out / p_out ) ** ( (gamma(f_out.fluid) - 1) / gamma(f_out.fluid) )
			
			f_in1.W * cp(f_in1.fluid) * f_in1.Tt + f_in2.W * cp(f_in2.fluid) * f_in2.Tt = f_out.W * cp(f_out.fluid) * f_out.Tt
			pt_out = f_out.pt
			
		END IF
		
		IF (Type == Incompressible_g) INSERT
		
			M_in1 = 0
			M_in2 = 0
			M_out = 0
			
			T_in1 = 0
			T_in2 = 0
			T_out = 0
			
			f_in1.W = V_in1 * f_in1.pt / R(f_in1.fluid) / f_in1.Tt * A_in1
			f_in2.W = V_in2 * f_in2.pt / R(f_in2.fluid) / f_in2.Tt * A_in2
			f_out.W = V_out * f_out.pt / R(f_out.fluid) / f_out.Tt * A_out
			
			f_in1.pt = p_in1 + 0.5 * f_in1.pt / ( R(f_in1.fluid) * f_in1.Tt) * V_in1**2
			f_in2.pt = p_in2 + 0.5 * f_in2.pt / ( R(f_in2.fluid) * f_in2.Tt) * V_in2**2
			f_out.pt = p_out + 0.5 * f_out.pt / ( R(f_out.fluid) * f_out.Tt) * V_out**2
			
			p_in1 = p_in2
			
			f_out.W * V_out + p_out * A_out = f_in1.W * V_in1 + p_in1 * A_in1 + f_in2.W * V_in2 + p_in2 * A_in2
			
			f_in1.W * cp(f_in1.fluid) * f_in1.Tt + f_in2.W * cp(f_in2.fluid) * f_in2.Tt = f_out.W * cp(f_out.fluid) * f_out.Tt
			
			pt_out * TPL = f_out.pt
		
		END IF
		
		IF (Type == Incompressible_l) INSERT
		
			M_in1 = 0
			M_in2 = 0
			M_out = 0
			
			T_in1 = 0
			T_in2 = 0
			T_out = 0
			
			f_in1.W = V_in1 * rho(f_in1.fluid) * A_in1
			f_in2.W = V_in2 * rho(f_in2.fluid) * A_in2
			f_out.W = V_out * rho(f_out.fluid) * A_out
			
			f_in1.pt = p_in1 + 0.5 * rho(f_in1.fluid) * V_in1**2
			f_in2.pt = p_in2 + 0.5 * rho(f_in2.fluid) * V_in2**2
			f_out.pt = p_out + 0.5 * rho(f_out.fluid) * V_out**2
			
			p_in1 = p_in2
			
			f_out.W * V_out + p_out * A_out = f_in1.W * V_in1 + p_in1 * A_in1 + f_in2.W * V_in2 + p_in2 * A_in2
			
			f_in1.W * cp(f_in1.fluid) * f_in1.Tt + f_in2.W * cp(f_in2.fluid) * f_in2.Tt = f_out.W * cp(f_out.fluid) * f_out.Tt
			
			pt_out * TPL = f_out.pt
		
		END IF
		
		
		
	
END COMPONENT




COMPONENT Tank_real (ENUM Type_Deposit Type = Known_k)

	"Liquid tank"
	
	PORTS
		IN Fluid g						"Inlet fluid port"
		OUT Fluid l						"Outlet fluid port"
	DATA		
		ENUM Liquids fluid_l	= LOX									"Working liquid name"
		REAL T_d = 288.15					UNITS u_K				"Tank temperature"		
		REAL Ag = 0.001					UNITS u_m2				"Gas input area"
		REAL T_0 = 298.15					UNITS u_K 				"Gas initial temperature"
		REAL p_0 = 2900000				UNITS u_Pa
		REAL Vol_t	= 100					UNITS u_m3				"Tank volume"
		REAL IVR = 0.01					UNITS no_units 		"Initial volume ratio Gas volume / Tank volume"
		REAL Qp_data = 0
		REAL k_data = 1
		ENUM t_interp intMet= AKIMA						"Interpolation method"
		
	DECLS		
		DISCR REAL Rgas					UNITS u_J_kgK 
		DISCR REAL Cv						UNITS u_J_kgK 
		DISCR REAL Cp						UNITS u_J_kgK
		REAL p_g								UNITS u_Pa			"Pressurisation gas pressure"
		REAL d 								UNITS u_kg_m3
		REAL T_g								UNITS u_K			"Pressurisation gas temperature"
		REAL U 								UNITS u_J			"Pressurisation gas internal energy"
		DISCR REAL U_0						UNITS u_J			"Initial ressurisation gas internal energy"
		DISCR REAL d_0						UNITS u_kg_m3		"Initial ressurisation gas internal energy"
		REAL Vol_g							UNITS u_m3			"Pressurisation gas volume"
		REAL Qp 								UNITS u_W			"Heat per second added to the pressurisation gas"
		REAL Q								UNITS u_J			"Total heat added to pressurisation gas"
		ALG REAL M_g = 0.1				UNITS no_units		"Pressurisation gas Mach number"
		REAL p_d								UNITS u_Pa			"Design tank pressure"
		REAL m	 							UNITS u_kg			"Pressurisation gas mass"
		REAL flux 							UNITS u_J
		REAL W								UNITS u_J
		
		
		REAL k								UNITS no_units
		DISCR REAL Qp_num, xsaved[4], timePrev[4] = 0
		
	INIT PRIORITY 100
	
		Init_fluid(fluid_l, l.fluid)
		Rgas = R(g.fluid)
		Cp = cp(g.fluid)
		Cv = cv(g.fluid)
		--Cp = 5193
		--Cv = 3115.85015
		--Rgas = 2077.14985
		Vol_g = IVR * Vol_t
		m = p_0 * Vol_g /(Rgas * T_0)
		U = Cv * m * T_0
		U_0 = Cv * m * T_0
		d_0 = m / Vol_g
		p_g = p_0
		T_g = T_0
		xsaved[1]= 0
      xsaved[2]= Q			
		
	DISCRETE
		ASSERT (Gas == State(g.fluid)) FATAL "ONLY GASES CAN PRESSURISE THE TANK!"
		ASSERT (M_g <= 1.) KILLPOINT "THE PRESSURISATION GAS INLET CANNOT BE CHOKED!"
		ASSERT(calculateDerivAkima(Q, Qp_num, TIME, xsaved, timePrev, intMet) ) NOTE ""
		ASSERT(Vol_t > Vol_g) KILLPOINT "THE TANK IS EMPTY"	
	CONTINUOUS
	
		-- Inlet gas
		Init_fluid(fluid_l, l.fluid)
		l.Tt = T_d
		l.pt = p_g
		
		g.W * sqrt(R(g.fluid) * g.Tt) / Ag / g.pt = sqrt(gamma(g.fluid)) * M_g * (1. + (gamma(g.fluid) - 1.) / 2. * M_g**2)**(- (gamma(g.fluid) + 1.) / 2. / (gamma(g.fluid) - 1.))
		(g.pt / p_g)**((gamma(g.fluid) - 1.) / gamma(g.fluid)) = 1. + (gamma(g.fluid) - 1.) / 2. * M_g**2	
		
		IF (Type == Known_Qp) INSERT
		Qp_data = Qp
		Q' =  Qp
		END IF 
		
		IF (Type == Known_k) INSERT
		k_data = k
		Qp = Qp_num
		END IF
		
		p_g = p_d
		
		-- Expansion process
		--Q' =  Qp
		p_g = m * Rgas * T_g/Vol_g
		--T_g' = (Qp - p_g * l.W / rho(l.fluid) + g.W * Cp * g.Tt - g.W * Cv * T_g)/ m / Cv 
		U-U_0 = W + Q + flux
		--p_g' = (T_g' / T_g + g.W / m - l.W / rho(l.fluid)/Vol_g)*p_g
		m' = g.W
		U = m * Cv * T_g
		(p_g / p_0) = (d /d_0)**k
		Vol_g' = l.W / rho(l.fluid)
		W' = -p_g * l.W / rho(l.fluid)
		flux' = g.W * Cp * g.Tt
		d = p_g / (Rgas * T_g)

		
END COMPONENT



COMPONENT TankBD_ideal

	"Liquid tank for blowdown process"
	
	PORTS
		IN Fluid g						"Inlet fluid port"
		OUT Fluid l						"Outlet fluid port"
	DATA		
		ENUM Liquids fluid_l	= LOX									"Working liquid name"
		REAL T_d = 288.15					UNITS u_K				"Tank temperature"		
		REAL T_0 = 298.15					UNITS u_K 				"Gas initial temperature"
		REAL p_0 = 2900000				UNITS u_Pa
		REAL Vol_t	= 100					UNITS u_m3				"Tank volume"
		REAL IVR = 0.01					UNITS no_units 		"Initial volume ratio Gas volume / Tank volume"
		REAL Qp = 0							UNITS u_W
		
	DECLS		
		DISCR REAL Rgas					UNITS u_J_kgK 
		DISCR REAL Cv						UNITS u_J_kgK 
		DISCR REAL Cp						UNITS u_J_kgK
		REAL p_g								UNITS u_Pa			"Pressurisation gas pressure"
		REAL d 								UNITS u_kg_m3
		REAL T_g								UNITS u_K			"Pressurisation gas temperature"
		DISCR REAL d_0						UNITS u_kg_m3			"Initial ressurisation gas internal energy"
		REAL Vol_g							UNITS u_m3			"Pressurisation gas volume"
		REAL Q								UNITS u_J			"Total heat added to pressurisation gas"
		REAL p_d								UNITS u_Pa			"Design tank pressure"
		REAL m	 							UNITS u_kg			"Pressurisation gas mass"
		
		
		REAL k								UNITS no_units
		DISCR REAL Qp_num, xsaved[4], timePrev[4] = 0
		
	INIT PRIORITY 100
	
		Init_fluid(fluid_l, l.fluid)
		Rgas = R(g.fluid)
		Cp = cp(g.fluid)
		Cv = cv(g.fluid)
		--Cp = 5193
		--Cv = 3115.85015
		--Rgas = 2077.14985
		Vol_g = IVR * Vol_t
		m = p_0 * Vol_g /(Rgas * T_0)
		
		d_0 = m / Vol_g
		p_g = p_0
		T_g = T_0
	DISCRETE
		ASSERT (Gas == State(g.fluid)) FATAL "ONLY GASES CAN PRESSURISE THE TANK!"
		
	CONTINUOUS
	
		-- Inlet gas
		Init_fluid(fluid_l, l.fluid)
		l.Tt = T_d
		l.pt = p_g
		p_g = g.pt
		p_g = p_d
		
		-- Expansion process
		Q' = Qp
		p_g := m * Rgas * T_g/Vol_g
		T_g' = (Qp - p_g * l.W / rho(l.fluid) + g.W * Cp * g.Tt - g.W * Cv * T_g)/ m / Cv 
		p_g' = (T_g' / T_g + g.W / m - l.W / rho(l.fluid)/Vol_g)*p_g
		(p_g / p_0) = (d /d_0)**k
		Vol_g' = l.W / rho(l.fluid)
		d = p_g / (Rgas * T_g)

		
END COMPONENT



COMPONENT TankPF_ideal

	"Liquid tank for pressure fed process"
	
	PORTS
		IN Fluid g						"Inlet fluid port"
		OUT Fluid l						"Outlet fluid port"
	DATA		
		ENUM Liquids fluid_l	= LOX									"Working liquid name"
		REAL T_d = 288.15					UNITS u_K				"Tank temperature"		
		REAL T_0 = 298.15					UNITS u_K 				"Gas initial temperature"
		REAL p_0 = 2900000				UNITS u_Pa				"Gas initial pressure. Process pressure"
		REAL Vol_t	= 100					UNITS u_m3				"Tank volume"
		REAL IVR = 0.01					UNITS no_units 		"Initial volume ratio Gas volume / Tank volume"
		REAL Qp = 0							UNITS u_W
		
	DECLS		
		DISCR REAL Rgas					UNITS u_J_kgK 
		DISCR REAL Cv						UNITS u_J_kgK 
		DISCR REAL Cp						UNITS u_J_kgK
		REAL p_g								UNITS u_Pa			"Pressurisation gas pressure"
		REAL d 								UNITS u_kg_m3		"Pressurisation gas density"
		REAL T_g								UNITS u_K			"Pressurisation gas temperature"
		DISCR REAL d_0						UNITS u_kg_m3		"Initial ressurisation gas internal energy"
		REAL Vol_g							UNITS u_m3			"Pressurisation gas volume"
		REAL Q								UNITS u_J			"Total heat added to pressurisation gas"
		ALG REAL M_g = 0.1				UNITS no_units		"Pressurisation gas Mach number"
		REAL Ag 								UNITS u_m2			"Regulator area"
		REAL p_d								UNITS u_Pa			"Design tank pressure"
		REAL m	 							UNITS u_kg			"Pressurisation gas mass"
		REAL k								UNITS no_units		
		
	INIT PRIORITY 100
	
		Init_fluid(fluid_l, l.fluid)
		Rgas = R(g.fluid)
		Cp = cp(g.fluid)
		Cv = cv(g.fluid)
		--Cp = 5193
		--Cv = 3115.85015
		--Rgas = 2077.14985
		Vol_g = IVR * Vol_t
		m = p_0 * Vol_g /(Rgas * T_0)
		
		d_0 = m / Vol_g
		p_g = p_0
		T_g = T_0
	DISCRETE
		ASSERT (Gas == State(g.fluid)) FATAL "ONLY GASES CAN PRESSURISE THE TANK!"
		ASSERT (M_g <= 1.) KILLPOINT "REGULATOR GAS INLET CANNOT BE  CHOKED"
		
	CONTINUOUS
	
		-- Inlet gas
		Init_fluid(fluid_l, l.fluid)
		l.Tt = T_d
		l.pt = p_g
		p_g = p_d
	
		
		-- Choking process in the regulator
		g.W * sqrt(R(g.fluid) * g.Tt) / Ag / g.pt = sqrt(gamma(g.fluid)) * M_g * (1. + (gamma(g.fluid) - 1.) / 2. * M_g**2)**(- (gamma(g.fluid) + 1.) / 2. / (gamma(g.fluid) - 1.))
		(g.pt / p_g) = (1. + (gamma(g.fluid) - 1.) / 2. * M_g**2)**(gamma(g.fluid) / (gamma(g.fluid)-1))
		
		-- Expansion process
		Q' = Qp
		p_g := m * Rgas * T_g/Vol_g
		T_g' = (Qp - p_g * l.W / rho(l.fluid) + g.W * Cp * g.Tt - g.W * Cv * T_g)/ m / Cv 
		p_g' = (T_g' / T_g + g.W / m - l.W / rho(l.fluid)/Vol_g)*p_g
		(p_g / p_0) = (d /d_0)**k
		Vol_g' = l.W / rho(l.fluid)
		d = p_g / (Rgas * T_g)
		g.W = (l.W / rho(l.fluid)/Vol_g*T_g/g.Tt - Qp/(Cp*m*g.Tt))*m
		
END COMPONENT


COMPONENT DeadStart 

	"Conditions of a fluid inlet in the system"
	
	PORTS
		OUT Fluid f						"Outlet fluid port"
	
	DATA
		ENUM LiquidsGases fluid = LOX											"Working fluid name"

	INIT PRIORITY 100
		Init_fluid(fluid, f.fluid)

	CONTINUOUS
		Init_fluid(fluid, f.fluid)
		f.W = 0
		

END COMPONENT

COMPONENT DeadEnd 

	"Conditions of a fluid inlet in the system"
	
	PORTS
		OUT Fluid f						"Outlet fluid port"
	
	DATA

	INIT PRIORITY 100

	CONTINUOUS
		f.W = 0
		

END COMPONENT

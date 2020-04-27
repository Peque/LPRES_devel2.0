--------------------------------------------------------------------------------
-- EcosimPro Simulation Source Code
-- FILE NAME: LPRES_Ports.el
-- DESCRIPTION: Defines common port types
-- NOTES:
-- AUTHOR: Pablo Sierra Heras
-- CREATION DATE: 05/12/2014
--------------------------------------------------------------------------------

-- Libraries
USE MATH VERSION "3.1"


-- Ports
--------------------------------------------------------------------------------
-- Port for fluid flow connections
--------------------------------------------------------------------------------
PORT Fluid SINGLE IN

	"Fluid port"

	EQUAL 			REAL Tt						UNITS u_K			RANGE ZERO,Inf		"Total temperature" 
	EQUAL 			REAL pt						UNITS u_Pa			RANGE ZERO,Inf		"Total pressure"
	SUM				REAL W						UNITS u_kg_s		RANGE ZERO,Inf		"Mass flow"
	HIDDEN EQUAL	REAL fluid[ChemName]	UNITS no_units		RANGE ZERO,Inf		"Working fluid"
	
END PORT

--------------------------------------------------------------------------------
-- Port for injection in a combustion chamber
--------------------------------------------------------------------------------
PORT FluidInj SINGLE

	"Fluid injection port"

	EQUAL 			REAL T						UNITS u_K			RANGE ZERO,Inf		"Temperature" 
	EQUAL 			REAL p						UNITS u_Pa			RANGE ZERO,Inf		"Pressure"
	SUM				REAL W						UNITS u_kg_s		RANGE ZERO,Inf		"Mass flow"
	HIDDEN EQUAL 	REAL fluid[ChemName]		UNITS no_units		RANGE ZERO,Inf		"Working fluid"
	EQUAL 			REAL p_c						UNITS u_Pa			RANGE ZERO,Inf		"Chamber pressure"
	
END PORT

--------------------------------------------------------------------------------
-- Port for gas flow connections through a nozzle
--------------------------------------------------------------------------------
PORT GasNozzle SINGLE

	 "Gas port through a nozzle"

	EQUAL 			REAL Tt						UNITS u_K			RANGE ZERO,Inf		"Total temperature" 
	EQUAL 			REAL pt						UNITS u_Pa			RANGE ZERO,Inf		"Total pressure"
	SUM				REAL W						UNITS u_kg_s		RANGE ZERO,Inf		"Mass flow"
	HIDDEN EQUAL	REAL fluid[ChemName]		UNITS no_units		RANGE ZERO,Inf		"Working fluid"
	EQUAL 			REAL A_out					UNITS u_m2			RANGE ZERO,Inf		"Nozzle output area"
	EQUAL 			REAL c_star					UNITS u_m_s			RANGE ZERO,Inf		"Characteristic velocity"
	EQUAL				REAL A_th					UNITS u_m2			RANGE ZERO,Inf		"Throat area"
		
END PORT

--------------------------------------------------------------------------------
-- Port for mechanical connections
--------------------------------------------------------------------------------
PORT Mechanical

	"Mechanical port"
	
   SUM   REAL Power 		UNITS u_W										"Mechanical power"
   EQUAL REAL N			UNITS u_rad_s		RANGE ZERO,Inf			"Rotational speed"

END PORT

--------------------------------------------------------------------------------
-- Port for heat exchanges
--------------------------------------------------------------------------------
PORT Heat SINGLE IN

	"Heat exchange port"
	
	SUM	REAL Q		UNITS u_W 									"Heat flux"
	EQUAL REAL T		UNITS u_K			RANGE ZERO,Inf		"Temperature" 
	EQUAL	REAL A		UNITS u_m2			RANGE ZERO,Inf		"Contact area"
	
END PORT

--------------------------------------------------------------------------------
-- Port for heat exchanges
--------------------------------------------------------------------------------
PORT Heat_2 SINGLE IN

	"Heat exchange port"
	
	SUM	REAL Q		UNITS u_W 									"Heat flux"
	EQUAL REAL T		UNITS u_K			RANGE ZERO,Inf		"Temperature" 
	EQUAL	REAL A		UNITS u_m2			RANGE ZERO,Inf		"Contact area"
	
END PORT
	
--------------------------------------------------------------------------------
-- Port for information exchanges between components
--------------------------------------------------------------------------------
PORT Info (INTEGER n)

	"Information exchange port"
	
	SUM REAL Data[n]			UNITS no_units 			"Information array"
	
END PORT
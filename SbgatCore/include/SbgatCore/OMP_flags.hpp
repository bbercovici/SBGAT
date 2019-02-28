/**
\file OMP_flags.hpp
\author Benjamin Bercovici
\author Jay McMahon
\date October 2018
\brief Defines a series of flags enabling/disabling OMP
\details Defines a series of flags enabling/disabling OMP
		Enable/Disable the corresponding functionality by setting the flag 
		- 1 (enabled)
		- 0 (disabled)
 \copyright MIT License, Benjamin Bercovici and Jay McMahon
*/

// Use OMP multithreading in ShapeModel methods
#define USE_OMP_SHAPE_MODEL 1

// Use OMP multithreading in DynamicAnalysis methods
#define USE_OMP_DYNAMIC_ANALYSIS 1

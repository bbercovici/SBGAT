
if (EXISTS /home/bebe0705/.am_fortuna)
	set(IS_FORTUNA ON)
	message("-- This is Fortuna")

else()
	set(IS_FORTUNA OFF)
endif()


if(${IS_FORTUNA})
	set(SBGATCORE_INCLUDE_HEADER /home/bebe0705/libs/local/include/SbgatCore/)
	set(SBGATCORE_LIBRARY /home/bebe0705/libs/local/lib/libSbgatCore.so)

else()
	set(SBGATCORE_INCLUDE_HEADER /usr/local/include/SbgatCore/)

	if (APPLE)
		set(SBGATCORE_LIBRARY /usr/local/lib/libSbgatCore.dylib)
	elseif(UNIX AND NOT APPLE)
		set(SBGATCORE_LIBRARY /usr/local/lib/libSbgatCore.so)
	else()
		message(FATAL_ERROR "Unsupported platform")
	endif()
endif()


message("-- Found SbgatCore: " ${SBGATCORE_LIBRARY})

set(SBGATCORE_INCLUDE_HEADER /usr/local/include/SbgatCore/)

if (APPLE)
	set(SBGATCORE_LIBRARY /usr/local/lib/libSbgatCore.dylib)
elseif(UNIX AND NOT APPLE)
	set(SBGATCORE_LIBRARY /usr/local/lib/libSbgatCore.so)
else()
	message(FATAL_ERROR "Unsupported platform")
endif()

message("-- Found SbgatCore: "${SBGATCORE_LIBRARY})


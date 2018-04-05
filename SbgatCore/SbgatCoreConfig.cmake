# MIT License

# Copyright (c) 2018 Benjamin Bercovici

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

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

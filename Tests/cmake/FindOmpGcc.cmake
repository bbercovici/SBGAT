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
#


# If running on a MAC, this will look for an OMP compliant compiler installed through Homebrew
# in /usr/local/Cellar
if(APPLE)

	# Checking if a built-from-source GCC lives in Homebrew's Cellar
	if(EXISTS /usr/local/Cellar/gcc)

		# Creating a glob storing the potential directories holding the compiler we want to use
		file(GLOB compiler_dirs /usr/local/Cellar/gcc/*)
		
		# Number of potential compilers
		list(LENGTH compiler_dirs len)

		# If len == 0, nothing to do here
		if(${len} EQUAL 0)
			message("No OMP-compliant compiler was found on this Mac.")
		else()
			# Looping over each directory to extract major/intermediate versions
			foreach(dir ${compiler_dirs})
				get_filename_component(name ${dir} NAME)

				string(REPLACE "." ";" split_name ${name})

				# major version
				list(GET split_name 0 major_version)
				list(APPEND major_version_list ${major_version})

				# intermediate version
				list(GET split_name 1 intermediate_version)
				list(APPEND intermediate_version_list ${intermediate_version})
			endforeach()

			# Finding the greatest major version 
			list(LENGTH compiler_dirs n_compilers)

			# If len == 1, there's only one compiler choice
			if (${n_compilers} EQUAL 1)
				list(GET compiler_dirs 0 OMP_FRIENDLY_GCC_PATH)
				list(GET major_version_list 0 compiler_major_version)
				set(OMP_FRIENDLY_GCC_PATH ${OMP_FRIENDLY_GCC_PATH}/bin/)
				set(OMP_FRIENDLY_GCC_MAJOR_VERSION ${compiler_major_version})

				message("Found OMP-compliant compiler: ${OMP_FRIENDLY_GCC_PATH}")
				
			else()
				# Sorting the compilers to find the most recent major version
				set(OMP_FRIENDLY_GCC_MAJOR_VERSION -1)
				set(major_index -1)

				foreach(major_inner ${major_version_list})
					MATH(EXPR major_index "${major_index}+1")
					if(major_inner GREATER OMP_FRIENDLY_GCC_MAJOR_VERSION)
						set(OMP_FRIENDLY_GCC_MAJOR_VERSION ${major_inner})
						set(OMP_FRIENDLY_GCC_MAJOR_VERSION_index ${major_inner})
					endif()
				endforeach()

				list(GET compiler_dirs ${major_index} OMP_FRIENDLY_GCC_PATH)
				list(GET major_version_list ${major_index} compiler_major_version)
				set(OMP_FRIENDLY_GCC_PATH ${OMP_FRIENDLY_GCC_PATH}/bin/)
				message("Found OMP-compliant compiler: ${OMP_FRIENDLY_GCC_PATH}")
				
			endif()	
			set(CMAKE_C_COMPILER ${OMP_FRIENDLY_GCC_PATH}gcc-${OMP_FRIENDLY_GCC_MAJOR_VERSION} CACHE STRING "C Compiler" FORCE)
			set(CMAKE_CXX_COMPILER ${OMP_FRIENDLY_GCC_PATH}g++-${OMP_FRIENDLY_GCC_MAJOR_VERSION} CACHE STRING "C++ Compiler" FORCE)

		endif()
	else()
		message("No OMP-compliant compiler was found on this Mac.")
	endif()

else() 
	# Running on Linux. Will switch back to compiler in /usr/local/bin
	message("Switching to /usr/local/gcc ")
	set(CMAKE_C_COMPILER "/usr/local/bin/gcc" CACHE STRING "C Compiler" FORCE)
	set(CMAKE_CXX_COMPILER "/usr/local/bin/g++" CACHE STRING "C++ Compiler" FORCE)
endif()

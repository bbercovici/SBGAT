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

# Proceeds to find the Qt install at some current locations

# if (${HOME_DIR})

# 	# Path to Qt5 cmake config file (author's mac configuration)
# 	set (QT_CMAKE_CONFIG_MAC "${HOME_DIR}/Qt/5.10.0/clang_64/lib/cmake")

# 	# Path to Qt5 cmake config file (author's ubuntu configuration)
# 	set (QT_CMAKE_CONFIG_LINUX "${HOME_DIR}/Qt/5.10.0/gcc_64/lib/cmake/")

# else()
	# Path to Qt5 cmake config file (author's mac configuration)
	set (QT_CMAKE_CONFIG_MAC "/usr/local/opt/qt")

	# Path to Qt5 cmake config file (author's ubuntu configuration)
	# set (QT_CMAKE_CONFIG_LINUX "~/Qt/5.10.0/gcc_64/lib/cmake/")

# endif()


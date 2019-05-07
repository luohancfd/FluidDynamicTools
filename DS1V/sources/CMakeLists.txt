cmake_minimum_required (VERSION 2.6)
project (DS1V)
enable_language (Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
# get_filename_component (default_prefix ".." ABSOLUTE)
# set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
#      "Choose the installation directory; by default it installs in the NORMA directory."
#      FORCE)

# FFLAGS depend on the compiler
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-fopenmp -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # ifort
  if (CMAKE_HOST_SYSTEM_NAME MATCHES "Windows")
    set (CMAKE_Fortran_FLAGS_RELEASE "/Qopenmp  /Qinit:zero /O3")
    set (CMAKE_Fortran_FLAGS_DEBUG   "/Qinit:zero /Od /debug:full")
  else ()
    set (CMAKE_Fortran_FLAGS_RELEASE "-qopenmp -init=zero -O3")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-init=zero -O0 -g")
  endif()
else()
  message("Untested compiler" ${CMAKE_Fortran_COMPILER_ID})
  set (CMAKE_Fortran_FLAGS_RELEASE "-fopenmp -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fopenmp -O0 -g")
endif()


# build executables
file (GLOB_RECURSE sources src/*.f90 src/*.h)
add_executable ("DS1_openmp_reac" ${sources})

# install executables and scripts
# install (TARGETS ${EXECUTABLES}
#          RUNTIME DESTINATION "bin")
# install (PROGRAMS ${SCRIPTS}
#          DESTINATION "bin")

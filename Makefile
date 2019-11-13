## Which compiler: (g++, clang++) [no spaces]
CXX=g++
#CXX=clang++

## Use OpenMP (parallelisation) yes/no:
UseOpenMP=yes
#UseOpenMP=no

## Build mode (changes warnings + optimisation level)
#Build=release
Build=dev
#Build=debug

## Optional: set directory for executables (by default: current directory)
XD=.

################################################################################
## None of the below options should need changing
################################################################################
## Set directories for source files (SD), and output object files (OD)
SD=./src
OD=./obj

## c++ standard. must be at least c++14
CXXSTD=-std=c++14
#CXXSTD=-std=c++17

## Build config + options:
include $(SD)/buildOptions.mk

## Build targets (must update if new programs/files are added):
include $(SD)/buildTargets.mk

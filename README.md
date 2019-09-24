# Relativistic, self-consistent atomic structure code.

Solves the Dirac equation for atomic systems using the Hartree-Fock method.
Fully relativistic, includes finite-nuclear size, and can solve for continuum states (energy normalisation).

 * With reasonable choices for the integration grids, typically converges to better than a few parts in 10^16
 * Includes an option to vary the effective speed of light - allowing non-relativistic approximation.
 * Orbitals are in form psi = (1/r) [f,ig], (using Dirac basis)

### Compiling and use:

 * Uses c++14
 * All programs compiled using the Makefile
 (run _$make_ or _$make programName_)
 * Must have LAPACK, GSL libraries installed already (see below)
 * Must create a directory called _./obj/_
   * code places object files inside here, but doesn't create the directory
 * All executables run like _$./programName_. By default, placed in same
  directory as makefile (you may change this in the Makefile)
 * All programs have input options, stored and read from a file.
By default, the input file name is: 'programName.in' (but can also give a different input file on runtime, e.g., to read input from
    `otherFile.txt`: _$./programName otherFile.txt_)
 * Note: below just tells how to use existing programs, to see how they work,
 see the comments/instructions inside the source code (all in /src/)
 * Tested with g++, clang++, and icpc.

Note: makes use of GSL libraries: https://www.gnu.org/software/gsl/

  * For example, with ubuntu: _$sudo apt-get install libgsl-dev_
  * Also needs LAPACK/BLAS libraries:
  _$sudo apt-get install libatlas-base-dev_ [and _liblapack-dev_ , but not _libblas-dev_ ?]

The above instructions are for linux (ubuntu). For windows, the easiest way (for me, anyway) is to make use of the recent 'windows subsystem for linux'. Instructions on installation/use here: https://www.roberts999.com/posts/2018/11/wsl-coding-windows-ubuntu
Then, the compilation + use can proceed as per above.


## hartreeFock (main program)

 * Solves relativistic Hartree Fock potential for core + valence states
 * Input taken from a plain text file.
 * See "hartreeFock.in" for minimal input example.
  May re-name this file (e.g., to "filename.txt"), then run as:
    * _$ ./hartreeFock filename.txt_
    * (Otherwise, program will assume file name is 'hartreeFock.in')
 * see README_input for a full list of input options + descriptions

--------------------------------------------------------------------------------

## Other programs:

### periodicTable

 * Gives info regarding particular element, including Z, default A, and electron configuration. Takes input in one line from command line.
 * Also looks up nuclear data from table and writes to screen.
 * Input is Z A; Z may be int or string. e.g.:
   * _$./periodicTable Cs 133_
 * Leave A blank (or put 0) to get default A value.
 * Put 'all' to list all available A values.
 * Note: numbers come from online database, and have some errors, so should be checked if needed.


### dmeXSection

 * Calculates the cross-section and event rates for ionisation of atoms
 by scattering of DM particle.
 * Takes in the output of "atomicKernal" (see README_input for details)
 * Also calculates "observable" event rates, accounting for detector thresholds
 and resolutions. For now, just for DAMA detector. Will add for XENON
 * For definitions/details, see:
   * B.M. Roberts, V.V. Flambaum, arXiv:1904.07127 (2019) [arXiv:1904.07127](https://arxiv.org/abs/1904.07127 "free download").
   * B.M.Roberts, V.A.Dzuba, V.V.Flambaum, M.Pospelov, Y.V.Stadnik,
 [Phys.Rev.D 93, 115037 (2016)](https://link.aps.org/doi/10.1103/PhysRevD.93.115037 "pay-walled");
 [arXiv:1604.04559](https://arxiv.org/abs/1604.04559 "free download").


### wigner

 * Small routine to calculate 3,6,9-j symbols, and Clebsch Gordan coeficients
 * Either give input via command line directly (quote marks required)
   * e.g., _./wigner '<0.5 -0.5, 0.5 0.5| 1 0>'_
   * or e.g., _./wigner '(0.5 1 0.5, -0.5 0 0.5)'_ etc.
 * Or, give an input file, that contains any number of symbols, all on new line
   * e.g., _./wigner -f myInputFile.in_
   * nb: the '-f' flag can be dropped in the '.in' file extension is used
   * Do not use quote marks in input file. Lines marked '!' or '#' are comments
 * 3j symbols must start with '('; 6,9j with '{', and CG with '<' (this is how code knows which symbol to calculate).
 * but, each number can be seperated by any symbol (space, comma etc.)

# Dependencies for Wavefunction

$(BD)/ContinuumOrbitals.o: $(SD)/Wavefunction/ContinuumOrbitals.cpp \
$(SD)/Wavefunction/ContinuumOrbitals.hpp $(SD)/DiracODE/DiracODE.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp $(SD)/Wavefunction/Wavefunction.hpp \
$(SD)/Maths/Grid.hpp \
$(SD)/Physics/AtomData.hpp $(SD)/Physics/PhysConst_constants.hpp
	$(COMP)

$(BD)/DiracSpinor.o: $(SD)/Wavefunction/DiracSpinor.cpp $(SD)/Wavefunction/DiracSpinor.hpp \
$(SD)/Maths/Grid.hpp $(SD)/Maths/NumCalc_quadIntegrate.hpp \
$(SD)/Physics/AtomData.hpp
	$(COMP)

$(BD)/BSplineBasis.o: $(SD)/Wavefunction/BSplineBasis.cpp $(SD)/Wavefunction/BSplineBasis.hpp\
$(SD)/Maths/LinAlg_MatrixVector.hpp $(SD)/Maths/BSplines.hpp \
$(SD)/Wavefunction/Hamiltonian.hpp \
$(SD)/Wavefunction/Wavefunction.hpp $(SD)/Wavefunction/DiracSpinor.hpp $(SD)/Maths/Grid.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp $(SD)/Physics/AtomData.hpp
	$(COMP)

$(BD)/Wavefunction.o: $(SD)/Wavefunction/Wavefunction.cpp $(SD)/Wavefunction/Wavefunction.hpp\
$(SD)/DiracODE/DiracODE.hpp $(SD)/Wavefunction/DiracSpinor.hpp $(SD)/Maths/Grid.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp $(SD)/Physics/AtomData.hpp \
$(SD)/Wavefunction/BSplineBasis.hpp \
$(SD)/Physics/NuclearPotentials.hpp $(SD)/Physics/PhysConst_constants.hpp
	$(COMP)
# Dependencies for MBPT

$(BD)/CorrelationPotential.o: $(SD)/MBPT/CorrelationPotential.cpp \
$(SD)/MBPT/CorrelationPotential.hpp \
$(SD)/MBPT/GreenMatrix.hpp \
$(SD)/Angular/Angular_tables.hpp \
$(SD)/Coulomb/Coulomb.hpp \
$(SD)/Coulomb/YkTable.hpp \
$(SD)/IO/FRW_fileReadWrite.hpp \
$(SD)/IO/SafeProfiler.hpp \
$(SD)/Maths/Grid.hpp \
$(SD)/Maths/Interpolator.hpp \
$(SD)/Maths/LinAlg_MatrixVector.hpp \
$(SD)/Physics/AtomData.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp \
$(SD)/Wavefunction/Wavefunction.hpp
	$(COMP)

$(BD)/DiagramRPA.o: $(SD)/MBPT/DiagramRPA.cpp \
$(SD)/MBPT/DiagramRPA.hpp \
$(SD)/Angular/Angular_369j.hpp \
$(SD)/Angular/Angular_tables.hpp \
$(SD)/Coulomb/Coulomb.hpp \
$(SD)/Coulomb/YkTable.hpp \
$(SD)/IO/FRW_fileReadWrite.hpp \
$(SD)/DiracOperator/DiracOperator.hpp \
$(SD)/IO/SafeProfiler.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp
	$(COMP)

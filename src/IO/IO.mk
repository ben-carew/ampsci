# Dependencies for IO

$(OD)/UserInput.o: $(SD)/IO/UserInput.cpp $(SD)/IO/UserInput.hpp \
$(SD)/IO/FileIO_fileReadWrite.hpp
	$(COMP)

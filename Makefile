alignPep.o: main.cu
	nvcc -std=c++11 main.cu device_functions.cu AlignmentMatrix.cpp FastaFile.cpp -o alignPep

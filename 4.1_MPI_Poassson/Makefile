CXX = mpic++
CXXFLAGS = -O3 -fopenmp

main:main.cpp Field.h 

all:main
	$(CXX) $(CXXFLAGS) 

clean:
	touch U.vtk F.vtk
	rm *.vtk main main
run:
	mpirun -n 8 ./main

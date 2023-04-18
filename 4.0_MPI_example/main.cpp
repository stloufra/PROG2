#include <mpi.h>
#include <iostream>
#include <math.h>

int main(int argc, char **argv) {

    MPI_Init(&argc,&argv);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	std::cout << "Procesor " << rank << "/" << size << "\n";

    double moje_pi = 0.0;
    const int n_max = 100000000/size;

    double pi_lokalni_soucet = 0.0;
	int j1 = rank*n_max;
	int j2 = (rank+1)*n_max;

    for (int j=j1;j<j2;++j) {
         int index = j;
         double znamenko = 1 - 2*(index % 2); // Timto nahradime vypocet mocniny (-1)^k pro stridani znamenka
         pi_lokalni_soucet += znamenko*4.0/double(2*index + 1);
    }
    std::cout << "lokalni soucet = " << pi_lokalni_soucet << "\n";

	MPI_Reduce(&pi_lokalni_soucet,&moje_pi,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

	if (rank == 0) {
    	std::cout << "chyba = " << M_PI - moje_pi << "\n";
	}

    MPI_Finalize();

	system("sleep 2");

	return 0;
}

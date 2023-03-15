// Program pro paralelni numericky vypocet hodnoty pi
// Aproximujeme pi jako soucet Taylorovy rady funkce 4*arctg(1)
// Jednotlive cleny sumy odpovidaji 4 * (-1)^i / (2*i + 1) pro i >= 0

// Paralelni vypocet pro obecny pocet vlaken pomoci OpenMP

// Verze 2 - automaticke rozdeleni prace mezi vlakna, update sdilene promenne pres atomickou operaci

// Rozdeleni iteraci mezi jednotliva vlakna provede openmp v cyklu nasledujicim direktivu #pragma omp for. 
// Nakonec castecne soucty secteme pomoci atomicke operace do sdilene promenne, cimz ziskame vysledek.

// Kod je nutne prelozit s parametrem prekladace -fopenmp

#include <omp.h>
#include <iostream>

int main() 
{
    double moje_pi = 0.0;
    const int pocet_iteraci = 1000000000;

    #pragma omp parallel
    {
	double s = 0.0;

        #pragma omp for
        for (int i=0; i < pocet_iteraci; ++i) {
            int znamenko = 1 - 2*(i % 2);
            s += znamenko * 4.0/double(2*i + 1);
        }

	// Soucet do sdilene promenne
        #pragma omp atomic  // rychlejsi nez zamek
        moje_pi += s;
    }

    std::cout << "pi = " << moje_pi << "\n";

    return 0;

} 
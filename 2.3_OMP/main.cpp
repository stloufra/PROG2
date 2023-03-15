// Program pro paralelni numericky vypocet hodnoty pi
// Aproximujeme pi jako soucet Taylorovy rady funkce 4*arctg(1)
// Jednotlive cleny sumy odpovidaji 4 * (-1)^i / (2*i + 1) pro i >= 0

// Paralelni vypocet pro obecny pocet vlaken pomoci OpenMP

// Verze 3 - pomoci paralelni redukce

// Kod je nutne prelozit s parametrem prekladace -fopenmp

#include <omp.h>
#include <iostream>

int main() 
{
    double moje_pi = 0.0;
    const int pocet_iteraci = 1000000000;

    #pragma omp parallel for reduction(+ : moje_pi)
    for (int i=0; i < pocet_iteraci; ++i) {
        int znamenko = 1 - 2*(i % 2);
        moje_pi += znamenko * 4.0/double(2*i + 1);
    }

    std::cout << "pi = " << moje_pi << "\n";

    return 0;

} 
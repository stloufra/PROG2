#include <omp.h>
#include <iostream>
#include <cmath>

// Program pro paralelni numericky vypocet hodnoty pi
// Aproximujeme pi jako soucet Taylorovy rady funkce 4*arctg(1)
// Jednotlive cleny sumy odpovidaji 4 * (-1)^i / (2*i + 1) pro i >= 0

// Paralelni vypocet pro obecny pocet vlaken pomoci OpenMP

// Verze 1 - manualni rozdeleni prace mezi vlakna, update sdilene promenne pres kritickou sekci

// Rozdelime sumu na bloky iteraci pro jednotliva vlakna, kazde s vlastni horni a dolni mezi indexu. 
// Castecne soucty secteme do sdilene promenne v kriticke sekci (pristup pouze jednoho vlakna najednou), 
// takze nedojde k desynchronizaci dat.

// Kod je nutne prelozit s parametrem prekladace -fopenmp

// systemova promena export "OMP_NUM_THREADS=2"

#include <omp.h>
#include <iostream>

int main() 
{
    double moje_pi = 0.0;
    const int pocet_iteraci = 1000000000;

    #pragma omp parallel
    {
        const int pocet_vlaken = omp_get_num_threads();
        const int cislo_vlakna = omp_get_thread_num();

        double s = 0.0;

        // Rozdelime iterace mezi jednotliva vlakna po blocich, definujeme dolni a horni mez pro kazde vlakno
    	int lokalni_pocet_iteraci = pocet_iteraci/pocet_vlaken;
        int dolni_mez = cislo_vlakna*lokalni_pocet_iteraci;
        int horni_mez;

        // V pripade, ze pocet iteraci neni delitelny poctem vlaken, posledni vlakno musi dopocitat vzdy az do posledni iterace
        if (cislo_vlakna == pocet_vlaken - 1) horni_mez = pocet_iteraci; 
        else horni_mez = (cislo_vlakna + 1)*lokalni_pocet_iteraci;

        // Zde probiha vypocet castecneho souctu Taylorovy rady funkce 4*arctg(1)
        for (int i=dolni_mez; i < horni_mez; ++i) {
            int znamenko = 1 - 2*(i % 2);
            s += znamenko * 4.0/double(2*i + 1);
        }

	// Soucet do sdilene promenne
        #pragma omp critical 
        {
            moje_pi += s;
        }
    }

    std::cout << "pi = " << moje_pi << "\n";

    return 0;

} 
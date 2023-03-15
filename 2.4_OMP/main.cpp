// Program pro paralelni numericky vypocet hodnoty pi
// Aproximujeme pi jako soucet Taylorovy rady funkce 4*arctg(1)
// Jednotlive cleny sumy odpovidaji 4 * (-1)^i / (2*i + 1) pro i >= 0

// Paralelni vypocet pro obecny pocet vlaken pomoci OpenMP

// Verze 4 - pomoci task-based paralelismu

// Nejprve definujeme dva ukoly - soucet lichych a sudych clenu. Tyto spustime paralelne.

// Pokud mame vice nez 2 vlakna k dispozici, vypocet lichych a sudych clenu muze sam o sobe probehnout paralelne pomoci N/2 vlaken.

// Kod je nutne prelozit s parametrem prekladace -fopenmp

#include <omp.h>
#include <iostream>
#include <exception>

// Funkce pro castecny soucet clenu rady, podle sablonoveho parametru vybereme kladne nebo zaporne
// Promenna n_vlaken definuje, kolik vlaken vyuzije funkce k vypoctu castecneho souctu
template <int ZNAMENKO>
double castecny_soucet_clenu(int n_iteraci, int n_vlaken) {
        double s = 0.0;

        #pragma omp parallel for reduction(+ : s) num_threads(n_vlaken)
        for (int i=0; i < n_iteraci; i++) {
            s += ZNAMENKO * 4.0/double(4*i + 2 - ZNAMENKO);
        }
 
        return s;
}

int main() 
{
    omp_set_nested(true); // Zapneme nested paralelismus

    double moje_pi = 0.0;
    const int pocet_iteraci = 1000000000;
    const int max_pocet_vlaken = omp_get_max_threads();
    if (max_pocet_vlaken < 2) throw std::runtime_error("Nelze spustit paralelni vypocet");

    const int pocet_vlaken = 2 * (max_pocet_vlaken / 2); // Vynutime sudy pocet vlaken

    // Rozdelime na dve vlakna, ale v tuto chvili spustime pouze ridici vlakno (master thread), proto #pragma omp single
    #pragma omp parallel num_threads(2)
    #pragma omp single
    {
        // Definujeme dva ukoly, ktere mohou byt spusteny paralelne - vypocet castecneho souctu kladnych a zapornych clenu
        #pragma omp task
        #pragma omp atomic
        moje_pi += castecny_soucet_clenu<1>(pocet_iteraci / 2, pocet_vlaken / 2);

        #pragma omp task
        #pragma omp atomic
        moje_pi += castecny_soucet_clenu<-1>(pocet_iteraci / 2, pocet_vlaken / 2);
    }

    std::cout << "pi = " << moje_pi << "\n";

    return 0;

} 
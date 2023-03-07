// Program pro paralelni numericky vypocet hodnoty pi
// Aproximujeme pi jako soucet Taylorovy rady funkce 4*arctg(1)
// Jednotlive cleny sumy odpovidaji 4 * (-1)^i / (2*i + 1) pro i >= 0

// Rozdelime sumu na castecne soucty sudych/lichych clenu a tyto poscitame paralelne ve dvou
// nezavislych vlaknech. Nakonec tyto castecne soucty secteme do globalni promenne, cimz ziskame
// vysledek. Tento krok musime provest za pomoci tzv. mutexu, ktery zaruci, ze vlakna budou pristupovat
// ke sdilene promenne vzdy postupne a nedojde k desynchronizaci dat.

// Kod je nutne prelozit s parametrem prekladace -lpthread

#include <iostream>
#include <thread>
#include <mutex>
#include <iomanip>

const int pocet_vlaken = 2;
double moje_pi = 0.0; 

// Zamek, ktery dovoli vstup do zamcene casti kodu pouze 1 vlaknu soucasne, ostatni musi pockat
std::mutex m; //mutual exclusion s

void spocitejSoucet(int cislo_vlakna) 
{

    // Zde probiha vypocet castecneho souctu Taylorovy rady funkce 4*arctg(1)
    // Vlakno 0 zpracovava pouze cleny se sudym indexem, vlakno 1 pouze liche cleny

    double s = 0.0;
    
    const int znamenko = 1 - 2*cislo_vlakna;

    for (int i=0; i < 10000; ++i) {

        // Castecny soucet, lokalni pro jedno vlakno
        s += znamenko * 4.0/double(2*(i*pocet_vlaken + cislo_vlakna) + 1);

    }

    // Zde zacina zamcena cast
    m.lock();

    // Nemuze tak zde dojit k soucasnemu pristupu vice vlaken k jedne globalni promenne
    // Castecne soucty se prictou do sdilene globalni promenne
    moje_pi += s;

    // Zde se odemkne
    m.unlock();

}

int main() 
{
    // Spustime funkci spocitejSoucet dvakrat jako 2 samostatna vlakna (cislo vlakna predame jako parametr)
    std::thread vlakno1(spocitejSoucet, 0);
    std::thread vlakno2(spocitejSoucet, 1);

    // Ukonceni samostatnych vlaken, pokracuje standardni seriovy program
    vlakno1.join();
    vlakno2.join();

    std::cout << std::setprecision(12) << "pi = " << moje_pi << "\n";

    return 0;

} 

/* 
Work Sharing 
    rozdělení práce 
                    -1 procesor 1-1000 
                    -2 procesor 1000-2000
    přístup k datům
                    -dle případu 
                    -nemusí si lézt do dat
                    -sdílení části paměti 
    !!! na vytváření vláken nejlépe pouze jednou !!! 
                    -vlakno vždy použivat blok
                    -

*/
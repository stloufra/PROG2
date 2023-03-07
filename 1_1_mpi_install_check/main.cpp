#include <iostream>
#include <omp.h>

int main() {
    #pragma omp parallel for
    for (int i = 0; i < 10; ++i) {
        std::cout << "Thread " << omp_get_thread_num() << " executing iteration " << i << std::endl;
    }
    return 0;
}
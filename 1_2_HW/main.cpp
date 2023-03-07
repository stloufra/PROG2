#include "Field.h"
#include <iostream>
#include <iomanip>
#include <chrono>


int main() {
	constexpr int N = 10000;
	constexpr int M = 10000;
	
	Field<double> f = Field<double>(N,M);

	for (int i=0; i<N; ++i) {
		for (int j=0; j<M; ++j) {
			f(i,j) = 1.0 + i/double(N) + j/double(M);
		}
	}

	std::cout << "Field created \n";

    auto stop1 = std::chrono::high_resolution_clock::now();

	std::cout << std::fixed << f.mean() << "\n";
	
	auto stop2 = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(stop2 - stop1).count() << " ms\n";

	auto stop3 = std::chrono::high_resolution_clock::now();

	std::cout << std::fixed << f.mean_paralel() << "\n";
	
	auto stop4 = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(stop4 - stop3).count() << " ms\n";


	return 0;
}
#include <iostream>
#include <fstream>
#include <chrono>

#include <TNL/Timer.h>
#include <TNL/Logger.h>
#include <TNL/Containers/Array.h>
#include <TNL/Containers/StaticArray.h>
#include <TNL/Algorithms/parallelFor.h>

using namespace TNL;
using namespace TNL::Containers;
using namespace TNL::Algorithms;

using DeviceType = Devices::Host;
using RealType = float;

using ArrayType = typename TNL::Containers::Array< RealType, DeviceType >;

template< typename Dev >






int main()
{
    const int N = 5000;
    const int it = 1000;

    
    
}


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

using DeviceType = Devices::Cuda;
using RealType = float;

using ArrayType = typename TNL::Containers::Array< RealType, DeviceType >;



//-------------------------FUNKCE-----------------------------
template< typename Dev >
void init_U(Array< RealType, Dev > &U, RealType &in, const int &N)
{
    auto Uview = U.getView();
    auto init = [=] __cuda_callable__ ( int i ) mutable 
    {
        Uview[i] = in;
    };

    parallelFor< Dev > (0, U.getSize(), init);
}

template< typename Dev >
void init_F(Array< RealType, Dev > &F, RealType &in1, RealType &in2, RealType &r, const int &N)
{
    auto Fview = F.getView();
    auto init = [=] __cuda_callable__ ( const StaticArray< 2, int >& i ) mutable
    {
        RealType lev = (RealType(i.x())/RealType(N)-0.5)*(RealType(i.x())/RealType(N)-0.5) + (RealType(i.y())/RealType(N)-0.5)*(RealType(i.y())/RealType(N)-0.5);
        if(lev < r)
        {
            Fview[ i.y()*N + i.x() ] = in2;
        }
        
        else
        {
            Fview[ i.y()*N + i.x() ] = in1;
        }
    };
   StaticArray< 2, int > begin{ 0, 0};
   StaticArray< 2, int > end{ N, N};
   parallelFor< Dev >( begin, end, init );
}

template< typename Dev >
void VTK_out_U(Array< RealType, Dev > &Un, const int &N)
{   
    if(std::is_same_v< Dev, TNL::Devices::Cuda > )
    {
        std::cout<<"\nCuda -> Host\n"<<std::endl;
        Array< RealType, TNL::Devices::Host > U( N*N );
        
        U = Un;

        std::ofstream out_file("U.vtk");

        out_file << "# vtk DataFile Version 2.0\n" ;
        out_file << "LBE mesh\n" ;
        out_file << "ASCII\n";
        out_file << "DATASET STRUCTURED_POINTS\n" ;
        out_file << "DIMENSIONS "<<N<<" "<<N<<" 1\n" ;
        out_file << "ASPECT_RATIO 1 1 1\n";
        out_file << "ORIGIN 0 0 0\n";
        out_file << "POINT_DATA "<< N*N <<"\n";
    
        out_file << "SCALARS "<<"U "<<"double 1\n";
        out_file << "LOOKUP_TABLE default\n";

        std::cout<<"Host -> .vtk\n"<<std::endl;
        for(int j=0; j < N; j++)
        {
            for(int i=0; i < N; i++)
            {
                out_file <<U[N*j + i]<<"\n";
            }
        }
    }
    else if (std::is_same<Dev, TNL::Devices::Host >::value )
    {
        std::ofstream out_file("U.vtk");

        out_file << "# vtk DataFile Version 2.0\n" ;
        out_file << "LBE mesh\n" ;
        out_file << "ASCII\n";
        out_file << "DATASET STRUCTURED_POINTS\n" ;
        out_file << "DIMENSIONS "<<N<<" "<<N<<" 1\n" ;
        out_file << "ASPECT_RATIO 1 1 1\n";
        out_file << "ORIGIN 0 0 0\n";
        out_file << "POINT_DATA "<< N*N <<"\n";
    
        out_file << "SCALARS "<<"U "<<"double 1\n";
        out_file << "LOOKUP_TABLE default\n";

        std::cout<<"Host -> .vtk\n"<<std::endl;
        for(int j=0; j < N; j++)
        {
            for(int i=0; i < N; i++)
            {
                out_file <<Un[N*j + i]<<"\n";
            }
        }
    }
    else
    {
        std::cout<<"Chybka se vloudila na line 120\n"<<std::endl;
    }
    
}

template< typename Dev >
void solver_Jacobi(Array< RealType, Dev > &Un,Array< RealType, Dev > &Un_1,Array< RealType, Dev > &F,RealType &h2, const int &N)
{   

    auto Un_view = Un.getView();
    auto Un_1_view = Un_1.getView();
    auto F_view = F.getView();
    
    auto iteration = [=] __cuda_callable__ ( const StaticArray< 2, int >& i ) mutable
    {
        Un_1_view[ i.y()*N + i.x() ] = 0.25*(h2*F_view[ i.y()*N + i.x() ] + Un_view[ i.y()*N + i.x() + 1 ]+ Un_view[ i.y()*N + i.x() - 1 ]+ Un_view[ (i.y()+1)*N + i.x() ]+ Un_view[ (i.y()-1)*N + i.x() ]);
    };
   StaticArray< 2, int > begin1{1, 1};
   StaticArray< 2, int > end1{N-1, N-1};
   parallelFor< Dev >( begin1, end1, iteration );

    auto copy_step = [=] __cuda_callable__ ( const StaticArray< 2, int >& i ) mutable
    {
        Un_view[ i.y()*N + i.x() ] =Un_1_view[ i.y()*N + i.x() ];
    };
   StaticArray< 2, int > begin2{ 1, 1};
   StaticArray< 2, int > end2{ N-1, N-1};
   parallelFor< Dev >( begin2, end2, copy_step );
}

//---------------------MAIN------------------------
// cca 13s

int main()
{
    const int N = 5000;
    const int it = 1000;

    RealType init_1 = 0.f;
    RealType init_2 = 50.f;
    RealType r= 0.04f;

    RealType h = 1/RealType(N);
    RealType h2 = h*h;


    ArrayType F( N*N );
    ArrayType Un( N*N );
    ArrayType Un_1( N*N );
    
    
    init_U(Un, init_1, N);
    init_U(Un_1, init_1, N);
    init_F(F, init_1, init_2, r, N);
    
    
    Timer timer;
    Logger logger(50, std::cout);

    timer.start();

    for(int n=0; n<it; n++)
    {
        solver_Jacobi(Un,Un_1,F,h2,N);
    }

    timer.stop();
    
    timer.writeLog( logger, 0 );

    VTK_out_U(Un, N);

    return  0;
    
}




/*
template< typename Dev >
void VTK_out_U(Array< double, Dev > &U, const int &N)
{   
    auto Uview = U.getView();
    std::ofstream out_file("U.vtk");

    out_file << "# vtk DataFile Version 2.0\n" ;
    out_file << "LBE mesh\n" ;
    out_file << "ASCII\n";
    out_file << "DATASET STRUCTURED_POINTS\n" ;
    out_file << "DIMENSIONS "<<N<<" "<<N<<" 1\n" ;
    out_file << "ASPECT_RATIO 1 1 1\n";
    out_file << "ORIGIN 0 0 0\n";
    out_file << "POINT_DATA "<< N*N <<"\n";
 
    out_file << "SCALARS "<<"U "<<"double 1\n";
    out_file << "LOOKUP_TABLE default\n";

    
    for(int j=0; j < N; j++)
    {
        for(int i=0; i < N; i++)
        {
            out_file <<Uview[N*j + i]<<"\n";
        }
    }
}
*/
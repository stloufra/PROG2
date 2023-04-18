#include <iostream>
#include <fstream>
#include <chrono>
#include <TNL/Containers/Array.h>
#include <TNL/Containers/StaticArray.h>
#include <TNL/Algorithms/parallelFor.h>

using namespace TNL;
using namespace TNL::Containers;
using namespace TNL::Algorithms;
using hrc = std::chrono::high_resolution_clock;

using device = Devices::Host;


void init_U(Array< double, device > &U, double &in, const int &N);
void VTK_out_U(Array< double, device > &U, const int &N);
void init_F(Array< double, device > &F, double &in1, double &in, double &r, const int &N);
void solver_Jacobi(Array< double, device > &Un,Array< double, device > &Un_1,Array< double, device > &F,double &h2, const int &N);


// cca 193s
int main()
{
   const int N = 10000;
    const int it = 1000;

    double init_1 = 0.0;
    double init_2 = 50.0;
    double r= 0.04;

    double h = 1/double(N);
    double h2 = h*h;


    Array< double, device > F( N*N );
    Array< double, device > Un( N*N );
    Array< double, device > Un_1( N*N );

    init_U(Un, init_1, N);
    init_U(Un_1, init_1, N);
    init_F(F, init_1, init_2, r, N);

    auto start = hrc::now();

    for(int n=0; n<it; n++)
    {
        solver_Jacobi(Un,Un_1,F,h2,N);
    }

    auto stop = hrc::now();
 
    VTK_out_U(Un, N);

    std::chrono::duration<double, std::micro> duration = std::chrono::duration_cast<std::chrono::duration<double, std::micro>>(stop - start);
    std::cout << "\nTime taken by function: " << duration.count()/1000000 << " seconds" << std::endl;
    
    return  0;
}

//-------------------------FUNKCE-----------------------------
void solver_Jacobi(Array< double, device > &Un,Array< double, device > &Un_1,Array< double, device > &F,double &h2, const int &N)
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
   parallelFor< device >( begin1, end1, iteration );

    auto copy_step = [=] __cuda_callable__ ( const StaticArray< 2, int >& i ) mutable
    {
        Un_view[ i.y()*N + i.x() ] =Un_1_view[ i.y()*N + i.x() ];
    };
   StaticArray< 2, int > begin2{ 1, 1};
   StaticArray< 2, int > end2{ N-1, N-1};
   parallelFor< device >( begin2, end2, copy_step );
}


void init_U(Array< double, device > &U, double &in, const int &N)
{
    auto Uview = U.getView();
    auto init = [=] __cuda_callable__ ( int i ) mutable 
    {
        Uview[i] = in;
    };

    parallelFor< device > (0, U.getSize(), init);
}

void init_F(Array< double, device > &F, double &in1, double &in2, double &r, const int &N)
{
    auto Fview = F.getView();
    auto init = [=] __cuda_callable__ ( const StaticArray< 2, int >& i ) mutable
    {
        double lev = (double(i.x())/double(N)-0.5)*(double(i.x())/double(N)-0.5) + (double(i.y())/double(N)-0.5)*(double(i.y())/double(N)-0.5);
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
   parallelFor< device >( begin, end, init );
}

void VTK_out_U(Array< double, device > &U, const int &N)
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

    
    for(int j=0; j < N; j++)
    {
        for(int i=0; i < N; i++)
        {
            out_file <<U[N*j + i]<<"\n";
        }
    }
}

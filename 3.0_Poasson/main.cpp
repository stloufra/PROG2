#include <cmath>
#include"Field.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>

using hrc = std::chrono::high_resolution_clock;

const int N = 1000;
const int iter = 1000;
double h = 1/double(N);
double h2 = h*h;

void init_prava(Field<double> &F);
void init_U(Field<double> &U);
bool vkruh(double j, double i);
void solver_Jacobi(Field<double> &Un,Field<double> &Un_1,Field<double> &F);
void solver_Jacobi_parallel(Field<double> &Un,Field<double> &Un_1,Field<double> &F);
void VTK_out_U(Field<double> &U);
void VTK_out_F(Field<double> &U);

int main()
{
    Field<double> F =Field<double>(N,N);
    Field<double> Un =Field<double>(N,N);
    Field<double> Un_1 =Field<double>(N,N);

    init_prava(F);
    init_U(Un);
    init_U(Un_1);

    VTK_out_F(F);

    //mainloop

    auto start = hrc::now();

    for(int n; n<N; n++)
    {
        solver_Jacobi_parallel(Un,Un_1,F);
    }

    auto stop = hrc::now();

    VTK_out_U(Un);

    std::cout << h << "\n";

    std::chrono::duration<double, std::micro> duration = std::chrono::duration_cast<std::chrono::duration<double, std::micro>>(stop - start);
    std::cout << "\nTime taken by function: " << duration.count()/1000000 << " seconds" << std::endl;

    return  0;
}




//--------------------FUNKCE-----------------------------

void solver_Jacobi(Field<double> &Un,Field<double> &Un_1,Field<double> &F)
{
    for(int j =1; j<N-1; j++)
    {
        for(int i =1; i<N-1; i++)
        {
            Un_1(j,i)=0.25*(h2*F(j,i) + Un(j-1,i) + Un(j+1,i)+ Un(j,i-1)+ Un(j,i+1));
            
        }
    }

    for(int j =1; j<N-1; j++)
    {
        for(int i =1; i<N-1; i++)
        {
            Un(j,i)=Un_1(j,i);
            
        }
    }

}


void solver_Jacobi_parallel(Field<double> &Un,Field<double> &Un_1,Field<double> &F)
{   
    #pragma omp parallel for collapse(2) 
    for(int j =1; j<N-1; j++)
    {
        for(int i =1; i<N-1; i++)
        {
            Un_1(j,i)=0.25*(h2*F(j,i) + Un(j-1,i) + Un(j+1,i)+ Un(j,i-1)+ Un(j,i+1));
            
        }
    }

    #pragma omp parallel for collapse(2) 
    for(int j =1; j<N-1; j++)
    {
        for(int i =1; i<N-1; i++)
        {
            Un(j,i)=Un_1(j,i);
            
        }
    }

}

void init_prava(Field<double> &F)
{
    for(int j =0; j<N; j++)
    {
        for(int i =0; i<N; i++)
        {
            if(vkruh(j,i))
            {
                F(j,i)=50;
            }
            else
            {
                F(j,i)=0;
            }
        }
    }

}

void init_U(Field<double> &U)
{
    for(int j =0; j<N; j++)
    {
        for(int i =0; i<N; i++)
        {
            U(j,i) = 0;
        }
    }

}

bool vkruh(double j, double i)
{
    double lev=(i/N-0.5)*(i/N-0.5) + (j/N-0.5)*(j/N-0.5);
    bool vevnitr;

    if(lev<0.04)
    {
        vevnitr=1;
    }
    else
    {
        vevnitr=0;
    }

    return vevnitr;

}

void VTK_out_U(Field<double> &U)
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
            out_file <<U(j,i)<<"\n";
        }
    }
}

void VTK_out_F(Field<double> &U)
{
    std::ofstream out_file("F.vtk");

    out_file << "# vtk DataFile Version 2.0\n" ;
    out_file << "LBE mesh\n" ;
    out_file << "ASCII\n";
    out_file << "DATASET STRUCTURED_POINTS\n" ;
    out_file << "DIMENSIONS "<<N<<" "<<N<<" 1\n" ;
    out_file << "ASPECT_RATIO 1 1 1\n";
    out_file << "ORIGIN 0 0 0\n";
    out_file << "POINT_DATA "<<N*N<<"\n";

    out_file << "SCALARS "<<"U "<<"double 1\n";
    out_file << "LOOKUP_TABLE default\n";

    
    for(int j=0; j < N; j++)
    {
        for(int i=0; i < N; i++)
        {
            out_file <<U(j,i)<<"\n";
        }
    }
}

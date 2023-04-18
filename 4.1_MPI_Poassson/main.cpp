#include <cmath>
#include"Field.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <cassert>
#include <vector> 
#include <mpi.h>

using hrc = std::chrono::high_resolution_clock;

const int N = 800;
const int iter = 1000;
double h = 1/double(N);
double h2 = h*h;
const bool message_debug = 0;
const bool mess_compose =0;

//--------------------FUNKCE-----------------------------

void solver_Jacobi(Field<double> &Un,Field<double> &Un_1,Field<double> &F, int &rank_, int &size_)
{
    int rank = rank_;
    int size = size_;
    
    int y1; int y2;

    //bounds
    if(rank ==0)
    {
        y1 = int(1);
        y2 = int(N/size);
    }
    else if(rank == size-1)
    {
        y1 = int(rank*N/size);
        y2 = int(N-1);
    }
    else
    {
        y1 = int(rank*N/size);
        y2 = int((rank+1)*N/size);
    }

    //iteration step
    for(int j = y1; j< y2; j++)
    {
        for(int i =1; i<N-1; i++)
        {
            Un_1(j,i)=0.25*(h2*F(j,i) + Un(j-1,i) + Un(j+1,i)+ Un(j,i-1)+ Un(j,i+1));
        }
    }

    //copy step
    for(int j =y1; j<y2; j++)
    {
        for(int i =1; i<N-1; i++)
        {
            Un(j,i)=Un_1(j,i);  
        }
    }

}

void solver_Jacobi_parallel(Field<double> &Un,Field<double> &Un_1,Field<double> &F, int &rank_, int &size_)
{   
    int rank = rank_;
    int size = size_;
    
    int y1; int y2;

    //bounds
    if(rank ==0)
    {
        y1 = int(1);
        y2 = int(N/size);
    }
    else if(rank == size-1)
    {
        y1 = int(rank*N/size);
        y2 = int(N-1);
    }
    else
    {
        y1 = int(rank*N/size);
        y2 = int((rank+1)*N/size);
    }

    //iteration step
    #pragma omp parallel for
    for(int j = y1; j< y2; j++)
    {
        for(int i =1; i<N-1; i++)
        {
            Un_1(j,i)=0.25*(h2*F(j,i) + Un(j-1,i) + Un(j+1,i)+ Un(j,i-1)+ Un(j,i+1));
        }
    }

    //copy step
    #pragma omp parallel for
    for(int j =y1; j<y2; j++)
    {
        for(int i =1; i<N-1; i++)
        {
            Un(j,i)=Un_1(j,i);  
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

bool v_rec(double j, double i)
{
    bool vevnitr;

    if(i>0.4*N && i < 0.6*N)
    {
        vevnitr = 1;
    }
    else
    {
        vevnitr = 0;
    }


    return vevnitr;
}

void init_prava(Field<double> &F,int &rank_, int &size_)
{
    int rank = rank_;
    int size = size_;

    int y1; int y2;

    if(rank ==0)
    {
        y1 = 0;
        y2 = int(N/size)+1;
    }
    else if(rank == size-1)
    {
        y1 = int(rank*N/size -1);
        y2 = N;
    }
    else
    {
        y1 = int(rank*N/size -1);
        y2 = int((rank+1)*N/size +1);
    }

    for(int j = y1 ; j < y2; j++)
    {
        for(int i =0; i<N ;i++)
        {
            if(vkruh(j,i))
            //if(v_rec(j,i))
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

void init(Field<double> &U)
{
    for(int j =0; j<N; j++)
    {
        for(int i =0; i<N; i++)
        {
            U(j,i) = 0;
        }
    }

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

void VTK_out_U_i(Field<double> &U, int rank)
{
    std::string rank_string = std::to_string(rank);
    std::ofstream out_file("U"+rank_string+".vtk");

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

void VTK_out_F(Field<double> &U, int rank)
{
    std::string rank_string = std::to_string(rank);
    std::ofstream out_file("F"+rank_string+".vtk");

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

void send_and_recieve(int &rank_, int &size_, Field<double> &U)
{
    int rank = rank_;
    int size = size_;
    
    std::vector<double> send_l(N);
    std::vector<double> send_u(N);

    std::vector<double> recv_l(N);
    std::vector<double> recv_u(N);

    int y_u; int y_l;

    if(rank ==0)
    {   y_u = int(N/size - 1);

        if(message_debug)
        {
            std::cout<< y_u <<" this is y_u\n";
        }

        for(int x = 0; x < N; x++)
        {
            send_u[x] = U(y_u,x);
        }
    }
    else if(rank == size-1)
    {
        y_l = int(rank*N/size );
        for(int x = 0; x < N; x++)
        {
            send_l[x] = U(y_l, x);
        }
    }
    else 
    {
        y_l = int(rank*N/size);
        y_u = int((rank+1)*N/size-1);

        for(int x = 0; x < N; x++)
        {
            send_l[x] = U(y_l, x);
            send_u[x] = U(y_u,x);
        }  
    }


    if(rank == 0)
    {  
        MPI_Status status_u;

        MPI_Send(send_u.data(), send_u.size() , MPI_DOUBLE, 1, 0, MPI_COMM_WORLD); 
        if(message_debug)
        {
            std::cout<<"Rank: " << rank << " send its uppper to: 1.\n";
            for(int i = 0; i<N; i++)
            {
                std::cout<< send_u[i] <<", ";
            }
            std::cout<<"<- this was send. \n";

            for(int i = 0; i<N; i++)
            {
                std::cout<< U(y_u,i) <<", ";
            }
            std::cout<<"<- this should have been send. \n";
        }

        MPI_Recv(recv_u.data(), recv_u.size(), MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status_u);
        if(message_debug)
        {
            std::cout<<"Rank: "<< rank <<" recieved its upper from: 1 sucessfully.\n";
            for(int i = 0; i<N; i++)
            {
                std::cout<< recv_u[i] <<", ";
            }
            std::cout<<"<- this was recieved. \n";
        } 
    }

    else if(rank == size-1)
    {
        MPI_Status status_l;
        MPI_Recv(recv_l.data(), recv_l.size(), MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status_l);
        if(message_debug)
        {
            std::cout<<"Rank: "<< rank << " recieved its lower from: "<< rank -1 <<" sucessfully.\n";
        }

        MPI_Send(send_l.data(), send_l.size() , MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD); 
        if(message_debug)
        {
            std::cout<<"Rank: "<<rank<<" send its lower to: "<< rank -1 <<".\n";
        } 
    }
    else 
    {   
        MPI_Status status_l;
        MPI_Status status_u;

        MPI_Recv(recv_l.data(), recv_l.size(), MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status_l);
        if(message_debug)
        {
            std::cout<<"Rank: "<<rank<<" recieved its lower from: "<< rank -1 <<" sucessfully.\n";
        }

        MPI_Send(send_u.data(), send_u.size(), MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
        if(message_debug)
        {
            std::cout<<"Rank: "<<rank<<" send its upper to: "<< rank +1 <<".\n";
        } 

        MPI_Recv(recv_u.data(), recv_u.size(), MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &status_u);
        if(message_debug)
        {
            std::cout<<"Rank: "<<rank<<" recieved its upper from: "<< rank +1 <<" sucessfully.\n";
        } 

        MPI_Send(send_l.data(), send_l.size(), MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
        if(message_debug)
        {
            std::cout<<"Rank: "<<rank<<" send its lower to: "<< rank -1 <<".\n";
        } 
    }

    if(rank ==0)
    {   int y_u = int(N/size);
        for(int x = 0; x < N; x++)
        {
            U(y_u,x) = recv_u[x];
        }
    }
    else if(rank == size-1)
    {
        int y_l = int(rank*N/size -1);
        for(int x = 0; x < N; x++)
        {
            U(y_l, x) = recv_l[x];
        }
    }
    else 
    {
        int y_l = int(rank*N/size -1);
        int y_u = int((rank+1)*N/size );
        for(int x = 0; x < N; x++)
        {
            U(y_l, x) = recv_l[x];
            U(y_u,x) = recv_u[x];
        }  
    }

}

void send_and_recieve_not_blocking(int &rank_, int &size_, Field<double> &U)
{
    int rank = rank_;
    int size = size_;
    
    std::vector<double> send_l(N);
    std::vector<double> send_u(N);

    std::vector<double> recv_l(N);
    std::vector<double> recv_u(N);

    int y_u; int y_l; int n;

    if(rank ==0)
    {   y_u = int(N/size - 1);

        if(message_debug)
        {
            std::cout<< y_u <<" this is y_u\n";
        }

        for(int x = 0; x < N; x++)
        {
            send_u[x] = U(y_u,x);
        }

        n = 2;
    }
    else if(rank == size-1)
    {
        y_l = int(rank*N/size );
        for(int x = 0; x < N; x++)
        {
            send_l[x] = U(y_l, x);
        }

        n = 2;
    }
    else 
    {
        y_l = int(rank*N/size);
        y_u = int((rank+1)*N/size-1);

        for(int x = 0; x < N; x++)
        {
            send_l[x] = U(y_l, x);
            send_u[x] = U(y_u,x);
        } 

        n = 4; 
    }

    MPI_Request requests[n];

    if(rank == 0)
    {  
        
        MPI_Isend(send_u.data(), send_u.size() , MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &requests[0]); 
        if(message_debug)
        {
            std::cout<<"Rank: " << rank << " send its uppper to: 1.\n";
            for(int i = 0; i<N; i++)
            {
                std::cout<< send_u[i] <<", ";
            }
            std::cout<<"<- this was send. \n";

            for(int i = 0; i<N; i++)
            {
                std::cout<< U(y_u,i) <<", ";
            }
            std::cout<<"<- this should have been send. \n";
        }

        MPI_Irecv(recv_u.data(), recv_u.size(), MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &requests[1]);
        if(message_debug)
        {
            std::cout<<"Rank: "<< rank <<" recieved its upper from: 1 sucessfully.\n";
            for(int i = 0; i<N; i++)
            {
                std::cout<< recv_u[i] <<", ";
            }
            std::cout<<"<- this was recieved. \n";
        } 
    }

    else if(rank == size-1)
    {
        
        MPI_Irecv(recv_l.data(), recv_l.size(), MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &requests[0]);
        if(message_debug)
        {
            std::cout<<"Rank: "<< rank << " recieved its lower from: "<< rank -1 <<" sucessfully.\n";
        }

        MPI_Isend(send_l.data(), send_l.size() , MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &requests[1]); 
        if(message_debug)
        {
            std::cout<<"Rank: "<<rank<<" send its lower to: "<< rank -1 <<".\n";
        } 
    }
    else 
    {   
        

        MPI_Irecv(recv_l.data(), recv_l.size(), MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &requests[0]);
        if(message_debug)
        {
            std::cout<<"Rank: "<<rank<<" recieved its lower from: "<< rank -1 <<" sucessfully.\n";
        }

        MPI_Isend(send_u.data(), send_u.size(), MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &requests[1]);
        if(message_debug)
        {
            std::cout<<"Rank: "<<rank<<" send its upper to: "<< rank +1 <<".\n";
        } 

        MPI_Irecv(recv_u.data(), recv_u.size(), MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &requests[2]);
        if(message_debug)
        {
            std::cout<<"Rank: "<<rank<<" recieved its upper from: "<< rank +1 <<" sucessfully.\n";
        } 

        MPI_Isend(send_l.data(), send_l.size(), MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &requests[3]);
        if(message_debug)
        {
            std::cout<<"Rank: "<<rank<<" send its lower to: "<< rank -1 <<".\n";
        } 
    }

    MPI_Waitall(n, requests, MPI_STATUS_IGNORE);

    if(rank ==0)
    {   int y_u = int(N/size);
        for(int x = 0; x < N; x++)
        {
            U(y_u,x) = recv_u[x];
        }
    }
    else if(rank == size-1)
    {
        int y_l = int(rank*N/size -1);
        for(int x = 0; x < N; x++)
        {
            U(y_l, x) = recv_l[x];
        }
    }
    else 
    {
        int y_l = int(rank*N/size -1);
        int y_u = int((rank+1)*N/size );
        for(int x = 0; x < N; x++)
        {
            U(y_l, x) = recv_l[x];
            U(y_u,x) = recv_u[x];
        }  
    }

}


void compose(int &rank_, int &size_, Field<double> &U)
{
    int rank = rank_;
    int size = size_;

    if(rank == 0)
    {   
        const int vel = int(N*N/size);
        std::vector<double> rec(vel);

        for(int s = 1; s < size; s++)
        {   
            MPI_Status status;
            MPI_Recv(rec.data(), rec.size(), MPI_DOUBLE, s, 0, MPI_COMM_WORLD, &status);

            if(mess_compose)
            {
                std::cout<<"Rank: "<< rank <<" recieved from: "<< s<<" sucessfully.\n";
            }

            int y1; int y2;
            y1 = int(s*N/size);
            y2 = int((s+1)*N/size);

            for(int j = y1 ; j < y2; j++)
            {
                for(int i =0; i<N ;i++)
                {
                    U(j,i) = rec[N*(j-y1)+i];
                }
            }   
        
        }
    } 

    else if(rank == size-1)
    {
        const int vel = int(N*N/size);

        int y1; int y2;
        y1 = int(rank*N/size);
        y2 = int(N);

        if(mess_compose)
        {
            std::cout<<"Velikost: "<<vel<< "; rank: "<< rank << "; y1, y2: " << y1 <<","<< y2 <<", velikost z nich: "<<  (y2-y1)*N <<"\n"; 
        }

        assert(vel == (y2-y1)*N);

        std::vector<double> send(vel);

        for(int j = y1 ; j < y2; j++)
        {
            for(int i =0; i<N ;i++)
            {
                send[N*(j-y1)+i]=U(j,i);
            }
        }

        MPI_Send(send.data(), send.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    else
    {   
        const int vel = int(N*N/size);

        int y1; int y2;
        y1 = int(rank*N/size);
        y2 = int((rank+1)*N/size);

        if(mess_compose)
        {
            std::cout<<"Velikost: "<<vel<< "; rank: "<< rank << "; y1, y2: " << y1 <<","<< y2 <<", velikost z nich: "<<  (y2-y1)*N <<"\n"; 
        }

        assert(vel == (y2-y1)*N);

        std::vector<double> send(vel);

        for(int j = y1 ; j < y2; j++)
        {
            for(int i =0; i<N ;i++)
                {
                    send[N*(j-y1)+i]=U(j,i);
                }
        }

        MPI_Send(send.data(), send.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
}

int main(int argc, char **argv) {

    MPI_Init(&argc,&argv);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

    assert(N%size == 0);

    std::cout<<"I am rank/size :"<< rank << "/" << size << std::endl;
    
    Field<double> F =Field<double>(N,N);
    Field<double> Un =Field<double>(N,N);
    Field<double> Un_1 =Field<double>(N,N);

  
    
    //------------here------------------
    init(Un);
    init(Un_1);
    init(F);

    init_prava(F, rank, size);

    VTK_out_F(F, rank);
    
    //mainloop

    // poiter to start, nuber of element, data type, where rank?, tag - defined by me, comm_world)

    
    auto start = hrc::now();

    for(int n =0; n < iter; n++)
    {
        solver_Jacobi_parallel(Un,Un_1,F,rank,size);
        send_and_recieve_not_blocking(rank, size, Un);
        if(rank==0)
        {
            std::cout<<"Iteration "<<n<<"\n";
        }
    }

    auto stop = hrc::now();
    
    VTK_out_U_i(Un, rank);

    compose(rank, size,Un);

    if(rank == 0)
    {
        VTK_out_U(Un);
    }
    
    std::chrono::duration<double, std::micro> duration = std::chrono::duration_cast<std::chrono::duration<double, std::micro>>(stop - start);
    std::cout << "\nTime taken by function: " << duration.count()/1000000 << " seconds" << std::endl;
    
    MPI_Finalize();

    return  0;
}



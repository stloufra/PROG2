#ifndef VARIABLES_H
#define VARIABLES_H
#include <vector>
#include <cassert>

#pragma once
template <typename T>
class Variables
{
public:
    Variables()=delete;
    Variables(int r, int c)
    {
        rows_pr = r; 
        collums_pr = c;
        data.resize(rows_pr*collums_pr);
        collums = c;
        rows = r;
    }
    ~Variables() {} 

    T &operator() (int i, int j)
    {   
        assert(i >= 0 && i <= rows_pr && j >= 0 && j <= collums_pr );

        return  data[i*collums_pr + j ];
    }

    const T &operator() (int i, int j) const
    {   
        assert(i >= 0 && i <= rows_pr && j >= 0 && j <= collums_pr );

        return data[i*collums_pr + j ];
    }
    
    operator + (Variables<T> A, Variables<T> B)
    {
        assert(A.rows == B.rows && A.collums == B.colums);        

        for(int i=0; i<prvni.collums; i++)
        {
            for(int j=0; j<prvni.rows; j++)
            {
                C(j,i) = A(j,i) + B(j,i)
            }
        }
    }

    T mean()
    {
        T prum;
        T sum = data(1,1);
        
       for(int i=0; i<prvni.collums; i++)
        {
            for(int j=0; j<prvni.rows; j++)
            {
                sum =+ data(j,i);
            }
        }

        sum = sum - data(1,1);

        prum = sum / (rows_pr*collums_pr);

    }

    int rows;
    int collums;

private:
    
    std::vector<T> data;
    int rows_pr;
    int collums_pr;
};


#endif

//zapsat do nej jakoukoliv promenou typ
//operace pres cely field
//prumer pres cele pole 
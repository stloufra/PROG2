#ifndef FIELD_H
#define FIELD_H
#include <vector>
#include <cassert>

#pragma once
template <typename T>
class Field
{
public:
    Field()=delete;
    Field(int r, int c)
    {
        rows_pr = r; 
        collums_pr = c;
        data.resize(rows_pr*collums_pr);
        collums = c;
        rows = r;
    }
    ~Field() {} 

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
    
    Field<T> operator+(const Field<T>& other) const 
    {
        Field<T> result(*this);
        result += other;
        return result;
    }   

    Field<T>& operator+=(const Field<T>& other) {
        
        assert(data.size()==other.data.size());

        for (int i = 0; i < data.size(); ++i) 
        {
            data[i] += other.data[i];
        }

        return *this;
    }

    Field<T> operator*(double k) {

        Field<T> result(*this);

        for (int i = 0; i < data.size(); i++)
        {
            result.data[i]=k * data[i];
        }

        return result;
    }



    Field<T>& operator=(const Field& other) {

        
        for (int i = 0; i < data.size(); i++)
        {
           data[i]= other.data[i];
        }

        return *this;
    }
    
    
    T mean()
    {
        T prum;
        T sum = data[0];
        T j = 1;
        
       for(int i=1; i< collums_pr*rows_pr; i++)
        {
            sum += data[i];
            j +=1;
        }

        return prum = sum / j;

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
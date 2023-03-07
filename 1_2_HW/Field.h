#ifndef FIELD_H
#define FIELD_H
#include <vector>
#include <cassert>
#include <omp.h>
#include <mutex>
#include <thread>

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
        pocet_vlaken = 1;
                
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

    T mean_paralel()
    {   
       T prum;
       Sum = data[0] - data[0];

       pocet_vlaken = omp_get_max_threads();

       std::vector<std::thread> threads;

       int data_size = rows_pr*collums_pr;
       int zbytek_rozdeleni = (data_size)%pocet_vlaken;
       int delka_deleni = (data_size - zbytek_rozdeleni)/pocet_vlaken;

       for (int i=0; i < pocet_vlaken; i++)
       {
            threads.push_back(std::thread(&Field<T>::spocitej_soucet,this,i*delka_deleni,(i+1)*delka_deleni-1 ));
       }

        for (auto& t : threads)
        {
            t.join();
        }


        for (int i=delka_deleni*pocet_vlaken; i < delka_deleni*pocet_vlaken + zbytek_rozdeleni;i++)
        {
                Sum += data[i];
        }

        return prum = Sum/data_size;

    }

    void spocitej_soucet(int index_prvni, int index_posledni)
    {
        T Sum_thrd = data[index_prvni];

        for (int i=index_prvni+1; i < index_posledni;i++)
        {
            Sum_thrd += data[i];
        }

        m.lock();
        Sum += Sum_thrd;
        m.unlock();

    }

    int rows;
    int collums;

private:
    
    std::vector<T> data;
    int rows_pr;
    int collums_pr;
    int pocet_vlaken;

    std::mutex m;
    T Sum;
};


#endif

//zapsat do nej jakoukoliv promenou typ
//operace pres cely field
//prumer pres cele pole 
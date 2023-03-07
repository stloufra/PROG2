#include <iostream>
#include <vector>
#include "Field.h"
#include "vector3d.h"



int main ()
{   

    Field<double> typ1_1(7, 5);
    Field<double> typ1_2(7, 5);
    Field<double> typ1_3(7, 5);
    Field<double> typ1_sum(7, 5);

    Field<vector3d> typ2(7,5);
    vector3d h;
    
    h(0) = 1;
    h(1) = 2;
    h(2) = 3;
    
    for(int i=0; i<typ1_1.collums; i++)
    {
        for(int j=0; j<typ1_1.rows; j++)
        {
            typ1_1(j,i) = 0.45;
            typ1_2(j,i) = 1.1;

            typ2(j,i) = h;
        }
    }
    //Field<std::vector<3>> druha(7, 5);

    typ1_sum = typ1_1 + typ1_2;
    typ1_3 = typ1_2 * 5.0;
    auto typ1_4 = typ1_3;

    std::cout<< typ1_1.mean() <<"\n";
    std::cout<< typ1_3(2,3) <<"\n";
    std::cout<< typ1_4(2,3) <<"\n";
    std::cout<< typ1_4.rows <<"\n";
    std::cout<< typ1_sum.mean() <<"\n";
    std::cout<< typ2(2,3)(2) <<"\n";
    return 0;
}

// poiter 
//reference - rychlejsi 
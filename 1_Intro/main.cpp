#include <iostream>
#include <vector>
#include "Variables.h"


int main ()
{   

    Variables<double> prvni(7, 5);
    Variables<double[3]> druhy(7,5);

    for(int i=0; i<prvni.collums; i++)
    {
        for(int j=0; j<prvni.rows; j++)
        {
            prvni(j,i) = 0.45;
            druhy(j,i) = {1,2,3};
        }
    }
    //Variables<std::vector<3>> druha(7, 5);


    std::cout<<prvni(2,3)  <<"\n";
    return 0;
}
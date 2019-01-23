#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

#include <Eigen/Dense>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include "Common.h"
#include "Watch.h"
#include "marray.h"
#include "LR.h"
#include "invx.h"
#include <unistd.h> 
#include <pthread.h> 


//-------------------------------------------------------------------------------
int main(int argc, char const *argv[]){
    main_invx(argc, argv);
    //test1();
}

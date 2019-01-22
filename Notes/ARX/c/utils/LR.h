#ifndef LINEARREGRESSION_HPP
#define LINEARREGRESSION_HPP
#include <Eigen/Dense>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include "Common.h"
#include "marray.h"
#include "CSV.h"

void SetVal(Eigen::MatrixXd& m, int column, double val);
void SetVal(Eigen::MatrixXd& m, int column, const Eigen::VectorXd& val);
void p(const Eigen::VectorXd& m, const char* str="", int all=0);
Eigen::MatrixXd& RemoveColumn(Eigen::MatrixXd& m, int colToRemove);
Eigen::VectorXd& RemoveColumn(Eigen::VectorXd& m, int colToRemove);

class LinearRegression{
public:
    char             *fileName;
    CSV              *csvp;
    Eigen::VectorXd  theta;
    Eigen::VectorXd  Jhistory;
    Eigen::VectorXd  res;
    double           intercept;
    Eigen::VectorXd  mean;
    Eigen::VectorXd  std;
        
    LinearRegression();
    LinearRegression(const char* file);
    ~LinearRegression() { 
        if ( fileName) 
            free( fileName);
        if ( csvp)
            delete csvp;
    }

    void Dump();
    
    
    static Eigen::VectorXd Mean( Eigen::MatrixXd & m );
    static Eigen::VectorXd Std( Eigen::MatrixXd & m, const Eigen::VectorXd & mean );
    Eigen::VectorXd& CalculateIntercept();
    static double CalculateIntercept(Eigen::VectorXd& m, Eigen::VectorXd& t);
    
    Eigen::VectorXd& Normalize( Eigen::MatrixXd & x, int doMean = 1, int doStd=0 );

    double Cost( const Eigen::MatrixXd & x, const Eigen::VectorXd & y, 
                           const Eigen::VectorXd & theta );
    Eigen::VectorXd GD( const Eigen::MatrixXd & x, const Eigen::VectorXd & y, 
                   int num_iters = 1024 * 50 , int batch=10, 
                   double alpha=0.005001, double eps=1e-7);
    
    Eigen::VectorXd& SVD( const Eigen::MatrixXd & x, const Eigen::VectorXd & y);
    Eigen::VectorXd& QR( const Eigen::MatrixXd & x, const Eigen::VectorXd & y);
    Eigen::VectorXd& Norm( const Eigen::MatrixXd & x, const Eigen::VectorXd & y);
};

#endif


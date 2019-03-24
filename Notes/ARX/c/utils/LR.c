#include <stdio.h>
#include <stdlib.h>
#include "marray.h"
#include "CSV.h" 
#include "LR.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

void GetMatrix(const CSV& csv, Eigen::MatrixXd & x, int getall=0, int uniq=3);


//EigenUtilities:
void p(const Eigen::VectorXd& m, const char* str, int all) {
    printf("%s #%ld :", str, m.size());
    int i;
    int ht = 6; //print 6 before and from last

    if(all || m.size() <= 2 * ht ) {
        for (i=0; i < m.size(); i++)
            printf("%.4f ", m[i]);
        printf("\n");
        return;
    }
        
    for (i=0; i < ht; i++)
        printf("%.4f ", m[i]);
    printf(" ... ");
    for (i=m.size()-ht; i < m.size(); i++)
        printf("%.4f ", m[i]);
    printf("\n");
}

void SetVal(Eigen::MatrixXd& m, int column, const Eigen::VectorXd& val) {
    if ( column < 0 ) {
        column = m.cols() + column ;
    }
    for( int i = 0; i < m.rows(); ++i ) {
        m(i, column) = val[i];
        printf("==> %f %f ", m(i, column) , val[i]);
    }
}

void SetVal(Eigen::MatrixXd& m, int column, double val) {
    if ( column < 0 ) {
        column = m.cols() + column ;
    }
    for( int i = 0; i < m.rows(); ++i ) {
        m(i, column) = val;
    }
}

Eigen::MatrixXd& RemoveColumn(Eigen::MatrixXd& m, int colToRemove) {
    unsigned int nRows = m.rows();
    unsigned int nCols = m.cols()-1;

    int c = colToRemove; 
    if ( colToRemove < 0 ) {
        colToRemove = c = nCols+ c+1;
    }
    //cout << "Removing Column " << colToRemove <<endl;    
    if( colToRemove < nCols )
        m.block(0,c,nRows,nCols-c) = m.block(0,c+1,nRows,nCols-c);
    m.conservativeResize(nRows,nCols);
    return m;
}

// Start columns is zero
Eigen::VectorXd& RemoveColumn(Eigen::VectorXd& m, int colToRemove) {
    Eigen::VectorXd t(m);
    
    int c = colToRemove; 
    if ( colToRemove < 0 ) {
        colToRemove = c = m.size()+ c+1;
    }
    //printf("Will remove column %d\n", c);
    
    m.resize(m.rows()-1);
    if ( c == 0)
        m << t.segment(1, t.rows()-1);
    else if ( c == t.rows() )
        m << t.segment(0, t.rows()-1);
    else 
        m << t.segment(0, c), t.segment(c+1, t.rows()-c-1);
    return m;
}

//-------------------------------------------------------------------------
LinearRegression::LinearRegression(): fileName(NULL), csvp(NULL){
}
LinearRegression::LinearRegression(const char* f){
    fileName = strdup(f);
    csvp = new CSV;
    csvp->Read(fileName);
    csvp->Dump();
}


Eigen::VectorXd LinearRegression::Mean( Eigen::MatrixXd & m ){
    Eigen::VectorXd mean( m.cols() );
    for( int i = 0; i < m.cols(); ++i ){
        mean[ i ] = m.col( i ).mean();
    }
    return mean;
}

Eigen::VectorXd LinearRegression::Std( Eigen::MatrixXd & m, const Eigen::VectorXd & mean ){
    double variance; 
    Eigen::VectorXd std( m.cols() );
    for( int i = 1; i < m.cols(); ++i ){
        variance = ( ( m.array().col( i ) - mean[ i ] ).array().pow( 2 ) ).sum();
        std[ i ] = sqrt( variance / ( m.rows() - 1 ) );
        variance = 0;
    }
    return std;
}

 
Eigen::VectorXd& LinearRegression::Normalize( Eigen::MatrixXd & x, int doMean, int doStd ){
    mean = Mean( x );
    std  = Std( x, mean );
    
    for( int i = 1; i < x.cols(); ++i )
        x.array().col( i ) -= mean[ i ];

    if ( doStd ) {
        for( int k = 1; k < x.cols(); ++k ){
            if (std[k] != 0)
                x.col( k ) /= std[ k ];
        }
        for( int k = x.cols()-1; k > 0 ; --k ){
            if (std[k] == 0) {
                RemoveColumn(x, k);
            }
        }
    }
    return mean;
}

double LinearRegression::Cost( const Eigen::MatrixXd & x, const Eigen::VectorXd & y, const Eigen::VectorXd & theta ){
    int m = y.size();
    Eigen::VectorXd result = ( x * theta - y );
    result = result.array().pow( 2 );
    double J = result.sum();
    J /= ( 2 * m );
    return J; 
}

Eigen::VectorXd LinearRegression::GD( const Eigen::MatrixXd& x, const Eigen::VectorXd& y, 
                   int num_iters  , int batch, 
                   double alpha, double eps){
    
    theta.resize(x.cols());
    theta.setZero();
    p(theta, "Theta:");
    Jhistory.resize(num_iters);
    
    int m = y.size(); //Length of results 
    int n = x.cols(); // Number of features 
    double result;
    Eigen::MatrixXd temp;
    
    Eigen::VectorXd h;
    Eigen::VectorXd g;
    Eigen::VectorXd l;
    //cout << "START COST:" << Cost( x, y, theta ) << " => " ;
    for( int i = 0; i < num_iters; ++i ){
        h = x * theta; 
        l = h - y;
        g= (x.transpose() *l)/m;
        theta = theta - ( alpha * g );
        
        Jhistory[i] = Cost( x, y, theta );
        
        if (Jhistory[i] < eps){
            Jhistory = Jhistory.segment(0,i);
            break;
        }
    }
    CalculateIntercept();
    return Jhistory;
}

Eigen::VectorXd& LinearRegression::CalculateIntercept(){
    if ( mean.size() > 0 ) {
        intercept = CalculateIntercept (mean, theta);
        theta[theta.size()-1] = intercept;
    }
    return theta;
}
 
double LinearRegression::CalculateIntercept(Eigen::VectorXd& mean, Eigen::VectorXd& theta){
    double intercept = mean[mean.rows()-1] - mean.dot(theta);
    return intercept;
}

void LinearRegression::Dump() {
    printf("Linear Regression: Dump() => ");
    if ( Jhistory.size() > 0 )
        p(Jhistory, " Error History");
    p(theta, " Theta");
}


const Eigen::VectorXd LinearRegression::SVD( const Eigen::MatrixXd & x, const Eigen::VectorXd & y){
    const Eigen::VectorXd theta  = x.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(y);
    //CalculateIntercept();
    return theta;
}

const Eigen::VectorXd LinearRegression::QR( const Eigen::MatrixXd & x, const Eigen::VectorXd & y){
    const Eigen::VectorXd theta = x.colPivHouseholderQr().solve(y);
    //CalculateIntercept();
    return theta;
}

const Eigen::VectorXd LinearRegression::Norm( const Eigen::MatrixXd & x, const Eigen::VectorXd & y){
    const Eigen::VectorXd theta  = (x.transpose() * x).ldlt().solve(x.transpose() * y);
    //CalculateIntercept();
    return theta;
}
    
    
/*
Python code:
    import pandas as pd
    from sklearn.linear_model import LinearRegression
    
    df=pd.read_csv("data/test.csv")
    xx = df.values;
    x=xx[:,1:-1]
    y=xx[:,-1]
    arx = LinearRegression(fit_intercept=True).fit(x, y)
    print(arx.coef_, arx.intercept_)
*/
void testLR ( int method = 2) {
    LinearRegression lr = LinearRegression("../data/test1.csv");
    
    Eigen::MatrixXd  x;
    GetMatrix(*lr.csvp, x);
    Eigen::VectorXd& m = lr.Normalize( x, 0 ,0 );
    
    Eigen::VectorXd y = x.col(x.cols() -1);   //Extract last column as y
    SetVal(x, -1, 1); //In our case easy to set last column to 1
    RemoveColumn(x, 0); //first column is not used 
    RemoveColumn(m, 0); // remove it from the coreponding mean value

    switch(method) {
        case 0: lr.GD( x, y, 1024*15, 10, 0.005001); break;
        case 1: lr.SVD( x, y); break;
        case 2: lr.QR( x, y); break;
        default:
        case 3: lr.Norm( x, y); break;
    }
    lr.Dump();
}

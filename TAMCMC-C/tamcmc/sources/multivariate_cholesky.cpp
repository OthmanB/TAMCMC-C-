# include <iostream>
# include <iomanip>
# include <Eigen/Dense>

using namespace std;

using Eigen::VectorXd;
using Eigen::MatrixXd;

VectorXd multivariate_eigen(VectorXd x0, MatrixXd M0, int seed);
double *r8vec_normal_01 ( int n, int *seed );
double *r8vec_uniform_01 ( int n, int *seed );
//****************************************************************************80


VectorXd multivariate_eigen(VectorXd mu, MatrixXd r, int seed){
	
   double *y;

// Do a Cholesky decomposition
   Eigen::MatrixXd Lchol( r.llt().matrixL() );
//
//  Y = MxN matrix of samples of the 1D normal distribution with mean 0
//  and variance 1.  
//
    y = r8vec_normal_01 ( r.rows(), &seed );

/*    std::cout << "Vector of samples y: " << std::endl << Eigen::Map<VectorXd>(y, r.rows()) << std::endl;
    std::cout << "Matrix r : " << std::endl << Lchol.transpose() << std::endl;
    std::cout << " mu : " << std::endl << mu.transpose() << std::endl;

    std::cout << " ----- " << std::endl;
    std::cout << "R'Y from the new code:" << std::endl;
    std::cout << Lchol*Eigen::Map<VectorXd>(y, r.rows()) << std::endl;
*/

return mu + Lchol*Eigen::Map<VectorXd>(y, r.rows());
}

//****************************************************************************80





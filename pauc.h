#ifndef PAUC_H
#define PAUC_H
#include <eigen3/Eigen/Dense>
#include <complex>

class pauli
{
    private:
    Eigen::Matrix2cd emptyp;

    public:
    pauli( int pauli_type_ );
    Eigen::Matrix2cd output() { return emptyp;}
};



#endif

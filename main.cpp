#include <eigen3/Eigen/Dense>
#include <complex>
#include <random>
#include <ctime>
#include <iostream>
#include <fstream>
#include "pauc.h"

#define PI 3.141592

// Creates the 4x4 symmetry matrix from the pauli
// matrices. 
Eigen::Matrix4cd kron( int a, int b, bool anti )
{
	Eigen::Matrix2cd A = pauli(a).output();
	Eigen::Matrix2cd B = pauli(b).output();
	Eigen::Matrix4cd C;
	std::complex<double> I(0.0,1.0);
	for( int i = 0; i < 2; i++ )
	{
		for( int j = 0; j < 2; j++ )
		{
			C.block( i * 2, j * 2, 2, 2 ) = A(i,j) * B; 
		}
	}

	if( anti )
	{
		return I * C;
	}
	else
	{
		return C;
	}
}

// Generates the Hamiltonian
Eigen::Matrix4cd GenHam( std::vector<double> parameters, double p )
{
	Eigen::Matrix4cd Ham;
	std::complex<double> I(0.0,1.0);
	
	std::complex<double> a = parameters[2] + parameters[3] * exp( I * p );
	std::complex<double> b = parameters[4] - parameters[5] * exp( I * p );

	Ham(0,0) = parameters[0];
	Ham(0,1) = std::conj(a);
	Ham(0,2) = 0.0;
	Ham(0,3) = -std::conj(b);
	Ham(1,0) = std::conj(Ham(0,1));
	Ham(1,1) = parameters[1];
	Ham(1,2) = b;
	Ham(1,3) = 0.0;
	Ham(2,0) = std::conj(Ham(0,2));
	Ham(2,1) = std::conj(Ham(1,2));
	Ham(2,2) = -parameters[0];
	Ham(2,3) = -std::conj(a);
	Ham(3,0) = std::conj(Ham(0,3));
	Ham(3,1) = std::conj(Ham(1,3));
	Ham(3,2) = std::conj(Ham(2,3));
	Ham(3,3) = -parameters[1];
	
	return Ham;
}

// Computes the symmetry equation
// C(H(-p)^*C^{\dagger}+-H(p)=\alpha(p)
Eigen::Matrix4cd GenAlpha( std::vector<double> parameters, double p, int am, int bm, bool antiflag, bool PHorTR ) //function that computes the 
{
	
	Eigen::Matrix4cd Ham = GenHam( parameters, p );
	Eigen::Matrix4cd HamM = GenHam( parameters, -p );;
	Eigen::Matrix4cd alpha;
	Eigen::Matrix4cd C = kron(am,bm,antiflag);
	
	if( PHorTR )
	{
		alpha = C * Ham.conjugate() * C.adjoint() + HamM;
	}
	else
	{
		alpha = C * Ham.conjugate() * C.adjoint() - HamM;
	}

	return alpha;
}

// Checks of symmetry equation is obeyed
bool checkSym( Eigen::Matrix4cd alpha )
{
	bool foundflag = true;
	for( int tik1 = 0; tik1 < 4; tik1++ ) {
		for( int tik2 = 0; tik2 < 4; tik2++ ) {
			if( std::abs( alpha(tik1,tik2) ) < 1.0e-8 ) {
				foundflag = foundflag & true;
			}
			else {
				foundflag = foundflag & false;
			}
		}
	}
	return foundflag;
}

// Print to file
void printToFile( bool foundPH, bool foundTR, std::vector<bool> printCond, std::vector<double> parameters, int PHsigAtik, int PHsigBtik, int TRsigAtik, int TRsigBtik, std::ofstream& output_file_SC )
{		
	bool condition = true;
	for( int i = 0; i < 6; i++ ) {
		condition = condition * printCond[i];
	}
	condition = condition * foundPH * foundTR;

	if( condition ) {
		output_file_SC << "mu_a: " << parameters[0] << ", ";
		output_file_SC << "mu_b: " << parameters[1] << ", ";
		output_file_SC << "t_1: " << parameters[2] << ", ";
		output_file_SC << "t_2: " << parameters[3] << ", ";
		output_file_SC << "del_1: " << parameters[4] << ", ";
		output_file_SC << "del_2: " << parameters[5] << ", ";
		output_file_SC << "PH_OP: kron(" << PHsigAtik << "," <<  PHsigBtik << ") ";
		output_file_SC << "TR_OP: kron(" << TRsigAtik << "," << TRsigBtik << ")" << std::endl;			
	}
}

int main()
{
	std::vector<double> parameters(6);
	std::vector<bool> printCond(6);
	Eigen::Matrix4cd alphaPH, alphaTR;
	std::vector<double> values(5);
	bool foundPH, foundTR;

	double p = PI/2.4;
	bool PHantiflag = false, TRantiflag = false;	

	for( int i = 0; i < 5; i++ ) { values[i] = i-2.0; }  

	std::ofstream output_file_SC;
	output_file_SC.open( "symm.txt" );
	
	for( int p0tik = 0; p0tik < 5; p0tik++ ) {
		for( int p1tik = 0; p1tik < 5; p1tik++ ) {
			for( int p2tik = 0; p2tik < 5; p2tik++ ) {
				for( int p3tik = 0; p3tik < 5; p3tik++ ) {
					for( int p4tik = 0; p4tik < 5; p4tik++ ) {
						for( int p5tik = 0; p5tik < 5; p5tik++ ) {
							
							parameters[0] = values[p0tik];
							parameters[1] = values[p1tik];
							parameters[2] = values[p2tik];
							parameters[3] = values[p3tik];
							parameters[4] = values[p4tik];
							parameters[5] = values[p5tik];
							
							printCond[0] = parameters[0] != 0;
							printCond[1] = parameters[1] != 0;
							printCond[2] = parameters[2] != 0;
							printCond[3] = parameters[3] != 0;
							printCond[4] = parameters[4] == 0;
							printCond[5] = parameters[5] == 0;
							
							/*
							printCond[0] = true;
							printCond[1] = true;
							printCond[2] = true;
							printCond[3] = true;
							printCond[4] = true;
							printCond[5] = true;
							*/

							for( int PHsigAtik = 0; PHsigAtik < 4; PHsigAtik++ ) {
								for( int PHsigBtik = 0; PHsigBtik < 4; PHsigBtik++ ) {
									alphaPH = GenAlpha( parameters, p, PHsigAtik, PHsigBtik, PHantiflag, true );
									foundPH = checkSym( alphaPH );	
									if ( foundPH ) {
										for( int TRsigAtik = 0; TRsigAtik < 4; TRsigAtik++ ) {
											for( int TRsigBtik = 0; TRsigBtik < 4; TRsigBtik++ ) {
												alphaTR = GenAlpha( parameters, p, TRsigAtik, TRsigBtik, TRantiflag, false );
												foundTR = checkSym( alphaTR );	
												printToFile( foundPH, foundTR, printCond, parameters, PHsigAtik, PHsigBtik, TRsigAtik, TRsigBtik, output_file_SC );
											}
										}
									}
								}
							}
						
						}
					}
				}
			}
		}
	}
	output_file_SC.close();	
	
	return 0;	
}





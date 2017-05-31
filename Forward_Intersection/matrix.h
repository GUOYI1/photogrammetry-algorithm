#pragma once
#include"common.h"
class Matrix
{
public:
	Matrix();
	Matrix(int, int);
	Matrix(Image &c1);
	Matrix(Coordinate&);
	friend Matrix operator+(Matrix &m1, Matrix &m2);
	friend Matrix operator*(Matrix &m1, Matrix &m2);
	friend Matrix operator*(double m, Matrix &m2);
	friend bool operator==(Matrix m1, Matrix m2);
	friend std::ostream& operator <<(std::ostream&, Matrix&);
	friend void Copy_Matrix(Matrix &m1, Matrix m2);
	Matrix Transpose();
	Matrix Inverse();
	double Del();
	double Cofactor(int, int);
	Matrix Submatrix(int, int);
	Matrix Adjoint();
	Matrix CrossStitching(Matrix m);
	Matrix VerticalStitching(Matrix m);
	Matrix* Subdeterminant(int order);
	int Rank();
	friend class Image;
private:
	double **data;
	int row;
	int arrange;
};

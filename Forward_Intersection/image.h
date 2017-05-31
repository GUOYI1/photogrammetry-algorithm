#pragma once
#include"common.h"
class Image
{
public:
	Image(int);
	static void parameter();
	void GetXY();
	void Resection();									
	static void ForwardIntersection(Image&, Image&);	
	void getphotoCo(Matrix &);
	void Derivatives(Matrix &, Matrix &);				//Form the partial derivatives matrix
	int Correlation(Matrix &, Matrix &, Matrix &);		//Error equation correction
	friend class Matrix;
private:
	double XS;
	double YS;
	double ZS;
	double phi;
	double omega;
	double kappa;
	double Co[4][3];
	double photorealCo[4][2];
	double photoCo[4][2];
	double X_[4];
	double Y_[4];
	double Z_[4];
	static double f;											//focal length
	static double Tolerance;
	double Iterations;
	Coordinate Coo;
};

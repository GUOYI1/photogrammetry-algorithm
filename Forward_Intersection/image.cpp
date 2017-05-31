#include"image.h"
#include"matrix.h"
using namespace std;
Image::Image(int a)
{
	Co[0][0] = 0;
	Co[0][1] = 0;
	Co[0][2] = 0;
	Co[1][0] = 1 * 1000;
	Co[1][1] = 0;
	Co[1][2] = 0;
	Co[2][0] = 0;
	Co[2][1] = 1 * 1000;
	Co[2][2] = 0;
	Co[3][0] = 1 * 1000;
	Co[3][1] = 1 * 1000;
	Co[3][2] = 0;
	XS = 0;
	YS = 0;
	ZS = 7 * 1000;
	//As an example, We set the control points like this.
	//IMG_LEFT: A(2783£¬1758),B(3116£¬1760),C(2778£¬2092),D(3116£¬2096)
	//IMG_RIGHT(Corresponding Points): A'(2477£¬1768),B'(2809£¬1763),C'(2474£¬2102),D'(2811£¬2100)
	if (a == 1)
	{
		photorealCo[0][0] = (2783 - 2880)*pixel;
		photorealCo[0][1] = (1758 - 1920)*pixel;
		photorealCo[1][0] = (3116 - 2880)*pixel;
		photorealCo[1][1] = (1760 - 1920)*pixel;
		photorealCo[2][0] = (2778 - 2880)*pixel;
		photorealCo[2][1] = (2092 - 1920)*pixel;
		photorealCo[3][0] = (3116 - 2880)*pixel;
		photorealCo[3][1] = (2096 - 1920)*pixel;
	}
	else
	{
		photorealCo[0][0] = (2477 - 2880)*pixel;
		photorealCo[0][1] = (1768 - 1920)*pixel;
		photorealCo[1][0] = (2809 - 2880)*pixel;
		photorealCo[1][1] = (1763 - 1920)*pixel;
		photorealCo[2][0] = (2474 - 2880)*pixel;
		photorealCo[2][1] = (2102 - 1920)*pixel;
		photorealCo[3][0] = (2811 - 2880)*pixel;
		photorealCo[3][1] = (2100 - 1920)*pixel;
	}
	for (int i = 0; i < 4; i++)
	{
		XS = XS + (Co[i][0] / 4);
		YS = YS + (Co[i][1] / 4);
		ZS = ZS + (Co[i][2] / 4);
	}
	phi = 0;
	omega = 0;
	kappa = 0;
	Iterations = 1;
}
double Image::f = 0;
double Image::Tolerance = 0;
void Image::parameter()
{
	cout << "Input Focal Length(mm£©£º";
	cin >> f;
	cout << endl;
	cout << "Input Tolerance:";
	cin >> Tolerance;
	cout << endl;
}
void Image::GetXY()
{
	system("start readpixel.exe");
	system("Pause");
	ifstream infile;
	infile.open("data.dat", ios::in);
	infile >> Coo.x >> Coo.y;
	//Here we assume that the default size for an image is 5616*3744
	Coo.x = (Coo.x - 2880)*pixel;
	Coo.y = (Coo.y - 1920)*pixel;
	infile.close();
}
void Image::Resection()
{
	int x;
	do
	{
		Matrix m(*this), A(8, 6), L(8, 1), X(6, 1);
		getphotoCo(m);
		Derivatives(A, m);
		x = Correlation(A, L, X);
	} while (x == 0);
}
void Image::getphotoCo(Matrix &m)
{

	for (int i = 0; i < 4; i++)
	{
		X_[i] = m.data[0][0] * (Co[i][0] - XS) + m.data[1][0] * (Co[i][1] - YS) + m.data[2][0] * (Co[i][2] - ZS);
		Y_[i] = m.data[0][1] * (Co[i][0] - XS) + m.data[1][1] * (Co[i][1] - YS) + m.data[2][1] * (Co[i][2] - ZS);
		Z_[i] = m.data[0][2] * (Co[i][0] - XS) + m.data[1][2] * (Co[i][1] - YS) + m.data[2][2] * (Co[i][2] - ZS);
	}

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 2; j++)
		{
			if (j == 0)
				photoCo[i][j] = -f*X_[i] / Z_[i];
			else
				photoCo[i][j] = -f*Y_[i] / Z_[i];
		}
}
void Image::Derivatives(Matrix &A, Matrix &m)
{
	for (int i = 0; i < A.row; i++)
	{
		if (i % 2 == 0)
		{
			A.data[i][0] = (m.data[0][0] * f + m.data[0][2] * photorealCo[i / 2][0]) / Z_[i / 2];
			A.data[i][1] = (m.data[1][0] * f + m.data[1][2] * photorealCo[i / 2][0]) / Z_[i / 2];
			A.data[i][2] = (m.data[2][0] * f + m.data[2][2] * photorealCo[i / 2][0]) / Z_[i / 2];
			A.data[i][3] = photorealCo[i / 2][1] * sin(omega) - (photorealCo[i / 2][0] / f*(photorealCo[i / 2][0] * cos(kappa) - photorealCo[i / 2][1] * sin(kappa)) + f*cos(kappa))*cos(omega);
			A.data[i][4] = -f*sin(kappa) - photorealCo[i / 2][0] / f*(photorealCo[i / 2][0] * sin(kappa) + photorealCo[i / 2][1] * cos(kappa));
			A.data[i][5] = photorealCo[i / 2][1];
		}
		else
		{
			int j;
			j = (i - 1) / 2;
			A.data[i][0] = (m.data[0][1] * f + m.data[0][2] * photorealCo[j][1]) / Z_[j];
			A.data[i][1] = (m.data[1][1] * f + m.data[1][2] * photorealCo[j][1]) / Z_[j];
			A.data[i][2] = (m.data[2][1] * f + m.data[2][2] * photorealCo[j][1]) / Z_[j];
			A.data[i][3] = -photorealCo[j][0] * sin(omega) - (photorealCo[j][1] / f*(photorealCo[j][0] * cos(kappa) - photorealCo[j][1] * sin(kappa)) - f*sin(kappa))*cos(omega);
			A.data[i][4] = -f*cos(kappa) - photorealCo[j][1] / f*(photorealCo[j][0] * sin(kappa) + photorealCo[j][1] * cos(kappa));
			A.data[i][5] = -photorealCo[j][0];
		}
	}
}
int Image::Correlation(Matrix &A, Matrix &L, Matrix &X)
{
	int i, j;
	for (i = 0; i < 8; i = i + 2)
	{
		j = i / 2;
		L.data[i][0] = photorealCo[j][0] - photoCo[j][0];
		L.data[i + 1][0] = photorealCo[j][1] - photoCo[j][1];
	}
	X = ((A.Transpose()*A).Inverse())*(A.Transpose())*L;
	for (i = 0; i < 6; i++)
	{
		if (fabs(X.data[i][0])>Tolerance)
			break;
	}
	if (i == 6)
	{
		cout << "Iterations:" << Iterations << endl;
		cout << "Central point XS:" << XS << endl;
		cout << "Central point YS:" << YS << endl;
		cout << "Central point ZS:" << ZS << endl;
		cout << "Rotation angle along the x axis(Unit: degree):" << phi * 180 / PI << endl;
		cout << "Rotation angle along the y axis(Unit: degree):" << omega * 180 / PI << endl;
		cout << "Rotation angle along the z axis(Unit: degree):" << kappa * 180 / PI << endl;
		return(1);
	}
	else
	{
		XS = XS + X.data[0][0];
		YS = YS + X.data[1][0];
		ZS = ZS + X.data[2][0];
		phi = phi + X.data[3][0];
		omega = omega + X.data[4][0];
		kappa = kappa + X.data[5][0];
		Iterations++;
		return(0);
	}
}
void Image::ForwardIntersection(Image& p1, Image& p2)
{
	double N1, N2, Bx, By, Bz;
	Matrix R1(p1), R2(p2), PS1(3, 1), PS2(3, 1), A1(3, 1), A2(3, 1), M1(3, 1), M2(3, 1), m1(p1.Coo), m2(p2.Coo);
	Bx = p2.XS - p1.XS;							//Computer the photogrammetric base, Bx, By, Bz
	By = p2.YS - p1.YS;
	Bz = p2.ZS - p1.ZS;
	PS1.data[0][0] = p1.XS;
	PS1.data[1][0] = p1.YS;
	PS1.data[2][0] = p1.ZS;
	PS2.data[0][0] = p2.XS;
	PS2.data[1][0] = p2.YS;
	PS2.data[2][0] = p2.ZS;
	M1 = R1*m1;									//Get the image space auxiliary coordinates. Save them in M1 and M2
	M2 = R2*m2;
	cout << "Matrix" << endl;
	cout << M1 << endl;
	cout << M2 << endl;
	N1 = (Bx*M2.data[2][0] - Bz * M2.data[0][0]) / (M1.data[0][0] * M2.data[2][0] - M2.data[0][0] * M1.data[2][0]);
	N2 = (Bx*M1.data[2][0] - Bz * M1.data[0][0]) / (M1.data[0][0] * M2.data[2][0] - M2.data[0][0] * M1.data[2][0]);
	A1 = PS1 + N1*M1;
	A2 = PS2 + N2*M2;
	cout << "Projection coefficient for Image1:" << N1 << endl;
	cout << "The real coordinate in Image1:" << endl;
	cout << 0.001*A1 << endl;
	cout << "Projection coefficient for Image2:" << N2 << endl;
	cout << "The real coordinate in Image2:" << endl;
	cout << 0.001*A2 << endl;
}
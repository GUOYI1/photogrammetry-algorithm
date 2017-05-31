#include<iostream >
#include<cstdlib>
#include<cmath>
#include<string>
#include<fstream >
#define pixel 0.006410256
#define PI 3.14159265358
using namespace std;
struct Coordinates
{
	double x;
	double y;
};
class Matrix;
class Calculation
{
public:
	Calculation();
	double Length(Matrix &m1);
	double Width(Matrix &m1);
	double Height(Matrix &m1);
	friend class Matrix;
private:
	Coordinates Co[4];
	double f;
	double H;
	double a;
	double w;
	double k;
};
class Matrix													//Matrix
{
public:
	Matrix(int, int);
	Matrix(Calculation &c1);
	Matrix(double, double, Calculation &c1);
	friend Matrix operator*(Matrix &m1, Matrix &m2);
	friend Matrix operator*(double m, Matrix &m2);
	friend ostream& operator <<(ostream&, Matrix&);
	friend double Calculation::Length(Matrix &m1);
	friend double Calculation::Width(Matrix &m1);
	friend double Calculation::Height(Matrix &m1);
private:
	double **temp;
	int row;
	int arrange;

};
Calculation::Calculation()
{
	int i;
	ifstream infile;
	infile.open("temp.dat", ios::in);
	for (i = 0; i<4; i++)
	{
		infile >> Co[i].x >> Co[i].y;
	}
	infile.close();
	cout << "Input focal distance£¨Unit: mm£©:";
	cin >> f;
	cout << endl;
	cout << "Input the distance between lens and ground£¨Unit: m£©:";
	cin >> H;
	cout << endl;
	cout << "Input deflection angle a£¨Unit: degree£©:";
	cin >> a;
	a = a / 360 * PI * 2;
	cout << endl;
	cout << "Input deflection angle w£¨Unit: degree£©:";
	cin >> w;
	w = w / 360 * PI * 2;
	cout << endl;
	cout << "Input deflection angle k£¨Unit: degree£©:";
	cin >> k;
	k = k / 360 * PI * 2;
	cout << endl;
}
double Calculation::Length(Matrix &m1)
{
	double p, X, Y, X1, Y1, SY;
	X = Co[0].x*pixel;
	Y = Co[0].y*pixel;
	Matrix m2(18, 12, *this);
	Matrix m((m1*m2).row, (m1*m2).arrange);
	X1 = -f*(m1.temp[0][0] * X + m1.temp[0][1] * Y - f*m1.temp[0][2]) / (m1.temp[2][0] * X + m1.temp[2][1] * Y - f*m1.temp[2][2]);
	Y1 = -f*(m1.temp[1][0] * X + m1.temp[1][1] * Y - f*m1.temp[1][2]) / (m1.temp[2][0] * X + m1.temp[2][1] * Y - f*m1.temp[2][2]);
	SY = -f*(m1.temp[1][0] * 18 + m1.temp[1][1] * 12 - f*m1.temp[1][2]) / (m1.temp[2][0] * 18 + m1.temp[2][1] * 12 - f*m1.temp[2][2]);
	p = H * 1000 / (Y1 - SY);
	m = (p / 1000)*m1*m2;
	if (m.temp[0][2]<0)
		m.temp[0][2] = -m.temp[0][2];
	return(m.temp[0][2]);
};
double Calculation::Width(Matrix &m1)
{
	double p, Xa, Ya, Xb, Yb, Xb1, Yb1, SY;
	Xa = Co[1].x*pixel;
	Ya = Co[1].y*pixel;
	Xb = Co[2].x*pixel;
	Yb = Co[2].y*pixel;
	Matrix m4(Xa, Ya, *this);
	Matrix m5(Xb, Yb, *this);
	Matrix m2((m1*m4).row, (m1*m4).arrange), m3((m1*m5).row, (m1*m5).arrange);
	Xb1 = -f*(m1.temp[0][0] * Xb + m1.temp[0][1] * Yb - f*m1.temp[0][2]) / (m1.temp[2][0] * Xb + m1.temp[2][1] * Yb - f*m1.temp[2][2]);
	Yb1 = -f*(m1.temp[1][0] * Xb + m1.temp[1][1] * Yb - f*m1.temp[1][2]) / (m1.temp[2][0] * Xb + m1.temp[2][1] * Yb - f*m1.temp[2][2]);
	SY = -f*(m1.temp[1][0] * 18 + m1.temp[1][1] * 12 - f*m1.temp[1][2]) / (m1.temp[2][0] * 18 + m1.temp[2][1] * 12 - f*m1.temp[2][2]);
	p = H * 1000 / (Yb1 - SY);
	m2 = (p / 1000)*m1*m4;
	m3 = (p / 1000)*m1*m5;
	return(m3.temp[0][0] - m2.temp[0][0]);
};
double Calculation::Height(Matrix &m1)
{
	double p, Xa, Ya, Xb, Yb, Xa1, Ya1, SY;

	Xa = Co[1].x*pixel;
	Ya = Co[1].y*pixel;
	Xb = Co[3].x*pixel;
	Yb = Co[3].y*pixel;
	Matrix m4(Xa, Ya, *this);
	Matrix m5(Xb, Yb, *this);
	Matrix m2((m1*m4).row, (m1*m4).arrange), m3((m1*m5).row, (m1*m5).arrange);
	Xa1 = -f*(m1.temp[0][0] * Xa + m1.temp[0][1] * Ya - f*m1.temp[0][2]) / (m1.temp[2][0] * Xa + m1.temp[2][1] * Ya - f*m1.temp[2][2]);
	Ya1 = -f*(m1.temp[1][0] * Xa + m1.temp[1][1] * Ya - f*m1.temp[1][2]) / (m1.temp[2][0] * Xa + m1.temp[2][1] * Ya - f*m1.temp[2][2]);
	SY = -f*(m1.temp[1][0] * 18 + m1.temp[1][1] * 12 - f*m1.temp[1][2]) / (m1.temp[2][0] * 18 + m1.temp[2][1] * 12 - f*m1.temp[2][2]);
	p = H * 1000 / (Ya1 - SY);
	m2 = (p / 1000)*m1*m4;
	m3 = (p / 1000)*m1*m5;
	return(m2.temp[0][1] - m3.temp[0][1]);
};
Matrix::Matrix(int r, int a)
{
	row = r;
	arrange = a;
	temp = new double*[row];
	temp[0] = new double[row*arrange];
	for (int i = 1; i<row; i++)
		temp[i] = temp[i - 1] + arrange;
}
Matrix::Matrix(Calculation &c1)										//Initialization
{
	row = 3;
	arrange = 3;
	temp = new double*[row];
	temp[0] = new double[row*arrange];
	for (int i = 1; i<row; i++)
		temp[i] = temp[i - 1] + arrange;
	temp[0][0] = cos(c1.a)*cos(c1.k) - sin(c1.a)*sin(c1.w)*sin(c1.k);
	temp[0][1] = -cos(c1.a)*sin(c1.k) - sin(c1.a)*sin(c1.w)*cos(c1.k);
	temp[0][2] = -sin(c1.a)*cos(c1.w);
	temp[1][0] = cos(c1.w)*sin(c1.k);
	temp[1][1] = cos(c1.w)*cos(c1.k);
	temp[1][2] = -sin(c1.w);
	temp[2][0] = sin(c1.a)*cos(c1.k) + cos(c1.a)*sin(c1.w)*sin(c1.k);
	temp[2][1] = -sin(c1.a)*sin(c1.k) + cos(c1.a)*sin(c1.w)*cos(c1.k);
	temp[2][2] = cos(c1.a)*cos(c1.w);
}
Matrix::Matrix(double x, double y, Calculation &c1)
{
	row = 3;
	arrange = 1;
	temp = new double*[row];
	temp[0] = new double[row*arrange];
	for (int i = 1; i<row; i++)
		temp[i] = temp[i - 1] + arrange;
	temp[0][0] = x;
	temp[0][1] = y;
	temp[0][2] = -c1.f;
}
ostream& operator <<(ostream &output, Matrix &m)					//Matrix Output
{
	for (int i = 0; i<m.row; i++)
	{
		for (int j = 0; j<m.arrange; j++)
			cout << m.temp[i][j] << " ";
		cout << endl;
	}
	return(output);
}
Matrix operator*(double m, Matrix &m2)
{
	Matrix m1(m2.row, m2.arrange);
	for (int i = 0; i<m2.row; i++)
	for (int j = 0; j<m2.arrange; j++)
	{
		m1.temp[i][j] = m*m2.temp[i][j];
	}
	return(m1);
}
Matrix operator*(Matrix &m1, Matrix &m2)							//Matrix Multiplication
{
	Matrix m(m1.row, m2.arrange);
	if (m1.arrange != m2.row)
	{
		cout << "Math Error" << endl;
		exit(0);
	}

	for (int i = 0; i<m1.row; i++)
	for (int j = 0; j<m2.arrange; j++)
	{
		m.temp[i][j] = 0;
		for (int k = 0; k<m1.arrange; k++)
			m.temp[i][j] = m.temp[i][j] + m1.temp[i][k] * m2.temp[k][j];
	}
	return(m);
}
int main()
{
	system("start readpixels.exe");
	system("Pause");
	Calculation c1;
	Matrix m1(c1);
	cout << "Length of the corridor:	" << c1.Length(m1) << endl;
	cout << "Width of the corridor:		" << c1.Width(m1) << endl;
	cout << "Height of the corridor:	" << c1.Height(m1) << endl;
	system("Pause");
	return(0);
}

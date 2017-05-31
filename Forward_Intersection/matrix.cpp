#include"matrix.h"
#include"image.h"
using namespace std;
Matrix::Matrix() :data(NULL), row(0), arrange(0)
{}
Matrix::Matrix(int r, int a)
{
	row = r;
	arrange = a;
	data = new double*[row];
	data[0] = new double[row*arrange];
	for (int i = 1; i<row; i++)
		data[i] = data[i - 1] + arrange;
}
Matrix::Matrix(Image& c1)
{
	row = 3;
	arrange = 3;
	data = new double*[row];
	data[0] = new double[row*arrange];
	for (int i = 1; i<row; i++)
		data[i] = data[i - 1] + arrange;
	data[0][0] = cos(c1.phi)*cos(c1.kappa) - sin(c1.phi)*sin(c1.omega)*sin(c1.kappa);
	data[0][1] = -cos(c1.phi)*sin(c1.kappa) - sin(c1.phi)*sin(c1.omega)*cos(c1.kappa);
	data[0][2] = -sin(c1.phi)*cos(c1.omega);
	data[1][0] = cos(c1.omega)*sin(c1.kappa);
	data[1][1] = cos(c1.omega)*cos(c1.kappa);
	data[1][2] = -sin(c1.omega);
	data[2][0] = sin(c1.phi)*cos(c1.kappa) + cos(c1.phi)*sin(c1.omega)*sin(c1.kappa);
	data[2][1] = -sin(c1.phi)*sin(c1.kappa) + cos(c1.phi)*sin(c1.omega)*cos(c1.kappa);
	data[2][2] = cos(c1.phi)*cos(c1.omega);
}
Matrix::Matrix(Coordinate& Coo)
{
	row = 3;
	arrange = 1;
	data = new double*[row];
	data[0] = new double[row*arrange];
	for (int i = 1; i<row; i++)
		data[i] = data[i - 1] + arrange;
	data[0][0] = Coo.x;
	data[1][0] = Coo.y;
	data[2][0] = -(Image::f);
}
ostream& operator <<(ostream &output, Matrix &m)
{
	for (int i = 0; i<m.row; i++)
	{
		for (int j = 0; j<m.arrange; j++)
			cout << m.data[i][j] << " ";
		cout << endl;
	}
	return(output);
}
Matrix operator+(Matrix &m1, Matrix &m2)
{
	Matrix m(m1.row, m1.arrange);
	if (m1.arrange != m2.arrange || m1.row != m2.row)
	{
		cout << "Math Error" << endl;
		exit(0);
	}
	for (int i = 0; i<m1.row; i++)
		for (int j = 0; j<m1.arrange; j++)
			m.data[i][j] = m1.data[i][j] + m2.data[i][j];
	return(m);
}
Matrix operator*(Matrix &m1, Matrix &m2)
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
			m.data[i][j] = 0;
			for (int k = 0; k<m1.arrange; k++)
				m.data[i][j] = m.data[i][j] + m1.data[i][k] * m2.data[k][j];
		}
	return(m);
}
Matrix operator*(double m, Matrix &m2)
{
	Matrix m1(m2.row, m2.arrange);
	for (int i = 0; i<m2.row; i++)
		for (int j = 0; j<m2.arrange; j++)
		{
			m1.data[i][j] = m*m2.data[i][j];
		}
	return(m1);
}
bool operator==(Matrix m1, Matrix m2)
{
	if (m1.row != m2.row || m1.arrange != m2.arrange)
		return false;
	else
	{
		int i = 0, j = 0;
		for (i = 0; i < m1.row; i++)
		{
			for (j = 0; j < m1.arrange; j++)
				if (m1.data[i][j] != m2.data[i][j])
					break;
			if (j != m1.arrange)break;
		}
		if (i == m1.row) return(true);
		else return false;
	}
}
Matrix Matrix::Transpose()
{
	Matrix m(arrange, row);
	for (int i = 0; i < m.row; i++)
		for (int j = 0; j < m.arrange; j++)
			m.data[i][j] = data[j][i];
	return(m);
}
Matrix Matrix::Submatrix(int r, int a)											
{
	r--; a--;
	Matrix m(row - 1, arrange - 1);
	if (row != arrange)
		cout << "Math Error" << endl;
	else
	{
		for (int i = 0; i < m.row; i++)
			for (int j = 0; j < m.arrange; j++)
			{
				if (i < r&&j < a)
					m.data[i][j] = data[i][j];
				else
				{
					if (i >= r&&j < a)
						m.data[i][j] = data[i + 1][j];
					else
					{
						if (i < r&&j >= a)
							m.data[i][j] = data[i][j + 1];
						else
							m.data[i][j] = data[i + 1][j + 1];
					}
				}
			}
	}
	return(m);
}
Matrix Matrix::Adjoint()														
{
	Matrix m(row, arrange);
	for (int i = 0; i < row; i++)
		for (int j = 0; j < arrange; j++)
			m.data[i][j] = Cofactor(i + 1, j + 1);
	m = m.Transpose();
	return (m);
}
Matrix Matrix::Inverse()
{
	Matrix R(row, arrange);
	Matrix E(row, arrange);
	for (int i = 0; i < E.row; i++)
	{
		for (int j = 0; j < E.arrange; j++)
		{
			if (i == j)
				E.data[i][j] = 1;
			else
				E.data[i][j] = 0;
		}
	}
	for (int k = 0; k < E.row - 1; k++)
	{
		for (int i = k + 1; i < E.row; i++)
		{
			double m = data[i][k] / data[k][k];
			for (int j = 0; j < E.arrange; j++)
			{
				data[i][j] = data[i][j] - m*data[k][j];
				E.data[i][j] = E.data[i][j] - m*E.data[k][j];
			}
		}
	}
	for (int k = E.row - 1; k > 0; k--)
	{
		for (int i = k - 1; i >= 0; i--)
		{
			double m = data[i][k] / data[k][k];
			for (int j = 0; j < arrange; j++)
			{
				data[i][j] = data[i][j] - m*data[k][j];
				E.data[i][j] = E.data[i][j] - m*E.data[k][j];
			}
		}
	}
	for (int i = 0; i < E.row; i++)
	{
		for (int j = 0; j < E.arrange; j++)
		{
			R.data[i][j] = E.data[i][j] / data[i][i];
		}
	}
	return R;
}

Matrix Matrix::CrossStitching(Matrix m)									
{
	Matrix m1(row, arrange + m.arrange);
	if (row != m.row)cout << "Stitichin Error" << endl;
	else
	{
		for (int i = 0; i < m1.row; i++)
		{
			for (int j = 0; j < m1.arrange; j++)
			{
				if (j < arrange)
					m1.data[i][j] = data[i][j];
				else
					m1.data[i][j] = m.data[i][j - arrange];
			}
		}
	}
	return(m1);
}
Matrix Matrix::VerticalStitching(Matrix m)									
{
	Matrix m1(row + m.row, arrange);
	if (arrange != m.arrange)cout << "Stitichin Error" << endl;
	else
	{
		for (int i = 0; i < m1.arrange; i++)
		{
			for (int j = 0; j < m1.row; j++)
			{
				if (j < row)
					m1.data[j][i] = data[j][i];
				else
					m1.data[j][i] = m.data[j - row][i];
			}
		}
	}
	return(m1);
}
double Matrix::Del()													
{
	double del = 0;
	if (row != arrange)
		cout << "Math Error" << endl;
	else
	{
		if (row == 1)return data[0][0];
		if (row == 2)
			return(data[0][0] * data[1][1] - data[0][1] * data[1][0]);
		else
		{
			for (int k = 0; k < arrange; k++)
			{
				Matrix m(row - 1, arrange - 1);
				m = Submatrix(1, k + 1);
				if (k % 2 != 0)
					del = del + (-1)*data[0][k] * m.Del();
				else
					del = del + data[0][k] * m.Del();
			}
		}
	}
	return(del);
}
double Matrix::Cofactor(int r, int a)										
{
	Matrix m(row - 1, arrange - 1);
	if (row != arrange)
		cout << "Math Error" << endl;
	else
	{
		m = Submatrix(r, a);
		r--; a--;
		if ((r + a) % 2 != 0)
			return ((-1) * m.Del());
		else
			return(m.Del());
	}
}

Matrix* Matrix::Subdeterminant(int order)
{
	int array[2] = { row,arrange };
	if (order > Min(array, 2))
		cout << order << "th subdeterminant does not exist" << endl;
	int num = Combine_Num(row, order)*Combine_Num(arrange, order);
	int** row_select = Combine_Result(row, order);
	int** arrange_select = Combine_Result(arrange, order);
	
	int pos = 0;
	Matrix* result = new Matrix[num];
	for (int i = 0; i < num; i++)
		result[i] = Matrix(order, order);


	for (int row_select_pos = 0; row_select_pos < Combine_Num(row, order); row_select_pos++)
		for (int arrange_select_pos = 0; arrange_select_pos < Combine_Num(arrange, order); arrange_select_pos++)
		{
			for (int i = 0; i < order; i++)
				for (int j = 0; j < order; j++)
				{
					result[pos].data[i][j] = data[row_select[row_select_pos][i]]\
						[arrange_select[arrange_select_pos][j]];

				}
			pos = pos++;
		}

	delete row_select[0];
	delete arrange_select[0];
	delete row_select;
	delete arrange_select;
	return(result);
}

int Matrix::Rank()
{
	int array[2] = { row,arrange };
	Matrix m1(row, arrange);
	if (*this == m1)
		return 0;
	int i = 0;
	for (i = Min(array, 2); i >0; i--)
	{
		Matrix* data = Subdeterminant(i);
		int j = 0;
		for (j = 0; j < Combine_Num(row, i)*Combine_Num(arrange, i); j++)
			if (data[j].Del() != 0) break;
		delete data;
		if (j != Combine_Num(row, i)*Combine_Num(arrange, i))
			break;
	}
	return(i);
}
void Copy_Matrix(Matrix &m1, Matrix m2)
{
	m1.row = m2.row;
	m1.arrange = m2.arrange;
	m1.data = new double*[m1.row];
	m1.data[0] = new double[m1.row*m1.arrange];
	for (int i = 1; i < m1.row; i++)
		m1.data[i] = m1.data[i - 1] + m1.arrange;
	for (int i = 0; i < m1.row; i++)
		for (int j = 0; j < m1.arrange; j++)
			m1.data[i][j] = m2.data[i][j];
}
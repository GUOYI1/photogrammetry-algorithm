#pragma once
#include<fstream >
#include<cstdlib>
#include<iostream>
#include<math.h>
#define No_Repeat 0
#define Repeat 1
#define PI 3.14159265358
#define pixel 0.00625

struct Coordinate
{
	double x;
	double y;
	double z;
};
class Image;
class Matrix;

template<typename GY_TYPE>
inline GY_TYPE Min(GY_TYPE array[], int length)
{
	GY_TYPE min = array[0];
	for (int i = 0; i < length; i++)
		if (abs(array[i] - min)>1e-6)
			if (array[i] < min) min = array[i];
	return(min);
}

template<typename GY_TYPE>
inline GY_TYPE Max(GY_TYPE array[], int length)
{
	GY_TYPE max = array[0];
	for (int i = 0; i < int length; i++)
		if (abs(array[i] - min)>1e-6)
			if (array[i] >min) max = array[i];
	return(max);
}

template<typename GY_TYPE>
inline bool Array_in(GY_TYPE data, GY_TYPE* array, int length)
{
	int i = 0
		for (i = 0; i < length; i++)
			if (abs(data - array[i])<1e-6) break;
	if (i == length)return false;
	else return true;

}

inline int Factorial(int num)
{
	int Product = 1;
	for (int i = num; i >0; i--)
	{
		Product = Product*i;
	}
	return Product;
}

inline int Array_Num(int ele_num, int select_num)                     //Permutation
{
	return(Factorial(ele_num) / Factorial(ele_num - select_num));
}

inline int Combine_Num(int ele_num, int select_num)                   //组合
{
	return(Array_Num(ele_num, select_num) / Array_Num(select_num, select_num));
}

inline void Array_Iteration(int**result, int* data, int dimension, int &num, int ele_num, int select_num, bool repeat_flag)//排列迭代
{
	if (repeat_flag)
	{
		for (int i = 0; i < ele_num; i++)
		{
			data[dimension - 1] = i;
			if (dimension != select_num)
				Array_Iteration(result, data, dimension + 1, num, ele_num, select_num, repeat_flag);
			else
			{
				for (int j = 0; j < dimension; j++)
					result[num][j] = data[j];
				num++;
			}
		}
	}
	else
	{
		for (int i = 0; i < ele_num; i++)
		{
			if (dimension != 1)
			{
				int j = 0;
				for (j = 0; j < dimension - 1; j++)
					if (data[j] == i)
						break;
				if (j != dimension - 1) continue;
			}
			data[dimension - 1] = i;
			if (dimension != select_num)
				Array_Iteration(result, data, dimension + 1, num, ele_num, select_num, repeat_flag);
			else
			{
				for (int j = 0; j < select_num; j++)
					result[num][j] = data[j];
				num++;
			}

		}
	}
}
inline int** Array_Result(int ele_num, int select_num, bool repeat_flag)                                               //排列结果
{
	int num = 0;
	int* data = new int[select_num];
	int** result = NULL;
	memset(data, 0, select_num*sizeof(int));
	if (repeat_flag)
	{
		result = new int*[pow(ele_num, select_num)];
		result[0] = new int[pow(ele_num, select_num)*select_num];
		memset(result[0], 0, pow(ele_num, select_num)*select_num*sizeof(int));
		for (int i = 1; i < pow(ele_num, select_num); i++)
			result[i] = result[i - 1] + select_num;
		Array_Iteration(result, data, 1, num, ele_num, select_num, repeat_flag);

	}
	else
	{
		result = new int*[Array_Num(ele_num, select_num)];
		result[0] = new int[Array_Num(ele_num, select_num)*select_num];
		memset(result[0], 0, Array_Num(ele_num, select_num)*sizeof(int));
		for (int i = 1; i < Array_Num(ele_num, select_num); i++)
			result[i] = result[i - 1] + select_num;
		Array_Iteration(result, data, 1, num, ele_num, select_num, repeat_flag);

	}
	delete data;
	return(result);

}

inline void Combine_Iteration(int**result, int* data, int dimension, int &num, int ele_num, int select_num, int start)//排列迭代
{
	for (int i = start; i < ele_num; i++)
	{
		data[dimension - 1] = i;
		if (dimension != select_num)
			Combine_Iteration(result, data, dimension + 1, num, ele_num, select_num, i + 1);
		else
		{
			for (int j = 0; j < select_num; j++)
				result[num][j] = data[j];
			num++;
		}
	}
}

inline int** Combine_Result(int ele_num, int select_num)
{
	int num = 0;
	int* data = new int[select_num];
	int** result = NULL;
	memset(data, 0, select_num*sizeof(int));
	result = new int*[Combine_Num(ele_num, select_num)];
	result[0] = new int[Combine_Num(ele_num, select_num)*select_num];
	memset(result[0], 0, Combine_Num(ele_num, select_num)*sizeof(int));
	for (int i = 1; i < Combine_Num(ele_num, select_num); i++)
		result[i] = result[i - 1] + select_num;
	Combine_Iteration(result, data, 1, num, ele_num, select_num, 0);
	delete data;
	return(result);

}
inline double* least_square_method(double(*data)[2], int n)
{
	double sum_x = 0, sum_y = 0, avg_x = 0, avg_y = 0, \
		Sxy = 0, Sx = 0, Sy = 0, ceigma = 0;
	double* result = new double[4];
	for (int i = 0; i < n; i++)
	{
		sum_x += data[i][0];
		sum_y += data[i][1];
	}
	avg_x = sum_x / n;
	avg_y = sum_y / n;
	for (int i = 0; i < n; i++)
	{
		Sx += pow(data[i][0] - avg_x, 2);
		Sy += pow(data[i][1] - avg_y, 2);
		Sxy += (data[i][0] - avg_x)*(data[i][1] - avg_y);
	}
	Sx = sqrt(Sx / (n - 1));
	Sy = sqrt(Sy / (n - 1));
	Sxy = Sxy / (n - 1);
	result[0] = Sxy / (Sx*Sy);
	result[1] = Sxy / pow(Sx, 2);
	result[2] = avg_y - result[1] * avg_x;
	for (int i = 0; i < n; i++)
	{
		double y_i = data[i][0] * result[1] + result[2];
		ceigma += pow(data[i][1] - y_i, 2);
	}
	ceigma = sqrt(ceigma);
	result[3] = ceigma;
	return(result);
}

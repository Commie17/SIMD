#include <immintrin.h>
#include <iostream>
#include <chrono>
#include <random>

using namespace std;
int len = 4096;
//дл€ доступа к функци€м умножени€ объ€вим их перед функцией main
void Cmp(__int64** matrix_v, __int64** matrix_c);
void Scl(__int64** matrix_a, __int64** matrix_b,__int64** matrix_c);
void Vec(__int64** matrix_a, __int64** matrix_b, __int64** matrix_c);

int main() 
{
	//создаЄм матрицы
	__int64** Matrix_vec= new __int64* [len];
	__int64** Matrix_sc=new __int64* [len];
	__int64** Matrix_a=new __int64* [len];
	__int64** Matrix_b= new __int64* [len];

	//заполн€ем матрицы случайными числами
	for (int i = 0;i < len;i++)
	{
		Matrix_vec[i] = new __int64[len];
		Matrix_sc[i] = new __int64[len];
		Matrix_a[i] = new __int64[len];
		Matrix_b[i] = new __int64[len];

		for (int j = 0;j < len;j++)
		{
			Matrix_vec[i][j]=0;
			Matrix_sc[i][j]=0;
			Matrix_a[i][j] = 1 + rand() % 400;
			Matrix_b[i][j] = 1 + rand() % 400;
		}
	}

	//считаем врем€ дл€ скал€рного метода
	auto scl_start= chrono::high_resolution_clock::now();
	Scl(Matrix_a, Matrix_b, Matrix_sc);
	auto scl_end= chrono::high_resolution_clock::now();
	chrono::duration<double> scl_dur = (scl_end - scl_start)*1000;

	//транспониурем матрицу Matrix_b
	__int64** TemporaryMatrix = new __int64* [len];
	for (int i = 0;i < len;i++)
	{
		TemporaryMatrix[i] = new __int64[len];
		for (int j = 0;j < len;j++)
		{
			TemporaryMatrix[i][j] = Matrix_b[j][i];
		}
	}

	//мен€ем матрицу 
	for (int i = 0;i < len;i++)
	{
		for (int j = 0;j < len;j++)
		{
			Matrix_b[i][j]=TemporaryMatrix[i][j];
		}
	}

	//очищаем пам€ть
	for (int i=0; i<len;i++)
	{
		delete[] TemporaryMatrix[i];
	}
	delete[] TemporaryMatrix;

	//начинаем считать врем€ дл€ векторного рассчЄта
	auto vec_start=chrono::high_resolution_clock::now();
	Vec(Matrix_a, Matrix_b, Matrix_vec);
	auto vec_end=chrono::high_resolution_clock::now();
	chrono::duration<double> vec_dur = (vec_end - vec_start) * 1000;

	Cmp(Matrix_sc,Matrix_vec);

	cout << "Scalar time: " << scl_dur.count()<<endl;
	cout << "Vector time: " << vec_dur.count() << endl;

	for (int i = 0; i < len;i++)
	{
		delete[] Matrix_vec[i];
		delete[] Matrix_sc[i];
		delete[] Matrix_a[i];
		delete[] Matrix_b[i];
	}
	delete[] Matrix_vec;
	delete[] Matrix_sc;
	delete[] Matrix_a;
	delete[] Matrix_b;
}

//функци€ дл€ скал€рного(стандартного) умножени€
void Scl(__int64** matrix_a, __int64** matrix_b, __int64** matrix_c)
{
	for (int i=0;i< len;i++)
	{
		for (int j=0;j< len;j++)
		{
			for (int k=0;k< len;k++)
			{
				matrix_c[i][j] += matrix_a[i][k] * matrix_b[k][j];
			}
		}
	}
}

//функци€ векторного вычислени€ произведени€
void Vec(__int64** matrix_a, __int64** matrix_b, __int64** matrix_c)
{
	__m128i temp_1;
	__m128i sum;
	__m128i vector_a_1;
	__m128i vector_b_1;
	__m128i temp_2;
	__m128i vector_a_2;
	__m128i vector_b_2;
	for (int i = 0;i < len;i++)
	{
		temp_1 = _mm_setzero_si128();
		temp_2 = _mm_setzero_si128();
		sum = _mm_setzero_si128();
		for (int j = 0;j < len;j++)
		{
			for (int k = 0;k < len;k+=4)
			{
				vector_a_1 = _mm_load_si128((__m128i*) &matrix_a[i][k]);
				vector_b_1 = _mm_load_si128((__m128i*) & matrix_b[j][k]);
				vector_a_2 = _mm_load_si128((__m128i*) & matrix_a[i][k+2]);
				vector_b_2 = _mm_load_si128((__m128i*) & matrix_b[j][k+2]);
				temp_1 = _mm_mul_epu32(vector_a_1, vector_b_1);
				temp_2 = _mm_mul_epu32(vector_a_2, vector_b_2);
				sum = _mm_add_epi64(temp_1, sum);
				sum = _mm_add_epi64(temp_2, sum);
			}
			matrix_c[i][j] += _mm_extract_epi64(sum, 1) + _mm_extract_epi64(sum, 0);
			temp_1 = _mm_setzero_si128();
			sum = _mm_setzero_si128();
			temp_2 = _mm_setzero_si128();
		}
	}
}

//функци€ сравнени€ матриц дл€ проверки результатов
void Cmp(__int64** matrix_v, __int64** matrix_c)
{
	bool isEqual = true;
	for (int i = 0;i < len;i++)
	{
		for (int j = 0;j < len;j++)
		{
			if (matrix_v[i][j] != matrix_c[i][j])
			{
				isEqual=false;
				break;
			}
		}
		if (isEqual == false)
			break;
	}
	if (isEqual)
	{
		cout << "All right, matrixs are equal \n";
	}
	else
	{
		cout << "Oh my god, matrixs aren't equal \n";
	}
}
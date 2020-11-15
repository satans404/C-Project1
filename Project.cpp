#pragma GCC optimize(3)
#pragma GCC optimize("Ofast")
#pragma GCC optimize("inline")
#include <bits/stdc++.h>
#include <random>
#include <omp.h>
#include <immintrin.h>
#include <avx2intrin.h>
#include <cblas.h>
using namespace std;
const int Bolcksize = 128;
long long n, k;
struct Matrix
{
    float *matrix;
    int n, m;
    Matrix(int n, int m)
    {
        //matrix = new float[n * m + 10];
        matrix = (float *)aligned_alloc(32, (n * m + n * m % 32 + 32) * sizeof(float));
        this->n = n;
        this->m = m;
    }
    void Random()
    {
        default_random_engine e;
        uniform_real_distribution<double> u(-1.2, 3.5);
        int tmp;
        for (int i = 0; i < n; ++i)
        {
            tmp = i * m;
            for (int j = 0; j < m; ++j)
            {
                matrix[i * n + j] = 1;//u(e);
            }
        }
    }
    Matrix Brute_force1(const Matrix &b) //暴力1
    {
        Matrix ans(this->n, b.m);
        if (this->m != b.n)
        {
            puts("Error:These two matrices can not be product!");
            return ans;
        }
        int tmp = -n, tmpb, x;
        for (int i = 0; i < n; ++i)
        {
            //tmp = i * n;
            tmp += n;
            tmpb = -b.n;
            for (int k = 0; k < n; ++k)
            {
                //tmpb = k * b.n;
                tmpb += b.n;
                x = this->matrix[tmp + k];
                for (int j = 0; j < b.m; ++j)
                {
                    ans.matrix[tmp + j] += (x * b.matrix[tmpb + j]);
                    //ans.matrix[tmp + j] += (matrix[tmp + k] * b.matrix[tmpb + j]);
                }
            }
        }
        return ans;
    }
    Matrix Brute_force12(const Matrix &b) //暴力1
    {
        Matrix ans(this->n, b.m);
        if (this->m != b.n)
        {
            puts("Error:These two matrices can not be product!");
            return ans;
        }
        int tmp, tmpb, x;
        for (int i = 0; i < n; ++i)
        {
            tmp = i * n;
            for (int k = 0; k < n; ++k)
            {
                tmpb = k * b.n;
                x = this->matrix[tmp + k];
                for (int j = 0; j < b.m; ++j)
                {
                    ans.matrix[i * n + j] += (x * b.matrix[k * b.n + j]);
                }
            }
        }
        return ans;
    }
    Matrix Brute_force13(const Matrix &b) //暴力1
    {
        Matrix ans(this->n, b.m);
        if (this->m != b.n)
        {
            puts("Error:These two matrices can not be product!");
            return ans;
        }
        int tmp, tmpb, x;
#pragma omp parallel for num_threads(32)
        for (int i = 0; i < n; ++i)
        {
            tmp = i * n;
            for (int k = 0; k < n; ++k)
            {
                tmpb = k * b.n;
                x = this->matrix[tmp + k];
                for (int j = 0; j < b.m; ++j)
                {
                    ans.matrix[tmp + j] += (x * b.matrix[tmpb + j]);
                }
            }
        }
        return ans;
    }
    Matrix Brute_force2(const Matrix &b) //暴力1
    {
        Matrix ans(this->n, b.m);
        if (this->m != b.n)
        {
            puts("Error:These two matrices can not be product!");
            return ans;
        }
        int tmp, tmpb;
        for (int i = 0; i < n; ++i)
        {
            tmp = i * n;
            for (int j = 0; j < b.m; ++j)
            {
                for (int k = 0; k < n; ++k)
                {
                    ans.matrix[tmp + j] += (this->matrix[tmp + k] * b.matrix[tmpb + j]);
                }
            }
        }
        return ans;
    }
    Matrix Brute_force3(const Matrix &b) //暴力1
    {
        Matrix ans(this->n, b.m);
        if (this->m != b.n)
        {
            puts("Error:These two matrices can not be product!");
            return ans;
        }
        int tmp, tmpb;
        for (int k = 0; k < n; ++k)
        {
            tmp = k * n;
            for (int i = 0; i < n; ++i)
            {
                tmpb = i * b.n;
                for (int j = 0; j < b.m; ++j)
                {
                    ans.matrix[tmp + j] += (this->matrix[tmp + k] * b.matrix[tmpb + j]);
                }
            }
        }
        return ans;
    }

    inline void BlockSum(Matrix &ans, const Matrix &b, int I, int J, int K)
    {
        int tmp, tmpb, ii = min(n, I + Bolcksize), kk = min(n, K + Bolcksize), jj = min(b.m, J + Bolcksize);
        float x = 0;
        for (int i = I; i < ii; ++i)
        {
            tmp = i * n;
            for (int k = K; k < kk; ++k)
            {
                tmpb = k * b.n;
                x = matrix[tmp + k];
                for (int j = J; j < jj; ++j)
                {
                    ans.matrix[tmp + j] += (x * b.matrix[tmpb + j]);
                }
            }
        }
    }
    Matrix First(const Matrix &b)
    {
        Matrix ans(this->n, b.m);
        if (this->m != b.n)
        {
            puts("Error:These two matrices can not be product!");
            return ans;
        }
        int tmp, tmpb;
        for (int i = 0; i < n; i += Bolcksize)
        {
            tmp = i * n;
            for (int k = 0; k < n; k += Bolcksize)
            {
                tmpb = k * b.n;
                for (int j = 0; j < b.m; j += Bolcksize)
                {
                    BlockSum(ans, b, i, j, k);
                }
            }
        }
        return ans;
    }
    Matrix Second(Matrix &b)
    {
        Matrix ans(this->n, b.m);
        if (this->m != b.n)
        {
            puts("Error:These two matrices can not be product!");
            return ans;
        }
        register int tmp = -n, tmpb = -b.n, x = n, xx = b.n;
#pragma omp parallel for num_threads(32)
        for (int i = 0; i < n; i += Bolcksize)
        {
            tmp += x;
            tmpb = -b.n;
            for (int k = 0; k < n; k += Bolcksize)
            {
                tmpb += xx;
                for (int j = 0; j < b.m; j += Bolcksize)
                {
                    BlockSum(ans, b, i, j, k);
                }
            }
        }
        return ans;
    }
    inline void _BlockSum(Matrix &ans, Matrix &b, int I, int J, int K)
    {

        register int tmp = (I - 1) * n, tmpb = K * b.n, ii = min(n, I + Bolcksize), kk = min(n, K + Bolcksize), jj = min(b.m, J + Bolcksize), xx = 8 * b.n, tmpc = (K - 8) * b.n, tmpd = (K - 8) * b.n;
        __m256 A, B, C, D;
        for (register int i = I; i < ii; ++i)
        {
            tmp += n;
            tmpb = (K - 8) * b.n;
            for (register int k = K; k < kk; k += 8)
            {
                tmpb += xx;
                B = _mm256_load_ps(matrix + tmp + k);
                //B=_mm256_set1_ps(1.0);
                for (register int j = J; j < jj; j += 8)
                {
                    /*A = _mm256_loadu_ps(b.matrix + tmpb + j);
                    C = _mm256_loadu_ps(ans.matrix + tmp + j);
                    C = _mm256_add_ps(_mm256_mul_ps(A, B), C);
                    _mm256_storeu_ps(ans.matrix + tmp + j, C);
                    //A = _mm256_load_ps(b.matrix + tmpb + j);
                    //C = _mm256_load_ps(ans.matrix + tmp + j);
                    //C = _mm256_add_ps(_mm256_mul_ps(A, B), C);*/
                    C = _mm256_fmadd_ps(_mm256_load_ps(b.matrix + tmpb + j), B, _mm256_load_ps(ans.matrix + tmp + j));
                    _mm256_store_ps(ans.matrix + tmp + j, C);
                }
            }
        }
    }
    inline void _BlockSum2(Matrix &ans, Matrix &b, int I, int J, int K)
    {
        register int tmp = (I - 1) * n, tmpb = K * b.n, ii = min(n, I + Bolcksize), kk = min(n, K + Bolcksize), jj = min(b.m, J + Bolcksize);
        __m256 A, B, C, D;
        for (register int i = I; i < ii; ++i)
        {
            tmp = i * n;

            for (register int k = K; k < kk; k += 8)
            {
                tmpb = k * b.n;
                B = _mm256_load_ps(matrix + tmp + k);
                for (register int j = J; j < jj; j += 8)
                {
                    C = _mm256_fmadd_ps(_mm256_load_ps(b.matrix + tmpb + j), B, _mm256_load_ps(ans.matrix + tmp + j));
                    _mm256_store_ps(ans.matrix + tmp + j, C);
                }
            }
        }
    }
    inline void _BlockSum3(Matrix &ans, Matrix &b, int I, int J, int K)
    {
        register int tmp = (I - 1) * n, tmpb = K * b.n, ii = min(n, I + Bolcksize), kk = min(n, K + Bolcksize), jj = min(b.m, J + Bolcksize);
        __m256 A, B, C, D;
        for (register int i = I; i < ii; ++i)
        {
            tmp = i * n;

            for (register int k = K; k < kk; k += 8)
            {
                tmpb = k * b.n;
                B = _mm256_loadu_ps(matrix + tmp + k);
                for (register int j = J; j < jj; j += 8)
                {
                    A = _mm256_loadu_ps(b.matrix + tmpb + j);
                    C = _mm256_loadu_ps(ans.matrix + tmp + j);
                    C = _mm256_add_ps(_mm256_mul_ps(A, B), C);
                    _mm256_storeu_ps(ans.matrix + tmp + j, C);
                    //A = _mm256_load_ps(b.matrix + tmpb + j);
                    //C = _mm256_load_ps(ans.matrix + tmp + j);
                    //C = _mm256_add_ps(_mm256_mul_ps(A, B), C);
                    //C = _mm256_fmadd_ps(_mm256_load_ps(b.matrix + tmpb + j), B, _mm256_load_ps(ans.matrix + tmp + j));
                    //_mm256_store_ps(ans.matrix + tmp + j, C);
                }
            }
        }
    }
    Matrix Third(Matrix &b)
    {
        Matrix ans(this->n, b.m);
        if (this->m != b.n)
        {
            puts("Error:These two matrices can not be product!");
            return ans;
        }
        if (m * b.n % 32)
        {
            return Third3(b);
        }
        register int tmp = 0, tmpb = 0, x = Bolcksize * n, xx = Bolcksize * b.n;
#pragma omp parallel for num_threads(32)
        for (int i = 0; i < n; i += Bolcksize)
        {
            //tmp = i * n;
            for (int k = 0; k < n; k += Bolcksize)
            {
                //tmpb = k * b.n;
                for (int j = 0; j < b.m; j += Bolcksize)
                {
                    _BlockSum(ans, b, i, j, k);
                    //ans.out();
                }
                tmpb += xx;
            }
            tmp += x;
        }
        return ans;
    }
    Matrix Third2(Matrix &b)
    {
        Matrix ans(this->n, b.m);
        if (this->m != b.n)
        {
            puts("Error:These two matrices can not be product!");
            return ans;
        }
        if (m * b.n % 32)
        {
            return Third3(b);
        }
        register int tmp = 0, tmpb = 0, x = Bolcksize * n, xx = Bolcksize * b.n;
#pragma omp parallel for num_threads(32)
        for (int i = 0; i < n; i += Bolcksize)
        {
            //tmp = i * n;
            for (int k = 0; k < n; k += Bolcksize)
            {
                //tmpb = k * b.n;
                for (int j = 0; j < b.m; j += Bolcksize)
                {
                    _BlockSum2(ans, b, i, j, k);
                    //ans.out();
                }
                tmpb += xx;
            }
            tmp += x;
        }
        return ans;
    }
    Matrix Third3(Matrix &b)
    {
        Matrix ans(this->n, b.m);
        if (this->m != b.n)
        {
            puts("Error:These two matrices can not be product!");
            return ans;
        }
        register int tmp = 0, tmpb = 0, x = Bolcksize * n, xx = Bolcksize * b.n;
#pragma omp parallel for num_threads(32)
        for (int i = 0; i < n; i += Bolcksize)
        {
            //tmp = i * n;
            for (int k = 0; k < n; k += Bolcksize)
            {
                //tmpb = k * b.n;
                for (int j = 0; j < b.m; j += Bolcksize)
                {
                    _BlockSum3(ans, b, i, j, k);
                }
                tmpb += xx;
            }
            tmp += x;
        }
        return ans;
    }
    void out()
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                printf("%.2lf ", matrix[i * n + j]);
            }
            puts("");
        }
    }
};
void Brute_force_test(int na, int ma, int nb, int mb) //暴力1
{

    float *c[na + 10], *b[nb + 10], *a[na + 10];
    for (int i = 0; i <= na; i++)
    {
        c[i] = new float[mb + 10];
        a[i] = new float[ma + 10];
    }
    for (int i = 0; i <= nb; i++)
    {
        b[i] = new float[mb + 10];
    }
    clock_t startTime, endTime;
    startTime = clock();
    for (int i = 0; i < na; ++i)
    {
        for (int k = 0; k < na; ++k)
        {
            for (int j = 0; j < mb; ++j)
            {
                c[i][j] += (a[i][k] * b[k][j]);
            }
        }
    }
    endTime = clock(); 
    cout << "The run time is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
}
void openblas_test(int M,int N,int K)
{
    float *A = new float[M * K];
	float *B = new float[K * N];
	float *C = new float[M * N];
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1, A, K, B, N, 0, C, N);
}
signed main()
{

    clock_t startTime, endTime;
    startTime = clock();
    Matrix A(15000, 15000);
    A.Random();
    endTime = clock(); 
    cout << "The run time is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    startTime = clock();
    A.Second(A);
    endTime = clock();
    cout << "The run time is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
}

#include <CppUTest/TestHarness.h>
#include <CppUTestExt/MockSupport.h>

extern "C"
{
#include "matrix.h"
#include <stdio.h>
}

#define EPSILON (0.000001)
#define M_  5
#define N_  5

TEST_GROUP( MATRIX )
{
    matrix * dut;

    void setup()
    {
        dut = matrix_alloc(M_, N_);
    }

    void teardown()
    {
        matrix_free(dut);
    }
};

TEST(MATRIX, alloc_and_free)
{
    matrix * a = matrix_alloc(5, 5);
    matrix * b = matrix_calloc(5, 5);

    matrix_free(a);
    matrix_free(b);
}

TEST(MATRIX, set_zero)
{
    matrix * a = matrix_alloc(5, 5);

    matrix_set_zero(a);

    for (int i=0; i<25; i++)
    {
        DOUBLES_EQUAL(0, a->data[i], EPSILON);
    }

    matrix_free(a);
}

TEST(MATRIX, set_all)
{
    matrix * a = matrix_alloc(5, 5);

    matrix_set_all(a, 42);

    for (size_t i=0; i<25; i++)
    {
        DOUBLES_EQUAL(42, a->data[i], EPSILON);
    }

    matrix_free(a);
}

TEST(MATRIX, get_set)
{
    const size_t m = 5, n = 5;
    matrix * a = matrix_alloc(m, n);
    double set, result;

    for (size_t i=0; i<m; i++)
    {
        for (size_t j=0; j<n; j++)
        {
            set = i*n + j;
            matrix_set(a, i, j, set);
            result = matrix_get(a, i, j);
            DOUBLES_EQUAL(set, result, EPSILON);
        }
    }

    matrix_free(a);
}

TEST(MATRIX, set_identity)
{
    const size_t m = 5, n = 5;
    double expected, result;

    matrix * a = matrix_alloc(m, n);

    matrix_set_identity(a);

    for (size_t i=0; i<m; i++)
    {
        for (size_t j=0; j<n; j++)
        {
            expected = (i==j) ? 1 : 0;
            result = matrix_get(a, i, j);
            DOUBLES_EQUAL(expected, result, EPSILON);
        }
    }

    matrix_free(a);
}

TEST(MATRIX, mem_cpy)
{
    const size_t m = 5, n = 5;
    double expected, result;
    int error;

    matrix * a = matrix_alloc(m, n);
    matrix * b = matrix_calloc(m, n);

    error = matrix_memcpy(b, a);

    for (size_t i=0; i<m; i++)
    {
        for (size_t j=0; j<n; j++)
        {
            expected = matrix_get(a, i, j);
            result = matrix_get(b, i, j);
            DOUBLES_EQUAL(expected, result, EPSILON);
        }
    }

    CHECK_FALSE(error);

    matrix_free(a);
    matrix_free(b);
}

TEST(MATRIX, swap)
{
    const size_t m = 5, n = 5;
    double expected, result, set;
    int error;

    matrix * a = matrix_alloc(m, n);
    matrix * b = matrix_alloc(m, n);

    for (size_t i=0; i<m; i++)
    {
        for (size_t j=0; j<n; j++)
        {
            set = i*n + j;
            matrix_set(a, i, j, set);
            matrix_set(b, i, j, 24-set);
        }
    }

    error = matrix_swap(a, b);

    for (size_t i=0; i<m; i++)
    {
        for (size_t j=0; j<n; j++)
        {
            expected = 24 - (i*n + j);
            result = matrix_get(a, i, j);
            DOUBLES_EQUAL(expected, result, EPSILON);

            expected = i*n + j;
            result = matrix_get(b, i, j);
            DOUBLES_EQUAL(expected, result, EPSILON);
        }
    }

    CHECK_FALSE(error);

    matrix_free(a);
    matrix_free(b);
}

TEST(MATRIX, swap2)
{
    const size_t m = 5, n = 5;
    double result;
    int error, error2, error3;
    size_t i = 0, j = 0;

    matrix * a = matrix_alloc(m, n);
    matrix * b = matrix_alloc(m, n);
    matrix * c = matrix_calloc(m, n);

    for (i=0; i<m; i++)
    {
        for (j=0; j<n; j++)
        {
            matrix_set(a, i, j, i);
            matrix_set(b, i, j, j);
            matrix_set(c, 0, j, 1);
        }
    }

    error = matrix_swap_rows(a, 0, 1);
    error2 = matrix_swap_columns(b, 0, 1);
    error3 = matrix_swap_rowcol(c, 0, 0);

    for (j=0; j<n; j++)
    {
        result = matrix_get(a, 0, j);
        DOUBLES_EQUAL(1, result, EPSILON);
        result = matrix_get(a, 1, j);
        DOUBLES_EQUAL(0, result, EPSILON);
    }

    for (i=0; i<m; i++)
    {
        result = matrix_get(b, i, 0);
        DOUBLES_EQUAL(1, result, EPSILON);
        result = matrix_get(b, i, 1);
        DOUBLES_EQUAL(0, result, EPSILON);

        result = matrix_get(c, i, 0);
        DOUBLES_EQUAL(1, result, EPSILON);
    }

    CHECK_FALSE(error);
    CHECK_FALSE(error2);
    CHECK_FALSE(error3);

    matrix_free(a);
    matrix_free(b);
    matrix_free(c);
}

TEST(MATRIX, equal)
{
    const size_t m = 5, n = 5;
    int isequal;

    matrix * a = matrix_calloc(m, n);
    matrix * b = matrix_calloc(m, n);

    for (size_t i=0; i<m; i++)
    {
        for (size_t j=0; j<n; j++)
        {
            matrix_set(a, i, j, i*n + j);
            matrix_set(b, j, i, i*n + j);
        }
    }

    isequal = matrix_equal(a, b);
    CHECK_FALSE(isequal);

    matrix_memcpy(a, b);
    isequal = matrix_equal(a, b);
    CHECK_TRUE(isequal);

    matrix_free(a);
    matrix_free(b);
}

TEST(MATRIX, transpose)
{
    const size_t n = 3;
    int error, isequal;
    size_t i = 0, j = 0;

    matrix * a = matrix_alloc(n, n);
    matrix * b = matrix_alloc(n, n);
    matrix * c = matrix_alloc(n, n);

    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
        {
            matrix_set(a, i, j, i*n + j);
            matrix_set(b, j, i, i*n + j);
        }
    }

    error = matrix_transpose_memcpy(c, a);
    CHECK_FALSE(error);
    isequal = matrix_equal(b, c);
    CHECK_TRUE(isequal);

    error = matrix_transpose(a);
    CHECK_FALSE(error);
    isequal = matrix_equal(b, c);
    CHECK_TRUE(isequal);

    matrix_free(a);
    matrix_free(b);
    matrix_free(c);
}

TEST(MATRIX, minmax)
{
    size_t i = 0, j = 0;
    double max = 0, min = 100;
    size_t imax = 100, jmax = 100, imin = 100, jmin = 100;

    for (i=0; i<M_; i++)
    {
        for (j=0; j<N_; j++)
        {
            matrix_set(dut, i, j, i*M_ + j);
        }
    }

    max = matrix_max(dut);
    matrix_max_index(dut, &imax, &jmax);
    min = matrix_min(dut);
    matrix_min_index(dut, &imin, &jmin);

    DOUBLES_EQUAL((M_*N_ -1), max, EPSILON);
    DOUBLES_EQUAL(0, min, EPSILON);
    CHECK_EQUAL(0, imin);
    CHECK_EQUAL(0, jmin);
    CHECK_EQUAL(M_ - 1, imax);
    CHECK_EQUAL(N_ - 1, jmax);

    max = 0;
    min = 100;
    imax = 100;
    jmax = 100;
    imin = 100;
    jmin = 100;

    matrix_minmax(dut, &min, &max);
    matrix_minmax_index(dut, &imin, &jmin, &imax, &jmax);

    DOUBLES_EQUAL((M_*N_ -1), max, EPSILON);
    DOUBLES_EQUAL(0, min, EPSILON);
    CHECK_EQUAL(0, imin);
    CHECK_EQUAL(0, jmin);
    CHECK_EQUAL(M_ - 1, imax);
    CHECK_EQUAL(N_ - 1, jmax);
}

TEST(MATRIX, add)
{
    size_t i, j;
    matrix * b = matrix_alloc(M_, N_);
    int error = 0;

    for (i=0; i<M_; i++)
    {
        for (j=0; j<N_; j++)
        {
            matrix_set(dut, i, j, 42);
            matrix_set(b, i, j, 43);
        }
    }

    error = matrix_add(dut, b);
    CHECK_FALSE(error);

    for (i=0; i<M_; i++)
    {
        for (j=0; j<N_; j++)
        {
            DOUBLES_EQUAL(85, matrix_get(dut, i, j), EPSILON);
        }
    }

    matrix_free(b);
}

TEST(MATRIX, sub)
{
    size_t i, j;
    matrix * b = matrix_alloc(M_, N_);
    int error = 0;

    for (i=0; i<M_; i++)
    {
        for (j=0; j<N_; j++)
        {
            matrix_set(dut, i, j, 43);
            matrix_set(b, i, j, 42);
        }
    }

    error = matrix_sub(dut, b);
    CHECK_FALSE(error);

    for (i=0; i<M_; i++)
    {
        for (j=0; j<N_; j++)
        {
            DOUBLES_EQUAL(1, matrix_get(dut, i, j), EPSILON);
        }
    }

    matrix_free(b);
}

TEST(MATRIX, mul_elem)
{
    size_t i, j;
    matrix * b = matrix_alloc(M_, N_);
    int error = 0;

    for (i=0; i<M_; i++)
    {
        for (j=0; j<N_; j++)
        {
            matrix_set(dut, i, j, 43);
            matrix_set(b, i, j, 42);
        }
    }

    error = matrix_mul_elements(dut, b);
    CHECK_FALSE(error);

    for (i=0; i<M_; i++)
    {
        for (j=0; j<N_; j++)
        {
            DOUBLES_EQUAL(1806, matrix_get(dut, i, j), EPSILON);
        }
    }

    matrix_free(b);
}

TEST(MATRIX, div_elem)
{
    size_t i, j;
    matrix * b = matrix_alloc(M_, N_);
    int error = 0;

    for (i=0; i<M_; i++)
    {
        for (j=0; j<N_; j++)
        {
            matrix_set(dut, i, j, 120);
            matrix_set(b, i, j, 42);
        }
    }

    error = matrix_div_elements(dut, b);
    CHECK_FALSE(error);

    for (i=0; i<M_; i++)
    {
        for (j=0; j<N_; j++)
        {
            DOUBLES_EQUAL(2.857142857142857, matrix_get(dut, i, j), EPSILON);
        }
    }

    matrix_free(b);
}

TEST(MATRIX, scale)
{
    size_t i, j;
    int error = 0;

    for (i=0; i<M_; i++)
    {
        for (j=0; j<N_; j++)
        {
            matrix_set(dut, i, j, 120);
        }
    }

    error = matrix_scale(dut, 0.1);
    CHECK_FALSE(error);

    for (i=0; i<M_; i++)
    {
        for (j=0; j<N_; j++)
        {
            DOUBLES_EQUAL(12, matrix_get(dut, i, j), EPSILON);
        }
    }
}

TEST(MATRIX, add_constant)
{
    size_t i, j;
    int error = 0;

    for (i=0; i<M_; i++)
    {
        for (j=0; j<N_; j++)
        {
            matrix_set(dut, i, j, 43);
        }
    }

    error = matrix_add_constant(dut, -1);
    CHECK_FALSE(error);

    for (i=0; i<M_; i++)
    {
        for (j=0; j<N_; j++)
        {
            DOUBLES_EQUAL(42, matrix_get(dut, i, j), EPSILON);
        }
    }
}

TEST(MATRIX, add_diag)
{
    size_t i, j;
    int error = 0;

    for (i=0; i<M_; i++)
    {
        for (j=0; j<N_; j++)
        {
            matrix_set(dut, i, j, 43);
        }
    }

    error = matrix_add_diagonal(dut, -1);
    CHECK_FALSE(error);

    for (i=0; i<M_; i++)
    {
        for (j=0; j<N_; j++)
        {
            if (i == j)
            {
                DOUBLES_EQUAL(42, matrix_get(dut, i, j), EPSILON);
            }
            else
            {
                DOUBLES_EQUAL(43, matrix_get(dut, i, j), EPSILON);
            }
        }
    }
}

TEST(MATRIX, mul)
{
    int error = 0;
    int isequal = 0;
    const size_t m1 = 3, n1 = 2, m2 = 2, n2 = 4;
    matrix * a = matrix_alloc(m1, n1);
    matrix * b = matrix_alloc(m2, n2);
    matrix * c = matrix_calloc(m1, n2);
    double arr[12] = {  11, 14, 17, 20,\
                        23, 30, 37, 44,\
                        35, 46, 57, 68};
    matrix * res = matrix_alloc_array(m1, n2, arr, sizeof(arr)/sizeof(arr[0]));

    for (int i=0; i<m1; i++)
    {
        for(int j=0; j<n1; j++)
        {
            matrix_set(a, i, j, i*n1 + j + 1);
        }
    }

    for (int i=0; i<m2; i++)
    {
        for(int j=0; j<n2; j++)
        {
            matrix_set(b, i, j, i*n2 + j + 1);
        }
    }

    error = matrix_mul(a, b, c);
    CHECK_FALSE(error);

    isequal = matrix_equal(c, res);
    CHECK(isequal);
}
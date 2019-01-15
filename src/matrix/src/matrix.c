/***********************************************************************************
*	File Name: matrix.c
*	Version 1.00
*
*	Description:
*
*	Note:
*
************************************************************************************
*	MIT License
*
*	Copyright (c) 2019 Frank Sodari
*
*	Permission is hereby granted, free of charge, to any person obtaining a copy
*	of this software and associated documentation files (the "Software"), to deal
*	in the Software without restriction, including without limitation the rights
*	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*	copies of the Software, and to permit persons to whom the Software is
*	furnished to do so, subject to the following conditions:
*
*	The above copyright notice and this permission notice shall be included in all
*	copies or substantial portions of the Software.
*
*	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
*	SOFTWARE.
***********************************************************************************/

#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define ACCESS(m, i, j) m->data[i*m->n_cols + j]

matrix * matrix_alloc (const size_t n_rows, const size_t n_cols)
{
    matrix * m = NULL;

    // 0 rows or columns is undefined
    if (n_rows && n_cols)
    {
        m = malloc(sizeof(*m));

        if (m)
        {
            m->n_rows = n_rows;
            m->n_cols = n_cols;
            m->data = malloc(n_rows*n_cols*sizeof(double));

            if (m->data == NULL)
            {
                free(m);
                m = NULL;
            }
        }
    }

    return (m);
}

matrix * matrix_calloc (const size_t n_rows, const size_t n_cols)
{
    matrix * m = NULL;

    // 0 rows or columns is undefined
    if (n_rows && n_cols)
    {
        m = malloc(sizeof(*m));

        if (m)
        {
            m->n_rows = n_rows;
            m->n_cols = n_cols;
            m->data = calloc(n_rows*n_cols*sizeof(double), sizeof(double));

            if (m->data == NULL)
            {
                free(m);
                m = NULL;
            }
        }
    }

    return (m);
}

matrix * matrix_alloc_array(const size_t n_rows, const size_t n_cols, const double arr[], size_t array_length)
{
    matrix * m = NULL;
    size_t i, j;

    // 0 rows or columns is undefined. Array length must equal rows*cols
    if (n_rows && n_cols && (n_rows*n_cols == array_length))
    {
        m = malloc(sizeof(*m));

        if (m)
        {
            m->n_rows = n_rows;
            m->n_cols = n_cols;
            m->data = malloc(n_rows*n_cols*sizeof(double));

            if (m->data == NULL)
            {
                free(m);
                m = NULL;
            }

            else
            {
                for (i=0; i < m->n_rows; i++)
                {
                    for (j=0; j < m->n_cols; j++)
                    {
                        ACCESS(m, i, j) = arr[i*m->n_cols + j];
                    }
                }
            }
        }
    }

    return m;
}

void matrix_free (matrix * m)
{
    free(m->data);
    free(m);
}

void matrix_set_zero (matrix * m)
{
    size_t i;
    size_t n = m->n_rows * m->n_cols;

    for (i=0; i<n; i++)
    {
        m->data[i] = 0;
    }
}

void matrix_set_identity (matrix * m)
{
    size_t i, j;

    for (i=0; i<m->n_rows; i++)
    {
        for (j=0; j<m->n_cols; j++)
        {
            ACCESS(m, i, j) = (i==j) ? 1 : 0;
        }
    }
}

void matrix_set_all (matrix * m, const double x)
{
    size_t i;
    size_t n = m->n_rows * m->n_cols;

    for (i=0; i<n; i++)
    {
        m->data[i] = x;
    }
}

int matrix_set_array(matrix *m, const double arr[], size_t length)
{
    int error = 0;
    size_t i, j;

    if (length != (m->n_rows * m->n_cols))
    {
        error = -1;
    }

    else
    {
        for (i = 0; i < m->n_rows; i++)
        {
            for (j = 0; j < m->n_cols; j++)
            {
                ACCESS(m, i, j) = arr[i*m->n_cols + j];
            }
        }
    }

    return error;
}

int matrix_memcpy (matrix * dest, const matrix * src)
{
    int error = 0;

    if (dest->n_rows != src->n_rows || dest->n_cols != src->n_cols)
    {
        error = -1;
    }

    else
    {
        memcpy(dest->data, src->data, src->n_rows * src->n_cols * sizeof(double));
    }

    return (error);
}

int matrix_swap (matrix * a, matrix * b)
{
    int error = 0;
    size_t  i, j;
    double temp;

    if (a->n_rows != b->n_rows || a->n_cols != b->n_cols)
    {
        error = -1;
    }

    else
    {
        for (i=0; i<a->n_rows; i++)
        {
            for(j=0; j<a->n_cols; j++)
            {
                temp = ACCESS(a, i, j);
                ACCESS(a, i, j) = ACCESS(b, i, j);
                ACCESS(b, i, j) = temp;
            }
        }
    }

    return (error);
}

int matrix_swap_rows (matrix * m, const size_t i, const size_t  j)
{
    size_t jj;
    double temp;
    int error = 0;

    if (i >= m->n_rows || j >= m->n_rows)
    {
        error = -1;
    }

    else
    {
        for (jj = 0; jj < m->n_cols; jj++)
        {
            temp = ACCESS(m, i, jj);
            ACCESS(m, i, jj) = ACCESS(m, j, jj);
            ACCESS(m, j, jj) = temp;
        }
    }

    return (error);
}

int matrix_swap_columns (matrix * m, const size_t i, const size_t j)
{
    size_t ii;
    double temp;
    int error = 0;

    if (i >= m->n_rows || j >= m->n_rows)
    {
        error = -1;
    }

    else
    {
        for (ii = 0; ii < m->n_cols; ii++)
        {
            temp = ACCESS(m, ii, i);
            ACCESS(m, ii, i) = ACCESS(m, ii, j);
            ACCESS(m, ii, j) = temp;
        }
    }

    return (error);
}

int matrix_swap_rowcol (matrix *m, const size_t i, const size_t j)
{
    size_t it;
    int error = 0;
    double temp;

    if (i >= m->n_rows || j >= m->n_rows)
    {
        error = -1;
    }

    else if (m->n_rows != m->n_cols)
    {
        error = -1;
    }

    else
    {
        for (it=0; it<m->n_rows; it++)
        {
            temp = ACCESS(m, i, it);
            ACCESS(m, i, it) = ACCESS(m, it, j);
            ACCESS(m, it, j) = temp;
        }
    }

    return (error);
}

int matrix_transpose (matrix * m)
{
    size_t i, j;
    int error = 0;
    double swap;

    // Transpose only allowed for square matrices
    if (m->n_rows != m->n_cols)
    {
        error = -1;
    }

    else
    {
        for (i=0; i<m->n_rows; i++)
        {
            for (j=i+1; j<m->n_cols; j++)
            {
                swap = ACCESS(m, i, j);
                ACCESS(m, i, j) = ACCESS(m, j, i);
                ACCESS(m, j, i) = swap;
            }
        }
    }

    return error;
}

int matrix_transpose_memcpy (matrix * dest, const matrix * src)
{
    size_t i, j;
    int error = 0;

    //
    if (src->n_rows != dest->n_cols || src->n_cols != dest->n_rows)
    {
        error = -1;
    }

    else
    {
        for (i=0; i<dest->n_rows; i++)
        {
            for (j=0; j<dest->n_cols; j++)
            {
                ACCESS(dest, i, j) = ACCESS(src, j, i);
            }
        }
    }

    return error;
}

double matrix_max (const matrix * m)
{
    size_t i, j;
    double max = ACCESS(m, 0, 0);
    double temp;

    for (i=0; i<m->n_rows; i++)
    {
        for (j=0; j<m->n_cols; j++)
        {
            temp = ACCESS(m, i, j);
            if (temp > max)
            {
                max = temp;
            }
        }
    }

    return max;
}

double matrix_min (const matrix * m)
{
    size_t i, j;
    double min = ACCESS(m, 0, 0);
    double temp;

    for (i=0; i<m->n_rows; i++)
    {
        for (j=0; j<m->n_cols; j++)
        {
            temp = ACCESS(m, i, j);
            if (temp < min)
            {
                min = temp;
            }
        }
    }

    return min;
}

void matrix_minmax (const matrix * m, double * min_out, double * max_out)
{
    size_t i, j;
    double temp;

    *min_out = ACCESS(m, 0, 0);
    *max_out = *min_out;

    for (i=0; i<m->n_rows; i++)
    {
        for (j=0; j<m->n_cols; j++)
        {
            temp = ACCESS(m, i, j);
            if (temp < *min_out)
            {
                *min_out = temp;
            }

            if (temp > *max_out)
            {
                *max_out = temp;
            }
        }
    }
}

void matrix_max_index (const matrix * m, size_t * imax, size_t * jmax)
{
    size_t i, j;
    double max = ACCESS(m, 0, 0);
    double temp;
    *imax = 0;
    *jmax = 0;

    for (i=0; i<m->n_rows; i++)
    {
        for (j=0; j<m->n_cols; j++)
        {
            temp = ACCESS(m, i, j);
            if (temp > max)
            {
                max = temp;
                *imax = i;
                *jmax = j;
            }
        }
    }
}

void matrix_min_index (const matrix * m, size_t * imin, size_t * jmin)
{
    size_t i, j;
    double min = ACCESS(m, 0, 0);
    double temp;

    *imin = 0;
    *jmin = 0;

    for (i=0; i<m->n_rows; i++)
    {
        for (j=0; j<m->n_cols; j++)
        {
            temp = ACCESS(m, i, j);
            if (temp < min)
            {
                min = temp;
                *imin = i;
                *jmin = j;
            }
        }
    }
}

void matrix_minmax_index (const matrix * m, size_t * imin, size_t * jmin, size_t * imax, size_t * jmax)
{
    size_t i, j;
    double min = ACCESS(m, 0, 0);
    double max = min;
    double temp;

    *imin = 0;
    *jmin = 0;
    *imax = 0;
    *jmax = 0;

    for (i=0; i<m->n_rows; i++)
    {
        for (j=0; j<m->n_cols; j++)
        {
            temp = ACCESS(m, i, j);
            if (temp < min)
            {
                min = temp;
                *imin = i;
                *jmin = j;
            }

            if (temp > max)
            {
                max = temp;
                *imax = i;
                *jmax = j;
            }
        }
    }
}

int matrix_equal (const matrix * a, const matrix * b)
{
    size_t  i, j;

    if (a->n_rows != b->n_rows || a->n_cols != b->n_cols)
    {
        return 0;
    }

    else
    {
        for (i=0; i<a->n_rows; i++)
        {
            for (j=0; j<a->n_cols; j++)
            {
                if (ACCESS(a, i, j) != ACCESS(b, i, j))
                {
                    return 0;
                }
            }
        }
    }

    return (1);
}

int matrix_add (matrix * a, const matrix * b)
{
    size_t i, j;
    int error = 0;

    if (a->n_rows != b->n_rows || a->n_cols != b->n_cols)
    {
        error = -1;
    }

    else
    {
        for (i=0; i<a->n_rows; i++)
        {
            for (j=0; j<a->n_cols; j++)
            {
                ACCESS(a, i, j) += ACCESS(b, i, j);
            }
        }
    }

    return error;
}

int matrix_sub (matrix * a, const matrix * b)
{
    size_t i, j;
    int error = 0;

    if (a->n_rows != b->n_rows || a->n_cols != b->n_cols)
    {
        error = -1;
    }

    else
    {
        for (i=0; i<a->n_rows; i++)
        {
            for (j=0; j<a->n_cols; j++)
            {
                ACCESS(a, i, j) -= ACCESS(b, i, j);
            }
        }
    }

    return error;
}

int matrix_mul_elements (matrix * a, const matrix * b)
{
    size_t i, j;
    int error = 0;

    if (a->n_rows != b->n_rows || a->n_cols != b->n_cols)
    {
        error = -1;
    }

    else
    {
        for (i=0; i<a->n_rows; i++)
        {
            for (j=0; j<a->n_cols; j++)
            {
                ACCESS(a, i, j) *= ACCESS(b, i, j);
            }
        }
    }

    return error;
}

int matrix_div_elements (matrix * a, const matrix * b)
{
    size_t i, j;
    int error = 0;

    if (a->n_rows != b->n_rows || a->n_cols != b->n_cols)
    {
        error = -1;
    }

    else
    {
        for (i=0; i<a->n_rows; i++)
        {
            for (j=0; j<a->n_cols; j++)
            {
                ACCESS(a, i, j) /= ACCESS(b, i, j);
            }
        }
    }

    return error;
}

int matrix_scale (matrix * a, const double x)
{
    size_t i, j;

    for (i=0; i<a->n_rows; i++)
    {
        for (j=0; j<a->n_cols; j++)
        {
            ACCESS(a, i, j) *= x;
        }
    }

    return 0;
}

int matrix_add_constant (matrix * a, const double x)
{
    size_t i, j;

    for (i=0; i<a->n_rows; i++)
    {
        for (j=0; j<a->n_cols; j++)
        {
            ACCESS(a, i, j) += x;
        }
    }

    return 0;
}

int matrix_add_diagonal (matrix * a, const double x)
{
    size_t i;
    int error = 0;

    if (a->n_rows != a->n_cols)
    {
        error = -1;
    }

    for (i=0; i<a->n_rows; i++)
    {

        ACCESS(a, i, i) += x;
    }

    return error;
}

int matrix_mul (const matrix * a, const matrix * b, matrix * c)
{
    size_t i, j, k;
    int error = 0;
    double sum;

    if (a->n_cols != b->n_rows || a->n_rows != c->n_rows || b->n_cols != c->n_cols)
    {
        error = -1;
    }
    else
    {
        // Iterate through result matrix
        for (i=0; i<c->n_rows; i++)
        {
            for (j=0; j<c->n_cols; j++)
            {
                sum = 0;
                for (k=0; k<a->n_cols; k++)
                {
                    sum += ACCESS(a, i, k) * ACCESS(b, k, j);
                }
                ACCESS(c, i, j) = sum;
            }
        }
    }

    return error;
}

double matrix_get (const matrix * m, const size_t i, const size_t j)
{
    if (i >= m->n_rows || j >= m->n_cols)
    {
        return -1;
    }

    return ACCESS(m, i, j);
}

void matrix_set (matrix * m, const size_t i, const size_t j, const double x)
{
    if (i >= m->n_rows || j >= m->n_cols)
    {
        return;
    }

    ACCESS(m, i, j) = x;
}

void matrix_print(matrix * m)
{
    size_t i, j;
    double element;

    printf("Matrix Contents: \n");
    for (i=0; i<m->n_rows; i++)
    {
        for (j=0; j<m->n_cols; j++)
        {
            element = ACCESS(m, i, j);
            printf("%lf, ", element);
        }
        printf("\n");
    }
}

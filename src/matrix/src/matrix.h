/***********************************************************************************
*	File Name: matrix.h
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

#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h>

typedef struct
{
    size_t n_rows;
    size_t n_cols;
    double * data;
} matrix;

matrix * matrix_alloc (const size_t n_rows, const size_t n_cols);
matrix * matrix_calloc (const size_t n_rows, const size_t n_cols);
matrix * matrix_alloc_array(const size_t n_rows, const size_t n_cols, const double arr[], size_t array_length);
void matrix_free (matrix * m);

void matrix_set_zero (matrix * m);
void matrix_set_identity (matrix * m);
void matrix_set_all (matrix * m, const double x);
int matrix_set_array(matrix *m, const double arr[], size_t length);

int matrix_memcpy (matrix * dest, const matrix * src);
int matrix_swap (matrix * a, matrix * b);

int matrix_swap_rows (matrix * m, const size_t i, const size_t  j);
int matrix_swap_columns (matrix * m, const size_t i, const size_t j);
int matrix_swap_rowcol (matrix *m, const size_t i, const size_t j);
int matrix_transpose (matrix * m);
int matrix_transpose_memcpy (matrix * dest, const matrix * src);

double matrix_max (const matrix * m);
double matrix_min (const matrix * m);
void matrix_minmax (const matrix * m, double * min_out, double * max_out);

void matrix_max_index (const matrix * m, size_t * imax, size_t * jmax);
void matrix_min_index (const matrix * m, size_t * imin, size_t * jmin);
void matrix_minmax_index (const matrix * m, size_t * imin, size_t * jmin, size_t * imax, size_t * jmax);

int matrix_equal (const matrix * a, const matrix * b);

int matrix_add (matrix * a, const matrix * b);
int matrix_sub (matrix * a, const matrix * b);
int matrix_mul_elements (matrix * a, const matrix * b);
int matrix_div_elements (matrix * a, const matrix * b);
int matrix_scale (matrix * a, const double x);
int matrix_add_constant (matrix * a, const double x);
int matrix_add_diagonal (matrix * a, const double x);

int matrix_mul (const matrix * a, const matrix * b, matrix * c);

double matrix_get (const matrix * m, const size_t i, const size_t j);
void matrix_set (matrix * m, const size_t i, const size_t j, const double x);
void matrix_print(matrix * m);

#endif

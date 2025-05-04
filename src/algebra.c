#include "algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
    if (a.rows != b.rows || a.cols != b.cols) {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
    
    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[i][j] = a.data[i][j] + b.data[i][j];
        }
    }
    return result;
}

Matrix sub_matrix(Matrix a, Matrix b) {
    if (a.rows != b.rows || a.cols != b.cols) {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
    
    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[i][j] = a.data[i][j] - b.data[i][j];
        }
    }
    return result;
}

Matrix mul_matrix(Matrix a, Matrix b) {
    if (a.cols != b.rows) {
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
        return create_matrix(0, 0);
    }
    
    Matrix result = create_matrix(a.rows, b.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < b.cols; j++) {
            result.data[i][j] = 0;
            for (int k = 0; k < a.cols; k++) {
                result.data[i][j] += a.data[i][k] * b.data[k][j];
            }
        }
    }
    return result;
}

Matrix scale_matrix(Matrix a, double k) {
    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[i][j] = a.data[i][j] * k;
        }
    }
    return result;
}

Matrix transpose_matrix(Matrix a) {
    Matrix result = create_matrix(a.cols, a.rows);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[j][i] = a.data[i][j];
        }
    }
    return result;
}

double det_matrix(Matrix a) {
    if (a.rows != a.cols) {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    
    int n = a.rows;
    if (n == 1) return a.data[0][0];
    
    double det = 0;
    for (int k = 0; k < n; k++) {
        Matrix submat = create_matrix(n-1, n-1);
        for (int i = 1; i < n; i++) {
            int subcol = 0;
            for (int j = 0; j < n; j++) {
                if (j == k) continue;
                submat.data[i-1][subcol++] = a.data[i][j];
            }
        }
        det += (k % 2 == 0 ? 1 : -1) * a.data[0][k] * det_matrix(submat);
    }
    return det;
}

Matrix inv_matrix(Matrix a) {
    if (a.rows != a.cols) {
        printf("Error: The matrix must be a square matrix.\n");
        return create_matrix(0, 0);
    }
    
    double det = det_matrix(a);
    if (fabs(det) < 1e-6) {
        printf("Error: The matrix is singular.\n");
        return create_matrix(0, 0);
    }
    
    int n = a.rows;
    Matrix inv = create_matrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Matrix minor = create_matrix(n-1, n-1);
            for (int x = 0, mi = 0; x < n; x++) {
                if (x == i) continue;
                for (int y = 0, mj = 0; y < n; y++) {
                    if (y == j) continue;
                    minor.data[mi][mj++] = a.data[x][y];
                }
                mi++;
            }
            double cofactor = ((i+j) % 2 == 0 ? 1 : -1) * det_matrix(minor);
            inv.data[j][i] = cofactor / det;
        }
    }
    return inv;
}

int rank_matrix(Matrix a) {
    Matrix mat = a;
    int rank = 0;
    int rows = mat.rows, cols = mat.cols;
    
    for (int col = 0; col < cols && rank < rows; col++) {
        int pivot = rank;
        for (int i = rank; i < rows; i++) {
            if (fabs(mat.data[i][col]) > fabs(mat.data[pivot][col])) {
                pivot = i;
            }
        }
        if (fabs(mat.data[pivot][col]) < 1e-6) continue;
        if (pivot != rank) {
            for (int j = col; j < cols; j++) {
                double temp = mat.data[rank][j];
                mat.data[rank][j] = mat.data[pivot][j];
                mat.data[pivot][j] = temp;
            }
        }
        for (int i = rank + 1; i < rows; i++) {
            double factor = mat.data[i][col] / mat.data[rank][col];
            for (int j = col; j < cols; j++) {
                mat.data[i][j] -= factor * mat.data[rank][j];
            }
        }
        rank++;
    }
    return rank;
}

double trace_matrix(Matrix a) {
    if (a.rows != a.cols) {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    
    double trace = 0;
    for (int i = 0; i < a.rows; i++) {
        trace += a.data[i][i];
    }
    return trace;
}

void print_matrix(Matrix a)
{
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}
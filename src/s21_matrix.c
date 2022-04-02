#include "s21_matrix.h"

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
    int res = SUCCESS;
    if (A->rows == B->rows && A->columns == B->columns &&
        A->matrix_type != INCORRECT_MATRIX &&
        B->matrix_type != INCORRECT_MATRIX) {
        for (int i = 0; i < A->rows && res != FAILURE; i++) {
            for (int j = 0; j < A->columns && res != FAILURE; j++) {
                double temp = fabs(A->matrix[i][j] - B->matrix[i][j]);
                if (temp >= 1e-7) {
                    res = FAILURE;
                }
            }
        }
    } else {
        res = FAILURE;
    }
    return res;
}

matrix_t s21_create_matrix(int rows, int columns) {
    matrix_t matrix;
    if (rows > 0 && columns > 0) {
        matrix.rows = rows;
        matrix.columns = columns;
        matrix.matrix = (double **)malloc(rows * sizeof(double *));
        for (int i = 0; i < rows; i++) {
            matrix.matrix[i] = (double *)calloc(columns, sizeof(double));
        }
        matrix.matrix_type = ZERO_MATRIX;
    } else {
        matrix.rows = 0;
        matrix.columns = 0;
        matrix.matrix_type = INCORRECT_MATRIX;
    }
    return matrix;
}

void s21_remove_matrix(matrix_t *A) {
    if (A->matrix_type != INCORRECT_MATRIX) {
        for (int i = 0; i < A->rows; i++) {
            free(A->matrix[i]);
        }
        free(A->matrix);
        A->rows = 0;
        A->columns = 0;
        A->matrix_type = INCORRECT_MATRIX;
    }
}

matrix_t s21_sum_matrix(matrix_t *A, matrix_t *B) {
    matrix_t matrix;
    if (A->rows == B->rows && A->columns == B->columns &&
        A->matrix_type != INCORRECT_MATRIX &&
        B->matrix_type != INCORRECT_MATRIX) {
        matrix = s21_create_matrix(A->rows, A->columns);
        for (int i = 0; i < matrix.rows; i++) {
            for (int j = 0; j < matrix.columns; j++) {
                matrix.matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
            }
        }
        if (check_matrix_zero(&matrix) == 1) {
            matrix.matrix_type = ZERO_MATRIX;
        } else if (check_matrix_identity(&matrix) == 1) {
            matrix.matrix_type = IDENTITY_MATRIX;
        } else {
            matrix.matrix_type = CORRECT_MATRIX;
        }
    } else {
        matrix.rows = 0;
        matrix.columns = 0;
        matrix.matrix_type = INCORRECT_MATRIX;
    }
    return matrix;
}

matrix_t s21_sub_matrix(matrix_t *A, matrix_t *B) {
    matrix_t matrix;
    if (A->rows == B->rows && A->columns == B->columns &&
        A->matrix_type != INCORRECT_MATRIX &&
        B->matrix_type != INCORRECT_MATRIX) {
        matrix = s21_create_matrix(A->rows, A->columns);
        for (int i = 0; i < matrix.rows; i++) {
            for (int j = 0; j < matrix.columns; j++) {
                matrix.matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
            }
        }
        if (check_matrix_zero(&matrix) == 1) {
            matrix.matrix_type = ZERO_MATRIX;
        } else if (check_matrix_identity(&matrix) == 1) {
            matrix.matrix_type = IDENTITY_MATRIX;
        } else {
            matrix.matrix_type = CORRECT_MATRIX;
        }
    } else {
        matrix = s21_create_matrix(0, 0);
        matrix.matrix_type = INCORRECT_MATRIX;
    }
    return matrix;
}

matrix_t s21_mult_number(matrix_t *A, double number) {
    matrix_t matrix;
    if (A->matrix_type != INCORRECT_MATRIX && number == number) {
        matrix = s21_create_matrix(A->rows, A->columns);
        for (int i = 0; i < matrix.rows; i++) {
            for (int j = 0; j < matrix.columns; j++) {
                matrix.matrix[i][j] = A->matrix[i][j] * number;
            }
        }
        if (check_matrix_zero(&matrix) == 1) {
            matrix.matrix_type = ZERO_MATRIX;
        } else if (check_matrix_identity(&matrix) == 1) {
            matrix.matrix_type = IDENTITY_MATRIX;
        } else {
            matrix.matrix_type = CORRECT_MATRIX;
        }
    } else {
        matrix = s21_create_matrix(0, 0);
        matrix.matrix_type = INCORRECT_MATRIX;
    }
    return matrix;
}

matrix_t s21_mult_matrix(matrix_t *A, matrix_t *B) {
    matrix_t matrix;
    if (A->rows == B->columns && A->matrix_type != INCORRECT_MATRIX &&
        B->matrix_type != INCORRECT_MATRIX) {
        matrix = s21_create_matrix(A->rows, B->columns);
        for (int i = 0; i < matrix.rows; ++i) {
            for (int j = 0; j < matrix.columns; ++j) {
                for (int k = 0; k < A->columns;) {
                    matrix.matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
                    k++;
                }
            }
        }
        if (check_matrix_zero(&matrix) == 1) {
            matrix.matrix_type = ZERO_MATRIX;
        } else if (check_matrix_identity(&matrix) == 1) {
            matrix.matrix_type = IDENTITY_MATRIX;
        } else {
            matrix.matrix_type = CORRECT_MATRIX;
        }
    } else {
        matrix = s21_create_matrix(0, 0);
        matrix.matrix_type = INCORRECT_MATRIX;
    }
    return matrix;
}

matrix_t s21_transpose(matrix_t *A) {
    matrix_t matrix;
    if (A->matrix_type != INCORRECT_MATRIX) {
        matrix = s21_create_matrix(A->columns, A->rows);
        for (int i = 0; i < matrix.rows; ++i) {
            for (int j = 0; j < matrix.columns; ++j) {
                matrix.matrix[i][j] = A->matrix[j][i];
            }
        }
        if (check_matrix_zero(&matrix) == 1) {
            matrix.matrix_type = ZERO_MATRIX;
        } else if (check_matrix_identity(&matrix) == 1) {
            matrix.matrix_type = IDENTITY_MATRIX;
        } else {
            matrix.matrix_type = CORRECT_MATRIX;
        }
    } else {
        matrix = s21_create_matrix(0, 0);
        matrix.matrix_type = INCORRECT_MATRIX;
    }
    return matrix;
}

matrix_t s21_mini_matrix(matrix_t *A, int x, int y) {
    matrix_t matrix;
    int columns = A->columns;
    int rows = A->rows;
    if (A->matrix_type != INCORRECT_MATRIX && A->columns == A->rows) {
        matrix = s21_create_matrix(A->columns - 1, A->rows - 1);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                if (i != x && j != y) {
                    if (i > x && j < y) {
                        matrix.matrix[i - 1][j] = A->matrix[i][j];
                    } else if (i < x && j > y) {
                        matrix.matrix[i][j - 1] = A->matrix[i][j];
                    } else if (i > x && j > y) {
                        matrix.matrix[i - 1][j - 1] = A->matrix[i][j];
                    } else {
                        matrix.matrix[i][j] = A->matrix[i][j];
                    }
                }
            }
        }
        if (check_matrix_zero(&matrix) == 1) {
            matrix.matrix_type = ZERO_MATRIX;
        } else if (check_matrix_identity(&matrix) == 1) {
            matrix.matrix_type = IDENTITY_MATRIX;
        } else {
            matrix.matrix_type = CORRECT_MATRIX;
        }
    } else {
        matrix = s21_create_matrix(0, 0);
        matrix.matrix_type = INCORRECT_MATRIX;
    }
    return matrix;
}

double s21_determinant(matrix_t *A) {
    double result = 0;
    if (eq_sqr(A) == 1 && A->matrix_type != INCORRECT_MATRIX) {
        if (A->columns == 1 && A->rows == 1) {
            result = A->matrix[0][0];
        } else if (A->columns == 2 && A->rows == 2) {
            result = A->matrix[0][0] * A->matrix[1][1] -
                     A->matrix[1][0] * A->matrix[0][1];
        } else if (A->rows > 2 && A->columns > 2) {
            for (int i = 1; i <= A->columns; i++) {
                matrix_t lil = s21_mini_matrix(A, 0, i - 1);
                result += pow(-1, 1 + (double)i) * A->matrix[0][i - 1] *
                          s21_determinant(&lil);

                s21_remove_matrix(&lil);
            }
        }
    } else {
        result = NAN;
    }
    return result;
}

matrix_t s21_calc_complements(matrix_t *A) {
    matrix_t matrix;
    matrix_t temp;
    if (A->matrix_type != INCORRECT_MATRIX && eq_sqr(A) == 1) {
        if (A->rows == 1 && A->columns == 1) {
            matrix = s21_create_matrix(1, 1);
            matrix.matrix[0][0] = 1;
        } else {
            matrix = s21_create_matrix(A->columns, A->rows);
            for (int i = 0; i < A->rows; i++) {
                for (int j = 0; j < A->columns; j++) {
                    temp = s21_mini_matrix(A, i, j);
                    matrix.matrix[i][j] =
                        s21_determinant(&temp) * pow(-1, (i + 1) + (j + 1));
                    s21_remove_matrix(&temp);
                }
            }
        }
        if (check_matrix_zero(&matrix) == 1) {
            matrix.matrix_type = ZERO_MATRIX;
        } else if (check_matrix_identity(&matrix) == 1) {
            matrix.matrix_type = IDENTITY_MATRIX;
        } else {
            matrix.matrix_type = CORRECT_MATRIX;
        }
    } else {
        matrix = s21_create_matrix(0, 0);
        matrix.matrix_type = INCORRECT_MATRIX;
    }
    return matrix;
}

int check_matrix_zero(matrix_t *A) {
    matrix_t matrix = s21_create_matrix(A->columns, A->rows);
    int result = FAILURE;
    if (s21_eq_matrix(A, &matrix) == SUCCESS) {
        A->matrix_type = ZERO_MATRIX;
        result = SUCCESS;
    }
    s21_remove_matrix(&matrix);
    return result;
}

int check_matrix_identity(matrix_t *A) {
    matrix_t matrix = s21_create_matrix(A->columns, A->rows);
    for (int i = 0, j = 0; i < matrix.rows && j < matrix.columns; i++, j++) {
        matrix.matrix[i][j] = 1;
    }
    int result = FAILURE;
    if (s21_eq_matrix(A, &matrix) == SUCCESS) {
        A->matrix_type = IDENTITY_MATRIX;
        result = SUCCESS;
    }
    s21_remove_matrix(&matrix);

    return result;
}

matrix_t s21_inverse_matrix(matrix_t *A) {
    matrix_t matrix;
    if (A->matrix_type != INCORRECT_MATRIX && eq_sqr(A) == 1) {
        double determinant = s21_determinant(A);
        if (determinant != 0) {
            matrix = s21_create_matrix(A->rows, A->columns);
            s21_remove_matrix(&matrix);
            matrix_t temp = s21_calc_complements(A);
            matrix_t transpose = s21_transpose(&temp);
            matrix = s21_mult_number(&transpose, 1 / determinant);
            s21_remove_matrix(&temp);
            s21_remove_matrix(&transpose);
            if (check_matrix_zero(&matrix) == 1) {
                matrix.matrix_type = ZERO_MATRIX;
            } else if (check_matrix_identity(&matrix) == 1) {
                matrix.matrix_type = IDENTITY_MATRIX;
            } else {
                matrix.matrix_type = CORRECT_MATRIX;
            }
        } else {
            matrix.columns = 0;
            matrix.rows = 0;
            matrix.matrix_type = INCORRECT_MATRIX;
        }
    } else {
        matrix = s21_create_matrix(0, 0);
        matrix.matrix_type = INCORRECT_MATRIX;
    }
    return matrix;
}

int eq_sqr(matrix_t *A) {
    int suck = FAILURE;
    if (A->rows == A->columns) {
        suck = SUCCESS;
    }
    return suck;
}

/*
 * DO NOT EDIT THIS FILE. Generated by checkmk.
 * Edit the original source file "s21_matrix-tests.check" instead.
 */

#include <check.h>

#line 1 "s21_matrix-tests.check"
#include <math.h>

#include "s21_matrix.h"

START_TEST(s21_create_matrix_1) {
#line 5
    matrix_t one = s21_create_matrix(4, 4);
    ck_assert_int_eq(one.matrix_type, 3);
    s21_remove_matrix(&one);
}
END_TEST

START_TEST(s21_create_matrix_2) {
#line 11
    matrix_t two = s21_create_matrix(-1, 4);
    ck_assert_int_eq(two.matrix_type, 1);
    s21_remove_matrix(&two);
}
END_TEST

START_TEST(s21_remove_matrix_1) {
#line 17
    matrix_t one = s21_create_matrix(4, 4);
    s21_remove_matrix(&one);
    ck_assert_int_eq(one.matrix_type, 1);
}
END_TEST

START_TEST(s21_remove_matrix_2) {
#line 23
    matrix_t one = s21_create_matrix(0, -1);
    s21_remove_matrix(&one);
    ck_assert_int_eq(one.matrix_type, 1);
}
END_TEST

START_TEST(s21_eq_matrix_1) {
#line 29
    matrix_t one = s21_create_matrix(4, 4);
    matrix_t two = s21_create_matrix(4, 4);
    one.matrix[0][0] = 1;
    one.matrix[0][1] = 0;
    one.matrix[0][2] = 0;
    one.matrix[0][3] = 0;
    one.matrix[1][0] = 0;
    one.matrix[1][1] = 1;
    one.matrix[1][2] = 0;
    one.matrix[1][3] = 0;
    one.matrix[2][0] = 0;
    one.matrix[2][1] = 0;
    one.matrix[2][2] = 1;
    one.matrix[2][3] = 0;
    one.matrix[3][0] = 0;
    one.matrix[3][1] = 0;
    one.matrix[3][2] = 0;
    one.matrix[3][3] = 1;
    int free = s21_eq_matrix(&one, &two);
    ck_assert_int_eq(one.matrix_type, 3);
    ck_assert_int_eq(two.matrix_type, 3);
    ck_assert_int_eq(free, 0);
    s21_remove_matrix(&one);
    s21_remove_matrix(&two);
}
END_TEST

START_TEST(s21_eq_matrix_2) {
#line 56
    matrix_t one = s21_create_matrix(4, 4);
    matrix_t two = s21_create_matrix(4, 4);
    int free = s21_eq_matrix(&one, &two);
    ck_assert_int_eq(one.matrix_type, 3);
    ck_assert_int_eq(two.matrix_type, 3);
    ck_assert_int_eq(free, 1);
    s21_remove_matrix(&one);
    s21_remove_matrix(&two);
}
END_TEST

START_TEST(s21_eq_matrix_3) {
#line 67
    matrix_t one = s21_create_matrix(3, 4);
    matrix_t two = s21_create_matrix(4, 4);
    int free = s21_eq_matrix(&one, &two);
    ck_assert_int_eq(one.matrix_type, 3);
    ck_assert_int_eq(two.matrix_type, 3);
    ck_assert_int_eq(free, 0);
    s21_remove_matrix(&one);
    s21_remove_matrix(&two);
}
END_TEST

START_TEST(s21_sum_matrix_1) {
#line 78
    matrix_t one = s21_create_matrix(4, 4);
    matrix_t two = s21_create_matrix(4, 4);
    one.matrix[0][0] = 1;
    one.matrix[0][1] = 0;
    one.matrix[0][2] = 0;
    one.matrix[0][3] = 0;
    one.matrix[1][0] = 0;
    one.matrix[1][1] = 1;
    one.matrix[1][2] = 0;
    one.matrix[1][3] = 0;
    one.matrix[2][0] = 0;
    one.matrix[2][1] = 0;
    one.matrix[2][2] = 1;
    one.matrix[2][3] = 0;
    one.matrix[3][0] = 0;
    one.matrix[3][1] = 0;
    one.matrix[3][2] = 0;
    one.matrix[3][3] = 1;
    int free = s21_eq_matrix(&one, &two);
    ck_assert_int_eq(one.matrix_type, 3);
    ck_assert_int_eq(two.matrix_type, 3);
    ck_assert_int_eq(free, 0);
    matrix_t four = s21_sum_matrix(&one, &two);
    ck_assert_int_eq(four.matrix_type, 2);
    s21_remove_matrix(&one);
    s21_remove_matrix(&two);
    s21_remove_matrix(&four);
}
END_TEST

START_TEST(s21_sum_matrix_2) {
#line 108
    matrix_t one = s21_create_matrix(4, 4);
    matrix_t two = s21_create_matrix(4, 4);
    int free = s21_eq_matrix(&one, &two);
    ck_assert_int_eq(one.matrix_type, 3);
    ck_assert_int_eq(two.matrix_type, 3);
    ck_assert_int_eq(free, 1);
    matrix_t four = s21_sum_matrix(&one, &two);
    ck_assert_int_eq(four.matrix_type, 3);
    s21_remove_matrix(&one);
    s21_remove_matrix(&two);
    s21_remove_matrix(&four);
}
END_TEST

START_TEST(s21_sum_matrix_3) {
#line 122
    matrix_t one = s21_create_matrix(4, 4);
    matrix_t two = s21_create_matrix(3, 4);
    int free = s21_eq_matrix(&one, &two);
    ck_assert_int_eq(free, 0);
    ck_assert_int_eq(one.matrix_type, 3);
    ck_assert_int_eq(two.matrix_type, 3);
    matrix_t four = s21_sum_matrix(&one, &two);
    ck_assert_int_eq(four.matrix_type, 1);
    s21_remove_matrix(&one);
    s21_remove_matrix(&two);
    s21_remove_matrix(&four);
}
END_TEST

START_TEST(s21_sum_matrix_4) {
#line 136
    matrix_t one = s21_create_matrix(4, 4);
    matrix_t two = s21_create_matrix(4, 4);
    one.matrix[0][0] = 1;
    one.matrix[0][1] = 0;
    one.matrix[0][2] = 3;
    one.matrix[0][3] = 0;
    one.matrix[1][0] = 0;
    one.matrix[1][1] = 1;
    one.matrix[1][2] = 0;
    one.matrix[1][3] = 4;
    one.matrix[2][0] = 0;
    one.matrix[2][1] = 0;
    one.matrix[2][2] = 1;
    one.matrix[2][3] = 0;
    one.matrix[3][0] = 2.123;
    one.matrix[3][1] = 0;
    one.matrix[3][2] = 0;
    one.matrix[3][3] = 1;
    int free = s21_eq_matrix(&one, &two);
    ck_assert_int_eq(one.matrix_type, 3);
    ck_assert_int_eq(two.matrix_type, 3);
    ck_assert_int_eq(free, 0);
    matrix_t four = s21_sum_matrix(&one, &two);
    ck_assert_int_eq(four.matrix_type, 0);
    s21_remove_matrix(&one);
    s21_remove_matrix(&two);
    s21_remove_matrix(&four);
}
END_TEST

START_TEST(s21_sub_matrix_1) {
#line 166
    matrix_t one = s21_create_matrix(4, 4);
    one.matrix[0][0] = 1;
    one.matrix[0][1] = 0;
    one.matrix[0][2] = 3;
    one.matrix[0][3] = 0;
    one.matrix[1][0] = 0;
    one.matrix[1][1] = 1;
    one.matrix[1][2] = 0;
    one.matrix[1][3] = 4;
    one.matrix[2][0] = 0;
    one.matrix[2][1] = 0;
    one.matrix[2][2] = 1;
    one.matrix[2][3] = 0;
    one.matrix[3][0] = 2.123;
    one.matrix[3][1] = 0;
    one.matrix[3][2] = 0;
    one.matrix[3][3] = 1;
    matrix_t two = s21_create_matrix(4, 4);
    two.matrix[0][0] = 1;
    two.matrix[0][1] = 0;
    two.matrix[0][2] = 3;
    two.matrix[0][3] = 0;
    two.matrix[1][0] = 0;
    two.matrix[1][1] = 1;
    two.matrix[1][2] = 0;
    two.matrix[1][3] = 4;
    two.matrix[2][0] = 0;
    two.matrix[2][1] = 0;
    two.matrix[2][2] = 1;
    two.matrix[2][3] = 0;
    two.matrix[3][0] = 2.123;
    two.matrix[3][1] = 0;
    two.matrix[3][2] = 0;
    two.matrix[3][3] = 1;
    int free = s21_eq_matrix(&one, &two);
    ck_assert_int_eq(one.matrix_type, 3);
    ck_assert_int_eq(two.matrix_type, 3);
    ck_assert_int_eq(free, 1);
    matrix_t four = s21_sub_matrix(&one, &two);
    ck_assert_int_eq(four.matrix_type, 3);
    s21_remove_matrix(&one);
    s21_remove_matrix(&two);
    s21_remove_matrix(&four);
}
END_TEST

START_TEST(s21_sub_matrix_2) {
#line 212
    matrix_t one = s21_create_matrix(4, 4);
    one.matrix[0][0] = 1;
    one.matrix[0][1] = 32;
    one.matrix[0][2] = 3;
    one.matrix[0][3] = 0;
    one.matrix[1][0] = 112;
    one.matrix[1][1] = 1;
    one.matrix[1][2] = 0;
    one.matrix[1][3] = 4;
    one.matrix[2][0] = 312;
    one.matrix[2][1] = 112;
    one.matrix[2][2] = 1;
    one.matrix[2][3] = 0.123;
    one.matrix[3][0] = 2.123;
    one.matrix[3][1] = 0;
    one.matrix[3][2] = 0.123123;
    one.matrix[3][3] = 1;
    matrix_t two = s21_create_matrix(4, 4);
    two.matrix[0][0] = 1;
    two.matrix[0][1] = 0;
    two.matrix[0][2] = 3;
    two.matrix[0][3] = 0;
    two.matrix[1][0] = 0;
    two.matrix[1][1] = 1;
    two.matrix[1][2] = 0;
    two.matrix[1][3] = 4;
    two.matrix[2][0] = 0;
    two.matrix[2][1] = 0;
    two.matrix[2][2] = 1;
    two.matrix[2][3] = 0;
    two.matrix[3][0] = 2.123;
    two.matrix[3][1] = 0;
    two.matrix[3][2] = 0;
    two.matrix[3][3] = 1;
    int free = s21_eq_matrix(&one, &two);
    ck_assert_int_eq(one.matrix_type, 3);
    ck_assert_int_eq(two.matrix_type, 3);
    ck_assert_int_eq(free, 0);
    matrix_t four = s21_sub_matrix(&one, &two);
    ck_assert_int_eq(four.matrix_type, 0);
    s21_remove_matrix(&one);
    s21_remove_matrix(&two);
    s21_remove_matrix(&four);
}
END_TEST

START_TEST(s21_sub_matrix_3) {
#line 257
    matrix_t one = s21_create_matrix(3, 4);
    matrix_t two = s21_create_matrix(4, 4);
    two.matrix[0][0] = 1;
    two.matrix[0][1] = 0;
    two.matrix[0][2] = 3;
    two.matrix[0][3] = 0;
    two.matrix[1][0] = 0;
    two.matrix[1][1] = 1;
    two.matrix[1][2] = 0;
    two.matrix[1][3] = 4;
    two.matrix[2][0] = 0;
    two.matrix[2][1] = 0;
    two.matrix[2][2] = 1;
    two.matrix[2][3] = 0;
    two.matrix[3][0] = 2.123;
    two.matrix[3][1] = 0;
    two.matrix[3][2] = 0;
    two.matrix[3][3] = 1;
    int free = s21_eq_matrix(&one, &two);
    ck_assert_int_eq(one.matrix_type, 3);
    ck_assert_int_eq(two.matrix_type, 3);
    ck_assert_int_eq(free, 0);
    matrix_t four = s21_sub_matrix(&one, &two);
    ck_assert_int_eq(four.matrix_type, 1);
    s21_remove_matrix(&one);
    s21_remove_matrix(&two);
    s21_remove_matrix(&four);
}
END_TEST

START_TEST(s21_sub_matrix_4) {
#line 287
    matrix_t one = s21_create_matrix(4, 4);
    one.matrix[0][0] = 2;
    one.matrix[0][1] = 0;
    one.matrix[0][2] = 0;
    one.matrix[0][3] = 0;
    one.matrix[1][0] = 0;
    one.matrix[1][1] = 2;
    one.matrix[1][2] = 0;
    one.matrix[1][3] = 0;
    one.matrix[2][0] = 0;
    one.matrix[2][1] = 0;
    one.matrix[2][2] = 2;
    one.matrix[2][3] = 0;
    one.matrix[3][0] = 0;
    one.matrix[3][1] = 0;
    one.matrix[3][2] = 0;
    one.matrix[3][3] = 2;
    matrix_t two = s21_create_matrix(4, 4);
    two.matrix[0][0] = 1;
    two.matrix[0][1] = 0;
    two.matrix[0][2] = 0;
    two.matrix[0][3] = 0;
    two.matrix[1][0] = 0;
    two.matrix[1][1] = 1;
    two.matrix[1][2] = 0;
    two.matrix[1][3] = 0;
    two.matrix[2][0] = 0;
    two.matrix[2][1] = 0;
    two.matrix[2][2] = 1;
    two.matrix[2][3] = 0;
    two.matrix[3][0] = 0;
    two.matrix[3][1] = 0;
    two.matrix[3][2] = 0;
    two.matrix[3][3] = 1;
    int free = s21_eq_matrix(&one, &two);
    ck_assert_int_eq(one.matrix_type, ZERO_MATRIX);
    ck_assert_int_eq(two.matrix_type, ZERO_MATRIX);
    ck_assert_int_eq(free, FAILURE);
    matrix_t four = s21_sub_matrix(&one, &two);
    ck_assert_int_eq(four.matrix_type, IDENTITY_MATRIX);
    s21_remove_matrix(&one);
    s21_remove_matrix(&two);
    s21_remove_matrix(&four);
}
END_TEST

START_TEST(s21_mult_number_1) {
#line 332
    matrix_t one = s21_create_matrix(4, 4);
    one.matrix[0][0] = 1;
    one.matrix[0][1] = 32;
    one.matrix[0][2] = 3;
    one.matrix[0][3] = 0;
    one.matrix[1][0] = 112;
    one.matrix[1][1] = 1;
    one.matrix[1][2] = 0;
    one.matrix[1][3] = 4;
    one.matrix[2][0] = 312;
    one.matrix[2][1] = 112;
    one.matrix[2][2] = 1;
    one.matrix[2][3] = 0.123;
    one.matrix[3][0] = 2.123;
    one.matrix[3][1] = 0;
    one.matrix[3][2] = 0.123123;
    one.matrix[3][3] = 1;
    ck_assert_int_eq(one.matrix_type, ZERO_MATRIX);
    matrix_t two = s21_mult_number(&one, 12);
    ck_assert_int_eq(two.matrix_type, CORRECT_MATRIX);
    s21_remove_matrix(&one);
    s21_remove_matrix(&two);
}
END_TEST

START_TEST(s21_mult_number_2) {
#line 356
    matrix_t one = s21_create_matrix(4, 4);
    one.matrix[0][0] = 1;
    one.matrix[0][1] = 32;
    one.matrix[0][2] = 3;
    one.matrix[0][3] = 0;
    one.matrix[1][0] = 112;
    one.matrix[1][1] = 1;
    one.matrix[1][2] = 0;
    one.matrix[1][3] = 4;
    one.matrix[2][0] = 312;
    one.matrix[2][1] = 112;
    one.matrix[2][2] = 1;
    one.matrix[2][3] = 0.123;
    one.matrix[3][0] = 2.123;
    one.matrix[3][1] = 0;
    one.matrix[3][2] = 0.123123;
    one.matrix[3][3] = 1;
    ck_assert_int_eq(one.matrix_type, ZERO_MATRIX);
    matrix_t two = s21_mult_number(&one, NAN);
    ck_assert_int_eq(two.matrix_type, INCORRECT_MATRIX);
    s21_remove_matrix(&one);
    s21_remove_matrix(&two);
}
END_TEST

START_TEST(s21_mult_number_3) {
#line 380
    matrix_t one = s21_create_matrix(4, 4);
    one.matrix[0][0] = 1;
    one.matrix[0][1] = 32;
    one.matrix[0][2] = 3;
    one.matrix[0][3] = 0;
    one.matrix[1][0] = 112;
    one.matrix[1][1] = 1;
    one.matrix[1][2] = 0;
    one.matrix[1][3] = 4;
    one.matrix[2][0] = 312;
    one.matrix[2][1] = 112;
    one.matrix[2][2] = 1;
    one.matrix[2][3] = 0.123;
    one.matrix[3][0] = 2.123;
    one.matrix[3][1] = 0;
    one.matrix[3][2] = 0.123123;
    one.matrix[3][3] = 1;
    ck_assert_int_eq(one.matrix_type, ZERO_MATRIX);
    matrix_t two = s21_mult_number(&one, 0);
    ck_assert_int_eq(two.matrix_type, ZERO_MATRIX);
    s21_remove_matrix(&one);
    s21_remove_matrix(&two);
}
END_TEST

START_TEST(s21_mult_matrix_1) {
#line 404
    int rows = 3;
    int columns = 3;
    matrix_t matrix = s21_create_matrix(rows, columns);
    matrix.matrix[0][0] = 1;
    matrix.matrix[0][1] = 2;
    matrix.matrix[0][2] = 0;
    matrix.matrix[1][0] = 0;
    matrix.matrix[1][1] = 4;
    matrix.matrix[1][2] = 2.55;
    matrix.matrix[2][0] = 5;
    matrix.matrix[2][1] = 2;
    matrix.matrix[2][2] = 1;

    int rows1 = 3;
    int columns1 = 3;
    matrix_t matrix2 = s21_create_matrix(rows1, columns1);
    matrix2.matrix[0][0] = -1;
    matrix2.matrix[0][1] = 2;
    matrix2.matrix[0][2] = -3;
    matrix2.matrix[1][0] = 0;
    matrix2.matrix[1][1] = 4;
    matrix2.matrix[1][2] = 2.56;
    matrix2.matrix[2][0] = 5;
    matrix2.matrix[2][1] = 0.1231;
    matrix2.matrix[2][2] = 1;

    int rows2 = 3;
    int columns2 = 3;
    matrix_t matrix_pr = s21_create_matrix(rows2, columns2);
    matrix_pr.matrix[0][0] = -1;
    matrix_pr.matrix[0][1] = 10;
    matrix_pr.matrix[0][2] = 2.12;
    matrix_pr.matrix[1][0] = 12.75;
    matrix_pr.matrix[1][1] = 16.313905;
    matrix_pr.matrix[1][2] = 12.79;
    matrix_pr.matrix[2][0] = 0;
    matrix_pr.matrix[2][1] = 18.1231;
    matrix_pr.matrix[2][2] = -8.88;

    matrix_t matrix_res = s21_mult_matrix(&matrix, &matrix2);

    ck_assert_int_eq(s21_eq_matrix(&matrix_res, &matrix_pr), 1);

    s21_remove_matrix(&matrix);
    s21_remove_matrix(&matrix2);
    s21_remove_matrix(&matrix_pr);
    s21_remove_matrix(&matrix_res);
}
END_TEST

START_TEST(s21_mult_matrix_2) {
#line 454
    int rows = 3;
    int columns = 3;
    matrix_t matrix = s21_create_matrix(rows, columns);
    matrix.matrix[0][0] = 0;
    matrix.matrix[0][1] = 0;
    matrix.matrix[0][2] = 0;
    matrix.matrix[1][0] = 0;
    matrix.matrix[1][1] = 0;
    matrix.matrix[1][2] = 0;
    matrix.matrix[2][0] = 0;
    matrix.matrix[2][1] = 0;
    matrix.matrix[2][2] = 0;

    int rows1 = 3;
    int columns1 = 3;
    matrix_t matrix2 = s21_create_matrix(rows1, columns1);
    matrix2.matrix[0][0] = 0;
    matrix2.matrix[0][1] = 0;
    matrix2.matrix[0][2] = 0;
    matrix2.matrix[1][0] = 0;
    matrix2.matrix[1][1] = 0;
    matrix2.matrix[1][2] = 0;
    matrix2.matrix[2][0] = 0;
    matrix2.matrix[2][1] = 0;
    matrix2.matrix[2][2] = 0;

    int rows2 = 3;
    int columns2 = 3;
    matrix_t matrix_pr = s21_create_matrix(rows2, columns2);
    matrix_pr.matrix[0][0] = 0;
    matrix_pr.matrix[0][1] = 0;
    matrix_pr.matrix[0][2] = 0;
    matrix_pr.matrix[1][0] = 0;
    matrix_pr.matrix[1][1] = 0;
    matrix_pr.matrix[1][2] = 0;
    matrix_pr.matrix[2][0] = 0;
    matrix_pr.matrix[2][1] = 0;
    matrix_pr.matrix[2][2] = 0;

    matrix_t matrix_res = s21_mult_matrix(&matrix, &matrix2);

    ck_assert_int_eq(s21_eq_matrix(&matrix_res, &matrix_pr), 1);

    s21_remove_matrix(&matrix);
    s21_remove_matrix(&matrix2);
    s21_remove_matrix(&matrix_pr);
    s21_remove_matrix(&matrix_res);
}
END_TEST

START_TEST(s21_mult_matrix_3) {
#line 503
    int rows = 3;
    int columns = 3;
    matrix_t matrix = s21_create_matrix(rows, columns);
    matrix.matrix[0][0] = 1;
    matrix.matrix[0][1] = 0;
    matrix.matrix[0][2] = 0;
    matrix.matrix[1][0] = 0;
    matrix.matrix[1][1] = 1;
    matrix.matrix[1][2] = 0;
    matrix.matrix[2][0] = 0;
    matrix.matrix[2][1] = 0;
    matrix.matrix[2][2] = 1;

    int rows1 = 3;
    int columns1 = 3;
    matrix_t matrix2 = s21_create_matrix(rows1, columns1);
    matrix2.matrix[0][0] = 1;
    matrix2.matrix[0][1] = 0;
    matrix2.matrix[0][2] = 0;
    matrix2.matrix[1][0] = 0;
    matrix2.matrix[1][1] = 1;
    matrix2.matrix[1][2] = 0;
    matrix2.matrix[2][0] = 0;
    matrix2.matrix[2][1] = 0;
    matrix2.matrix[2][2] = 1;

    int rows2 = 3;
    int columns2 = 3;
    matrix_t matrix_pr = s21_create_matrix(rows2, columns2);
    matrix_pr.matrix[0][0] = 1;
    matrix_pr.matrix[0][1] = 0;
    matrix_pr.matrix[0][2] = 0;
    matrix_pr.matrix[1][0] = 0;
    matrix_pr.matrix[1][1] = 1;
    matrix_pr.matrix[1][2] = 0;
    matrix_pr.matrix[2][0] = 0;
    matrix_pr.matrix[2][1] = 0;
    matrix_pr.matrix[2][2] = 1;

    matrix_t matrix_res = s21_mult_matrix(&matrix, &matrix2);

    ck_assert_int_eq(s21_eq_matrix(&matrix_res, &matrix_pr), 1);

    s21_remove_matrix(&matrix);
    s21_remove_matrix(&matrix2);
    s21_remove_matrix(&matrix_pr);
    s21_remove_matrix(&matrix_res);
}
END_TEST

START_TEST(s21_transpose_1) {
#line 552
    matrix_t one = s21_create_matrix(4, 4);
    one.matrix[0][0] = 1;
    one.matrix[0][1] = 32;
    one.matrix[0][2] = 3;
    one.matrix[0][3] = 0;
    one.matrix[1][0] = 112;
    one.matrix[1][1] = 1;
    one.matrix[1][2] = 0;
    one.matrix[1][3] = 4;
    one.matrix[2][0] = 312;
    one.matrix[2][1] = 112;
    one.matrix[2][2] = 1;
    one.matrix[2][3] = 0.123;
    one.matrix[3][0] = 2.123;
    one.matrix[3][1] = 0;
    one.matrix[3][2] = 0.123123;
    one.matrix[3][3] = 1;
    ck_assert_int_eq(one.matrix_type, ZERO_MATRIX);
    matrix_t two = s21_transpose(&one);
    ck_assert_int_eq(two.matrix_type, CORRECT_MATRIX);
    s21_remove_matrix(&one);
    s21_remove_matrix(&two);
}
END_TEST

START_TEST(s21_transpose_2) {
#line 576
    matrix_t one;
    one.rows = 4;
    one.columns = 6;
    one.matrix_type = INCORRECT_MATRIX;
    ck_assert_int_eq(one.matrix_type, INCORRECT_MATRIX);
    matrix_t two = s21_transpose(&one);
    ck_assert_int_eq(two.matrix_type, INCORRECT_MATRIX);
    s21_remove_matrix(&one);
    s21_remove_matrix(&two);
}
END_TEST

START_TEST(s21_transpose_3) {
#line 587
    matrix_t one = s21_create_matrix(4, 4);
    one.matrix[0][0] = 0;
    one.matrix[0][1] = 0;
    one.matrix[0][2] = 0;
    one.matrix[0][3] = 0;
    one.matrix[1][0] = 0;
    one.matrix[1][1] = 0;
    one.matrix[1][2] = 0;
    one.matrix[1][3] = 0;
    one.matrix[2][0] = 0;
    one.matrix[2][1] = 0;
    one.matrix[2][2] = 0;
    one.matrix[2][3] = 0;
    one.matrix[3][0] = 0;
    one.matrix[3][1] = 0;
    one.matrix[3][2] = 0;
    one.matrix[3][3] = 0;
    ck_assert_int_eq(one.matrix_type, ZERO_MATRIX);
    matrix_t two = s21_transpose(&one);
    ck_assert_int_eq(two.matrix_type, ZERO_MATRIX);
    s21_remove_matrix(&one);
    s21_remove_matrix(&two);
}
END_TEST

START_TEST(s21_calc_complements_1) {
#line 611
    matrix_t one = s21_create_matrix(3, 4);
    one.matrix[0][0] = 1;
    one.matrix[0][1] = 32;
    one.matrix[0][2] = 3;
    one.matrix[0][3] = 0;
    one.matrix[1][0] = 112;
    one.matrix[1][1] = 1;
    one.matrix[1][2] = 0;
    one.matrix[1][3] = 4;
    one.matrix[2][0] = 312;
    one.matrix[2][1] = 112;
    one.matrix[2][2] = 1;
    one.matrix[2][3] = 0.123;
    ck_assert_int_eq(one.matrix_type, ZERO_MATRIX);
    matrix_t two = s21_calc_complements(&one);
    ck_assert_int_eq(two.matrix_type, INCORRECT_MATRIX);
    s21_remove_matrix(&one);
    s21_remove_matrix(&two);
}
END_TEST

START_TEST(s21_calc_complements_2) {
#line 631
    matrix_t one = s21_create_matrix(4, 4);
    one.matrix[0][0] = 1;
    one.matrix[0][1] = 32;
    one.matrix[0][2] = 3;
    one.matrix[0][3] = 0;
    one.matrix[1][0] = 112;
    one.matrix[1][1] = 1;
    one.matrix[1][2] = 0;
    one.matrix[1][3] = 4;
    one.matrix[2][0] = 312;
    one.matrix[2][1] = 112;
    one.matrix[2][2] = 1;
    one.matrix[2][3] = 0.123;
    ck_assert_int_eq(one.matrix_type, ZERO_MATRIX);
    matrix_t two = s21_calc_complements(&one);
    ck_assert_int_eq(two.matrix_type, CORRECT_MATRIX);
    s21_remove_matrix(&one);
    s21_remove_matrix(&two);
}
END_TEST

START_TEST(s21_calc_complements_3) {
#line 651
    matrix_t one = s21_create_matrix(1, 1);
    one.matrix[0][0] = 0.23;
    ck_assert_int_eq(one.matrix_type, ZERO_MATRIX);
    matrix_t two = s21_calc_complements(&one);
    ck_assert_int_eq(two.matrix_type, IDENTITY_MATRIX);
    s21_remove_matrix(&one);
    s21_remove_matrix(&two);
}
END_TEST

START_TEST(s21_calc_complements_4) {
#line 660
    matrix_t one = s21_create_matrix(4, 4);
    ck_assert_int_eq(one.matrix_type, ZERO_MATRIX);
    matrix_t two = s21_calc_complements(&one);
    ck_assert_int_eq(two.matrix_type, ZERO_MATRIX);
    s21_remove_matrix(&one);
    s21_remove_matrix(&two);
}
END_TEST

START_TEST(s21_calc_complements_5) {
#line 668
    matrix_t one = s21_create_matrix(3, 3);
    one.matrix[0][0] = 1;
    one.matrix[0][1] = 2;
    one.matrix[0][2] = 3;
    one.matrix[1][0] = 0;
    one.matrix[1][1] = 4;
    one.matrix[1][2] = 2;
    one.matrix[2][0] = 5;
    one.matrix[2][1] = 2;
    one.matrix[2][2] = 1;
    matrix_t two = s21_create_matrix(3, 3);
    two.matrix[0][0] = 0;
    two.matrix[0][1] = 10;
    two.matrix[0][2] = -20;
    two.matrix[1][0] = 4;
    two.matrix[1][1] = -14;
    two.matrix[1][2] = 8;
    two.matrix[2][0] = -8;
    two.matrix[2][1] = -2;
    two.matrix[2][2] = 4;
    ck_assert_int_eq(one.matrix_type, ZERO_MATRIX);
    ck_assert_int_eq(two.matrix_type, ZERO_MATRIX);
    matrix_t three = s21_calc_complements(&one);
    for (int i = 0; i < three.rows; ++i) {
        for (int j = 0; j < three.columns; ++j) {
            ck_assert_double_eq(two.matrix[i][j], three.matrix[i][j]);
        }
    }
    ck_assert_int_eq(three.matrix_type, CORRECT_MATRIX);
    s21_remove_matrix(&one);
    s21_remove_matrix(&two);
    s21_remove_matrix(&three);
}
END_TEST

START_TEST(s21_determinant_1) {
#line 702
    matrix_t one = s21_create_matrix(3, 4);
    ck_assert_int_eq(one.matrix_type, ZERO_MATRIX);
    double two = s21_determinant(&one);
    ck_assert_double_nan(two);
    s21_remove_matrix(&one);
}
END_TEST

START_TEST(s21_determinant_2) {
#line 709
    matrix_t one = s21_create_matrix(3, 3);
    one.matrix[0][0] = 1;
    one.matrix[0][1] = 32;
    one.matrix[0][2] = 3;
    one.matrix[1][0] = 112;
    one.matrix[1][1] = 1;
    one.matrix[1][2] = 0;
    one.matrix[2][0] = 312;
    one.matrix[2][1] = 112;
    one.matrix[2][2] = 1;
    ck_assert_int_eq(one.matrix_type, ZERO_MATRIX);
    double two = s21_determinant(&one);
    ck_assert_double_eq(two, 33113);
    s21_remove_matrix(&one);
}
END_TEST

START_TEST(s21_determinant_3) {
#line 725
    matrix_t one = s21_create_matrix(3, 3);
    one.matrix[0][0] = 1;
    one.matrix[0][1] = 2;
    one.matrix[0][2] = 3;
    one.matrix[1][0] = 4;
    one.matrix[1][1] = 5;
    one.matrix[1][2] = 6;
    one.matrix[2][0] = 7;
    one.matrix[2][1] = 8;
    one.matrix[2][2] = 9;
    ck_assert_int_eq(one.matrix_type, ZERO_MATRIX);
    double two = s21_determinant(&one);
    ck_assert_double_eq(two, 0);
    s21_remove_matrix(&one);
}
END_TEST

START_TEST(s21_inverse_matrix_1) {
#line 741
    matrix_t one = s21_create_matrix(3, 3);
    one.matrix[0][0] = 1;
    one.matrix[0][1] = 2;
    one.matrix[0][2] = 3;
    one.matrix[1][0] = 4;
    one.matrix[1][1] = 5;
    one.matrix[1][2] = 6;
    one.matrix[2][0] = 7;
    one.matrix[2][1] = 8;
    one.matrix[2][2] = 9;
    ck_assert_int_eq(one.matrix_type, ZERO_MATRIX);
    double two = s21_determinant(&one);
    ck_assert_double_eq(two, 0);
    matrix_t three = s21_inverse_matrix(&one);
    ck_assert_int_eq(three.matrix_type, INCORRECT_MATRIX);
    s21_remove_matrix(&one);
    s21_remove_matrix(&three);
}
END_TEST

START_TEST(s21_inverse_matrix_2) {
#line 760
    matrix_t one = s21_create_matrix(2, 3);
    ck_assert_int_eq(one.matrix_type, ZERO_MATRIX);
    double two = s21_determinant(&one);
    ck_assert_double_nan(two);
    matrix_t three = s21_inverse_matrix(&one);
    ck_assert_int_eq(three.matrix_type, INCORRECT_MATRIX);
    s21_remove_matrix(&one);
    s21_remove_matrix(&three);
}
END_TEST

START_TEST(s21_inverse_matrix_3) {
#line 770
    matrix_t one = s21_create_matrix(3, 3);
    one.matrix[0][0] = 5;
    one.matrix[0][1] = -6;
    one.matrix[0][2] = 12;
    one.matrix[1][0] = 2;
    one.matrix[1][1] = 3;
    one.matrix[1][2] = -2;
    one.matrix[2][0] = 1;
    one.matrix[2][1] = 3;
    one.matrix[2][2] = 9;
    ck_assert_int_eq(one.matrix_type, ZERO_MATRIX);
    double two = s21_determinant(&one);
    ck_assert_double_eq(two, 321);
    matrix_t three = s21_inverse_matrix(&one);
    ck_assert_int_eq(three.matrix_type, CORRECT_MATRIX);
    s21_remove_matrix(&one);
    s21_remove_matrix(&three);
}
END_TEST

int main(void) {
    Suite *s1 = suite_create("Core");
    TCase *tc1_1 = tcase_create("Core");
    SRunner *sr = srunner_create(s1);
    int nf;

    suite_add_tcase(s1, tc1_1);
    tcase_add_test(tc1_1, s21_create_matrix_1);
    tcase_add_test(tc1_1, s21_create_matrix_2);
    tcase_add_test(tc1_1, s21_remove_matrix_1);
    tcase_add_test(tc1_1, s21_remove_matrix_2);
    tcase_add_test(tc1_1, s21_eq_matrix_1);
    tcase_add_test(tc1_1, s21_eq_matrix_2);
    tcase_add_test(tc1_1, s21_eq_matrix_3);
    tcase_add_test(tc1_1, s21_sum_matrix_1);
    tcase_add_test(tc1_1, s21_sum_matrix_2);
    tcase_add_test(tc1_1, s21_sum_matrix_3);
    tcase_add_test(tc1_1, s21_sum_matrix_4);
    tcase_add_test(tc1_1, s21_sub_matrix_1);
    tcase_add_test(tc1_1, s21_sub_matrix_2);
    tcase_add_test(tc1_1, s21_sub_matrix_3);
    tcase_add_test(tc1_1, s21_sub_matrix_4);
    tcase_add_test(tc1_1, s21_mult_number_1);
    tcase_add_test(tc1_1, s21_mult_number_2);
    tcase_add_test(tc1_1, s21_mult_number_3);
    tcase_add_test(tc1_1, s21_mult_matrix_1);
    tcase_add_test(tc1_1, s21_mult_matrix_2);
    tcase_add_test(tc1_1, s21_mult_matrix_3);
    tcase_add_test(tc1_1, s21_transpose_1);
    tcase_add_test(tc1_1, s21_transpose_2);
    tcase_add_test(tc1_1, s21_transpose_3);
    tcase_add_test(tc1_1, s21_calc_complements_1);
    tcase_add_test(tc1_1, s21_calc_complements_2);
    tcase_add_test(tc1_1, s21_calc_complements_3);
    tcase_add_test(tc1_1, s21_calc_complements_4);
    tcase_add_test(tc1_1, s21_calc_complements_5);
    tcase_add_test(tc1_1, s21_determinant_1);
    tcase_add_test(tc1_1, s21_determinant_2);
    tcase_add_test(tc1_1, s21_determinant_3);
    tcase_add_test(tc1_1, s21_inverse_matrix_1);
    tcase_add_test(tc1_1, s21_inverse_matrix_2);
    tcase_add_test(tc1_1, s21_inverse_matrix_3);

    srunner_run_all(sr, CK_ENV);
    nf = srunner_ntests_failed(sr);
    srunner_free(sr);

    return nf == 0 ? 0 : 1;
}

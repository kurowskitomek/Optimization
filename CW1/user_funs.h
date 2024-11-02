#pragma once

#include"ode_solver.h"

matrix ff2T(matrix x1, matrix x2, matrix = NAN);
matrix ff1T(matrix, matrix = NAN, matrix = NAN);
matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);
matrix ff1R(matrix x, matrix ud1, matrix ud2);
matrix df1(double t, matrix Y, matrix ud1, matrix ud2);


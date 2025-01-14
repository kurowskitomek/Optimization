#include"user_funs.h"
#include <cmath>

matrix ff1T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = -cos(0.1 * x(0)) * exp(-pow(0.1 * x(0) - 2 * 3.14159265, 2)) + 0.002 * pow((0.1 * x), 2);

	return y;
}

matrix ff2T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * 3.14159265 * x(0)) - cos(2.5 * 3.14159265 * x(1)) + 2;

	return y;
}

matrix gradient_example(matrix x, matrix ud1, matrix ud2)
{
	// Let's assume b is a predefined vector
	matrix b(2, 1);  // 2D column vector as an example
	b(0, 0) = 1.0;
	b(1, 0) = 2.0;

	// Calculate gradient: A * x - b (assuming A is the identity matrix)
	return x - b;  // Since A = I, gradient is just x - b
}

matrix ff3T_outside(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = (sin(3.14 * sqrt(pow(x(0) / 3.14, 2) + pow(x(1) / 3.14, 2)))) / (3.14 * sqrt(pow(x(0) / 3.14, 2) + pow(x(1) / 3.14, 2)));

	double a = m2d(ud1);
	double c = m2d(ud2);

	//g1
	if (-x(0) + 1 > 0)
		y = y + c * pow(-x(0) + 1, 2);
	//g2
	if (-x(1) + 1 > 0)
		y = y + c * pow(-x(1) + 1, 2);
	//g3
	if (norm(x) - a > 0)
		y = y + c * pow(norm(x) - a, 2);

	return y;
}

matrix ff3T_inside(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = (sin(3.14 * sqrt(pow(x(0) / 3.14, 2) + pow(x(1) / 3.14, 2)))) / (3.14 * sqrt(pow(x(0) / 3.14, 2) + pow(x(1) / 3.14, 2)));

	double a = m2d(ud1);
	double c = m2d(ud2);

	//g1
	if (-x(0) + 1 > 0)
		y = 1E10;
	else
		y = y - c / (-x(0) + 1);
	//g2
	if (-x(1) + 1 > 0)
		y = 1E10;
	else
		y = y - c / (-x(1) + 1);
	//g3
	if (norm(x) - a > 0)
		y = 1E10;
	else
		y = y - c / (norm(x) - a);

	return y;
}

matrix ff4T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0) + 2 * x(1) - 7, 2) + pow(2 * x(0) + x(1) - 5, 2);

	return y;
}

matrix df3(double t, matrix Y, matrix ud1, matrix ud2)
{
	//Wektor zmian po czasie
	matrix dY(4, 1);

	//Zmienne dane
	double x = Y(0);
	double v_x = Y(1);
	double y = Y(2);
	double v_y = Y(3);

	//Dane zadania
	double C = ud1(0);
	double rho = ud1(1);
	double r = ud1(2);
	double m = ud1(3);
	double g = ud1(4);

	double s = 3.14 * pow(r, 2);

	double D_x = (1.0 / 2.0) * C * rho * s * v_x * abs(v_x);
	double D_y = (1.0 / 2.0) * C * rho * s * v_y * abs(v_y);
	double FM_x = rho * v_y * m2d(ud2) * 3.14 * pow(r, 3);
	double FM_y = rho * v_x * m2d(ud2) * 3.14 * pow(r, 3);

	dY(0) = v_x;
	dY(1) = (-D_x - FM_x) / m;
	dY(2) = v_y;
	dY(3) = ((-m * g) - D_y - FM_y) / m;

	return dY;
}

matrix ff3R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0(4, new double[4] { 0.0, x(0), 100, 0 });
	matrix* Y = solve_ode(df3, 0.0, 0.01, 7.0, Y0, ud1, x(1));
	int n = get_len(Y[0]);
	int i50 = 0;
	int i_0 = 0;
	for (int i = 0; i < n; ++i)
	{
		if (abs(Y[1](i, 2) - 50.0) < abs(Y[1](i50, 2) - 50.0))
			i50 = i;
		if (abs(Y[1](i, 2)) < abs(Y[1](i_0, 2)))
			i_0 = i;
	}

	y = -Y[1](i_0, 0);

	if (abs(x(0)) - 10 > 0)
		y = y + ud2 * pow(abs(x(0)) - 10, 2);
	if (abs(x(1)) - 15 > 0)
		y = y + ud2 * pow(abs(x(1)) - 15, 2);
	if (abs(Y[1](i50, 0) - 5.0) - 0.5 > 0)
		y = y + ud2 * pow(abs(Y[1](i50, 0) - 5.0) - 0.5, 2);

	return y;
}

matrix ff1R(matrix x, matrix ud1, matrix ud2)
{
	matrix Y0 = matrix(3, new double[3] { 5, 1, 20 });
	matrix* y = solve_ode(df1, 0, 1, 2000, Y0, ud1, x);
	int n = get_len(y[0]);
	double max = y[1](0, 2);
	for (int i = 0; i < n; i++)
		if (max < y[1](i, 2))
			max = y[1](i, 2);

	return abs(max - 50);
}

matrix ff2R(matrix x, matrix ud1, matrix ud2)
{
	matrix y(1);
	y(0) = 0;
	//matrix x
	matrix Y0(2, 1), Yref(2, new double[2] { 3.14, 0 });
	//cout << "bbbb" << endl;
	matrix* Y = solve_ode(df2, 0, 0.1, 100, Y0, Yref, x);
	//cout << "bbbb" << endl;
	int n = get_len(Y[0]);
	for (int i = 0; i < n; i++)
	{
		y = y + 10.0 * pow(Yref(0) - Y[1](i, 0), 2)
			+ pow(Yref(1) - Y[1](i, 0), 2)
			+ pow(x(0) * (Yref(0) - Y[1](i, 0)) + x(1) * (Yref(1) - Y[1](i, 0)), 2);
		y = 0.1 * y;
	}
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

matrix df1(double t, matrix Y, matrix ud1, matrix ud2)
{
	double a = 0.98, b = 0.63, g = 9.81, PA = 0.5, TA = 90, PB = 1, DB = 0.00365665, Fin = 0.01, Tin = 20;

	matrix dY(3, 1);
	double FAout = Y(0) > 0 ? a * b * m2d(ud2) * sqrt(2 * g * Y(0) / PA) : 0;
	double FBout = Y(1) > 0 ? a * b * DB * sqrt(2 * g * Y(1) / PB) : 0;

	dY(0) = -FAout;
	dY(1) = FAout + Fin - FBout;
	dY(2) = Fin / Y(1) * (Tin - Y(2)) + FAout / Y(1) * (TA - Y(2));
	// cout << Y(0) - dY(0) << " " << Y(1) - dY(1) << " " << Y(2) - dY(2) << endl;
	return dY;
}

matrix df2(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);

	dY(0) = Y(1);
	double l = 1;
	double mr = 1;
	double mc = 5;
	double I = mr * l * l / 3 + mc * l * l;
	double b = 0.5;
	double Mt = ud2(0) * (ud1(0) - Y(0)) + ud2(1) * (ud1(1) - Y(1));

	dY(1) = (Mt - b * Y(1)) / I;

	return dY;
}

matrix ff0T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
	int n = get_len(Y[0]);
	double teta_max = Y[1](0, 0);
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);
	double m = 1, l = 0.5, b = 0.5, g = 9.81;
	double I = m*pow(l, 2);
	dY(0) = Y(1);
	dY(1) = ((t <= ud2(1))*ud2(0) - m*g*l*sin(Y(0)) - b*Y(1)) / I;
	return dY;
}



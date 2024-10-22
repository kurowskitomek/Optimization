#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try {
			double *p = new double[2]{0, 0};
			// Tu wpisz kod funkcji
			int i = 0;
			double x1 = x0 + d;
			if (ff(x1, 0, 0) == ff(x0, 0, 0)) {
					p[0] = x0;
					p[1] = x1;
					return p;
			}
			if (ff(x1, 0, 0) > ff(x0, 0, 0)) {
					d = -d;
					x1 = x0 + d;
					if (ff(x1, 0, 0) >= ff(x0, 0, 0)) {
							p[0] = x1;
							p[1] = x0 - d;
							return p;
					}
			}
			double last = 0;
			do {
					if (i > Nmax)
							throw;
					i++;
					last = x1;
					x1 = x0 + pow(alpha, i) * d;
			} while (ff(last, 0, 0) <= ff(x1, 0, 0));
			if (d > 0) {
					p[0] = x0 + pow(alpha, i - 2) * d;
					p[1] = x1;
					return p;
			}
			p[1] = x0 + pow(alpha, i - 2) * d;
			p[0] = x1;

			return p;
	} catch (string ex_info) {
			throw("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		int i = 0;
		double c = (a + b) / 2.0;

		solution A(a);
		solution B(b);
		solution C(c);
		solution D(0);
		solution Dold(0);

		A.fit_fun(ff, ud1, ud2);
		B.fit_fun(ff, ud1, ud2);
		C.fit_fun(ff, ud1, ud2);

		while (true)
		{
			Dold.x = D.x;
			Dold.fit_fun(ff, ud1, ud2);

			matrix l = 
				A.y * (pow(B.x, 2) - pow(C.x, 2)) +
				B.y * (pow(C.x, 2) - pow(A.x, 2)) +
				C.y * (pow(A.x, 2) - pow(B.x, 2));

			matrix m = 
				A.y * (B.x - C.x) +
				B.y * (C.x - A.x) +
				C.y * (A.x - B.x);

			if (m <= 0)
				throw;

			D.x = 0.5 * l / m;
			D.fit_fun(ff, ud1, ud2);

			double test_a = A.x(0);
			double test_b = B.x(0);
			double test_c = C.x(0);
			double test_d = D.x(0);


			double test_ay = A.y(0);
			double test_by = B.y(0);
			double test_cy = C.y(0);
			double test_dy = D.y(0);

			if (A.x < D.x && D.x < C.x)
				if (D.y < C.y)
				{
					//a_tab[2] = a_tab[1];
					B.x = C.x;
					B.fit_fun(ff, ud1, ud2);

					C.x = D.x;
					C.fit_fun(ff, ud1, ud2);
				}
				else
				{
					A.x = D.x;
					A.fit_fun(ff, ud1, ud2);
				}
			else
				if (C.x < D.x && D.x < B.x)
				{
					if (ff(D.x, NAN, NAN) < ff(C.x, NAN, NAN))
					{
						A.x = C.x;
						A.fit_fun(ff, ud1, ud2);

						C.x = D.x;
						C.fit_fun(ff, ud1, ud2);
					}
					else
					{
						B.x = D.x;
						B.fit_fun(ff, ud1, ud2);
					}
				}
				else
					throw;

			i++;

			if (solution::f_calls > Nmax)
				throw;

			if (B.x - A.x < epsilon || abs(D.x(0) - Dold.x(0)) < gamma)
				break;
		}

		Xopt = D;

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}

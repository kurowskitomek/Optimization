#include"opt_alg.h"
#include <iostream>





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

double* expansion(matrix(ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{

		double* p = new double[2] { 0, 0 }; //przedzial 
		int i = 0; //i = 0 
		solution X0(x0);
		solution X1(x0 + d);

		if (X1.fit_fun(ff, ud1, ud2) == X0.fit_fun(ff, ud1, ud2))
		{
			p[0] = x0;
			p[1] = x0 + d;
			return p;
		}

		if (X1.fit_fun(ff, ud1, ud2) > X0.fit_fun(ff, ud1, ud2))
		{
			d = -d;
			X1.x = X0.x + d;
			if (X1.fit_fun(ff, ud1, ud2) >= X0.fit_fun(ff, ud1, ud2))
			{
				p[0] = m2d(X1.x);
				p[1] = m2d(X0.x - d);
				return p;
			}
		}

		solution XH(x0);
		solution XI(x0);
		solution XJ(x0);

		do
		{
			if (solution::f_calls > Nmax)
				throw ("Max iter!\n");			

			i++;
			XH = XI;
			XI = XJ;
			XJ.x = X0.x + pow(alpha, i) * d;
		}
		while (XI.fit_fun(ff, ud1, ud2) > XJ.fit_fun(ff, ud1, ud2));

		if (d > 0)
		{
			p[0] = m2d(XH.x);
			p[1] = m2d(XJ.x);
			return p;
		}

		p[0] = m2d(XJ.x);
		p[1] = m2d(XH.x);
		
		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}
//
//
//double* expansion2(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
//{
//	try
//	{
//		double* p = new double[2] { 0, 0 };
//		// Tu wpisz kod funkcji
//		int i = 0;
//		double x1 = x0 + d;
//
//		if (ff(x1, NAN, NAN) == ff(x0, NAN, NAN))
//		{
//			p[0] = x0;
//			p[1] = x1;
//			return p;
//		}
//
//		if (ff(x1, NAN, NAN) > ff(x0, NAN, NAN))
//		{
//			d = -d;
//			x1 = x0 + d;
//
//			if (ff(x1, NAN, NAN) >= ff(x0, NAN, NAN))
//			{
//				p[0] = x1;
//				p[1] = x0 - d;
//				return p;
//			}
//		}
//
//		double last = 0;
//
//		while (true)
//		{
//			if (solution::f_calls > Nmax)
//				throw;
//
//			i++;
//			last = x1;
//			x1 = x0 + pow(alpha, i) * d;
//
//			if (ff(last, NAN, NAN) <= ff(x1, NAN, NAN))
//				break;
//		}
//
//		if (d > 0)
//		{
//			p[0] = last;//x0 + pow(alpha, i - 2) * d;
//			p[1] = x1;
//			return p;
//		}
//
//		p[0] = x1;//x0 + pow(alpha, i - 2) * d;
//		p[1] = last;
//
//		return p;
//	}
//	catch (string ex_info)
//	{
//		throw("double* expansion(...):\n" + ex_info);
//	}
//}

//returns index of fib number greater than x
double fibi(double x, double _prev = 0, double _curr = 0, double _index = 1)
{
    if(_curr > x)
    {
      return _index;
    }
    else {
    if(_curr != 0)
        return fibi(x, _curr, _prev+_curr, _index+1);
    else {
      return fibi(x, 0, 1, _index+1);
    }
    }
}
//returns fib number with index of index
double fibo(double index, double _prev = 0, double _curr = 0)
{
  if(index  == 1)
  {
    return _curr;
  }
  else {
    if(_curr == 0)
      return fibo(index - 1, 0, 1);
    else 
      return fibo(index - 1, _curr, _prev + _curr);
  }
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji
		double k = fibi((b - a) / epsilon);

		int intK = k;

		double A, B, C, D;

		A = a;
		B = b;
		C = b - fibo(k - 1) / fibo(k) * (b - a);
		D = a + b - C;

		solution A_sol, B_sol, C_sol, D_sol;

		for (int i = 0; i < k - 3; i++)
		{
			A_sol.x = A;
			C_sol.x = C;
			D_sol.x = D;
			C_sol.fit_fun(ff);
			D_sol.fit_fun(ff);

      cout<<B-A<<endl;
			if (C_sol.y < D_sol.y)
			{
				A = A;
				B = D;
			}
			else
			{
				B = B;
				A = C;
			}
			C = B - fibo(k - i - 2) / fibo(k - i - 1) * (B - A);
			D = A + B - C;
		}

		Xopt = C;
		Xopt.fit_fun(ff, ud1, ud2);

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
      cout<<B.x-A.x<<endl;
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
				return NULL;

			D.x = 0.5 * l / m;
			D.fit_fun(ff, ud1, ud2);

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
					if (D.y < C.y)
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
					return NULL;

			i++;

			if (solution::f_calls > Nmax)
				return NULL;

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

		solution x = x0;
		solution xb = x;
		solution xb2 = x;
		while (s > epsilon)
		{
			xb.fit_fun(ff, ud1, ud2);
			x = HJ_trial(ff, xb, s, ud1, ud2);
			x.fit_fun(ff, ud1, ud2);

			if (x.y < xb.y)
			{
				while (x.y < xb.y)
				{
					xb2.x = xb.x;
					xb.x = x.x;
					x.x = 2.0 * xb.x - xb2.x;
					x = HJ_trial(ff, x, s, ud1, ud2);
					xb.fit_fun(ff, ud1, ud2);
					x.fit_fun(ff, ud1, ud2);

					if (solution::f_calls > Nmax)
					{
						Xopt.flag = 0;
						return Xopt;
					}
				}
			}
			else
			{
				s = alpha * s;
			}
		}
		Xopt = xb;


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
		solution xps;
		solution xms;
		for (int i = 0; i < get_len(XB.x); i++)
		{
			xps.x = XB.x + s * exp(i + 1);
			xms.x = XB.x - s * exp(i + 1);
			xps.fit_fun(ff, ud1, ud2);
			xms.fit_fun(ff, ud1, ud2);
			if (xps.y < XB.y)
			{
				XB = xps;
			}
			else if (xms.y < XB.y)
			{
				XB = xms;
			}
		}

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	const int n = 2; // Amount of dimensions
	try
	{
		matrix D(n, n), l(n, 1), p(n, 1), s(s0);
		for (int i = 0; i < n; i++)
			D(i, i) = 1;

		solution X, Xt;
		X.x = x0;
		X.fit_fun(ff);
		while (true)
		{
			for (int j = 0; j < n; j++)
			{
				Xt.x = X.x + s(j) * D[j];
				Xt.fit_fun(ff);
				if (Xt.y < X.y)
				{
					X = Xt;
					l(j) += s(j);
					s(j) *= alpha;
				}
				else
				{
					p(j) = p(j) + 1;
					s(j) *= -beta;
				}
			}

			bool basisChange = true;
			for (int j = 0; j < n; j++)
				if (p(j) == 0 || l(j) == 0)
				{
					basisChange = false;
					break;
				}

			if (basisChange)
			{
				matrix Q(n, n), v(n, 1);
				for (int i = 0; i < n; i++)
					for (int j = 0; j <= i; j++)
						Q(i, j) = l(i);
				Q = D * Q;
				v = Q[0] / norm(Q[0]);
				D.set_col(v, 0);
				for (int i = 1; i < n; i++)
				{
					matrix temp(n, 1);
					for (int j = 0; j < i; j++)
						temp = temp + (trans(Q[i]) * D[j]) * D[j];
					v = Q[i] - temp;
					D.set_col(v, i);
				}
				s = s0;
				l = matrix(n, 1);
				p = matrix(n, 1);
			}
			double max = abs(s(0));
			for (int i = 1; i < n; i++)
				if (max < abs(s(i)))
					max = abs(s(i));
			if (max < epsilon || solution::f_calls > Nmax)
				return X;
		}
	}
	catch (string ex_info)
	{
		throw("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution XB;
		XB.x = x0;
		XB.fit_fun(ff, ud1, ud2);

		solution XT;
		XT = XB;

		double s = 0.5; //D³ugoœæ boku trójk¹ta
		double alpha = 1.0; //Wspó³czynnik odbicia
		double beta = 0.5; //Wspó³czynnik zwê¿enia
		double gamma = 2.0; //Wspó³czynnik ekspansji
		double delta = 0.5; //Wspó³czynnik redukcji

		do
		{
			XT.x = XB.x;
			XT = sym_NM(ff, XB.x, s, alpha, beta, gamma, delta, epsilon, Nmax, ud1, c);
			c *= dc;

			if (solution::f_calls > Nmax)
			{
				XT.flag = 0;
				throw std::string("Maximum amount of f_calls reached!");
			}

			if (norm(XT.x - XB.x) < epsilon)
				break;

			XB = XT;
		}
		while (true);

		return XT;
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
		//Funkcja pomocnicza do znajdywania maksymum normy
		auto max = [&](std::vector<solution> sim, int i_min) -> double
			{
				double result = 0.0;
				for (int i = 0; i < sim.size(); ++i)
				{
					double normal = norm(sim[i_min].x - sim[i].x);
					if (result < normal)
						result = normal;
				}
				return result;
			};

		int n = get_len(x0);

		//Tworzenie bazy ortogonalnej
		matrix d = matrix(n, n);
		for (int i = 0; i < n; ++i)
			d(i, i) = 1.0;

		//Tworzenie simplexu i uzupe³nianie go danymi
		std::vector<solution> simplex;
		simplex.resize(n + 1);
		simplex[0].x = x0;
		simplex[0].fit_fun(ff, ud1, ud2);
		for (int i = 1; i < simplex.size(); ++i)
		{
			simplex[i].x = simplex[0].x + s * d[i - 1];
			simplex[i].fit_fun(ff, ud1, ud2);
		}

		//Indeks najmniejszej wartoœci wierzcho³ka simplexu
		int i_min{};
		//Indeks najwiêkszej wartoœci wierzcho³ka simplexu
		int i_max{};

		while (max(simplex, i_min) >= epsilon)
		{
			//Wyznaczanie maksymalnego i minimalnego indeksu
			i_min = 0;
			i_max = 0;
			for (int i = 1; i < simplex.size(); ++i)
			{
				if (simplex[i].y < simplex[i_min].y)
					i_min = i;
				if (simplex[i].y > simplex[i_max].y)
					i_max = i;
			}

			//Wyznaczenie œrodka ciê¿koœci
			matrix simplex_CoG{};
			for (int i = 0; i < simplex.size(); ++i)
			{
				if (i == i_max)
					continue;
				simplex_CoG = simplex_CoG + simplex[i].x;
			}
			simplex_CoG = simplex_CoG / simplex.size();

			//Obliczanie wartoœci funkcji odbitego simplexu
			solution simplex_reflected{};
			simplex_reflected.x = simplex_CoG + alpha * (simplex_CoG - simplex[i_max].x);
			simplex_reflected.fit_fun(ff, ud1, ud2);

			if (simplex_reflected.y < simplex[i_min].y)
			{
				//Obliczanie wartoœci funkcji powiêkszonego simplexu
				solution simplex_expansion{};
				simplex_expansion.x = simplex_CoG + gamma * (simplex_reflected.x - simplex_CoG);
				simplex_expansion.fit_fun(ff, ud1, ud2);
				if (simplex_expansion.y < simplex_reflected.y)
					simplex[i_max] = simplex_expansion;
				else
					simplex[i_max] = simplex_reflected;
			}
			else
			{
				if (simplex[i_min].y <= simplex_reflected.y && simplex_reflected.y < simplex[i_max].y)
					simplex[i_max] = simplex_reflected;
				else
				{
					//Obliczanie wartoœci funkcji pomniejszonego simplexu
					solution simplex_narrowed{};
					simplex_narrowed.x = simplex_CoG + beta * (simplex[i_max].x - simplex_CoG);
					simplex_narrowed.fit_fun(ff, ud1, ud2);
					if (simplex_narrowed.y >= simplex[i_max].y)
					{
						for (int i = 0; i < simplex.size(); ++i)
						{
							if (i == i_min)
								continue;
							simplex[i].x = delta * (simplex[i].x + simplex[i_min].x);
							simplex[i].fit_fun(ff, ud1, ud2);
						}
					}
					else
						simplex[i_max] = simplex_narrowed;
				}
			}

			if (solution::f_calls > Nmax)
			{
				simplex[i_min].flag = 0;
				throw std::string("Maximum amount of f_calls reached!");
			}

		}

		return simplex[i_min];
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

double dotProduct(const matrix& v1, const matrix& v2)
{
	// Ensure both vectors have the same dimension
	int* size1 = get_size(v1);
	int* size2 = get_size(v2);

	if (size1[0] != size2[0] || size1[1] != 1 || size2[1] != 1)
	{
		throw std::string("dotProduct: Both inputs must be column vectors of the same size.");
	}

	// Calculate the dot product
	double result = 0.0;
	for (int i = 0; i < size1[0]; ++i)
	{
		result += v1(i, 0) * v2(i, 0); // Sum the element-wise products
	}

	delete[] size1;
	delete[] size2;
	return result;
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		// Initialize the solution object with the initial guess x0
		solution Xopt(x0);

		// Compute the initial gradient (g = A*x - b)
		matrix g = gf(Xopt.x, ud1, ud2);  // Gradient at x0

		// Store the gradient in the solution object
		Xopt.grad(gf, ud1, ud2);

		// Calculate the initial residual r = -g
		matrix r = -g;

		// The initial search direction p is the same as the residual
		matrix p = r;

		// Compute the squared residual (rsOld = r^T * r)
		double rsOld = dotProduct(r, r);

		// Iterate for up to Nmax iterations or until convergence
		for (int k = 0; k < Nmax; ++k)
		{
			// Compute A * p (where ff is the matrix-vector multiplication function)
			matrix Ap = ff(p, ud1, ud2);

			// Calculate the step size alpha
			double alpha = rsOld / dotProduct(p, Ap);  // Step size for the current iteration

			// Update the solution x
			Xopt.x = Xopt.x + (p * alpha);  // x = x + alpha * p

			// Update the gradient g
			g = g + (Ap * alpha);  // g = g + alpha * A * p (gradient update)

			// Update the residual r
			r = r - (Ap * alpha);  // r = r - alpha * A * p (residual update)

			// Compute the new squared residual (rsNew = r^T * r)
			double rsNew = dotProduct(r, r);

			// Check for convergence
			if (sqrt(rsNew) < epsilon)  // If the residual is small enough, we stop
			{
				std::cout << "Converged in " << k + 1 << " iterations." << std::endl;
				break;
			}

			// Update the search direction p
			p = r + (rsNew / rsOld) * p;  // p = r + (rsNew / rsOld) * p

			// Update the old residual squared value for the next iteration
			rsOld = rsNew;
		}

		// Return the solution object containing the final result for x
		return Xopt;
	}
	catch (const std::string& ex_info)
	{
		// Exception handling
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

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
		// Initialize variables
		int n = get_dim(x0);           // Number of dimensions
		solution Xopt;                 // To store the optimal solution
		int i = 0;                     // Iteration counter
		matrix xb = x0;                // Current point
		matrix sj = s0;                // Current step size vector
		matrix lambdaj(n);             // Lambda vector for direction length
		matrix pj(n);                  // Count of steps taken in each direction
		matrix dj(n);                  // Direction vector for current iteration

		// Initialize direction vectors
		for (int j = 0; j < n; j++)
		{
			dj(j) = 1; // Initialize with unit vectors
			lambdaj(j) = 0; // Initialize λj(0) = 0
			pj(j) = 0; // Initialize pj(0) = 0
		}

		// Start of the loop
		while (true)
		{
			// Loop over each dimension
			for (int j = 0; j < n; j++)
			{
				// Check if the function value decreases
				if (ff(xb + sj * dj, ud1, ud2) < ff(xb, ud1, ud2))
				{
					xb = xb + sj * dj;  // Update xb
					lambdaj(j) += sj;    // Update λj
					sj = alpha * sj;      // Expand step size
				}
				else
				{
					sj = -beta * sj;     // Contract step size
					pj(j) += 1;          // Increment the step count
				}
			}

			// Update iteration counter
			i++;
			// Store the current position
			Xopt = solution(xb);
			Xopt.y = ff(xb, ud1, ud2); // Store the function value at the current point


			bool basis_change_needed = false;
			for (int j = 0; j < n; j++)
			{
				if (lambdaj[j] != 0 && pj[j] != 0)
				{
					basis_change_needed = true;
					break;
				}
			}

			// Check for direction basis change
			if (basis_change_needed) // Check if there were steps taken
			{
				// Change basis for directions `dj` based on your method's logic
				//for (int j = 0; j < n; j++)
				//{
				//	dj(j) = /* new direction logic here based on your requirements */;
				//}

				// Reset lambda and step count
				lambdaj = 0; // λj(i) = 0
				pj = 0;      // pj(i) = 0
				sj = s0;     // Reset sj to initial step size

				for (int j = 0; j < n; j++)
				{
					if (pj[j] > some_threshold)
					{ // Adjust this condition based on your needs
						dj[j] = (dj[j] == 1) ? -1 : 1; // Simple flip, can also be reinitialized
						// Alternatively, you might want to calculate a new direction based on your specific requirements.
						// Example: dj[j] = calculate_new_direction(xb, j); // Custom function to calculate new direction.
					}
					else
					{
						dj[j] = 1; // Reset to a unit direction if it hasn't failed
					}
				}
			}

			// Check if maximum number of function calls exceeded
			if (solution::f_calls > Nmax)
			{
				throw std::runtime_error("Exceeded maximum function calls");
			}

			// Check termination condition
			double maxStepSize = sj.max_abs(); // Assuming max_abs gets the maximum absolute value of sj
			if (maxStepSize < epsilon) break; // Terminate if max step size is less than ε
		}

		return Xopt; // Return the optimal solution
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

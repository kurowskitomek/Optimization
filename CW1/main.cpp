/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include"opt_alg.h"

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		//lab0();
		lab1();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	system("pause");
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;
	int Nmax = 10000;
	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
	solution opt;
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Wahadlo
	Nmax = 1000;
	epsilon = 1e-2;
	lb = 0;
	ub = 5;
	double teta_opt = 1;
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
	ofstream Sout("symulacja_lab0.csv");
	Sout << hcat(Y[0], Y[1]);
	Sout.close();
	Y[0].~matrix();
	Y[1].~matrix();
}

void real_problem_lab1()
{
	solution::clear_calls();
	int Nmax = 1000;
	double epsilon = 1e-8;
	double gamma = 1e-10;
	solution optLag, optFib;

	cout << "Optymalizacja Lagrang'e\n";

	optLag = lag(ff1R, 0.0001, 0.01, epsilon, gamma, Nmax, NAN, NAN);
	cout << optLag << endl << endl;
	solution::clear_calls();

	cout << "\nOptymalizacja Fibonacci'e\n";

	optFib = fib(ff1R, 0.0001, 0.01, epsilon, NAN, NAN);
	cout << optFib << endl << endl;
	solution::clear_calls();

	std::ofstream symLag("symulacja_lab1_LAG.csv");
	double Y_m[3] = { 5, 1, 10 };

	if (symLag.is_open())
	{
		symLag << "t,Volume_A,Volume_B,Temperature_B\n";
		matrix Y0 = matrix(3, Y_m);
		matrix* Y = solve_ode(df1, 0, 1, 2000, Y0, NULL, optLag.x);
		for (int i = 0; i < get_len(Y[0]); i++)
		{
			symLag << Y[0](i) << "," << Y[1](i, 0) << "," << Y[1](i, 1) << "," << Y[1](i, 2) << "\n";
		}
		symLag.close();
		std::cout << "Results saved to symulacja_lab1_LAG.csv" << std::endl;
	}
	else
	{
		std::cerr << "Failed to open symulacja_lab1_LAG.csv for writing." << std::endl;
	}

	std::ofstream symFib("symulacja_lab1_FIB.csv");

	if (symFib.is_open())
	{
		symFib << "t,Volume_A,Volume_B,Temperature_B\n";
		matrix Y0 = matrix(3, Y_m);
		matrix* Y = solve_ode(df1, 0, 1, 2000, Y0, NULL, optFib.x);
		for (int i = 0; i < get_len(Y[0]); i++)
		{
			symFib << Y[0](i) << "," << Y[1](i, 0) << "," << Y[1](i, 1) << "," << Y[1](i, 2) << "\n";
		}
		symFib.close();
		std::cout << "Results saved to symulacja_lab1_FIB.csv" << std::endl;
	}
	else
	{
		std::cerr << "Failed to open symulacja_lab1_FIB.csv for writing." << std::endl;
	}
}

void test_problem_lab1()
{
	solution y;
	solution::clear_calls();

	cout << "EXP:\n";

	double* y_ex = expansion(lab1f, 30, 1.0, 2.25, 10000, NAN, NAN);
	cout << y_ex[0] << ", " << y_ex[1] << "\n";

	cout << "\nFIB:\n";

	solution::clear_calls();
	y = fib(lab1f, 50, 100, 0.1, NAN, NAN);

	cout << y << endl << endl;

	cout << "\nLAG:\n";

	solution::clear_calls();
	y = lag(lab1f, 0, 100, 0.1, 0.00001, 10000, NAN, NAN);

	cout << y << endl << endl;
}

void lab1()
{
	test_problem_lab1();
	real_problem_lab1();
}

void lab2()
{

}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}

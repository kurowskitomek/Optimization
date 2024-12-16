/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/
#include <iostream>
#include <random>
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
		//lab1();
		lab3();
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

/*void test_problem_lab1()
{
	for(int i =0; i < 100; i++){
		std::random_device rd;
    		std::mt19937 gen(rd());
		solution y;
		solution::clear_calls();

		cout << "EXP:\n";

		std::uniform_real_distribution<double> dist(1.1, 2.5);
		double random = dist(gen);
	
		double* y_ex = expansion(lab1f, random, 1.0, 2.25, 10000, NAN, NAN);
		cout << y_ex[0] << ", " << y_ex[1] << "\n";

		cout << "\nFIB:\n";

		solution::clear_calls();
		y = fib(lab1f, y_ex[0], y_ex[1], 0.1, NAN, NAN);

		cout << y << endl << endl;

		cout << "\nLAG:\n";

		solution::clear_calls();
		y = lag(lab1f, y_ex[0], y_ex[1], 0.1, 0.00001, 10000, NAN, NAN);

		cout << y << endl << endl;
	}
}*/

void test_problem_lab1()
{
    // Open a CSV file for writing
    std::ofstream csv_file("results.csv");
    if (!csv_file.is_open()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }

    // Write headers to the CSV file
    //csv_file << "expansion_start ; expansion_end ; fib_x ; fib_y ; lag_x ; lag_y\n";

    std::uniform_real_distribution<double> dist(-100.0, 100.0);
    double x0 = 0.25;
    for(int j = 0; j< 3;j++){
    x0+=1.0;
    for(int i = 0; i < 100; i++){
        std::random_device rd;
        std::mt19937 gen(rd());
        solution y, x;
        solution::clear_calls();

        // Expansion phase
        double random = dist(gen);
        double* y_ex = expansion(ff1T, random, 1.0, x0, 10000, NAN, NAN);
        double expansion_start = y_ex[0];
        double expansion_end = y_ex[1];
        //std::cout << "EXP:\n" << expansion_start << ", " << expansion_end << "\n";

        csv_file <<x0<<";"<<random<<";"<<expansion_start<<";"<<expansion_end<<";"<< solution::f_calls<<";";
        // Fibonacci phase
        solution::clear_calls();
        x = fib(ff1T, expansion_start, expansion_end, 0.1, NAN, NAN);
        //std::cout << "FIB:\n" << x<< std::endl;
        std::cout << "-------------------------------------------" << std::endl;
        std::string loc_glb;
        if(x.x > 50){
        loc_glb = "glob";
        }
        else {
          loc_glb = "loc";
        }

        csv_file <<x.x<<";"<<x.y<<";"<< solution::f_calls<<";"<<loc_glb<<";";

        // Lagrange phase
        solution::clear_calls();
        y = lag(ff1T, expansion_start, expansion_end, 0.1, 0.00001, 10000, NAN, NAN);
        //std::cout << "LAG:\n" << y<< std::endl << std::endl;
        std::cout << "-------------------------------------------next" << std::endl;

        // Write the results to the CSV file
        csv_file <<y.x<<";"<<y.y<<";"<< solution::f_calls<<";"<<loc_glb<<"\n";
    }
    }

    // Close the CSV file
    csv_file.close();
}


void lab1()
{
	test_problem_lab1();
	//real_problem_lab1();
}

void lab2()
{
	solution y;
	solution::clear_calls();

	cout << "HJ:\n";
	double s = 0.001;
	matrix s0(s);
	matrix x0(0);
	s0.add_row(s);
	x0.add_row(0);

	matrix s1(s);
	matrix x1(0);
	s1.add_row(s);
	x1.add_row(0);
	//s0.add_col(s);

	solution y_ex = HJ(ff2T, x1, 0.001, 0.2, 0.001, 1000, NAN);
	cout << y_ex.x << "\n";

	cout << "\nRosen:\n";

	solution::clear_calls();

	y = Rosen(ff2T, x0, s0, 2, 0.0001, 0.001, 100000000);

	cout << y.x << endl << endl;
}

void lab4()
{
	solution y;
	solution::clear_calls();

	cout << "CG:\n";
	double s = 0.001;
	matrix s0(s);
	matrix x0(0);
	s0.add_row(s);
	x0.add_row(0);

	matrix s1(s);
	matrix x1(0);
	s1.add_row(s);
	x1.add_row(0);
	//s0.add_col(s);

	solution y_ex = CG(ff4T, gradient_example, x1, 0.001, 0.2, 0.001, 1000, NAN);
	cout << y_ex.x << "\n";

}

void lab3()
{
	double epsilon = 1E-3;
	int Nmax = 10000;
	double c_inside = 100;
	double dc_inside = 0.2;
	double c_outside = 1.0;
	double dc_outside = 1.5;

	std::ofstream sout_test("test_lab3.csv");

	sout_test << "x0_1;x0_2;x1_out;x2_out;norm_out;y_out;f_calls_out;x1_in;x2_in;norm_in;y_in;f_calls_in\n";

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> x0_dist(1.5, 5.5);
	std::stringstream test_ss;
	solution test_sol;
	matrix a = matrix(4.0);
	matrix test_x0{};

	for (int i = 0; i < 3; ++i)
	{
		if (i == 0)
			a = matrix(4.0);
		else if (i == 1)
			a = matrix(4.4934);
		else
			a = matrix(5.0);

		for (int j = 0; j < 100; ++j)
		{
			test_x0 = matrix(2, new double[2] { x0_dist(gen), x0_dist(gen) });
			test_ss << test_x0(0) << "; " << test_x0(1) << "; ";

			// Zewnêtrzne rozwi¹zanie
			test_sol = pen(ff3T_outside, test_x0, c_outside, dc_outside, epsilon, Nmax, a);
			//cout << test_sol;
			test_ss << test_sol.x(0) << "; " << test_sol.x(1) << "; "
				<< sqrt(pow(test_sol.x(0), 2) + pow(test_sol.x(1), 2)) << "; "
				<< test_sol.y << "; " << test_sol.f_calls << "; ";
			solution::clear_calls();

			// Wewnêtrzne rozwi¹zanie
			test_sol = pen(ff3T_inside, test_x0, c_inside, dc_inside, epsilon, Nmax, a);
			//cout << test_sol;
			test_ss << test_sol.x(0) << "; " << test_sol.x(1) << "; "
				<< sqrt(pow(test_sol.x(0), 2) + pow(test_sol.x(1), 2)) << "; "
				<< test_sol.y << "; " << test_sol.f_calls << "\n";
			solution::clear_calls();
		}
	}

	sout_test << test_ss.str();
	sout_test.close();

	std::ofstream Sout("sym_lab3.csv");

	Sout << "time;x_position;y_position\n";

	matrix ud1 = matrix(5, new double[5]
	{
			0.47, //Wspó³czynnik oporu (C) [-]
			1.2, //Gêstoœæ powietrza (rho) [kg/m^3]
			0.12, //Promieñ pi³ki (r) [m]
			0.6, //Masa pi³ki (m) [kg]
			9.81 //Przyœpieszenie ziemskie (g) [m/s^2]
	});

	//Pocz¹tkowe wartoœci szukania minimum
	matrix x0 = matrix(2, new double[2] { -4.0, 2.0 });

	//Szukanie optymalnej prêdkoœci pocz¹tkowej po osi x i pocz¹tkowej prêdkoœci obrotowej
	solution opt = pen(ff3R, x0, c_outside, dc_outside, epsilon, Nmax, ud1);
	std::cout << opt << "\n";

	std::cout << opt.x(0) << std::endl;
	std::cout << opt.x(1) << std::endl;

	//Symulacja lotu pi³ki dla wyznaczonych ograniczeñ
	matrix Y0(4, new double[4] { 0.0, opt.x(0), 100, 0 });
	matrix* Y = solve_ode(df3, 0.0, 0.01, 7.0, Y0, ud1, opt.x(1));

	// Zapis wyników symulacji do pliku CSV
	int num_rows = get_len(Y[0]);
	for (int i = 0; i < num_rows; ++i)
	{
		Sout << Y[0](i) << "; "  // Czas
			<< Y[1](i, 0) << "; "  // Pozycja w osi x
			<< Y[1](i, 2) << "\n";  // Pozycja w osi y
	}

	// Zamkniêcie pliku i usuniêcie dynamicznie alokowanej pamiêci
	Sout.close();
	delete[] Y;
}

void lab5()
{

}

void lab6()
{

}

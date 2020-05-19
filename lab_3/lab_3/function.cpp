#pragma once
#include "vector_new.h"


double	f(double x, Data my_data)
{
//	return sin(x * PI);//test1
//	return x * (1 - x);//test2

//	return x * cos(PI * x);	//variant2
	return (x + 0.2) * sin(PI * x / 2);//variant6
}
double	g(double x, Data my_data)
{
//	return 0.0;//test1 test2

//	return x * (2 - x);//variant2
	return pow(1 + x, 2);//variant6
}
double	phi(double t, Data my_data)
{
//	return 0.0;//test1 test2

//	return 2 * t;//variant2
	return 0.0;//variant6
}
double	psi(double t, Data my_data)
{
//	return 0.0;//test1 test2

//	return -1;//variant2
	return 1.2 * (t + 1);//variant6
}

double	diff_f(double x, std::function<double(double, Data)> &f, Data my_data)
{
	double	diff;
	double	eps = 1.e-3;

	diff = (f(x + eps, my_data) - 2 * f(x, my_data) + f(x - eps, my_data)) / pow(eps, 2);

	return diff;
}

double integrate(std::vector<double> y1, std::vector<double> y2, double h, double tau, Data my_data)
{
	double sum = 0.0;
	int n = y1.size();

	sum += 0.5 * (pow((y2[1] - y1[1]) / tau, 2) + pow(my_data.a, 4) * pow((y1[2] - y1[0]) / (h * 2), 2)) * h;
	for (int i = 1; i < n - 3; i++)
		sum += (pow((y2[i + 1] - y1[i + 1]) / tau, 2) + pow(my_data.a, 4) * pow((y1[i + 2] - y1[i]) / (h * 2), 2)) * h;
	sum += 0.5 * (pow((y2[n - 2] - y1[n - 2]) / tau, 2) + pow(my_data.a, 4) * pow((y1[n - 1] - y1[n - 3]) / (h * 2), 2)) * h;

	return sum;
}

void	scheme_cross(double h, double tau, Data my_data, std::string order)
{
	int					n = my_data.L / h;
	std::ofstream		fout;
	std::ofstream		fout_enrgy;
	std::vector<double> y1(n + 1);
	std::vector<double> y2(n + 1);
	std::vector<double> y3(n + 1);
	std::string			str = "Scheme_cross";

	str += order + ".txt";
	fout.open(str);
	fout_enrgy.open("Energy.txt");
	// инициализация начальными данными
	for (int i = 0; i <= n; i++)
		y1[i] = my_data.initial_condition(i*h, my_data);
	for (int i = 0; i != y1.size(); i++)
		fout << 0 << '\t' << i * h << '\t' << y1[i] << '\n';


	//вычисление второго слоя
	for (int i = 0; i <= n; i++)
		y2[i] = y1[i] + tau * my_data.diff_initial_cond(i * h, my_data) + pow(my_data.a * tau, 2) / 2 * diff_f(i * h, my_data.initial_condition, my_data);
	for (int i = 0; i != y2.size(); i++)
		fout << tau << '\t' << i * h << '\t' << y2[i] << '\n';
//	fout_enrgy << tau << '\t' << integrate(y1, y2, h, tau, my_data) << '\n';


	//вычисление по временным слоям схемой "крест"
	for (double j = 2 * tau; j <= my_data.T; j += tau)
	{
		y3[0] = my_data.left_boarder(j, my_data);

		for (int i = 1; i < n; i++)
			y3[i] = 2 * y2[i] - y1[i] + pow(my_data.a * tau / h, 2) * (y2[i + 1] - 2 * y2[i] + y2[i - 1]);

		y3[n] = my_data.right_boarder(j, my_data);

		for (int i = 0; i != y3.size(); i++)
			fout << j << '\t' << i * h << '\t' << y3[i] << '\n';
//		fout_enrgy << j << '\t' << integrate(y2, y3, h, tau, my_data) << '\n';
		y1 = y2;
		y2 = y3;
	}
	fout.close();
	fout_enrgy.close();
	std::cout << "Your file is ready" << std::endl;
}

#pragma once
#ifndef VECTOR_NEW_H
#define VECTOR_NEW_H

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <functional>
#include <string>

extern double	EPS;
extern double	PI;

struct Data
{
	double	L = 2;//test1
//	double	L = 1;//test2
	double	a = 1;
	double	T = 1;

	std::function<double(double, Data)>	left_boarder;
	std::function<double(double, Data)>	right_boarder;
	std::function<double(double, Data)>	initial_condition;
	std::function<double(double, Data)>	diff_initial_cond;
};

void print(std::vector<double> x);

// переодпределение операций под вектора
std::vector<double> operator * (double a, std::vector<double> b);
std::vector<double> operator + (std::vector<double> a, std::vector<double> b);
std::vector<double> operator - (std::vector<double> a, std::vector<double> b);
std::vector<double> operator / (std::vector<double> a, double b);



double	f(double x, Data my_data);
double	g(double x, Data my_data);
double	phi(double t, Data my_data);
double	psi(double t, Data my_data);

void	scheme_cross(double h, double tau, Data my_data, std::string order);


void	error_check(double h, double tau, Data DATA);


#endif //VECTOR_NEW_H
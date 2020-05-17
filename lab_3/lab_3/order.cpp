#pragma once
#include "vector_new.h"

void error_check(double h, double tau, Data DATA)
{
	double	q = 2;

	if (true)//(fabs(sigma - 0.5) >= EPS)
	{
		for (int i = 1; i <= 5; i++)
		{
			scheme_cross(h, tau, DATA, "_" + std::to_string(pow(q, i)));
			h /= q;
			tau /= q;
		}
	}
	std::cout << "errors is ready" << std::endl;
}
 
#include "vector_new.h"

double	EPS = 1.e-3;
double	PI = 3.1415926535;

int main()
{
	Data	DATA;
	double	h = 0.01;
	double	tau = 0.01;

	DATA.left_boarder = &phi;
	DATA.right_boarder = &psi;
	DATA.initial_condition = &f; 
	DATA.diff_initial_cond = &g;

	scheme_cross(h, tau, DATA, "");
// 	error_check(0.1, 0.1, DATA);


	system("pause");
	return (0);
}
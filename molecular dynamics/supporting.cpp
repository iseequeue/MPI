#include "supporting.h"

double pbc(const double& x, const double& dx, const double& L_calc_cell)
{
	if (x + dx >= L_calc_cell) return x + dx - L_calc_cell;
	if (x + dx < -0.0)		   return x + dx + L_calc_cell;
	return x + dx;
}

double length(const std::vector<double>& r)
{
	double d = 0.0;  //counter
	for (auto i = 0u; i < r.size(); i++)
		d += r[i]*r[i];
	//return std::pow(d, 0.5);
	return sqrt(d);
}

double scalar(const std::vector<double>& r1, const std::vector<double>& r2)
{
	double w = 0.0; //counter
	for (auto i = 0u; i < r1.size(); i++)
		w += r1[i] * r2[i];
	return w;
}

void print(const std::vector<std::vector<double> >& r)
{
	for (auto i = 0u; i < r.size(); i++)
	{
		std::cout << i << ' ' << r[i][0] << ' ' << r[i][1] << ' ' << r[i][2] << '\n';
	}
}

void print(const std::vector<double>& v)
{
	for (const auto& i : v)
		std::cout << i << ' ';
	std::cout << std::endl;
}
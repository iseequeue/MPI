#include "file_work.h"

std::vector<std::vector<double>> read_config(std::string file_name, int n_particles)
{
	std::vector<std::vector<double>> ions;

	std::string line;

	std::ifstream in(file_name);
	if (in.is_open())
	{
		std::vector<double> coord;
		double input;
		for (auto i = 0u; i < n_particles; i++)
		{
			coord.clear();
			for (auto j = 0u; j < 3; j++)
			{
				in >> input;
				coord.push_back(input);
			}
			ions.push_back(coord);
		}
	}
	in.close();
	return ions;
}

void write_config(const std::string& file_name, const std::vector<std::vector<double>>& ions)
{
	std::ofstream out;
	out.open(file_name);
	out.precision(9);
	for (auto i : ions)
	{
		out << i[0] << ' ' << i[1] << ' ' << i[2] << std::endl;
	}

	out.close();

	return;
}

void write_vector(const std::string& file_name, const std::vector<double>& v)
{
	std::ofstream out;
	out.open(file_name);
	out.precision(9);

	for (const auto& i : v)
	{
		out << i << std::endl;
	}
	out.close();
	std::cout << "Writing to file completed" << std::endl;
}

void dump_step(const std::string& filename, const std::vector<std::vector<double>>& ions, const int& step, const double& L_calc_cell)
{
	std::ofstream out;
	out.open(filename, std::ios::app);
	out.precision(9);

	out << "ITEM: TIMESTEP\n";
	out << step << std::endl;
	out << "ITEM: NUMBER OF ATOMS\n";
	out << ions.size() << std::endl;
	out << "ITEM: BOX BOUNDS pp pp pp\n";
	out << 0 << ' ' << L_calc_cell << std::endl;
	out << 0 << ' ' << L_calc_cell << std::endl;
	out << 0 << ' ' << L_calc_cell << std::endl;
	out << "ITEM: ATOMS id x y z\n";
	for (auto i = 0u; i < ions.size(); i++)
	{
		out << i << ' ' << ions[i][0] << ' ' << ions[i][1] << ' ' << ions[i][2] << ' ' << std::endl;
	}

	out.close();
}
#pragma once
#include <vector>



#include "constants.h"
#include "supporting.h"

#include <cstdio>
#include "mpi.h"





double v1(const std::vector<double>& vector_ij, const double& L);

double v2(const std::vector<double>& vector_ij, const double& L_calc_cell);

double V0(const std::vector<std::vector<double>>& ions, const double& GAMMA);

double u(const int& i, const std::vector<std::vector<double>>& ions, const double& L_calc_cell);


double full_energy(const std::vector<std::vector<double>>& ions, const double& L_calc_cell, const double& GAMMA);

double tricky_energy(const std::vector<std::vector<double>>& ions_origin, const double& L_calc_cell, const double& GAMMA);


double spher_av_pot(const double& r_ij, const double& r_max);

std::vector<double> u_spher(const int& i, const std::vector<std::vector<double>>& ions_origin, const double& L_calc_cell);

std::vector<double> spherical_energy(const std::vector<std::vector<double>>& ions_origin, const double& L_calc_cell, const double& GAMMA);

std::vector<double> ewald_i(const int& i, const std::vector<double>& coord_i, const std::vector<std::vector<double>>& ions, const double& L_calc_cell);
std::vector<double> ewald_energy(const std::vector<std::vector<double>>& ions_origin, const double& L_calc_cell, const double& GAMMA);




double w_ewald(const std::vector<double>& r_ij, const double& L);

double w_i_ewald(const int& i, const std::vector<std::vector<double >>& ions, const double& L);


double W_ewald(const std::vector < std::vector<double>>& ions, const double& L);
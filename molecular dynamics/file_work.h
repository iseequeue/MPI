#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <fstream>



std::vector<std::vector<double>> read_config(std::string file_name, int n_particles);

void write_config(const std::string& file_name, const std::vector<std::vector<double>>& ions);

void write_vector(const std::string& file_name, const std::vector<double>& v);

void dump_step(const std::string& filename, const std::vector<std::vector<double>>& ions, const int& step, const double& L_calc_cell);
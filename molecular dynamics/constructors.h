#pragma once
#include <vector>
#include <cmath>
#include <random>

#include "constants.h"
#include "mpi.h"
#include <iostream>
#include "supporting.h"


double simple_cubic(std::vector<std::vector<double>>& ions, int elementary_max);

double base_cubic(std::vector<std::vector<double>>& ions, int elementary_max);

double face_cubic(std::vector<std::vector<double>>& ions, int elementary_max);

double random_state(std::vector<std::vector<double>>& ions, int n_particles);
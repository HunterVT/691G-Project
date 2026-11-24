//
// Created by hunter-vt on 11/19/26.
//

#include <iostream>
#include <vector>
#include <algorithm>
#include <array>

enum class State {
    A_PLUS, C_PLUS, G_PLUS, D_PLUS,
    A_MINUS, C_MINUS, G_MINUS, D_MINUS,
};

int main() {

    double in_island_prob = 0.02;
    double not_island = 0.98;

    // State ordering: A+, C+, G+, T+. A-, C-, G-, T-
    std::array<double, 8> initial_state_probabilities = {
        0.005, 0.005, 0.005, 0.005,
        0.245, 0.245, 0.245, 0.245
    };


    std::array<std::array<double, 8>, 8> transition_matrix{};


    return 0;
}
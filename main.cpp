//
// Created by hunter-vt on 11/19/26.
//

#include <iostream>
#include <vector>
#include <algorithm>
#include <array>
#include <cmath>
#include <string>

enum class State {
    A_PLUS, C_PLUS, G_PLUS, D_PLUS,
    A_MINUS, C_MINUS, G_MINUS, D_MINUS,
};

static void print_viterbi_results(
    const std::string& sequence,
    const std::vector<int>& optimal_path,
    const double final_max_score
) {
    constexpr std::array<char, 4> nucleotides = {
        'A', 'C', 'G', 'T'
    };

    // Print observed sequence.
    std::cout << "Observed sequence: ";
    for (const char nucleotide : sequence) {
        std::cout << nucleotide << "  ";
    }
    std::cout << '\n';

    // Print optimal hidden state sequence.
    std::cout << "Optimal path:      ";
    for (const int state : optimal_path) {
        const int nucleotide_idx = state % 4;
        const int region_idx = state / 4;

        std::cout << nucleotides[nucleotide_idx]
                  << (region_idx == 0 ? '+' : '-')
                  << ' ';
    }
    std::cout << '\n';

    // Print predicted CpG island membership.
    std::cout << "Island membership: ";
    for (const int state : optimal_path) {
        std::cout << (state / 4 == 0 ? '+' : '-')
                  << "  ";
    }
    std::cout << '\n';

    // Print probability of the optimal path.
    std::cout << "Final probability: "
              << std::scientific
              << final_max_score << '\n';
}

static int nucleotide_to_idx (const char nucleotide) {
    switch (nucleotide) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        default:
            return -1;
    }
}

int main() {

    // Probabilities of transitioning/staying to in/out of a CpG island
    const std::array<std::array<double, 2>, 2> region_transition_probs = {
        0.99, 0.01,
        0.001, 0.999
    };

    // Nucleotide ordering: A, C, G, T
    // island_transitions is probability of observing next nucleotide within an island
    // conditional on the current nucleotide
    const std::array<std::array<double, 4>, 4> island_transitions = {
        0.15, 0.35, 0.35, 0.15,
        0.15, 0.20, 0.50, 0.15,
        0.15, 0.35, 0.35, 0.15,
        0.15, 0.35, 0.35, 0.15
    };

    // background_transitions is probability of observing next nucleotide NOT within an island
    // conditional on the current nucleotide
    const std::array<std::array<double, 4>, 4> background_transitions = {
        0.30, 0.20, 0.20, 0.30,
        0.40, 0.12, 0.08, 0.40,
        0.30, 0.20, 0.20, 0.30,
        0.30, 0.20, 0.20, 0.30,
    };

    // State ordering: A+, C+, G+, T+. A-, C-, G-, T-
    std::array<double, 8> initial_state_probabilities = {
        0.005, 0.005, 0.005, 0.005,
        0.245, 0.245, 0.245, 0.245
    };

    std::array<std::array<double, 8>, 8> transition_matrix{};

    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {

            // Get current region and nucleotide indices
            const int current_region = i / 4;
            const int current_nucleotide = i % 4;

            // Get next (destination) region and nucleotide indices
            const int destination_region = j / 4;
            const int destination_nucleotide = j % 4;

            const double region_transition_prob = region_transition_probs[current_region][destination_region];
            double nucleotide_transition_prob = 0;


            // Calculate and store transition probability
            if (destination_region == 0) {
                nucleotide_transition_prob = island_transitions[current_nucleotide][destination_nucleotide];
            }
            else if (destination_region == 1) {
                nucleotide_transition_prob = background_transitions[current_nucleotide][destination_nucleotide];
            }
            transition_matrix[i][j] = region_transition_prob * nucleotide_transition_prob;
        }
    }

    // Verify that our transition matrix is sound
    constexpr double epsilon = 1e-9;
    for (int i = 0; i < 8; ++i) {
        double running_sum = 0;
        for (int j = 0; j < 8; ++j) {
            if (!std::isfinite(transition_matrix[i][j]) || transition_matrix[i][j] > 1 || transition_matrix[i][j] < 0) {
                std::cerr << "Invalid transition probability matrix!" << std::endl;
                std::cerr << "Entry/Entries are not finite!" << std::endl;
                return 1;
            }
            running_sum += transition_matrix[i][j];
        }
        if (!(std::abs(running_sum - 1) < epsilon)) {
            std::cerr << "Invalid transition probability matrix!" << std::endl;
            return 1;
        }
    }

    // Emission matrix
    std::array<std::array<double, 4>, 8> emission_matrix{};

    for (int i = 0; i < 8; ++i) {
        // Determine the nucleotide index associated with state i.
        const int nucleotide_idx = i % 4;

        // Set the corresponding emission probability to 1.
        emission_matrix[i][nucleotide_idx] = 1;

    }

    const std::string sequence =
        "ATATATATATATATATATAT"
        "CGCGCGCGCGCGCGCGCGCG"
        "ATATATATATATATATATAT";
    std::vector<std::array<double, 8>> viterbi_matrix(sequence.size());

    std::array<int, 8> initial_row{};
    initial_row.fill(-1);
    std::vector<std::array<int, 8>> traceback_matrix(sequence.size(), initial_row);
    if (sequence.empty()) {
        return 1;
    }
    const int first_nucleotide_idx = nucleotide_to_idx(sequence[0]);
    if (first_nucleotide_idx == -1) {
        return 1;
    }
    for (int j = 0; j < 8; ++j) {
        viterbi_matrix[0][j] = initial_state_probabilities[j] * emission_matrix[j][first_nucleotide_idx];
    }


    for (std::size_t t = 1; t < sequence.size(); ++t) {
        const int curr_nucleotide_idx = nucleotide_to_idx(sequence[t]);
        if (curr_nucleotide_idx == -1) {
            return 1;
        }
        for (int j = 0; j < 8; ++j) {

            double max_score = 0.0;
            int best_predecessor = -1;

            for (int i = 0; i < 8; ++i) {
                const double candidate_score = viterbi_matrix[t-1][i] * transition_matrix[i][j];
                if (candidate_score > max_score) {
                    max_score = candidate_score;
                    best_predecessor = i;
                }
            }

            viterbi_matrix[t][j] = emission_matrix[j][curr_nucleotide_idx] * max_score;

            if (viterbi_matrix[t][j] > 0) {
                traceback_matrix[t][j] = best_predecessor;
            }
        }
    }

    double final_max_score = 0.0;
    int final_state = -1;

    for (int j = 0; j < 8; ++j) {
        const double candidate_score = viterbi_matrix.back()[j];

        if (candidate_score > final_max_score) {
            // Update final maximum score
            final_max_score = candidate_score;
            // Update final state.
            final_state = j;
        }
    }

    if (final_state < 0 || final_state > 7) {
        std::cerr << "No valid final state found. "
                  << "Possible numerical underflow.\n";
        return 1;
    }

    std::vector<int> optimal_path(sequence.size(), -1);
    optimal_path.back() = final_state;
    for (std::size_t t = sequence.size() - 1; t > 0; --t) {
        // Obtain the optimal state at position t.
        const int optimal_state = optimal_path[t];

        // Look up its predecessor in traceback_matrix.
        int predecessor = traceback_matrix[t][optimal_state];

        if (predecessor < 0 || predecessor > 7) {
            return 1;
        }

        // Store the predecessor at position t - 1
        // in optimal_path.
        optimal_path[t-1] = predecessor;

    }

    print_viterbi_results(sequence, optimal_path, final_max_score);

    return 0;
}
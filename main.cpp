//
// Created by hunter-vt on 11/19/26.
//

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <iomanip>
#include <string>
#include <limits>
#include <optional>
#include <algorithm>
#include<fstream>
#include <stdexcept>
#include <chrono>

static std::string read_sequence(const std::string& filepath) {

    // Open the file
    std::ifstream file(filepath);
    // Verify that the file opened successfully
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    std::string sequence;

    // Read the nucleotide sequence
    // Reading nucleotide failed
    if (!std::getline(file, sequence) || sequence.empty()) {
        throw std::runtime_error("Could not read nucleotide");
    }

    // Return the sequence
    return sequence;
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


template <typename T>
struct HMM {
    // Initial state probabilities
    // State ordering: A+, C+, G+, T+. A-, C-, G-, T-
    std::array<T, 8> initial_state_probabilities{};

    // Transition probability matrix
    std::array<std::array<T, 8>, 8> transition_matrix{};

    // Emission probability matrix
    std::array<std::array<T, 4>, 8> emission_matrix{};
};

template <typename T>
static HMM<T> construct_synthetic_hmm() {
    HMM<T> hmm{};

    // Define constituent transition matrices.
    // Probabilities of transitioning/staying to in/out of a CpG island
    const std::array<std::array<T, 2>, 2> region_transition_probs = {
        static_cast<T>(0.99L), static_cast<T>(0.01L),
        static_cast<T>(0.001L), static_cast<T>(0.999L)
    };

    const std::array<std::array<T, 4>, 4> island_transitions = {
        static_cast<T>(0.15L), static_cast<T>(0.35L), static_cast<T>(0.35L), static_cast<T>(0.15L),
        static_cast<T>(0.15L), static_cast<T>(0.20L), static_cast<T>(0.50L), static_cast<T>(0.15L),
        static_cast<T>(0.15L), static_cast<T>(0.35L), static_cast<T>(0.35L), static_cast<T>(0.15L),
        static_cast<T>(0.15L), static_cast<T>(0.35L), static_cast<T>(0.35L), static_cast<T>(0.15L)
    };

    // Background_transitions is probability of observing next nucleotide NOT within an island
    // conditional on the current nucleotide
    const std::array<std::array<T, 4>, 4> background_transitions = {
        static_cast<T>(0.30L), static_cast<T>(0.20L), static_cast<T>(0.20L), static_cast<T>(0.30L),
        static_cast<T>(0.40L), static_cast<T>(0.12L), static_cast<T>(0.08L), static_cast<T>(0.40L),
        static_cast<T>(0.30L), static_cast<T>(0.20L), static_cast<T>(0.20L), static_cast<T>(0.30L),
        static_cast<T>(0.30L), static_cast<T>(0.20L), static_cast<T>(0.20L), static_cast<T>(0.30L),
    };

    // Initialize state probabilities.
    // State ordering: A+, C+, G+, T+. A-, C-, G-, T-
    hmm.initial_state_probabilities = {
        static_cast<T>(0.005L), static_cast<T>(0.005L), static_cast<T>(0.005L), static_cast<T>(0.005L),
        static_cast<T>(0.245L), static_cast<T>(0.245L), static_cast<T>(0.245L), static_cast<T>(0.245L)
    };

    // Construct the full transition matrix.
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {

            // Get current region and nucleotide indices
            const int current_region = i / 4;
            const int current_nucleotide = i % 4;

            // Get next (destination) region and nucleotide indices
            const int destination_region = j / 4;
            const int destination_nucleotide = j % 4;

            const T region_transition_prob = region_transition_probs[current_region][destination_region];
            T nucleotide_transition_prob = 0;


            // Calculate and store transition probability
            if (destination_region == 0) {
                nucleotide_transition_prob = island_transitions[current_nucleotide][destination_nucleotide];
            }
            else if (destination_region == 1) {
                nucleotide_transition_prob = background_transitions[current_nucleotide][destination_nucleotide];
            }
            hmm.transition_matrix[i][j] = region_transition_prob * nucleotide_transition_prob;
        }
    }

    // Construct the emission matrix.
    for (int i = 0; i < 8; ++i) {
        // Determine the nucleotide index associated with state i.
        const int nucleotide_idx = i % 4;

        // Set the corresponding emission probability to 1.
        hmm.emission_matrix[i][nucleotide_idx] = 1;

    }
    return hmm;
}

template <typename T>
static bool validate_hmm(const HMM<T>& hmm) {
    // Verify that our transition matrix is sound
    constexpr T epsilon = 16 * std::numeric_limits<T>::epsilon();

    // Check transition matrix
    for (int i = 0; i < 8; ++i) {
        T running_sum = 0;
        for (int j = 0; j < 8; ++j) {
            if (!std::isfinite(hmm.transition_matrix[i][j]) || hmm.transition_matrix[i][j] > 1 || hmm.transition_matrix[i][j] < 0) {
                std::cerr << "Invalid transition probability matrix!" << std::endl;
                std::cerr << "Entry/Entries are not finite!" << std::endl;
                return false;
            }
            running_sum += hmm.transition_matrix[i][j];
        }
        if (!(std::abs(running_sum - 1) < epsilon)) {
            std::cerr << "Invalid transition probability matrix!" << std::endl;
            return false;
        }
    }

    // Check initial matrix
    T init_sum = 0;
    for (int i = 0; i < 8; ++i) {
        if (!std::isfinite(hmm.initial_state_probabilities[i]) || hmm.initial_state_probabilities[i] > 1 || hmm.initial_state_probabilities[i] < 0) {
            std::cerr << "Invalid transition probability matrix!" << std::endl;
            std::cerr << "Entry/Entries are not finite!" << std::endl;
            return false;
        }
        init_sum += hmm.initial_state_probabilities[i];
    }
    if (!(std::abs(init_sum - 1) < epsilon)) {
        std::cerr << "Invalid transition probability matrix!" << std::endl;
        return false;
    }
    // Check emission matrix
    for (int i = 0; i < 8; ++i) {
        T emit_row_sum = 0.0;
        for (int k = 0; k < 4; ++k) {
            if (!std::isfinite(hmm.emission_matrix[i][k]) || hmm.emission_matrix[i][k] > 1 || hmm.emission_matrix[i][k] < 0) {
                std::cerr << "Invalid transition probability matrix!" << std::endl;
                std::cerr << "Entry/Entries are not finite!" << std::endl;
                return false;
            }
            emit_row_sum += hmm.emission_matrix[i][k];
        }
        if (!(std::abs(emit_row_sum - 1) < epsilon)) {
            std::cerr << "Invalid transition probability matrix!" << std::endl;
            return false;
        }
    }

    return true;
}

struct ObserveUnderflow {
    // First position containing a subnormal candidate score
    std::optional<std::size_t> first_subnormal_position{};

    // First position where a positive multiplication rounds to zero
    std::optional<std::size_t> first_underflow_position{};

    // First position where all Viterbi scores are zero
    std::optional<std::size_t> first_all_zero_position{};

    // Total number of subnormal candidate scores
    std::size_t subnormal_count = 0;

    // Total number of multiplications that underflow to zero
    std::size_t underflow_count = 0;
};

template <typename T>
struct ViterbiResult {
    // Final maximum score
    T final_score = 0;
    // Final state
    int final_state = -1;
    // Optimal path
    std::vector<int> optimal_path{};
    // Whether decoding was successful
    bool success = false;

    ObserveUnderflow underflow_observations;
};

template <typename T>
static ViterbiResult<T> viterbi(
    const HMM<T>& hmm,
    const std::string& sequence
) {
    ViterbiResult<T> result{};

    // Initialization
    std::vector<std::array<T, 8>> viterbi_matrix(sequence.size());
    std::array<int, 8> initial_row{};
    initial_row.fill(-1);
    std::vector<std::array<int, 8>> traceback_matrix(sequence.size(), initial_row);
    if (sequence.empty()) {
        return result;
    }
    const int first_nucleotide_idx = nucleotide_to_idx(sequence[0]);
    if (first_nucleotide_idx == -1) {
        return result;
    }
    for (int j = 0; j < 8; ++j) {
        viterbi_matrix[0][j] = hmm.initial_state_probabilities[j] * hmm.emission_matrix[j][first_nucleotide_idx];
    }

    // Recursion
    for (std::size_t t = 1; t < sequence.size(); ++t) {
        const int curr_nucleotide_idx = nucleotide_to_idx(sequence[t]);
        if (curr_nucleotide_idx == -1) {
            return result;
        }
        for (int j = 0; j < 8; ++j) {

            T max_score = 0.0;
            int best_predecessor = -1;

            for (int i = 0; i < 8; ++i) {
                const T previous_score = viterbi_matrix[t-1][i];
                const T transition_probability = hmm.transition_matrix[i][j];
                const T candidate_score = previous_score * transition_probability;

                // Check for underflow
                if (hmm.emission_matrix[j][curr_nucleotide_idx] > 0) {

                    // First check for and detect subnormal candidate scores
                    if (std::fpclassify(candidate_score) == FP_SUBNORMAL) {
                        ++result.underflow_observations.subnormal_count;

                        // Record first occurrence only
                        if (!result.underflow_observations.first_subnormal_position.has_value()) {
                            result.underflow_observations.first_subnormal_position = t;
                        }
                    }

                    // Detect underflow to zero
                    if (previous_score > 0 && transition_probability > 0 && std::fpclassify(candidate_score) == FP_ZERO) {
                        ++result.underflow_observations.underflow_count;

                        if (!result.underflow_observations.first_underflow_position.has_value()) {
                            result.underflow_observations.first_underflow_position = t;
                        }
                    }
                }

                if (candidate_score > max_score) {
                    max_score = candidate_score;
                    best_predecessor = i;
                }
            }
            viterbi_matrix[t][j] = hmm.emission_matrix[j][curr_nucleotide_idx] * max_score;

            if (viterbi_matrix[t][j] > 0) {
                traceback_matrix[t][j] = best_predecessor;
            }
        }
        const T row_max = *std::max_element(
viterbi_matrix[t].begin(),
viterbi_matrix[t].end()
);
        if (row_max == 0) {
            if (!result.underflow_observations.first_all_zero_position.has_value()) {
                result.underflow_observations.first_all_zero_position = t;
            }
        }
    }
    // Termination
    for (int j = 0; j < 8; ++j) {
        const T candidate_score = viterbi_matrix.back()[j];

        if (candidate_score > result.final_score) {
            // Update final maximum score
            result.final_score = candidate_score;
            // Update final state.
            result.final_state = j;
        }
    }

    if (result.final_state < 0 || result.final_state > 7) {
        std::cerr << "No valid final state found. "
                  << "Possible numerical underflow.\n";
        return result;
    }

    // Traceback
    result.optimal_path.resize(sequence.size(), -1);
    result.optimal_path.back() = result.final_state;
    for (std::size_t t = sequence.size() - 1; t > 0; --t) {
        // Obtain the optimal state at position t.
        const int optimal_state = result.optimal_path[t];

        // Look up its predecessor in traceback_matrix.
        int predecessor = traceback_matrix[t][optimal_state];

        if (predecessor < 0 || predecessor > 7) {
            result.optimal_path.clear();
            return result;
        }

        // Store the predecessor at position t - 1
        // in optimal_path.
        result.optimal_path[t-1] = predecessor;

    }
    result.success = true;
    return result;
}

// Helper function to print Viterbi run results
template <typename T>
static void print_viterbi_results(
    const std::string& sequence,
    const ViterbiResult<T>& result
) {
    if (!result.success) {
        std::cout << "Decoding unsuccessful.\n";
        return;
    }
    std::setprecision(std::numeric_limits<T>::max_digits10);
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
    for (const int state : result.optimal_path) {
        const int nucleotide_idx = state % 4;
        const int region_idx = state / 4;

        std::cout << nucleotides[nucleotide_idx]
                  << (region_idx == 0 ? '+' : '-')
                  << ' ';
    }
    std::cout << '\n';

    // Print predicted CpG island membership
    std::cout << "Island membership: ";
    for (const int state : result.optimal_path) {
        std::cout << (state / 4 == 0 ? '+' : '-')
                  << "  ";
    }
    std::cout << '\n';

    // Print probability of the optimal path
    std::cout << "Final probability: "
              << std::scientific
              << result.final_score << '\n';
}

// Helper function to print underflow occurrences
static void print_underflow_report(
    const ObserveUnderflow& observations
) {
    std::cout << "First subnormal position: ";
    if (observations.first_subnormal_position.has_value()) {
        std::cout << observations.first_subnormal_position.value();
    }
    else {
        std::cout << "None";
    }
    std::cout << '\n';

    std::cout << "First underflow position: ";
    if (observations.first_underflow_position.has_value()) {
        std::cout << observations.first_underflow_position.value();
    }
    else {
        std::cout << "None";
    }
    std::cout << '\n';

    std::cout << "First all-zero position: ";
    if (observations.first_all_zero_position.has_value()) {
        std::cout << observations.first_all_zero_position.value();
    }
    else {
        std::cout << "None";
    }
    std::cout << '\n';

    std::cout << "Subnormal count: "
              << observations.subnormal_count << '\n';

    std::cout << "Underflow count: "
              << observations.underflow_count << '\n';
}

int main(int argc, char* argv[]) {

    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <path_to_sequence> <float, double, or 'long_double'> \n";
        return 1;
    }

    const std::string sequence = read_sequence(argv[1]);

    const std::string data_type = argv[2];

    if (data_type != "float" && data_type != "double" && data_type != "long_double") {
        std::cerr << "Invalid data type: Use float, double, or long_double" << std::endl;
        return 1;
    }

    const auto hmm_float = construct_synthetic_hmm<float>();
    const auto hmm_double = construct_synthetic_hmm<double>();
    const auto hmm_long_double = construct_synthetic_hmm<long double>();

    if (!validate_hmm(hmm_float) ||
        !validate_hmm(hmm_double) ||
        !validate_hmm(hmm_long_double)) {
        return 1;
        }
/**
    const std::string sequence =
        "ATATATATATATATATATAT"
        "CGCGCGCGCGCGCGCGCGCG"
        "ATATATATATATATATATAT";
    */

    // const std::string sequence(20000, 'A');

    if (data_type == "float") {

        // Get duration and run viterbi
        const auto float_start = std::chrono::steady_clock::now();
        const auto result_float = viterbi(hmm_float, sequence);
        const auto float_end = std::chrono::steady_clock::now();
        std::chrono::duration<double, std::milli> float_duration_ms = float_end - float_start;

        // Print results
        std::cout << "\n--- Float ---\n";
        std::cout << "\nInput Size: " + std::to_string(sequence.size()) + " base pairs\n";
        print_viterbi_results(sequence, result_float);
        print_underflow_report(result_float.underflow_observations);
        std::cout << "Duration: " + std::to_string(float_duration_ms.count()) + "ms" << std::endl;
    }

    if (data_type == "double") {

        // Get duration and run viterbi
        const auto double_start = std::chrono::steady_clock::now();
        const auto result_double = viterbi(hmm_double, sequence);
        const auto double_end = std::chrono::steady_clock::now();
        std::chrono::duration<double, std::milli> double_duration_ms = double_end - double_start;

        // Print results
        std::cout << "\n--- Double ---\n";
        std::cout << "\nInput Size: " + std::to_string(sequence.size()) + " base pairs\n";
        print_viterbi_results(sequence, result_double);
        print_underflow_report(result_double.underflow_observations);
        std::cout << "Duration: " + std::to_string(double_duration_ms.count()) + "ms" << std::endl;



    }

    if (data_type == "long_double") {

        // Get duration and run viterbi
        const auto long_double_start = std::chrono::steady_clock::now();
        const auto result_long_double = viterbi(hmm_long_double, sequence);
        const auto long_double_end = std::chrono::steady_clock::now();
        std::chrono::duration<double, std::milli> long_double_duration_ms = long_double_end - long_double_start;

        // Print results
        std::cout << "\n--- Long Double ---\n";
        std::cout << "\nInput Size: " + std::to_string(sequence.size()) + + " base pairs\n";
        print_viterbi_results(sequence, result_long_double);
        print_underflow_report(result_long_double.underflow_observations);
        std::cout << "Duration: " + std::to_string(long_double_duration_ms.count()) + "ms" << std::endl;
    }








    return 0;
}
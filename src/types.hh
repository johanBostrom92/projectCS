#pragma once
#include "parameters.hh"
#include <vector>
#include <atomic>
#include <string>
#include <array>
#include <random>


/**
 * The status one agent can have
 * S = Susceptible - Not currently infected but can be by an Infected or Asymptotic agent
 * A = Asymptotic - Currently infected but not showing symptoms, can infect Susceptible agents so they become either Infected or Asymptotic
 * I = Infected - Currently infected and showing symptoms, can infect Susceptible agents so they become either Infected or Asymptotic
 * V = Vaccinated - Currently vaccinated, cant infect nor be infected
 * R = Recovered - Has been infected but recovered, cant infect nor be infected.
 *
 * WARNING: It is necessary for the code to work that all entries are given their default values (do not assign e.g. S=2),
 * and that the STATES_COUNT entry is always last.
 */
enum agent_status {
    S, A, I, V, R,
    STATES_COUNT // Keep this last, so it equals the number of states there are
};

/**
 * Defines a single agent or person in the model
 */
struct agent {
    agent_status status = S;
    int infection_radius = INFECTION_RADIUS;
    float infection_probability = INFECTION_PROBABILITY;
    int recovery_rate = RECOVERY_RATE;
    bool vaccination_progress = false;
    int vaccination_rate = VACCINATION_RATE;
    int DAY_INFECTED = 0;
};

/**
 * A group of agents interacting in a grid
 */
struct board {
    unsigned int dim;
    std::string name;

    std::vector<agent> agents;
    std::array<std::atomic_int, STATES_COUNT> status_counts;
    std::atomic_uint total_infections;
    std::vector<unsigned int> vaccination_weights;
    std::atomic_uint64_t vaccination_weight_sum;
    std::atomic_uint vaccinations_started;
    std::tuple<double, double> lat_long;
    std::vector<double> weights;

    /**
     * Creates a new square board
     * @param dim The size of the board along each axis
     * @param initial_infections The number of infected agents at creation
     * @param agent_types Specifications for what agents to populate the board with
     * @param name A name describing the board
     * @param A tuple of latitude and longitude of the board
     * @param gen Generator to use when generating initial infections
     */
    board(unsigned int dim, unsigned int initial_infections, const std::vector<agent_type>& agent_types, const std::string& name, std::tuple<double, double> lat_long, std::vector<double> weights,std::mt19937_64& gen);

    // std::atomic deletes these, so we need to redefine them.
    // Note that these are not atomic
    board(const board &);
    board(board &&);
    board &operator=(const board &);
    board &operator=(board &&);
};
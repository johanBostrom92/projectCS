#pragma once
#include "parameters.hh"
#include <vector>
#include <atomic>
#include <string>


/**
 * The status one agent can have
 * S = Susceptible - Not currently infected but can be by an Infected or Asymptotic agent
 * A = Asymptotic - Currently infected but not showing symptoms, can infect Susceptible agents so they become either Infected or Asymptotic
 * I = Infected - Currently infected and showing symptoms, can infect Susceptible agents so they become either Infected or Asymptotic
 * V = Vaccinated - Currently vaccinated, cant infect nor be infected
 * R = Recovered - Has been infected but recovered, cant infect nor be infected.
 */
enum agent_status {
    S, A, I, V, R
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
};

/**
 * A group of agents interacting in a grid
 */
struct board {
    unsigned int dim;
    double weight;
    std::string name;

    std::vector<agent> agents;
    std::atomic_int sus;
    std::atomic_int rem;
    std::atomic_int inf;
    std::atomic_int asymp;
    std::atomic_int vacc;
    std::vector<unsigned int> vaccination_weights;
    std::atomic_uint64_t vaccination_weight_sum;
    std::atomic_uint vaccinations_started;

    /**
     * Creates a new square board
     * @param dim The size of the board along each axis
     * @param initial_infections The number of infected agents at creation
     * @param agent_types Specifications for what agents to populate the board with
     */
    board(unsigned int dim, unsigned int initial_infections, const std::vector<agent_type> agent_types, double weight, std::string name);

    // std::atomic deletes these, so we need to redefine them.
    // Note that these are not atomic
    board(const board &);
    board(board &&);
    board &operator=(const board &);
    board &operator=(board &&);
};
#pragma once
#include "parameters.hh"
#include <vector>
#include <atomic>
#include <string>

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
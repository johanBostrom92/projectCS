#pragma once
#include <vector>

// TODO: add some descriptions
#define INFECTION_RADIUS 2
#define RECOVERY_RATE 5
#define INFECTION_PROBABILITY 100
#define MAKE_ASYM 50 //The probability that an infected individual is asymptomatic
#define ASYM_INF_PROB 10 //The probability for asymptomatic carriers to infect others
#define VACCINATION_RATE 28
#define VACCINATION_EFFICACY 75
#define DIM 3
#define MAX_TIME 100
#define STARTER_AGENTS 1
#define QUARANTINE_START 5
#define ENABLE_QUARANTINE false
#define LAMBDA 2.5
#define ONLY_ELIGIBLE false // If TRUE, it chooses only an eligible target to infect. If FALSE, any target can be chosen/tried.

/**
 * Defines a type of agent, by its infection radius and how often it should occur.
 * The infection radius is defined by a normal distribution.
 */
struct agent_type {
    /** The mean of the infection radius */
    float radius_mean;
    /** The standard deviation of the infection radius. Can be set to 0 to have all agents take radius_mean */
    float radius_stddev;
    /** The weight of the type, determining its relative frequency. A type with weight 2 will be twice as common as one with weight 1, etc. */
    double weight;
};

/**
 * Gives all agent types to use for the simulation
 */
static const std::vector<agent_type> AGENT_TYPES = {
    agent_type{ INFECTION_RADIUS, 0, 1 }
};
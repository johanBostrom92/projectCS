#pragma once
#include <vector>

// TODO: add some descriptions
#define INFECTION_RADIUS 50
#define RECOVERY_RATE 14
#define INFECTION_PROBABILITY 25
#define MAKE_ASYM 10 //The probability that an infected individual is asymptomatic
#define ASYM_INF_PROB 10 //The probability for asymptomatic carriers to infect others
#define VACCINATION_RATE 28
#define VACCINATION_EFFICACY 75
#define VACCINATION_START 10     // The time at which to begin vaccinating people
#define VACCINATIONS_PER_DAY 5000  // The number of people we can vaccinate per day
#define DIM 1000
#define MAX_TIME 200
#define STARTER_AGENTS 4
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

static_assert(DIM >= INFECTION_RADIUS);
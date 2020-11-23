#pragma once
#include <vector>

// TODO: add some descriptions
#define INFECTION_RADIUS 50
#define RECOVERY_RATE 14
#define INFECTION_PROBABILITY 25
#define DIM 1000
#define MAX_TIME 140
#define STARTER_AGENTS 4
#define QUARANTINE 5
#define Q_FLAG false
#define LAMBDA 2.5
#define ONE_CHANCE false

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
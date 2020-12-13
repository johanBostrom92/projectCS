#pragma once
#include <vector>


#define INFECTION_RADIUS 2 //The radius in which an agent can infect another
#define RECOVERY_RATE 5 //The time it takes for an agent to become Recovered
#define INFECTION_PROBABILITY 25 //Probability for an agent to infect another
#define MAKE_ASYM 50 //The probability that an infected individual is asymptomatic
#define ASYM_INF_PROB 10 //The probability for asymptomatic carriers to infect others

#define DIM 7 //The population of the community
#define MAX_TIME 100 //The amount of timesteps to run the simulation
#define STARTER_AGENTS 1 //The amount of agents which starts the simulation infected
#define QUARANTINE_START 5 //The timestep until agents starts quarantining
#define ENABLE_QUARANTINE false //The flag to enable quarantine
#define LAMBDA 2.5  // How fast the INFECTION_RADIUS decrease in case of a lock down
#define ONLY_ELIGIBLE false // If TRUE, it chooses only an eligible target to infect. If FALSE, any target can be chosen/tried.
#define PLOT true   // if we want to do plots (increases execution time drastically)
#define SCALE 1 //scale for geobubble plot. 

#define VACCINATION_RATE 28 //The time it takes for the vaccine to work
#define VACCINATION_EFFICACY 75 //The efficiency of the vaccine
#define VACCINATION_START 10     // The time at which to begin vaccinating people
#define VACCINATIONS_PER_DAY 5000  // The number of people we can vaccinate per day

enum class vaccination_strategy {
    UNIFORM,       // Every agents has equal change of being vaccinated
    HIGH_DENSITY,  // Agents in areas with a high density of infected agents have a *higher* chance of being vaccinated
    LOW_DENSITY    // Agents in areas with a high density of infected agents have a *lower* chance of being vaccinated
};

constexpr vaccination_strategy VACC_STRAT = vaccination_strategy::LOW_DENSITY;


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

static_assert(DIM >= INFECTION_RADIUS, "The infection radius cannot be greater than the board dimension");
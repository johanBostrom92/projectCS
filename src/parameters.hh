#pragma once
#include <vector>


#define INFECTION_RADIUS 50 //The radius in which an agent can infect another
#define RECOVERY_RATE 14 //The time it takes for an agent to become Recovered
#define INFECTION_PROBABILITY 25 //Probability for an agent to infect another
#define MAKE_ASYM 50 //The probability that an infected individual is asymptomatic
#define ASYM_INF_PROB 10 //The probability for asymptomatic carriers to infect others

#define CITIES_INPUT "cities_swe.csv"

#define MAX_TIME 365 //The amount of timesteps to run the simulation
#define STARTER_AGENTS 10 //The amount of agents which starts the simulation infected

//The quarantine variables
#define ENABLE_QUARANTINE true //The flag to enable quarantine
#define QUARANTINE_START 30 //The specific timestep when quarantining begins
#define QUARANTINE_END 60 //The specific timestep when quarantining stops
#define QUARANTINE_EFFICACY 100 //The efficiency of the quarantine (how much of the population that actually follows the quarantine)
#define LAMBDA 2.5  // How fast the INFECTION_RADIUS decrease in case of a lock down

#define ONLY_ELIGIBLE false // If TRUE, it chooses only an eligible target to infect. If FALSE, any target can be chosen/tried.
#define PLOT true   // if we want to do plots (increases execution time drastically)
#define SCALE 1 //scale for geobubble plot. 
#define SWAP_AMOUNT 1500 // sets the amount of agents to swap each time unit. 

//The two variables that determines reinfection times. 
//To disable reinfection, set min & max value to max simulation time plus one (MAX_TIME+1).
#define RECOVERED_MIN_THRESHOLD 20 // The minimum amount of immunity-days after a recovered agent can be reinfected again.
#define RECOVERED_MAX_THRESHOLD MAX_TIME // The maxmimum amount of immunity-days after a recovered agent can be reinfected again. Default is the max simulation time

#define VACCINATION_RATE 28 //The time it takes for the vaccine to work
#define VACCINATION_EFFICACY 75 //The efficiency of the vaccine
#define VACCINATION_START 180    // The time at which to begin vaccinating people
#define VACCINATIONS_PER_DAY 1000  // The number of people we can vaccinate per day

enum class vaccination_strategy {
    UNIFORM,       // Every agents has equal change of being vaccinated
    HIGH_DENSITY,  // Agents in areas with a high density of infected agents have a *higher* chance of being vaccinated
    LOW_DENSITY    // Agents in areas with a high density of infected agents have a *lower* chance of being vaccinated
};

constexpr vaccination_strategy VACC_STRAT = vaccination_strategy::HIGH_DENSITY;


/**
 * Defines a type of agent, by its infection radius, infection probability and how often it should occur.
 * The infection radius is defined by a normal distribution.
 */
struct agent_type {
    /** How far away agents may interact with other agents */
    float infection_radius;
    /** The probability of an infected agent infecting another agent it interacts with */
    unsigned int infection_probability;
    /** The weight of the type, determining its relative frequency. A type with weight 2 will be twice as common as one with weight 1, etc. */
    double weight;
};

/**
 * Gives all agent types to use for the simulation
 */
static const std::vector<agent_type> AGENT_TYPES = {
    agent_type{ INFECTION_RADIUS, INFECTION_PROBABILITY, 1 }
};
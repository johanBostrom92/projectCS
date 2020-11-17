#include <vector>
#include <omp.h>
#include <cstdlib>
#include <time.h>
#include <iostream>
#include <atomic>
#include <random>
#include <numeric>

#define INFECTION_RADIUS 1
#define RECOVERY_RATE 14
#define INFECTION_PROBABILITY 100
#define DIM 3
#define MAX_TIME 3
#define N_PATIENT_ZERO 4

// Defines a type of agent, by its infection radius and how often it should occur.
struct agent_type {
    float radius_mean;
    float radius_stddev;
    double weight;
};

static const std::vector<agent_type> agent_types = {
    agent_type{ 50, 0, 1 }
};

enum agent_status {
        S, I, R
};

struct agent {
    agent_status status = S;
    int infection_radius = INFECTION_RADIUS;
    float infection_probability = INFECTION_PROBABILITY;
    //int social_distancing;
    //bool inuarantine;
    //bool symptoms;
    int recovery_rate = RECOVERY_RATE;
};

struct board {
    unsigned int dim = DIM;
    std::vector<agent> agents;
};

void print_board(board& b){
    for(int i = 0; i < DIM*DIM; i++){
        if(i % DIM == 0){
            std::endl (std::cout);
        }
        std::cout <<  b.agents[i].status << " ";
    }
}

board generate_board() {
    board b = {
        DIM,
        std::vector<agent>(DIM*DIM)
    };
    // Set initial infections
    for (int i = 0; i < N_PATIENT_ZERO; i++) {
        int pZ = std::rand() % DIM*DIM;
        b.agents[pZ].status = I;
    }

    // Get the sum of all type weights
    double weight_sum = 0.0f;
    for (int i = 0; i < agent_types.size(); i++) {
        auto& type = agent_types[i];
        weight_sum += type.weight;
    }

    // Randomize an agent type for every agent, then generate and assign a radius using that type
    std::default_random_engine rand_generator;
    std::uniform_real_distribution type_dist(0.0, weight_sum);
    for(auto& agent : b.agents) {
        double type_val = type_dist(rand_generator);
        double cumulative_weight = 0.0f;
        for (int i = 0; i < agent_types.size(); i++) {
            auto& type = agent_types[i];
            if (cumulative_weight + type.weight >= type_val) {
                if (type.radius_stddev == 0) {  // Uses a fixed radius for agents of this type
                    agent.infection_radius = type.radius_mean;
                } else {  // Generates a radius from a normal distribution
                    agent.infection_radius = std::normal_distribution<float>(type.radius_mean, type.radius_stddev)(rand_generator);
                }
                break;
            }
            cumulative_weight += type.weight;
        }
    }
    return b;
}

int main() {
    srand(time(NULL));

    board previous = generate_board();

    board current = previous;
	std::atomic_int sus = DIM*DIM-1;
	std::atomic_int rem = 0;
	std::atomic_int inf = 1;
	for (unsigned int t = 0; t < MAX_TIME; t++) { //Loop tracking time
        #pragma omp parallel for collapse(2)
        for (unsigned int y = 0; y < DIM; y++) {
            for (unsigned int x = 0; x < DIM; x++) {

                //{ Update agent
                agent& self = previous.agents[y * DIM + x];
                agent& currentSelf = current.agents[y * DIM + x];
                if (self.status == I) { //The infected checks for susceptible neighbors within the box.
                    currentSelf.recovery_rate--;
                    if(self.recovery_rate == 0){
                        currentSelf.status = R;
						rem++;
						inf--;
						continue;
                    }
                    //{ Check neighbours
                    for (int y_box = -self.infection_radius; y_box <= self.infection_radius; y_box++) {
                        for (int x_box = -self.infection_radius; x_box <= self.infection_radius; x_box++) {
                            int y_other = y + y_box;
                            int x_other = x + x_box;
                            if (y_other < 0 || y_other >= DIM)  { //Checks for Bounds
                                break;
                            } else if (x_other < 0 || x_other >= DIM) {
                                continue;
                            }

                            agent& other = previous.agents[y_other * DIM + x_other];

                            if (other.status == S) { //If neighbor is susceptible
                                float inf = std::rand() % 100;
                                if(inf <= INFECTION_PROBABILITY){
                                    current.agents[y_other * DIM + x_other].status = I;
									inf++;
								}
                            }
                        }

                    }
                    //}
                }
                //}

            }
        }
        // TODO: optimize
        previous = current;
        if(t % 10 == 0) {
            std::cout << std::endl << "---- t: " << t + 1;
            print_board(current);
        }
		sus = sus - rem - inf;
	} // /for t
}

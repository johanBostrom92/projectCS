#include <vector>
#include <omp.h>
#include <cstdlib>
#include <time.h>
#include <iostream>
#include <atomic>
#include <random>
#include <matplot/matplot.h>
#include <numeric>

#define INFECTION_RADIUS 1
#define RECOVERY_RATE 5
#define INFECTION_PROBABILITY 25
#define DIM 10
#define MAX_TIME 10
#define STARTER_AGENTS 1

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
    //bool inQuarantine;
    //bool showingSymptoms;
    int recovery_rate = RECOVERY_RATE;
};

struct board {
    unsigned int dim = DIM;
    std::vector<agent> agents;
};

void print_board(board& b){
    //const std::string status[3] = { "S", "I", "R" };
    for(int i = 0; i < DIM*DIM; i++){
        if(i % DIM == 0){
            std::endl (std::cout);
        }
            std::cout << b.agents[i].status << " ";
    }
}

board generate_board() {
    board b = {
        DIM,
        std::vector<agent>(DIM*DIM)
    };
    // set initial infections
    //for (int i = 0; i < n_patient_zero; i++) {
    //    int pz = std::rand() % dim*dim;
    //    b.agents[pz].status = i;
    //}

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

//          move(uppsala_prev, sthlm_prev, 0);
void swap(board& previous, board& current, int idx) {
    agent swap_agent = current.agents[idx];
    current.agents[idx] = previous.agents[idx];
    previous.agents[idx] = swap_agent;
}

void step(board& previous, board& current, std::mt19937_64 gen) {
#ifdef _WIN32
#pragma omp parallel for
#else
#pragma omp parallel for collapse(2)
#endif
    for (int y = 0; y < DIM; y++) {
        for (int x = 0; x < DIM; x++) {

            //{ Update agent
            agent& self = previous.agents[y * DIM + x];
            agent& currentSelf = current.agents[y * DIM + x];
            if (self.status == I) { //The infected checks for susceptible neighbors within the box.
                currentSelf.recovery_rate--;
                if (currentSelf.recovery_rate == 0) {
                    currentSelf.status = R;
                    //rem++;
                    //inf--;
                    continue;
                }
                //{ Check neighbours
                for (int y_box = -self.infection_radius; y_box <= self.infection_radius; y_box++) {
                    for (int x_box = -self.infection_radius; x_box <= self.infection_radius; x_box++) {
                        int y_other = y + y_box;
                        int x_other = x + x_box;
                        if (y_other < 0 || y_other >= DIM) { //Checks for Bounds
                            break;
                        }
                        else if (x_other < 0 || x_other >= DIM || (x_other == x && y_other == y)) {
                            continue;
                        }

                        agent& other = previous.agents[y_other * DIM + x_other];
                        agent& otherCurr = current.agents[y_other * DIM + x_other];
                        if (other.status == S && otherCurr.status != I) { //If neighbour is susceptible
                            std::uniform_int_distribution<int> dis2(0, 100);
                            int prob = dis2(gen);
                            //int prob = std::rand() % 100;
                            if (prob <= INFECTION_PROBABILITY) {
                                otherCurr.status = I;
                                //inf++;
                            }
                        }
                    }

                }
                //}
            }
            //}

        }
    }
    previous = current;
}


int main() {
    board uppsala_prev= {
        DIM,
        std::vector<agent>(DIM*DIM)
    };
    board sthlm_prev = {
    DIM,
    std::vector<agent>(DIM*DIM)
    };

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<int> dis(0, (DIM*DIM-1));

    int seeded = 0;
    // Func to spawn starter infecitous agents
    while (seeded != STARTER_AGENTS) {
        int pz = dis(gen);
        if (uppsala_prev.agents[pz].status != I) {
            uppsala_prev.agents[pz].status = I;
            std::cout << "pz is: " << pz;
            seeded++;
        }

    }


    board uppsala_curr = uppsala_prev;
    board sthlm_curr = sthlm_prev;
	std::atomic_int sus = DIM*DIM-1;
	std::atomic_int rem = 0;
	std::atomic_int inf = 1;
	std::vector<double> susp = {};
    std::vector<double> infe = {};
    std::vector<double> remo = {};

    for (unsigned int t = 0; t < MAX_TIME; t++)
    { //Loop tracking time
        // TODO: optimize
        step(uppsala_prev, uppsala_curr, gen);
        step(sthlm_prev, sthlm_curr, gen);

        /*if(t % 10 == 0) {
            std::cout << std::endl << "---- t: " << t;
            print_board(current);
        }*/
        //sus = DIM * DIM - rem - inf;
        //susp.push_back(sus);
        //remo.push_back(rem);
        //infe.push_back(inf);

        std::cout <<  std::endl <<"TImestep : " << t << std::endl;
        std::cout << std::endl << " ----Uppsala---- " << std::endl;
        //std::cout << "seeded: " << seeded;
        print_board(uppsala_prev);
        std::cout << std::endl << "----- Sthlm----" << std::endl;
        print_board(sthlm_prev);
        std::cout << std::endl << "<-------------------------------------------------->" << std::endl << std::endl;
	} // /for t
	{
        using namespace matplot;

		std::vector<std::vector<double>> Y = {susp, remo, infe};
        plot(Y);
        title("infected people");
        xlabel("t (days)");
        ylabel("population");
#ifndef _WIN32
        legend({"s", "r", "i"});
#endif

        show();
    }
    return 0;
}

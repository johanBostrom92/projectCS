#include <vector>
#include <omp.h>
#include <cstdlib>
#include <time.h>
#include <iostream>
#include <atomic>
#include <random>
#include <matplot/matplot.h>
#include <numeric>
#include <cmath>
#include <chrono>
#include <tuple>
#include <algorithm>

#define INFECTION_RADIUS 50
#define RECOVERY_RATE 14
#define INFECTION_PROBABILITY 25
#define MAKE_ASYM 50 //The probability that an infected individual is asymptomatic
#define ASYM_INF_PROB 10 //The probability for asymptomatic carriers to infect others
#define DIM 1000
#define MAX_TIME 140
#define STARTER_AGENTS 4
#define QUARANTINE 5
#define Q_FLAG false
#define LAMBDA 2.5
#define ONE_CHANCE false // If TRUE, it chooses only an eligible target to infect. If FALSE, any target can be chosen/tried.

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
        S, A, I, R
};

struct agent {
    agent_status status = S;
    int infection_radius = INFECTION_RADIUS;
    float infection_probability = INFECTION_PROBABILITY;
    int recovery_rate = RECOVERY_RATE;
};

struct board {
    unsigned int dim = DIM;
    std::vector<agent> agents;
};

void print_board(board& b){
    //const std::string status[4] = { "S", "A", "I", "R" };
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

void infect(agent& self, agent& to_infect, std::mt19937_64& gen, std::atomic_int& inf) {
    std::uniform_int_distribution<int> dis2(1, 100);
    int infect_prob;
    if (self.status == I) {
        infect_prob = self.infection_probability;
    }
    else {
        infect_prob = ASYM_INF_PROB;
    }
    int prob = dis2(gen);
    if (prob <= infect_prob) {
        int prob2 = dis2(gen);
        if (prob2 <= MAKE_ASYM) {
            to_infect.status = A;
        }
        else {
            to_infect.status = I;
            inf++;
        }
    }
}

void step(board& previous, board& current, std::mt19937_64& gen, std::atomic_int& rem, std::atomic_int& inf, int t) {
#ifdef _WIN32
#pragma omp parallel for
#else
#pragma omp parallel for collapse(2)
#endif
    for (int y = 0; y < DIM; y++) {
        for (int x = 0; x < DIM; x++) {
            agent& self = previous.agents[y * DIM + x];
            agent& currentSelf = current.agents[y * DIM + x];
            if (self.status == I || self.status == A) { //The infected checks for susceptible neighbors within the box.
                currentSelf.recovery_rate--;
                if (currentSelf.recovery_rate == 0) {
                    if (self.status == I) {//Infected only decreased if carrier was not asymptomatic
                        inf--;
                    }
                    currentSelf.status = R;
                    rem++;
                    continue;
                }
                int rad = self.infection_radius;
                if (t >= QUARANTINE && Q_FLAG) {
                    rad = static_cast<double>(self.infection_radius) * exp(((-(t - static_cast<double>(QUARANTINE))) / LAMBDA));
                    //std::cout << "current radious is " << rad << std::endl;

                }
                if (rad <= 0) {
                    continue;
                }

                if (ONE_CHANCE) {
                    //std::random_device rd;
                    //std::mt19937_64 gen(rd());
                    //std::uniform_int_distribution<int> dis(-rad, rad);

                    std::vector<std::tuple<int, int>> checked;
                    for (int y_box = -rad; y_box <= rad; y_box++) {
                        for (int x_box = -rad; x_box <= rad; x_box++) {
                            int y_other = y + y_box;
                            int x_other = x + x_box;
                            if (y_other < 0 || y_other >= DIM) { //Checks for Bounds
                                break;
                            }
                            else if (x_other < 0 || x_other >= DIM || (x_other == x && y_other == y)) {
                                continue;
                            }
                            agent& other = previous.agents[y_other * DIM + x_other];
                            if (other.status == S) {
                                checked.push_back(std::make_tuple(x_other, y_other));
                            }
                        }
                    }

                Repeat:
                    int size = checked.size();
                    std::uniform_int_distribution<int> dis3(0, size - 1);
                    int rand_S = dis3(gen);
                    if (size <= 0) {
                        continue;
                    }
                    int y_other = std::get<1>(checked.at(rand_S));
                    int x_other = std::get<0>(checked.at(rand_S));
                    agent& other = previous.agents[y_other * DIM + x_other];
                    agent& otherCurr = current.agents[y_other * DIM + x_other];

                    if (other.status == S && otherCurr.status != I && otherCurr.status != A) { //If neighbour is susceptible
                        //int prob = std::rand() % 100;
                        infect(self, otherCurr, gen, inf);
                    } 

                    else {

                        checked.erase(checked.begin() + rand_S);

                        if (checked.size() == 0) { continue; }

                        goto Repeat;

                    }

                } 
                else { // Doesnt care if target can be infected or not.
                    //{ Update agent
                    std::uniform_int_distribution<int> dis(-rad, rad);
                Repeating:
                    int x_rand = dis(gen);
                    int y_rand = dis(gen);
                    int y_other = y_rand + y;
                    int x_other = x_rand + x;
                            if (y_other < 0 || y_other >= DIM || x_other < 0 || x_other >= DIM || (x_other == x && y_other == y)) { //Checks for Bounds
                                goto Repeating; //Repeat failures
                                //continue; //Skip failures
                            }
     
                    agent& other = previous.agents[y_other * DIM + x_other];
                    agent& otherCurr = current.agents[y_other * DIM + x_other];
                    if (other.status == S && otherCurr.status != I && otherCurr.status != A) { //If neighbour is susceptible
                        infect(self, otherCurr, gen, inf);
                    }

                    
                        
                    
                }
            }
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
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
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
	std::atomic_int sus = DIM*DIM-STARTER_AGENTS;
	std::atomic_int rem = 0;
	std::atomic_int inf = STARTER_AGENTS;
	std::vector<double> susp = {};
    std::vector<double> infe = {};
    std::vector<double> remo = {};

    for (unsigned int t = 0; t < MAX_TIME; t++)
    { //Loop tracking time
        // TODO: optimize
        step(uppsala_prev, uppsala_curr, gen,rem,inf,t);
        //step(sthlm_prev, sthlm_curr, gen,rem,inf);

        /*if(t % 10 == 0) {
            std::cout << std::endl << "---- t: " << t;
            print_board(current);
        }*/
        //sus = DIM * DIM - rem - inf;
        //susp.push_back(sus);
        //remo.push_back(rem);
        //infe.push_back(inf);

        //std::cout <<  std::endl <<"Timestep : " << t << std::endl;
        //std::cout << std::endl << " ----Uppsala---- " << std::endl;
        //std::cout << "seeded: " << seeded;
        //print_board(uppsala_prev);
        //std::cout << std::endl << "----- Sthlm----" << std::endl;
        //print_board(sthlm_prev);
        //std::cout << std::endl << "<-------------------------------------------------->" << std::endl << std::endl;
        sus = DIM * DIM - rem - inf;
       // std::cout << "printing sus: " << sus << " Printing Rem: " << rem << " Printing inf: " << inf << std::endl;
        susp.push_back(sus);
        remo.push_back(rem);
        infe.push_back(inf);
    } // /for t
	{
        using namespace matplot;

		std::vector<std::vector<double>> y = {susp, remo, infe};
        plot(y);
        title("infected people");
        xlabel("t (days)");
        ylabel("population");
/*#ifndef _WIN32
        legend({"s", "r", "i"});
#endif
*/
        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << std::endl << "It took  " << time_span.count() << " seconds." << std::endl;
        show(); //TODO - uncomment
    }
    return 0;
}

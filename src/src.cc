#include "types.hh"
#include "parameters.hh"
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
#include <cassert>



void print_board(board& b, std::string name, int t) {
    std::cout << std::endl << "Susceptible = 0 -- Asymptomatic = 1 -- Infected = 2 -- Vaccinated = 3 -- Recovered = 4 " << std::endl;
    std::cout << std::endl << "Day: " << t << std::endl;
    std::cout << std::endl << " ----------" << name << "----------" << std::endl;
    for (int i = 0; i < DIM * DIM; i++) {
        if (i % DIM == 0) {
            std::endl(std::cout);
        }
            std::cout << b.agents[i].status << " ";

    }
    std::cout << std::endl << std::endl << "Susceptible: " << b.sus << std::endl << "Recovered: " << b.rem << std::endl << "Infected: " << b.inf << std::endl << "Asymptomatic: " << b.asymp << std::endl << "Vaccinated: " << b.vacc;
    std::cout << std::endl << "-----------------------------" << std::endl << std::endl;
}


void swap(board& previous, board& current, int idx) {
    agent swap_agent = current.agents[idx];
    current.agents[idx] = previous.agents[idx];
    previous.agents[idx] = swap_agent;
}

void vaccinate(agent& agent) {
    agent.vaccination_progress = true;
}

void vaccinate(board& b, int idx) {
    b.agents[idx].vaccination_progress = true;
}

void infect(agent& self, agent& to_infect, std::mt19937_64& gen, board& b) {
    std::uniform_int_distribution<int> dis2(0, 99);
    int infect_prob;
    if (self.status == I) {
        infect_prob = self.infection_probability;
    }
    else {
        infect_prob = ASYM_INF_PROB;
    }
    int prob = dis2(gen);
    if (prob < infect_prob) {
        int prob2 = dis2(gen);
        if (prob2 < MAKE_ASYM) {
            to_infect.status = A;
            b.asymp++;
        }
        else {
            to_infect.status = I;
            b.inf++;
        }
    }
}

void step(board& previous, board& current, std::mt19937_64& gen, int t) {
#ifdef _WIN32
#pragma omp parallel for
#else
#pragma omp parallel for collapse(2)
#endif
    for (int y = 0; y < DIM; y++) {
        for (int x = 0; x < DIM; x++) {
            agent& self = previous.agents[y * DIM + x];
            agent& currentSelf = current.agents[y * DIM + x];
            std::uniform_int_distribution<int> vacc_dis(0, 99);
            if (self.vaccination_progress) {
                currentSelf.vaccination_rate--;
                if (currentSelf.vaccination_rate == 0) {
                    int vacc_check = vacc_dis(gen);
                    if (vacc_check < VACCINATION_EFFICACY) {
                        //Congrats! You're vaccinated!
                        // TODO: handle A and I specially?
                        switch (currentSelf.status) {
                            case S:
                                current.sus--;
                                break;
                            case A:
                                current.asymp--;
                                break;
                            case I:
                                current.inf--;
                                break;
                            case V:
                                current.vacc--;
                                break;
                            case R:
                                current.rem--;
                                break;
                            default:
                                std::cout << "WARN: Unknown agent status: " << currentSelf.status << std::endl;
                                break;
                        }
                        currentSelf.status = V;
                        current.vacc++;
                        continue;
                    }
                    else {
                        //Better luck next time!
                        //currentSelf.vaccination_progress = false;
                    }
                }
            }
            if (self.status == I || self.status == A) { //The infected checks for susceptible neighbors within the box.
                currentSelf.recovery_rate--;
                if (currentSelf.recovery_rate == 0) {
                    if (self.status == I) {//Infected only decreased if carrier was not asymptomatic
                        current.inf--;
                    } else if (self.status == A) {
                        current.asymp--;
                    }
                    currentSelf.status = R;
                    current.rem++;
                    continue;
                }
                int rad = self.infection_radius;
                if (t >= QUARANTINE_START && ENABLE_QUARANTINE) {
                    rad = static_cast<double>(self.infection_radius) * exp(((-(t - static_cast<double>(QUARANTINE_START))) / LAMBDA));

                }
                if (rad <= 0) {
                    continue;
                }

                if (ONLY_ELIGIBLE) {
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

                    if (other.status == S && otherCurr.status == S) { //If neighbour is susceptible
                        infect(self, otherCurr, gen, current);
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
                    if (other.status == S && otherCurr.status == S) { //If neighbour is susceptible
                        infect(self, otherCurr, gen, current);
                    }




                }
            }
        }
    }
    current.sus = current.dim * current.dim - current.rem - current.inf - current.asymp - current.vacc;
    previous = current;
}

#define STATUS_TO_WEIGHT(s) ((s) == I ? 1 : 0)

/**
 * Updates all vaccination weights for a board, to determine who should be more
 * or less likely to receive a vaccine. This is based on the number of infected agents within
 * INFECTION_RADIUS of each agent.
 * TODO: make it possible to customize the strategy for how the weights are chosen
 */
void update_vaccination_weights(board& b) {
    b.vaccination_weight_sum = 0;
    std::vector<unsigned int> tmp(b.vaccination_weights.size());

    // Essentially, this is a box blur, which we can separate into one pass for each axis
    // (we first blur on the x axis and store it to tmp, then blur tmp on the y axis and store it back to the board)

    #pragma omp parallel for
    for (int y = 0; y < b.dim; y++) {
        unsigned int sum = 0;
        for (int x_box = 0; x_box < INFECTION_RADIUS; x_box++) {
            sum += STATUS_TO_WEIGHT(b.agents[y*b.dim + x_box].status);
        }
        for (int x = 0; x < b.dim; x++) {
            if (x + INFECTION_RADIUS < b.dim) {
                sum += STATUS_TO_WEIGHT(b.agents[y*b.dim + x + INFECTION_RADIUS].status);
            }
            tmp[y*b.dim + x] = sum;
            if (x - INFECTION_RADIUS >= 0) {
                sum -= STATUS_TO_WEIGHT(b.agents[int(y*b.dim) + x - INFECTION_RADIUS].status);
            }
        }
    }
    #pragma omp parallel for
    for (int x = 0; x < b.dim; x++) {
        unsigned int sum = 0;
        for (int y_box = 0; y_box < INFECTION_RADIUS; y_box++) {
            sum += tmp[y_box*b.dim + x];
        }
        for (int y = 0; y < b.dim; y++) {
            if (y + INFECTION_RADIUS < b.dim) {
                sum += tmp[(y + INFECTION_RADIUS)*b.dim + x];
            }
            int idx = y*b.dim + x;
            int weight = 1 + sum;
            if (b.agents[idx].status == V || b.agents[idx].vaccination_progress) {
                weight = 0;
            }
            b.vaccination_weights[idx] = weight;
            b.vaccination_weight_sum.fetch_add(weight, std::memory_order_acq_rel);
            if (y - INFECTION_RADIUS >= 0) {
                sum -= tmp[(y - INFECTION_RADIUS)*int(b.dim) + x];
            }
        }
    }
}

/**
 * Vaccinates a number of agents in the board, based on the board's vaccination weights
 * @param board The board to perform vaccinations in
 * @param n_vaccinations The number of vaccinations to perform
 */
void vaccinate(board& b, unsigned int n_vaccinations, std::mt19937_64& gen) {
    // Make sure we don't perform more vaccinations than there are unvaccinated agents
    n_vaccinations = std::min((size_t) n_vaccinations, b.agents.size() - b.vaccinations_started);
    if (n_vaccinations == 0 || b.vaccination_weight_sum == 0) {
        return;
    }

    // For each agent, compute the combined weight of it and all agents with a lower index
    // This lets us later use lower_bound to do a binary search for the agent to vaccinate
    std::vector<unsigned int> cumulative_weights(b.agents.size());
    unsigned int weight_sum = 0;
    for(int i = 0; i < b.agents.size(); i++) {
        weight_sum += b.vaccination_weights[i];
        cumulative_weights[i] = weight_sum;
    }
    assert(cumulative_weights[b.agents.size() - 1] = b.vaccination_weight_sum);

    std::uniform_int_distribution<int> dis(1, (int) b.vaccination_weight_sum); // Generates (cumulative) weight values determining who gets vaccinated
    std::vector<std::atomic_flag> vaccinated(b.agents.size());  // Stores which agents have been vaccinated this round
    unsigned int collisions = 0;

    #pragma omp parallel for
    for (int i = 0; i < n_vaccinations; i++) {
        int n = dis(gen);
        // Finds the first agent with cumulative weight greater than n
        auto it = std::lower_bound(cumulative_weights.begin(), cumulative_weights.end(), n);
        unsigned int agent = it - cumulative_weights.begin();
        // Atomically check if this agent has been vaccinated this round
        if (!vaccinated[agent].test_and_set()) {
            vaccinate(b.agents[agent]);
            b.vaccinations_started.fetch_add(1, std::memory_order_acq_rel);
            b.vaccination_weight_sum.fetch_sub(b.vaccination_weights[agent], std::memory_order_acq_rel);
            b.vaccination_weights[agent] = 0;
        } else {
            // If it's been vaccinated, we will need to do it over
            collisions++;
        }
    }
    if (collisions > 0) {
        // Redo the vaccinations that collided. It might be better to use a collision-free algorithm here.
        vaccinate(b, collisions, gen);
    }
}

int main() {
    board uppsala_prev(DIM, STARTER_AGENTS, AGENT_TYPES);
    board sthlm_prev(DIM, STARTER_AGENTS, AGENT_TYPES);
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    std::random_device rd;
    std::mt19937_64 gen(rd());

    board uppsala_curr = uppsala_prev;
    board sthlm_curr = sthlm_prev;

    std::vector<double> uppsala_susp = {};
    std::vector<double> uppsala_infe = {};
    std::vector<double> uppsala_remo = {};
    std::vector<double> uppsala_asymp = {};
    std::vector<double> uppsala_vacc = {};
    std::vector<double> sthlm_susp = {};
    std::vector<double> sthlm_infe = {};
    std::vector<double> sthlm_remo = {};
    std::vector<double> sthlm_asymp = {};
    std::vector<double> sthlm_vacc = {};

    for (unsigned int t = 0; t < MAX_TIME; t++)
    { //Loop tracking time

        step(uppsala_prev, uppsala_curr, gen, t);
        if (t > VACCINATION_START) {
            update_vaccination_weights(uppsala_curr);
            vaccinate(uppsala_curr, VACCINATIONS_PER_DAY, gen);

        }
        //step(sthlm_prev, sthlm_curr, gen, t);
        // print_board(uppsala_prev, "Uppsala", t);
        uppsala_susp.push_back(uppsala_curr.sus);
        uppsala_remo.push_back(uppsala_curr.inf);
        uppsala_infe.push_back(uppsala_curr.rem);
        uppsala_asymp.push_back(uppsala_curr.asymp);
        uppsala_vacc.push_back(uppsala_curr.vacc);
        sthlm_susp.push_back(sthlm_curr.sus);
        sthlm_remo.push_back(sthlm_curr.inf);
        sthlm_infe.push_back(sthlm_curr.rem);
        sthlm_asymp.push_back(sthlm_curr.asymp);
        sthlm_vacc.push_back(sthlm_curr.vacc);

    } // /for t
    {
        using namespace matplot;

        std::vector<std::vector<std::vector<double>>> plot_data{
            { uppsala_susp, uppsala_remo, uppsala_infe, uppsala_asymp, uppsala_vacc },
            { sthlm_susp, sthlm_remo, sthlm_infe }
        };

        std::vector<std::string> comm_names = { "Uppsala", "Stockholm" };

        for (int i = 0; i < plot_data.size(); i++) {
            auto f = figure();
            auto ax = f->current_axes();
            plot(ax, plot_data[i]);
            title(ax, comm_names[i]);
            xlabel(ax, "t (days)");
            ylabel(ax, "population");
#ifndef _WIN32
            legend(ax, {"s", "r", "i"});
#endif
        }

        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << std::endl << "It took  " << time_span.count() << " seconds." << std::endl;
    }
    std::cout << "Press Enter to exit..." << std::endl;
    std::cin.get();
    return 0;
}

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

void print_board(board& b){
    //const std::string status[3] = { "S", "I", "R" };
    for(int i = 0; i < DIM*DIM; i++){
        if(i % DIM == 0){
            std::endl (std::cout);
        }
        std::cout << b.agents[i].status << " ";
    }
}

//          move(uppsala_prev, sthlm_prev, 0);
void swap(board& previous, board& current, int idx) {
    agent swap_agent = current.agents[idx];
    current.agents[idx] = previous.agents[idx];
    previous.agents[idx] = swap_agent;
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
            if (self.status == I) { //The infected checks for susceptible neighbors within the box.
                currentSelf.recovery_rate--;
                if (currentSelf.recovery_rate == 0) {
                    currentSelf.status = R;
                    current.rem++;
                    current.inf--;
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

                    if (other.status == S && otherCurr.status != I) { //If neighbour is susceptible
                        //int prob = std::rand() % 100;
                        std::uniform_int_distribution<int> dis2(0, 99);
                        int prob = dis2(gen);
                        if (prob < INFECTION_PROBABILITY) {
                            otherCurr.status = I;
                            current.inf++;
                        }

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
                    if (other.status == S && otherCurr.status != I) { //If neighbour is susceptible

                        std::uniform_int_distribution<int> dis2(0, 99);
                        int prob = dis2(gen);
                        if (prob < INFECTION_PROBABILITY) {
                            // std::cout << "Infected" << std::endl;
                            otherCurr.status = I;
                            current.inf++;
                        }

                    }




                }
            }
        }
    }
    current.sus = current.dim * current.dim - current.rem - current.inf;
    previous = current;
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
    std::vector<double> sthlm_susp = {};
    std::vector<double> sthlm_infe = {};
    std::vector<double> sthlm_remo = {};

    for (unsigned int t = 0; t < MAX_TIME; t++)
    { //Loop tracking time
        // TODO: optimize
        step(uppsala_prev, uppsala_curr, gen, t);
        //step(sthlm_prev, sthlm_curr, gen,rem,inf);

        /*if(t % 10 == 0) {
            std::cout << std::endl << "---- t: " << t;
            print_board(current);
        }*/
        //sus = DIM * DIM - rem - inf;
        //susp.push_back(sus);
        //remo.push_back(rem);
        //infe.push_back(inf);

        //std::cout <<  std::endl <<"TImestep : " << t << std::endl;
        //std::cout << std::endl << " ----Uppsala---- " << std::endl;
        //std::cout << "seeded: " << seeded;
        //print_board(uppsala_prev);
        //std::cout << std::endl << "----- Sthlm----" << std::endl;
        //print_board(sthlm_prev);
        //std::cout << std::endl << "<-------------------------------------------------->" << std::endl << std::endl;
       // std::cout << "printing sus: " << sus << " Printing Rem: " << rem << " Printing inf: " << inf << std::endl;
        uppsala_susp.push_back(uppsala_curr.sus);
        uppsala_remo.push_back(uppsala_curr.inf);
        uppsala_infe.push_back(uppsala_curr.rem);
        sthlm_susp.push_back(sthlm_curr.sus);
        sthlm_remo.push_back(sthlm_curr.inf);
        sthlm_infe.push_back(sthlm_curr.rem);
    } // /for t
    {
        using namespace matplot;

        std::vector<std::vector<std::vector<double>>> plot_data{
            { uppsala_susp, uppsala_remo, uppsala_infe },
            { sthlm_susp, sthlm_remo, sthlm_infe }
        };

        for (int i = 0; i < plot_data.size(); i++) {
            auto f = figure();
            auto ax = f->current_axes();
            plot(ax, plot_data[i]);
            title(ax, "Community " + std::to_string(i+1));
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

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
#include <math.h>


void print_board(board& b, int t) {
    if(t== 0)
    std::cout << std::endl << "Susceptible = 0 -- Asymptomatic = 1 -- Infected = 2 -- Vaccinated = 3 -- Recovered = 4 " << std::endl;
    std::cout << std::endl << "Day: " << t << std::endl;
    std::cout << std::endl << " ----------" << b.name << "----------" << std::endl;
    for (int i = 0; i < b.dim * b.dim; i++) {
        if (i % b.dim == 0) {
            std::endl(std::cout);
        }
            std::cout << b.agents[i].status << " ";

    }
    std::cout << std::endl << std::endl << "Susceptible: " << b.sus << std::endl << "Recovered: " << b.rem << std::endl << "Infected: " << b.inf << std::endl << "Asymptomatic: " << b.asymp << std::endl << "Vaccinated: " << b.vacc;
    std::cout << std::endl << "-----------------------------" << std::endl << std::endl;
}


void swap(board& fromBoard, board& toBoard, int fromIndex, int toIndex) {
    agent swap_agent = toBoard.agents[toIndex];
    toBoard.agents[toIndex] = fromBoard.agents[fromIndex];
    fromBoard.agents[fromIndex] = swap_agent;
    std::cout << std::endl << fromBoard.name << ":" << fromIndex << " <==> " << toBoard.name << ":" << toIndex << std::endl;
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
std::vector<std::atomic_flag> infected(previous.agents.size());
#ifdef _WIN32
#pragma omp parallel for
#else
#pragma omp parallel for collapse(2)
#endif
    for (int y = 0; y < current.dim; y++) {
        for (int x = 0; x < current.dim; x++) {

            agent& self = previous.agents[y * current.dim + x];
            agent& currentSelf = current.agents[y * current.dim + x];
            std::uniform_int_distribution<int> vacc_dis(0, 99);
            if (self.vaccination_progress) {
                currentSelf.vaccination_rate--;
                if (currentSelf.vaccination_rate == 0) {
                    int vacc_check = vacc_dis(gen);
                    if (vacc_check < VACCINATION_EFFICACY) {
                        //Congrats! You're vaccinated!
                        currentSelf.status = V;
                        current.vacc++;
                        continue;
                    }
                    else {
                        //Better luck next time!
                    }
                }
            }
            if (self.status == I || self.status == A) { //The infected checks for susceptible neighbors within the box.
                currentSelf.recovery_rate--;
                if (currentSelf.recovery_rate == 0) {
                    if (self.status == I) {//Infected only decreased if carrier was not asymptomatic
                        current.inf--;
                    }
                    else if(self.status == A) {
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
                            if (y_other < 0 || y_other >= current.dim) { //Checks for Bounds
                                break;
                            }
                            else if (x_other < 0 || x_other >= current.dim || (x_other == x && y_other == y)) {
                                continue;
                            }
                            agent& other = previous.agents[y_other * current.dim + x_other];
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
                    agent& other = previous.agents[y_other * current.dim + x_other];
                    agent& otherCurr = current.agents[y_other * current.dim + x_other];

                    if (other.status == S && !infected[y_other * current.dim + x_other].test_and_set()) { //If neighbour is susceptible
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
                            if (y_other < 0 || y_other >= current.dim || x_other < 0 || x_other >= current.dim || (x_other == x && y_other == y)) { //Checks for Bounds
                        goto Repeating; //Repeat failures
                                //continue; //Skip failures
                    }

                    agent& other = previous.agents[y_other * current.dim + x_other];
                    agent& otherCurr = current.agents[y_other * current.dim + x_other];
                    if (other.status == S && !infected[y_other * current.dim + x_other].test_and_set()) { //If neighbour is susceptible
                        infect(self, otherCurr, gen, current);
                    }




                }
            }
        }
    }
    current.sus = current.dim * current.dim - current.rem - current.inf - current.asymp - current.vacc;
    previous = current;
}

int whereToMove(std::vector<double> items, std::mt19937_64& gen) {
    std::uniform_real_distribution move_dist(0.0, 1.0);
    double cumulative = 1;
    double rand = move_dist(gen);
    for (int i = 0; i < items.size(); i++) {
        cumulative -= items[i];
        if (rand >= cumulative) {
            return i;
        }
    }

}






void updateBoard(board& from, board& to, agent agentFrom, agent agentTo) {
    int agentFromA, agentFromS, agentFromV, agentFromI, agentFromR;
    agentFromA = agentFromS = agentFromV = agentFromI = agentFromR = 0;
    int agentToA, agentToS, agentToV, agentToI, agentToR;
    agentToA = agentToS = agentToV = agentToI = agentToR = 0;
    switch (agentFrom.status) {
        case S:
            agentFromS++;
            break;
        case A:
            agentFromA++;
            break;
        case I:
            agentFromI++;
            break;
        case V:
            agentFromV++;
            break;
        case R:
            agentFromR++;
            break;
        default:
            //Should not execute
            std::cout << "Unknown type!";
   }

    switch (agentTo.status) {
        case S:
            agentToS++;
            break;
        case A:
            agentToA++;
            break;
        case I:
            agentToI++;
            break;
        case V:
            agentToV++;
            break;
        case R:
            agentToR++;
            break;
        default:
            //Should not execute
            std::cout << "Unknown type!";
    }

    if (&from != &to) {
        //from.sus -= agentFromS; Not needed
        from.asymp -= agentFromA;
        from.vacc -= agentFromV;
        from.inf -= agentFromI;
        from.rem -= agentFromR;

        //to.sus -= agentToS;
        to.asymp -= agentToA;
        to.vacc -= agentToV;
        to.inf -= agentToI;
        to.rem -= agentToR;

        //from.sus += agentToS;
        from.asymp += agentToA;
        from.vacc += agentToV;
        from.inf += agentToI;
        from.rem += agentToR;

        //to.sus += agentFromS;
        to.asymp += agentFromA;
        to.vacc += agentFromV;
        to.inf += agentFromI;
        to.rem += agentFromR;
    }

}

int main() {

    std::vector<std::string> comm_names = { "Uppsala", "Stockholm" };
    std::vector<double> weight = {             0.50,       0.50 };
    std::vector<int> dimensions = {};
    int population = 10;
    std::vector<board> prev_board = {};
    std::vector<board> curr_board = {};

    for (int i = 0; i < comm_names.size(); i++)
    {
        int calc_dim = ceil(sqrt(weight[i] * population));
        dimensions.push_back(calc_dim);

        board new_board(calc_dim, STARTER_AGENTS, AGENT_TYPES, weight[i], comm_names[i]);


        prev_board.push_back(new_board);
        board curr = new_board;
        curr_board.push_back(curr);
    }

    //std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::vector<double> uppsala_susp = {};
    std::vector<double> uppsala_infe = {};
    std::vector<double> uppsala_remo = {};
    std::vector<double> uppsala_asymp = {};
    std::vector<double> uppsala_vacc = {};
   /* std::vector<double> sthlm_susp = {};
    std::vector<double> sthlm_infe = {};
    std::vector<double> sthlm_remo = {};
    std::vector<double> sthlm_asymp = {};
    std::vector<double> sthlm_vacc = {};*/

    std::uniform_int_distribution<int> board_dis(0, comm_names.size()-1);

    for (unsigned int t = 0; t < MAX_TIME; t++)
    { //Loop tracking time
        // TODO: optimize

        for (int i = 0; i < 5; i++) {
            int fromBoardIdx = board_dis(gen);
            board& fromBoard = curr_board[fromBoardIdx];
            int fromDim = fromBoard.dim;

            std::uniform_int_distribution<int> fromAgent_dis(0, (fromDim * fromDim - 1));
            int agentfromBoardIdx = fromAgent_dis(gen);
            agent agentFromBoard = fromBoard.agents[agentfromBoardIdx];
            int toBoardIdx = whereToMove(weight, gen);
            board& toBoard = curr_board[toBoardIdx];
            int targetDim = toBoard.dim;

            std::uniform_int_distribution<int> targetAgent_dis(0, (targetDim * targetDim - 1));
            int agentToBoardIdx = targetAgent_dis(gen);
            agent agentToBoard = toBoard.agents[agentToBoardIdx];
            swap(fromBoard, toBoard, agentfromBoardIdx, agentToBoardIdx);
            updateBoard(fromBoard, toBoard, agentFromBoard, agentToBoard);
        }

        for (int i = 0; i < comm_names.size(); i++)
        {
            prev_board[i] = curr_board[i];
            step(prev_board[i], curr_board[i], gen, t);
            print_board(prev_board[i], t);
        }


        /*uppsala_susp.push_back(curr_board[0].sus);
        uppsala_remo.push_back(curr_board[0].inf);
        uppsala_infe.push_back(curr_board[0].rem);
        uppsala_asymp.push_back(curr_board[0].asymp);
        uppsala_vacc.push_back(curr_board[0].vacc);
       
        */
        /* sthlm_susp.push_back(sthlm_curr.sus);
        sthlm_remo.push_back(sthlm_curr.inf);
        sthlm_infe.push_back(sthlm_curr.rem);
        sthlm_asymp.push_back(sthlm_curr.asymp);
        sthlm_vacc.push_back(sthlm_curr.vacc);*/

    } // /for t
    {
//        using namespace matplot;
//
//        std::vector<std::vector<std::vector<double>>> plot_data{
//            { uppsala_susp, uppsala_remo, uppsala_infe },
//            //{ sthlm_susp, sthlm_remo, sthlm_infe }
//        };
//
//        for (int i = 0; i < plot_data.size(); i++) {
//            auto f = figure();
//            auto ax = f->current_axes();
//            plot(ax, plot_data[i]);
//            title(ax, comm_names[i]);
//            xlabel(ax, "t (days)");
//            ylabel(ax, "population");
//#ifndef _WIN32
//            legend(ax, {"s", "r", "i"});
//#endif
//        }

       /* std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << std::endl << "It took  " << time_span.count() << " seconds." << std::endl;*/
    }
    std::cout << "Press Enter to exit..." << std::endl;
    std::cin.get();
    return 0;
}

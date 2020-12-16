#define _USE_MATH_DEFINES
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
#include <visualization.hh>
#include <cassert>
#include <iomanip> 




// Holds one random number generator for each thread
static std::vector<std::mt19937_64> generators;

/**
 * Prints the current board
 * @param b The board to be printed
 * @param name The name of the board/community
 * @param t The current timestep
 */
void print_board(board& b, std::string name, int t) {
    std::cout << std::endl << "Susceptible = 0 -- Asymptomatic = 1 -- Infected = 2 -- Vaccinated = 3 -- Recovered = 4 " << std::endl;
    std::cout << std::endl << "Day: " << t << std::endl;
    std::cout << std::endl << " ----------" << b.name << "----------" << std::endl;
    for (int i = 0; i < b.dim * b.dim; i++) {
        if (i % b.dim == 0) {
            std::endl(std::cout);
        }
        std::cout << b.agents[i].status << " ";

    }
    std::cout << std::endl << std::endl << "Susceptible: " << b.status_counts[S] << std::endl << "Recovered: " << b.status_counts[R] << std::endl << "Infected: " << b.status_counts[I] << std::endl << "Asymptomatic: " << b.status_counts[A] << std::endl << "Vaccinated: " << b.status_counts[V];
    std::cout << std::endl << "-----------------------------" << std::endl << std::endl;
}

//          move(uppsala_prev, sthlm_prev, 0);


/**
 * Swaps the position of an agent from one board to another
 * @param previous Board to swap from
 * @param current Board to swap to
 * @param idx Index of the agent to swap
 */
void swap(board& fromBoard, board& toBoard, int fromIndex, int toIndex) {
    agent swap_agent = toBoard.agents[toIndex];
    toBoard.agents[toIndex] = fromBoard.agents[fromIndex];
    fromBoard.agents[fromIndex] = swap_agent;
    std::cout << std::endl << fromBoard.name << ":" << fromIndex << " <==> " << toBoard.name << ":" << toIndex << std::endl;
}


/**
 * Start the vaccination process of an agent
 * @param agent The agent to vaccinate
 */
void vaccinate(agent& agent) {
    agent.vaccination_progress = true;
}

/**
 * One agent which infects another agent
 * @param self The agent which will infect
 * @param to_infect The agent to infect
 * @param gen The current random number generator
 * @param b The board b in which both ag ents reside
 */
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
        }
        else {
            to_infect.status = I;
        }
        b.status_counts[to_infect.status]++;
        b.total_infections++;
    }
}

/**
 * Take one step in the simulation
 * @param previous The previous state of the board
 * @param current The current state of the board
 * @param t The current timestep
 */
void step(board& previous, board& current, int t) {
    std::vector<std::atomic_flag> infected(previous.agents.size());
#ifdef _WIN32
#pragma omp parallel for
#else
#pragma omp parallel for collapse(2)
#endif
    for (int y = 0; y < current.dim; y++) {
        for (int x = 0; x < current.dim; x++) {
            std::mt19937_64& gen = generators[omp_get_thread_num()];

            agent& self = previous.agents[y * current.dim + x];
            agent& currentSelf = current.agents[y * current.dim + x];
            std::uniform_int_distribution<int> vacc_dis(0, 99);
            if (self.vaccination_progress) { //Check first if agent is vaccinated
                currentSelf.vaccination_rate--;
                if (currentSelf.vaccination_rate == 0) {
                    int vacc_check = vacc_dis(gen);
                    if (vacc_check < VACCINATION_EFFICACY) {
                        //Congrats! You're vaccinated!
                        // TODO: handle A and I specially? i.e., does the vaccinate 'beat' the infection if it's already there?
                        current.status_counts[currentSelf.status]--;
                        currentSelf.status = V;
                        current.status_counts[V]++;
                        continue;
                    }
                }
            }
            if (self.status == I || self.status == A) { //The infected checks for susceptible neighbors within the box.
                currentSelf.recovery_rate--;
                if (currentSelf.recovery_rate == 0) {
                    current.status_counts[self.status]--;
                    currentSelf.status = R;
                    current.status_counts[R]++;
                    continue;
                }
                int rad = self.infection_radius;
                if (t >= QUARANTINE_START && ENABLE_QUARANTINE) {
                    rad = static_cast<double>(self.infection_radius) * exp(((-(t - static_cast<double>(QUARANTINE_START))) / LAMBDA));

                }
                if (rad <= 0) {
                    continue;
                }

                if constexpr (ONLY_ELIGIBLE) {
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
    current.status_counts[S] = current.dim * current.dim - current.status_counts[R] - current.status_counts[I] - current.status_counts[A] - current.status_counts[V];
    previous = current;
}

#define IS_INFECTED(s) ((s) == I ? 1 : 0)

/**
 * Updates all vaccination weights for a board, to determine who should be more
 * or less likely to receive a vaccine. This is based on the vaccination strategy
 * denoted by the VACC_STRAT parameter.
 */
void update_vaccination_weights(board& b) {
    b.vaccination_weight_sum = 0;
    std::vector<unsigned int> tmp(b.vaccination_weights.size());

    if constexpr (VACC_STRAT == vaccination_strategy::UNIFORM) {
        std::fill(b.vaccination_weights.begin(), b.vaccination_weights.end(), 1);
        b.vaccination_weight_sum = b.vaccination_weights.size();
    }
    else if constexpr (VACC_STRAT == vaccination_strategy::HIGH_DENSITY || VACC_STRAT == vaccination_strategy::LOW_DENSITY) {
        // Essentially, this is a box blur, which we can separate into one pass for each axis
        // (we first blur on the x axis and store it to tmp, then blur tmp on the y axis and store it back to the board)
#pragma omp parallel for
        for (int y = 0; y < b.dim; y++) {
            unsigned int sum = 0;
            for (int x_box = 0; x_box < std::min(INFECTION_RADIUS, (int)b.dim - 1); x_box++) {
                sum += IS_INFECTED(b.agents[y * b.dim + x_box].status);
            }
            for (int x = 0; x < b.dim; x++) {
                if (x + INFECTION_RADIUS < b.dim) {
                    sum += IS_INFECTED(b.agents[y * b.dim + x + INFECTION_RADIUS].status);
                }
                tmp[y * b.dim + x] = sum;
                if (x - INFECTION_RADIUS >= 0) {
                    sum -= IS_INFECTED(b.agents[int(y * b.dim) + x - INFECTION_RADIUS].status);
                }
            }
        }
        unsigned int maxSum = 0;
        for (int x = 0; x < b.dim; x++) {
            unsigned int sum = 0;
            for (int y_box = 0; y_box < std::min(INFECTION_RADIUS, (int)b.dim - 1); y_box++) {
                sum += tmp[y_box * b.dim + x];
            }
            for (int y = 0; y < b.dim; y++) {
                if (y + INFECTION_RADIUS < b.dim) {
                    sum += tmp[(y + INFECTION_RADIUS) * b.dim + x];
                }
                int idx = y * b.dim + x;
                maxSum = std::max(sum, maxSum);
                b.vaccination_weights[idx] = sum;
                if (y - INFECTION_RADIUS >= 0) {
                    sum -= tmp[(y - INFECTION_RADIUS) * int(b.dim) + x];
                }
            }
        }
        for (int i = 0; i < b.agents.size(); i++) {
            if constexpr (VACC_STRAT == vaccination_strategy::HIGH_DENSITY) {
                b.vaccination_weights[i] += 1;
            }
            else {
                b.vaccination_weights[i] = 1 + maxSum - b.vaccination_weights[i];
            }
            b.vaccination_weight_sum.fetch_add(b.vaccination_weights[i], std::memory_order_acq_rel);
        }
    }
    // Regardless of vaccination strategy, perform a final pass to remove weight from those already vaccinated
    for (int i = 0; i < b.agents.size(); i++) {
        if (b.agents[i].status == V || b.agents[i].vaccination_progress || b.agents[i].status == R) {
            b.vaccination_weight_sum -= b.vaccination_weights[i];
            b.vaccination_weights[i] = 0;
        }
    }
}

/**
* Converts the latitude and longitude to Radians for Distance Calculation.
* @param degree Takes the Coordinate to translate. 
*/

long double toRadians(const long double degree)
{
    //M_PI is cmath predefined value for pi.
    long double one_deg = (M_PI) / 180;
    return (one_deg * degree);
}
/**
* Calculates the distance between two sets of coordinates.
* @param loc1 Takes a tuple consisting of latitude and longitude values
* @param loc1 Takes a tuple consisting of latitude and longitude values
*/
int toDistance(std::tuple<double,double> loc1, std::tuple<double,double> loc2)
{
    /*double lat1 = 59.85882;
    double long1 = 17.63889;
    double lat2 = 59.4439;
    double long2 = 18.06872;
    */
    double lat1 = std::get<0>(loc1);
    double long1 = std::get<1>(loc1);
    double lat2 = std::get<0>(loc2);
    double long2 = std::get<1>(loc2);

    lat1 = toRadians(lat1);
    long1 = toRadians(long1);
    lat2 = toRadians(lat2);
    long2 = toRadians(long2);

    // Haversine Formula for calculating distance
    double difflong = long2 - long1;
    double difflat = lat2 - lat1;

    double ans = pow( sin(difflat / 2), 2) + cos(lat1) * cos(lat2) * pow( sin(difflong / 2), 2);

    ans = 2 * asin(sqrt(ans));

    // Sets radius of the earth, 6371 for KM and 3956 for miles  
    double R = 6371;
    ans = ans * R; //Calculate final distance.

    return static_cast<int>(ans); //Returns int
}

/**
* Calculates weights for a community based on size and distance to other communities. 
* @param comm The vector containing all the communities. 
* @param curr The community to calculate weights for. 
*/
std::vector<double> calculateWeight(std::vector<board>& comm,board& curr) {
    std::vector<double> weights = {};
    #pragma omp parallel for ordered
    for (int i = 0; i < comm.size(); i++) {
        double population = comm[i].dim* comm[i].dim;
        std::tuple<double, double> current_latlong = curr.lat_long;
        std::tuple<double, double> target_latlong = comm[i].lat_long;
        int distance = toDistance(current_latlong, target_latlong);
        if (distance == 0) {
            #pragma omp ordered
            weights.push_back(1);
        }
        else {
            /*std::cout << std::setprecision(15) << std::fixed; //To check values with higher precision
            std::cout << std::endl << "pop " << population << " distance " << distance << std::endl;
            std::cout << std::endl << "The weight is " << population / distance << std::endl; */
            #pragma omp ordered
            weights.push_back(population / distance);
        }
    }
    return weights;
}







/**
 * Vaccinates a number of agents in the board, based on the board's vaccination weights
 * @param board The board to perform vaccinations in
 * @param n_vaccinations The number of vaccinations to perform
 */
void vaccinate(board& b, unsigned int n_vaccinations) {
    // Make sure we don't perform more vaccinations than there are unvaccinated agents
    n_vaccinations = std::min((size_t)n_vaccinations, b.agents.size() - b.vaccinations_started);
    if (n_vaccinations == 0 || b.vaccination_weight_sum == 0) {
        return;
    }

    // For each agent, compute the combined weight of it and all agents with a lower index
    // This lets us later use lower_bound to do a binary search for the agent to vaccinate
    std::vector<uint64_t> cumulative_weights(b.agents.size());
    unsigned int weight_sum = 0;
    for (int i = 0; i < b.agents.size(); i++) {
        weight_sum += b.vaccination_weights[i];
        cumulative_weights[i] = weight_sum;
    }
    assert(cumulative_weights[b.agents.size() - 1] = b.vaccination_weight_sum);

    std::uniform_int_distribution<unsigned int> dis(1, (int)b.vaccination_weight_sum); // Generates (cumulative) weight values determining who gets vaccinated
    std::vector<std::atomic_flag> vaccinated(b.agents.size());  // Stores which agents have been vaccinated this round
    unsigned int collisions = 0;

#pragma omp parallel for
    for (int i = 0; i < n_vaccinations; i++) {
        std::mt19937_64& gen = generators[omp_get_thread_num()];
        uint64_t n = dis(gen);
        // Finds the first agent with cumulative weight greater than n
        auto it = std::lower_bound(cumulative_weights.begin(), cumulative_weights.end(), n);
        unsigned int agent = it - cumulative_weights.begin();
        // Atomically check if this agent has been vaccinated this round
        if (!vaccinated[agent].test_and_set()) {
            vaccinate(b.agents[agent]);
            b.vaccinations_started.fetch_add(1, std::memory_order_acq_rel);
            b.vaccination_weight_sum.fetch_sub(b.vaccination_weights[agent], std::memory_order_acq_rel);
            b.vaccination_weights[agent] = 0;
        }
        else {
            // If it's already been vaccinated, we will need to do this vaccination over
            collisions++;
        }
    }
    if (collisions > 0) {
        // Redo the vaccinations that collided. It might be better to use a collision-free algorithm here.
        vaccinate(b, collisions);
    }
}


/**
* Picks a random community to travel to based on weights.
* @param items The list of weights from the given community to another.
* @param gen The generator object that is used to randomize numbers. 
*/
int weightRand(std::vector<double> items, std::mt19937_64& gen) {

    double cumulative = 0;
    for (int i = 0; i < items.size(); i++)
    {
        cumulative += items[i];
    }
    std::uniform_int_distribution<int> move_dist(0, cumulative);
    //std::uniform_real_distribution move_dist(0.0, cumulative);
    int rand = move_dist(gen);
    for (int i = 0; i < items.size(); i++) {
        cumulative -= items[i];
        if (rand >= cumulative) {
            return i;
        }
    }
    std::cout << std::endl << "The system failed to find a weight in the double version of weightRand" << std::endl;
    exit(1);
}
/**
* Picks a random community based on population size.
* @param items The list of community sizes.
* @param gen The generator object that is used to randomize numbers.
*/
int weightRand(std::vector<int> items, std::mt19937_64& gen) {

    double cumulative = 0;
    for (int i = 0; i < items.size(); i++)
    {
        cumulative += items[i];
    }

    std::uniform_real_distribution move_dist(0.0, cumulative);
    double rand = move_dist(gen);
    for (int i = 0; i < items.size(); i++) {
        cumulative -= items[i];
        if (rand >= cumulative) {
            return i;
        }
    }
    std::cout << std::endl << "The system failed to find a weight in the int version of weightRand" << std::endl;
    exit(1);
}






void updateBoard(board& from, board& to, agent agentFrom, agent agentTo) {
    from.status_counts[agentTo.status]++;
    from.status_counts[agentFrom.status]--;
    to.status_counts[agentFrom.status]++;
    to.status_counts[agentTo.status]--;
}


/**
* Move functionality, allowing agents to be moved between communities at random based of weights and population sizes.
* @param curr_board The list of current communities.
* @param agents The amount of agents to be swapped.
* @param weight The list of community population sizes. 
*/
void moveAgents(std::vector<board>& curr_board, int agents, std::vector<int> weight) {
    std::mt19937_64 gen = generators[0];

    for (int i = 0; i < agents; i++) {
        int fromBoardIdx = weightRand(weight, gen);
        board& fromBoard = curr_board[fromBoardIdx];
        int fromDim = fromBoard.dim;

        std::uniform_int_distribution<int> fromAgent_dis(0, (fromDim * fromDim - 1));
        int agentfromBoardIdx = fromAgent_dis(gen);
        agent agentFromBoard = fromBoard.agents[agentfromBoardIdx];
        int toBoardIdx = weightRand(curr_board[fromBoardIdx].weights, gen);
        board& toBoard = curr_board[toBoardIdx];
        int targetDim = toBoard.dim;

        std::uniform_int_distribution<int> targetAgent_dis(0, (targetDim * targetDim - 1));
        int agentToBoardIdx = targetAgent_dis(gen);
        agent agentToBoard = toBoard.agents[agentToBoardIdx];
        swap(fromBoard, toBoard, agentfromBoardIdx, agentToBoardIdx);
        updateBoard(fromBoard, toBoard, agentFromBoard, agentToBoard);
    }
}

int main() {
    
    // Create random generators
    unsigned int n_threads;
    // seems like we have to run omp_get_num_threads in a parallel region to get the actual number of threads
#pragma omp parallel
    n_threads = omp_get_num_threads();
    std::random_device rd;
    for (int i = 0; i < n_threads; i++) {
        generators.push_back(std::mt19937_64(rd() + i));
    }



    //std::vector<std::string> comm_names = { "Uppsala", "Stockholm", "Eskilstuna" }; //Provided by user
    //std::vector<double> weight = {             0.34,       0.33,       0.33 }; //Provided by user
    //std::vector<std::vector<double>> inter_weight = { {1.0, 0.7, 0.3}, {0.7, 1.0, 0.3}, {0.4, 0.6, 1.0} }; //Provided by user
    //int population = 100000; //Provided by user


    //read_data_long = (cities, coordinates, population)
    std::tuple<std::vector<std::string>, std::vector<std::tuple<double, double>>, std::vector<int>> csv_data = read_data_from_csv();

    std::vector<int> dimensions = {};

    std::vector<board> prev_board = {};
    std::vector<board> curr_board = {};
    std::vector<int> populations = {};
    // Stores overall counts each step. This data is plotted after the simulation.
    std::vector<std::vector<unsigned int>> status_history(STATES_COUNT, std::vector<unsigned int>(MAX_TIME + 1, 0));
    std::vector<unsigned int> infection_history(MAX_TIME + 1, 0);

    std::vector<std::string> comm_names = std::get<0>(csv_data);
    std::vector<std::tuple<double, double>> coordinates = std::get<1>(csv_data);

    for (int i = 0; i < comm_names.size(); i++)
    {
        int population = std::get<2>(csv_data)[i];
        int calc_dim = ceil(sqrt(population));
        dimensions.push_back(calc_dim);
        std::vector<double> weights = {};
        board new_board(calc_dim, STARTER_AGENTS, AGENT_TYPES, comm_names[i], coordinates[i], weights, generators[0]);
        prev_board.push_back(new_board);
        board curr = new_board;
        curr_board.push_back(curr);
        populations.push_back(calc_dim * calc_dim);
        for (int s = 0; s < STATES_COUNT; s++) {
            status_history[s][0] += curr_board[i].status_counts[s];
        }
        infection_history[0] += curr_board[i].total_infections;
    }
    for (int i = 0; i < curr_board.size(); i++)
    {
        std::vector<double> weights = calculateWeight(curr_board, curr_board[i]);
        curr_board[i].weights = weights;
        prev_board[i].weights = weights;
    }
    

    /*for (int i = 0; i < curr_board.size(); i++)
    {
        std::cout << "Hello name : " << std::get<0>(curr_board[i].lat_long) << " " << std::get<1>(curr_board[i].lat_long) << std::endl;
    }*/
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

    if (remove("..\\lib\\built_covid_data\\coviddata.csv") != 0) {
        perror("Error deleting file");
    }
    else {
        puts("File successfully deleted");
        std::ofstream file;

        file.open("..\\lib\\built_covid_data\\coviddata.csv", std::ios::app);

        // Initialize first row
        file << "city" << ";" << "popu" << ";" << "lat" << ";" << "long" << ";" << "month" << ";" << "agent" << std::endl;
        file.close();
    }

    for (unsigned int t = 0; t < MAX_TIME; t++)
    { //Loop tracking
        if (t % 10 == 0) {
            for (int i = 0; i < curr_board.size(); i++)
            {
                visualization_of_board(curr_board[i], t);
            }

        }
        //moveAgents(curr_board, SWAP_AMOUNT, populations);
        //The magic number is how many agents should swap each timestep.
        moveAgents(curr_board, SWAP_AMOUNT, populations);
        prev_board = curr_board;

        for (int i = 0; i < comm_names.size(); i++)
        {
            prev_board[i] = curr_board[i];

            step(prev_board[i], curr_board[i], t);
            // print_board(prev_board[i], comm_names[i], t);
            if (t > VACCINATION_START) {
                // TODO: we can optimize this by avoiding to update weights when there are no vaccinations left to do
                update_vaccination_weights(curr_board[i]);
                vaccinate(curr_board[i], VACCINATIONS_PER_DAY);
            }
        }


        // Save status counts
        for (int i = 0; i < comm_names.size(); i++) {
            for (int s = 0; s < STATES_COUNT; s++) {
                status_history[s][t + 1] += curr_board[i].status_counts[s];
            }
            infection_history[t + 1] += curr_board[i].total_infections;
        }

    } // /for t
    {
        using namespace matplot;

        std::vector<std::vector<unsigned int>> plot_data{
            status_history[I], status_history[A], status_history[S], status_history[V], status_history[R]
        };
        auto f = figure();
        auto ax = f->current_axes();
        area(ax, plot_data);
        title(ax, "Agent count per status");
        xlabel(ax, "t (days)");
        ylabel(ax, "population");
#ifndef _WIN32 //Must be set in allcaps to work
        legend(ax, { "i", "a", "s", "v", "r" });
#endif

        // Plot the total number of infections
        auto f2 = figure();
        auto ax2 = f2->current_axes();
        plot(ax2, infection_history);
        title(ax2, "Cumulative infections (incl. asymptomatic)");
        xlabel(ax2, "t (days)");
        ylabel(ax2, "Number of infections");
        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << std::endl << "It took  " << time_span.count() << " seconds." << std::endl;
    }
    std::cout << "Press Enter to exit..." << std::endl;
    std::cin.get();
    return 0;
}

#include <vector>
#include <omp.h>
#include <cstdlib>
#include <time.h>
#include <iostream>
#include <atomic>
#include <matplot/matplot.h>
#include <cmath>
#include <random>
#include <chrono>
#include <tuple>
#include <algorithm>

#define INFECTION_RADIUS 50
#define RECOVERY_RATE 14
#define INFECTION_PROBABILITY 25
#define DIM 1000
#define MAX_TIME 140
#define QUARANTINE 5
#define Q_FLAG false
#define LAMBDA 2.5
#define ONE_CHANCE true

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
    for(int i = 0; i < DIM*DIM; i++){
        if(i % DIM == 0){
            std::endl (std::cout);
        }
        std::cout <<  b.agents[i].status << " ";
    }
}

int main() {
    board previous = {
        DIM,
        std::vector<agent>(DIM*DIM)
    };
    srand(time(NULL));
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<int> dis(0, (DIM*DIM-1));
    int seeded = 0;
    // Func to spawn starter infecitous agents
    while (seeded != 4) {
        int pz = dis(gen);
        previous.agents[pz].status = I;
        std::cout << "pz is: " << pz;
        seeded++;
        

    }
    /*int pZ = std::rand() % DIM * (DIM - 1);
    previous.agents[pZ].status = I;
    pZ = std::rand() % DIM*(DIM-1);
    previous.agents[pZ].status = I;
    pZ = std::rand() % DIM * (DIM - 1);
    previous.agents[pZ].status = I;
    pZ = std::rand() % DIM * (DIM - 1);
    previous.agents[pZ].status = I;*/
    board current = previous;
	std::atomic_int sus = DIM*DIM-4;
	std::atomic_int rem = 0;
	std::atomic_int inf = 4;
	std::vector<double> susp = {};
    std::vector<double> infe = {};
    std::vector<double> remo = {};

    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

	for (unsigned int t = 0; t < MAX_TIME; t++)
	{ //Loop tracking time
       
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
                        rem++;
                        inf--;
                        continue;
                    }
                    //{ Check neighbours
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
                        std::uniform_int_distribution<int> dis3(0, size-1);
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
                            std::uniform_int_distribution<int> dis2(0, 100);
                            int prob = dis2(gen);
                            if (prob <= INFECTION_PROBABILITY) {
                                otherCurr.status = I;
                                inf++;
                            }
                            
                        }
                        else {
                            
                            checked.erase(checked.begin() + rand_S);
                            
                            if (checked.size() == 0 ) { continue; }
                            
                            goto Repeat;
                            
                        }

                    }
                    else{
                        for (int y_box = -rad; y_box <= rad; y_box++) {
                            for (int x_box = -rad; x_box <= rad; x_box++) {
                                int y_other = y + y_box;
                                int x_other = x + x_box;
                                if (y_other < 0 || y_other >= DIM)  { //Checks for Bounds
                                    break;
                                } else if (x_other < 0 || x_other >= DIM || (x_other == x && y_other == y)) {
                                    continue;
                                }
                                agent& other = previous.agents[y_other * DIM + x_other];
                                agent& otherCurr = current.agents[y_other * DIM + x_other];
                                if (other.status == S && otherCurr.status != I) { //If neighbour is susceptible
                                    int prob = std::rand() % 100;
                                    if(prob < INFECTION_PROBABILITY){
                                        otherCurr.status = I;
 									    inf++;
                                    
								    }
                                    //goto Exit; //Testing theory
                                }
                            }

                        }
                    
                }
                }
            //Exit: //Testing theory regarding beta compared to paper.
                //int testing = 0; //Required for line above
            }
        }
        // TODO: optimize
        previous = current;
        /*if(t % 10 == 0) {
            std::cout << std::endl << "---- t: " << t;
            print_board(current);
        }*/
		sus = DIM*DIM - rem - inf;
		susp.push_back(sus);
		remo.push_back(rem);
		infe.push_back(inf);
        //std::cout << t << std::endl;
        // print_board(previous);
        // std::cout << " --------- "  << std::endl;
	} // /for t
	{
        using namespace matplot;

		std::vector<std::vector<double>> Y = {susp, remo, infe};
        std::cout << "sus " << sus << "rem " << rem << "infe " << inf<<std::endl;
            plot(Y);
            title("infected people");
            xlabel("t (days)");
            ylabel("population");
//#ifndef _WIN32
//    legend({ "s", "r", "i" });
//#endif
        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << std::endl<<"It took  " << time_span.count() << " seconds."<<std::endl;
        show();
    }
    
    return 0;
}


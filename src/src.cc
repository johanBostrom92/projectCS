#include <vector>
#include <omp.h>
#include <cstdlib>
#include <time.h>
#include <iostream>
#include <atomic>
#include <matplot/matplot.h>
#include <cmath>
#include <random>

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
    int pZ = std::rand() % DIM * (DIM - 1);
    previous.agents[pZ].status = I;
    pZ = std::rand() % DIM*(DIM-1);
    previous.agents[pZ].status = I;
    pZ = std::rand() % DIM * (DIM - 1);
    previous.agents[pZ].status = I;
    pZ = std::rand() % DIM * (DIM - 1);
    previous.agents[pZ].status = I;
    board current = previous;
	std::atomic_int sus = DIM*DIM-1;
	std::atomic_int rem = 0;
	std::atomic_int inf = 1;
	std::vector<double> susp = {};
    std::vector<double> infe = {};
    std::vector<double> remo = {};

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
                    Repeat:
                        std::random_device rd;
                        std::mt19937_64 gen(rd());
                        std::uniform_int_distribution<int> dis(-rad, rad);
                        int x_rand = dis(gen);
                        int y_rand = dis(gen);
                        int y_other = y_rand + y;
                        int x_other = x_rand + x;
                        if (y_other < 0 || y_other >= DIM || x_other < 0 || x_other >= DIM) {
                            goto Repeat;
                        }
                        agent& other = previous.agents[y_other * DIM + x_other];
                        agent& otherCurr = current.agents[y_other * DIM + x_other];
                        if (other.status == S && otherCurr.status != I) { //If neighbour is susceptible
                            int prob = std::rand() % 100;
                            if (prob < INFECTION_PROBABILITY) {
                                otherCurr.status = I;
                                inf++;
                            }
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
		
            plot(Y);
            title("infected people");
            xlabel("t (days)");
            ylabel("population");
#ifndef _WIN32
    legend({ "s", "r", "i" });
#endif
        show();
    }
    return 0;
}


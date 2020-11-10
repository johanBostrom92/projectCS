#include <vector>
#include <omp.h>
#include <cstdlib>
#include <time.h>
#include <iostream>
#include <atomic>

#define INFECTION_RADIUS 1
#define RECOVERY_RATE 14
#define INFECTION_PROBABILITY 5
#define DIM 100
#define MAX_TIME 200

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

int main() {
    board previous = {
        DIM,
        std::vector<agent>(DIM*DIM)
    };
    srand(time(NULL));
    int pZ = std::rand() % DIM*DIM-1;
    previous.agents[pZ].status = I;
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
            std::cout << std::endl << "---- t: " << t;
            print_board(current);
        }
		sus = sus - rem - inf;
	} // /for t
}


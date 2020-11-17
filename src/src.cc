#include <vector>
#include <omp.h>
#include <cstdlib>
#include <time.h>
#include <iostream>
#include <atomic>
#include <random>
#include <matplot/matplot.h>

#define INFECTION_RADIUS 1
#define RECOVERY_RATE 5
#define INFECTION_PROBABILITY 25
#define DIM 5
#define MAX_TIME 10
#define STARTER_AGENTS 1

enum agent_status {
        S, I, R, X
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
        if (b.agents[i].status == X) {
            std::cout << "X" << " ";
        }
        else {
            std::cout << b.agents[i].status << " ";
        }
    }
}


//          move(uppsala_prev, sthlm_prev, 0);
void move(board& previous, board& current, int idx) {
    agent sthlm_agent = current.agents[idx];
    current.agents[idx] = previous.agents[idx];
    previous.agents[idx] = sthlm_agent;
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
                            if (prob < INFECTION_PROBABILITY) {
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
    //while (seeded != STARTER_AGENTS) {
    //    int pZ = dis(gen);
    //    if (uppsala_prev.agents[pZ].status != I) {
    //        uppsala_prev.agents[pZ].status = I;
    //        std::cout << "pz is: " << pZ;
    //        seeded++;
    //    }

    //}
    uppsala_prev.agents[0].status = X;

    //previous.agents[pZ].status = I;
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
        if (t == 1) {
            move(uppsala_curr, sthlm_curr, 0);
        }
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
		/*    {1, 3, 1, 2}, {5, 2, 5, 6}, {3, 7, 3, 1}};
            */
		// auto f = gcf();
        // f->width(f->width() * 2);
        /*
        subplot(1, 2, 0);
        area(Y);
        title("Stacked");
        legend({"S", "R", "I"});
        */
        //subplot(1, 2, 1);
        //auto handles = plot(Y);
        //handles[0]->marker(line_spec::marker_style::point);
        //handles[1]->marker(line_spec::marker_style::point);
        //handles[2]->marker(line_spec::marker_style::point);
        //title("infected people");
        //xlabel("t (days)");
        //ylabel("population");
        //legend({"s", "r", "i"});


        //show();
    }
    return 0;
}


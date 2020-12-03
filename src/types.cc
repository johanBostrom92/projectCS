#include "types.hh"
#include <random>

board::board(unsigned int dim, unsigned int initial_infections, const std::vector<agent_type> agent_types, double weight, std::string name)
    : dim(dim),
    agents(dim*dim),
    sus(dim*dim - initial_infections),
    rem(0),
    inf(initial_infections),
    asymp(0),
    vacc(0),
    weight(weight),
    name(name){

    std::default_random_engine rand_generator;

    // Generate initial infections
    std::uniform_int_distribution<int> dis(0, (dim*dim-1));
    int seeded = 0;
        while (seeded != initial_infections) {
            int pz = dis(rand_generator);
            if (agents[pz].status != I) {
                agents[pz].status = I;
                seeded++;
            }

        }


    // Get the sum of all type weights
    double weight_sum = 0.0f;
    for (int i = 0; i < agent_types.size(); i++) {
        auto& type = agent_types[i];
        weight_sum += type.weight;
    }

    // Randomize an agent type for every agent, then generate and assign a radius using that type
    std::uniform_real_distribution type_dist(0.0, weight_sum);
    for(auto& agent : agents) {
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
}

board::board(board&& other) {
   *this = std::move(other);
}

board::board(const board& other) {
   *this = other;
}

board& board::operator=(board&& other) {
    if (this != &other) {
        this->agents = std::move(other.agents);
        this->dim = other.dim;
        this->name = other.name;
        this->weight = other.weight;
        other.dim = 0;

        this->sus = other.sus.load();
        this->rem = other.rem.load();
        this->inf = other.inf.load();
        this->asymp = other.asymp.load();
        this->vacc = other.vacc.load();
    }
    return *this;
}

board& board::operator=(const board& other) {
    if (&other != this) {
        this->agents = other.agents;
        this->dim = other.dim;
        this->name = other.name;
        this->weight = other.weight;

        this->sus = other.sus.load();
        this->rem = other.rem.load();
        this->inf = other.inf.load();
        this->asymp = other.asymp.load();
        this->vacc = other.vacc.load();
    }
    return *this;
}
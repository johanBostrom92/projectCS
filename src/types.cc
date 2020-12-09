#include "types.hh"
#include <random>

board::board(unsigned int dim, unsigned int initial_infections, const std::vector<agent_type> agent_types, double weight, std::string name)
    : dim(dim),
    agents(dim*dim),
    total_infections(initial_infections),
    weight(weight),
    name(name),
    vaccination_weights(dim*dim, 1),
    vaccination_weight_sum(dim*dim),
    vaccinations_started(0) {

    for (int s = 0; s < STATES_COUNT; s++) {
        status_counts[s] = 0;
    }
    status_counts[S] = dim*dim - initial_infections;
    status_counts[I] = initial_infections;

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
        this->total_infections = other.total_infections.load();
        this->vaccination_weights = std::move(other.vaccination_weights);
        this->vaccination_weight_sum = other.vaccination_weight_sum.load();
        this->vaccinations_started = other.vaccinations_started.load();
        this->dim = other.dim;
        this->name = other.name;
        this->weight = other.weight;
        other.dim = 0;

        for (int s = 0; s < STATES_COUNT; s++) {
            this->status_counts[s] = other.status_counts[s].load();
        }
    }
    return *this;
}

board& board::operator=(const board& other) {
    if (&other != this) {
        this->agents = other.agents;
        this->total_infections = other.total_infections.load();
        this->vaccination_weights = other.vaccination_weights;
        this->vaccination_weight_sum = other.vaccination_weight_sum.load();
        this->vaccinations_started = other.vaccinations_started.load();
        this->dim = other.dim;
        this->name = other.name;
        this->weight = other.weight;

        for (int s = 0; s < STATES_COUNT; s++) {
            this->status_counts[s] = other.status_counts[s].load();
        }
    }
    return *this;
}
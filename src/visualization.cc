#include <vector>
#include "types.hh"
#include "parameters.hh"
#include <matplot/matplot.h>
#include <omp.h>
#include <fstream>
#include <iostream>
#include <stdio.h>



std::tuple<std::vector<std::string>, std::vector<std::tuple<double, double>>, std::vector<int>>  read_data_from_csv() {
    std::vector<std::string> cities = {};
    std::vector<std::tuple<double, double>> coordinates = {};
    std::vector<int> populations = {};


    
    
    std::ifstream data_file;
    std::string file_name = CITIES_INPUT;
    std::string path = "..//lib//cities_data//" + file_name;
    data_file.open(path);
    std::cout << "Successfully read file " << path << std::endl;
    std::string name;
    std::string lat;
    std::string lat_long;
    std::string country;
    std::string population;

    //Skip first line on purpose to access data
    std::getline(data_file, name);

    while (data_file.good()) {

        std::getline(data_file, name, ';');
        std::getline(data_file, lat, ';');
        std::getline(data_file, lat_long, ';');
        std::getline(data_file, country, ';');
        std::getline(data_file, population, '\n');
       

        std::tuple tmp = std::make_tuple(stod(lat), stod(lat_long));
        if (!name.empty()){
              cities.push_back(name);
            coordinates.push_back(tmp);
            populations.push_back(std::stoi(population));
        }
    }
    data_file.close();

    std::tuple<std::vector<std::string>, std::vector<std::tuple<double, double>>, std::vector<int>> res = std::make_tuple(cities, coordinates, populations);
    return res;
}


/**
 * Translates a agent status single character to a string
 * @param a agent_status
 * @return a string of the translated status
 */
std::string translate_agent_status(agent_status status) {
    if (status == S) {
        return "Susceptible";
    }     if (status == A) {
        return "Asymptomatic";
    }     if (status == V) {
        return "Vaccinated";
    }    if (status == I) {
        return "Infected";
    }    if (status == R) {
        return "Recovered";
    }
}


void write_data_to_csv(std::string city, agent_status status[], std::vector<unsigned int> size, std::tuple<double,double> lat_long, unsigned int time_unit) {

    // create our file
    std::ofstream file;

    file.open("..//lib//built_covid_data//coviddata.csv", std::ios::app);


    for (int i = 0; i < STATES_COUNT; i++) {
   
        
        //TODO: send tuple!
        std::string trans_status = translate_agent_status(status[i]);
  
        file << city << ";" << size[i]/SCALE << ";" << std::get<0>(lat_long) << ";" << std::get<1>(lat_long) << ";" << time_unit << ";" << trans_status << std::endl;

    }

    //closing 
    file.close();
}


void visualization_of_board(board b, unsigned int t) {
    if (PLOT) {

        std::vector<unsigned int> population = {};
        agent_status status_name[5] = { S, A, I, V, R };

        for (int i = 0; i < STATES_COUNT; i++)
        {
            population.push_back(b.status_counts[i]);
        }
        write_data_to_csv(b.name, status_name, population, b.lat_long, t);

    }
}


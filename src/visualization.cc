#include <vector>
#include "types.hh"
#include "parameters.hh"
#include <matplot/matplot.h>
#include <omp.h>
#include <fstream>
#include <iostream>


void scatterplot(board b){
     //Transform array index to X,Y cordinates as col,row

        // X,Y for S
        std::vector<double> col_s;
        std::vector<double> row_s;

        // X,Y for A
        std::vector<double> col_a;
        std::vector<double> row_a;

        // X,Y for V
        std::vector<double> col_v;
        std::vector<double> row_v;

        // X,Y for I
        std::vector<double> col_i;
        std::vector<double> row_i;

        // X,Y for R
        std::vector<double> col_r;
        std::vector<double> row_r;

        #pragma omp parallel for ordered
        for(int i = 0; i < DIM*DIM; i++){
            if (b.agents[i].status == S) {
                #pragma omp ordered
                {
                    row_s.push_back(i / DIM);
                    col_s.push_back(i % DIM);
                }
            }
            else if (b.agents[i].status == A) {
                #pragma omp ordered
                {
                    row_a.push_back(i / DIM);
                    col_a.push_back(i % DIM);
                }
            }
            else if (b.agents[i].status == V) {
                #pragma omp ordered
                {
                    row_v.push_back(i / DIM);
                    col_v.push_back(i % DIM);
                }
            }
            else if (b.agents[i].status == I) {
                #pragma omp ordered
                {
                    row_i.push_back(i / DIM);
                    col_i.push_back(i % DIM);
                }
            }
            else if (b.agents[i].status == R) {
                #pragma omp ordered
                {
                    row_r.push_back(i / DIM);
                    col_r.push_back(i % DIM);
                }
            }
            else {
                std::cout <<  "Error didn't find any match ";
            }

        }
        {
            //scatterplots using matplot++ on converted array -> XY cordinates. 
            using namespace matplot;
            
            int size = 8; //TODO fix dynamical size depending on table DIM size. 
        
            //dummy algoritm
            // set size of graph window
            //based on DIM calculate amount of circles given a radius 
            // set size as calculated radius 

            hold(on);
            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    if (!col_s.empty()) {
                        auto scat_s = scatter(col_s, row_s, size);
                        scat_s->marker_color({ 0, 0, 0 });
                        scat_s->marker_face_color({ 0.2, 0.4, 1 });
                    } 
                }
                 #pragma omp section
                {
                    if (!col_a.empty()) {
                        auto scat_a = scatter(col_a, row_a, size);
                        scat_a->marker_color({ 0, 0, 0 });
                        scat_a->marker_face_color({ 1, 0.3, 0.3 });
                    }
                }

                #pragma omp section
                {
                    if (!col_v.empty()) {
                        auto scat_v = scatter(col_v, row_v, size);
                        scat_v->marker_color({ 0, 0, 0 });
                        scat_v->marker_face_color({ 0.84, 0.733, 0.36 });
                    }
                }
                #pragma omp section
                {
                    if (!col_i.empty()){
                        auto scat_i = scatter(col_i, row_i, size);
                        scat_i->marker_color({ 0, 0, 0 });
                        scat_i->marker_face_color({ 1, 0, 0 });
                    }
                }
                #pragma omp section
                {
                    if (!col_r.empty()) {
                        auto scat_r = scatter(col_r, row_r, size);
                        scat_r->marker_color({ 0, 0, 0 });
                        scat_r->marker_face_color({ 0, 1, 0 });
                    }
                }
            }
           show();
        }
}


void heatmap_plot(board b){
        using namespace matplot;
/*
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z; 
   contourf(x, y, z, {2, 3}); //->contour_text(true)
   */ 
  // /* 
   std::vector<std::vector<double>> x = {{1,2},{1,2}};
    std::vector<std::vector<double>> y = {{1,2},{1,2}};
    std::vector<std::vector<double>> z = {{5,2},{1,2}};

  
    contourf(x,y,z);
   // */
    show();
}



void geo_plot(){
    using namespace matplot; 
    std::vector<double> lat = {60};
    std::vector<double> lon = {18}; 
    std::vector<double> size = {(1)}; 
    std::vector<double> type = {0}; 


    geolimits({55, 75}, {10, 25});
   auto log_size = transform(size, [](double x) { return log(x)/(10000000000000); });
   std::cout << "size is: " << log_size[0]; 
   hold(on);
   auto gb = geobubble(lat, lon, log_size); 
    

    
    // geolimits({45, 62}, {-155, -120});

    show();
}

std::tuple<double, double> get_lat_long(std::string city){
    using namespace std; 
    
    transform(city.begin(), city.end(), city.begin(), ::tolower);
   
    ifstream data_file;
    data_file.open("..\\lib\\cities_data\\cities_swe.csv");

    string name;
    string lat;
    string lat_long;


    while(data_file.good()){

        getline(data_file, name, ';');
        getline(data_file, lat, ';');
        getline(data_file, lat_long, '\n');
        
        if (name.compare(city) == 0){
            tuple res = make_tuple(stod(lat),stod(lat_long));
            data_file.close();
            return res;
        }
    }
    data_file.close();

return make_tuple(0,0);
}


std::string translate_agent_status(agent_status status){
    if (status == S){
        return "Susceptable";
    }     if (status == A){
        return "Asymptotic";
    }     if (status == V){
        return "Vaccinated";
    }    if (status == I){
        return "Infected";
    }    if (status == R){
        return "Recovered";
    }
}


void write_data_to_csv(std::string city[], agent_status status[], int size[], int time_unit[], int size_of_array){
    
 
    // create our file
    std::ofstream file;

    file.open("..\\lib\\built_covid_data\\coviddata.csv",std::ios::app);				

    // Initialize first row
    //file << "city" << ";" << "popu" << ";" << "lat" << ";" << "long"  << ";" << "month" << ";" << "agent" << std::endl;

    for (int i = 0; i < size_of_array; i++){
        std::cout << "kom hit " << i << std::endl;
        std::tuple values = get_lat_long(city[i]);
        std::string trans_status = translate_agent_status(status[i]);
        file << city[i] << ";" << size[i] << ";" << std::get<0>(values) << ";" << std::get<1>(values)  << ";" << time_unit[i] << ";" << trans_status << std::endl;

    }
  
    //closing 
    file.close();
}

void visualization_of_board(board b){
    if (PLOT){
    //    heatmap_plot(b);
    
      //  geo_plot();
    //scatterplot( b);


    std::string test[4] = {"stockholm", "Uppsala", "Malmoe", "Goeteborg"};
    agent_status test2[4] = {A,S,S,A}; 
    int test_pop[4] = {10000, 5000, 666, 9999}; 
    int time_unit[4] = {1,2,3,4}; 
    write_data_to_csv(test, test2, test_pop, time_unit, 4);
  
    }
}


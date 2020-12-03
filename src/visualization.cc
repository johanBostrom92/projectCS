#include <vector>
#include "types.hh"
#include "parameters.hh"
#include <matplot/matplot.h>
#include <omp.h>


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

        #pragma omp parallel for 
        for(int i = 0; i < DIM*DIM; i++){
            if (b.agents[i].status == S) {
                row_s.push_back(i / DIM);
                col_s.push_back(i % DIM);
            }
            else if (b.agents[i].status == A) {
                row_a.push_back(i / DIM);
                col_a.push_back(i % DIM);
            }
            else if (b.agents[i].status == V) {
                row_v.push_back(i / DIM);
                col_v.push_back(i % DIM);
            }
            else if (b.agents[i].status == I) {
                row_i.push_back(i / DIM);
                col_i.push_back(i % DIM);
            }
            else if (b.agents[i].status == R) {
                row_r.push_back(i / DIM);
                col_r.push_back(i % DIM);
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


void heatmap_plot(){
        using namespace matplot;

   // auto [X, Y, Z] = peaks();
   // contourf(X, Y, Z, {2, 3}); //->contour_text(true)
    std::vector<std::vector<double>> x;
    std::vector<std::vector<double>> y;
    std::vector<std::vector<double>> z;

   x = {{1,2,3,4,5},{1,2,3,4,5}};
   y = {{1,2,3,4,5},{1,2,3,4,5}};
   z = {{1,2,3,4,5},{1,2,3,4,5}};
    contourf(x,y,z);
    show();
}



void geo_plot(){
    using namespace matplot; 
    std::vector<double> lat = {60};
    std::vector<double> lon = {18}; 
    std::vector<double> size = {(1)}; 
    std::vector<double> type = {0}; 


    geolimits({55, 75}, {10, 25});
   auto log_size = transform(size, [](double x) { return log(x + 2); });
    auto gb = geobubble(lat, lon, log_size); 
    

    
    // geolimits({45, 62}, {-155, -120});

    show();
}

void visualization_of_board(board b){
    if (PLOT){
       // heatmap_plot();
       // tsunami_example();
       // geo_plot();
    scatterplot( b);
    }
}

#include "parameters.hh"
#include <vector>
#include "types.hh"

/**
 * Used to call visualisation
 * @param a board b
 */
void visualization_of_board(board b); 

/**
 * The function translates our array to vectors of x,y cordinates and categories them based on agent
 * @param a board
 * prints the board using a scatter plot visualizing neighbours that can infect each others
 */ 
void scatterplot(board b);

/**
 * Writes covid data to a csv file
 * @param a string array with name of affected cities
 * @param a agent_status array with the status of agents
 * @param a int array of the amount of affected agents given a specific status
 * @param a array of the current time unit 
 * @param a int for the size of the current arrays, all arrays needs to be the same size as this int
 * writes data to a CSV file which will be read by of geoplotter.py
 */
void write_data_to_csv(std::string city[], agent_status status[], int size[], int time_unit[], int size_of_array);
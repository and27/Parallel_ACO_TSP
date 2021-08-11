#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>

#include "ants.h"

using namespace std;

//create an ant structure to store ant related data (visited, path)
struct Ant {
  float tourLength;  
  bool tabu[MAX_NODES];
  int path[MAX_NODES];
  int curNode, nextNode, pathIndex;
};

Ant antColony[MAX_NODES];
float pheromone[MAX_NODES][MAX_NODES];
float best = (float)MAX_TOUR;
int bestIndex;

//this can be implemented as a function member or constructor if we construct the ant class
void init(){
  for (int from = 0; from < MAX_NODES; from ++){
   for (int to=0; to < MAX_NODES; to ++){
     pheromone[from][to] = INIT_PHER;
   }
  }

  //Initialize the departure node on all the antColony
  for (int ant=0; ant<MAX_ANTS; ant++){
    antColony[ant].curNode = rand()%MAX_NODES;
    //Put all the cities available (tabu = false (0))
    for (int to = 0; to <MAX_NODES; to++){
      antColony[ant].tabu[to]=0;
      }
    antColony[ant].pathIndex = 1; 
    antColony[ant].path[0] = antColony[ant].curNode;
    antColony[ant].nextNode = -1;
    antColony[ant].tourLength=0;
    //put the selected city as visited (tabu = true (1))
    antColony[ant].tabu[antColony[ant].curNode] = 1; 
  }
}


void restartAnts() {
  for (int ant = 0; ant < MAX_ANTS; ant++) {
    //At the begining we put the nextNode as disabled (-1)	  
    antColony[ant].nextNode = -1;
    antColony[ant].tourLength = 0.0;
    for (int i = 0; i < MAX_NODES; i++) {
      antColony[ant].tabu[i] = 0;
      antColony[ant].path[i] = -1;
    }
    //Select the first node randomly
    antColony[ant].curNode = rand() % MAX_NODES;
    antColony[ant].pathIndex = 1;
    antColony[ant].path[0] = antColony[ant].curNode;
    antColony[ant].tabu[antColony[ant].curNode] = 1;
  }
}

//Calculate node selection probability
float calculateProbability(int from, int to, EdgeMatrix *dist){
   return (pow(pheromone[from][to], ALPHA)* pow((1.0/(*dist)[from][to]), BETA));
}

//Select the next node using the calculateProbability function (p)
int selectNextCity(int ant, EdgeMatrix *dist){
  int from = antColony[ant].curNode;
  float sum = 0.0;
  for (int to = 0; to <MAX_NODES; to++){
    if (antColony[ant].tabu[to]==0){
      sum+=calculateProbability(from, to, dist);
    }
   }

  int lastBestIndex = 0.0; 
  float luckyNumber = (float)rand()/RAND_MAX;
  //cout << "random: " << luckyNumber<<endl; 

  for (int to = 0; to < MAX_NODES; to++){
    if (antColony[ant].tabu[to]==0){
      float product = calculateProbability(from, to, dist) /sum;
      sleep(0.5);
      if(product > 0 ){
      // cout << "probability: "  << product <<endl; 
       luckyNumber-= product; 
       lastBestIndex=to;
       	
       if(luckyNumber<=0.0){
	  //cout << "Selected node "<< to<<endl;
        return to;
       }
      }
    }
   }
  cout << "not reached" << endl; 
  return lastBestIndex;
 }


void updatePheromone(){
 int from, to, i, ant;
 for (from = 0; from <MAX_NODES; from ++){
  for (to=0; to<from; to++){
   pheromone[from][to]*=1.0-RHO;
   
   if (pheromone[from][to]<0.0){
    pheromone[from][to]=INIT_PHER;
   }
   pheromone[to][from] = pheromone[from][to];
  }
 }

  for (ant = 0; ant < MAX_ANTS; ant++) {
    for (i = 0; i < MAX_NODES; i++) {
     from = antColony[ant].path[i];
      if (i < MAX_NODES - 1) {
      to = antColony[ant].path[i+1];
      } else {
      to = antColony[ant].path[0];
       }
      pheromone[from][to] += (QVAL / antColony[ant].tourLength);
      pheromone[to][from] = pheromone[from][to];
    }
     }
} 

float euclideanDistance(int x1, int x2, int y1, int y2) {
  int xd = x1 - x2;
  int yd = y1 - y2;
  return (int) (sqrt(xd * xd + yd * yd) + 0.5);
}

float pseudoEuclideanDistance(int x1, int x2, int y1, int y2) {
  int xd = x1 - x2;
  int yd = y1 - y2;
  float rij = sqrt((xd * xd + yd * yd) / 10.0);
  return ceil(rij);
}

void constructTSP(std::string graph, cityType *cities, EdgeMatrix *dist) {
  // Load cities from file
  std::ifstream infile((graph + ".tsp").c_str());
  std::string line;
  bool euclidean = true; // whether to run EUC_2D or ATT distance metric

  int city, x, y;
  bool reading_nodes = false;
  while (std::getline(infile, line)) {
    istringstream iss(line);
    string word;
    if (!reading_nodes) {
      iss >> word;
      if (word.compare("EDGE_WEIGHT_TYPE") == 0) {
        iss >> word >> word;
        std::cout << "edge type: " << word << std::endl;
        euclidean = !word.compare("EUC_2D");
      } else if (word.compare("NODE_COORD_SECTION") == 0) {
        reading_nodes = true;
      }
    } else if (iss >> city >> x >> y) {
      cities[city-1].x = x;
      cities[city-1].y = y;
    }
  }
  infile.close();

  // Compute distances between cities (edge weights)
  for (int from = 0; from < MAX_NODES; from++) {
    (*dist)[from][from] = 0.0;

    for (int to = from + 1; to < MAX_NODES; to++) {
      float edge_dist;
      if (euclidean) {
        edge_dist = euclideanDistance(cities[from].x, cities[to].x, cities[from].y, cities[to].y);
      } else {
        edge_dist = pseudoEuclideanDistance(cities[from].x, cities[to].x, cities[from].y, cities[to].y);
      }
      if (edge_dist == 0) {
        edge_dist = 1.0;
      }
      printf("edge[%d][%d] = %f\n", from, to, edge_dist);
      (*dist)[from][to] = edge_dist;
      (*dist)[to][from] = edge_dist;
    }
  }
}

   
//Select the best found solution
int bestSolution(){
  for (int ant = 0; ant < MAX_ANTS; ant++) {
    if (antColony[ant].tourLength < best) {
      best = antColony[ant].tourLength;
      printf("new best: %1.f\n", best);
      bestIndex = ant;
    }
  }
  return best;
}
 


int main(int argc, char *argv[]){
  
  cityType cities[MAX_NODES];
  string graph = argv[1];
  EdgeMatrix *dist = new EdgeMatrix();
  //float *dist = new float[MAX_NODES * MAX_NODES];
  constructTSP(graph, cities, dist);
  init();
  for (int i =0; i<MAX_ITERATIONS; i++){
    //For each tour we run the following actions and -- > Update Pheromone to improve the next tour	  
    for (int ant = 0; ant < MAX_ANTS; ant++){
      while (antColony[ant].pathIndex < MAX_NODES){
      //cout << "Lets select the next city " << endl;
        antColony[ant].nextNode = selectNextCity(ant, dist);
        antColony[ant].tabu[antColony[ant].nextNode]=1;
        //cout << "current -> " << antColony[ant].curNode << endl;
        //cout << "next -> " << antColony[ant].nextNode << endl;
        antColony[ant].path[antColony[ant].pathIndex++] = antColony[ant].nextNode;
        antColony[ant].tourLength += (*dist)[antColony[ant].curNode][antColony[ant].nextNode];   
        //Now the current node is the node selected (next node)
        antColony[ant].curNode = antColony[ant].nextNode;
      }
     antColony[ant].tourLength += (*dist)[antColony[ant].path[MAX_NODES-1]][antColony[ant].path[0]]; 
    }
    bestSolution();
    updatePheromone();
    if (i != MAX_ITERATIONS - 1) {
      restartAnts();
    }// end ants for
  }//End for of TOUR ITERATIONS
}

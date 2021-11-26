/***************************************************************************
    local_search.cpp
    (C) 2021 by C.Blum & M.Blesa
    
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../HEADERS/Timer.h"
#include "../HEADERS/Random.h"
#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <vector>
#include <set>
#include <limits>
#include <iomanip>

// global variables concerning the random number generator (in case needed)
time_t t;
Random* rnd;

// Data structures for the problem data
int n_of_nodes;
int n_of_arcs;
vector< set<int> > neighbors;

// string for keeping the name of the input file
string inputFile;

// number of applications of local search
int n_apps = 1;

// dummy parameters as examples for creating command line parameters -> 
// see function read_parameters(...)
int dummy_integer_parameter = 0;
int dummy_double_parameter = 0.0;


inline int stoi(string &s) {

  return atoi(s.c_str());
}

inline double stof(string &s) {

  return atof(s.c_str());
}

void read_parameters(int argc, char **argv) {

    int iarg = 1;

    while (iarg < argc) {
        if (strcmp(argv[iarg],"-i")==0) inputFile = argv[++iarg];
        
        // reading the number of applications of local search 
        // from the command line (if provided)
        else if (strcmp(argv[iarg],"-n_apps")==0) n_apps = atoi(argv[++iarg]); 
        
        // example for creating a command line parameter param1 -> 
        // integer value is stored in dummy_integer_parameter
        else if (strcmp(argv[iarg],"-param1")==0) {
            dummy_integer_parameter = atoi(argv[++iarg]); 
        }
        // example for creating a command line parameter param2 -> 
        // double value is stored in dummy_double_parameter
        else if (strcmp(argv[iarg],"-param2")==0) {
            dummy_double_parameter = atof(argv[++iarg]);  
        }
        iarg++;
    }
}

void MergeSortedIntervals(vector<set<int> >& v, int s, int m, int e) {

    // temp is used to temporary store the vector obtained by merging
    // elements from [s to m] and [m+1 to e] in v
    vector<set<int> > temp;

    int i, j;
    i = s;
    j = m + 1;

    while (i <= m && j <= e) {

        if (v[i].size() <= v[j].size()) {
            temp.push_back(v[i]);
            ++i;
        }
        else {
            temp.push_back(v[j]);
            ++j;
        }

    }

    while (i <= m) {
        temp.push_back(v[i]);
        ++i;
    }

    while (j <= e) {
        temp.push_back(v[j]);
        ++j;
    }

    for (int i = s; i <= e; ++i)
        v[i] = temp[i - s];

}

// the MergeSort function
// Sorts the array in the range [s to e] in v using
// merge sort algorithm
void MergeSort(vector<set<int> >& v, int s, int e) {
    if (s < e) {
        int m = (s + e) / 2;
        MergeSort(v, s, m);
        MergeSort(v, m + 1, e);
        MergeSortedIntervals(v, s, m, e);
    }
}

int computeH(const set<int>& v, set<int>& S) {

    int compt = 0;
    for (int a : v) {                                   // per tot vei de a
        std::set<int>::iterator it = S.find(a);
        if (it != S.end()) ++compt;                     // Veins de v a S
    }

    int n = ceil(float(v.size())/2);                    // upper bound de grau(v)/2
    return n-compt;
}

int cover_degree(const vector<set<int> >& G, set<int>& S, set<int>& SC, int i) {

    int argMax = 0;
    int aux = 0;
    int nCoverMax = -1;

    for (int a : G[i]) {

        std::set<int>::iterator it = SC.find(a);
        int covered = 0;
        if (it != SC.end()) {
            for (int aresta : G[a]) {
                if (computeH(G[aresta], S) > 0) {
                    ++covered;
                    aux = aresta;
                }
            }

            if (covered > nCoverMax) {
                nCoverMax = covered;
                argMax = a;
            }
        }
    }
    return argMax;
}

void greedy(const vector<set<int> >& G, set<int>& S) {

    int n = G.size();
    set<int> SC;
    for (int i = 0; i < n; ++i) {
        SC.insert(i);
    }

    for (int i = 0; i < n; ++i) {

        int p = computeH(G[i], S);                  // mirem si el vertex i esta cobert
        if (p > 0) {
        
            for (int j = 0; j < p; ++j) {
                int argmax = cover_degree(G, S, SC, i);             // cover degree dels vertexs adjacents a i que estan en SC
                S.insert(argmax);
                SC.erase(argmax);
            }
        }
    }
}

double getHeuristic(const set<int> &solution, const vector<set<int>> &neighbors, int& weight, int s, int e) {
    int aux = weight;
    for(set<int>::iterator it = neighbors[s].begin(); it != neighbors[s].end(); ++it) {
        if (solution.find(*it) == solution.end()) aux += 1;
    }
    for(set<int>::iterator it = neighbors[e].begin(); it != neighbors[e].end(); ++it) {
        if (solution.find(*it) == solution.end()) aux -= 1;
    }
    return aux;
}

bool is_Solution(const set<int> &solution, const vector<set<int>> &neighbours, const int &j) {
    //cout << "Deleted from solution vertex " << j << endl;
    for (int i : neighbours[j]) {
        int needed_neighbors = round(neighbours[i].size() / 2.0);
        //cout << "Seeing if neighbor " << i << " has enough vertex. (Needed: " << neighbours[i].size() / 2.0 << ")." << endl;
        //cout << i << "'s first neighbor is " << *(neighbours[i].begin()) << endl;
        int counted_neighbors = 0;
        for (set<int>::iterator it = neighbours[i].begin(); it != neighbours[i].end(); it++) {
            if (solution.find(*it) != solution.end()) {
                //cout << "FOUND A NEIGHBOR: " << *it << endl;
                counted_neighbors++;
            }
        }
        //cout << "Found " << counted_neighbors << " counted neighbors." << endl << endl;
        if (counted_neighbors < needed_neighbors) return false;
    }
    return true;
}


void swap(set<int>& aux,int i, int j) {
    aux.erase(i);
    aux.insert(j);
}

bool eliminateVertex(set<int> &actual, vector< set<int>>& neighbors, set<int>& newsol, set<int>::iterator& last, set<int>& not_sol, int& weight) {
    vector<set<int>> successors;
    for (set<int>::iterator it = last; it != actual.end(); ++it) {
       set<int> aux = actual;
       aux.erase(*it);
       if (is_Solution(aux, neighbors, *it )) {
           newsol = aux;
           last = it;
           not_sol.insert(*it);
            for(set<int>::iterator it2 = neighbors[*it].begin(); it2 != neighbors[*it].end(); ++it2) {
                if (actual.find(*it2) == actual.end()) weight -= 1;
            }
            cout << *it << endl;
           return true;
       }
    }
    return false;
}

bool nextSwap(set<int> &actual, vector< set<int>>& neighbors, set<int>& newsol, set<int>& not_sol, int& weight) {

    double heurActual = weight;
    double heurNewSol;
    pair<int, int> bestSwap;
    bool primer = true;
    int n = neighbors.size();
    for (set<int>::iterator it = actual.begin(); it != actual.end(); ++it) {
        for (set<int>::iterator it2 = not_sol.begin(); it2 != not_sol.end(); ++it2) {
            set<int> aux = actual;
            if (actual.find(*it2) == actual.end()) {
                swap(aux, *it,  *it2);
                if (is_Solution(aux, neighbors, *it )) {
                    cout << "Swapping vertex " << *it << " for vertex " << *it2 << endl;
                    if (primer) {
                        primer = false;
                        bestSwap.first = *it;
                        bestSwap.second = *it2;
                        heurNewSol = getHeuristic(aux, neighbors, weight, *it, *it2);
                    } else {
                        double heurAux = getHeuristic(aux, neighbors, weight, *it, *it2);
                        if ( heurAux > heurNewSol) {
                            bestSwap.first = *it;
                            bestSwap.second = *it2;
                            heurNewSol = heurAux;
                            cout << *it << " " << *it2 << endl;
                        }
                    }
                }
            }
        }
    }

    if (heurNewSol > heurActual) {
        swap(actual, bestSwap.first, bestSwap.second);
        not_sol.insert(bestSwap.first);
        not_sol.erase(bestSwap.second);
        newsol = actual;
        for(set<int>::iterator it = neighbors[bestSwap.first].begin(); it != neighbors[bestSwap.first].end(); ++it) {
            if (actual.find(*it) == actual.end()) weight += 1;
        }
        for(set<int>::iterator it = neighbors[bestSwap.second].begin(); it != neighbors[bestSwap.second].end(); ++it) {
            if (actual.find(*it) == actual.end()) weight -= 1;
        }
        cout << weight <<endl;
        return true;
    }else {
        return false;
    }
    
}

/**********
Main function
**********/

int main( int argc, char **argv ) {

        read_parameters(argc,argv);
        
        // setting the output format for doubles to 2 decimals after the comma
        std::cout << std::setprecision(2) << std::fixed;

        // initializing the random number generator.
        // A random number in (0,1) is obtained with: double rnum = rnd->next();
        rnd = new Random((unsigned) time(&t));
        rnd->next();

        // vectors for storing the result and the computation time 
        // obtained by the <n_apps> applications of local search
        vector<double> results(n_apps, std::numeric_limits<int>::max());
        vector<double> times(n_apps, 0.0);

        // opening the corresponding input file and reading the problem data
        ifstream indata;
        indata.open(inputFile.c_str());
        if(not indata) { // file couldn't be opened
            cout << "Error: file could not be opened" << endl;
        }

        indata >> n_of_nodes;
        indata >> n_of_arcs;
        neighbors = vector< set<int> >(n_of_nodes);
        int u, v;
        while(indata >> u >> v) {
            neighbors[u - 1].insert(v - 1);
            neighbors[v - 1].insert(u - 1);
        }
        indata.close();

        // main loop over all applications of local search
        for (int na = 0; na < n_apps; ++na) {

            // the computation time starts now
            Timer timer;

            // Example for requesting the elapsed computation time at any moment: 
            // double ct = timer.elapsed_time(Timer::VIRTUAL);

            cout << "start application " << na + 1 << endl;
            

            // Solució inicial (tots els vertex)

            set<int> S;                 // S will contain the final solution
            vector<set<int>> initial_neighbors = neighbors;
            /*
            MergeSort(neighbors, 0, neighbors.size()-1);            // ordenem primer els vertexs respecte el seu grau de manera ascendent
            greedy(neighbors, S);
            */

            for (int i = 0; i < neighbors.size(); i++) S.insert(i);  // Solucio inicial plena

            set<int> actual = S;

            bool end = false;

            /*
            for(int i = 0; i < initial_neighbors.size(); ++i) {
                cout << "Vertex :" << i << "  Veins: ";
                for (set<int>::iterator it = initial_neighbors[i].begin(); it != initial_neighbors[i].end(); ++it) {
                    cout << *it << " ";
                }
                cout << endl;
            }
            */
           int weight = 0;

           for(int i = 0; i < neighbors.size(); ++i) {
                if (actual.find(i) == actual.end()) {
                    for (set<int>::iterator it = neighbors[i].begin(); it != neighbors[i].end(); ++it) {
                        if (actual.find(*it) == actual.end()) weight += 1;
                    }
                }
            }

        set<int>::iterator last = actual.begin();
            
        cout << actual.size() << endl;

        set<int> solaux;

        set<int> not_sol;
        for (int i = 0; i < initial_neighbors.size(); i++) {
            if (actual.find(i) == actual.end()) not_sol.insert(i);
        }


        while (not end) {
            
            if (eliminateVertex(actual, initial_neighbors, solaux, last, not_sol, weight)) {
                cout << "Escollit erase" << endl;
                actual = solaux;
            } else {
                //cout << "Starting swap" << endl;
                if (nextSwap(actual, initial_neighbors, solaux, not_sol, weight)) {
                    cout << "Escollit swap" << endl;
                    actual = solaux;
                } else {
                    end = true;
                }
                last = actual.begin();
            }
            //cout << actual.size() << endl;
        }

        cout << actual.size() << endl;

        
        
        for (set<int>::iterator it = actual.begin(); it != actual.end(); ++it) { 
            cout << *it << " ";
        }
        cout << endl;
        
        
        // HERE GOES YOUR LOCAL SEARCH METHOD

        // The starting solution for local search may be randomly generated, 
        // or you may incorporate your greedy heuristic in order to produce 
        // the starting solution.
        
        // Empty initial solution
        // Whenever you move to a new solution, first take the computation 
        // time as explained above. Say you store it in variable ct.
        // Then, write the following to the screen: 
        // cout << "value " << <value of the current solution>;
        // cout << "\ttime " << ct << endl;

        // When a local minimum is reached, store the value of the 
        // corresponding solution in vector results: 
        // results[na] = <value of the local minimum>;
        
        // Finally store the needed computation time (that is, the time 
        // measured once the local minimum is reached) in vector times: 
        // times[na] = ct;

        cout << "end application " << na + 1 << endl;
    }

    // calculating the average of the results and computation times, and 
    // their standard deviations, and write them to the screen
    double r_mean = 0.0;
    int r_best = std::numeric_limits<int>::max();
    double t_mean = 0.0;
    for (int i = 0; i < int(results.size()); i++) {
        r_mean = r_mean + results[i];
        if (int(results[i]) < r_best) r_best = int(results[i]);
        t_mean = t_mean + times[i];
    }
    r_mean = r_mean/double(results.size());
    t_mean = t_mean/double(times.size());
    double rsd = 0.0;
    double tsd = 0.0;
    for (int i = 0; i < int(results.size()); i++) {
        rsd = rsd + pow(results[i]-r_mean,2.0);
        tsd = tsd + pow(times[i]-t_mean,2.0);
    }
    rsd = rsd/double(results.size());
    if (rsd > 0.0) {
        rsd = sqrt(rsd);
    }
    tsd = tsd/double(results.size());
    if (tsd > 0.0) {
        tsd = sqrt(tsd);
    }
    cout << r_best << "\t" << r_mean << "\t" << rsd << "\t";
    cout << t_mean << "\t" << tsd << endl;
}
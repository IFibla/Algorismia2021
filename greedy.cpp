/***************************************************************************
    greedy.cpp 
    (C) 2021 by C. Blum & M. Blesa
    
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

#include "./HEADERS/Timer.h"
#include "./HEADERS/Random.h"
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

// dummy parameters as examples for creating command line parameters 
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
        
        // example for creating a command line parameter param1 
        //-> integer value is stored in dummy_integer_parameter
        else if (strcmp(argv[iarg],"-param1")==0) 
            dummy_integer_parameter = atoi(argv[++iarg]); 
            
        // example for creating a command line parameter param2 
        //-> double value is stored in dummy_double_parameter
        else if (strcmp(argv[iarg],"-param2")==0) 
            dummy_double_parameter = atof(argv[++iarg]);  
            
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

int computeH(const set<int> v, set<int>& S) {
    int compt = 0;

    for (int a : v) {

        std::set<int>::iterator it = S.find(a);
        if (it != S.end()) ++compt;

    }
     int n = v.size()/2; 

    return n-compt;

}

int cover_deegree(const vector<set<int> >& G, const set<int>& veins, set<int>& S) {
    int covered = 0;
    for (int a : veins) {
        if (computeH(G[a], S) > 0) ++covered;
    }
    return covered;
}

void greedy(const vector<set<int> >& G, set<int>& S) {

    int n = G.size();
    set<int> SC;

    for (int i = 1; i <= n; ++i) {
        SC.insert(i);
    }

    for (int i = 0; i < n; ++i) {

        int p = computeH(G[i], S);

        if (p > 0) {
            
            for (int j = 0; j < p; ++j) {

                int nCoverMax = -1;

                for (int a : G[i]) {

                    int nCover = 0;
                    std::set<int>::iterator it = SC.find(a);
                    if (it != SC.end()) nCover = cover_deegree(G, G[a], S);
                    if (nCover > nCoverMax) nCoverMax = nCover;
                }

                S.insert(nCoverMax);
                SC.erase(nCoverMax);
            }

        }
        
    }
    cout << S.size() << " " << SC.size() << endl;
}

/************
Main function
*************/

int main( int argc, char **argv ) {

    read_parameters(argc,argv);
    
    // setting the output format for doubles to 2 decimals after the comma
    std::cout << std::setprecision(2) << std::fixed;

    // initializing the random number generator. 
    // A random number in (0,1) is obtained with: double rnum = rnd->next();
    rnd = new Random((unsigned) time(&t));
    rnd->next();

    // variables for storing the result and the computation time 
    // obtained by the greedy heuristic
    double results = std::numeric_limits<int>::max();
    double time = 0.0;

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

    // the computation time starts now
    Timer timer;

    // Example for requesting the elapsed computation time at any moment: 
    // double ct = timer.elapsed_time(Timer::VIRTUAL);

    // HERE GOES YOUR GREEDY HEURISTIC
    // When finished with generating a solution, first take the computation 
    // time as explained above. Say you store it in variable ct.
    // Then write the following to the screen: 
    

    set<int> S;                 // S will contain the final solution

    MergeSort(neighbors, 0, neighbors.size()-1);            // O(nlg(n))
    greedy(neighbors, S);

    cout << "VERTEXS a D" << endl;
    for (int i : S) {
        cout << "Vertex " << i << endl;
    }



    // cout << "value " << <value of your solution> << "\ttime " << ct << endl;

}


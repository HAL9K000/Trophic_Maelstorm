/**
 * The cpp file does the following:
 * 
 * Input a 1D array of size L*L, representing a 2D grid of size L*L.
 * Input a grid size L and a grid distance dx.
 * Input an array of radii R = [r1, r2, ..., rN].
 * Let r = max(R).
 * Then for each lattice site i*L + j, we want to find the number and indices of lattice sites within a norm 2 radius of r from i*L + j..
 * This output will be stored in a 2D array output[] of size (L*L)*(pi*r*r/(dx*dx)).
 * For each row in this output array, we will store the indices of the lattice sites within the norm 2 radius of r from the lattice site i*L + j.
 * This output will be stored in ascending order of the norm 2 distance from the lattice site i*L + j (with the lattice site i*L + j itself being the first entry).
 * Also note this norm 2 distance calculation will obey periodic boundary conditions.
 * Finally, in a seperate array count_r of size N, we will store the number of lattice sites within a norm 2 radius of ri from the lattice site i*L + j.
 * More specifically count_r[k] will store the number of lattice sites within a norm 2 radius of rk from the lattice site i*L + j (1 <= k <= N).
 * output[] and count_r[] will be passed by reference.
 * output[] will be a 2D array of size (L*L)*(pi*r*r/(dx*dx)) and count_r[] will be an array of size N.
 * output will be initialised to -1, and count_r will be initialised to 0.
 * 
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <omp.h>
#include <chrono>
#include <random>
#include <map>
#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstdlib>
#include <bits/stdc++.h>
using namespace std::chrono;
using namespace std;

typedef std::vector<std::vector <double>> D2Vec_Double;
typedef std::vector<std::vector <int>> D2Vec_Int;

typedef std::vector<std::vector <std::vector <double>>> D3Vec_Double;
typedef std::vector<std::vector <std::vector <int>>> D3Vec_Int;

typedef std::vector<std::vector <std::vector <std::vector <double>>>> D4Vec_Double;
typedef std::vector<std::vector <std::vector <std::vector <int>>>> D4Vec_Int;

struct coordinates {
   int x;
   int y;
};

struct coord_dist {
   int x;
   int y;
   double dist;
};

// Custom comparator for sorting pairs based on the squared distance
bool comparePairs(const pair<int, int>& a, const pair<int, int>& b) {
    return (a.first * a.first + a.second * a.second) < (b.first * b.first + b.second * b.second);
}

// Function to calculate the norm 2 distance between two lattice sites i*L + j and k*L + l, obeying periodic boundary conditions.

double norm2_dist(int i, int j, int k, int l, int L)
{
    int x_dist = abs(i - k);
    int y_dist = abs(j - l);
    x_dist = min(x_dist, L - x_dist);
    y_dist = min(y_dist, L - y_dist);
    return sqrt(x_dist*x_dist + y_dist*y_dist);
}

vector<pair<int, int>> computeNeighboringSitesCentral(int R) 
{
    vector<pair<int, int>> centralNeighboringSites;
    // Compute neighboring sites for the central lattice site (0, 0)
    for (int x = -R; x <= R; ++x) 
    {
        for (int y = -R; y <= R; ++y) 
        {
            int distanceSquared = x * x + y * y;
            if (distanceSquared <= R * R) 
            {
                centralNeighboringSites.emplace_back(x, y);
            }
        }
    }
    // Sort the vector of neighboring sites based on squared distance
    // Sort the centralNeighboringSites based on squared distance
    sort(centralNeighboringSites.begin(), centralNeighboringSites.end(), comparePairs);

    return centralNeighboringSites;
}

// Function to get neighbors within a ball of radius fR
vector<int> getNeighborsInBall(const vector<int>& sortedNeighboringSites, double f) {
    int numNeighbors = sortedNeighboringSites.size();
    int maxDistanceIndex = int(f*f*numNeighbors);
    return vector<int>(sortedNeighboringSites.begin(), sortedNeighboringSites.begin() + maxDistanceIndex);
}

vector<int> generateNeighboringSitesFromCentral(const vector<pair<int, int>>& centralNeighboringSites, int i, int j, int L) {
    vector<int> neighboringSites;

    // Generate neighboring sites for the lattice site (i, j)
    for (const auto& offset : centralNeighboringSites) {
        int dx = offset.first;
        int dy = offset.second;

        int nx = (i + dx + L) % L;
        int ny = (j + dy + L) % L;
        neighboringSites.push_back(nx * L + ny);
    }

    return neighboringSites;
}

int main() 
{
    // Set the grid size (L) and the radius of the circle (R)
    int L = 20;
    int R = 8;

    // Compute neighboring sites for the central lattice site (0, 0)
    vector<pair<int, int>> centralNeighboringSites = computeNeighboringSitesCentral(R);

    // Example: Generate neighboring sites for lattice site (2, 2)
    vector<int> resultFor2_2 = generateNeighboringSitesFromCentral(centralNeighboringSites, 2, 2, L);

    // Print the result
    cout << "Neighboring sites for (2, 2): ";
  
    //Format the output such that the elements of the vector are outputted as concentric circles centered at (2, 2).
    double dist = 0;
    for (int i = 0; i < resultFor2_2.size(); i++)
    {
        if (i == 0)
        {
            cout << "Centre: (" << resultFor2_2[i] / L << ", " << resultFor2_2[i] % L << ")";
        }
        else
        {   
            double dist_i = norm2_dist(2, 2, resultFor2_2[i] / L, resultFor2_2[i] % L, L);
            dist_i *= dist_i;
            if (dist_i == dist)
            {
                
                cout << " (" << resultFor2_2[i] / L << ", " << resultFor2_2[i] % L << ")";
            }
            else
            {
                dist = dist_i;
                cout << endl;
                cout << "Norm 2 Distance Sq " << dist << ": ";
                cout << " (" << resultFor2_2[i] / L << ", " << resultFor2_2[i] % L << ")";
            }
        }
    }
    cout << endl;

    return 0;
}
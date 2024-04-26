/**
 * @file test_par_hdf5.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2024-04-22
 * 
 * @copyright Copyright (c) 2024
 * 
 * This is a test file for the parallel HDF5 implementation.
 * In this script, the following will occur:
 * 1. Using OpenMP, 10 threads will be created.
 * 2. Each thread will have multiple "replicates" given by n (achieved by a for loop), and in each replicate, a 2D vector of fixed size "L" will be created.
 * 3. In each replicate, The 2D vector will be filled with random numbers between 0 and 1.
 * 4. The 2D vector will be written to an HDF5 file, and will be preceded by a header.
 * 5. Each thread will correspond to a different subfolder (or group) in the HDF5 file, named "/thread_i" where i is the thread number.
 * 6 The HDF5 file will be chunked and compressed, and will be initiallised with "n" number of chunks (corresponding to the number of replicates).
 * 7. In each replicate i, the 2D vector will be written to the  i-th chunk.
 * 8. The HDF5 file will be closed.
 * 9. Mutex locks will be used to ensure that the HDF5 file is not written to by multiple threads at the same time.
 * 
 */

//#include <hdf5.h>


#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <omp.h>
#include <mutex>
#include "H5Cpp.h"

using namespace std;
using namespace std::chrono;

// Define mutex for synchronization
std::mutex mtx;
std::mutex hypemtx;

// HDF5 dimensions
const int L = 10; // Size of the 2D vector
const int n = 5;  // Number of replicates

// Function to generate random numbers between 0 and 1
double getRandomNumber() 
{
    static std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    return distribution(generator);
}

// Function to write 2D vector to HDF5 file
void writeDataToHDF5_parallel(H5::H5File& file, int threadNum) 
{
    std::vector<std::vector<double>> data(L, std::vector<double>(L));

    // Fill the 2D vector with random numbers
    for (int i = 0; i < L; ++i) 
    {
        for (int j = 0; j < L; ++j) 
        {
            data[i][j] = getRandomNumber();
        }
    }
    
    //H5::H5File file(filename, H5F_ACC_RDWR | H5F_ACC_CREAT);
    


    // Create group for the current thread if it doesn't exist

    std::string groupName = "/thread_" + std::to_string(threadNum);
    H5::Group group;
    {
        std::lock_guard<std::mutex> lock(mtx);
        if (!file.exists(groupName)) 
            group = file.createGroup(groupName);
        
        else 
            group = file.openGroup(groupName);
    }

    // Define dataset properties
    hsize_t dims[2] = {L, L};
    H5::DataSpace dataspace(2, dims);
    H5::DSetCreatPropList plist;
    plist.setChunk(2, dims);
    plist.setDeflate(6);

    // Create or open the dataset
    std::string datasetName = groupName + "/data";
    H5::DataSet dataset;
    if (!group.exists("data")) {
        std::lock_guard<std::mutex> lock(mtx);
        dataset = group.createDataSet("data", H5::PredType::NATIVE_DOUBLE, dataspace, plist);
    } else {
        dataset = group.openDataSet("data");
    }

    


    
    // Write data to dataset
    {
        std::lock_guard<std::mutex> lock(hypemtx);
        // Extend dataset to accommodate new chunk
        hsize_t currentDims[2];
        dataset.getSpace().getSimpleExtentDims(currentDims);
        hsize_t newDims[2] = {currentDims[0] + 1, L};
        dataset.extend(newDims);
        hsize_t start[2] = {currentDims[0], 0};
        hsize_t count[2] = {1, L};
        H5::DataSpace memspace(2, count);
        // Select hyperslab for writing data
        dataset.getSpace().selectHyperslab(H5S_SELECT_SET, count, start);
        // Write data to dataset
        dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE, memspace, dataspace);
    }

    // Close the HDF5 file
    file.close();
}

int main() {
    // Set number of threads for OpenMP
    int numThreads = 10;
    omp_set_num_threads(numThreads);

    // Start clock
    auto start = std::chrono::high_resolution_clock::now();

    H5::H5File file("test_data.h5", H5F_ACC_RDWR | H5F_ACC_CREAT);
    cout << "File opened" << endl;

    // Execute parallel region
#pragma omp parallel
    {
        // Get thread number
        int threadNum = omp_get_thread_num();

        // Write data to HDF5 file
        writeDataToHDF5_parallel(file, threadNum);
    }

    // Stop clock
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = duration_cast<seconds>(end - start);
    cout<< "Time taken: " << elapsed.count() << " s" << endl;

    return 0;
}
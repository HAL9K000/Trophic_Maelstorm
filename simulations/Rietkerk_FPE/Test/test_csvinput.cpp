#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;
int g= 128;
int Sp= 4;

// Function to extract a specific column from a CSV file and store it as doubles in a vector
void extractColumn(std::vector <std::vector<double>> &array, const std::string& filename, const std::vector<int> &columns, std::vector <std::vector<double>> &const_ind_val) 
{
    // Create an empty 2D vector with Sp rows (and no columns) called result.
    std::vector<std::vector <double>> result(3, vector<double>(1, 0.0));
    /**
    for(int s=0; s< Sp; s++)
    {
        result.push_back({0});
    }*/

    // Open the CSV file
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) 
    {
        std::cerr << "Error opening file: " << filename << std::endl;
    }



    // Check for a header line and discard it

    std::string line;
    int columnIndex = columns[columns.size()-1];
    while (std::getline(inputFile, line)) 
    {
        std::stringstream ss(line);
        std::string cell;

        if(line.find("a_c") != std::string::npos)
		{   cout << "FUCK YOU" <<endl;   continue;}

        //Store last column index from columns vector
        
        // Extract the desired columns
        for (int i = 0; i <= columnIndex; ++i) 
        {
            if (!std::getline(ss, cell, ',')) 
            {
                std::cerr << "Error: Column index " << columnIndex << " out of bounds." << std::endl;
                inputFile.close();
                
            }
            //Check if i is in columns vector
            if (std::find(columns.begin(), columns.end(), i) != columns.end()) 
            {
                // Store the cell value in the result vector as double
                double value;
                try 
                {
                    value = std::stod(cell);
                    result[i-2].push_back(value);
                } 
                catch (const std::invalid_argument& e) 
                {
                    std::cerr << "Error: Unable to convert cell value to double. Column index: " << i << std::endl;
                } 
                catch (const std::out_of_range& e) 
                {
                    std::cerr << "Error: Cell value out of range for double. Column index: " << i << std::endl;
                }
            }
        }

    }

    // First store result in array.

    for(int s=0; s< 1; s++)
    {
        for(int i=0; i< g*g; i++)
        {
            array[s][i] = result[s][i+1];
        }
    }

    for(int i=0; i< g*g; i++)
    {
            array[Sp-2][i] = result[1][i+1];    // Sotring W in the second last row of array
            array[Sp-1][i] = result[2][i+1];    // Sotring O in the last row of array
    }
    

    //Print result
    for(int s=0; s< 3; s++)
    {
        cout << "For species in RESULT:\t" << s <<   " with size:\t" << result[s].size() <<endl;
        for(int i=0; i< g*g; i++)
        {
            //cout << result[s][i+1] << " ";
        }
        cout << endl;
    }

    // const_ind_val is a 2D vector of size (Sp - columns.size()) x 2, which stores the constant index 
    //for each column not being extracted (and thus not in result) and the constant value to be substituted in the array at that index value.

    // These values are pre-assigned in the main function.

    // Now, we need to substitute the constant values in the array.

    for(int i=0; i< const_ind_val.size(); i++)
    {
        int index = const_ind_val[i][0];
        double value = const_ind_val[i][1];

        for(int j=0; j< g*g; j++)
        {
            array[index][j] = value;
        }
    }

    // Close the file
    inputFile.close();

    //Free up memory allocated to result
    vector <vector<double>>().swap(result);
    for(int s=0; s< Sp; s++)
    {
        //result[s].clear();
        //vector <double>().swap(result[s]);
    }


}

int main() 
{
    stringstream t, a, R; t << 83176.2; a  << 0.045; R << 1;
    std::string filename= "Test_CSVs/Rogue_T_" + t.str() + "_A_" + a.str() + "_R_" + R.str() + ".csv";
    int columnIndex;

    std::vector <std::vector<double>> array(Sp, vector<double> (g*g, 0.0));

    // Get user input for filename and column index
    
    vector <int> cols ={2,3,4};
    vector <vector<double>> const_ind_val = {{1, 0.65}};
    // Extract the specified column as doubles
    extractColumn(array, filename, cols, const_ind_val);

    // Dislay all values in array
    for(int s=0; s< Sp; s++)
    {
        cout << "For species:\t" << s << endl;
        for(int i=g*g-10; i< g*g; i++)
        {
            cout << array[s][i] << " ";
        }
        cout << endl;
    }

    return 0;
}
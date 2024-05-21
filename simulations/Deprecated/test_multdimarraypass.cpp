#include "test_multdimarraypass.h"


//const int L =tempL;
/**
void determine_neighbours_R2(int g, int gm, int neighbours_R2[][L*L][2][3])
{
    for(int x = 0; x < L*L; x++)
    {
        for(int y=0; y < 2; y++)
        {
            for(int z= 0; z< 3; z++)
            {
                neighbours_R2[0][x][y][z] = 1; neighbours_R2[1][x][y][z] = 1;
            }
        }

    }
}
*/
void determine_neighbours_R2(int g, int gm, auto &neighbours_R2)
{
    for(int x = 0; x < L*L; x++)
    {
        for(int y=0; y < 2; y++)
        {
            for(int z= 0; z< 2; z++)
            {
                neighbours_R2[0][x][y][z] = 1; neighbours_R2[1][x][y][z] = 2;
            }
        }

    }
}

int main()
{
    int gm = 3; int g2 =L*L;
    //int nR2[2][L*L][2][3] ={0};
    createVec<4, int> nR2(2, g2, 2, 3, 7);
    determine_neighbours_R2(2, gm, nR2);

    for(int x = 0; x < L*L; x++)
    {
        for(int y=0; y < 2; y++)
        {
            for(int z= 0; z< 3; z++)
            {
                cout << nR2[0][x][y][z] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

    return 0;

}
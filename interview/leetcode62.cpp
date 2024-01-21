#include <iostream>
#include <vector>
using namespace std;

class Solution {
public:

    int numOfPath(int i, int j, int m, int n)
    {
        if (i == m-1)
        {
            return 1;
        }
        if (j == n-1)
        {
            return 1;
        }
        else
        {
            return numOfPath(i+1, j, m, n) + numOfPath(i, j+1, m, n);
        }
    }


    int uniquePaths(int m, int n) {
        vector<vector<int>> mat;
        mat.resize(m, vector<int> (n));
        for (int i = 0; i != m; i++)
        {
            for (int j = 0; j != n; j++)
            {
                mat[i][j] = 0;
            }
        }
        for (int i = 0; i != m; i++)
        {
            mat[i][n-1] = 1;
        }

        for (int j = 0; j != n; j++)
        {
            mat[m-1][j] = 1;
        }
        
        for (int i = m-2; i != -1; i--)
        {
            for(int j = n-2; j != -1; j--)
            {
                mat[i][j] = mat[i+1][j] + mat[i][j+1];
            }
        }

        return mat[0][0];
        
    }
};


int main()
{
    Solution S;
    cout << S.uniquePaths(51,9) << endl;
    return 0;
}
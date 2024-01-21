#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

class Solution {
public:
    int maxSubArray(vector<int>& nums) {
        vector<int> dp;
        dp.resize(nums.size());
        dp[0] = nums[0];
        for (int i = 1; i != nums.size(); i++)
        {
            dp[i] = max(nums[i], dp[i-1] + nums[i]);
        }
        sort(dp.begin(), dp.end());
        return dp[dp.size()-1];
    }
};


int main()
{
    Solution S;
    vector<int> nums = {-2,1,-3,4,-1,2,1,-5,4};
    cout << S.maxSubArray(nums) << endl;

    return 0;

}
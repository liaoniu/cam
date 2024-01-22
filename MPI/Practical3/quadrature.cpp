#include <iostream>
#include <iomanip>
#include <math.h>

#define STEPS 27720

int main(int argc, char *argv[])
{
  double total = 0.0;
  
  for(int i = 0; i < STEPS; i++)
    {
      total += exp( -pow(10.0 * ((i + 0.5) / (STEPS + 1)), 2));
    }

  std::cout << std::fixed << std::setprecision(16);
  std::cout << "Grand total = " << total << ", integral = " << 10.0 * total / STEPS << std::endl;

  return 0;
}

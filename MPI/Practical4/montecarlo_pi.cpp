#include<iostream>
#include<iomanip>
#include<cstdlib>
#include <cmath>

int main(int argc, char **argv){
  int i, count, hit = 0;
  double x, y;

  count = atoi(argv[1]);

  for(i = 0; i < count; i++){
    x=rand() / (double) RAND_MAX;
    y=rand() / (double) RAND_MAX;
    if (x*x + y*y < 1.0)
      {
	hit++;
      }
  }
  std::cout << std::fixed << std::setprecision(16);
  std::cout << "Ratio: " << (double)(hit / (double) count) << std::endl;
  std::cout << "Pi/4 = " << M_PI / 4.0 << "\tError = " << M_PI / 4.0 - (double)(hit / (double) count) << std::endl;
  return 0;
}

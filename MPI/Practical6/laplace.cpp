#include <fstream>
#include <math.h>

#define SIZE 400

using namespace std;

int main(){
  int i, j, n, ht, niter;
  double *g1, *g2, *p;
  ofstream img;
  
  ht=SIZE+2;
  niter=40000;

  g1 = new double[SIZE * ht];
  g2 = new double[SIZE * ht];

  for(i = 0; i < SIZE * ht; i++) {
    g1[i] = 0.0;
  }
  
  for(i = SIZE / 8; i < 7 * SIZE / 8; i++) {
    g1[i] = 1.0;
  }
  
  /* Top */
  for(i = 0; i < SIZE; i++) {
    g2[i] = g1[i];
  }

  /* Bottom */
  for(i = 0; i < SIZE; i++) {
    g2[i + (ht - 1) * SIZE] = g1[i + (ht - 1) * SIZE];
  }
  
  for(n = 0; n < niter; n++){

    for(j = 1; j < ht - 1 ;j++){
      /* Left */
      g2[j * SIZE]=( g1[j * SIZE + 1] + 
		  g1[(j - 1) * SIZE] + g1[(j + 1) * SIZE] ) / 3.0;

      /* Body */
      for(i = 1; i < SIZE - 1; i++){
        g2[j * SIZE + i] = 0.25 * ( g1[j * SIZE + i - 1] + g1[j * SIZE + i + 1] + 
			   g1[(j - 1) * SIZE + i] + g1[(j + 1) * SIZE + i]);
      }
      /* Right */
      g2[j * SIZE + SIZE - 1] = (g1[j * SIZE + SIZE - 2] + 
			 g1[(j - 1) * SIZE + SIZE - 1] + g1[(j + 1) * SIZE + SIZE - 1])/3.0;
    }

    p = g1;g1 = g2;g2 = p;
  }
  img.open("laplace.pam");
  img << "P2\n" << SIZE << " " << SIZE << " 255\n";
  for(i = SIZE; i<(SIZE + 1) * SIZE; i++)
    {
      img << (int)(255*g1[i]) << "\n";
    }

  img.close();
  return 0;
}

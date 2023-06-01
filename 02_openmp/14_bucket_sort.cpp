#include <cstdio>
#include <cstdlib>
#include <vector>
#include <omp.h>

int main() {
  int n = 50;
  int range = 5;
  std::vector<int> key(n);
#pragma omp parallel for
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");

  std::vector<int> bucket(range,0); 
#pragma omp parallel for shared(bucket)
  for (int i=0; i<n; i++)
#pragma omp atomic update
    bucket[key[i]]++;
  std::vector<int> offset(range,0);
  std::vector<int> offsetout(range,0);
#pragma omp parallel
  for(int j=1; j<range; j<<=1){
#pragma omp for
  for(int i=0; i<range; i++)
    offsetout[i] = offset[i];
#pragma omp for
  for(int i=j; i<range;i++)
    offset[i] += offsetout[i-j] + bucket[i-j];
  }

#pragma omp parallel for
  for (int i=0; i<range; i++) {
    int j = offset[i];
    for (; bucket[i]>0; bucket[i]--) {
      key[j++] = i;
    }
  }

#pragma omp parallel
  for (int i=0; i<n; i++) {
#pragma omp single
    printf("%d ",key[i]);
  }
  printf("\n");
}

#include <cstdio>
#include <cstdlib>

int main() {
  const int N=8;
  int a[N], b[N];
  for(int i=0; i<N; i++) {
    a[i] = rand() & 3;
    printf("%*d ",2,a[i]);
  }
  printf("\n");
#pragma omp parallel
  for(int j=1; j<N; j<<=1) {
//    printf("j: %d ", j);
#pragma omp for
    for(int i=0; i<N; i++){
      b[i] = a[i];
      printf("a[%d]=%d \n", i, a[i]);
    }
#pragma omp for
    for(int i=j; i<N; i++){
      a[i] += b[i-j];
      printf("a[%d] = %d b[%d] = %d \n\n", i, a[i], i-j, b[i-j]);
    }
  }
  for(int i=0; i<N; i++)
    printf("%*d ",2,a[i]);
  printf("\n");
}

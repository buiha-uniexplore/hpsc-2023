#include <cstdio>
#include <cstdlib>
#include <vector>

__global__ void initiate(int *bucket){
  bucket[threadIdx.x] = 0;
}

__global__ void increment(int *bucket, int *key){

  int i = key[threadIdx.x];
 
  atomicAdd(&bucket[i], 1);

}

__device__ __managed__ int jj=0;

__global__ void sort(int * bucket, int *key, int &jj){
  i = threadIdx.x;
  extern __shared__ int jjj;
  for (int t = bucket[i]; bucket[i] < 0; bucket[i]--){
    
  }

}


int main() {
  int n = 50;
  int range = 5;
  //std::vector<int> key(n);
  int *key;
  cudaMallocManaged(&key, n*sizeof(int));
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");

  /*
  std::vector<int> bucket(range); 
  for (int i=0; i<range; i++) {
    bucket[i] = 0;
  }
  */

  int *bucket;
  cudaMallocManaged(&bucket, range*sizeof(int));
  initiate<<<1,range>>>(bucket);
  cudaDeviceSynchronize();
  for (int i=0; i<range; i++) {
    printf("%d ", bucket[i]);
  }
  printf("\n");

  /*
  for (int i=0; i<n; i++) {
    bucket[key[i]]++;
  }
  */

  increment<<<1,n>>>(bucket, key);
  cudaDeviceSynchronize();
  for (int i=0; i<range; i++) {
    printf("%d ", bucket[i]);
  }
  printf("\n");

  // for (int i=0, j=0; i<range; i++) {
  //   for (; bucket[i]>0; bucket[i]--) {
  //     key[j++] = i;
  //   }
  // }

  sort<<<1,range>>>(bucket, key, jj);
  cudaDeviceSynchronize();

  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");

  cudaFree(bucket);
}

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

__global__ void offsetcalc(int *a, int *b, int n){
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  for(int j=1; j<n; j<<=1){
    b[i] = a[i];
    __syncthreads();
    a[i] += b[i-j];
    __syncthreads();
  }
  for (int i=0; i<n; i++) {
  //  printf("%d ",a[i]);
  }
}

__global__ void sort(int *bucket, int *key, int *offset){
  int i = threadIdx.x;

  int offsetnum = offset[i-1];

  for (int jj=0; jj<bucket[i]; jj++){
    key[jj+offsetnum] = i;
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
  //  printf("%d ", bucket[i]);
  }
  //printf("\n");

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

  int *offset, *offsetout;
  cudaMallocManaged(&offset, range*sizeof(int));
  cudaMallocManaged(&offsetout, range*sizeof(int));

  for (int i = 0; i<range; i++){
    offset[i] = bucket[i];
  }


  offsetcalc<<<1, range>>>(offset, offsetout, range);
  cudaDeviceSynchronize();
  for (int i=0; i<range; i++) {
    printf("%d ",offsetout[i]);
  }
  printf("\n");

  sort<<<1,range>>>(bucket, key, offset);
  cudaDeviceSynchronize();

  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");

  cudaFree(bucket);
  cudaFree(key);
  cudaFree(offset);
  cudaFree(offsetout);
}

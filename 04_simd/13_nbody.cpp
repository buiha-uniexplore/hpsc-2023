#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <immintrin.h>
#include <pmmintrin.h>

int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
   // printf("%d %g %g\n",i,x[i],y[i]);
  }
  
   __m256 xvec = _mm256_load_ps(x);
   __m256 yvec = _mm256_load_ps(y);
   __m256 mvec = _mm256_load_ps(m);

  float iindex[N];
  for(int i=0; i<N; i++){
    iindex[i] = i;
  }  
  float jindex[N];
  for(int j=0; j<N; j++){
    jindex[j] = j;
  }  



  for(int i=0; i<N; i++) {
        __m256 xivec = _mm256_set1_ps(x[i]);
        __m256 xjvec = _mm256_load_ps(x);
        __m256 rxvec = _mm256_sub_ps(xivec, xjvec);
        
        
        __m256 ivec = _mm256_set1_ps(i);
        __m256 jvec = _mm256_load_ps(jindex);
        __m256 mask = _mm256_cmp_ps(ivec, jvec, _CMP_EQ_OQ);
        rxvec = _mm256_blendv_ps(rxvec, xjvec, mask);

        __m256 yivec = _mm256_set1_ps(y[i]);
        __m256 yjvec = _mm256_load_ps(y);
        __m256 ryvec = _mm256_sub_ps(yivec, yjvec);
        
  
        ryvec = _mm256_blendv_ps(ryvec, yjvec, mask);
      

        float fxi[N];
        float r[N];
	__m256 rxsqr = _mm256_mul_ps(rxvec, rxvec);  
	__m256 rysqr = _mm256_mul_ps(ryvec, ryvec);
        __m256 rsqr = _mm256_add_ps(rxsqr, rysqr);  

        __m256 rrvec = _mm256_rsqrt_ps(rsqr);
        __m256 fxivec = _mm256_mul_ps(rxvec, mvec);
        fxivec = _mm256_mul_ps(fxivec, rrvec);
        fxivec = _mm256_mul_ps(fxivec, rrvec);
        fxivec = _mm256_mul_ps(fxivec, rrvec);
  
        /*
        float testvec[N];
        _mm256_store_ps(testvec, xjvec);
        */
	float rxdb[N];
        float rydb[N];
        float rsqrdb[N];
        float fxidb[N];
	_mm256_store_ps(rxdb, rxvec);
	_mm256_store_ps(rydb, ryvec);
	_mm256_store_ps(rsqrdb, rsqr);
        _mm256_store_ps(fxidb, fxivec);

	for(int k=0; k<N; k++){
          printf("i=%d rx[%d]=%g ry[%d]=%g\n ", i, k, rxdb[k], k, rydb[k]);
        }
	for(int k=0; k<N; k++){
          printf("i=%d rsqr[%d]=%g fxi[%d]=%g\n", i, k, rsqrdb[k], k, fxidb[k]);
        }


   for(int j=0; j<N; j++) {
        


      if(i != j) {
             
        float rx = x[i] - x[j];
        float ry = y[i] - y[j];

        printf("i=%d j=%d rx=%g ry=%g\n",i,j,rx,ry);
        float r = std::sqrt(rx * rx + ry * ry);
        printf("i=%d j=%d r=%g\n",i,j,r);

        fx[i] -= rx * m[j] / (r * r * r);
        printf("i=%d j=%d fx[%d]=%g\n",i,j,i, fx[i]);

        fy[i] -= ry * m[j] / (r * r * r);
      }
    }

    printf("%d %g %g\n",i,fx[i],fy[i]);
  }
}

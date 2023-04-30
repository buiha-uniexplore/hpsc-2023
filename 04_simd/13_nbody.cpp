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
 
  }
  
   __m256 xvec = _mm256_load_ps(x);
   __m256 yvec = _mm256_load_ps(y);
   __m256 mvec = _mm256_load_ps(m);

  
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
        __m256 x0vec = _mm256_setzero_ps();
        __m256 mask = _mm256_cmp_ps(ivec, jvec, _CMP_EQ_OQ);
        rxvec = _mm256_blendv_ps(rxvec, xivec, mask);

        __m256 yivec = _mm256_set1_ps(y[i]);
        __m256 yjvec = _mm256_load_ps(y);
        __m256 ryvec = _mm256_sub_ps(yivec, yjvec);
        __m256 y0vec = _mm256_setzero_ps();
  
        ryvec = _mm256_blendv_ps(ryvec, yjvec, mask);
      
        float fxidb[N];
        float fyidb[N];

        float fxi[N];
        float r[N];
	__m256 rxsqr = _mm256_mul_ps(rxvec, rxvec);  
	__m256 rysqr = _mm256_mul_ps(ryvec, ryvec);
        __m256 rsqr = _mm256_add_ps(rxsqr, rysqr);  
    
        __m256 rrvec = _mm256_rsqrt_ps(rsqr);
    
        __m256 fxivec = _mm256_mul_ps(rxvec, mvec);       
        fxivec = _mm256_blendv_ps(fxivec, x0vec, mask);
        fxivec = _mm256_mul_ps(fxivec, rrvec);    
        fxivec = _mm256_mul_ps(fxivec, rrvec);
	fxivec = _mm256_mul_ps(fxivec, rrvec);
        fxivec = _mm256_sub_ps(x0vec, fxivec);
 
	__m256 fxibvec = _mm256_permute2f128_ps(fxivec, fxivec, 1);
	fxibvec = _mm256_add_ps(fxibvec, fxivec);
        fxibvec = _mm256_hadd_ps(fxibvec, fxibvec);
        fxibvec = _mm256_hadd_ps(fxibvec, fxibvec);
        _mm256_store_ps(fxidb, fxibvec); 
 
        __m256 fyivec = _mm256_mul_ps(ryvec, mvec);       
        fyivec = _mm256_blendv_ps(fyivec, y0vec, mask);
        fyivec = _mm256_mul_ps(fyivec, rrvec);    
        fyivec = _mm256_mul_ps(fyivec, rrvec);
	fyivec = _mm256_mul_ps(fyivec, rrvec);
        fyivec = _mm256_sub_ps(y0vec, fyivec);

	__m256 fyibvec = _mm256_permute2f128_ps(fyivec, fyivec, 1);
	fyibvec = _mm256_add_ps(fyibvec, fyivec);
        fyibvec = _mm256_hadd_ps(fyibvec, fyibvec);
        fyibvec = _mm256_hadd_ps(fyibvec, fyibvec);
        _mm256_store_ps(fyidb, fyibvec); 
      
	float rxdb[N];
        float rydb[N];
        float rsqrdb[N];
        float fxibvecdb[N];
	_mm256_store_ps(rxdb, rxvec);
	_mm256_store_ps(rydb, ryvec);
	_mm256_store_ps(rsqrdb, rsqr);


        printf("%d %g %g \n", i, fxidb[0],fyidb[0]);
        

/*
   for(int j=0; j<N; j++) {
        


      if(i != j) {
             
        float rx = x[i] - x[j];
        float ry = y[i] - y[j];

        //printf("i=%d j=%d rx=%g ry=%g\n",i,j,rx,ry);
        float r = std::sqrt(rx * rx + ry * ry);
        //printf("i=%d j=%d r=%g\n",i,j,r);

        float test1 = rx * m[j];
        fx[i] -= rx * m[j] / (r * r * r);
//        printf("i=%d j=%d rx*m[%d]=%g fx[%d]=%g\n",i,j,j,test1,i, fx[i]);

        fy[i] -= ry * m[j] / (r * r * r);
      }
    }*/

   // printf("%d %g %g\n",i,fx[i],fy[i]);
  }

}

#include <iostream>
#include <chrono>
#include <vector>
using namespace std;
typedef vector<vector<float>> matrix;
#include <immintrin.h>

int main(){
    const int nx = 41;
    const int ny = 41;
    int nt = 500;
    int nit = 50;
    float dx = 2 / (nx - 1);
    float dy = 2 / (ny - 1);
    float dt = .01;
    float rho = 1;
    float nu = .02;
    vector<float> x[nx];
    vector<float> y[ny];
    matrix u(ny, vector<float> (nx));
    matrix v(ny, vector<float> (nx));
    matrix b(ny, vector<float> (nx));
    matrix p(ny, vector<float> (nx));
    matrix un(ny, vector<float> (nx));
    matrix vn(ny, vector<float> (nx));
    matrix bn(ny, vector<float> (nx));
    matrix pn(ny, vector<float> (nx));
    chrono::steady_clock::time_point tic, toc;
    double time;

    // __m256 uvec = _mm256_load_ps(u);
    // __m256 vvec = _mm256_load_ps(v);
    // __m256 bvec = _mm256_load_ps(b);
    // __m256 pvec = _mm256_load_ps(p);
    // __m256 unvec = _mm256_load_ps(un);
    // __m256 vnvec = _mm256_load_ps(vn);
    // __m256 bnvec = _mm256_load_ps(bn);
    // __m256 pnvec = _mm256_load_ps(pn);

    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            u[i][j] = 0;
            v[i][j] = 0;
            p[i][j] = 0;
            b[i][j] = 0;
        }
    }

    
    for (int n=0; n<nt;n++){
        toc = chrono::steady_clock::now();
        for (int j = 1; j<ny-1; j++){
            for (int i = 1; i<nx-1; i++){

                __m256 bvec;
                __m256 uim1 = _mm256_loadu_ps(&u[j][i-1]);
                __m256 uip1 = _mm256_loadu_ps(&u[j][i+1]);
                __m256 ujm1 = _mm256_loadu_ps(&u[j-1][i]);
                __m256 ujp1 = _mm256_loadu_ps(&u[j+1][i]);
                __m256 vim1 = _mm256_loadu_ps(&v[j][i-1]);
                __m256 vip1 = _mm256_loadu_ps(&v[j][i+1]);
                __m256 vjm1 = _mm256_loadu_ps(&v[j-1][i]);
                __m256 vjp1 = _mm256_loadu_ps(&v[j+1][i]);

                bvec = rho * (1 / dt *\
                        ((uip1 - uim1) / (2 * dx) + (vjp1 - vjm1) / (2 * dy)) -\
                        ((uip1 - uim1) / (2 * dx))*((uip1 - uim1) / (2 * dx)) - 2 * ((ujp1 - ujm1) / (2 * dy) *\
                        (vip1 - vim1) / (2 * dx)) - ((vjp1 - vjm1) / (2 * dy))*((vjp1 - vjm1) / (2 * dy)));
            
                _mm256_storeu_ps(&b[i][j], bvec);
            }
        }
        for (int it=0; it<nit; it++){
            for (int j = 1; j<ny-1; j++){
                for (int i = 1; i<nx-1; i++){
                    pn[j][i] = p[j][i];
                }
            }

            for (int j = 1; j<ny-1; j++){
                for (int i = 1; i<nx-1; i++){
                    p[j][i] = (dy*dy * (pn[j][i+1] + pn[j][i-1]) +\
                            dx*dx * (pn[j+1][i] + pn[j-1][i]) -\
                            b[j][i] * dx*dx * dy*dy)\
                            / (2 * (dx*dx + dy*dy));
                }
            }
            
            for(int j = 1; j<ny-1; j++){
                p[j][nx-1] = p[j][nx-2];
                p[j][0] = p[j][1];
            }

            for(int i = 1 ;i<nx-1; i++){
                p[0][i] = p[1][i];
                p[ny-1][i] = p[ny-2][i]; //0?
            }
        }
        for (int j = 1; j<ny-1; j++){
            for (int i = 1; i<nx-1; i++){
                un[j][i] = u[j][i];
                vn[j][i] = v[j][i];
            }
        }
        tic = chrono::steady_clock::now();
        time = chrono::duration<double>(tic-toc).count();
        printf("step =%d : %lf s\n", n, time);
        for (int j = 1; j<ny-1; j++){
            for (int i = 1; i<nx-1; i++){
                __m256 unvec = _mm256_loadu_ps(&un[j][i]);
                __m256 unim1 = _mm256_loadu_ps(&un[j][i-1]);
                __m256 unip1 = _mm256_loadu_ps(&un[j][i+1]);
                __m256 unjm1 = _mm256_loadu_ps(&un[j-1][i]);
                __m256 unjp1 = _mm256_loadu_ps(&un[j+1][i]);
                __m256 pip1 = _mm256_loadu_ps(&p[j][i+1]);
                __m256 pim1 = _mm256_loadu_ps(&p[j][i-1]);
                __m256 uvec;
                uvec = unvec - unvec * dt / dx * (unvec - unim1)\
                                - unvec * dt / dy * (unvec - unjm1)\
                                - dt / (2 * rho * dx) * (pip1 - pim1)\
                                + nu * dt / dx*dx * (unip1 - 2 * unvec + unim1)\
                                + nu * dt / dy*dy * (unjp1 - 2 * unvec + unjm1);

                // _mm256_storeu_ps(&un[j][i], unvec);
                // _mm256_storeu_ps(&un[j][i-1], unim1);
                // _mm256_storeu_ps(&un[j][i+1], unip1);
                // _mm256_storeu_ps(&un[j-1][i], unjm1);
                // _mm256_storeu_ps(&un[j+1][i], unjp1);
                // _mm256_storeu_ps(&p[j][i+1], pip1);
                // _mm256_storeu_ps(&p[j][i-1], pim1);
                _mm256_storeu_ps(&u[j][i], uvec);

                __m256 vnvec = _mm256_loadu_ps(&vn[j][i]);
                __m256 vnim1 = _mm256_loadu_ps(&vn[j][i-1]);
                __m256 vnip1 = _mm256_loadu_ps(&vn[j][i+1]);
                __m256 vnjm1 = _mm256_loadu_ps(&vn[j-1][i]);
                __m256 vnjp1 = _mm256_loadu_ps(&vn[j+1][i]);
                __m256 pjp1 = _mm256_loadu_ps(&p[j+1][i]);
                __m256 pjm1 = _mm256_loadu_ps(&p[j-1][i]);
                __m256 vvec;
                vvec = vnvec - vnvec * dt / dx * (vnvec - vnim1)\
                                - vnvec * dt / dy * (vnvec - vnjm1)\
                                - dt / (2 * rho * dx) * (pjp1 - pjm1)\
                                + nu * dt / dx*dx * (vnip1 - 2 * vnvec + vnim1)\
                                + nu * dt / dy*dy * (vnjp1 - 2 * vnvec + vnjm1);

                // _mm256_storeu_ps(&p[j+1][i], pjp1);
                // _mm256_storeu_ps(&p[j-1][i], pjm1);
                _mm256_storeu_ps(&v[j][i], vvec);
               
            }
        }

        for(int j = 1; j<ny-1; j++){
            u[j][0] = 0;
            p[j][nx-1] = 0;
            v[j][0] = 0;
            v[j][nx-1] = 0;
        }

        for(int i = 1 ;i<nx-1; i++){
            u[0][i] = 0;
            u[ny-1][i] = 1;
            v[0][i] = 0;
            v[ny-1][i] = 0;
        }
    }
}



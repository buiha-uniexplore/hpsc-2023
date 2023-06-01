#include <iostream>
#include <chrono>
#include <vector>
using namespace std;
typedef vector<vector<float>> matrix;
#include <openacc.h>

__global__ void initialization(){

}


int main(){
    const int nx = 41;
    const int ny = 41;
    int nt = 500;
    int nit = 50;
    double dx = 2 / (nx - 1);
    double dy = 2 / (ny - 1);
    double dt = .01;
    double rho = 1;
    double nu = .02;
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

    cudaMallocManaged(&u, nx*nx*sizeof(float));

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
                b[j][i] = rho * (1 / dt *\
                        ((u[j][i+1] - u[j][i-1]) / (2 * dx) + (v[j+1][i] - v[j-1][i]) / (2 * dy)) -\
                        ((u[j][i+1] - u[j][i-1]) / (2 * dx))*((u[j][i+1] - u[j][i-1]) / (2 * dx)) - 2 * ((u[j+1][i] - u[j-1][i]) / (2 * dy) *\
                        (v[j][i+1] - v[j][i-1]) / (2 * dx)) - ((v[j+1][i] - v[j-1][i]) / (2 * dy))*((v[j+1][i] - v[j-1][i]) / (2 * dy)));
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
                u[j][i] = un[j][i] - un[j][i] * dt / dx * (un[j][i] - un[j][i - 1])\
                                - un[j][i] * dt / dy * (un[j][i] - un[j - 1][i])\
                                - dt / (2 * rho * dx) * (p[j][i+1] - p[j][i-1])\
                                + nu * dt / dx*dx * (un[j][i+1] - 2 * un[j][i] + un[j][i-1])\
                                + nu * dt / dy*dy * (un[j+1][i] - 2 * un[j][i] + un[j-1][i]);
                v[j][i] = vn[j][i] - vn[j][i] * dt / dx * (vn[j][i] - vn[j][i - 1])\
                                - vn[j][i] * dt / dy * (vn[j][i] - vn[j - 1][i])\
                                - dt / (2 * rho * dx) * (p[j+1][i] - p[j-1][i])\
                                + nu * dt / dx*dx * (vn[j][i+1] - 2 * vn[j][i] + vn[j][i-1])\
                                + nu * dt / dy*dy * (vn[j+1][i] - 2 * vn[j][i] + vn[j-1][i]);
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
    cudaDeviceSynchronize();
}



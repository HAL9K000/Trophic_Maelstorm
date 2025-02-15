#define DEFINE_CUDA_CONSTANTS // Define the CUDA constant arrays only once through this .cu file in the header file

#include "rietkerk_bjork_basic.h"
#include <cstdio>



// AtomicAdd function analogue for double precision floating point numbers for devices with compute capability < 6.0
#if __CUDA_ARCH__ < 600
__device__ double atomicAdd_Double(double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}
#endif

cudaError_t copyToDeviceConstantMemory_AdvTerms(const int* sig_D_ScaledBounds, const int2* sigma_vD_ScBounds, 
        const double* mu_vel_prefactor,  const int* size_gauss_D, const int* size_gauss_VXY)
{
    cudaError_t err;

    err = cudaMemcpyToSymbol(d_sigD_ScaleBounds, sig_D_ScaledBounds, CuSpNV * sizeof(int), 0, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
    printf("cudaMemcpyToSymbol failed for d_sigD_ScaleBounds: %s\n", cudaGetErrorString(err)); return err;
    }
    else{
        printf("d_sigD_ScaleBounds copied to constant memory with values:\n");
        for(int i=0; i<CuSpNV; i++)
            printf("%d\t", sig_D_ScaledBounds[i]);
        }

    err = cudaMemcpyToSymbol(d_sigvD_ScaleBounds, sigma_vD_ScBounds, CuSpNV * sizeof(int2), 0, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
    printf("cudaMemcpyToSymbol failed for d_sigvD_ScaleBounds: %s\n", cudaGetErrorString(err)); return err;
    }
    else{
        printf("\nd_sigvD_ScaleBounds copied to constant memory with values:\n");
        for(int i=0; i<CuSpNV; i++)
            printf("[ %d , %d , %d ]\t", i, sigma_vD_ScBounds[i].x, sigma_vD_ScBounds[i].y);
        
    }

    err = cudaMemcpyToSymbol(d_mu_vel_prefactor, mu_vel_prefactor, CuSpNV * sizeof(double), 0, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
    printf("cudaMemcpyToSymbol failed for d_mu_vel_prefactor: %s\n", cudaGetErrorString(err));
    return err;
    }
    else{
            printf("\nd_mu_vel_prefactor copied to constant memory with values:\n");
            for(int i=0; i<CuSpNV; i++)
                printf("%f\t", mu_vel_prefactor[i]);
        }

    err = cudaMemcpyToSymbol(d_size_gauss_D, size_gauss_D, CuSpNV * sizeof(int), 0, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
    printf("cudaMemcpyToSymbol failed for d_size_gauss_D: %s\n", cudaGetErrorString(err));
    return err;
    }
    else{
            printf("\nd_size_gauss_D copied to constant memory with values:\n");
            for(int i=0; i<CuSpNV; i++)
                printf("%d\t", size_gauss_D[i]);
        }

    err = cudaMemcpyToSymbol(d_size_gauss_VXY, size_gauss_VXY, CuSpNV * sizeof(int), 0, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        printf("cudaMemcpyToSymbol failed for d_size_gauss_VXY: %s\n", cudaGetErrorString(err)); return err;
    }
    else{
            printf("\nd_size_gauss_VXY copied to constant memory with values:\n");
            for(int i=0; i<CuSpNV; i++)
                printf("%d\t", size_gauss_VXY[i]);
        }

    return cudaSuccess;
}


// CUDA kernel to update the concentration field
__global__ void updateMovement( double* Rho_t, double* Rho_tsar, double* gamma, std::pair<double, double> *v_eff,
    double* gaussian_stencil_D,  double* gaussian_stencil_VXY, int L) 
{
    int L2 = L*L;
    // Thread indices
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = gridDim.x * blockDim.x;

    for(int s=0; s< CuSpNV; s++)
	{
        // IMPORTANT: s = 0 is GRAZER, s = 1 is PREDATOR
        // s has been shifted by -1 compared to CPU code to save memory.

        // Iterate over lattice sites in a grid-stride loop
        for (int i = idx; i < L*L; i += stride) 
        {
            int c_i = int(i/L);  // Row index
            int c_j = i%L;  // Column index

            if (Rho_t[s*L2 + i] < epsilon) 
                continue;

            // Diffusion stencil bounds
            int min_p = c_i - d_sigD_ScaleBounds[s]; int max_p = c_i + d_sigD_ScaleBounds[s];
            int min_q = c_j - d_sigD_ScaleBounds[s]; int max_q = c_j + d_sigD_ScaleBounds[s];

            // Advection mean coordinates
            int mu_x = (int)(round(c_i + v_eff[s*L2 + i].first*(d_mu_vel_prefactor[s])) + L) % L; // Assuming initial condition v_x(t) = v_y(t) =0
            
            int mu_y = (int)(round(c_j + v_eff[s*L2 + i].second*(d_mu_vel_prefactor[s])) + L) % L;

            // Advection stencil bounds
            /**  RECALL: d_sigvD_Scalebound is __constant__ int2 array of size CuSpNV
             and corresponds to sig_vD_ScaledBounds[s] from CPU code **/
            int min_p_vx = mu_x - d_sigvD_ScaleBounds[s].x; int max_p_vx = mu_x + d_sigvD_ScaleBounds[s].x;
            int min_q_vy = mu_y - d_sigvD_ScaleBounds[s].y; int max_q_vy = mu_y + d_sigvD_ScaleBounds[s].y;
            
            // Update diffusion term
            for (int p = min_p; p <= max_p; ++p) 
            {
                // Apply periodic boundary conditions
				int dx = std::abs(p - c_i); // Calculate distance from mean //int x = (p + L) % L; 
                for (int q = min_q; q <= max_q; ++q) 
                {
                    //int y = (q + L) % L;
                    int dy = abs(q - c_j);
                    //C_t[x * Nx + y] += (1 - a1) * (1 - a1) * C[i * Nx + j] * gaussian_stencil_D[dx] * gaussian_stencil_D[dy];
                    atomicAdd(&Rho_tsar[s*L2 + i], (gamma[s*L2 + i]) * (gamma[s*L2 + i]) * Rho_t[s*L2 + i] 
                    *(gaussian_stencil_D[d_size_gauss_D[s] + dx] *gaussian_stencil_D[d_size_gauss_D[s] + dy]));
                }
                // Next update d_C_t[] using joint probability of advection and diffusion terms
                for(int q = min_q_vy; q <= max_q_vy; ++q)
                {
                    //int y = (q + L) % L;
                    int dy = std::abs(q - mu_y);
                    //C_t[x][y] += (1-a1)*(a1)*C[i][j]*(gaussian_stencil_D[dx]*gaussian_stencil_vy[dy]);
                    atomicAdd(&Rho_tsar[s*L2 + i], (gamma[s*L2 + i]) *(1 - gamma[s*L2 + i])* Rho_t[s*L2 + i] *
                    (gaussian_stencil_D[d_size_gauss_D[s]+dx]*gaussian_stencil_VXY[d_size_gauss_VXY[s]+dy]));
                }
            }

            // Update advection term
            for (int p = min_p_vx; p <= max_p_vx; ++p) 
            {
                //int x = (p + L) % L;
                int dx = abs(p - mu_x);
                // Probability of advection in x-axis and diffusion in y-direction
                for (int q = min_q; q <= max_q; ++q) 
                {
                    //int y = (q + L) % L;
                    int dy = std::abs(q - c_j);
                    //C_t[x][y] += (a1)*(1-a1)*C[i][j]*(gaussian_stencil_vx[dx]*gaussian_stencil_D[dy]);
                    atomicAdd(&Rho_tsar[s*L2 + i], (1 - gamma[s*L2 + i])*(gamma[s*L2 + i])*Rho_t[s*L2 + i] 
                    *(gaussian_stencil_VXY[d_size_gauss_VXY[s]+dx]*gaussian_stencil_D[d_size_gauss_D[s]+dy]));
                }

                // Update pure advection term joint probability
                for (int q = min_q_vy; q <= max_q_vy; ++q) 
                {
                    //int y = (q + L) % L;
                    int dy = abs(q - mu_y);
                    //C_t[x * Nx + y] += a1 * a1 * C[i * Nx + j] * gaussian_stencil_vx[dx] * gaussian_stencil_vy[dy];
                    atomicAdd(&Rho_tsar[s*L2 + i], (1-gamma[s*L2 + i])*(1-gamma[s*L2 + i])* Rho_t[s*L2 + i] 
                    *(gaussian_stencil_VXY[d_size_gauss_VXY[s]+dx]*gaussian_stencil_VXY[d_size_gauss_VXY[s]+dy]));
                }
            }
        }
    }
}

// Function that reports all the values stored in __constant__ memory
__global__ void reportConstantMemory_AdvTerms()
{
    printf("Printing values stored in __constant__ memory \n");
    printf("(Sp, d_size_gauss_D): \t");
    for(int i=0; i<CuSpNV; i++){
        printf("( %d , %d )\t", i, d_size_gauss_D[i]);
    }
    printf("\n");

    printf("d_size_gauss_VXY: ");
    for(int i=0; i<CuSpNV; i++){
        printf("%d\t", d_size_gauss_VXY[i]);
    }
    printf("\n");

    printf("d_mu_vel_prefactor: ");
    for(int i=0; i<CuSpNV; i++){
        printf("%f\t", d_mu_vel_prefactor[i]);
    }
    printf("\n");

    printf("d_sigD_ScaleBounds: ");
    for(int i=0; i<CuSpNV; i++){
        printf("%d\t", d_sigD_ScaleBounds[i]);
    }
    printf("\n");

    printf("(Sp, d_sigvD_ScaleBounds.X, d_sigvD_ScaleBounds.Y): ");
    for(int i=0; i<CuSpNV; i++){
        printf("[ %d , %d , %d ]\t", i, d_sigvD_ScaleBounds[i].x, d_sigvD_ScaleBounds[i].y);
    }
    printf("\n");


}

void reportConstantMemory()
{
    reportConstantMemory_AdvTerms<<<1,1>>>();
}


// Function to launch the FKE- Advection Diffusion CUDA kernel
void launchAdvDiffKernel_MultiSp(double* Rho_t, double* Rho_tsar, double* gamma, std::pair<double, double> *v_eff,
    double* gaussian_stencil_D,  double* gaussian_stencil_VXY, int L) 
{
    // Calculate grid and block sizes
    int blockSize = 256;
    int numBlocks = (L * L + blockSize - 1) / blockSize;

    // Launch the CUDA kernel
    updateMovement<<<numBlocks, blockSize>>>(Rho_t, Rho_tsar, gamma, v_eff, gaussian_stencil_D,
        gaussian_stencil_VXY, L);
}
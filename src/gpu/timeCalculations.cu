

__global__ void buildUstar(float *Ustar, float *U, float *R, float *ShearSource, float dt, int m, int n)
{
	// Calculate the row and column of the thread within the thread block
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	
	// First check if the thread is operating on a cell inside of the block's two cell deep ghost cells
	if (col > 1 && row > 1 && col < n-2 && row < m-2)
	{
		// Calculate the index of this cell in the U, Ustar, and R matrices
		int i = row*n*3 + col*3;
		
		// Calculate the index of this cell's shear source value
		int shearIndex = row*n + col;
		
		// Build Ustar
		Ustar[i] = U[i] + dt * R[i];
		Ustar[i + 1] = (U[i + 1] + dt * R[i + 1]) / (1.0f + dt * ShearSource[shearIndex]);
		Ustar[i + 2] = (U[i + 2] + dt * R[i + 2]) / (1.0f + dt * ShearSource[shearIndex]);
	}
}

__global__ void buildUnext(float *Unext, float *U, float *Ustar, float *Rstar, float *ShearSourceStar, float dt, int m, int n)
{
	// Calculate the row and column of the thread within the thread block
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	
	// First check if the thread is operating on a cell inside of the block's two cell deep ghost cells
	if (col > 1 && row > 1 && col < n-2 && row < m-2)
	{
		// Calculate the index of this cell in the U, Unext, Ustar, and R matrices
		int i = row*n*3 + col*3;
		
		// Calculate the index of this cell's shear source value
		int shearIndex = row*n + col;
		
		// Build Unext
		Unext[i] = 0.5f * U[i] + 0.5f * (Ustar[i] + dt * Rstar[i]);
		Unext[i + 1] = (0.5f * U[i + 1] + 0.5f * (Ustar[i + 1] + dt * Rstar[i + 1])) / (1.0f + 0.5f * dt * ShearSourceStar[shearIndex]);
		Unext[i + 2] = (0.5f * U[i + 2] + 0.5f * (Ustar[i + 2] + dt * Rstar[i + 2])) / (1.0f + 0.5f * dt * ShearSourceStar[shearIndex]);
	}
}

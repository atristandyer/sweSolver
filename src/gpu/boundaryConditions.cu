

__global__ void applyWallBoundaries(float *U, int m, int n)
{
	// Wall boundaries are applied by mirroring the two cells nearest the
	// boundary and changing the sign of the normal discharge component
	
	// Calculate the row and column of the thread within the thread block
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	
	// Calculate this cell's location in the U matrix
	int i = row*n*3 + col*3;
	int rIndex;
	
	// The bottom row
	if (row < 2 && col < n)
	{
		// Calculate the index of the cell we will get reflected values from
		rIndex = (3-row)*n*3 + col*3;
		
		// Reflect the values and reverse the sign of the normal velocity
		U[i] = U[rIndex];
		U[i + 2] = -U[rIndex + 2];
	}
	
	// The top row
	else if (row > m-3 && row < m && col < n)
	{
		// Calculate the index of the cell we will get reflected values from
		rIndex = (2*m-5-row)*n*3 + col*3;
		
		// Reflect the values and reverse the sign of the normal velocity
		U[i] = U[rIndex];
		U[i + 2] = -U[rIndex + 2];
	}
	
	// The left column
	if (col < 2 && row < m)
	{
		// Calculate the index of the cell we will get reflected values from
		rIndex = row*n*3 + (3-col)*3;
		
		// Reflect the values and reverse the sign of the normal velocity
		U[i] = U[rIndex];
		U[i + 1] = -U[rIndex + 1];
	}
	
	// The right column
	else if (col > n-3 && col < n && row < m)
	{
		// Calculate the index of the cell we will get reflected values from
		rIndex = row*n*3 + (2*n-5-col)*3;
		
		// Reflect the values and reverse the sign of the normal velocity
		U[i] = U[rIndex];
		U[i + 1] = -U[rIndex + 1];
	}
}

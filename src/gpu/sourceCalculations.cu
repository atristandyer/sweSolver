__global__ void bedSlopeSourceSolver(float *BedSlopeSource, float *UIntPts, float *BottomIntPts, int m, int n, float dx, float dy)
{
	// Calculate the row and column of the thread within the thread block
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	
	// First check if the thread is operating on a cell inside of the block's one cell deep ghost cells
	if (col > 0 && row > 0 && col < n-1 && row < m-1)
	{
		// Calculate index of this cell in the BedSlopeSource matrix
		int slopeIndex = row*n*2 + col*2;
		
		// Calculate index of this cell in the UIntPts matrix
		int uIndex = row*n*4*3 + col*4*3;
		
		// Calculate index of this cell in the BottomIntPts matrix
		int bottomIndex = row*(n+1)*2 + col*2;
		
		//// Calculate the water depth at the center of the cell
		// Water depth 
	}
}

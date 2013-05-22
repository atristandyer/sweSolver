__global__ void bedSlopeSourceSolver(float *BedSlopeSource, float *U, float *BottomIntPts, int m, int n, float dx, float dy)
{
	// Calculate the row and column of the thread within the thread block
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	
	// First check if the thread is operating on a cell inside of the block's one cell deep ghost cells
	if (col > 0 && row > 0 && col < n-1 && row < m-1)
	{
		// Calculate index of this cell in the BedSlopeSource matrix
		int slopeIndex = row*n*2 + col*2;
		
		// Calculate index of this cell in the U matrix
		int uIndex = row*n*3 + col*3;
		
		// Calculate index of this cell in the BottomIntPts matrix
		int bottomIndex = row*(n+1)*2 + col*2;
		
		// Calculate the water depth at the center of the cell
		// Note: h = w - B;
		//	 w is the water column height at the cell center measured from z = 0;
		//	 B is the terrain elevation measured from z = 0, and is calculated as the average
		//	   elevation of the center of two opposing edges (in this case, left and right)
		float hCenter =  U[uIndex] - (BottomIntPts[bottomIndex+1] + BottomIntPts[bottomIndex + 1*2 + 1])/2.0f;
		
		// Calculate the slope of the terrain in both x- and y-directions
		float slopeX = (BottomIntPts[bottomIndex + 3] - BottomIntPts[bottomIndex + 1]) / dx;
		float slopeY = (BottomIntPts[bottomIndex + (n+1)*2] - BottomIntPts[bottomIndex]) / dy;
		
		// Calculate the bed slope source terms
		BedSlopeSource[slopeIndex] = -9.81f * slopeX * hCenter;
		BedSlopeSource[slopeIndex+1] = -9.81f * slopeY * hCenter;		 
	}
}

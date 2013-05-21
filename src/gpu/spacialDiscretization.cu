___device___ float minmod(float a, float b, float c)
{
	float ab = fminf(fabsf(a), fabs(b)) * (copysignf(1.0f, a) + copysignf(1.0f, b)) * 0.5f;
	return fminf(fabsf(ab, fabsf(c)) * (copysignf(1.0f, ab) + copysignf(1.0f, c)) * 0.5f;
}
    
    
__global__ void ReconstructFreeSurface(float *U, float *BottomIntPts, float *UIntPts, float *huvIntPts, int m, int n, float dx, float dy)
{
	// Calculate the row and column of the thread within the thread block
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	
	// Calculate the index of the up, down, left, right terrain elevation values in the BottomIntPts vector
	// 	Note: We only need elevation values from two point on opposite sides of the cell
	//	      because the bottom is piecewise-bilinear, which guarantees that the average
	// 	      of both sets of points on opposing sides of the cell is equal to the value
	//	      at the center of the cell.
	int upTerrainIndex = (row+1)*(n+1)*2 + col*2;
	int downTerrainIndex = row*(n+1)*2 + col*2;
	int leftTerrainIndex = row*(n+1)*2 + col*2 + 1;
	int rightTerrainIndex = row*(n+1)*2 + (col+1)*2 + 1;
	
	// Before performing calculations for a timestep, we need to reset any residual values left over
	// from previous calculations in case a cell goes from wet to dry. In the case of a cell going
	// from wet to dry, the values of huvIntPts and UIntPts for that cell would otherwise remain at the 
	// values from the previous timestep when the cell was wet, causing errors.
	// 	huvIntPts: set h = u = v = 0
	//	UIntPts: set w = ground elevation, u = v = 0
	int cellIntPtIndex = row*n*4*2 + col*4*2;
	for (int i=0; i<4; i++)
	{
		for (int j=0; i<2; j++)
		{
			huvIntPts[cellIntPtIndex + i*2 + j] = 0.0f;
			if (j > 0)
			{
				UIntPts[cellIntPtIndex + i*2 + j] = 0.0f;
			}
		}
	}
	UIntPts[cellIntPtIndex + 0*2] = BottomIntPts[upTerrainIndex];
	UIntPts[cellIntPtIndex + 1*2] = BottomIntPts[downTerrainIndex];
	UIntPts[cellIntPtIndex + 2*2] = BottomIntPts[rightTerrainIndex];
	UIntPts[cellIntPtIndex + 3*2] = BottomIntPts[leftTerrainIndex];

	// First check if the thread is operating on a cell inside of the block's one cell deep ghost cells
	if (col > 0 && row > 0 && col < n-1 && row < m-1)
	{
	
		// Calculate the index of the [w, hu, hv] values in the U vector
		int wIndex = row*n*3 + col*3 + 0;
		int huIndex = row*n*3 + col*3 + 1;
		int hvIndex = row*n*3 + col*3 + 2;
	
		// Now check if the cell has water (evaluates to true if it does)
		if (U[wIndex] > (BottomIntPts[leftTerrainIndex] + BottomIntPts[rightTerrainIndex])/2.0f) 
		{
			float N[2], S[2], E[2], W[2]; // These are the [w, hu] and [w, hv] vectors
			float north[2], south[2], east[2], west[2]; // These are the [h, u] and [h, v] vectors
			
			float forward, central, backward, slope;
			float Kappa = 0.01f * fmaxf(1.0f, fminf(dx, dy));
			float sqrt2 = sqrtf(2.0f);
			float g = 9.81f;
			
			// Calculate indices of cells surrounding the current cell in the U vector
			int upIndex = (row+1)*n*3 + col*3;
			int downIndex = (row-1)*n*3 + col*3;
			int leftIndex = row*n*3 + (col-1)*3;
			int rightIndex = row*n*3 + (col+1)*3;
			
			// Reconstruct the free surface so that it is sloped based on the free surface height of adjacent cells
			for (int i=0; i<2; i++)
			{
				// This value is used to make N/S values used [w, hv] instead of [w, hu]
				int extraIndex = 0;
				if (i == 1)
				{
					extraIndex += 1;
				}
				
				// North and South
				forward = (U[upIndex+i+extraIndex] - U[wIndex+i+extraIndex])/dy;
				central = (U[upIndex+i+extraIndex] - U[downIndex+i+extraIndex])/(2.0f*dy);
				backward = (U[wIndex+i+extraIndex] - U[downIndex+i+extraIndex])/dy;
				slope = minmod(1.3f*forward, central, 1.3f*backward);
				
				N[i] = U[wIndex+i+extraIndex] + (dy/2.0f)*slope;
				S[i] = U[wIndex+i+extraIndex] - (dy/2.0f)*slope;
				
				// East and West
				forward = (U[rightIndex+i] - U[wIndex+i])/dx;
				central = (U[rightIndex+i] - U[leftIndex+i)/(2.0f*dx);
				backward = (U[wIndex+i] - U[leftIndex+i])/dx;
				slope = minmod(1.3f*forward, central, 1.3f*backward);
				
				E[i] = U[wIndex+i] + (dx/2.0f)*slope;
				W[i] = U[wIndex+i] - (dx/2.0f)*slope;
			}
			
			// Check the water free surface at the cell interfaces for positivity and correct if necessary
			if (N[0] < BottomIntPts[upTerrainIndex])
			{
				N[0] = BottomIntPts[upTerrainIndex];
				S[0] = 2*U[wIndex] - BottomIntPts[upTerrainIndex];
			}
			else if (S[0] < BottomIntPts[downTerrainIndex])
			{
				S[0] = BottomIntPts[downTerrainIndex];
				N[0] = 2*U[wIndex] - BottomIntPts[downTerrainIndex];
			}
			if (E[0] < BottomIntPts[rightTerrainIndex])
			{
				E[0] = BottomIntPts[rightTerrainIndex];
				W[0] = 2*U[wIndex] - BottomIntPts[rightTerrainIndex];
			}
			else if (W[0] < BottomIntPts[leftTerrainIndex])
			{
				W[0] = BottomIntPts[leftTerrainIndex];
				E[0] = 2*U[wIndex] - BottomIntPts[leftTerrainIndex];
			}
			
			// Calculate the values of h, u, and v using the damping calculations for u and v
			north[0] = N[0] - BottomIntPts[upTerrainIndex];
			south[0] = S[0] - BottomIntPts[downTerrainIndex];
			east[0] = E[0] - BottomIntPts[rightTerrainIndex];
			west[0] = W[0] - BottomIntPts[leftTerrainIndex];
			
			north[1] = (sqrt2 * north[0] * N[1]) / sqrtf(powf(north[0], 4.0f) + fmaxf(powf(north[0], 4.0f), Kappa));
			south[1] = (sqrt2 * south[0] * S[1]) / sqrtf(powf(south[0], 4.0f) + fmaxf(powf(south[0], 4.0f), Kappa));
			east[1] = (sqrt2 * east[0] * E[1]) / sqrtf(powf(east[0], 4.0f) + fmaxf(powf(east[0], 4.0f), Kappa));
			west[1] = (sqrt2 * west[0] * W[1]) / sqrtf(powf(west[0], 4.0f) + fmaxf(powf(west[0], 4.0f), Kappa));
			
			// Update the values of hu and hv based on new values of u and v
			N[1] = north[0] * north[1];
			S[1] = south[0] * south[1];
			E[1] = east[0] * east[1];
			W[1] = west[0] * west[1];
			
			// Put the calculated interface values into global memory
			for (int i=0; i<2; i++)
			{
				UIntPts[cellIntPtIndex + 0*2 + i] = N[i];
				UIntPts[cellIntPtIndex + 1*2 + i] = S[i];
				UIntPts[cellIntPtIndex + 2*2 + i] = E[i];
				UIntPts[cellIntPtIndex + 3*2 + i] = W[i];
				
				huvIntPts[cellIntPtIndex + 0*2 + i] = north[i];
				huvIntPts[cellIntPtIndex + 1*2 + i] = south[i];
				huvIntPts[cellIntPtIndex + 2*2 + i] = east[i];
				huvIntPts[cellIntPtIndex + 3*2 + i] = west[i];
			}
			
			// End the kernel here. The values of the conserved variable [h, u, v] need to be stored
			// in global memory for the entire domain before we can begin calculating propagation speeds.
		}
	}
}


__global__ void CalculatePropSpeeds(float *UIntPts, float *huvIntPts, float *propSpeeds, int m, int n)
{
	// Calculate the row and column of the thread within the thread block
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	
	// Constants
	float g = 9.81f;
	
	// Get rid of any residual values in the propSpeeds matrix in case
	// a cell has gone from wet to dry in the previous timestep
	for (int i=0; i<4; i++)
	{
		propSpeeds[row*n*4 + col*4 + i] = 0.0f;
	}
	
	// First check if the thread is operating on a cell inside of the block's one cell deep ghost cells
	if (col > 0 && row > 0 && col < n-1 && row < m-1)
	{
		
		// Make sure cell is wet by making sure at least one value of h is greater than 0
		float hNorth = huvIntPts[row*n*4*2 + col*4*2 + 0*2];		// value of h at the north interface of the current cell
		float hSouth = huvIntPts[(row+1)*n*4*2 + col*4*2 + 1*2];	// value of h at the south interface of the cell above
		float hEast = huvIntPts[row*n*4*2 + col*4*2 + 2*2];		// value of h at the east interface of the current cell
		float hWest = huvIntPts[row*n*4*2 + (col+1)*4*2 + 3*2];		// value of h at the west interface of the cell to the right
		if (hNorth > 0.0f || hSouth > 0.0f || hEast > 0.0f || hWest > 0.0f)
		{
			// Get the rest of the values needed from huvIntPts
			float vNorth = huvIntPts[row*n*4*2 + col*4*2 + 0*2 + 1];	// value of v at the north interface of the current cell
			float vSouth = huvIntPts[(row+1)*n*4*2 + col*4*2 + 1*2 + 1];	// value of v at the south interface of the cell above
			float uEast = huvIntPts[row*n*4*2 + col*4*2 + 2*2 + 1];		// value of u at the east interface of the current cell
			float uWest = huvIntPts[row*n*4*2 + (col+1)*4*2 + 3*2 + 1];	// value of u at the west interface of the cell to the right
		
			// Each cell in propSpeeds contains four values [N, S, E, W]
			// Each thread will calculate the N, E values of it's own cell, the S of the cell above, and the W of
			// the cell to the right.
			int N = row*n*4 + col*4 + 0;	// North value of this cell
			int S = (row+1)*n*4 + col*4 + 1;// South value of cell above
			int E = row*n*4 + col*4 + 2;	// East value of this cell
			int W = row*n*4 + (col+1)*4 + 3;// West value of cell to the right
		
			// Calculate north propagation speed of the current cell
			propSpeeds[N] = fminf(fminf(vNorth - sqrtf(g*hNorth), vSouth - sqrtf(g*hSouth)), 0.0f);
		
			// Calculate south propagation speed of the cell above
			propSpeeds[S] = fmaxf(fmaxf(vNorth + sqrtf(g*hNorth), vSouth + sqrtf(g*hSouth)), 0.0f);
		
			// Calculate east propagation speed of the current cell
			propSpeeds[E] = fminf(fminf(uEast - sqrtf(g*hEast), uWest - sqrtf(g*hWest)), 0.0f);
		
			// Calculate west propagation speed of the cell to the right
			propSpeeds[W] = fmaxf(fmaxf(uEast + sqrtf(g*hEast), uWest + sqrtf(g*hWest)), 0.0f);
		}
	}
}

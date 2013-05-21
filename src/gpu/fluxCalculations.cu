
// This function calculates the flux terms in the x-direction
__device__ void F(float F_vec[3], float U_vec[3], float bottomElevation)
{
	float h = U_vec[0] - bottomElevation;
	if (h <= 0.0f)
	{
		F_vec[0] = 0.0f;
		F_vec[1] = 0.0f;
		F_vec[2] = 0.0f;
	} else {
		F_vec[0] = U_vec[1];
		F_vec[1] = (powf(U_vec[1], 2.0f)/h) + 0.5f * 9.81f * powf(h, 2.0f);
		F_vec[2] = (U_vec[1] * U_vec[2]) / h;
	}
}

// This function calculates the flux terms in the y-direction
__device__ void G(float G_vec[3], float U_vec[3], float bottomElevation)
{
	float h = U_vec[0] - bottomElevation;
	if (h <= 0.0f)
	{
		G_vec[0] = 0.0f;
		G_vec[1] = 0.0f;
		G_vec[2] = 0.0f;
	} else {
		G_vec[0] = U_vec[2];
		G_vec[1] = (U_vec[1] * U_vec[2]) / h;
		G_vec[2] = (powf(U_vec[2], 2.0f)/h) + 0.5f * 9.81f * powf(h, 2.0f);
	}
}


__global__ void FluxSolver(float *Fluxes, float *UIntPts, float *BottomIntPts, float *propSpeeds, int m, int n)
{
	// Calculate the row and column of the thread within the thread block
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	
	// First check if the thread is operating on a cell inside of the block's one cell deep ghost cells
	if (col > 0 && row > 0 && col < n-1 && row < m-1)
	{
		//// Retrieve data and calculate indices
		// Calculate the starting index of this thread's cell in the Fluxes array
		int fluxCellIndex = row*n*2*3 + col*2*3;
	
		// Calculate the starting index of each interface point in UIntPts
		int northIndex = row*n*4*3 + col*4*3 + 0*3;	// The north index of U for the current cell
		int southIndex = (row+1)*n*4*3 + col*4*3 + 1*3;	// The south index of U for the cell above
		int eastIndex = row*n*4*3 + col*4*3 + 2*3;	// The east index of U for the current cell
		int westIndex = row*n*4*3 + (col+1)*4*3 + 3*3;	// The west index of U for the current cell
	
		// Get the two bottom elevations we need
		float northBottomElevation = BottomIntPts[(row+1)*(n+1)*2 + col*2 + 0];
		float eastBottomElevation = BottomIntPts[row*(n+1)*2 + (col+1)*2 + 1];


		//// East flux calculations
		// Get the propagation speeds to the left and right of the right side interface
		float aLeft = propSpeeds[(row+1)*n*4 + col*4 + 3];
		float aRight = propSpeeds[row*n*4 + col*4 + 2];
		
		// At least one propagation speed needs to be nonzero for flux calculation
		if (aLeft != 0.0f || aRight != 0.0f)
		{
			// Get the U vectors to the left and right of the right side interface
			float Uleft[3] = {UIntPts[westIndex], UIntPts[westIndex + 1], UIntPts[westIndex + 2]};
			float Uright[3] = {UIntPts[eastIndex], UIntPts[eastIndex + 1], UIntPts[eastIndex + 2]};
			
			// Create the F vectors and calculate them
			float Fright[3], Fleft[3];
			F(Fright, Uright, eastBottomElevation);
			F(Fleft, Uleft, eastBottomElevation);
			
			// Calculate the fluxes
			for (int i=0; i<3; i++)
			{
				Fluxes[fluxCellIndex + 1*3 + i] = ((aRight*Fright[i] - aLeft*Fleft[i]) / (aRight - aLeft)) + ((aRight*aLeft)/(aRight-aLeft))*(Uleft[i] - Uright[i]);
			}
		} else {
			for (int i=0; i<3; i++)
			{
				Fluxes[fluxCellIndex + 1*3 + i] = 0.0f;
			}
		}
		
		
		//// North flux calculations
	}
}




























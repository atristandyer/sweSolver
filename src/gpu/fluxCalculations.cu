
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
		float aLeft = propSpeeds[row*n*4 + (col+1)*4 + 3];
		float aRight = propSpeeds[row*n*4 + col*4 + 2];
		
		// At least one propagation speed needs to be nonzero for flux calculation
		if (aLeft != 0.0f || aRight != 0.0f)
		{
			// Get the U vectors to the left and right of the right side interface
			float Uleft[3] = {UIntPts[eastIndex], UIntPts[eastIndex + 1], UIntPts[eastIndex + 2]}; // U just to the left of the interface (ie. East value of current cell)
			float Uright[3] = {UIntPts[westIndex], UIntPts[westIndex + 1], UIntPts[westIndex + 2]};// U just to the right of the interface (ie. West value of cell to the right)
			
			// Create the F vectors and calculate them
			float Fright[3], Fleft[3];
			F(Fright, Uright, eastBottomElevation); // Calculate F based on values of U just to the right of the cell interface
			F(Fleft, Uleft, eastBottomElevation);   // Calculate F based on values of U just to the left of the cell interface
			
			// Calculate the flux across the right side interface
			for (int i=0; i<3; i++)
			{
				Fluxes[fluxCellIndex + 1*3 + i] = ((aRight*Fleft[i] - aLeft*Fright[i]) / (aRight - aLeft)) + ((aRight*aLeft)/(aRight-aLeft))*(Uright[i] - Uleft[i]);
			}
		} else {
			for (int i=0; i<3; i++)
			{
				// If neither this cell nor the cell to the right has a propagation speed (ie. is dry), there is no flux
				Fluxes[fluxCellIndex + 1*3 + i] = 0.0f;
			}
		}
		
		
		//// North flux calculations
		// Get the propagation speeds above and below the upper interface
		float bUp = propSpeeds[(row+1)*n*4 + col*4 + 1];
		float bDown = propSpeeds[row*n*4 + col*4 + 0];
		
		// At least one propagation speed needs to be nonzero for flux calculation
		if (bUp != 0.0f || bDown != 0.0f)
		{
			// Get the U vectors above and below the upper interface
			float Uup[3] = {UIntPts[southIndex], UIntPts[southIndex + 1], UIntPts[southIndex + 2]};  // U just above the interface (ie. South value of the current cell)
			float Udown[3] = {UIntPts[northIndex], UIntPts[northIndex + 1], UIntPts[northIndex + 2]};// U just below the interface (ie. North value of the current cell)
			
			// Create the G vectors and calculate them
			float Gup[3], Gdown[3];
			G(Gup, Uup, northBottomElevation);     // Calculate G based on the values just above the interface
			G(Gdown, Udown, northBottomElevation); // Calculate G based on the values just below the interface
			
			// Calculate the flux across the upper interface
			for (int i=0; i<3; i++)
			{
				Fluxes[fluxCellIndex + i] = ((bUp*Gdown[i] - bDown*Gup[i]) / (bUp - bDown)) + ((bUp*bDown)/(bUp-bDown))*(Uup[i] - Udown[i]);
			}
		} else {
			for (int i=0; i<3; i++)
			{
				// If neither this cell nor the cell above has a propagation speed (ie. is dry), there is no flux
				Fluxes[fluxCellIndex + i] = 0.0f;
			}
		}
	}
}


__global__ void buildRValues(float *RValues, float *Fluxes, float *SlopeSource, float *WindSource, int m, int n)
{
	// Calculate the row and column of the thread within the thread block
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	
	// Calculate the cell's index in the RValues vector
	int RCellIndex = row*n*3 + col*3;
	
	// Actually, we don't need this because all three R values will evaluate to zero for a dry cell
	// Clear residual values from the R matrix in case a cell goes from wet to dry
	//for (int i=0; i<3; i++)
	//{
	//	RValues[RCellIndex+i] = 0.0f;
	//}
	
	// First check if the thread is operating on a cell inside of the block's two cell deep ghost cells
	if (col > 1 && row > 1 && col < n-2 && row < m-2)
	{
		// Calculate the indices of the flux arrays for the current cell (each cell stores it's own north and east values)
		int north = row*n*2*3 + col*2*3;		// Index of the flux array at the north interface
		int south = (row-1)*n*2*3 + col*2*3; 		// Index of the flux array at the south interface
		int east = row*n*2*3 + col*2*3 + 1*3;		// Index of the flux array at the east interface
		int west = row*n*2*3 + (col-1)*2*3 + 1*3;	// Index of the flux array at the west interface
		
		// Calculate the source index (slope and wind source matrices are the same size so we can use the same indices)
		int source = row*n*2 + col*2;
		
		// Calculate values of R
		RValues[RCellIndex] = 0.0f - (Fluxes[east] - Fluxes[west]) - (Fluxes[north] - Fluxes[south]);
		RValues[RCellIndex+1] = WindSource[source] + SlopeSource[source] - (Fluxes[east+1] - Fluxes[west+1]) - (Fluxes[north+1] - Fluxes[south+1]);
		RValues[RCellIndex+2] = WindSource[source+1] + SlopeSource[source+1] - (Fluxes[east+2] - Fluxes[west+2]) - (Fluxes[north+2] - Fluxes[south+2]);
	}
}



























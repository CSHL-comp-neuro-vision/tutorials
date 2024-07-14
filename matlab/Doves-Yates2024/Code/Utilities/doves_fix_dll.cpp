/**************************************************************************/
/**************************************************************************/

// This function is used to compute fixations given the x and y eye 
// tracker co-ordinates in VISUAL DEGREES. The algorithm is what is used by 
// the ASL routine. I had to write this since I did not want to save my 
// data in ASL's format and call their routine......

// PROGRAMMER - Umesh Rajashekar (umesh@ece.utexas.edu)
// DATE - May 19,2001 : Updated 13.05.04

/**************************************************************************/
/*** Include Libraries ****************************************************/
/**************************************************************************/

#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <malloc.h>

//function in MATLAB is to be called something like 
// [x_fix y_fix] = asl_fix(x_data,y_data,MIN_SAMPLES);
// THE DATA IS ASSUMED TO BE FED IN 1 ROW AND X COLUMNS (I give an error msg else)

/**************************************************************************/
/*** mexFunction **********************************************************/
/**************************************************************************/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//FUNCTION PROTOTYPES
int asl_fix(double x[],int x_size, double y[], int y_size,double x_fix[], double y_fix[], double fixt_start_index[], double fix_time[], double saccade_time[], int MIN_SAMPLES);
//umesh modified Feb 28 to add fix and saccade time params

int x_size				= 0;		//size of x eye data 
int y_size				= 0;		//size of y eye data (Yeah!! redundant.. but i used it like an idiot as a variable)

int num_x_cols			= mxGetM(prhs[0]);
int num_y_cols			= mxGetM(prhs[1]);

if ((num_x_cols !=1) || (num_y_cols!=1))
	mexErrMsgTxt("Sorry but this function accepts data only in 1 row and N cols");

x_size = mxGetN(prhs[0]);	//find the size of the input
y_size = mxGetN(prhs[1]);	//find the size of the input

if (x_size != y_size)
	mexErrMsgTxt("x and y data should be of same length");

//TEMPORARY OUTPUTS FROM FUNCTION (since only NFIX of these are valid)
double* x_fix			= new double [x_size];
double* y_fix			= new double [y_size];
double* fix_start_index = new double [x_size]; // Start index of the fixation - Umesh Feb 28,2002; Mar 28, 2003
double* fix_time		= new double [x_size]; // Fixation time in terms of num of samples - Umesh Feb 28,2002
double* saccade_time	= new double [x_size]; // Saccade time in terms of num of samples - Umesh Feb 28,2002

double* x = mxGetPr(prhs[0]);
double* y = mxGetPr(prhs[1]);
int min_samples=4;	//Minimum number of samples for fixation computation

if (nrhs == 3)
	min_samples = mxGetScalar(prhs[2]);


int num_fix = asl_fix(x,x_size,y,y_size,x_fix,y_fix,fix_start_index, fix_time,saccade_time,min_samples);

//int num_fix =1;
plhs[0] = mxCreateDoubleMatrix(1,num_fix,mxREAL);
plhs[1] = mxCreateDoubleMatrix(1,num_fix,mxREAL);
plhs[2] = mxCreateDoubleMatrix(1,num_fix,mxREAL);
plhs[3] = mxCreateDoubleMatrix(1,num_fix,mxREAL); // Return variable for the fixation time parameter - Umesh Feb 28,2002
plhs[4] = mxCreateDoubleMatrix(1,num_fix,mxREAL); // Return variable for the saccade time parameter - Umesh Feb 28,2002

//*(mxGetPr(plhs[0]))=1;

for ( int i =0;i<num_fix;i++)
{
	*((mxGetPr(plhs[0]))+i)= x_fix[i];
	*((mxGetPr(plhs[1]))+i)= y_fix[i];
	*((mxGetPr(plhs[2]))+i)= fix_start_index[i];
	*((mxGetPr(plhs[3]))+i)= fix_time[i];
	*((mxGetPr(plhs[4]))+i)= saccade_time[i];
}

//I allocated this memory and am responsible for deleting it
if (x_fix) 
	delete [] x_fix;
if (y_fix)
	delete [] y_fix;
}

/**************************************************************************/
/* This functions computes the fixations given the x and y eye data points /
/**************************************************************************/

int asl_fix(double x[],int x_size, double y[], int y_size,double x_fix[], double y_fix[], double fix_start_index [], double fix_time[], double saccade_time[], int MIN_SAMPLES)
{
	//CONSTANTS
	//int MIN_SAMPLES = 4; 
	double my_const = 0.0;	// be more generous in allowing for fixations
	double CR_1 = 0.5+my_const, CR_2 = 1+my_const, CR_3 = 1.5+my_const;
	int MAX_CR2POINTS = int(MIN_SAMPLES/2); //Max number of points that can exceed CR2 (3 in ASL)

	//FUNCTION PROTOTYPES
	double my_std(double* x,int begin_index, int NUM_SAMPLES);
	double my_mean(double* x,  int begin_index, int NUM_SAMPLES);

	int NFIX  =0 ;	//total number of fixations
	int n =0;		//initialize index

    while (1)
    {
    	//FIND THE TEMPORARY FIXATION CO-ORDINATES
        //	for (int n = 0;n < x_size-MIN_SAMPLES+1; n++)


    	int SAC_DUR = 0; // Duration of saccade			//Umesh added Feb 28, 2002

    	while( (n+MIN_SAMPLES-1) < x_size )
    	{		
    		double std_x = my_std(x,n,MIN_SAMPLES);
    		double std_y = my_std(y,n,MIN_SAMPLES);

    		if ( (std_x < CR_1) && (std_y < CR_1))
    			break;
    		saccade_time [NFIX]= SAC_DUR ++;
    		n=n+1;
    	}
			
    	if ( (n+MIN_SAMPLES-1) >= x_size )
    		break;			//break out of all loops if out of array

    	// "n" is the index of the first temporary fixation point
    	NFIX+=1	;// increase num of fixations found
    	fix_start_index [NFIX-1] = double(n); //Set the n as the start value of the fixation
    	int FIX_DUR = MIN_SAMPLES;	//Duration of the fixation
    	n = n + MIN_SAMPLES -1 ;	//Update index to next data point

    	double	X_t = my_mean(x,n-(MIN_SAMPLES-1),MIN_SAMPLES);	//temporary fixation point
    	double	Y_t = my_mean(y,n-(MIN_SAMPLES-1),MIN_SAMPLES);	//temporary fixation point
    
    	int count = 0;
    	double sum_x = X_t*MIN_SAMPLES; double sum_x_temp =X_t*MIN_SAMPLES;double mean_x_temp=0;
    	double sum_y = Y_t*MIN_SAMPLES; double sum_y_temp =Y_t*MIN_SAMPLES;double mean_y_temp=0;
    
    	//while (count <3)
    	while (count < MAX_CR2POINTS)
    	{
    		n = n+1;
    
    		if (n >= x_size)	//break out of while loop if exceeding array dimension
    			break;
    
    		FIX_DUR = FIX_DUR+1;	//original point of increasing fixation durations
    
    		double delta_x = fabs(X_t - x[n]);
    		double delta_y = fabs(Y_t - y[n]);
    
    
    		if ( (delta_x <= CR_2) && (delta_y<= CR_2) )
    		{
    
    			//FIX_DUR = FIX_DUR+1;	//Umesh put this
    
    			if (count > 0)
    			{
    				 sum_x = sum_x_temp;
    				 sum_y = sum_y_temp;
    				 if (sum_x_temp == X_t*MIN_SAMPLES)
    					 FIX_DUR-=count;
    			}
			
			sum_x = sum_x + x[n]; sum_x_temp = sum_x;
			sum_y = sum_y + y[n]; sum_y_temp = sum_y;

			count =0;
		}
		

		if ( (delta_x > CR_2) || (delta_y > CR_2) )
		{
			count = count + 1;
			
			if ( (delta_x <= CR_3) && (delta_y <= CR_3) )
			{
				sum_x_temp = sum_x_temp + x[n];
				sum_y_temp = sum_y_temp + y[n];
			}
				
		}

		// 3 consecutive points "outside" region but still mean is acceptable
		//if (count == 3)
		if (count == MAX_CR2POINTS)
		{
			 //mean_x_temp = my_mean(x,n-2,3); mean_y_temp = my_mean(y,n-2,3);
			
			mean_x_temp = my_mean(x,n-2,MAX_CR2POINTS); mean_y_temp = my_mean(y,n-2,MAX_CR2POINTS);
			if ( (fabs(mean_x_temp-X_t) <= CR_2) && (fabs(mean_y_temp-Y_t) <= CR_2) )
			{	
				//FIX_DUR = FIX_DUR+3;	//Umesh put this here
				sum_x = sum_x + x[n-2] + x[n-1]+x[n];sum_x_temp = sum_x;
				sum_y = sum_y + y[n-2] + y[n-1]+y[n];sum_y_temp = sum_y;
				count = 0;
			}

		}
	
	}	// loop for count <3
	
	if (n >= x_size)		// break out of all loops if n exceeds index dimensions
	{
		FIX_DUR-=count;
		x_fix[NFIX-1] = sum_x/FIX_DUR;
		y_fix[NFIX-1] = sum_y/FIX_DUR;
	    fix_time[NFIX-1]=FIX_DUR;
		break;

	}
	
	// 3 consecutive points outside region and mean diff failure
	if ( (fabs(mean_x_temp-X_t) > CR_2) || (fabs(mean_y_temp-Y_t) > CR_2) )
	{
		//FIX_DUR-=3;
		FIX_DUR-=MAX_CR2POINTS;
		x_fix[NFIX-1] = sum_x/FIX_DUR;;
		y_fix[NFIX-1] = sum_y/FIX_DUR;
		fix_time[NFIX-1]=FIX_DUR;		
	}
				
}	//end of while(1)

//for (int i = 0;i<NFIX;i++)
//	cout<<x_fix[i] << "\t" << y_fix[i] << endl;

return (NFIX);
}

//*************************************************************************************************

double my_std(double* x, int begin_index, int NUM_SAMPLES)
{
	double temp_x1  = 0;
	double temp_x2  = 0;

	for (int i= begin_index; i < begin_index+NUM_SAMPLES; i++)
	{
		temp_x1 = temp_x1 + x[i]*x[i];
		temp_x2 = temp_x2 + x[i];
	}
	
	temp_x1         = temp_x1/NUM_SAMPLES;
	temp_x2         = temp_x2/NUM_SAMPLES;
	double std_x    = temp_x1 - temp_x2*temp_x2;
	std_x           = sqrt(std_x);
	
	return(std_x);
}

//*************************************************************************************************

double my_mean(double* x,  int begin_index, int NUM_SAMPLES)
{
	double mean_x   = 0;
	for (int i=begin_index;i < begin_index + NUM_SAMPLES ;i++)
	{
		mean_x+= x[i];
	}
	mean_x/=NUM_SAMPLES;	
	return(mean_x);
}












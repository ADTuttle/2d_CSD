#include <stdlib.h>
#include <math.h>

double array_max(double *array,size_t size)
{
	double max=0;
	for(size_t ind=0;ind<size;ind++)
	{
		if(fabs(array[ind])>max)
		{
			max=fabs(array[ind]);
		}
	}
	return max;
}
double array_diff_max(double *array1,double *array2,size_t size)
{
	double max=0;
	for(size_t ind=0;ind<size;ind++)
	{
		if(fabs(array1[ind]-array2[ind])>max)
		{
			max=fabs(array1[ind]-array2[ind]);
		}
	}
	return max;
}

double l2_norm(double *array1,size_t size)
{
	double max=0;
	for(size_t ind=0;ind<size;ind++)
	{
		max += array1[ind]*array1[ind];
	}
	return sqrt(max);
}

#include <stdio.h>
#include <math.h>
#define ERROR 1e-5

static double cor_coeff_rect[7][2]={
       0,1.1,
       0.1,1.08,
       0.2,1.06,
       0.4,1.04,
       0.6,1.02,
       0.8,1.01,
       1.0,1.0,
};

static double cor_coeff_tri[8][2]={
       0,0.75,
       10,0.84,
       20,0.89,
       30,0.93,
       40,0.96,
       60,0.98,
       80,0.99,
       90,1.0,
};

static double cor_coeff_ecllipse[10][2]={
       0.1,1.21,
       0.2,1.16,
       0.3,1.11,
       0.4,1.08,
       0.5,1.05,
       0.6,1.03,
       0.7,1.02,
       0.8,1.01,
       0.9,1.01,
       1.0,1.0,
};

int inter_rect(double width_length_ratio, double* cor_coeff){
       int index;
       double k;
       if(width_length_ratio<1e-5||width_length_ratio>1)
       {
       		return -1;
       }
       else if(width_length_ratio>0&&width_length_ratio<0.1)
       {
       		index=0;
       }	
       else if(width_length_ratio<0.2)
       {
       		index=1;
       }
       else if(width_length_ratio<0.4)
       {
       		index=2;
       }		
       else if(width_length_ratio<0.6)
       {
       		index=3;
       }
	   else if(width_length_ratio<0.8)
	   {
   			index=4;
	   }
	   else if(width_length_ratio<1||fabs(width_length_ratio-1)<ERROR)
	   {
   			index=5;
	   }
	   k=(cor_coeff_rect[index+1][1]-cor_coeff_rect[index][1])/(cor_coeff_rect[index+1][0]-cor_coeff_rect[index][0]);
	   *cor_coeff=k*(width_length_ratio-cor_coeff_rect[index][0])+cor_coeff_rect[index][1];
	   return 0;
}

int inter_tri(double half_ang, double* cor_coeff)
{
	int index;
	double k;
	if(half_ang<ERROR||half_ang>90||fabs(half_ang-90)<ERROR)
	{
		return -1;
	}
	else if(half_ang>ERROR&&half_ang<10)
	{
		index=0;
	}
	else if(half_ang<20)
	{
		index=1;
	}
	else if(half_ang<30)
	{
		index=2;
	}
	else if(half_ang<40)
	{
		index=3;
	}
	else if(half_ang<60)
	{
		index=4;
	}
	else if(half_ang<80)
	{
		index=5;
	}
	else if(half_ang<90)
	{
		index=6;
	}
	k=(cor_coeff_tri[index+1][1]-cor_coeff_tri[index][1])/(cor_coeff_tri[index+1][0]-cor_coeff_tri[index][0]);
 	*cor_coeff=k*(half_ang-cor_coeff_tri[index][0])+cor_coeff_tri[index][1];
	return 0;	
}

int inter_ec(double radius_ratio, double* cor_coeff)
{
	int index;
	double k;
	if(radius_ratio<ERROR||radius_ratio>1+ERROR)
	{
		return -1;
	}
	else if(radius_ratio>ERROR&&radius_ratio<0.1)
	{
		index=0;
	}
	else if(radius_ratio<0.2)
	{
		index=1;
	}	
	else if(radius_ratio<0.3)
	{
		index=2;
	}
	else if(radius_ratio<0.4)
	{
		index=3;
	}
	else if(radius_ratio<0.5)
	{
		index=4;
	}
	else if(radius_ratio<0.6)
	{
		index=5;
	}
	else if(radius_ratio<0.7)
	{
		index=6;
	}
	else if(radius_ratio<0.8)
	{
		index=7;
	}
	else if(radius_ratio<0.9)
	{
		index=8;
	}
	else if(radius_ratio<1||fabs(radius_ratio-1)<ERROR)
	{
		index=9;
	}
	k=(cor_coeff_ecllipse[index+1][1]-cor_coeff_ecllipse[index][1])/(cor_coeff_ecllipse[index+1][0]-cor_coeff_ecllipse[index][0]);
 	*cor_coeff=k*(radius_ratio-cor_coeff_ecllipse[index][0])+cor_coeff_ecllipse[index][1];
	return 0;	
}      

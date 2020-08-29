#include<math.h>
#include"data_structure.h"

double karman_nikuradse(double re_num,double x)
{
	return 1/x+0.4-1.737*log(re_num*x);
}

double karman_nikuradse_diff(double x)
{
	return -1/(x*x)-1.737/x;
}

int karman_nikuradse_iter(double* x_final,double re_num,double x_ini,int max_fric_iter_times)
{
	double x_this,x_prev;
	int i=0;
	x_prev=x_ini-karman_nikuradse(re_num,x_ini)/karman_nikuradse_diff(x_ini);
	do
	{
		x_this=x_prev-karman_nikuradse(re_num,x_prev)/karman_nikuradse_diff(x_prev);
//		printf("%d,%lf,%lf,%lf,%.10lf\n",i,x_prev,x_this,fabs(x_this-x_prev),karman_nikuradse_diff(x_prev));
		if(fabs(x_this-x_prev)<ERR)
		{
			break;
		}
		else
		{
			x_prev=x_this;
		}
		i++;
	}
	while(i<max_fric_iter_times);
	if(i==max_fric_iter_times)
	{
		return -1;
	}
	else
	{
		*x_final=x_this;
	}
	return 0;
}

double round_rib_fric(double hydra_dia, double rib_height, double rib_space)
{
	double cf_temp;
	cf_temp=2.5*log(hydra_dia/(2*rib_height))+0.95*pow((rib_space/rib_height),0.53)-3.75;
	return 2/pow(cf_temp,2);
}

double rec_rib_fric(double re_num, double hydra_dia, double rib_height, double rib_space, double h_w_ratio)
{
	double r_plus,cf0;
	r_plus=1.429*pow((rib_space/rib_height),0.35);
	cf0=0.046*pow(re_num,-0.2);
	return (h_w_ratio*cf0+2/pow((r_plus-2.5-2.5*log((2*rib_height/hydra_dia)*2/(1+h_w_ratio))),2))/(1+h_w_ratio);
}

//for C3X calculation only
double C3X_nu(double re_num, double prandtl_num, double cr_num)
{
	return cr_num*0.022*pow(prandtl_num,0.5)*pow(re_num,0.8);
}
//for C3X calculation only


double dittus_boelter(double re_num, double prandtl_num)
{
	return 0.023*pow(prandtl_num,0.5)*pow(re_num,0.8);
}

double petukhov_kirillov(double fric_coeff, double re_num, double prandtl_num)
{
	double k1,k2;
	k1=fric_coeff*re_num*prandtl_num/2;
	k2=1.07+12.7*sqrt(fric_coeff/2)*(pow(prandtl_num,2.0/3.0)-1);
	return k1/k2;	
}

double gnielinski(double fric_coeff, double re_num, double prandtl_num)
{
	double k1,k2;
	k1=fric_coeff*(re_num-1000)*prandtl_num/2;
	k2=1+12.7*sqrt(fric_coeff/2)*(pow(prandtl_num,2.0/3.0)-1);
	return k1/k2;
}

double rec_smooth_nu(double re_num, double prandtl_num, double h_w_ratio)
{
	return 0.023*pow(re_num,0.8)*pow(prandtl_num,0.4)*pow((1+0.156*h_w_ratio),-0.1);
}


double round_rib_nu(double fric_coeff, double re_num, double hydr_dia, double rib_height, double prandtl_num, double rib_space)
{
	double h_plus;
	h_plus=rib_height*re_num*sqrt(fric_coeff/2)/hydr_dia;
	return re_num*prandtl_num*(fric_coeff/2)/(1+sqrt(fric_coeff/2)*(4.5*pow(h_plus,0.28)*pow(prandtl_num,0.57)-0.95*pow((rib_space/rib_height),0.53)));
}

double rec_rib_nu(double fric_coeff, double re_num, double prandtl_num, double hydr_dia, double rib_height, double rib_space, double h_w_ratio)
{
	double r_plus,g_plus,h_plus,cf0;
	h_plus=rib_height*re_num*sqrt(fric_coeff/2)/hydr_dia;
	g_plus=3.7*pow(h_plus,0.28);
	r_plus=1.429*pow((rib_space/rib_height),0.35);
	return re_num*prandtl_num*0.5*(fric_coeff+h_w_ratio*(fric_coeff-cf0))/((g_plus-r_plus)*sqrt((fric_coeff+h_w_ratio*(fric_coeff-cf0))/2)+1);
}

double para_wall_nu(double re_num, double prandtl_num, double wall_dist, double stream_length, double visc_cool, double visc_wall)
{
	return pow(prandtl_num,1.0/3.0)*pow((visc_cool/visc_wall),0.14)*0.02*pow(re_num,0.8)*(1+4.6*wall_dist/stream_length);
}

double impinge_stag_nu(double re_num, double hole_dia, double hole_pitch, double target_dia, double target_dist)
{
	return 0.44*pow(re_num,0.7)*pow(hole_dia/hole_pitch,0.8)*exp(-0.85*(target_dist/hole_dia)*(hole_dia/hole_pitch)*pow(hole_dia/target_dia,0.4));
}

double impinge_head_nu(double re_num, double hole_dia, double hole_pitch, double target_dia, double target_dist)
{
	return 0.63*pow(re_num,0.7)*pow(hole_dia/hole_pitch,0.5)*pow(hole_dia/target_dia,0.6)*exp(-1.27*(target_dist/hole_dia)*pow(hole_dia/hole_pitch,0.5)*pow(hole_dia/target_dia,1.2));
}

double impinge_body_h(double ave_velo, double hole_pitch, double dynamic_visc, double ave_dens)
{
	double Kc=4;
	double re_cool;
	re_cool=ave_velo*hole_pitch*ave_dens/dynamic_visc;
	return 0.286*pow(re_cool,0.625)*Kc/hole_pitch;
}

double pinfin_fric(double re_num, double pitch, double diameter, double pinfin_row)
{
	double fp;
	fp=2.06*pow(pitch/diameter,1.1)*pow(re_num,-0.16);
	return 4*fp*(pinfin_row-1);
}

double pinfin_nu(double re_num, double prandtl_num)
{
	return 0.248*pow(re_num,0.594)*pow(prandtl_num,0.333);
}

double colebrook(double roughness_factor,double re_num,double x)
{
    double A,B;
    A=roughness_factor/3.7;
    B=2.51/re_num;
    return x+2*log10(A+B*x);
}
double colebrook_diff(double roughness_factor,double re_num,double x)
{
    double A,B;
    A=roughness_factor/3.7;
    B=2.51/re_num;
    return 1+2*B*log10(2.718281829)/(A+B*x);
}


int colebrook_iter(double* x_final,double roughness_factor,double re_num,double x_ini,int max_fric_iter_times)
{
    double x_this,x_prev;
    int i=0;
    
    x_prev=x_ini-colebrook(roughness_factor,re_num,x_ini)/colebrook_diff(roughness_factor,re_num,x_ini);
    do
    {
        x_this=x_prev-colebrook(roughness_factor,re_num,x_prev)/colebrook_diff(roughness_factor,re_num,x_prev);
        if(fabs(x_this-x_prev)<ERR)
        { 
            break;
        }
        else
        {
            x_prev=x_this;
			i++;
        }
        
        
    }
    while(i<max_fric_iter_times);
    if(i==max_fric_iter_times)
    {
        return -1;
    }
    else
    {
        *x_final=x_this;
    }
    return 0;
}

//int main()
//{
//	double x;
//	colebrook_iter(&x,0.005,95976.76370,0.1,100);
//	printf(",,,,,,%lf\n",x);
//}


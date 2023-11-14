/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* This is the solver for FLOWNET ANALYSIS.                                                              */
/* VERSION 2.0                                                                                           */
/*                                                                                                       */
/* WANG ZHENRONG (Edison. WANG), an independent C program developer, reserves all rights of this program.*/
/* Contacts: +8613661536219(cel.) zhenrong_w@163.com(email).                                             */
/* No license has been applied for this code.                                                            */
/*                                                                                                       */
/* Main function:  solving flow network problem.                                                         */
/* Fluid Type:     WATER, STEAM, AIR (ideal gas).                                                        */
/* Problem type:   [STEADY STATE!] flow under reference temperature/ flow with heat effects.             */
/* Geom Structure: round/rect/valve/pipe/pump...                                                         */
/* Important NOTE: For more detailed information, please read the help doc.                              */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include"data_structure.h"
#ifdef _WIN32
#include"..\\_fluid_air\\property_calculation.c"
#include"..\\_fluid_steam\\steam_property_calc.c"
#include"..\\_fluid_water\\water_property.c"
#else
#include"../_fluid_air/property_calculation.c"
#include"../_fluid_steam/steam_property_calc.c"
#include"../_fluid_water/water_property.c"
#endif
#include"flow_thermal_coeffs.c"
#include"cor_coeff_cal.c"


void welcome(FILE* output_console)
{
    time_t rtime;
    struct tm* timeinfo=NULL;
    time(&rtime);
    timeinfo=localtime(&rtime);

    printf("\n# FLOW NETWORK CALCULATION PROGRAM\n");
    printf("# VERSION 0.1\n");
    printf("# WangZhr. ALL RIGHTS RESERVED.\n\n");
    printf("# W E L C O M E !\n");
    printf("# CURRENT DATE AND TIME: %s\n",asctime(timeinfo));
    if(output_console!=NULL)
    {
        fprintf(output_console,"\n\n**********************************************************************\t%s",asctime(timeinfo));
        fprintf(output_console,"\n# FLOW NETWORK CALCULATION PROGRAM\n# VERSION 0.1\n");
        fprintf(output_console,"# WangZhr. ALL RIGHTS RESERVED.\n\n");
        fprintf(output_console,"# W E L C O M E !\n");
        fprintf(output_console,"# CURRENT DATE AND TIME: %s\n",asctime(timeinfo));
    }
}

int read_setup(FILE* setup_file,setup* p_setup,FILE* output_console)
{
    char choice;

    fscanf(setup_file,"%d%d%d",&p_setup->snode_num_max,&p_setup->selem_num_max,&p_setup->sbound_num_max);
    fscanf(setup_file,"%d",&p_setup->sfluid_type);
    fscanf(setup_file,"%d%lf",&p_setup->senergy_flag,&p_setup->sref_temp);
    fscanf(setup_file,"%d%d",&p_setup->shtc_flag,&p_setup->scorrelation_flag);
    fscanf(setup_file,"%d%d",&p_setup->sfriction_flag,&p_setup->sfric_coeff_flag);
    fscanf(setup_file,"%lf%lf%lf",&p_setup->sini_pres,&p_setup->sini_flow_rate,&p_setup->sini_temp);
    fscanf(setup_file,"%d%d%d%lf%lf",&p_setup->max_fric_iter_times,&p_setup->max_flow_iter_times,&p_setup->max_total_iter_times,&p_setup->flow_converge_flag,&p_setup->energy_converge_flag);

    if(p_setup->snode_num_max<1||p_setup->selem_num_max<1||p_setup->sbound_num_max<1)
    {
        printf("FATAL ERROR: NODE_NUM/ELEM_NUM/BOUNDARY_NODE_NUM SHOULD BE A POSITIVE NUMBER.\nPLEASE CHECK SETUP FILE(LINE 1).\nPROGRAM ABORTED.\n");
        if(output_console!=NULL)
        {
            fprintf(output_console,"FATAL ERROR: NODE_NUM/ELEM_NUM/BOUNDARY_NODE_NUM SHOULD BE A POSITIVE NUMBER.\nPLEASE CHECK SETUP FILE(LINE 1).\nPROGRAM ABORTED.\n");
        }
        return -1;
    }

    if(p_setup->sfluid_type<0||p_setup->sfluid_type>2)
    {
        printf("FATAL ERROR: FLUID TYPE SETUP ERROR.\nPLEASE CHECK SETUP FILE(LINE 2).\nPROGRAM ABORTED.\n");
        if(output_console!=NULL)
        {
            fprintf(output_console,"FATAL ERROR: FLUID TYPE SETUP ERROR.\nPLEASE CHECK SETUP FILE(LINE 2).\nPROGRAM ABORTED.\n");
        }
        return -1;
    }

    if(p_setup->senergy_flag!=0&&p_setup->senergy_flag!=1)
    {
        printf("FATAL ERROR: ENERGY EQUATION SETUP ERROR.\nPLEASE CHECK SETUP FILE(LINE 3).\nPROGRAM ABORTED.\n");
        if(output_console!=NULL)
        {
            fprintf(output_console,"FATAL ERROR: ENERGY EQUATION SETUP ERROR.\nPLEASE CHECK SETUP FILE(LINE 3).\nPROGRAM ABORTED.\n");
        }
        return -1;
    }
	if(p_setup->senergy_flag==1)
	{
		if(p_setup->shtc_flag!=0&&p_setup->shtc_flag!=1)
    	{
        	printf("FATAL ERROR: HEAT TRANSFER COEFFICIENT SETUP ERROR.\nPLEASE CHECK SETUP FILE(LINE 4).\nPROGRAM ABORTED.\n");
        	if(output_console!=NULL)
        	{
            	fprintf(output_console,"FATAL ERROR: HEAT TRANSFER COEFFICIENT SETUP ERROR.\nPLEASE CHECK SETUP FILE(LINE 4).\nPROGRAM ABORTED.\n");
        	}
        	return -1;
    	}
    	
		if(p_setup->scorrelation_flag!=1&&p_setup->scorrelation_flag!=2&&p_setup->scorrelation_flag!=0)
    	{
        	printf("WARNING: HEAT TRANSFER CORRELATION CHOICE ERROR.\n PLEASE CHECK SETUP FILE(LINE 4).\nAUTO CHOOSE DITTUS-BOELTER CORRELATION.\n");
        	if(output_console!=NULL)
        	{
            	fprintf(output_console,"WARNING: HEAT TRANSFER CORRELATION CHOICE ERROR.\n PLEASE CHECK SETUP FILE(LINE 4).\nAUTO CHOOSE DITTUS-BOELTER CORRELATION.\n");
        	}
        	p_setup->scorrelation_flag==1;
    	}
    	
		if(p_setup->shtc_flag==1&&p_setup->scorrelation_flag==0&&p_setup->selem_num_max!=10)
    	{
        	printf("WARNING: HEAT TRANSFER CORRELATION CHOICE ERROR.\n PLEASE CHECK SETUP FILE(LINE 4).\nAUTO CHOOSE DITTUS-BOELTER CORRELATION.\n");
        	if(output_console!=NULL)
        	{
            	fprintf(output_console,"WARNING: HEAT TRANSFER CORRELATION CHOICE ERROR.\n PLEASE CHECK SETUP FILE(LINE 4).\nAUTO CHOOSE DITTUS-BOELTER CORRELATION.\n");
        	}
        	p_setup->scorrelation_flag==1;
    	}
	}
    
    if(p_setup->sfriction_flag!=0&&p_setup->sfriction_flag!=1)
    {
        printf("FATAL ERROR: FRICTION SETUP ERROR.\nPLEASE CHECK SETUP FILE(LINE 5).\n{PROGRAM ABORTED.\n");
        if(output_console!=NULL)
        {
            fprintf(output_console,"FATAL ERROR: FRICTION SETUP ERROR.\nPLEASE CHECK SETUP FILE(LINE 5).\n{PROGRAM ABORTED.\n");
        }
        return -1;
    }
    if(p_setup->sfric_coeff_flag!=0&&p_setup->sfric_coeff_flag!=1)
    {
    	printf("FATAL ERROR: FRICTION SETUP ERROR.\nPLEASE CHECK SETUP FILE(LINE 5).\n{PROGRAM ABORTED.\n");
        if(output_console!=NULL)
        {
            fprintf(output_console,"FATAL ERROR: FRICTION SETUP ERROR.\nPLEASE CHECK SETUP FILE(LINE 5).\n{PROGRAM ABORTED.\n");
        }
        return -1;
    }

    if(p_setup->senergy_flag==0&&p_setup->sref_temp<ERR)
    {
        printf("FATAL ERROR: REFERENCE TEMPERATURE SHOULD BE POSITIVE.\nPLEASE CHECK SETUP FILE(LINE 3).\nPROGRAM ABORTED.\n");
        if(output_console!=NULL)
        {
            fprintf(output_console,"FATAL ERROR: REFERENCE TEMPERATURE SHOULD BE POSITIVE.\nPLEASE CHECK SETUP FILE(LINE 3).\nPROGRAM ABORTED.\n");
        }
        return -1;
    }

    if(p_setup->max_fric_iter_times<1||p_setup->max_flow_iter_times<1||p_setup->max_total_iter_times<1)
    {
        printf("FATAL ERROR: ITERATION TIMES SHOULD BE POSITIVE.\nPLEASE CHECK SETUP FILE(LINE 7).\nPROGRAM ABORTED.\n");
        if(output_console!=NULL)
        {
            fprintf(output_console,"FATAL ERROR: ITERATION TIMES SHOULD BE POSITIVE.\nPLEASE CHECK SETUP FILE(LINE 7).\nPROGRAM ABORTED.\n");
        }
        return -1;
    }
    if(p_setup->flow_converge_flag<1e-9||p_setup->energy_converge_flag<1e-9)
    {
        printf("WARNING: CONVERGENCE CRITERIA SETUP ERROR. AUTO SET TO DEFAULT VALUE 1E-6.\n");
        if(output_console!=NULL)
        {
            fprintf(output_console,"WARNING: CONVERGENCE CRITERIA SETUP ERROR. AUTO SET TO DEFAULT VALUE 1E-6.\n");
        }
        p_setup->flow_converge_flag=1e-6;
        p_setup->energy_converge_flag=1e-6;
    }

    if(p_setup->senergy_flag==1&&fabs(p_setup->sini_temp)<1e-4)
    {
        printf("WARNING: INITIALIZATION NUMBER OF TEMPERATURE FIELD SHOULD NOT BE ZERO. AUTO SET TO DEFAULT VALUE 500.\n");
        if(output_console!=NULL)
        {
            fprintf(output_console,"WARNING: INITIALIZATION NUMBER OF TEMPERATURE FIELD SHOULD NOT BE ZERO. AUTO SET TO DEFAULT VALUE 500.\n");
        }
        p_setup->sini_temp=500;
    }

    if(fabs(p_setup->sini_pres)<1e-4||fabs(p_setup->sini_flow_rate)<1e-4)
    {
        printf("WARNING: INITIALIZATION NUMBER OF FLOW FIELD SHOULD NOT BE ZERO. AUTO SET TO DEFAULT VALUE (P:1E6, MFR:10).\n");
        if(output_console!=NULL)
        {
            fprintf(output_console,"WARNING: INITIALIZATION NUMBER OF FLOW FIELD SHOULD NOT BE ZERO. AUTO SET TO DEFAULT VALUE (P:1E6, MFR:10).");
        }
        p_setup->sini_pres=1e6;
        p_setup->sini_flow_rate=10;
    }
    return 0;
}

int create_elem_array(elem* elem_array,FILE* elem_file,setup* p_setup)
{
    int i=0;
    double width_length_ratio;
    double cor_coeff_temp;
    double a1,a2,h,phi1,phi2;
    double a0,beta;
    double b0;
    double k1,k2,k3,k4,k5;
    double pinfins,pinfinx,pinfind,pinfinx1;
		
    do
    {
  		(elem_array+i)->elem_type=0;
		(elem_array+i)->elem_num=-1;
		(elem_array+i)->inlet_node_num=-1;
		(elem_array+i)->outlet_node_num=-1;
		(elem_array+i)->length=0;
		(elem_array+i)->diameter=0;
		(elem_array+i)->section_length=0;
		(elem_array+i)->section_width=0;
		(elem_array+i)->wall_distance=0;
		(elem_array+i)->pitch=0;
		(elem_array+i)->target_dia=0;
		(elem_array+i)->target_dist=0;
		(elem_array+i)->cor_coeff=0;
		(elem_array+i)->pinfin_row=0;
		(elem_array+i)->incline_angle=0;
		(elem_array+i)->param_1=0;
		(elem_array+i)->param_2=0;
		(elem_array+i)->roughness_factor=0;
		(elem_array+i)->rib_height=0;
		(elem_array+i)->rib_space=0;
		(elem_array+i)->resist_factor=0;
		(elem_array+i)->flow_angle=0;
		(elem_array+i)->flow_rate=0;
		(elem_array+i)->ave_dens=0;
		(elem_array+i)->ave_temp=0;
		(elem_array+i)->re_num=0;
		(elem_array+i)->inlet_area=0;
		(elem_array+i)->inlet_radius=0;
		(elem_array+i)->outlet_area=0;
		(elem_array+i)->outlet_radius=0;
		(elem_array+i)->krot=1;
		(elem_array+i)->ave_area=0;
		(elem_array+i)->hydr_dia=0;
		(elem_array+i)->omega=0;
		(elem_array+i)->ht_coeff=0;
		(elem_array+i)->ht_area=0;
		(elem_array+i)->elem_specific_heat=0;
		(elem_array+i)->ht_temp=0;
		(elem_array+i)->fluid_cp=0;
		(elem_array+i)->fluid_visc=0;
		(elem_array+i)->fluid_thcond=0;
		(elem_array+i)->fluid_pranum=0;	

        fscanf(elem_file,"%d",&(elem_array+i)->elem_type);
        fscanf(elem_file,"%d%d%d",&(elem_array+i)->elem_num,&(elem_array+i)->inlet_node_num,&(elem_array+i)->outlet_node_num);

        if((elem_array+i)->elem_type<1)
        {
            return -1;
        }
//round_no_rib
        if((elem_array+i)->elem_type==1)
        {
            fscanf(elem_file,"%lf%lf%lf%lf",&(elem_array+i)->length,&(elem_array+i)->flow_angle,&(elem_array+i)->resist_factor,&(elem_array+i)->roughness_factor);
            fscanf(elem_file,"%lf%lf%lf%lf%lf",&(elem_array+i)->diameter,&(elem_array+i)->inlet_area,&(elem_array+i)->inlet_radius,&(elem_array+i)->outlet_area,&(elem_array+i)->outlet_radius);
            fscanf(elem_file,"%lf%lf%lf%lf%lf",&(elem_array+i)->omega,&(elem_array+i)->krot,&(elem_array+i)->ht_coeff,&(elem_array+i)->ht_area,&(elem_array+i)->ht_temp);
            (elem_array+i)->ave_area=0.5*((elem_array+i)->inlet_area+(elem_array+i)->outlet_area);
            (elem_array+i)->hydr_dia=(elem_array+i)->diameter;
            (elem_array+i)->cor_coeff=1;
        }
//rect_no_rib
        else if((elem_array+i)->elem_type==2)
        {
            fscanf(elem_file,"%lf%lf%lf%lf",&(elem_array+i)->length,&(elem_array+i)->flow_angle,&(elem_array+i)->resist_factor,&(elem_array+i)->roughness_factor);
            fscanf(elem_file,"%lf%lf%lf%lf%lf%lf",&(elem_array+i)->section_length,&(elem_array+i)->section_width,&(elem_array+i)->inlet_area,&(elem_array+i)->inlet_radius,&(elem_array+i)->outlet_area,&(elem_array+i)->outlet_radius);
            fscanf(elem_file,"%lf%lf%lf%lf%lf",&(elem_array+i)->omega,&(elem_array+i)->krot,&(elem_array+i)->ht_coeff,&(elem_array+i)->ht_area,&(elem_array+i)->ht_temp);
            (elem_array+i)->ave_area=0.5*((elem_array+i)->inlet_area+(elem_array+i)->outlet_area);
            width_length_ratio=(elem_array+i)->section_width/(elem_array+i)->section_length;
            if(inter_rect(width_length_ratio,&cor_coeff_temp)==-1)
            {
                return -1;
            }
            (elem_array+i)->cor_coeff=cor_coeff_temp;
            (elem_array+i)->hydr_dia=2*(elem_array+i)->section_length*(elem_array+i)->section_width/((elem_array+i)->section_length+(elem_array+i)->section_width);
        }
//ladder_no_rib
        else if((elem_array+i)->elem_type==3)
        {
            fscanf(elem_file,"%lf%lf%lf%lf",&(elem_array+i)->length,&(elem_array+i)->flow_angle,&(elem_array+i)->resist_factor,&(elem_array+i)->roughness_factor);
            fscanf(elem_file,"%lf%lf%lf%lf%lf",&a1,&a2,&h,&phi1,&phi2);
            k1=h/(a1+a2);
            k2=1/sin(phi1*PI/180)+1/sin(phi2*PI/180);
            k3=2*h/(1+k1*k2);
            (elem_array+i)->hydr_dia=k3;
            width_length_ratio=2*h/(a1+a2);
            if(inter_rect(width_length_ratio,&cor_coeff_temp)==-1)
            {
                return -1;
            }
            (elem_array+i)->cor_coeff=cor_coeff_temp;
            fscanf(elem_file,"%lf%lf%lf%lf",&(elem_array+i)->inlet_area,&(elem_array+i)->inlet_radius,&(elem_array+i)->outlet_area,&(elem_array+i)->outlet_radius);
            fscanf(elem_file,"%lf%lf%lf%lf%lf",&(elem_array+i)->omega,&(elem_array+i)->krot,&(elem_array+i)->ht_coeff,&(elem_array+i)->ht_area,&(elem_array+i)->ht_temp);
            (elem_array+i)->ave_area=0.5*((elem_array+i)->inlet_area+(elem_array+i)->outlet_area);
        }
//triangle_no_rib
        else if ((elem_array+i)->elem_type==4)
        {
            fscanf(elem_file,"%lf%lf%lf%lf",&(elem_array+i)->length,&(elem_array+i)->flow_angle,&(elem_array+i)->resist_factor,&(elem_array+i)->roughness_factor);
            fscanf(elem_file,"%lf%lf%lf",&a0,&h,&beta);
            k1=tan(beta*PI/180);
            k2=1/(k1*k1);
            k3=sqrt(1+k2);
            k4=2*h/(1+k3);
            (elem_array+i)->hydr_dia=k4;
            if(inter_tri(beta,&cor_coeff_temp)==-1)
            {
                return -1;
            }
            (elem_array+i)->cor_coeff=cor_coeff_temp;
            fscanf(elem_file,"%lf%lf%lf%lf",&(elem_array+i)->inlet_area,&(elem_array+i)->inlet_radius,&(elem_array+i)->outlet_area,&(elem_array+i)->outlet_radius);
            fscanf(elem_file,"%lf%lf%lf%lf%lf",&(elem_array+i)->omega,&(elem_array+i)->krot,&(elem_array+i)->ht_coeff,&(elem_array+i)->ht_area,&(elem_array+i)->ht_temp);
            (elem_array+i)->ave_area=0.5*((elem_array+i)->inlet_area+(elem_array+i)->outlet_area);
        }
//eclipse_no_rib
        else if((elem_array+i)->elem_type==5)
        {
            fscanf(elem_file,"%lf%lf%lf%lf",&(elem_array+i)->length,&(elem_array+i)->flow_angle,&(elem_array+i)->resist_factor,&(elem_array+i)->roughness_factor);
            fscanf(elem_file,"%lf%lf",&a0,&b0);
            k1=sqrt(a0*b0);
            k2=1.5*(a0+b0);
            k3=4*a0*b0/(k2-k1);
            (elem_array+i)->hydr_dia=k3;
            if(inter_ec(b0/a0,&cor_coeff_temp)==-1)
            {
                return -1;
            }
            (elem_array+i)->cor_coeff=cor_coeff_temp;
            fscanf(elem_file,"%lf%lf%lf%lf",&(elem_array+i)->inlet_area,&(elem_array+i)->inlet_radius,&(elem_array+i)->outlet_area,&(elem_array+i)->outlet_radius);
            fscanf(elem_file,"%lf%lf%lf%lf%lf",&(elem_array+i)->omega,&(elem_array+i)->krot,&(elem_array+i)->ht_coeff,&(elem_array+i)->ht_area,&(elem_array+i)->ht_temp);
            (elem_array+i)->ave_area=0.5*((elem_array+i)->inlet_area+(elem_array+i)->outlet_area);
        }
//film_cooling_hole
        else if((elem_array+i)->elem_type==6)
        {
            fscanf(elem_file,"%lf%lf%lf%lf",&(elem_array+i)->length,&(elem_array+i)->incline_angle,&(elem_array+i)->resist_factor,&(elem_array+i)->roughness_factor);
            fscanf(elem_file,"%lf%lf%lf",&(elem_array+i)->diameter,&(elem_array+i)->inlet_area,&(elem_array+i)->outlet_area);
            fscanf(elem_file,"%lf%lf%lf",&(elem_array+i)->ht_coeff,&(elem_array+i)->ht_area,&(elem_array+i)->ht_temp);
            (elem_array+i)->ave_area=0.5*((elem_array+i)->inlet_area+(elem_array+i)->outlet_area);
            (elem_array+i)->hydr_dia=(elem_array+i)->diameter;
            (elem_array+i)->cor_coeff=1;
        }
        else if((elem_array+i)->elem_type==7)
        {
        	fscanf(elem_file,"%lf%lf%lf",&(elem_array+i)->param_1,&(elem_array+i)->param_2,&(elem_array+i)->ave_area);
        	(elem_array+i)->hydr_dia=sqrt(4*(elem_array+i)->ave_area/PI);
        }
        
        else if((elem_array+i)->elem_type==8)
        {
        	fscanf(elem_file,"%lf%lf%lf",&(elem_array+i)->param_1,&(elem_array+i)->param_2,&(elem_array+i)->hydr_dia);
        	(elem_array+i)->ave_area=PI*(elem_array+i)->hydr_dia*(elem_array+i)->hydr_dia/4;
		}
		
		else if((elem_array+i)->elem_type==9)
		{
			fscanf(elem_file,"%lf%lf%lf%lf%lf%lf%lf",&(elem_array+i)->param_1,&(elem_array+i)->inlet_area,&(elem_array+i)->inlet_radius,&(elem_array+i)->outlet_area,&(elem_array+i)->outlet_radius,&(elem_array+i)->omega,&(elem_array+i)->krot);
			(elem_array+i)->ave_area=0.5*((elem_array+i)->inlet_area+(elem_array+i)->outlet_area);
			(elem_array+i)->hydr_dia=sqrt(4*(elem_array+i)->ave_area/PI);
		}
        
//round_with_ribs
        else if((elem_array+i)->elem_type==11)
        {
            fscanf(elem_file,"%lf%lf%lf",&(elem_array+i)->length,&(elem_array+i)->flow_angle,&(elem_array+i)->resist_factor);
            fscanf(elem_file,"%lf%lf%lf",&(elem_array+i)->diameter,&(elem_array+i)->rib_height,&(elem_array+i)->rib_space);
            fscanf(elem_file,"%lf%lf%lf%lf",&(elem_array+i)->inlet_area,&(elem_array+i)->inlet_radius,&(elem_array+i)->outlet_area,&(elem_array+i)->outlet_radius);
            fscanf(elem_file,"%lf%lf%lf%lf%lf",&(elem_array+i)->omega,&(elem_array+i)->krot,&(elem_array+i)->ht_coeff,&(elem_array+i)->ht_area,&(elem_array+i)->ht_temp);
            (elem_array+i)->ave_area=0.5*((elem_array+i)->inlet_area+(elem_array+i)->outlet_area);
            (elem_array+i)->hydr_dia=(elem_array+i)->diameter;
            (elem_array+i)->cor_coeff=1;
        }
//rect_with_ribs
        else if((elem_array+i)->elem_type==12)
        {
            fscanf(elem_file,"%lf%lf%lf",&(elem_array+i)->length,&(elem_array+i)->flow_angle,&(elem_array+i)->resist_factor);
            fscanf(elem_file,"%lf%lf%lf%lf",&(elem_array+i)->section_length,&(elem_array+i)->section_width,&(elem_array+i)->rib_height,&(elem_array+i)->rib_space);
            fscanf(elem_file,"%lf%lf%lf%lf",&(elem_array+i)->inlet_area,&(elem_array+i)->inlet_radius,&(elem_array+i)->outlet_area,&(elem_array+i)->outlet_radius);
            fscanf(elem_file,"%lf%lf%lf%lf%lf",&(elem_array+i)->omega,&(elem_array+i)->krot,&(elem_array+i)->ht_coeff,&(elem_array+i)->ht_area,&(elem_array+i)->ht_temp);
            (elem_array+i)->ave_area=0.5*((elem_array+i)->inlet_area+(elem_array+i)->outlet_area);
            width_length_ratio=(elem_array+i)->section_width/(elem_array+i)->section_length;
            (elem_array+i)->cor_coeff=1;
            (elem_array+i)->hydr_dia=2*(elem_array+i)->section_length*(elem_array+i)->section_width/((elem_array+i)->section_length+(elem_array+i)->section_width);
        }

//parallel_wall
        else if((elem_array+i)->elem_type==21)
        {
            fscanf(elem_file,"%lf%lf%lf%lf",&(elem_array+i)->length,&(elem_array+i)->flow_angle,&(elem_array+i)->resist_factor,&(elem_array+i)->roughness_factor);
            fscanf(elem_file,"%lf%lf%lf%lf%lf",&(elem_array+i)->wall_distance,&(elem_array+i)->inlet_area,&(elem_array+i)->inlet_radius,&(elem_array+i)->outlet_area,&(elem_array+i)->outlet_radius);
            fscanf(elem_file,"%lf%lf%lf%lf%lf",&(elem_array+i)->omega,&(elem_array+i)->krot,&(elem_array+i)->ht_coeff,&(elem_array+i)->ht_area,&(elem_array+i)->ht_temp);
            (elem_array+i)->ave_area=0.5*((elem_array+i)->inlet_area+(elem_array+i)->outlet_area);
            (elem_array+i)->hydr_dia=2*(elem_array+i)->wall_distance;
            (elem_array+i)->cor_coeff=1;
        }
//U_sharp_turn
        else if((elem_array+i)->elem_type==31)
        {
            fscanf(elem_file,"%lf%lf",&(elem_array+i)->section_length,&(elem_array+i)->section_width);
            fscanf(elem_file,"%lf%lf%lf",&(elem_array+i)->ht_coeff,&(elem_array+i)->ht_area,&(elem_array+i)->ht_temp);
            (elem_array+i)->ave_area=(elem_array+i)->section_length*(elem_array+i)->section_width;
            (elem_array+i)->hydr_dia=2*(elem_array+i)->section_length*(elem_array+i)->section_width/((elem_array+i)->section_length+(elem_array+i)->section_width);
            (elem_array+i)->length=(elem_array+i)->hydr_dia;
            (elem_array+i)->resist_factor=4.01;

            (elem_array+i)->inlet_area=(elem_array+i)->outlet_area=(elem_array+i)->ave_area;
        }
//U_round_turn
        else if((elem_array+i)->elem_type==32)
        {
            fscanf(elem_file,"%lf%lf",&(elem_array+i)->section_length,&(elem_array+i)->section_width);
            fscanf(elem_file,"%lf%lf%lf",&(elem_array+i)->ht_coeff,&(elem_array+i)->ht_area,&(elem_array+i)->ht_temp);
            (elem_array+i)->ave_area=(elem_array+i)->section_length*(elem_array+i)->section_width;
            (elem_array+i)->hydr_dia=2*(elem_array+i)->section_length*(elem_array+i)->section_width/((elem_array+i)->section_length+(elem_array+i)->section_width);
            (elem_array+i)->length=(elem_array+i)->hydr_dia;
            (elem_array+i)->resist_factor=4.51;

            (elem_array+i)->inlet_area=(elem_array+i)->outlet_area=(elem_array+i)->ave_area;
        }
//Stagnation_point_impingement
        else if((elem_array+i)->elem_type==41)
        {
            fscanf(elem_file,"%lf%lf%lf%lf",&(elem_array+i)->diameter,&(elem_array+i)->pitch,&(elem_array+i)->target_dia,&(elem_array+i)->target_dist);
            fscanf(elem_file,"%lf%lf%lf",&(elem_array+i)->ht_coeff,&(elem_array+i)->ht_area,&(elem_array+i)->ht_temp);
            (elem_array+i)->ave_area=PI*(elem_array+i)->diameter*(elem_array+i)->diameter/4;
            (elem_array+i)->hydr_dia=(elem_array+i)->diameter;
            (elem_array+i)->length=(elem_array+i)->hydr_dia;
            (elem_array+i)->resist_factor=3.7;
            (elem_array+i)->inlet_area=(elem_array+i)->outlet_area=(elem_array+i)->ave_area;
            //
            //
            //
        }
//head_area_impingement
        else if((elem_array+i)->elem_type==42)
        {
            fscanf(elem_file,"%lf%lf%lf%lf",&(elem_array+i)->diameter,&(elem_array+i)->pitch,&(elem_array+i)->target_dia,&(elem_array+i)->target_dist);
            fscanf(elem_file,"%lf%lf%lf",&(elem_array+i)->ht_coeff,&(elem_array+i)->ht_area,&(elem_array+i)->ht_temp);
            (elem_array+i)->flow_angle=0;
            (elem_array+i)->ave_area=PI*(elem_array+i)->diameter*(elem_array+i)->diameter/4;
            (elem_array+i)->hydr_dia=(elem_array+i)->diameter;
            (elem_array+i)->length=(elem_array+i)->hydr_dia;
            (elem_array+i)->resist_factor=3.7;

            (elem_array+i)->inlet_area=(elem_array+i)->outlet_area=(elem_array+i)->ave_area;
            //
            //
            //
        }
//body_impingement
        else if((elem_array+i)->elem_type==43)
        {
            fscanf(elem_file,"%lf%lf%lf%lf",&(elem_array+i)->diameter,&(elem_array+i)->pitch,&(elem_array+i)->target_dia,&(elem_array+i)->target_dist);
            fscanf(elem_file,"%lf%lf%lf",&(elem_array+i)->ht_coeff,&(elem_array+i)->ht_area,&(elem_array+i)->ht_temp);
            (elem_array+i)->ave_area=PI*(elem_array+i)->diameter*(elem_array+i)->diameter/4;
            (elem_array+i)->hydr_dia=(elem_array+i)->diameter;
            (elem_array+i)->length=(elem_array+i)->hydr_dia;
            (elem_array+i)->resist_factor=3.7;
            (elem_array+i)->inlet_area=(elem_array+i)->outlet_area=(elem_array+i)->ave_area;

            //
            //
            //
        }
//pin_fin_ribs        
        else if((elem_array+i)->elem_type==51)
        {
            fscanf(elem_file,"%lf%lf%lf%lf%lf",&(elem_array+i)->length,&(elem_array+i)->diameter,&(elem_array+i)->pitch,&(elem_array+i)->wall_distance,&(elem_array+i)->pinfin_row);
            fscanf(elem_file,"%lf%lf%lf",&(elem_array+i)->ht_coeff,&(elem_array+i)->ht_area,&(elem_array+i)->ht_temp);
            (elem_array+i)->ave_area=((elem_array+i)->pitch-(elem_array+i)->diameter)*(elem_array+i)->wall_distance;
            (elem_array+i)->hydr_dia=4*(elem_array+i)->ave_area*(elem_array+i)->length/(elem_array+i)->ht_area;
            (elem_array+i)->length=(elem_array+i)->hydr_dia;
            (elem_array+i)->inlet_area=(elem_array+i)->outlet_area=(elem_array+i)->ave_area;
        }
		else if((elem_array+i)->elem_type==61)
		{
			fscanf(elem_file,"%lf%lf%lf",&(elem_array+i)->resist_factor,&(elem_array+i)->inlet_area,&(elem_array+i)->outlet_area);
			if((elem_array+i)->outlet_area<1e-9)
			{
				return -1;
			}
			else if(fabs((elem_array+i)->inlet_area)<1e-9&&fabs((elem_array+i)->outlet_area)<1e-9)
			{
				return -1;
			}
			else if((elem_array+i)->inlet_area>1e-9&&(elem_array+i)->inlet_area<(elem_array+i)->outlet_area)
			{
				return -1;
			}
			(elem_array+i)->length=(elem_array+i)->hydr_dia=2*pow((elem_array+i)->outlet_area/PI,0.5);
			(elem_array+i)->ave_area=(elem_array+i)->outlet_area;
		}
		
		else if((elem_array+i)->elem_type==62)
		{
			fscanf(elem_file,"%lf%lf%lf",&(elem_array+i)->resist_factor,&(elem_array+i)->inlet_area,&(elem_array+i)->outlet_area);
			if((elem_array+i)->inlet_area<1e-9)
			{
				return -1;
			}
			else if(fabs((elem_array+i)->inlet_area)<1e-9&&fabs((elem_array+i)->outlet_area)<1e-9)
			{
				return -1;
			}
			else if((elem_array+i)->outlet_area>1e-9&&(elem_array+i)->inlet_area>(elem_array+i)->outlet_area)
			{
				return -1;
			}
			(elem_array+i)->length=(elem_array+i)->hydr_dia=2*pow((elem_array+i)->inlet_area/PI,0.5);
			(elem_array+i)->ave_area=(elem_array+i)->outlet_area;
		}
		
        if(p_setup->senergy_flag==0)
        {
            (elem_array+i)->ave_temp=p_setup->sref_temp;
        }
        if(p_setup->sfriction_flag==0)
        {
            (elem_array+i)->resist_factor=0;
        }
        i++;
    }
    while(!feof(elem_file)&&i<p_setup->selem_num_max);
    return 0;
}

void create_pipe_index(int* index_table,int node_num_max,elem* elem_array,int elem_num_max)
{
    int i,j;
    int i_index,j_index;
    for(i=0; i<node_num_max; i++)
    {
        for(j=0; j<node_num_max; j++)
        {
            *(index_table+i*node_num_max+j)=-1;
        }
    }
    for(i=0; i<elem_num_max; i++)
    {
        i_index=(elem_array+i)->inlet_node_num;
        j_index=(elem_array+i)->outlet_node_num;
        *(index_table+i_index*node_num_max+j_index)=i;
    }
//    for(i=0;i<node_num_max;i++)
 //  {
 //   	for(j=0;j<node_num_max;j++)
 //   	{
//	    	printf("%d,",*(index_table+i*node_num_max+j));
//	    }
//	    printf("\n");
 //   }
}

int create_bound_array(bound* bound_array,FILE* bound_file,int bound_num_max)
{
    int i=0;
    do
    {
        fscanf(bound_file,"%d%d%lf%lf",&(bound_array+i)->node_num,&(bound_array+i)->flow_bound_type,&(bound_array+i)->flow_bound_value,&(bound_array+i)->bound_temp);
        if((bound_array+i)->flow_bound_type!=1&&(bound_array+i)->flow_bound_type!=0)
        {
            return -1;
        }
        i++;
    }
    while(!feof(bound_file)&&i<bound_num_max);
    return 0;
//	for(i=0;i<bound_num_max;i++){
//		printf("%d\t%d\t%.2lf\t%.2lf\n",(bound_array+i)->node_num,(bound_array+i)->bound_type,(bound_array+i)->bound_pressure,(bound_array+i)->bound_density);
//	}
}

void create_node_array(node* node_array, elem* elem_array,bound* bound_array,int* index_table,int* node_not_pbound,setup* p_setup)
{
    int i,j;
    int n_n_b=0;
    int n_n_t_b=0;
    double temp_density;
    double no_use;
    for(i=0; i<p_setup->snode_num_max; i++)
    {
    	(node_array+i)->node_num=-1;
		(node_array+i)->node_temp=0;
		(node_array+i)->node_pressure=0;
		(node_array+i)->node_density=1;
		(node_array+i)->node_radius=0;
		(node_array+i)->node_omega=0;
		(node_array+i)->node_flow_rate=0;
		(node_array+i)->node_enth=0;
		(node_array+i)->flow_rate_bound_flag=-1;
		(node_array+i)->flow_bound_flag=-1;
		(node_array+i)->bound_index_pres=-1;
		(node_array+i)->bound_index_temp=-1;
		(node_array+i)->table_index_pres=-1;
		(node_array+i)->table_index_temp=-1;	
    	
        (node_array+i)->node_num=i;
        for(j=0; j<p_setup->sbound_num_max; j++)
        {
            if((bound_array+j)->node_num==i&&(bound_array+j)->flow_bound_type==0)
            {
                (node_array+i)->flow_rate_bound_flag=0;
                (node_array+i)->node_flow_rate=0;
                (node_array+i)->bound_index_pres=j;
                (node_array+i)->node_pressure=(bound_array+j)->flow_bound_value;
                (node_array+i)->flow_bound_flag=1;
                if(p_setup->senergy_flag==0)
                {
                    (node_array+i)->node_temp=p_setup->sref_temp;
                    (node_array+i)->bound_index_temp=-1;
                }
                else if(p_setup->senergy_flag==1&&(bound_array+j)->bound_temp!=-1)
                {
                    (node_array+i)->node_temp=(bound_array+j)->bound_temp;
                    (node_array+i)->bound_index_temp=j;
//         			steam_prop_calc((node_array+i)->node_pressure,(node_array+i)->node_temp,&(node_array+i)->node_density,&no_use);
                }
                else if(p_setup->senergy_flag==1&&(bound_array+j)->bound_temp==-1)
                {
                    (node_array+i)->node_temp=-1;
                    (node_array+i)->bound_index_temp=-1;
                }
                break;
            }
            else if((bound_array+j)->node_num==i&&(bound_array+j)->flow_bound_type==1)
            {
                (node_array+i)->flow_rate_bound_flag=1;
                (node_array+i)->node_flow_rate=(bound_array+j)->flow_bound_value;
                (node_array+i)->bound_index_pres=-1;
                (node_array+i)->node_pressure=-1;
                (node_array+i)->flow_bound_flag=1;
                if(p_setup->senergy_flag==0)
                {
                    (node_array+i)->node_temp=p_setup->sref_temp;
                    (node_array+i)->bound_index_temp=-1;
                }
                else
                {
                    (node_array+i)->node_temp=(bound_array+j)->bound_temp;
                    (node_array+i)->bound_index_temp=j;
                }
                break;
            }
            else
            {
                (node_array+i)->flow_rate_bound_flag=0;
                (node_array+i)->flow_bound_flag=0;
                (node_array+i)->node_flow_rate=0;
                (node_array+i)->node_pressure=-1;
                (node_array+i)->bound_index_pres=-1;
                (node_array+i)->bound_index_temp=-1;
                if(p_setup->senergy_flag==0)
                {
                    (node_array+i)->node_temp=p_setup->sref_temp;
                }
                else
                {
                    (node_array+i)->node_temp=-1;
                }
            }
        }
//        if(p_setup->sfluid_type==0&&p_setup->senergy_flag==0)
//        {
//            (node_array+i)->node_density=p_setup->WATER_DENS_23;
//        }
//        else
//        {
//            (node_array+i)->node_density=ideal_gas((node_array+i)->node_pressure,(node_array+i)->node_temp);

//        }
        for(j=0; j<p_setup->snode_num_max; j++)
        {
            if(*(index_table+i*p_setup->snode_num_max+j)!=-1)
            {
                (node_array+i)->node_radius=(elem_array+*(index_table+i*p_setup->snode_num_max+j))->inlet_radius;
                (node_array+i)->node_omega=(elem_array+*(index_table+i*p_setup->snode_num_max+j))->omega;
                break;
            }
        }
        for(j=0; j<p_setup->snode_num_max; j++)
        {
            if(*(index_table+j*p_setup->snode_num_max+i)!=-1)
            {
                (node_array+i)->node_radius=(elem_array+*(index_table+j*p_setup->snode_num_max+i))->outlet_radius;
                (node_array+i)->node_omega=(elem_array+*(index_table+j*p_setup->snode_num_max+i))->omega;
//				printf("%d\t%d\t%d\t%d\t%lf\n",(node_array+i)->node_num,i,j,*(index_table+j*node_num_max+i),(elem_array+*(index_table+j*node_num_max+i))->omega);
            }
        }
        if((node_array+i)->bound_index_pres==-1)
        {
            n_n_b++;
        }
//		printf("%d\t%d\t%d\t%lf\n",(node_array+i)->node_num,(node_array+i)->bound_index_pres,(node_array+i)->flow_rate_bound_flag,(node_array+i)->node_flow_rate);
    }

    *node_not_pbound=n_n_b;
}

void create_node_not_pbound_index(int* node_not_pbound_index,node* node_array,int node_num_max,int node_not_pbound)
{
    int i;
    int j=0;

    for(i=0; i<node_num_max; i++)
    {
        if((node_array+i)->bound_index_pres!=-1)
        {
            (node_array+i)->table_index_pres=-1;
            continue;
        }
        else
        {
            *(node_not_pbound_index+j)=i;
            (node_array+i)->table_index_pres=j;
            j++;
        }
    }
//    for(i=0; i<node_not_pbound; i++)
//    {
//        printf("%d,",*(node_not_pbound_index+i));
//    }
//	for(i=0;i<node_num_max;i++){
//		printf("%d\t%d\n",(node_array+i)->node_num,(node_array+i)->table_index_pres);
//	}
}

int calc_temp_bound(bound* bound_array, int* index_table, node* node_array, elem* elem_array, setup* p_setup)
{
	int i,j;
	int connect_elem;
	int node_num_max, bound_num_max;
	int bound_nnum=-1;
	int flow_to_elem, flow_from_elem;
	int nntb=0;
//	double no_use,no_use2;
	
	node_num_max=p_setup->snode_num_max;
	bound_num_max=p_setup->sbound_num_max;
	
	for(i=0;i<bound_num_max;i++)
	{
		bound_nnum=(bound_array+i)->node_num;
		if((bound_array+i)->flow_bound_type==1&&(bound_array+i)->flow_bound_value>0)
		{
			(node_array+bound_nnum)->node_temp=(bound_array+i)->bound_temp;
			(node_array+bound_nnum)->bound_index_temp=i;
//			steam_prop_calc((node_array+i)->node_pressure,(node_array+i)->node_temp,&(node_array+i)->node_density,&no_use,&no_use2);
			continue;	
		}
		else if((bound_array+i)->flow_bound_type==1&&(bound_array+i)->flow_bound_value<0)
		{
			(node_array+bound_nnum)->node_temp=-1;
			(node_array+bound_nnum)->bound_index_temp=-1;
			continue;
		}
		else
		{
			flow_to_elem=-1;
			flow_from_elem=-1;
			for(j=0;j<node_num_max;j++)
			{
				if(*(index_table+node_num_max*bound_nnum+j)!=-1)
				{
					flow_to_elem=*(index_table+node_num_max*bound_nnum+j);
					break;
				}
			}
			if(flow_to_elem==-1)
			{
				for(j=0;j<node_num_max;j++)
				{
					if(*(index_table+j*node_num_max+bound_nnum)!=-1)
					{
						flow_from_elem=*(index_table+j*node_num_max+bound_nnum);
						break;
					}
				}
			}
		}
		if(flow_to_elem!=-1&&(elem_array+flow_to_elem)->flow_rate>0)
		{
			
			(node_array+bound_nnum)->node_temp=(bound_array+i)->bound_temp;
			(node_array+bound_nnum)->bound_index_temp=i;
			continue;	
		}
		else if(flow_to_elem!=-1&&(elem_array+flow_to_elem)->flow_rate<0)
		{
			(node_array+bound_nnum)->node_temp=-1;
			(node_array+bound_nnum)->bound_index_temp=-1;
			continue;
		}
		else if(flow_from_elem!=-1&&(elem_array+flow_from_elem)->flow_rate>0)
		{
			(node_array+bound_nnum)->node_temp=-1;
			(node_array+bound_nnum)->bound_index_temp=-1;
			continue;	
		}
		else if(flow_from_elem!=-1&&(elem_array+flow_from_elem)->flow_rate<0)
		{
			(node_array+bound_nnum)->node_temp=(bound_array+i)->bound_temp;
			(node_array+bound_nnum)->bound_index_temp=i;
			continue;
		}			
	}
	
	for(i=0;i<node_num_max;i++)
	{
		if((node_array+i)->bound_index_temp==-1)
		{
			nntb++;
		}
	}
	return nntb;
}

void read_node_temp(double* prev_temp, node* node_array, int node_num_max)
{
	int i;
	for(i=0;i<node_num_max;i++)
	{
		*(prev_temp+i)=(node_array+i)->node_temp;
	}
}

void create_node_not_tbound_index(int* node_not_tbound_index,node* node_array,int node_num_max,int node_not_tbound)
{
    int i;
    int j=0;

    for(i=0; i<node_num_max; i++)
    {
        if((node_array+i)->bound_index_temp!=-1)
        {
            (node_array+i)->table_index_temp=-1;
            continue;
        }
        else
        {
            *(node_not_tbound_index+j)=i;
            (node_array+i)->table_index_temp=j;
            j++;
        }
    }
//	for(i=0;i<node_num_max;i++){
//		printf("%d\t%d\n",(node_array+i)->node_num,(node_array+i)->table_index_temp);
//	}
//    return index;
}


void initial_coeff(double* coeff_array,int coeff_num)
{
    int i;
    for(i=0; i<coeff_num; i++)
    {
        *(coeff_array+i)=0;
    }
}

void create_energy_equations_new(double** p_energy_eqns,double* right_hand_array,elem* elem_array,node* node_array,int* index_table,int node_not_tbound,int* node_not_tbound_index, int node_num_max)
{
	int i,node_index,k, upstream_node,linked_pipe_index;
	double coeff_sum,constant_sum,constant_temp,fric_coeff;
	double downstr_flow_rate;
	
	for(i=0;i<node_not_tbound;i++)
	{
		coeff_sum=0;
		constant_sum=0;
		downstr_flow_rate=0;
		node_index=*(node_not_tbound_index+i);
		initial_coeff(*(p_energy_eqns+i),node_not_tbound);
		for(k=0;k<node_num_max;k++)
		{
			linked_pipe_index=*(index_table+k*node_num_max+node_index);
			if(linked_pipe_index!=-1&&(elem_array+linked_pipe_index)->flow_rate>0)
			{
				upstream_node=(elem_array+linked_pipe_index)->inlet_node_num;
				fric_coeff=(elem_array+linked_pipe_index)->resist_factor*(elem_array+linked_pipe_index)->length/(2*(elem_array+linked_pipe_index)->ave_dens*(elem_array+linked_pipe_index)->hydr_dia*(elem_array+linked_pipe_index)->ave_area*(elem_array+linked_pipe_index)->ave_area);
				constant_temp=((node_array+upstream_node)->node_pressure-(node_array+node_index)->node_pressure-fric_coeff*pow((elem_array+linked_pipe_index)->flow_rate,2))*(elem_array+linked_pipe_index)->flow_rate/(elem_array+linked_pipe_index)->ave_dens;
				constant_sum+=constant_temp;
				constant_temp=(elem_array+linked_pipe_index)->ht_coeff*(elem_array+linked_pipe_index)->ht_area*((elem_array+linked_pipe_index)->ht_temp-0.5*((node_array+upstream_node)->node_temp+(node_array+node_index)->node_temp));
				constant_sum-=constant_temp;
				downstr_flow_rate+=(elem_array+linked_pipe_index)->flow_rate;
				if((node_array+upstream_node)->bound_index_temp!=-1)
				{
					constant_sum-=(elem_array+linked_pipe_index)->flow_rate*(node_array+upstream_node)->node_enth;
				}
				else
				{
					*(*(p_energy_eqns+i)+(node_array+upstream_node)->table_index_temp)=(elem_array+linked_pipe_index)->flow_rate;
				}
			}
			else if(linked_pipe_index!=-1&&(elem_array+linked_pipe_index)->flow_rate<0)
			{
				coeff_sum-=fabs((elem_array+linked_pipe_index)->flow_rate);
			}
		}
		for(k=0;k<node_num_max;k++)
		{
			linked_pipe_index=*(index_table+node_index*node_num_max+k);
			if(linked_pipe_index!=-1&&(elem_array+linked_pipe_index)->flow_rate<0)
			{
				upstream_node=(elem_array+linked_pipe_index)->outlet_node_num;
				fric_coeff=(elem_array+linked_pipe_index)->resist_factor*(elem_array+linked_pipe_index)->length/(2*(elem_array+linked_pipe_index)->ave_dens*(elem_array+linked_pipe_index)->hydr_dia*(elem_array+linked_pipe_index)->ave_area*(elem_array+linked_pipe_index)->ave_area);
				constant_temp=((node_array+upstream_node)->node_pressure-(node_array+node_index)->node_pressure-fric_coeff*pow((elem_array+linked_pipe_index)->flow_rate,2))*(-(elem_array+linked_pipe_index)->flow_rate)/(elem_array+linked_pipe_index)->ave_dens;
				constant_sum+=constant_temp;
				constant_temp=(elem_array+linked_pipe_index)->ht_coeff*(elem_array+linked_pipe_index)->ht_area*((elem_array+linked_pipe_index)->ht_temp-0.5*((node_array+upstream_node)->node_temp+(node_array+node_index)->node_temp));
				constant_sum-=constant_temp;
				downstr_flow_rate-=(elem_array+linked_pipe_index)->flow_rate;
				
				if((node_array+upstream_node)->bound_index_temp!=-1)
				{
					constant_sum-=fabs((elem_array+linked_pipe_index)->flow_rate)*(node_array+node_index)->node_enth;
				}
				else
				{
					*(*(p_energy_eqns+i)+(node_array+upstream_node)->table_index_temp)=-(elem_array+linked_pipe_index)->flow_rate;
				}
			}
			else if(linked_pipe_index!=-1&&(elem_array+linked_pipe_index)->flow_rate>0)
			{
				coeff_sum-=(elem_array+linked_pipe_index)->flow_rate;
			}
		}
		if((node_array+node_index)->flow_bound_flag==1)
		{
			*(*(p_energy_eqns+i)+i)=-downstr_flow_rate;
		}
		else
		{
			*(*(p_energy_eqns+i)+i)=coeff_sum;
		}	
		*(right_hand_array+i)=constant_sum;
	}	
}


/*void update_node_temp(node* node_array,double* enth_solu,int* node_not_tbound_index,int node_not_tbound)
{
    int i;
    double new_temp;
    
    for(i=0; i<node_not_tbound; i++)
    {
        (node_array+*(node_not_tbound_index+i))->node_temp=*(enth_solu+i);
    }
}*/

int calc_node_enth(node* node_array, int node_num_max, int fluid_type)
{
	int i,flag=0;
	double no_use1, no_use2, no_use3, no_use4;
	for(i=0;i<node_num_max;i++)
	{
		if(fluid_type==0||fluid_type==1)
		{
			flag=steam_prop_calc((node_array+i)->node_pressure,(node_array+i)->node_temp,&no_use1,&no_use2,&no_use3,&(node_array+i)->node_enth,&no_use4);
			if(flag<0)
			{
				printf("ERRRRRRRROR\n");
				return -1;
			}	
		}
		else if(fluid_type==2)
		{
			(node_array+i)->node_enth=spe_enth_air((node_array+i)->node_temp);
		}
		
		
//		printf("%d,%lf\n",i,(node_array+i)->node_enth);
	}
	return 0;
}

void update_node_enth(node* node_array,double* enth_solu,int* node_not_tbound_index,int node_not_tbound)
{
	int i;    
    for(i=0; i<node_not_tbound; i++)
    {
        (node_array+*(node_not_tbound_index+i))->node_enth=*(enth_solu+i);
    }
}

int update_node_temp_new(node* node_array,int* node_not_tbound_index,int node_not_tbound, int fluid_type)
{
	int i;
	int flag=0;
	for(i=0;i<node_not_tbound;i++)
	{
		if(fluid_type==0||fluid_type==1)
		{
			flag=steam_temp_calc_ph(&(node_array+*(node_not_tbound_index+i))->node_temp,(node_array+*(node_not_tbound_index+i))->node_pressure,(node_array+*(node_not_tbound_index+i))->node_enth);
			if(flag<0)
			{
				return -1;
			}
		}
		else if(fluid_type==2)
		{
			flag=air_temp_calc_h(&(node_array+*(node_not_tbound_index+i))->node_temp,(node_array+*(node_not_tbound_index+i))->node_enth);
			if(flag<0)
			{
				return -1;
			}
		}	
	}
	return 0;
}

void update_elem_temp(elem* elem_array, node* node_array,int elem_num_max)
{
    int i;
    double inlet_temp,outlet_temp;
    for(i=0; i<elem_num_max; i++)
    {
        inlet_temp=(node_array+(elem_array+i)->inlet_node_num)->node_temp;
        outlet_temp=(node_array+(elem_array+i)->outlet_node_num)->node_temp;
        (elem_array+i)->ave_temp=0.5*(inlet_temp+outlet_temp);
    }
}

int update_flow_thermal_coeffs(elem* p_elem, setup* p_setups)
{
    double viscous,ave_velo,resist_factor_temp;
    double nusselt_num;
    double* p_resist_factor_temp=&resist_factor_temp;
    double nw_pranum,nw_thcond,nw_dynamic_viscous;
    air_prop near_wall;
    int iter_flag;
//    air_prop near_wall;

//for C3X calculation ONLY!
    double cr_num[10]= {1.118,1.118,1.118,1.118,1.118,1.118,1.118,1.056,1.056,1.025};
//for C3X calculation ONLY!

    if(p_setups->sfluid_type==0)
    {
		iter_flag=colebrook_iter(&resist_factor_temp,p_elem->roughness_factor,p_elem->re_num,0.4,p_setups->max_fric_iter_times);
//        printf("\t\t%d,,%lf\n",iter_flag,resist_factor_temp);
		if(iter_flag==-1)
        {
            return -1;
        }
//        printf("\t\t%lf\n",resist_factor_temp);
        p_elem->resist_factor=1/pow(resist_factor_temp,2);
        
        return 0;
//need extension
//heat transfer coeff calculation.
    }
    else if(p_setups->sfluid_type==1||p_setups->sfluid_type==2)
    {
        if(p_elem->elem_type<11)
        {
            if(karman_nikuradse_iter(p_resist_factor_temp,p_elem->re_num,0.01,p_setups->max_fric_iter_times)==-1)
            {
                return -1;
            }
            p_elem->resist_factor=pow(resist_factor_temp,2)*p_elem->cor_coeff;
            if(p_elem->elem_type==6)
            {
                p_elem->resist_factor=p_elem->resist_factor+(1.5*p_elem->diameter/p_elem->length);
            }
        }

        else if(p_elem->elem_type==11)
        {
            p_elem->resist_factor=round_rib_fric(p_elem->hydr_dia,p_elem->rib_height,p_elem->rib_space);
        }

        else if(p_elem->elem_type==12)
        {
            p_elem->resist_factor=rec_rib_fric(p_elem->re_num,p_elem->hydr_dia,p_elem->rib_height,p_elem->rib_space,p_elem->section_width/p_elem->section_length);
        }

        else if(p_elem->elem_type==21)
        {
            if(karman_nikuradse_iter(p_resist_factor_temp,p_elem->re_num,0.01,p_setups->max_fric_iter_times)==-1)
            {
                return -1;
            }
            p_elem->resist_factor=pow(resist_factor_temp,2)*p_elem->cor_coeff;
        }

        else if(p_elem->elem_type==51)
        {
            p_elem->re_num=fabs(p_elem->ave_dens*ave_velo*p_elem->diameter/viscous);
            p_elem->resist_factor=pinfin_fric(p_elem->re_num,p_elem->pitch,p_elem->diameter,p_elem->pinfin_row);
        }
        
        else if(p_elem->elem_type==61)
        {
        	if(p_elem->inlet_area==-1)
        	{
	        	p_elem->resist_factor=0.5;
	        }
	        else
	        {
        		p_elem->resist_factor=0.5*(1-p_elem->outlet_area/p_elem->inlet_area);
        	}
        }
        
        else if(p_elem->elem_type==62)
        {
        	if(p_elem->outlet_area==-1)
        	{
	        	p_elem->resist_factor=1;
	        }
	        else
	        {
        		p_elem->resist_factor=(1-p_elem->inlet_area/p_elem->outlet_area)*(1-p_elem->inlet_area/p_elem->outlet_area);
        	}
       	}


        if(p_setups->senergy_flag==1&&p_setups->shtc_flag==1)
        {
            if(p_elem->elem_type==1||p_elem->elem_type==3||p_elem->elem_type==4||p_elem->elem_type==5||p_elem->elem_type==31||p_elem->elem_type==32)
            {
                if(p_setups->scorrelation_flag==0)
                {
                    nusselt_num=C3X_nu(p_elem->re_num,p_elem->fluid_pranum,cr_num[p_elem->elem_num]);
                }
                else if(p_setups->scorrelation_flag==1)
                {
                    nusselt_num=dittus_boelter(p_elem->re_num,p_elem->fluid_pranum);
                }
                else if(p_setups->scorrelation_flag==2)
                {
                    if(p_elem->re_num<10000)
                    {
                        nusselt_num=gnielinski(p_elem->resist_factor,p_elem->re_num,p_elem->fluid_pranum);
                    }
                    else if(p_elem->re_num>10000||p_elem->re_num==10000)
                    {
                        nusselt_num=petukhov_kirillov(p_elem->resist_factor,p_elem->re_num,p_elem->fluid_pranum);
                    }
                }
                p_elem->ht_coeff=nusselt_num*p_elem->fluid_thcond/p_elem->hydr_dia;
            }

            else if(p_elem->elem_type==2)
            {
                nusselt_num=rec_smooth_nu(p_elem->re_num,p_elem->fluid_pranum,p_elem->section_width/p_elem->section_length);
                p_elem->ht_coeff=nusselt_num*p_elem->fluid_thcond/p_elem->hydr_dia;
            }

            else if(p_elem->elem_type==11)
            {
                nusselt_num=round_rib_nu(p_elem->resist_factor,p_elem->re_num,p_elem->hydr_dia,p_elem->rib_height,p_elem->fluid_pranum,p_elem->rib_space);
                p_elem->ht_coeff=nusselt_num*p_elem->fluid_thcond/p_elem->hydr_dia;
            }

            else if(p_elem->elem_type==12)
            {
                nusselt_num=rec_rib_nu(p_elem->resist_factor,p_elem->re_num,p_elem->fluid_pranum,p_elem->hydr_dia,p_elem->rib_height,p_elem->rib_space,p_elem->section_width/p_elem->section_length);
                p_elem->ht_coeff=nusselt_num*p_elem->fluid_thcond/p_elem->hydr_dia;
            }

            else if(p_elem->elem_type==21)
            {
                cal_property(&near_wall,p_elem->ht_temp);
                nusselt_num=para_wall_nu(p_elem->re_num,p_elem->fluid_pranum,p_elem->wall_distance,p_elem->length,p_elem->fluid_visc,near_wall.dynamic_viscous);
                p_elem->ht_coeff=nusselt_num*p_elem->fluid_thcond/p_elem->hydr_dia;
            }

            else if(p_elem->elem_type==41)
            {
                nusselt_num=impinge_stag_nu(p_elem->re_num,p_elem->diameter,p_elem->pitch,p_elem->target_dia,p_elem->target_dist);
                p_elem->ht_coeff=nusselt_num*p_elem->fluid_thcond/p_elem->hydr_dia;
            }

            else if(p_elem->elem_type==42)
            {
                nusselt_num=impinge_head_nu(p_elem->re_num,p_elem->diameter,p_elem->pitch,p_elem->target_dia,p_elem->target_dist);
                p_elem->ht_coeff=nusselt_num*p_elem->fluid_thcond/p_elem->hydr_dia;
            }
            else if(p_elem->elem_type==43)
            {
                p_elem->ht_coeff=impinge_body_h(ave_velo,p_elem->pitch,viscous,p_elem->ave_dens);
            }

            else if(p_elem->elem_type==51)
            {
                nusselt_num=pinfin_nu(p_elem->re_num,p_elem->fluid_pranum);
                p_elem->ht_coeff=nusselt_num*p_elem->fluid_thcond/p_elem->hydr_dia;
            }
        }
        return 0;
    }
}

int update_node_vars(node* node_array,double* prev_solu,int* node_not_pbound_index,int node_not_pbound,setup* p_setup)
{
    int i,prop_calc_flag;
    double no_use,no_use2,no_use3,no_use4;
    for(i=0; i<node_not_pbound; i++)
    {
        (node_array+*(node_not_pbound_index+i))->node_pressure=*(prev_solu+i);
//		printf("%.4lf\t%.4lf\t%.12lf\n",*(prev_solu+i),(node_array+*(node_not_pbound_index+i))->node_temp,(node_array+*(node_not_pbound_index+i))->node_density);
    }
    
   	for(i=0;i<p_setup->snode_num_max;i++)
   	{
    	if(p_setup->sfluid_type==0||p_setup->sfluid_type==1)
    	{
    		prop_calc_flag=steam_prop_calc((node_array+i)->node_pressure,(node_array+i)->node_temp,&(node_array+i)->node_density,&no_use,&no_use2,&no_use3,&no_use4);
    		if(prop_calc_flag==-1)
    		{
		    	return -1;
		    }
    	}
    	else if(p_setup->sfluid_type==2)
    	{
    		(node_array+i)->node_density=ideal_gas((node_array+i)->node_pressure,(node_array+i)->node_temp);
    	}
    }
    return 0;
}

int update_elem_vars(elem* elem_array,double* prev_solu,node* node_array,int node_not_pbound,setup* p_setup)
{
    double inlet_pressure,outlet_pressure,ave_pressure,viscous,ave_velo,cp,cv,drdp;
    double no_use,no_use2,no_use3,no_use4;
    int i,elem_index,flag;
	air_prop prop;    

    if(p_setup->sfluid_type==0)
    {
        for(i=node_not_pbound; i<node_not_pbound+p_setup->selem_num_max; i++)
        {
            elem_index=i-node_not_pbound;
            (elem_array+elem_index)->flow_rate=*(prev_solu+i);
            inlet_pressure=(node_array+(elem_array+elem_index)->inlet_node_num)->node_pressure;
            outlet_pressure=(node_array+(elem_array+elem_index)->outlet_node_num)->node_pressure;
			ave_pressure=0.5*(inlet_pressure+outlet_pressure);
            
			if(steam_prop_calc(ave_pressure,(elem_array+elem_index)->ave_temp,&(elem_array+elem_index)->ave_dens,&(elem_array+elem_index)->fluid_cp,&no_use,&no_use2,&no_use3)==-1)
			{
				return -1;
			}
            
			viscous=water_viscous((elem_array+elem_index)->ave_temp);
			(elem_array+elem_index)->fluid_visc=viscous;
			water_thcond_pranum_calc(&(elem_array+elem_index)->fluid_thcond,&(elem_array+elem_index)->fluid_pranum,(elem_array+elem_index)->ave_temp);
//			printf("%lf,%lf,%lf\n",ave_velo,viscous,(elem_array+elem_index)->hydr_dia);
        	ave_velo=(elem_array+elem_index)->flow_rate/((elem_array+elem_index)->ave_dens*(elem_array+elem_index)->ave_area);
        	(elem_array+elem_index)->re_num=fabs((elem_array+elem_index)->ave_dens*ave_velo*(elem_array+elem_index)->hydr_dia/viscous);
        	
            if(p_setup->sfriction_flag==1&&p_setup->sfric_coeff_flag==1)
            {
            	if((elem_array+elem_index)->elem_type!=7&&(elem_array+elem_index)->elem_type!=8&&(elem_array+elem_index)->elem_type!=9)
            	{
//            		printf("\t\t\t%d,\n",(elem_array+elem_index)->elem_type);
	            	flag=update_flow_thermal_coeffs(elem_array+elem_index,p_setup);
	            	if(flag==-1)
                	{
                		printf(".....................");
                    	return -2;
                	}
	            }
            }
            if(p_setup->senergy_flag==1)
            {
                (elem_array+elem_index)->elem_specific_heat=(elem_array+elem_index)->fluid_cp;
            }
//            printf("\n\t\t\t%lf\n",(elem_array+elem_index)->resist_factor);
        }
    }
    else if(p_setup->sfluid_type==1||p_setup->sfluid_type==2)
    {
        for(i=node_not_pbound; i<node_not_pbound+p_setup->selem_num_max; i++)
        {
            elem_index=i-node_not_pbound;
/*            if((elem_array+elem_index)->elem_type==61&&*(prev_solu+i)<0)
            {
            	return -2;
            }
            else if((elem_array+elem_index)->elem_type==62&&*(prev_solu+i)<0)
            {
            	return -2;
            }*/
            (elem_array+elem_index)->flow_rate=*(prev_solu+i);
            inlet_pressure=(node_array+(elem_array+elem_index)->inlet_node_num)->node_pressure;
            outlet_pressure=(node_array+(elem_array+elem_index)->outlet_node_num)->node_pressure;
//			printf("%lf,%lf\n",(inlet_pressure+outlet_pressure)/2,(elem_array+elem_index)->ave_temp);
            if(p_setup->sfluid_type==1)
            {
            	if(steam_prop_calc((inlet_pressure+outlet_pressure)/2,(elem_array+elem_index)->ave_temp,&(elem_array+elem_index)->ave_dens,&(elem_array+elem_index)->fluid_cp,&cv,&no_use,&drdp)==-1)
            	{
	            	return -1;
	            }
	           	(elem_array+elem_index)->fluid_visc=steam_visc_calc((elem_array+elem_index)->ave_temp,(elem_array+elem_index)->ave_dens);
	           	(elem_array+elem_index)->fluid_thcond=steam_thcond_calc((elem_array+elem_index)->ave_temp,(elem_array+elem_index)->ave_dens,(elem_array+elem_index)->fluid_visc,(elem_array+elem_index)->fluid_cp,cv,drdp,(inlet_pressure+outlet_pressure)/2);
	            (elem_array+elem_index)->fluid_pranum=0.7;            	
            }
            else
            {
            	(elem_array+elem_index)->ave_dens=(inlet_pressure+outlet_pressure)/((2*(elem_array+elem_index)->ave_temp)*GAS_CONST_AIR);
            	cal_property(&prop,(elem_array+elem_index)->ave_temp);
            	viscous=prop.dynamic_viscous;
            	(elem_array+elem_index)->fluid_visc=viscous;
            	(elem_array+elem_index)->fluid_pranum=prop.prandtl_num;
            	(elem_array+elem_index)->fluid_thcond=prop.thermal_conduct;
            	(elem_array+elem_index)->fluid_cp=prop.specific_heat_p;
            	(elem_array+elem_index)->elem_specific_heat=prop.specific_heat_p;
        		ave_velo=(elem_array+elem_index)->flow_rate/((elem_array+elem_index)->ave_dens*(elem_array+elem_index)->ave_area);
        		(elem_array+elem_index)->re_num=fabs((elem_array+elem_index)->ave_dens*ave_velo*(elem_array+elem_index)->hydr_dia/viscous);
            }			
//			printf("%lf\n",(elem_array+elem_index)->ave_dens);
//            
//            viscous=prop.dynamic_viscous;
//        	ave_velo=(elem_array+elem_index)->flow_rate/((elem_array+elem_index)->ave_dens*(elem_array+elem_index)->ave_area);
//        	(elem_array+elem_index)->re_num=fabs((elem_array+elem_index)->ave_dens*ave_velo*(elem_array+elem_index)->hydr_dia/viscous);
            if(p_setup->sfriction_flag==1&&p_setup->sfric_coeff_flag==1)
            {
                if(update_flow_thermal_coeffs(elem_array+elem_index,p_setup)==-1)
                {
                    return -2;
                }
            }
        }
    }
    return 0;
}

double cal_coeff_m(elem* p_elem,node* node_array)
{
    double inlet_dens,outlet_dens,inlet_temp,outlet_temp,inlet_pres,outlet_pres;

    double inertia_coeff;
    double friction_coeff;

    inlet_dens=(node_array+p_elem->inlet_node_num)->node_density;
    outlet_dens=(node_array+p_elem->outlet_node_num)->node_density;
    
	if(p_elem->elem_type==7)
    {
    	if(p_elem->flow_rate>0)
		{
			return p_elem->param_2*p_elem->ave_area;
		} 
		else
		{
			return -(p_elem->param_2*p_elem->ave_area);
		}
    }
    
    else if(p_elem->elem_type==8)
    {
    	if(p_elem->flow_rate>0)
    	{
	    	return -(p_elem->param_1/p_elem->re_num+p_elem->param_2*(1+1/p_elem->hydr_dia))/(2*p_elem->ave_dens*p_elem->ave_area);
	    }
	    else
	    {
    		return (p_elem->param_1/p_elem->re_num+p_elem->param_2*(1+1/p_elem->hydr_dia))/(2*p_elem->ave_dens*p_elem->ave_area);
    	}
    }
    
    
    else if(p_elem->elem_type<31)
    {
        if(p_elem->flow_rate>0)
        {
        	inertia_coeff=1/(p_elem->ave_area*p_elem->ave_dens)-(p_elem->ave_area*p_elem->ave_dens)/(p_elem->outlet_area*p_elem->outlet_area*outlet_dens*outlet_dens);
//           inertia_coeff=p_elem->ave_area/(inlet_dens*p_elem->inlet_area*p_elem->inlet_area)-1/(p_elem->ave_dens*p_elem->ave_area);
//            inertia_coeff=p_elem->ave_area/(inlet_dens*p_elem->inlet_area*p_elem->inlet_area)-p_elem->ave_area/(outlet_dens*p_elem->outlet_area*p_elem->outlet_area);
//			inertia_coeff=0.5*(p_elem->ave_dens*p_elem->ave_area)*(1/(inlet_dens*inlet_dens*p_elem->inlet_area*p_elem->inlet_area)-1/(p_elem->ave_dens*p_elem->ave_dens*p_elem->ave_area*p_elem->ave_area));
//        	inertia_coeff=1/(inlet_dens*p_elem->ave_area)-1/(p_elem->ave_dens*p_elem->ave_area);
		}
        else
        {
        	inertia_coeff=p_elem->ave_area*p_elem->ave_dens/(inlet_dens*inlet_dens*p_elem->inlet_area*p_elem->inlet_area)-1/(p_elem->ave_dens*p_elem->ave_area);
//        	inertia_coeff=p_elem->ave_area/(inlet_dens*p_elem->inlet_area*p_elem->inlet_area)-p_elem->ave_area/(outlet_dens*p_elem->outlet_area*p_elem->outlet_area);
//            inertia_coeff=1/(p_elem->ave_area*p_elem->ave_dens)-(p_elem->ave_area)/(p_elem->outlet_area*p_elem->outlet_area*outlet_dens);
//            inertia_coeff=1/(p_elem->inlet_area*inlet_dens)-1/(p_elem->outlet_area*outlet_dens);
//        	inertia_coeff=0.5*(p_elem->ave_dens*p_elem->ave_area)*(1/(p_elem->ave_dens*p_elem->ave_dens*p_elem->ave_area*p_elem->ave_area)-1/(outlet_dens*outlet_dens*p_elem->outlet_area*p_elem->outlet_area));
//			inertia_coeff=1/(p_elem->ave_area*p_elem->ave_dens)-1/(p_elem->ave_area*outlet_dens);
		}
//		printf("\t\t%lf\n",inertia_coeff);
    }
    
    else if(p_elem->elem_type==61)
    {
    	if(p_elem->inlet_area==-1)
    	{
	    	inertia_coeff=-1/(p_elem->ave_dens*p_elem->ave_area);
	    }
	    else
	    {
    		inertia_coeff=1/(inlet_dens*p_elem->inlet_area)-1/(p_elem->ave_dens*p_elem->ave_area);
    	}    	
    }
    
    else if(p_elem->elem_type==62)
    {
    	if(p_elem->outlet_area==-1)
    	{
	    	inertia_coeff=-1/(p_elem->ave_dens*p_elem->ave_area);
	    }
	    else
	    {
    		inertia_coeff=1/(inlet_dens*p_elem->inlet_area)-1/(p_elem->ave_dens*p_elem->ave_area);
    	}    	
    }
    
    else
    {
        inertia_coeff=0;
    }
	
 	if(p_elem->elem_type==9)
    {
    	if(fabs(p_elem->param_1)>1e-6)
    	{
	    	friction_coeff=1/(2*(p_elem->ave_dens*p_elem->param_1*p_elem->param_1*p_elem->ave_area));
	    }
		else
		{
			friction_coeff=0;
		}	
    }
    else
    {
    	friction_coeff=p_elem->resist_factor*p_elem->length/(2*p_elem->ave_dens*p_elem->hydr_dia*p_elem->ave_area);
    }
    
    if(p_elem->flow_rate>0)
    {
//		printf("%lf\t%lf\t%lf\n",inertia_coeff,friction_coeff,inertia_coeff-friction_coeff);
        return	inertia_coeff-friction_coeff;
    }
    else
    {
//		printf("%lf\t%lf\t%lf\n",inertia_coeff,friction_coeff,inertia_coeff+friction_coeff);
        return  inertia_coeff+friction_coeff;
    }
}

void create_nonlinear_matrix(double** nonlinear_matrix,elem* elem_array,node* node_array,bound* bound_array,int* index_table,int elem_num_max,int node_num_max,int node_not_pbound)
{
    int i,j,k;
    int matrix_index=0;
    int elem_index;
    int inlet_node,outlet_node;
    double elem_coeff;
    int coeff_num=2*node_num_max+elem_num_max*3;
//    FILE* nm=fopen("..\\test_files\\nonlinear_matrix.txt","w");

    for(i=0; i<node_num_max; i++)
    {
        if((node_array+i)->bound_index_pres!=-1)
        {
            continue;
        }
        else
        {
            initial_coeff(*(nonlinear_matrix+matrix_index),coeff_num);
            for(k=0; k<node_num_max; k++)
            {
                if(*(index_table+i*node_num_max+k)!=-1)
                {
                    *(*(nonlinear_matrix+matrix_index)+2*node_num_max+*(index_table+i*node_num_max+k)*3+1)=-1;
                }
            }
            for(k=0; k<node_num_max; k++)
            {
                if(*(index_table+k*node_num_max+i)!=-1)
                {
                    *(*(nonlinear_matrix+matrix_index)+2*node_num_max+*(index_table+k*node_num_max+i)*3+1)=1;
                }
            }
            matrix_index++;
        }
    }

    for(i=matrix_index; i<node_not_pbound+elem_num_max; i++)
    {
        initial_coeff(*(nonlinear_matrix+i),coeff_num);
        elem_index=i-node_not_pbound;
        inlet_node=(elem_array+elem_index)->inlet_node_num;
        outlet_node=(elem_array+elem_index)->outlet_node_num;
		elem_coeff=cal_coeff_m(elem_array+elem_index,node_array);
//		printf("\n\t\t%d:::%lf",elem_index,elem_coeff);
//		printf(",%lf",(elem_array+elem_index)->ave_area);
//        printf("%.6lf\n",elem_coeff);
        if((node_array+inlet_node)->bound_index_pres==-1&&(node_array+outlet_node)->bound_index_pres==-1)
        {
            *(*(nonlinear_matrix+i)+inlet_node*2)=(elem_array+elem_index)->ave_area;
            *(*(nonlinear_matrix+i)+outlet_node*2)=-(elem_array+elem_index)->ave_area;
        }
        else if((node_array+inlet_node)->bound_index_pres!=-1&&(node_array+outlet_node)->bound_index_pres==-1)
        {
            *(*(nonlinear_matrix+i)+inlet_node*2+1)=(elem_array+elem_index)->ave_area*(node_array+inlet_node)->node_pressure;
            *(*(nonlinear_matrix+i)+outlet_node*2)=-(elem_array+elem_index)->ave_area;
        }
        else if((node_array+inlet_node)->bound_index_pres==-1&&(node_array+outlet_node)->bound_index_pres!=-1)
        {
            *(*(nonlinear_matrix+i)+inlet_node*2)=(elem_array+elem_index)->ave_area;
            *(*(nonlinear_matrix+i)+outlet_node*2+1)=-(elem_array+elem_index)->ave_area*(node_array+outlet_node)->node_pressure;
        }
        else
        {
            *(*(nonlinear_matrix+i)+inlet_node*2+1)=(elem_array+elem_index)->ave_area*(node_array+inlet_node)->node_pressure;
            *(*(nonlinear_matrix+i)+outlet_node*2+1)=-(elem_array+elem_index)->ave_area*(node_array+outlet_node)->node_pressure;
        }
        *(*(nonlinear_matrix+i)+2*node_num_max+3*elem_index)=elem_coeff;
        if((elem_array+elem_index)->elem_type==7)
        {
        	*(*(nonlinear_matrix+i)+2*node_num_max+3*elem_index+2)=(elem_array+elem_index)->param_1*(elem_array+elem_index)->ave_area;
        }
    }
//    printf("\n\n%d\n\n",node_not_pbound);
    //printf("%d\n",matrix_index);
/*  for(i=0;i<node_not_pbound+elem_num_max;i++){
  		for(j=0;j<coeff_num;j++){
	    	fprintf(nm,"%.5lf,",*(*(nonlinear_matrix+i)+j));
	    }
	    fprintf(nm,"\n");
	 }*/
}

void create_jacobi(double** jacobi_matrix,double** nonlinear_matrix,double* prev_solu,int* node_not_pbound_index,int elem_num_max,int node_num_max,int node_not_pbound)
{
    int i,j;
    int j_coeff_num=node_not_pbound+elem_num_max;

    for(i=0; i<node_not_pbound; i++)
    {
        initial_coeff(*(jacobi_matrix+i),j_coeff_num);
        for(j=0; j<node_not_pbound; j++)
        {
            *(*(jacobi_matrix+i)+j)=0;
        }
        for(j=node_not_pbound; j<j_coeff_num; j++)
        {
            *(*(jacobi_matrix+i)+j)=*(*(nonlinear_matrix+i)+node_num_max*2+(j-node_not_pbound)*3+1);
        }
    }
    for(i=node_not_pbound; i<j_coeff_num; i++)
    {
        initial_coeff(*(jacobi_matrix+i),j_coeff_num);
        for(j=0; j<node_not_pbound; j++)
        {
            *(*(jacobi_matrix+i)+j)=*(*(nonlinear_matrix+i)+*(node_not_pbound_index+j)*2);
        }
        for(j=node_not_pbound; j<j_coeff_num; j++)
        {
            *(*(jacobi_matrix+i)+j)=*(*(nonlinear_matrix+i)+node_num_max*2+(j-node_not_pbound)*3)*2*(*(prev_solu+i));
        }
    }
//	printf("\n");
//	for(i=0;i<j_coeff_num;i++){
//		for(j=0;j<j_coeff_num;j++){
//			printf("%lf,",*(*(jacobi_matrix+i)+j));
//		}
//		printf("\n");
//	}
}

void initialization_node_temp(node* node_array, double ini_temp, int node_num_max)
{
    int i;
    for(i=0;i<node_num_max;i++)
    {
   		if((node_array+i)->bound_index_temp==-1)
    	{
    		(node_array+i)->node_temp=ini_temp;
     	}
    } 
}

void initialization_flow(double* prev_solu, node* node_array, int* node_not_pbound_index, double ini_num1,double ini_num2,int node_not_bound,int elem_num_max)
{
    int i;
    for(i=0; i<node_not_bound; i++)
    {
        *(prev_solu+i)=ini_num1;
        (node_array+*(node_not_pbound_index+i))->node_pressure=ini_num1;
        
    }
    
    for(i=node_not_bound;i<node_not_bound+elem_num_max;i++)
    {
    	*(prev_solu+i)=ini_num2;
    }
}

void cal_constants(double* constants_array,int* node_not_pbound_index,elem* elem_array,node* node_array,int node_not_pbound,int elem_num_max)
{
    double inlet_rad,outlet_rad;
    double grav_force,rot_force;
    int i;
	
    for(i=0; i<node_not_pbound; i++)
    {
        if((node_array+*(node_not_pbound_index+i))->flow_rate_bound_flag==1)
        {
            *(constants_array+i)=(node_array+*(node_not_pbound_index+i))->node_flow_rate;
        }
        else
        {
            *(constants_array+i)=0;
        }
    }

    for(i=0; i<elem_num_max; i++)
    {
        inlet_rad=(elem_array+i)->inlet_radius;
        outlet_rad=(elem_array+i)->outlet_radius;
        grav_force=(elem_array+i)->ave_area*(elem_array+i)->ave_dens*GRAV_FACTOR*(elem_array+i)->length*(-sin((elem_array+i)->flow_angle*PI/180));
        rot_force=0.5*(elem_array+i)->ave_dens*pow((elem_array+i)->krot,2)*pow((elem_array+i)->omega,2)*(elem_array+i)->ave_area*(pow(outlet_rad,2)-pow(inlet_rad,2));
//        *(constants_array+node_not_pbound+i)=grav_force+0.5*pow((elem_array+i)->omega,2)*(elem_array+i)->ave_dens*(pow(outlet_rad,2)-pow(inlet_rad,2))*(elem_array+i)->ave_area;
        *(constants_array+node_not_pbound+i)=grav_force+rot_force;
//        printf("%\t\t%lf,%lf,\n",grav_force,rot_force);
//		printf("\n\t\t\t%lf",0.5*pow((elem_array+i)->omega,2)*(elem_array+i)->ave_dens*(pow(outlet_rad,2)-pow(inlet_rad,2))*(elem_array+i)->ave_area);
	}
}

void nonlinear_compute(double* compute_solu,double** nonlinear_matrix,double* prev_solu,double* constants_array,node* node_array,int node_num_max,int elem_num_max,int node_not_pbound)
{
    int i,j;
    int para_index;
    int para_num=node_not_pbound+elem_num_max;
    double temp;

    for(i=0; i<para_num; i++)
    {
        *(compute_solu+i)=0;
        para_index=0;
        for(j=0; j<node_num_max; j++)
        {
            if((node_array+j)->bound_index_pres!=-1)
            {
                temp=(*(*(nonlinear_matrix+i)+j*2+1));
                *(compute_solu+i)-=temp;
            }
            else
            {
                temp=(*(*(nonlinear_matrix+i)+j*2))*(*(prev_solu+para_index));
                para_index++;
//				printf("%d\t%d\t%lf\n",j,bound_or_not(node_not_bound_index,j,node_not_bound),temp);
                *(compute_solu+i)-=temp;
            }
        }
        for(j=node_num_max; j<node_num_max+elem_num_max; j++)
        {
            temp=*(*(nonlinear_matrix+i)+node_num_max*2+(j-node_num_max)*3)*pow(*(prev_solu+node_not_pbound+j-node_num_max),2);
            temp+=*(*(nonlinear_matrix+i)+node_num_max*2+(j-node_num_max)*3+1)*(*(prev_solu+node_not_pbound+j-node_num_max));
            temp+=*(*(nonlinear_matrix+i)+node_num_max*2+(j-node_num_max)*3+2);
            *(compute_solu+i)-=temp;
        }
        *(compute_solu+i)-=*(constants_array+i);
//        printf("\n\t\t\t%lf",*(compute_solu+i));
    }
}

void swap_line(double** left_matrix,double* right_array,int i,int j,int elem_num_max,int node_not_bound)
{
    double temp;
    int para_num=node_not_bound+elem_num_max;
    int k;
    temp=*(right_array+i);
    *(right_array+i)=*(right_array+j);
    *(right_array+j)=temp;
    for(k=0; k<para_num; k++)
    {
        temp=*(*(left_matrix+i)+k);
        *(*(left_matrix+i)+k)=*(*(left_matrix+j)+k);
        *(*(left_matrix+j)+k)=temp;
    }
}

int Gaussian_eliminate(double* solution,double** left_matrix,double* right_hand_array,int elem_num_max,int node_not_bound)
{
    int para_num=node_not_bound+elem_num_max;
    int i,j,k,kk;
    int max_index;
    double max_temp;
    double coeff_temp;
    double sigma_temp;
    int la,lb;
/*    FILE* Gus_test=fopen("..\\test_files\\g_test.txt","w");
    FILE* Gus_test2=fopen("..\\test_files\\g_test2.txt","w");
    
    printf("%d\n\n",para_num);
    for(k=0;k<para_num;k++)
    {
    	for(j=0;j<para_num;j++)
    	{
	    	fprintf(Gus_test2,"%lf\t",*(*(left_matrix+k)+j));
	    }
	    fprintf(Gus_test2,"%.6lf\n",*(right_hand_array+k));
    }*/

    for(i=0; i<para_num-1; i++)
    {
        max_temp=fabs(*(*(left_matrix+i)+i));
        max_index=i;
        for(j=i+1; j<para_num; j++)
        {
            if(fabs(*(*(left_matrix+j)+i))>max_temp)
            {
                max_temp=fabs(*(*(left_matrix+j)+i));
                max_index=j;
            }
        }

        if(fabs(max_temp-0)<=ERR)
        {
            return -1;
        }
        if(max_index!=i)
        {
            swap_line(left_matrix,right_hand_array,i,max_index,elem_num_max,node_not_bound);
        }

        for(j=i+1; j<para_num; j++)
        {
            coeff_temp=*(*(left_matrix+j)+i)/(*(*(left_matrix+i)+i));
            *(right_hand_array+j)=*(right_hand_array+j)-coeff_temp*(*(right_hand_array+i));
            for(k=i; k<para_num; k++)
            {
                *(*(left_matrix+j)+k)=*(*(left_matrix+j)+k)-coeff_temp*(*(*(left_matrix+i)+k));
            }
        }
    }
/*   for(k=0;k<para_num;k++)
    {
    	for(j=0;j<para_num;j++)
    	{
	    	fprintf(Gus_test,"%lf,",*(*(left_matrix+k)+j));
	    }
	    fprintf(Gus_test,"\n");
    }*/
    
    if(fabs(*(*(left_matrix+para_num-1)+para_num-1))<ERR)
    {
    	return -1;
    }
    *(solution+para_num-1)=*(right_hand_array+para_num-1)/(*(*(left_matrix+para_num-1)+para_num-1));
//    printf("%lf,%lf,\n",*(right_hand_array+para_num-1),*(*(left_matrix+para_num-1)+para_num-1));
    for(k=para_num-2; k>-1; k--)
    {
        sigma_temp=0;
        for(kk=k+1; kk<para_num; kk++)
        {
            sigma_temp+=*(*(left_matrix+k)+kk)*(*(solution+kk));
        }
        *(solution+k)=(*(right_hand_array+k)-sigma_temp)/(*(*(left_matrix+k)+k));
    }
    
//    printf("\n%d\n",node_not_bound);
//    if(node_not_bound==16)
//    {
//    	for(k=0;k<node_not_bound+elem_num_max;k++)
//   		{
//   			printf("%lf\n",*(solution+k));
//    	}
//    }
//    
    
//    fclose(Gus_test);
//    fclose(Gus_test2);
    return 0;
}

void correction(double* solu,double* correc,int elem_num_max,int node_not_pbound)
{
    int i;
    for(i=0; i<node_not_pbound+elem_num_max; i++)
    {
        *(solu+i)+=*(correc+i);
//       printf("%lf,",*(solu+i));
    }
}

void array_subtract(double* sub,double* prev_temp_solu,double* this_temp_solu,int para_num)
{
    int i;
    for(i=0; i<para_num; i++)
    {
//    	printf("%lf,,,%lf,,,\n",*(prev_temp_solu+i),*(this_temp_solu+i));
        *(sub+i)=*(prev_temp_solu+i)-*(this_temp_solu+i);
    }
}

void array_copy(double* prev_solu,double* this_solu,int para_num)
{
    int i;
    for(i=0; i<para_num; i++)
    {
        *(prev_solu+i)=*(this_solu+i);
    }
}

int converge_or_not(double* correc,double* p_residual,int elem_num_max,int node_not_bound,FILE* output_console,setup* p_setups,int equation_flag)
{
    int i;
    int flag;
    double temp;
    temp=fabs(*(correc));
    for(i=0; i<node_not_bound+elem_num_max; i++)
    {
        if(fabs(*(correc+i))>temp)
        {
            temp=fabs(*(correc+i));
        }
    }
    if(equation_flag==0)
    {
        printf("  MAX FLOW RESIDUAL VALUE=%.7lf\n",temp);
        if(output_console!=NULL)
        {
            fprintf(output_console,"  MAX FLOW RESIDUAL VALUE=%.7lf\n",temp);
        }
        *p_residual=temp;
        if(temp<p_setups->flow_converge_flag)
        {
            return 0;
        }
    }
    else if(equation_flag==1)
    {
        printf("  MAX ENERGY RESIDUAL VALUE=%.7lf\n",temp);
        if(output_console!=NULL)
        {
            fprintf(output_console,"  MAX ENERGY RESIDUAL VALUE=%.7lf\n",temp);
        }
        *p_residual=temp;
        if(temp<p_setups->energy_converge_flag)
        {
            return 0;
        }
    }
    return 1;
}

int solve_and_output(FILE* elem_file,FILE* bound_file,FILE* setup_file,FILE* output_elem_file,FILE* output_node_file,FILE* residual_file,FILE* output_console)
{
    elem* elem_array=NULL;
    int* index_table=NULL;
    node* node_array=NULL;
    bound* bound_array=NULL;
    setup setups;
    setup* p_setups=&setups;

    int memory_alloc_flag=0;
    int memory_alloc_index=0;

    double** p_energy_eqns=NULL;
    double* right_hand_array=NULL;

    double** nonlinear_matrix=NULL;
    double** jacobi_matrix=NULL;

    int node_not_pbound=0;
    int node_not_tbound=-1;

    int* node_not_pbound_index=NULL;
    int* node_not_tbound_index=NULL;
    int i,j,k;

    int node_num_max;
    int elem_num_max;
    int bound_num_max;
    int fluid_type;
    double density_in_file;
    int energy_flag;
    double ref_temp;
    int friction_flag;
    int fric_coeff_flag;
    double ini_pres;
    double ini_flow_rate;
    double ini_temp_num;

    int flow_converge_flag=-1;
    int temp_converge_flag=-1;

    int total_iteration_flag=0;
    int flow_iteration_flag;

    double flow_residual=0;
    double temp_residual=0;

    double* prev_flow_solu=NULL;
    double* prev_enth_solu=NULL;
    double* this_enth_solu=NULL;
    double* subtract=NULL;
    double* constants_array=NULL;
    double* correc=NULL;
    double* prev_compute=NULL;
    double* prev_temp=NULL;
    double* this_temp=NULL;
    int update_elem_flag=0;
    int update_node_flag=0;
    int failure_flag=0;
    int iii=0;

    if(read_setup(setup_file,p_setups,output_console)==-1)
    {
        return -1;
    }

    node_num_max=p_setups->snode_num_max;
    elem_num_max=p_setups->selem_num_max;
    bound_num_max=p_setups->sbound_num_max;
    fluid_type=p_setups->sfluid_type;
    energy_flag=p_setups->senergy_flag;
    ref_temp=p_setups->sref_temp;
    friction_flag=p_setups->sfriction_flag;
    fric_coeff_flag=p_setups->sfric_coeff_flag;
    ini_pres=p_setups->sini_pres;
    ini_flow_rate=p_setups->sini_flow_rate;
    ini_temp_num=p_setups->sini_temp;
//	printf("%d %d\n",compress_flag,energy_flag);

    node_array=(node*)malloc(node_num_max*sizeof(node));
    elem_array=(elem*)malloc(elem_num_max*sizeof(elem));
    bound_array=(bound*)malloc(bound_num_max*sizeof(bound));
    index_table=(int*)malloc(node_num_max*node_num_max*sizeof(int));

    if(create_elem_array(elem_array,elem_file,p_setups)==-1)
    {
        free(node_array);
        node_array=NULL;
        free(elem_array);
        elem_array=NULL;
        free(bound_array);
        bound_array=NULL;
        free(index_table);
        index_table=NULL;
        return -2;
    }
    create_pipe_index(index_table,node_num_max,elem_array,elem_num_max);
    if(create_bound_array(bound_array,bound_file,bound_num_max)==-1)
    {
        free(node_array);
        node_array=NULL;
        free(elem_array);
        elem_array=NULL;
        free(bound_array);
        bound_array=NULL;
        free(index_table);
        index_table=NULL;
        return -3;
    }
    create_node_array(node_array,elem_array,bound_array,index_table,&node_not_pbound,p_setups);

    node_not_pbound_index=(int*)malloc(node_not_pbound*sizeof(int));
    create_node_not_pbound_index(node_not_pbound_index,node_array,node_num_max,node_not_pbound);

    nonlinear_matrix=(double**)malloc((node_not_pbound+elem_num_max)*sizeof(double*));
    jacobi_matrix=(double**)malloc((node_not_pbound+elem_num_max)*sizeof(double*));
    prev_flow_solu=(double*)malloc((node_not_pbound+elem_num_max)*sizeof(double));
    prev_compute=(double*)malloc((node_not_pbound+elem_num_max)*sizeof(double));
    constants_array=(double*)malloc((node_not_pbound+elem_num_max)*sizeof(double));
    correc=(double*)malloc((node_not_pbound+elem_num_max)*sizeof(double));
    
    prev_temp=(double*)malloc((node_num_max)*sizeof(double));
    this_temp=(double*)malloc((node_num_max)*sizeof(double));
    subtract=(double*)malloc((node_num_max)*sizeof(double));
    
    if(node_array==NULL||elem_array==NULL||bound_array==NULL||index_table==NULL||node_not_pbound_index==NULL||nonlinear_matrix==NULL||jacobi_matrix==NULL||prev_flow_solu==NULL||prev_compute==NULL||constants_array==NULL||correc==NULL||prev_temp==NULL||this_temp==NULL||subtract==NULL)
    {
        memory_alloc_flag=-1;
    }

    if(memory_alloc_flag!=-1)
    {
        for(i=0; i<node_not_pbound+elem_num_max; i++)
        {
            *(nonlinear_matrix+i)=(double*)malloc((node_num_max*2+elem_num_max*3)*sizeof(double));
            *(jacobi_matrix+i)=(double*)malloc((node_not_pbound+elem_num_max)*sizeof(double));
            if(*(nonlinear_matrix+i)==NULL||*(jacobi_matrix+i)==NULL)
            {
                memory_alloc_flag=-1;
                memory_alloc_index=i;
                break;
            }
        }
        if(memory_alloc_flag==-1)
        {
            for(i=0; i<memory_alloc_index; i++)
            {
                free(*(nonlinear_matrix+i));
                *(nonlinear_matrix+i)=NULL;
                free(*(jacobi_matrix+i));
                *(jacobi_matrix+i)=NULL;
            }
        }
    }

    if(memory_alloc_flag==-1)
    {
        failure_flag=-4;
        goto free_memory;
    }

    initialization_flow(prev_flow_solu,node_array,node_not_pbound_index,ini_pres,ini_flow_rate,node_not_pbound,elem_num_max);
    if(energy_flag==1)
    {
         initialization_node_temp(node_array,ini_temp_num,node_num_max);
		 calc_node_enth(node_array,node_num_max,fluid_type);   
    } 
    
    printf("\n# STARTING SOLVING...\n");
    if(output_console!=NULL)
    {
        fprintf(output_console,"\n# STARTING SOLVING...\n");
    }

    do
    {
        if(energy_flag==1)
        {
            printf("\n# TOTAL ITERATION %d: \n",total_iteration_flag+1);
            if(output_console!=NULL)
            {
                fprintf(output_console,"\n# TOTAL ITERATION %d: \n",total_iteration_flag+1);
            }
            fprintf(residual_file,"%d\t\t",total_iteration_flag);
            read_node_temp(prev_temp,node_array,node_num_max);
			update_elem_temp(elem_array,node_array,elem_num_max);
        }
		
        flow_iteration_flag=0;

        do
        {
            printf("\tFLOW ITERATION %d:",flow_iteration_flag+1);
            if(output_console!=NULL)
            {
                fprintf(output_console,"\tFLOW ITERATION %d:",flow_iteration_flag+1);
            }
            update_node_flag=update_node_vars(node_array,prev_flow_solu,node_not_pbound_index,node_not_pbound,p_setups);
            if(update_node_flag==-1)
            {
            	printf("!WARNING: STEAM/WATER DENSITY CALCULATION ERROR (NODE). RESULTS MIGHT BE WRONG.\n");
            }
            update_elem_flag=update_elem_vars(elem_array,prev_flow_solu,node_array,node_not_pbound,p_setups);
            if(update_elem_flag==-1)
            {
            	printf("!WARNING: STEAM/WATER DENSITY CALCULATION ERROR (ELEM). RESULTS MIGHT BE WRONG.\n");
            }
            if(update_elem_flag==-2)
            {
                failure_flag=-5;
                goto free_memory;
            }
            
            create_nonlinear_matrix(nonlinear_matrix,elem_array,node_array,bound_array,index_table,elem_num_max,node_num_max,node_not_pbound);
            create_jacobi(jacobi_matrix,nonlinear_matrix,prev_flow_solu,node_not_pbound_index,elem_num_max,node_num_max,node_not_pbound);
            cal_constants(constants_array,node_not_pbound_index,elem_array,node_array,node_not_pbound,elem_num_max);
            nonlinear_compute(prev_compute,nonlinear_matrix,prev_flow_solu,constants_array,node_array,node_num_max,elem_num_max,node_not_pbound);
            if(Gaussian_eliminate(correc,jacobi_matrix,prev_compute,elem_num_max,node_not_pbound)==-1)
            {
                failure_flag=-6;
                goto free_memory;
            }
            correction(prev_flow_solu,correc,elem_num_max,node_not_pbound);
            flow_converge_flag=converge_or_not(correc,&flow_residual,elem_num_max,node_not_pbound,output_console,p_setups,0);
            flow_iteration_flag++;
            if(energy_flag==1)
            {
                fprintf(residual_file,"%.10lf\t",flow_residual);
            }
            else
            {
                fprintf(residual_file,"%d\t%.10lf\n",flow_iteration_flag,flow_residual);
            }
        }
        while(flow_converge_flag!=0&&flow_iteration_flag<p_setups->max_flow_iter_times);
        if(flow_iteration_flag<p_setups->max_flow_iter_times)
        {
            printf("\n\t! FLOW SOLUTION CONVERGED AT FLOW ITERATION STEP %d.\n\n",flow_iteration_flag);
            if(output_console!=NULL)
            {
                fprintf(output_console,"\n\t! FLOW SOLUTION CONVERGED AT FLOW ITERATION STEP %d.\n\n",flow_iteration_flag);
            }
            if(energy_flag==0)
            {
                break;
            }
        }
        else
        {
            printf("\n\t!! FLOW SOLUTION UNCONVERGED AT TOTAL ITERATION STEP %d.\n",total_iteration_flag+1);
            printf("!! PROGRAM ABORTED.\n");
            if(output_console!=NULL)
            {
                fprintf(output_console,"\n\t!! FLOW SOLUTION UNCONVERGED AT TOTAL ITERATION STEP %d.\n",total_iteration_flag+1);
                fprintf(output_console,"!! PROGRAM ABORTED.\n");
            }
            failure_flag=-7;
            goto free_memory;
        }
        
        node_not_tbound=calc_temp_bound(bound_array,index_table,node_array,elem_array,p_setups);
        node_not_tbound_index=(int*)malloc(node_not_tbound*sizeof(int));
        create_node_not_tbound_index(node_not_tbound_index,node_array,node_num_max,node_not_tbound);
        p_energy_eqns=(double**)malloc(node_not_tbound*sizeof(double*));
        right_hand_array=(double*)malloc(node_not_tbound*sizeof(double));
        
        prev_enth_solu=(double*)malloc(node_not_tbound*sizeof(double));
        this_enth_solu=(double*)malloc(node_not_tbound*sizeof(double));

        if(node_not_tbound_index==NULL||p_energy_eqns==NULL||right_hand_array==NULL||prev_enth_solu==NULL||this_enth_solu==NULL||subtract==NULL)
        {
            memory_alloc_flag=-1;
        }
        if(memory_alloc_flag!=-1)
        {
            for(i=0; i<node_not_tbound; i++)
            {
                *(p_energy_eqns+i)=(double*)malloc(node_not_tbound*sizeof(double));
                if(*(p_energy_eqns+i)==NULL)
                {
                    memory_alloc_flag=-1;
                    memory_alloc_index=i;
                    break;
                }
            }
            if(memory_alloc_flag==-1)
            {
                for(i=0; i<memory_alloc_index; i++)
                {
                    free(*(p_energy_eqns+i));
                    *(p_energy_eqns+i)=NULL;
                }
            }
        }  
		if(memory_alloc_flag==-1)
		{
			failure_flag=-4;
			goto free_memory;
		}
//		 printf("FLAG.\n");     
        create_energy_equations_new(p_energy_eqns,right_hand_array,elem_array,node_array,index_table,node_not_tbound,node_not_tbound_index,node_num_max);
//        printf("FLAG..\n");
		if(Gaussian_eliminate(this_enth_solu,p_energy_eqns,right_hand_array,0,node_not_tbound)==-1)
		{
			failure_flag=-8;
			goto free_memory;
		}
//		printf("FLAG...\n");
        update_node_enth(node_array,this_enth_solu,node_not_tbound_index,node_not_tbound);
        if(update_node_temp_new(node_array,node_not_tbound_index,node_not_tbound,fluid_type)==-1)
        {
        	failure_flag=-9;
        	goto free_memory;
        }
        read_node_temp(this_temp,node_array,node_num_max);
        array_subtract(subtract,prev_temp,this_temp,node_num_max);
//        array_subtract(subtract,prev_temp_solu,this_temp_solu,node_not_tbound+elem_num_max);
        temp_converge_flag=converge_or_not(subtract,&temp_residual,0,node_num_max,output_console,p_setups,1);
//        array_copy(prev_temp_solu,this_temp_solu,node_not_tbound+elem_num_max);
 //       update_node_temp(node_array,prev_temp_solu,node_not_tbound_index,node_not_tbound);
        fprintf(residual_file,"\t+%.10lf\n",temp_residual);
		
		for(iii=0;iii<node_not_tbound;iii++)
		{
			free(*(p_energy_eqns+iii));
			*(p_energy_eqns+iii)=NULL;
		}
		if(p_energy_eqns!=NULL)
		{
			free(p_energy_eqns);
			p_energy_eqns=NULL;
		}
		
//		printf("FLAG....\n");
		free(node_not_tbound_index);
		node_not_tbound_index=NULL;
//		printf("FLAG.....\n");
		free(right_hand_array);
		right_hand_array=NULL;
//		printf("FLAG......\n");
		free(prev_enth_solu);
		prev_enth_solu=NULL;
//		printf("FLAG..........\n");
		free(this_enth_solu);
		this_enth_solu=NULL;      
//        printf("FLAG...............\n");
        total_iteration_flag++;
    }
    while(temp_converge_flag!=0&&total_iteration_flag<p_setups->max_total_iter_times);
    if(energy_flag==1&&total_iteration_flag	<p_setups->max_total_iter_times)
    {
        printf("\n! ENERGY SOLUTION CONVERGED AT TOTAL ITERATION STEP %d.\n\n",total_iteration_flag);
        if(output_console!=NULL)
        {
            fprintf(output_console,"\n! ENERGY SOLUTION CONVERGED AT TOTAL ITERATION STEP %d.\n\n",total_iteration_flag);
        }
    }
    else if(energy_flag==1&&total_iteration_flag==p_setups->max_total_iter_times)
    {
        printf("\n!! ENERGY SOLUTION UNCONVERGED AT TOTAL ITERATION STEP %d.\n\n",p_setups->max_total_iter_times);
        printf("!! PROGRAM ABORTED.\n");
        if(output_console!=NULL)
        {
            fprintf(output_console,"\n!! ENERGY SOLUTION UNCONVERGED AT TOTAL ITERATION STEP %d.\n\n",p_setups->max_total_iter_times);
            fprintf(output_console,"!! PROGRAM ABORTED.\n");
        }
        failure_flag=-10;
        goto free_memory;
    }
    printf("# SOLUTION:\n");
    printf("# NODE\tPRESSURE\tND_TEMPERATURE\tND_DENSITY\n");
    fprintf(output_node_file,"# THIS IS AN NODE INFORMATION OUTPUT FILE FOR FLOW CALCULATION PROGRAM (VERSION 0.1).\n");
    fprintf(output_node_file,"# VARIABLES: NODE_NUM, PRESSURE, TEMPERATURE, DENSITY, ROT_RADIUS, ROT_VELOCITY.\n");
    for(i=0; i<node_num_max; i++)
    {
        printf("  %d\t%.6lf\t%.6lf\t%.6lf\n",i,(node_array+i)->node_pressure,(node_array+i)->node_temp,(node_array+i)->node_density);
        fprintf(output_node_file,"  %d\t%.6lf\t%.6lf\t%.6lf\t",i,(node_array+i)->node_pressure,(node_array+i)->node_temp,(node_array+i)->node_density);
        fprintf(output_node_file,"%.6lf\t%.6lf\n",(node_array+i)->node_radius,(node_array+i)->node_omega);
    }
    printf("\n");
    printf("# PIPE\tFLOW_RATE\tAVE_TEMPERATURE\tRE_NUM\tHTC\n");
    fprintf(output_elem_file,"# THIS IS AN ELEMEMT INFORMATION OUTPUT FILE FOR FLOW CALCULATION PROGRAM (VERSION 0.1).\n");
    fprintf(output_elem_file,"# VARIABLES: ELEM_NUM, FLOW_RATE, AVERAGE_TEMPERATURE, AVERAGE_DENS, RE_NUM, HEAT_TRANSFER_COEFFICIENT.\n");
    for(i=0; i<elem_num_max; i++)
    {
        printf("  %d\t%.6lf\t%.6lf\t%.6lf\t%.6lf\n",i,(elem_array+i)->flow_rate,(elem_array+i)->ave_temp,(elem_array+i)->re_num,(elem_array+i)->ht_coeff);
        fprintf(output_elem_file,"  %d\t%.6lf\t%.6lf\t%.6lf\t%.6lf\t%.6lf\t",i,(elem_array+i)->flow_rate,(elem_array+i)->ave_temp,(elem_array+i)->ave_dens,(elem_array+i)->re_num,(elem_array+i)->ht_coeff);
        fprintf(output_elem_file,"%.6lf\t%.6lf\n",(node_array+(elem_array+i)->inlet_node_num)->node_pressure-(node_array+(elem_array+i)->outlet_node_num)->node_pressure,(node_array+(elem_array+i)->outlet_node_num)->node_temp-(node_array+(elem_array+i)->inlet_node_num)->node_temp);
    }

free_memory:
    free(node_array);
    node_array=NULL;
    free(elem_array);
    elem_array=NULL;
    free(bound_array);
    bound_array=NULL;
    free(index_table);
    index_table=NULL;
    free(p_setups);
    p_setups=NULL;
    free(node_not_pbound_index);
    node_not_pbound_index=NULL;
    free(prev_flow_solu);
    prev_flow_solu=NULL;
    free(prev_temp);
    prev_temp=NULL;
    free(this_temp);
    this_temp=NULL;
    free(subtract);
    subtract=NULL;
    
    if(failure_flag!=-4)
    {
        for(i=0; i<node_not_pbound+elem_num_max; i++)
        {
            free(*(nonlinear_matrix+i));
            *(nonlinear_matrix+i)=NULL;
            free(*(jacobi_matrix+i));
            *(jacobi_matrix+i)=NULL;
        }
    }
    free(nonlinear_matrix);
    nonlinear_matrix=NULL;
    free(jacobi_matrix);
    jacobi_matrix=NULL;
    free(prev_compute);
    prev_flow_solu=NULL;
    free(constants_array);
    constants_array=NULL;
    free(correc);
    correc=NULL;
    if(energy_flag==1&&failure_flag!=-8)
    {
    	if(node_not_tbound_index!=NULL)
    	{
	    	free(node_not_tbound_index);
        	node_not_tbound_index=NULL;
	    }
        
        if(prev_enth_solu!=NULL)
    	{
	    	free(prev_enth_solu);
        	prev_enth_solu=NULL;
	    }
        
        if(this_enth_solu!=NULL)
        {
        	 free(this_enth_solu);
        	this_enth_solu=NULL;
        }
       	if(right_hand_array!=NULL)
       	{
       		free(right_hand_array);
        	right_hand_array=NULL;
       	}
        if(p_energy_eqns!=NULL)
        {
        	if(failure_flag!=-3)
        	{
            	for(i=0; i<node_not_tbound+elem_num_max; i++)
            	{
            		if(*(p_energy_eqns+i)!=NULL)
            		{
	            		free(*(p_energy_eqns+i));
                		*(p_energy_eqns+i)=NULL;
	            	}  
            	}
        	}
			free(p_energy_eqns);
        	p_energy_eqns=NULL;
  		}        
    }
    return failure_flag;
}

int main(void)
{
    FILE* elem_file=NULL;
    FILE* bound_file=NULL;
    FILE* setup_file=NULL;
    FILE* output_elem_file=NULL;
    FILE* output_node_file=NULL;
    FILE* residual_file=NULL;
    FILE* output_console=NULL;
    int failure_flag;
    char choice;

    output_console=fopen("..\\output_console.txt","a");
    if(output_console==NULL)
    {
        printf("\nWARNING: FILE I/O ERROR.\nCANNOT CREATE OR OPEN CONSOLE OUTPUT FILE.\n");
        printf("SOLVING PROCESS WILL NOT BE RECORDED.\nA FILE CHECK IS STRONGLY RECOMMENDED.\n");
        printf("PRESS 'Y' TO ABORT PROGRAM AND CHECK FILES, OTHERS TO START SOLVING:");
        fflush(stdin);
        choice=getchar();
        if(choice=='Y'&&choice=='y')
        {
            goto close_files;
        }
    }
    welcome(output_console);
    elem_file=fopen("..\\input_files\\input_elems.dat","r");
    bound_file=fopen("..\\input_files\\input_boundaries.dat","r");
    setup_file=fopen("..\\input_files\\input_setups.dat","r");

    if(elem_file==NULL||bound_file==NULL||setup_file==NULL)
    {
        printf("\nFATAL ERROR: FILE I/O ERROR.\nCHECK INPUT FILES.\nPROGRAM ABORTED.\n");
        fprintf(output_console,"\nFATAL ERROR: FILE I/O ERROR.\nCHECK INPUT FILES.\nPROGRAM ABORTED.\n");
        goto close_files;
    }

    output_elem_file=fopen("..\\output_files\\output_elems.dat","w");
    output_node_file=fopen("..\\output_files\\output_nodes.dat","w");
    residual_file=fopen("..\\output_files\\output_residuals.dat","w");

    if(output_elem_file==NULL||output_node_file==NULL||residual_file==NULL)
    {
        printf("\nFATAL ERROR: FILE I/O ERROR.\nCHECK OUTPUT FILES.\nPROGRAM ABORTED.\n");
        fprintf(output_console,"\nFATAL ERROR: FILE I/O ERROR.\nCHECK OUTPUT FILES.\nPROGRAM ABORTED.\n");
        goto close_files;
    }
    failure_flag=solve_and_output(elem_file,bound_file,setup_file,output_elem_file,output_node_file,residual_file,output_console);
    if(failure_flag==-1)
    {
        printf("\nFATAL ERROR: CALCULATION SETUP ERROR.\nPLEASE CHECK SETUP FILE.\nPROGRAM ABORTED.\n");
        fprintf(output_console,"\nFATAL ERROR: CALCULATION SETUP ERROR.\nPLEASE CHECK SETUP FILE.\nPROGRAM ABORTED.\n");
        goto close_files;
    }

    if(failure_flag==-2)
    {
        printf("FATAL ERROR: ELEMENT INFORMATION ERROR.\nPLEASE CHECK THE FILE.\nPROGRAM ABORTED.\n");
        fprintf(output_console,"FATAL ERROR: ELEMENT INFORMATION ERROR.\nPLEASE CHECK THE FILE.\nPROGRAM ABORTED.\n");
        goto close_files;
    }
    if(failure_flag==-3)
    {
        printf("\nFATAL ERROR: BOUNDARY CONDITION ERROR.\nPLEASE CHECK THE FILE.\nPROGRAM ABORTED.\n");
        fprintf(output_console,"\nFATAL ERROR: BOUNDARY CONDITION ERROR.\nPLEASE CHECK THE FILE.\nPROGRAM ABORTED.\n");
        goto close_files;
    }
    if(failure_flag==-4)
    {
        printf("\nFATAL ERROR: MEMORY ALLOCATION FAILURE.\nPROGRAM ABORTED.\n");
        fprintf(output_console,"\nFATAL ERROR: MEMORY ALLOCATION FAILURE.\nPROGRAM ABORTED.\n");
        goto close_files;
    }
    else if(failure_flag==-5)
    {
        printf("\nFATAL ERROR: RESIST COEFFICIENT CALCULATION UNCONVERGED.\nPROGRAM ABORTED.\n");
        fprintf(output_console,"\nFATAL ERROR: RESIST COEFFICIENT CALCULATION UNCONVERGED.\nPROGRAM ABORTED.\n");
        goto close_files;
    }
    else if (failure_flag==-6)
    {
    	printf("\nFATAL ERROR: GAUSSIAN ELIMINATION ERROR (FLOW CALCULATION).\nPROGRAM ABORTED.\n");
    	fprintf(output_console,"\nFATAL ERROR: GAUSSIAN ELIMINATION ERROR.\nPROGRAM ABORTED.\n");
    	goto close_files;
    }
    else if(failure_flag==-8)
    {
        printf("\nFATAL ERROR: GAUSSIAN ELIMINATION ERROR (ENERGY CALCULATION).\nPROGRAM ABORTED.\n");
        fprintf(output_console,"\nFATAL ERROR: GAUSSIAN ELIMINATION ERROR (ENERGY CALCULATION).\nPROGRAM ABORTED.\n");
        goto close_files;
    }
    else if(failure_flag==-10)
    {
    	printf("\nFATAL ERROR: NODE TEMPERATURE CALCULATION ERROR.\nPROGRAM ABORTED.\n");
        fprintf(output_console,"\nFATAL ERROR: NODE TEMPERATURE CALCULATION ERROR.\nPROGRAM ABORTED.\n");
        goto close_files;
    }
    
close_files:
    if(elem_file!=NULL)
    {
        fclose(elem_file);
        elem_file=NULL;
    }
    if(bound_file!=NULL)
    {
        fclose(bound_file);
        bound_file=NULL;
    }
    if(setup_file!=NULL)
    {
        fclose(setup_file);
        setup_file=NULL;
    }
    if(output_elem_file!=NULL)
    {
        fclose(output_elem_file);
        output_elem_file=NULL;
    }
    if(output_node_file!=NULL)
    {
        fclose(output_node_file);
        output_node_file=NULL;
    }
    if(residual_file!=NULL)
    {
        fclose(residual_file);
        residual_file=NULL;
    }
    if(output_console!=NULL)
    {
        fclose(output_console);
        output_console=NULL;
    }
    printf("\nPLEASE PRESS ENTER TO EXIT:");
    fflush(stdin);
    getchar();

    return failure_flag;
}


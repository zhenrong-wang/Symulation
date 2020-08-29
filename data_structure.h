#ifndef DATA_STRUCTURE_H_
#define DATA_STRUCTURE_H_

typedef struct{
	int node_num;
	double node_temp;
	double node_pressure;
	double node_density;
	double node_radius;
	double node_omega;
	double node_flow_rate;
	double node_enth;
	int flow_rate_bound_flag;
	int flow_bound_flag;
	int bound_index_pres;
	int bound_index_temp;
	int table_index_pres;
	int table_index_temp;
}node;

//1-round channel, 2-rect channel, 3-trapzium channel, 4-triangle channel, 5-ecllipse channel, 6-film cooling hole
//11-round rib channel, 12-rect rib channel
//21-parallel wall
//31-uturn sharp, 32-uturn round
//41-impingement stagger point, 42-impingement head area, 43-impingement hody.
//51-pin-fin
typedef struct{
	int elem_type;
	int elem_num;
	int inlet_node_num;
	int outlet_node_num;
	
	double length;
	double diameter; //for 1&31
	double section_length; //for 2&21&22
	double section_width; //for 2&21&22
	double wall_distance; //for 11 only
	double pitch; //for 41 only
//	double wall_thick; //for 41 only
	double target_dia; //for 41 only
	double target_dist; //for 41 only
	double cor_coeff;  //for 2~5
	double pinfin_row; //for 51
	double incline_angle; //for 6 only
	double param_1;
	double param_2;
	
	double roughness_factor;
	double rib_height;
	double rib_space;
	
	double resist_factor;//calculate_by Renolds_number
	
	double flow_angle;
	double flow_rate; //unknown parameter
	double ave_dens;//unknown parameter
	double ave_temp;//calculate in the program
	double re_num;
	
	double inlet_area;
	double inlet_radius; //inlet rotating radius
	double outlet_area; 
	double outlet_radius; //outlet rotating radius
	double omega; //rotating speed
	double krot; //rotation factor, normally equals to 1
	
	double ave_area; //calculate in the program
	double hydr_dia; //calculate in the program
	
	double ht_coeff;//!calculate in the program
	double ht_area;	//heat transfer area
	double elem_specific_heat;//!calculate in the program
	double ht_temp;	
	
	double fluid_cp;
	double fluid_visc;
	double fluid_thcond;
	double fluid_pranum;
		
}elem;

//boundry condition
typedef struct{
	int node_num; //boundry node
	int flow_bound_type; //boundry condition type for flow (mass flow rate or pressure)
	double flow_bound_value; //boundary condition value
	double bound_temp; //boundary temperature
}bound;

//calculation setups
typedef struct{
	int snode_num_max;
	int selem_num_max;
	int sbound_num_max;
	
	int sfluid_type;
	
	int senergy_flag;
	int shtc_flag;
	int scorrelation_flag;
	double sref_temp;
	
	int sfriction_flag;
	int sfric_coeff_flag;
	
	double sini_pres;
	double sini_flow_rate;
	double sini_temp;
	
	int max_fric_iter_times;
	int max_flow_iter_times;
	int max_total_iter_times;
	
	double flow_converge_flag;
	double energy_converge_flag; 
}setup;


#define GRAV_FACTOR 9.8
#define PI 3.14159265359
#define WATER_SPECIFIC_HEAT 4200
#define ERR 1e-6
#define WATER_DENS_23 997.77

#endif

/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "custom.h"
#include "../BioFVM/BioFVM.h"  
using namespace BioFVM;
// using namespace std;

// declare cell definitions here 

std::vector<bool> nodes;

void create_cell_types( void )
{
	// set the random seed 
	// SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 

	initialize_default_cell_definition(); 
	cell_defaults.functions.update_phenotype = NULL;
	// cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	// cell_defaults.functions.update_migration_bias = NULL; 

	if (parameters.bools("use_damage_model") == true)
	{
		cell_defaults.functions.pre_update_intracellular = pre_update_intracellular_drug_effect;
		std::cout<<"Using damage model"<<std::endl;
	}

	cell_defaults.functions.post_update_intracellular = post_update_intracellular_drug_effect; 
	// cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;

	// cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	// cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	cell_defaults.custom_data.add_variable(parameters.strings("node_to_visualize"), "dimensionless", 0.0 ); //for paraview visualization

	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 
	
	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/

	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries();

	/*
	   This summarizes the setup. 
	*/

	/*
       Cell rule definitions 
	*/

	setup_cell_rules(); 

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 

	// Cell_Definition* TLGL = find_cell_definition("TLGL");
	// Cell_Definition* TLGL_resistant = find_cell_definition("TLGL_resistant");
	// TLGL->functions.update_phenotype = transition_to_resistant_cell_type;
	// TLGL_resistant->functions.update_phenotype = transition_to_post_resistant_cell_type;
	display_cell_definitions( std::cout ); 


	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	// load cells from your CSV file
	load_cells_from_pugixml(); 	
}

void evaluateEffect(double impact, double half_max, double hill_power, Phenotype& phenotype, std::string effectNode)
{
    double random = UniformRandom();
    
    double hill_probability_response = Hill_response_function(impact, half_max, hill_power);
    if (hill_probability_response > random)
    {
        phenotype.intracellular->set_boolean_variable_value(effectNode, true);
		
    }
    else
    {
        phenotype.intracellular->set_boolean_variable_value(effectNode, false);
    }
}

void pre_update_intracellular_drug_effect(Cell* pCell, Phenotype& phenotype, double dt)
{

	double half_max = pCell->custom_data["half_max_effect_hill_function"];
	double hill_power = pCell->custom_data["exp_effect_hill_function"];
	double hill_probability_response;

	// ASSUMES TARGET INPUT NODES AND SUSBSTRATE NAMES ARE THE SAME!!!!!

	if (parameters.ints("num_substrates") == 1)
	{
		std::string damage_variable = parameters.strings("EffectNode1") + "_damage";
		double impact = pCell->custom_data[damage_variable];
		evaluateEffect(impact, half_max, hill_power, phenotype, parameters.strings("EffectNode1"));

	}

	else if (parameters.ints("num_substrates") == 2)
	{
		std::string damage1 = parameters.strings("EffectNode1") + "_damage";
		std::string damage2 = parameters.strings("EffectNode2") + "_damage";
		double impact1 = pCell->custom_data[damage1];
		double impact2 = pCell->custom_data[damage2];

		evaluateEffect(impact1, half_max, hill_power, phenotype, parameters.strings("EffectNode1"));
		evaluateEffect(impact2, half_max, hill_power, phenotype, parameters.strings("EffectNode2"));

	}

	else if (parameters.ints("num_substrates") == 3)
	{
		std::string damage1 = parameters.strings("EffectNode1") + "_damage";
		std::string damage2 = parameters.strings("EffectNode2") + "_damage";
		std::string damage3 = parameters.strings("EffectNode3") + "_damage";
		double impact1 = pCell->custom_data[damage1];
		double impact2 = pCell->custom_data[damage2];
		double impact3 = pCell->custom_data[damage3];

		evaluateEffect(impact1, half_max, hill_power, phenotype, parameters.strings("EffectNode1"));
		evaluateEffect(impact2, half_max, hill_power, phenotype, parameters.strings("EffectNode2"));
		evaluateEffect(impact3, half_max, hill_power, phenotype, parameters.strings("EffectNode3"));
	}

}

double update_value(double value, double previous_value, double smoothing) {
	double smoothed_value = (previous_value * smoothing + value)/(smoothing + 1);
	return smoothed_value;

	// I am going to have to go my own way on this ... And really try some kind of rolling average - I don't see how this will ever equal 1... 
	// Maybe multiple cell types IS easier. 
}

void evaluateResistance(Cell* pCell, Phenotype& phenotype, double dt)
{
	bool Apoptosis = pCell->phenotype.intracellular->get_boolean_variable_value("Apoptosis"); // ? 1.0 : 0.0;

	static int resistance_index = pCell->custom_data.find_variable_index("resistance");
	if (resistance_index < 0)
	{
		std::cout << "resistance variable not found" << std::endl;
		std::exit(-1);
	}

	double resitanceValue = pCell->custom_data["resistance"];
	
	// I think I can just autocast teh bool to a double
	pCell->custom_data["resistance"] = update_value(Apoptosis, resitanceValue, 1);
	std::cout<<"Resistance = "<<pCell->custom_data["resistance"]<<std::endl;



	// Then need to pass through hill funcction ... whats up there is more like 'update resistance likelyhood" - NEXT is evaluate (and thats where I'll have to figure out how to make it permanent)
	// How will I make the resistance permanent????? Some extxr int?? (more and more complicated ... )
	// double half_max = pCell->custom_data["half_max_effect_hill_function"];
	// double hill_power = pCell->custom_data["exp_effect_hill_function"];
	// double hill_probability_response;

	// std::string damage_variable = parameters.strings("EffectNode1") + "_damage";
	// double impact = pCell->custom_data[damage_variable];
	// double random = UniformRandom();
	// double hill_probability_response = Hill_response_function(impact, half_max, hill_power);
	// if (hill_probability_response > random)
	// {
	// 	return 1;
	// }
	// else
	// {
	// 	return 0;
	// }
}	

void post_update_intracellular_drug_effect(Cell* pCell, Phenotype& phenotype, double dt)
{
	bool Apoptosis = pCell->phenotype.intracellular->get_boolean_variable_value("Apoptosis"); // ? 1.0 : 0.0;
	bool Proliferation = pCell->phenotype.intracellular->get_boolean_variable_value("Proliferation"); // ? 1.0 : 0.0;

	
	double base_apoptosis = get_single_base_behavior(pCell, "apoptosis");
	// if(pCell->ID=1)
	// {
	// 	std::cout<<"Time: " << PhysiCell_globals.current_time <<std::endl;

	// 	std::cout<<"cell: " << pCell->type_name << " " << pCell->ID <<std::endl;
	// 	// std::getchar();
	// }
	// Evaluate to see if resistance likelihood should be increased
	// evaluateResistance(pCell, phenotype, dt);   
	
	// was using the logical test - not sure if I should use the value of the variable instead??? 
	// and why isn't the set single behavior working as expected?
	if(Apoptosis)
	{
		// std::cout<<"apoptosis = true"<<std::endl;
		// double base_cycle_entry = get_single_base_behavior(pCell, "cycle entry");
		// std::cout<<"base_cycle_entry = "<<base_cycle_entry<<std::endl;
		// set_single_behavior( pCell , "cycle entry" , base_cycle_entry  ); 

		// double base_apoptosis = get_single_base_behavior(pCell, "apoptosis");
		// std::cout<<"base_apoptosis = "<<base_apoptosis<<std::endl;

		double total_apoptosis_rate = base_apoptosis + pCell->custom_data["intervention_induced_apoptosis"];
		// std::cout<<"total_apoptosis_rate = "<<total_apoptosis_rate<<std::endl;

		set_single_behavior( pCell , "apoptosis" , total_apoptosis_rate  );
		// std::cout<<"adjusted_apoptosis = "<<get_single_behavior(pCell, "apoptosis") <<std::endl;
	}

	else
	{
		set_single_behavior( pCell , "apoptosis" , base_apoptosis  );
		// std::cout<<"apoptosis = false"<<std::endl;
	}

	// Resitance

	// from MaBoSS


	// if(Apoptosis == true && Proliferation == true)
	// {
	// 	double base_cycle_entry = get_single_base_behavior(pCell, "cycle entry");
	// 	std::cout<<"base_cycle_entry = "<<base_cycle_entry<<std::endl;
	// 	set_single_behavior( pCell , "cycle entry" , base_cycle_entry * pCell->custom_data["proliferation_multiplier"] ); 
		
	// 	double base_apoptosis = get_single_base_behavior(pCell, "apoptosis");
	// 	// std::cout<<"base_apoptosis = "<<base_apoptosis<<std::endl;
	// 	set_single_behavior( pCell , "apoptosis" , base_apoptosis * pCell->custom_data["apoptosis_multiplier"]  );
	// }

	// else if(Apoptosis == false && Proliferation == true)
	// {
	// 	double base_cycle_entry = get_single_base_behavior(pCell, "cycle entry");
	// 	std::cout<<"base_cycle_entry = "<<base_cycle_entry<<std::endl;
	// 	set_single_behavior( pCell , "cycle entry" , base_cycle_entry * pCell->custom_data["proliferation_multiplier"] ); 

	// 	set_single_behavior( pCell , "apoptosis" , 0  );
	// }

	// else if(Apoptosis == true && Proliferation == false)
	// {
	// 	set_single_behavior( pCell , "cycle entry" , 0  ); 

	// 	double base_apoptosis = get_single_base_behavior(pCell, "apoptosis");
	// 	std::cout<<"base_apoptosis = "<<base_apoptosis<<std::endl;
	// 	set_single_behavior( pCell , "apoptosis" , base_apoptosis * pCell->custom_data["apoptosis_multiplier"]  );
	// }

	// else
	// {
	// 	set_single_behavior( pCell , "cycle entry" , 0  ); 
	// 	set_single_behavior( pCell , "apoptosis" , 0  );
	// }
}

void add_compound( double drug_amount, double dose_interval, std::string substrate_name) 
{

	// Adds compounds uniformly to the microenvironment at concentration/amount "drug_amount" and every "drug_interval" minutes

    int number_of_voxels = microenvironment.mesh.voxels.size();

	// std::cout<<number_of_voxels<<" voxels"<<std::endl;
	
	int substrate_index = BioFVM::microenvironment.find_density_index( substrate_name ); 
	if (substrate_index < 0) 
    {
        // std::cout << "        static int << pro_GAP_index = " <<substrate_index << std::endl;
		std::cout << "        static int << " << substrate_name << "_index = " << substrate_index << std::endl;
        std::exit(-1);  //rwh: should really do these for each
    }

	std::cout<<"concentration in voxel 0 = "<<BioFVM::microenvironment.density_vector(0)[substrate_index]<<std::endl;

    for( int n=0; n < number_of_voxels ; n++ )
    {
		BioFVM::microenvironment.density_vector(n)[substrate_index] += drug_amount;
		// (*p_density_vectors)[n][ECM_density_index] = ecm.ecm_voxels[n].density;
		// std::cout<<BioFVM::microenvironment.density_vector(n)[ECM_anisotropy_index]<<std::endl;
		// std::cout<<&BioFVM::microenvironment.density_vector(n)[ECM_anisotropy_index]<<std::endl;
		// BioFVM::microenvironment.density_vector(n)[ECM_anisotropy_index] = ecm.ecm_voxels[n].anisotropy;
		// std::cout<<BioFVM::microenvironment.density_vector(n)[ECM_anisotropy_index]<<std::endl;
		// BioFVM::microenvironment.density_vector(n)[ECM_density_index] = ecm.ecm_voxels[n].density;
		// std::cout<<BioFVM::microenvironment.density_vector(n)[ECM_density_index]<<std::endl;
		// BioFVM::microenvironment.voxels[i].density_vector[ECM_density_index] = ecm.ecm_voxels[i].density;
		//  = ecm.ecm_voxels[i].anisotropy;
	
        // BioFVM::microenvironment.voxels[i].density_vector[ECM_density_index] = ecm.ecm_voxels[i].density;
    }


	std::cout<<"After addition concentration in voxel 0 = "<<BioFVM::microenvironment.density_vector(0)[substrate_index]<<std::endl;

    return;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Old - not using this anymore - 02.04.24 - replaced with XML (but probably chould have just used this instead of changing to XML....)


// Old - not using this anymore - 02.04.24
void transition_to_resistant_cell_type(Cell* pCell, Phenotype& phenotype, double dt)
{
	bool Apoptosis = pCell->phenotype.intracellular->get_boolean_variable_value("Apoptosis"); // ? 1.0 : 0.0;
	bool substrate = pCell->phenotype.intracellular->get_boolean_variable_value(parameters.strings("substrate_name")); // ? 1.0 : 0.0;

	// if(pro_GAP == true && Apoptosis == false)
	if(substrate == true)
	{
		// double base_apoptosis = get_single_base_behavior(pCell, "transform to TLGL_resistant");
		// std::cout<<"base_apoptosis = "<<base_apoptosis<<std::endl;
		set_single_behavior( pCell , "transform to TLGL_resistant" , 1E9  );
	}
	return;
}

// Old - not using this anymore - 02.04.24
void transition_to_post_resistant_cell_type( Cell* pCell, Phenotype& phenotype, double dt)
{
	bool Apoptosis = pCell->phenotype.intracellular->get_boolean_variable_value("Apoptosis"); // ? 1.0 : 0.0;
	bool substrate = pCell->phenotype.intracellular->get_boolean_variable_value(parameters.strings("substrate_name")); // ? 1.0 : 0.0;

	if(substrate == false && Apoptosis == false)
	{
		// double base_apoptosis = get_single_base_behavior(pCell, "transform to TLGL_post_resistant");
		// std::cout<<"base_apoptosis = "<<base_apoptosis<<std::endl;
		set_single_behavior( pCell , "transform to TLGL_post_resistant" , 1E9  );
	}
	return;
}

// Old - not using this anymore - 02.04.24
void from_nodes_to_cell(Cell* pCell, Phenotype& phenotype, double dt)
{
	double prosurvival_value = pCell->phenotype.intracellular->get_boolean_variable_value("Apoptosis") ? 1.0 : 0.0;
	std::cout<<"Apotosis = "<< pCell->phenotype.intracellular->get_boolean_variable_value("Apoptosis")<<std::endl;
	pCell->phenotype.intracellular->set_boolean_variable_value("Apoptosis", true);
	std::cout<<"Apotosis = "<<pCell->phenotype.intracellular->get_boolean_variable_value("Apoptosis")<<std::endl;
	static int start_phase_index; // Q_phase_index; 
	static int end_phase_index; // K_phase_index;
	double multiplier = 1.0;

	// live model 
			
	if( pCell->phenotype.cycle.model().code == PhysiCell_constants::live_cells_cycle_model )
	{
		start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
		end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );

		multiplier = ( ( prosurvival_value * 20 ) + 1 ); //[1, 21]
		phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = multiplier *	phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index);
		std::cout<< "prosurvival_value: " << prosurvival_value << " multiplier: " << multiplier << std::endl;
	}
	else
	{
		std::cout << "Warning: from_nodes_to_cell() is only implemented for live cells model. \n";
		std::cout << "         This function will not affect the transition rates for the current cell cycle model. \n";
	}

	pCell->set_internal_uptake_constants(dt); // why this? because we are not using the default uptake and secretion rates???
}

// Old - not using this anymore - 02.04.24
// void pre_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt )
// {
// 	if (PhysiCell::PhysiCell_globals.current_time >= 100.0 
// 		&& pCell->phenotype.intracellular->get_parameter_value("$time_scale") == 0.0
// 	){
// 		pCell->phenotype.intracellular->set_parameter_value("$time_scale", 0.1);
// 	}

// }
// Old - not using this anymore - 02.04.24
void post_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt )
{
	color_node(pCell);
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	std::vector< std::string > output( 4 , "rgb(0,0,0)" );
	
	if ( !pCell->phenotype.intracellular->get_boolean_variable_value( parameters.strings("node_to_visualize") ) )
	{
		output[0] = "rgb(255,0,0)"; // Red
		output[2] = "rgb(125,0,0)";
		
	}
	else{
		output[0] = "rgb(0, 255,0)"; // Green
		output[2] = "rgb(0, 125,0)";
	}
	
	return output;
}

void color_node(Cell* pCell){
	std::string node_name = parameters.strings("node_to_visualize");
	pCell->custom_data[node_name] = pCell->phenotype.intracellular->get_boolean_variable_value(node_name);
}
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
	
	// cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	// cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.pre_update_intracellular = pre_update_intracellular_drug_effect; 
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

	Cell_Definition* TLGL = find_cell_definition("TLGL");
	Cell_Definition* TLGL_resistant = find_cell_definition("TLGL_resistant");
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

void pre_update_intracellular_drug_effect(Cell* pCell, Phenotype& phenotype, double dt)
{
	Cell_Definition* pCD = find_cell_definition(pCell->type_name);	

	double drug_conc = get_single_signal(pCell, "pro_GAP");
	// std::cout<<"drug_conc = "<<drug_conc<<std::endl;
	if(drug_conc > pCell->custom_data["drug_threshold"])
	{
		pCell->phenotype.intracellular->set_boolean_variable_value("pro_GAP", true);
	}
	else
	{
		pCell->phenotype.intracellular->set_boolean_variable_value("pro_GAP", false);
	}
	
	// Without doing something with the accumulating substrate, internal values rises without end. So - skipping for now 
	// and assuming uptake is just proportional to the external concentration. I could do something fancy 
	// like slowing the rate of uptake as the internal concentration rises, but I don't think that's merited at this time.
}

void post_update_intracellular_drug_effect(Cell* pCell, Phenotype& phenotype, double dt)
{
	bool Apoptosis = pCell->phenotype.intracellular->get_boolean_variable_value("Apoptosis"); // ? 1.0 : 0.0;
	bool Proliferation = pCell->phenotype.intracellular->get_boolean_variable_value("Proliferation"); // ? 1.0 : 0.0;

	if(Apoptosis == true)
	{
		// double base_cycle_entry = get_single_base_behavior(pCell, "cycle entry");
		// std::cout<<"base_cycle_entry = "<<base_cycle_entry<<std::endl;
		// set_single_behavior( pCell , "cycle entry" , base_cycle_entry  ); 

		double base_apoptosis = get_single_base_behavior(pCell, "apoptosis");
		// std::cout<<"base_apoptosis = "<<base_apoptosis<<std::endl;
		set_single_behavior( pCell , "apoptosis" , base_apoptosis * pCell->custom_data["apoptosis_multiplier"]  );
	}

	else
	{
		set_single_behavior( pCell , "apoptosis" , 0  );
	}


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

void add_compound( double drug_amount, double dose_interval ) 
{

	// Adds compounds uniformly to the microenvironment at concentration/amount "drug_amount" and every "drug_interval" minutes

    int number_of_voxels = microenvironment.mesh.voxels.size();

	// std::cout<<number_of_voxels<<" voxels"<<std::endl;
	
	static int pro_GAP_index = BioFVM::microenvironment.find_density_index( "pro_GAP" ); 
	if (pro_GAP_index < 0) 
    {
        std::cout << "        static int pro_GAP_index = " <<pro_GAP_index << std::endl;
        std::exit(-1);  //rwh: should really do these for each
    }

    for( int n=0; n < number_of_voxels ; n++ )
    {
		BioFVM::microenvironment.density_vector(n)[pro_GAP_index] = drug_amount;
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
    return;

}

void transition_to_resistant_cell_type(Cell* pCell, Phenotype& phenotype, double dt)
{
	bool Apoptosis = pCell->phenotype.intracellular->get_boolean_variable_value("Apoptosis"); // ? 1.0 : 0.0;
	bool pro_GAP = pCell->phenotype.intracellular->get_boolean_variable_value("pro_GAP"); // ? 1.0 : 0.0;

	// if(pro_GAP == true && Apoptosis == false)
	if(pro_GAP == true)
	{
		// double base_apoptosis = get_single_base_behavior(pCell, "transform to TLGL_resistant");
		// std::cout<<"base_apoptosis = "<<base_apoptosis<<std::endl;
		set_single_behavior( pCell , "transform to TLGL_resistant" , 1E9  );
	}
	return;
}

void transition_to_post_resistant_cell_type( Cell* pCell, Phenotype& phenotype, double dt )
{
	bool Apoptosis = pCell->phenotype.intracellular->get_boolean_variable_value("Apoptosis"); // ? 1.0 : 0.0;
	bool pro_GAP = pCell->phenotype.intracellular->get_boolean_variable_value("pro_GAP"); // ? 1.0 : 0.0;

	if(pro_GAP == false && Apoptosis == false)
	{
		// double base_apoptosis = get_single_base_behavior(pCell, "transform to TLGL_post_resistant");
		// std::cout<<"base_apoptosis = "<<base_apoptosis<<std::endl;
		set_single_behavior( pCell , "transform to TLGL_post_resistant" , 1E9  );
	}
	return;
}

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

void pre_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt )
{
	if (PhysiCell::PhysiCell_globals.current_time >= 100.0 
		&& pCell->phenotype.intracellular->get_parameter_value("$time_scale") == 0.0
	){
		pCell->phenotype.intracellular->set_parameter_value("$time_scale", 0.1);
	}

}

void post_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt )
{
	color_node(pCell);
}

// std::vector<std::string> my_coloring_function( Cell* pCell )
// {
// 	std::vector< std::string > output( 4 , "rgb(0,0,0)" );
	
// 	if ( !pCell->phenotype.intracellular->get_boolean_variable_value( parameters.strings("node_to_visualize") ) )
// 	{
// 		output[0] = "rgb(255,0,0)"; // Red
// 		output[2] = "rgb(125,0,0)";
		
// 	}
// 	else{
// 		output[0] = "rgb(0, 255,0)"; // Green
// 		output[2] = "rgb(0, 125,0)";
// 	}
	
// 	return output;
// }

void color_node(Cell* pCell){
	std::string node_name = parameters.strings("node_to_visualize");
	pCell->custom_data[node_name] = pCell->phenotype.intracellular->get_boolean_variable_value(node_name);
}
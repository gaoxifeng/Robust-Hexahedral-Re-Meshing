//    This file is part of the implementation of

//    Robust Structure Simplification for Hex Re-meshing
//    Xifeng Gao, Daniele Panozzo, Wenping Wang, Zhigang Deng, Guoning Chen
//    In ACM Transactions on Graphics (Proceedings of SIGGRAPH ASIA 2017)
// 
// Copyright (C) 2017 Xifeng Gao<gxf.xisha@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.

#include "simplification.h"
char path_IOH[300];
char Choices[300]="SIM";
char Hex_NUM_Ratio[300]="1.0";
char Iteration_Base[300] = "2";
char ToBe_Removed_Cuboid_Ratio[300] = "0.9";
char Hard_Feature[300] = "1";
char Hausdorff_ratio_t[300] = "0.01";
char temp_string[300];
h_io io;
base_complex bc;
simplification sim;
int main( int argc, char* argv[] )
{
	//int nprocess = -1;
	//tbb::task_scheduler_init init(nprocess == -1 ? tbb::task_scheduler_init::automatic: nprocess);

	if (strcmp(Choices, "SIM") == 0 || strcmp(Choices, "OPT") == 0) {
		if (argc != 7) {
			cout << "#parameters are not exactly 7!" << endl;
		}
		sprintf(Hex_NUM_Ratio, "%s", argv[2]);
		sprintf(Iteration_Base, "%s", argv[3]);
		sprintf(ToBe_Removed_Cuboid_Ratio, "%s", argv[4]);
		sprintf(Hard_Feature, "%s", argv[5]);
		sprintf(path_IOH, "%s", argv[6]);
		if(argc == 8) sprintf(Hausdorff_ratio_t, "%s", argv[7]);
	}
	sprintf(path_out,"%s",path_IOH);

	sim.mesh.type = Mesh_type::Hex;
	io.read_hybrid_mesh_VTK(sim.mesh, path_IOH);
	build_connectivity(sim.mesh);
	bc.singularity_structure(sim.si, sim.mesh);
	bc.base_complex_extraction(sim.si, sim.frame, sim.mesh);

	double hex_num_ratiod = 0.9, tobe_removed_cuboid_ratio = 1.0, hausdorff_ratio_td = 0.01;
	int iteration_base = 3;
	bool hard_feature = false;

	std::string::size_type sz;
	if (strcmp(Hex_NUM_Ratio, " ") != 0)
		hex_num_ratiod = std::stod(Hex_NUM_Ratio, &sz);

	if (strcmp(Iteration_Base, " ") != 0)
		iteration_base = std::stoi(Iteration_Base, &sz);

	if (strcmp(ToBe_Removed_Cuboid_Ratio, " ") != 0)
		tobe_removed_cuboid_ratio = std::stod(ToBe_Removed_Cuboid_Ratio, &sz);

	if (strcmp(Hard_Feature, " ") != 0)
		hard_feature = std::stoi(Hard_Feature, &sz);

	if (strcmp(Hausdorff_ratio_t, " ") != 0)
		hausdorff_ratio_td = std::stod(Hausdorff_ratio_t, &sz);

	if (strcmp(Choices, "SIM") == 0){
		//simplification
		sim.set_sharp_feature(hard_feature);
		sim.set_slim_iteration_base(iteration_base);
		sim.set_cuboid_ratio(tobe_removed_cuboid_ratio);
		if (hex_num_ratiod > 10)sim.set_target_hex_num(hex_num_ratiod);
		else sim.set_target_hex_num(hex_num_ratiod * sim.mesh.Hs.size());
		sim.hausdorff_ratio_threshould = hausdorff_ratio_td;

		if (!sim.initialize()) return false;
		sim.pipeline();
	}
	else if (strcmp(Choices, "OPT") == 0) {
		//optimization
		sim.set_sharp_feature(hard_feature);
		sim.set_slim_iteration_base(iteration_base);
		sim.set_cuboid_ratio(tobe_removed_cuboid_ratio);
		if (hex_num_ratiod > 10)sim.set_target_hex_num(hex_num_ratiod);
		else sim.set_target_hex_num(hex_num_ratiod * sim.mesh.Hs.size());

		if (!sim.initialize()) return false;
		sim.optimization();
		char path[300];
		sprintf(path, "%s%s", path_IOH, "_optimized.vtk");
		io.write_hybrid_mesh_VTK(sim.mesh, path);
	}
	
	return 0;
}
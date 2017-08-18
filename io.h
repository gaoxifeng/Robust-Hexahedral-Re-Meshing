#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "global_types.h"
#include "global_functions.h"

using namespace std;
class h_io
{
private:
	int counter;
public:
	h_io(void) { counter = -1; }; ~h_io(void) {};

	void read_hybrid_mesh_OFF(Mesh &hmi, char *path);
	void write_hybrid_mesh_OFF(Mesh &hmi, char *path);

	void read_hybrid_mesh_VTK(Mesh &hmi, char * path);
	void write_hybrid_mesh_VTK(Mesh &hmi, char * path);
	void write_hybrid_mesh_VTK_ele_tag(Mesh &hmi, Eigen::VectorXi &ele_tag, char * path);
	void write_hybrid_mesh_VTK_ele_tag(Mesh &hmi, Eigen::VectorXd &ele_tag, char * path);

	void read_hybrid_mesh_MESH(Mesh &hmi, char * path);
	void write_hybrid_mesh_MESH(Mesh &hmi, char * path);

	void write_hybrid_mesh_OBJ(Mesh &hmi, char * path);
	//singularity structure
	void write_singularG_VTK(Singularity &si, Mesh &mesh, char *path);
	//base complex
	void write_Frame_VTK(Frame &frame, Mesh &mesh, char *path);
	//sheet
	void write_Sheet_VTK(Sheet &s, Mesh &mesh, Frame &frame, char *path);
	void write_Points_VTK(MatrixXd &V, char *path);
	void write_Vs_Groups_VTK(Mesh &mesh, vector<vector<uint32_t>> &Vs_Group, char *path, bool twoside);
	//chord
	void write_Chord_VTK(CHord &c, Mesh &mesh, Frame &frame, char *path);

};


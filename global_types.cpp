#pragma once
#include "global_types.h"

//meshes
Mesh_Feature mf;


//local mesh
vector<MatrixXi> LocalMeshTri_F;
vector<MatrixXd> LocalMeshTri_V;
vector<VectorXd> LocalMeshTri_C;
vector<MatrixXi> LocalMeshEs;
//local projection
vector<Feature_Constraints> SeLocalProjection;
//mesh
vector<MatrixXi> MeshTri_F;
vector<MatrixXd> MeshTri_V, MeshTri_N, MeshTri_VN;
vector<VectorXd> MeshTri_C;
vector<MatrixXi> MeshEs;

vector<MatrixXi> SkeletonTri_F;
vector<MatrixXd> SkeletonTri_V;
vector<MatrixXd> SkeletonTri_N;
vector<MatrixXd> SkeletonTri_C;

vector<MatrixXi> FrameTri_F;
vector<MatrixXd> FrameTri_V, FrameTri_N;
vector<VectorXd> FrameTri_C;

vector<MatrixXi> MeshSkeletonTri_F;
vector<MatrixXd> MeshSkeletonTri_V, MeshSkeletonTri_N, MeshSkeletonTri_VN;
vector<VectorXd> MeshSkeletonTri_C;
vector<VectorXd> AOSkeletons;
vector<VectorXd> AOs;

vector<MatrixXi> CollapseTri_F;
vector<MatrixXd> CollapseTri_V, CollapseTri_N;
vector<VectorXd> CollapseTri_C;



char path_out[300];
int32_t GRAIN_SIZE;

MatrixXd tenC, tenC_perm;

double diagonal_len;

Mesh mesh_sheet, mesh_sheetS;
bool HEXAHEDRAL_COLLASPE = false;
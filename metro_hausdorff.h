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

#pragma once
#define _SILENCE_STDEXT_HASH_DEPRECATION_WARNINGS 1
#define _CRT_SECURE_NO_WARNINGS 1
#include <time.h>
#include <vcg/math/histogram.h>
#include <vcg/complex/complex.h>
#include <vcg/simplex/face/component_ep.h>
//#include <wrap/io_trimesh/import.h>
//#include <wrap/io_trimesh/export.h>
#include <vcg/complex/algorithms/update/component_ep.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include "sampling.h"
#include "global_types.h"
//#include <wrap/io_trimesh/import_off.h>

using namespace std;
// project definitions.
// error messages




class CFace;
class CVertex;
struct UsedTypes :public vcg::UsedTypes< vcg::Use<CFace>::AsFaceType, vcg::Use<CVertex>::AsVertexType> {};
class CVertex : public vcg::Vertex<UsedTypes, vcg::vertex::Coord3d, vcg::vertex::Qualityf, vcg::vertex::Normal3d, vcg::vertex::Color4b, vcg::vertex::BitFlags> {};
class CFace : public vcg::Face< UsedTypes, vcg::face::VertexRef, vcg::face::Normal3d, vcg::face::EdgePlane, vcg::face::Color4b, vcg::face::Mark, vcg::face::BitFlags> {};
class CMesh : public vcg::tri::TriMesh< std::vector<CVertex>, std::vector<CFace> > {};


// -----------------------------------------------------------------------------------------------

using namespace vcg;

extern int compute(Mesh &mesh0, Mesh & mesh1, double &hausdorff_ratio, double &hausdorff_ratio_threshold);
#include "io.h"

void h_io::read_hybrid_mesh_OFF(Mesh &hmi, char *path)
{
	hmi.Vs.clear(); hmi.Fs.clear(); hmi.Hs.clear();

	char fpath[300];
	sprintf(fpath, "%s%s", path, ".off");
	std::fstream file(fpath, std::ios::in);
	char s[1024], sread[1024];
	uint32_t vnum, tnum;	float x, y, z;
	file.getline(s, 1023); file.getline(s, 1023);
	sscanf(s, "%d %d %s", &vnum, &tnum, &sread);
	hmi.V.resize(3, vnum);
	for (uint32_t i = 0; i<vnum; i++){
		file.getline(s, 1023);
		//sscanf(s, "%f %f %f", &x, &y, &z);
		int ind;
		sscanf(s, "%f %f %f %d", &x, &y, &z, &ind);
		Hybrid_V v; v.v.resize(3);
		v.v[0] = x;
		v.v[1] = y;
		v.v[2] = z;
		hmi.V.col(i) = Vector3d(x, y, z);
		v.id = i;
		v.boundary = false;
		hmi.Vs.push_back(v);
	}

	if (hmi.type == Mesh_type::Tri || hmi.type == Mesh_type::Qua){
		Hybrid_F f;
		if (hmi.type == Mesh_type::Tri) f.vs.resize(3); else f.vs.resize(4);
		uint32_t a, b, c, d, e;
		hmi.Fs.resize(tnum);
		for (uint32_t i = 0; i<tnum; i++) {
			file.getline(s, 1023);
			if (sscanf(s, "%d %d %d %d", &vnum, &a, &b, &c) == 4 && vnum == 3){
				//a--;b--;c--;//temporarily
				f.vs[0]=a;
				f.vs[1]=b;
				f.vs[2]=c;
			}
			else if (sscanf(s, "%d %d %d %d %d", &vnum, &a, &b, &c, &d) == 5 && vnum == 4){
				//a--;b--;c--;d--;//temporarily
				f.vs[0] = a;
				f.vs[1] = b;
				f.vs[2] = c;
				f.vs[2] = d;
			}
			else { std::cout << "Wrong format of input file!" << endl; system("PAUSE"); }

			f.id = i; hmi.Fs[i] = f;
			for (uint32_t i = 0; i < f.vs.size(); i++) hmi.Vs[f.vs[i]].neighbor_fs.push_back(f.id);
		}
	}
	else if(hmi.type == Mesh_type::Tet || hmi.type == Mesh_type::Hyb || hmi.type == Mesh_type::Hex)
	{
		uint32_t a, b, c, d, e, f, g, m,o,p,q;
		Hybrid h; 
		hmi.Hs.resize(tnum);
		for (uint32_t i = 0; i<tnum; i++) {
			file.getline(s, 1023);
		//while (file.getline(s, 1023)){
			if (sscanf(s, "%d %d %d %d %d", &vnum, &a, &b, &c, &d) == 5 && vnum == 4){
				h.vs.resize(4);
				//a--;b--;c--;d--;//temporarily
				h.vs[0]=a;
				h.vs[1]=b;
				h.vs[2]=c;
				h.vs[3]=d;
			}else if (sscanf(s, "%d %d %d %d %d %d", &vnum, &a, &b, &c, &d, &e) == 6 && vnum == 5){
				h.vs.resize(5);
				h.vs[0] = a;
				h.vs[1] = b;
				h.vs[2] = c;
				h.vs[3] = d;
				h.vs[4] = e;
			}
			else if (sscanf(s, "%d %d %d %d %d %d %d", &vnum, &a, &b, &c, &d, &e, &f) == 7 && vnum == 6)
			{
				h.vs.resize(6);
				h.vs[0] = a;
				h.vs[1] = b;
				h.vs[2] = c;
				h.vs[3] = d;
				h.vs[4] = e;
				h.vs[5] = f;
			}
			else if (sscanf(s, "%d %d %d %d %d %d %d %d %d", &vnum, &a, &b, &c, &d, &e, &f, &g, &m) == 9 && vnum == 8)
			{
				h.vs.resize(8);
				h.vs[0] = a;
				h.vs[1] = b;
				h.vs[2] = c;
				h.vs[3] = d;
				h.vs[4] = e;
				h.vs[5] = f;
				h.vs[6] = g;
				h.vs[7] = m;
			}
			else if (sscanf(s, "%d %d %d %d %d %d %d %d %d %d %d", &vnum, &a, &b, &c, &d, &e, &f, &g, &m, &o, &p) == 11 && vnum == 10)
			{
				h.vs.resize(8);
				h.vs[0] = a;
				h.vs[1] = b;
				h.vs[2] = c;
				h.vs[3] = d;
				h.vs[4] = e;
				h.vs[5] = f;
				h.vs[6] = g;
				h.vs[7] = m;
			}
			else if (sscanf(s, "%d %d %d %d %d %d %d %d %d %d %d %d", &vnum, &a, &b, &c, &d, &e, &f, &g, &m, &o, &p, &q) == 12 && vnum == 10)
			{
				h.vs.resize(8);
				h.vs[0] = a;
				h.vs[1] = b;
				h.vs[2] = c;
				h.vs[3] = d;
				h.vs[4] = e;
				h.vs[5] = f;
				h.vs[6] = g;
				h.vs[7] = m;
			}
			else if (sscanf(s, "%d %d %d %d %d %d %d %d %d", &vnum, &a, &b, &c, &d, &e, &f, &g, &m, &o, &p, &q) == 9 && vnum == 10)
			{
				h.vs.resize(8);
				h.vs[0] = a;
				h.vs[1] = b;
				h.vs[2] = c;
				h.vs[3] = d;
				h.vs[4] = e;
				h.vs[5] = f;
				h.vs[6] = g;
				h.vs[7] = m;
			}
			else { std::cout << "Wrong format of input file!" << endl; system("PAUSE"); }

			h.id = i; hmi.Hs[i] = h;
			for (uint32_t i = 0; i < h.vs.size(); i++) hmi.Vs[h.vs[i]].neighbor_hs.push_back(h.id);
		}
	}
	file.close();
}
void h_io::write_hybrid_mesh_OFF(Mesh &hmi, char *path)
{
	fstream f(path, ios::out);

	f << "OFF" << endl;
	if(hmi.type == Mesh_type::Tri || hmi.type == Mesh_type::Qua) 
		f << hmi.Vs.size() << " " << hmi.Fs.size() << " " << 0 << endl;
	else f << hmi.Vs.size() << " " << hmi.Hs.size() << " " << 0 << endl;
	for (uint32_t i = 0; i<hmi.Vs.size(); i++) 
		f << setprecision(10) << hmi.Vs[i].v[0] << " " << setprecision(10) << hmi.Vs[i].v[1] << " " << setprecision(10) << hmi.Vs[i].v[2] << endl;

	if (hmi.type == Mesh_type::Tri || hmi.type == Mesh_type::Qua) {
		for (uint32_t i = 0; i < hmi.Fs.size(); i++){
			f << hmi.Fs[i].vs.size() << " ";
			for (uint32_t j = 0; j < hmi.Fs[i].vs.size(); j++) f << hmi.Fs[i].vs[j] << " ";
			f << endl;
		}
	}
	else {
		for (uint32_t i = 0; i < hmi.Hs.size(); i++){
			f << hmi.Hs[i].vs.size() << " ";
			for (uint32_t j = 0; j < hmi.Hs[i].vs.size(); j++) f << hmi.Hs[i].vs[j] +1<< " ";
			f << endl;
		}
	}
	f.close();
}
void h_io::read_hybrid_mesh_VTK(Mesh &hmi, char * path)
{
	char file[300];
	sprintf(file, "%s%s", path, ".vtk");
	std::ifstream ff(file, std::ios::in);
	char s[1024], sread[1024], sread2[1024];
	uint32_t vnum, hnum;	double x, y, z;

	bool find = false; uint32_t lines = 0;
	while (!find)
	{
		ff.getline(s, 1023);
		if (sscanf(s, "%s %d %s", &sread, &vnum, &sread2) == 3 && (strcmp(sread, "POINTS") == 0))
			find = true;
		if (++lines>10)
			throw std::runtime_error("cannot find head of VTK!");
	}
	hmi.V.resize(3, vnum);
	hmi.V.setZero();
	hmi.Vs.resize(vnum);
	for (uint32_t i = 0; i<vnum; i++)
	{
		ff.getline(s, 1023);
		sscanf(s, "%lf %lf %lf", &x, &y, &z);

		hmi.V(0, i) = x;
		hmi.V(1, i) = y;
		hmi.V(2, i) = z;

		Hybrid_V v;
		v.id = i; v.boundary = false;
		hmi.Vs[i] = v;
	}

	find = false;
	while (!find)
	{
		uint32_t temp_int;
		ff.getline(s, 1023);
		if (sscanf(s, "%s %d %d", &sread, &hnum, &temp_int) == 3 && (strcmp(sread, "CELLS") == 0))
			find = true;
	}
	hmi.Hs.resize(hnum); 
	Hybrid h;
	uint32_t a, b, c, d, e, f, g, m, o, p, q;
	for (uint32_t i = 0; i<hnum; i++)
	//while (ff.getline(s, 1023))
	{
		ff.getline(s, 1023);
		if (sscanf(s, "%d %d %d %d %d", &vnum, &a, &b, &c, &d) == 5 && vnum == 4) {
			h.vs.resize(4);
			//a--;b--;c--;d--;//temporarily
			h.vs[0] = a;
			h.vs[1] = b;
			h.vs[2] = c;
			h.vs[3] = d;
		}
		else if (sscanf(s, "%d %d %d %d %d %d", &vnum, &a, &b, &c, &d, &e) == 6 && vnum == 5) {
			h.vs.resize(5);
			h.vs[0] = a;
			h.vs[1] = b;
			h.vs[2] = c;
			h.vs[3] = d;
			h.vs[4] = e;
		}
		else if (sscanf(s, "%d %d %d %d %d %d %d", &vnum, &a, &b, &c, &d, &e, &f) == 7 && vnum == 6)
		{
			h.vs.resize(6);
			h.vs[0] = a;
			h.vs[1] = b;
			h.vs[2] = c;
			h.vs[3] = d;
			h.vs[4] = e;
			h.vs[5] = f;
		}
		else if (sscanf(s, "%d %d %d %d %d %d %d %d %d", &vnum, &a, &b, &c, &d, &e, &f, &g, &m) == 9 && vnum == 8)
		{
			h.vs.resize(8);
			h.vs[0] = a;
			h.vs[1] = b;
			h.vs[2] = c;
			h.vs[3] = d;
			h.vs[4] = e;
			h.vs[5] = f;
			h.vs[6] = g;
			h.vs[7] = m;
			//std::reverse(h.vs.begin(), h.vs.begin() + 4);
			//std::reverse(h.vs.begin() + 4, h.vs.end());
		}
		else { 
			std::cout << "Wrong format of input file!" << endl; system("PAUSE"); 
		}

		h.id = i; hmi.Hs[h.id] = h;
		for (uint32_t i = 0; i < h.vs.size(); i++) hmi.Vs[h.vs[i]].neighbor_hs.push_back(h.id);
	}
	ff.close();
}
void h_io::write_hybrid_mesh_VTK(Mesh &hmi, char * path)
{
	std::fstream f (path, std::ios::out);

	f << "# vtk DataFile Version 2.0" << std::endl << "mesh vtk data - converted from .off" << std::endl;
	f << "ASCII" << std::endl;
	f << "DATASET UNSTRUCTURED_GRID" << std::endl;

	f << "POINTS " << hmi.V.cols() << " double" << std::endl;
	for (uint32_t i = 0; i<hmi.V.cols(); i++)
		f << hmi.V(0, i) << " " << hmi.V(1, i) << " " << hmi.V(2, i) << std::endl;

	if (hmi.type == Mesh_type::Tri || hmi.type == Mesh_type::Qua) {
		uint32_t vnum = hmi.Fs[0].vs.size();
		f << "CELLS " << hmi.Fs.size() << " " << hmi.Fs.size() * (vnum + 1) << std::endl;

		for (uint32_t i = 0; i < hmi.Fs.size(); i++) {
			f << " " << vnum << " ";
			for (uint32_t j = 0; j < vnum; j++) f << hmi.Fs[i].vs[j] << " ";
			f << std::endl;
		}
		f << "CELL_TYPES " << hmi.Fs.size() << std::endl;
		for (uint32_t i = 0; i < hmi.Fs.size(); i++)
			if (hmi.type == 0) f << 5 << std::endl; else f << 9 << std::endl;
	}
	else {
		uint32_t vnum = hmi.Hs[0].vs.size();
		f << "CELLS " << hmi.Hs.size() << " " << hmi.Hs.size() * (vnum + 1) << std::endl;

		for (uint32_t i = 0; i < hmi.Hs.size(); i++)
		{
			f << " " << vnum << " ";
			for (uint32_t j = 0; j < vnum; j++)
				f << hmi.Hs[i].vs[j] << " ";
			/*f << hmi.Hs[i].vs[3] << " ";
			f << hmi.Hs[i].vs[2] << " ";
			f << hmi.Hs[i].vs[1] << " ";
			f << hmi.Hs[i].vs[0] << " ";
			f << hmi.Hs[i].vs[7] << " ";
			f << hmi.Hs[i].vs[6] << " ";
			f << hmi.Hs[i].vs[5] << " ";
			f << hmi.Hs[i].vs[4] << " ";*/
			f<< std::endl;
		}
		f << "CELL_TYPES " << hmi.Hs.size() << std::endl;
		for (uint32_t i = 0; i < hmi.Hs.size(); i++)
			if(hmi.type== Mesh_type::Tet)
				f << 10 << std::endl;
			else
				f << 12 << std::endl;
	}


	f << "POINT_DATA " << hmi.Vs.size() << std::endl;
	f << "SCALARS fixed int" << std::endl;
	f << "LOOKUP_TABLE default" << std::endl;
	for (uint32_t i = 0; i<hmi.Vs.size(); i++) {
		if (hmi.Vs[i].boundary) f << "1" << std::endl; else f << "0" << std::endl;
	}
	f.close();
}
void h_io::write_hybrid_mesh_VTK_ele_tag(Mesh &hmi, Eigen::VectorXi &ele_tag, char * path)
{

}
void h_io::write_hybrid_mesh_VTK_ele_tag(Mesh &hmi, Eigen::VectorXd &ele_tag, char * path)
{

}
void h_io::read_hybrid_mesh_MESH(Mesh &hmi, char * path)
{
	hmi.Vs.clear(); hmi.Hs.clear();
	char file[300];
	sprintf(file, "%s%s", path, ".mesh");
	std::fstream f(file, std::ios::in);
	char s[1024], sread[1024];
	int vnum, hnum;	double x, y, z;

	int find = false;
	while (!find)
	{
		f.getline(s, 1023);
		if (sscanf(s, "%s%d", &sread, &vnum) == 2 && (strcmp(sread, "Vertices") == 0))
		{
			find = true;
		}
		else if(sscanf(s, "%s", &sread) == 1 && (strcmp(sread, "Vertices") == 0))
		{
			find = true;
			f.getline(s, 1023);
			sscanf(s, "%d", &vnum);
		}
	}
	hmi.V.resize(3, vnum);
	hmi.V.setZero();
	hmi.Vs.resize(vnum);

	for (int i = 0; i<vnum; i++)
	{
		f.getline(s, 1023);
		sscanf(s, "%lf %lf %lf", &x, &y, &z);

		hmi.V(0, i) = x;
		hmi.V(1, i) = y;
		hmi.V(2, i) = z;

		Hybrid_V v;
		v.id = i; v.boundary = false;
		hmi.Vs[i] = v;
	}
	find = false;
	while (!find)
	{
		int temp_int;
		f.getline(s, 1023);
		if (sscanf(s, "%s%d", &sread, &hnum) == 2 && (strcmp(sread, "Hexahedra") == 0))
		{
			find = true;
		}
		else if (sscanf(s, "%s", &sread) == 1 && (strcmp(sread, "Hexahedra") == 0))
		{
			f.getline(s, 1023);
			sscanf(s, "%d", &hnum);
			find = true;
		}
	}
	hmi.Hs.resize(hnum);
	Hybrid h;
	int32_t a, b, c, d, e, ff, g, m, o, p, q;

	int hid = 0;
	for (int i = 0; i<hnum; i++)
	{
		f.getline(s, 1023);
		int a, b, c, d, e, f, g, m;
		//sscanf(s,"%d %d %d %d %d %d %d %d %d",&vnum,&a,&b,&c,&d,&e,&f,&g,&m);
		sscanf(s, "%d %d %d %d %d %d %d %d %d", &a, &b, &c, &d, &e, &ff, &g, &m, &vnum);
		
		a--; b--; c--; d--; e--; ff--; g--; m--;

		h.vs.resize(8);
		h.vs[0] = a;
		h.vs[1] = b;
		h.vs[2] = c;
		h.vs[3] = d;
		h.vs[4] = e;
		h.vs[5] = ff;
		h.vs[6] = g;
		h.vs[7] = m;

		h.id = i; hmi.Hs[h.id] = h;
		for (uint32_t i = 0; i < h.vs.size(); i++) hmi.Vs[h.vs[i]].neighbor_hs.push_back(h.id);
	}

	f.close();
}
void h_io::write_hybrid_mesh_MESH(Mesh &hmi, char * path)
{
	std::fstream f(path, std::ios::out);

	f << "MeshVersionFormatted 1" << std::endl;
	f << "Dimension 3" << std::endl;
	f << "Vertices" << " " << hmi.V.cols() << std::endl;
	for (int i = 0; i<hmi.V.cols(); i++)
		f<< hmi.V(0, i) << " " << hmi.V(1, i) << " " << hmi.V(2, i) << " " << 0 << std::endl;

	if (hmi.type == Mesh_type::Tri) {
		f<< "Triangles" << endl;
		f<< hmi.Fs.size() << std::endl;

		for (int i = 0; i<hmi.Fs.size(); i++){
			f<< hmi.Fs[i].vs[0] + 1 << " " << hmi.Fs[i].vs[1] + 1 << " " << hmi.Fs[i].vs[2] + 1 << " " << 0 << std::endl;
		}
	}else if(hmi.type == Mesh_type::Hex) {
		f<< "Hexahedra" << endl;
		f<< hmi.Hs.size() << std::endl;

		for (int i = 0; i < hmi.Hs.size(); i++) {
			for (auto vid : hmi.Hs[i].vs)
				f<< vid + 1 << " ";
			f<<0<<  endl;
		}
	}
	f << "End";
	f.close();
}

void h_io::write_hybrid_mesh_OBJ(Mesh &hmi, char * path) {
	std::fstream f(path, std::ios::out);
	for (int i = 0; i<hmi.V.cols(); i++)
		f << "v " << hmi.V(0, i) << " " << hmi.V(1, i) << " " << hmi.V(2, i) << std::endl;

	for (int i = 0; i < hmi.Fs.size(); i++) {
		f << "f";
		for (int j = 0; j < hmi.Fs[i].vs.size(); j++)
			f << " " << hmi.Fs[i].vs[j] + 1;
		f << endl;
	}
	f.close();
}

void h_io::write_singularG_VTK(Singularity &si, Mesh &mesh, char *path)
{
	std::fstream ff(path, std::ios::out);
	ff << "# vtk DataFile Version 2.0" << std::endl << "mesh vtk data" << std::endl;
	ff << "ASCII" << std::endl;
	ff << "DATASET POLYDATA" << std::endl;


	ff << "POINTS " << mesh.Vs.size() << " double" << std::endl;
	for (uint32_t i = 0; i<mesh.Vs.size(); i++)
		ff << mesh.V(0, i) << " " << mesh.V(1, i) << " " << mesh.V(2, i) << std::endl;

	ff << "VERTICES " << si.SVs.size() << " " << si.SVs.size() * 2 << endl;
	for (uint32_t i = 0; i<si.SVs.size(); i++) ff << "1 " << si.SVs[i].hid << endl;

	int line_num = 0;
	for (uint32_t i = 0; i<si.SEs.size(); i++) line_num += si.SEs[i].vs_link.size();

	ff << "LINES " << si.SEs.size() << " " << line_num + si.SEs.size() << endl;
	for (uint32_t i = 0; i<si.SEs.size(); i++) {
		ff << si.SEs[i].vs_link.size() << " ";
		for (uint32_t j = 0; j<si.SEs[i].vs_link.size(); j++) ff << si.SEs[i].vs_link[j] << " ";
		ff << endl;
	}


	ff << "POINT_DATA " << mesh.Vs.size() << std::endl;
	ff << "SCALARS V_Scalars int" << std::endl;
	ff << "LOOKUP_TABLE V_Table" << std::endl;
	std::vector<short> V_tag(mesh.Vs.size(), 0);
	for (uint32_t i = 0; i < mesh.Vs.size(); i++) if (mesh.Vs[i].boundary) V_tag[i] = 1;
	for (uint32_t j = 0; j<si.SVs.size(); j++) V_tag[si.SVs[j].hid] = 3;
	for (uint32_t j = 0; j<si.SEs.size(); j++) for (uint32_t k = 0; k<si.SEs[j].vs_link.size(); k++)
		V_tag[si.SEs[j].vs_link[k]] = 3;

	for (uint32_t i = 0; i<V_tag.size(); i++) ff << V_tag[i] << std::endl;

	ff.close();
}
void h_io::write_Frame_VTK(Frame &frame, Mesh &mesh, char *path)
{
	std::fstream ff(path, std::ios::out);
	ff << "# vtk DataFile Version 2.0" << std::endl << "mesh vtk data" << std::endl;
	ff << "ASCII" << std::endl;
	ff << "DATASET POLYDATA" << std::endl;


	ff << "POINTS " << mesh.Vs.size() << " double" << std::endl;
	for (uint32_t i = 0; i<mesh.Vs.size(); i++)
		ff << mesh.V(0, i) << " " << mesh.V(1, i) << " " << mesh.V(2, i) << std::endl;

	ff << "VERTICES " << frame.FVs.size() << " " << frame.FVs.size() * 2 << endl;
	for (uint32_t i = 0; i<frame.FVs.size(); i++) ff << "1 " << frame.FVs[i].hid << endl;

	int line_num = 0;
	for (uint32_t i = 0; i<frame.FEs.size(); i++) line_num += frame.FEs[i].vs_link.size();

	ff << "LINES " << frame.FEs.size() << " " << line_num + frame.FEs.size() << endl;
	for (uint32_t i = 0; i<frame.FEs.size(); i++) {
		ff << frame.FEs[i].vs_link.size() << " ";
		for (uint32_t j = 0; j<frame.FEs[i].vs_link.size(); j++) ff << frame.FEs[i].vs_link[j] << " ";
		ff << endl;
	}

	int polygon_num = 0;
	for (int i = 0; i<frame.FFs.size(); i++)
		polygon_num += frame.FFs[i].ffs_net.size();

	ff << "POLYGONS " << polygon_num << " " << polygon_num * 5 << endl;
	for (int i = 0; i<frame.FFs.size(); i++)
	{
		for (int j = 0; j<frame.FFs[i].ffs_net.size(); j++)
		{
			ff << "4 " << mesh.Fs[frame.FFs[i].ffs_net[j]].vs[0] << " ";
			ff << mesh.Fs[frame.FFs[i].ffs_net[j]].vs[1] << " " << mesh.Fs[frame.FFs[i].ffs_net[j]].vs[2] << " ";
			ff << mesh.Fs[frame.FFs[i].ffs_net[j]].vs[3] << endl;
		}
	}



	ff << "POINT_DATA " << mesh.Vs.size() << std::endl;
	ff << "SCALARS V_Scalars int" << std::endl;
	ff << "LOOKUP_TABLE V_Table" << std::endl;
	std::vector<short> V_tag(mesh.Vs.size(), 0);
	for (uint32_t i = 0; i < mesh.Vs.size(); i++) if (mesh.Vs[i].boundary) V_tag[i] = 1;
	for (uint32_t j = 0; j<frame.FVs.size(); j++) V_tag[frame.FVs[j].hid] = 3;
	for (uint32_t j = 0; j<frame.FEs.size(); j++) for (uint32_t k = 0; k<frame.FEs[j].vs_link.size(); k++)
		V_tag[frame.FEs[j].vs_link[k]] = 3;

	for (uint32_t i = 0; i<V_tag.size(); i++) ff << V_tag[i] << std::endl;

	ff.close();
}
void h_io::write_Sheet_VTK(Sheet &s, Mesh &mesh, Frame &frame, char *path) {
	vector<bool> H_flag(mesh.Hs.size(), false);
	for (auto cid : s.cs) for (auto hid : frame.FHs[cid].hs_net) H_flag[hid] = true;
	vector<uint32_t> Hs;
	for (uint32_t i = 0; i < H_flag.size(); i++) if (H_flag[i]) Hs.push_back(i);

	std::fstream f(path, std::ios::out);

	f << "# vtk DataFile Version 2.0" << std::endl << "hex mesh vtk data" << std::endl;
	f << "ASCII" << std::endl;
	f << "DATASET UNSTRUCTURED_GRID" << std::endl;


	f << "POINTS " << mesh.V.cols() << " double" << std::endl;
	for (uint32_t i = 0; i<mesh.V.cols(); i++)
		f << mesh.V(0, i) << " " << mesh.V(1, i) << " " << mesh.V(2, i) << std::endl;

	f << "CELLS " << Hs.size() << " " << Hs.size() * 9 << std::endl;
	for (auto h:Hs){
		f << " " << 8 << " " << mesh.Hs[h].vs[0] << " " << mesh.Hs[h].vs[1] << " " << mesh.Hs[h].vs[2] << " " << mesh.Hs[h].vs[3];
		f << " " << mesh.Hs[h].vs[4] << " " << mesh.Hs[h].vs[5] << " " << mesh.Hs[h].vs[6] << " " << mesh.Hs[h].vs[7] << std::endl;
	}

	f << "CELL_TYPES " << Hs.size() << std::endl;
	for (int i = 0; i < Hs.size(); i++) f << 12 << std::endl;

	f << "POINT_DATA " << mesh.V.cols() << std::endl;
	f << "SCALARS fixed int" << std::endl;
	f << "LOOKUP_TABLE default" << std::endl;
	for (int i = 0; i<mesh.Vs.size(); i++) if (mesh.Vs[i].boundary) f << "1" << std::endl; else f << "0" << std::endl;

	f.close();
}
void h_io::write_Vs_Groups_VTK(Mesh &mesh, vector<vector<uint32_t>> &Vs_Group, char *path, bool twoside){
	std::fstream f(path, std::ios::out);

	f << "# vtk DataFile Version 2.0" << std::endl << "hex mesh vtk data" << std::endl;
	f << "ASCII" << std::endl;
	f << "DATASET UNSTRUCTURED_GRID" << std::endl;
	
	VectorXi b;
	MatrixXd bc;

	uint32_t constraint_Num = 0;
	for (auto vs : Vs_Group) constraint_Num += vs.size();
	b.resize(constraint_Num);
	bc.resize(constraint_Num, 3); bc.setZero();
	constraint_Num = 0;
	for (auto vs : Vs_Group) {
		Vector3d v; v.setZero();
		for (auto vid : vs) v = v + mesh.V.col(vid);
		v = v / vs.size();
		for (auto vid : vs) {
			b[constraint_Num] = vid;
			//bc.row(constraint_Num) = v;
			bc.row(constraint_Num) = mesh.V.col(vid);
			constraint_Num++;
		}
	}


	f << "POINTS " << constraint_Num << " double" << std::endl;
	for (uint32_t i = 0; i<bc.rows(); i++)
		f << bc(i, 0) << " " << bc(i, 1) << " " << bc(i, 2) << std::endl;

	f << "POINT_DATA " << bc.rows() << std::endl;
	f << "SCALARS fixed int" << std::endl;
	f << "LOOKUP_TABLE default" << std::endl;
	for (uint32_t i =0;i<Vs_Group.size();i++)
		for (auto vid : Vs_Group[i]) f << i << std::endl;

	f.close();
}

void h_io::write_Points_VTK(MatrixXd &V, char *path) {
	std::fstream f(path, std::ios::out);

	f << "# vtk DataFile Version 2.0" << std::endl << "hex mesh vtk data" << std::endl;
	f << "ASCII" << std::endl;
	f << "DATASET UNSTRUCTURED_GRID" << std::endl;

	f << "POINTS " << V.cols() << " double" << std::endl;
	for (uint32_t i = 0; i<V.cols(); i++)
		f << V(0, i) << " " << V(1, i) << " " << V(2, i) << std::endl;

	f << "POINT_DATA " << V.cols() << std::endl;
	f << "SCALARS fixed int" << std::endl;
	f << "LOOKUP_TABLE default" << std::endl;
	for (uint32_t i = 0; i<V.cols(); i++)
		f << i << std::endl;

	f.close();
}
void h_io::write_Chord_VTK(CHord &c, Mesh &mesh, Frame &frame, char *path) {
	vector<bool> H_flag(mesh.Hs.size(), false);
	for (auto cid : c.cs) for (auto hid : frame.FHs[cid].hs_net) H_flag[hid] = true;
	vector<uint32_t> Hs;
	for (uint32_t i = 0; i < H_flag.size(); i++) if (H_flag[i]) Hs.push_back(i);

	std::fstream f(path, std::ios::out);

	f << "# vtk DataFile Version 2.0" << std::endl << "hex mesh vtk data" << std::endl;
	f << "ASCII" << std::endl;
	f << "DATASET UNSTRUCTURED_GRID" << std::endl;


	f << "POINTS " << mesh.V.cols() << " double" << std::endl;
	for (uint32_t i = 0; i<mesh.V.cols(); i++)
		f << mesh.V(0, i) << " " << mesh.V(1, i) << " " << mesh.V(2, i) << std::endl;

	f << "CELLS " << Hs.size() << " " << Hs.size() * 9 << std::endl;
	for (auto h : Hs) {
		f << " " << 8 << " " << mesh.Hs[h].vs[0] << " " << mesh.Hs[h].vs[1] << " " << mesh.Hs[h].vs[2] << " " << mesh.Hs[h].vs[3];
		f << " " << mesh.Hs[h].vs[4] << " " << mesh.Hs[h].vs[5] << " " << mesh.Hs[h].vs[6] << " " << mesh.Hs[h].vs[7] << std::endl;
	}

	f << "CELL_TYPES " << Hs.size() << std::endl;
	for (int i = 0; i < Hs.size(); i++) f << 12 << std::endl;

	f << "POINT_DATA " << mesh.V.cols() << std::endl;
	f << "SCALARS fixed int" << std::endl;
	f << "LOOKUP_TABLE default" << std::endl;
	for (int i = 0; i<mesh.Vs.size(); i++) if (mesh.Vs[i].boundary) f << "1" << std::endl; else f << "0" << std::endl;

	f.close();
}

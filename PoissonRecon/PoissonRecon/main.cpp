#include "PoissonRecon.h"
#include "plyio.h"

//#include "stlio.h"

void main(int argc, char *argv[])
{
	PoissonReconstruction fsp;

	orth::MeshModel mm,mm2;
	tinyply::plyio io;
	//tinystl::stlio io2;
	io.read_ply_file(argv[1], mm);
	vector<orth::MeshModel> mms;
	mms.push_back(mm);
	//for (int i = 0; i < 10000; i++)
	{
		//std::cout << "------------------------------------ " << i << " ---------------------------------------" << std::endl;
		fsp.run(mm,  9,0.5);

	}
	
	//io2.write_stl_file(argv[2], mm, true);
	io.write_ply_file(argv[2], mm, true);
}
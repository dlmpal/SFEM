// Convert a gmsh file to native sfem mesh format
// .msh2 files are suported
// The program is executed as follows:
//   gmshToSfem $gmsh_path $mesh_dir

#include "sfem.h"

int main(int argc, char **argv)
{
    sfem::initialize(&argc, &argv, "GMSH_TO_SFEM");

    std::string gmsh_path = argv[1];
    std::string mesh_dir = argv[2];

    auto mesh = sfem::io::read_gmsh(gmsh_path);
    sfem::io::write_mesh(mesh_dir, mesh);

    sfem::finalize();
    return 0;
}

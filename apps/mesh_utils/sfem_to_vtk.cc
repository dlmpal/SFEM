// Convert a mesh in native format to VTK.
// The program is executed, inside the directory where the field values are located, as follows:
//   sfemToVTK $mesh_dir $time_steps $field_name_1 $field_name_2 ...
// It is expected that for each field and time step there exists a separate file,
// e.g. ($field_name_1)_0 for the first field and first time step.
// The program produces as many VTK files as there are timesteps
// The VTK files are placed in the current directory, and named as sfem_0.vtk, sfem_1.vtk, etc

#include "sfem.h"

int main(int argc, char **argv)
{
    sfem::initialize(&argc, &argv, "SFEM_TO_VTK");

    std::string mesh_dir = argv[1];
    auto mesh = sfem::io::read_mesh(mesh_dir);
    int time_steps = std::atoi(argv[2]);
    int n_fields = argc - 3;
    std::vector<sfem::mesh::Field> fields;

    // Initialize all fields
    for (auto i = 0; i < n_fields; i++)
    {
        std::string name = argv[i + 3];
        std::string path = name + "_" + std::to_string(0);
        std::ifstream file(path);
        if (!file.is_open())
        {
            sfem::error::invalid_filename_error(path, __FILE__, __LINE__);
        }
        int n_vars;
        file >> n_vars;
        fields.push_back(sfem::mesh::Field(name, n_vars, mesh));
    }

    // Create a .vtk file for each timestep
    for (auto time = 0; time < time_steps; time++)
    {
        // Update field values every timestep
        for (auto &field : fields)
        {
            std::string path = field.name() + "_" + std::to_string(time);
            sfem::io::read_field_values(path, field);
        }

        std::string path = "sfem_" + std::to_string(time) + ".vtk";
        sfem::io::write_vtk(path, mesh, fields);
    }

    sfem::finalize();
    return 0;
}

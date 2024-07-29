#include "field.h"
#include "../common/error.h"

namespace sfem::io
{
    //=============================================================================
    void read_field_values(const std::string &path, mesh::Field &field)
    {
        std::ifstream file(path);
        if (!file.is_open())
        {
            error::invalid_filename_error(path, __FILE__, __LINE__);
        }

        int n_vars;
        file >> n_vars;
        if (n_vars != field.n_vars())
        {
            error::invalid_size_error(field.n_vars(), n_vars, __FILE__, __LINE__);
        }

        auto im = field.mesh().node_im();

        std::vector<Scalar> values(field.n_dof_local());
        std::vector<Scalar> value(field.n_vars());
        for (int i = 0; i < im.n_global(); i++)
        {
            for (int j = 0; j < n_vars; j++)
            {
                file >> value[j];
            }

            // Only add the values for the local DOF
            int local_idx = im.global_to_local(i);
            if (local_idx >= 0)
            {
                for (int j = 0; j < n_vars; j++)
                {
                    values[local_idx * n_vars + j] = value[j];
                }
            }
        }
        field.set_values(values);
    }
    //=============================================================================
    void write_field_values(const std::string &path, const mesh::Field &field, bool assemble_global)
    {
        std::string _path = path;
        std::vector<Scalar> values;
        if (assemble_global)
        {
            values = mesh::gather_field_values(field);
        }
        else
        {
            values = field.values();
            if (Logger::instance().n_procs() > 1)
            {
                _path += "_proc_" + std::to_string(Logger::instance().proc_rank());
            }
        }
        if (Logger::instance().proc_rank() == SFEM_ROOT || assemble_global == false)
        {
            std::ofstream file(_path);
            if (!file.is_open())
            {
                error::invalid_filename_error(_path, __FILE__, __LINE__);
            }
            file << field.n_vars() << "\n";
            for (std::size_t i = 0; i < values.size(); i++)
            {
                file << values[i] << "\n";
            }
        }
    }
}
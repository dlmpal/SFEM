#include "field.h"
#include "../common/error.h"
#include <mpi.h>

namespace sfem::mesh
{
    //=============================================================================
    Field::Field(const std::string &name,
                 int n_vars,
                 mesh::Mesh &mesh,
                 const std::vector<std::string> &component_names)
        : name_(name),
          n_vars_(n_vars),
          mesh_(mesh),
          dof_im_(mesh.node_im().renumber())
    {
        if (n_vars <= 0)
        {
            Logger::instance().error("Field: " + name + " has an invalid number of variables: " + std::to_string(n_vars) + "\n",
                                     __FILE__, __LINE__);
        }

        if (component_names.size() != static_cast<std::size_t>(n_vars))
        {
            if (component_names.size() == 0)
            {
                for (int i = 0; i < n_vars; i++)
                {
                    comp_names_.push_back(name + "_" + std::to_string(i));
                }
            }
            else
            {
                error::invalid_size_error(n_vars, component_names.size(), __FILE__, __LINE__);
            }
        }
        else
        {
            for (int i = 0; i < n_vars; i++)
            {
                comp_names_.push_back(component_names[i]);
            }
        }

        values_.resize(n_vars * mesh.n_nodes_local());
    }
    //=============================================================================
    std::string Field::name() const
    {
        return name_;
    }
    //=============================================================================
    int Field::n_vars() const
    {
        return n_vars_;
    }
    //=============================================================================
    mesh::Mesh &Field::mesh() const
    {
        return mesh_;
    }
    //=============================================================================
    std::vector<std::string> Field::comp_names() const
    {
        return comp_names_;
    }
    //=============================================================================
    common::IndexMap Field::dof_im() const
    {
        return dof_im_;
    }
    //=============================================================================
    int Field::n_dof_local() const
    {
        return n_vars_ * dof_im_.n_local();
    }
    //=============================================================================
    int Field::n_dof_owned() const
    {
        return n_vars_ * dof_im_.n_owned();
    }
    //=============================================================================
    int Field::n_dof_ghost() const
    {
        return n_vars_ * dof_im_.n_ghost();
    }
    //=============================================================================
    int Field::n_dof_global() const
    {
        return n_vars_ * dof_im_.n_global();
    }
    //=============================================================================
    std::vector<int> Field::map_node_dof(const std::vector<int> &nodes) const
    {
        std::vector<int> dof(nodes.size() * n_vars_);
        for (std::size_t i = 0; i < nodes.size(); i++)
        {
            for (int j = 0; j < n_vars_; j++)
            {
                dof[i * n_vars_ + j] = nodes[i] * n_vars_ + j;
            }
        }
        return dof;
    }
    //=============================================================================
    std::vector<int> Field::get_owned_dof() const
    {
        auto owned_nodes = dof_im_.get_owned_idxs();
        return map_node_dof(owned_nodes);
    }
    //=============================================================================
    std::vector<int> Field::get_ghost_dof() const
    {
        auto ghost_nodes = dof_im_.get_ghost_idxs();
        return map_node_dof(ghost_nodes);
    }
    //=============================================================================
    std::vector<int> Field::get_cell_dof(const mesh::Cell &cell) const
    {
        auto cell_nodes = mesh_.get_cell_nodes(cell);
        for (std::size_t i = 0; i < cell_nodes.size(); i++)
        {
            cell_nodes[i] = dof_im_.local_to_global(cell_nodes[i]);
        }
        return map_node_dof(cell_nodes);
    }
    //=============================================================================
    std::vector<Scalar> Field::get_cell_values(const mesh::Cell &cell) const
    {
        auto cell_nodes = mesh_.get_cell_nodes(cell);
        std::vector<Scalar> cell_values(cell.n_nodes() * n_vars_);
        for (int i = 0; i < cell.n_nodes(); i++)
        {
            for (int j = 0; j < n_vars_; j++)
            {
                cell_values[i * n_vars_ + j] = values_[cell_nodes[i] * n_vars_ + j];
            }
        }
        return cell_values;
    }
    //=============================================================================
    void Field::add_fixed_dof(const std::string &region_name, int var, Scalar value)
    {
        auto region_nodes = mesh_.get_region_nodes(region_name);
        for (auto node : region_nodes)
        {
            node = dof_im_.local_to_global(node);
            fixed_dof_[node * n_vars_ + var] = value;
        }
    }
    //=============================================================================
    std::vector<int> Field::get_fixed_dof() const
    {
        std::vector<int> fdof;
        fdof.reserve(fixed_dof_.size());
        for (auto kv : fixed_dof_)
        {
            fdof.push_back(kv.first);
        }
        return fdof;
    }
    //=============================================================================
    std::vector<Scalar> Field::get_fixed_dof_values() const
    {
        std::vector<Scalar> fdof_values;
        fdof_values.reserve(fixed_dof_.size());
        for (auto kv : fixed_dof_)
        {
            fdof_values.push_back(kv.second);
        }
        return fdof_values;
    }
    //=============================================================================
    void Field::clear_fixed_dof()
    {
        fixed_dof_.clear();
    }
    //=============================================================================
    void Field::set_all(const std::vector<Scalar> &value)
    {
        if (value.size() != static_cast<std::size_t>(n_vars_))
        {
            error::invalid_size_error(n_vars_, value.size(), __FILE__, __LINE__);
        }

        for (int i = 0; i < mesh_.n_nodes_local(); i++)
        {
            for (int j = 0; j < n_vars_; j++)
            {
                values_[i * n_vars_ + j] = value[j];
            }
        }
    }
    //=============================================================================
    void Field::set_values(const std::vector<Scalar> &values)
    {
        if (values.size() != static_cast<std::size_t>(n_dof_local()))
        {
            error::invalid_size_error(n_dof_local(), values.size(), __FILE__, __LINE__);
        }
        values_ = values;
    }
    //=============================================================================
    const std::vector<Scalar> &Field::values() const
    {
        return values_;
    }
    //=============================================================================
    std::vector<Scalar> gather_field_values(const Field &field)
    {
        int n_procs = Logger::instance().n_procs();
        int proc_rank = Logger::instance().proc_rank();

        // For serial execution, return early
        if (n_procs == 1)
        {
            return field.values();
        }

        int n_dof_owned = field.n_dof_owned();
        int n_dof_global = field.n_dof_global();

        auto owned_nodes = field.mesh().node_im().get_owned_idxs();
        auto send_buffer_dof = field.map_node_dof(owned_nodes);
        const auto &send_buffer_values = field.values();

        std::vector<int> recv_buffer_dof;
        std::vector<Scalar> recv_buffer_values;
        std::vector<int> recv_counts(n_procs, 0);
        std::vector<int> recv_displs(n_procs, 0);

        if (proc_rank == SFEM_ROOT)
        {
            recv_buffer_dof.resize(n_dof_global);
            recv_buffer_values.resize(n_dof_global);
        }

        // Compute receive counts and displacements
        MPI_Gather(&n_dof_owned, 1, MPI_INT,
                   recv_counts.data(), 1, MPI_INT,
                   SFEM_ROOT, SFEM_COMM_WORLD);
        for (int i = 0; i < n_procs; i++)
        {
            if (i > 0)
            {
                recv_displs[i] = recv_displs[i - 1] + recv_counts[i - 1];
            }
        }

        // Gather all the DoF and their values at the root process
        MPI_Gatherv(send_buffer_dof.data(), n_dof_owned, MPI_INT,
                    recv_buffer_dof.data(), recv_counts.data(), recv_displs.data(), MPI_INT,
                    SFEM_ROOT, SFEM_COMM_WORLD);
        MPI_Gatherv(send_buffer_values.data(), n_dof_owned, SFEM_MPI_FLOAT,
                    recv_buffer_values.data(), recv_counts.data(), recv_displs.data(), SFEM_MPI_FLOAT,
                    SFEM_ROOT, SFEM_COMM_WORLD);

        // Re-order the values
        std::vector<Scalar> field_values;
        if (proc_rank == SFEM_ROOT)
        {
            field_values.resize(n_dof_global);
            for (int i = 0; i < n_dof_global; i++)
            {
                field_values[recv_buffer_dof[i]] = recv_buffer_values[i];
            }
        }

        return field_values;
    }
}
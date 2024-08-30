#pragma once

#include "mesh.h"

namespace sfem::mesh
{
    class Field
    {
    public:
        /// @brief Create a Field
        /// @param name Name with which the field is referenced
        /// @param n_vars Number of variables per node
        /// @param mesh Mesh on which the field is defined
        /// @param component_names The names of the field components. By default they are set to {$name}_i
        Field(const std::string &name,
              int n_vars,
              mesh::Mesh &mesh,
              const std::vector<std::string> &component_names = {});

        /// @brief Copy constructor
        Field(const Field &) = default;

        // Copy assignment (deleted)
        Field &operator=(const Field &) = delete;

        /// @brief Move constructor
        Field(Field &&) = default;

        /// @brief Move assignment
        Field &operator=(Field &&mesh) = default;

        /// @brief Destructor
        ~Field() = default;

        /// @brief Get the name
        std::string name() const;

        /// @brief Get the number of variables
        int n_vars() const;

        /// @brief Get a pointer to the Mesh
        mesh::Mesh &mesh() const;

        /// @brief Get the component names
        std::vector<std::string> comp_names() const;

        /// @brief Get the DoF index map
        common::IndexMap dof_im() const;

        /// @brief Get the number of owned DoF for this process
        int n_dof_owned() const;

        /// @brief Get the number of ghost DoF for this process
        int n_dof_ghost() const;

        /// @brief Get the number of local DoF
        int n_dof_local() const;

        /// @brief Get the global number of DoF
        int n_dof_global() const;

        /// @brief  Get the DoF corresponding to the given nodes
        std::vector<int> map_node_dof(const std::vector<int> &nodes) const;

        /// @brief Get the owned DoF for this process
        std::vector<int> get_owned_dof() const;

        /// @brief Get the ghsot DoF for this process
        std::vector<int> get_ghost_dof() const;

        /// @brief Get the DoF local to this process (owned + ghost)
        std::vector<int> get_local_dof() const;

        /// @brief Get the DoF belonging to a cell
        /// @note The DoF are returned in global indexing
        std::vector<int> get_cell_dof(const mesh::Cell &cell) const;

        /// @brief Get the values belonging to a cell
        std::vector<Scalar> get_cell_values(const mesh::Cell &cell) const;

        /// @brief Assign fixed values to desired DoF
        void add_fixed_dof(const std::string &region_name, int var, Scalar value);

        /// @brief Get the fixed DoF
        std::vector<int> get_fixed_dof() const;

        /// brief Get the values of the fixed DoF
        std::vector<Scalar> get_fixed_dof_values() const;

        /// @brief Clear all existing fixed DoF
        void clear_fixed_dof();

        /// @brief Set local DoF values uniformly
        /// @param value A vector of size n_vars
        void set_all(const std::vector<Scalar> &value);

        /// @brief Set local DoF values
        void set_values(const std::vector<Scalar> &values);

        /// @brief Get local DoF values
        const std::vector<Scalar> &values() const;

    private:
        /// @brief Field name
        std::string name_;

        /// @brief Number of variables per node
        int n_vars_;

        /// @brief Mesh on which the Field operates
        mesh::Mesh &mesh_;

        /// @brief Component names
        std::vector<std::string> comp_names_;

        /// @brief DoF index map
        common::IndexMap dof_im_;

        /// @brief Constrained DoF and their values
        std::unordered_map<int, Scalar> fixed_dof_;

        /// @brief Values corresponding to the local DoF
        std::vector<Scalar> values_;
    };

    /// @brief Assemble all Field values to the root process
    std::vector<Scalar> gather_field_values(const Field &field);
}
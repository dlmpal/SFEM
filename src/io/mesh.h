#pragma once

#include "../mesh/sfem_mesh.h"

namespace sfem::io
{
    /// @brief
    /// @param path
    /// @param distributed
    /// @param cell_im
    /// @param node_im
    /// @return
    std::pair<std::vector<mesh::Cell>, mesh::Connectivity>
    read_cells(const std::string &path, bool distributed, const common::IndexMap &cell_im, const common::IndexMap &node_im);

    /// @brief
    /// @param path
    /// @param distributed
    /// @param node_im
    /// @return
    std::vector<Scalar> read_xpts(const std::string &path, bool distributed, const common::IndexMap &node_im);

    /// @brief
    /// @param path
    /// @return
    std::vector<mesh::Region> read_regions(const std::string &path);

    /// @brief Read a mesh from file
    /// @note For distributed meshes, the mesh is first partitioned
    /// @param dir Directory in which the mesh files are located
    /// @param partitioner_type Type of MeshPartitioner to be used. Defaults to "METIS"
    /// @return The portion of the mesh corresponding to this process
    mesh::Mesh read_mesh(const std::string &dir, const std::string &partitioner_type = "METIS");

    /// @brief
    /// @param path
    /// @param mesh
    void write_cells(const std::string &path, const mesh::Mesh &mesh);

    /// @brief
    /// @param path
    /// @param mesh
    void write_xts(const std::string &path, const mesh::Mesh &mesh);

    /// @brief
    /// @param path
    /// @param mesh
    void write_regions(const std::string &path, const mesh::Mesh &mesh);

    /// @brief Write a Mesh to file
    /// @note Not designed for parallel execution
    void write_mesh(const std::string &dir, const mesh::Mesh &mesh);
}
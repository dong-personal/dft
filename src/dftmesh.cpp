

#include "dftmesh.h"

#include <algorithm>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <utility>

namespace dft
{

void DFTMesh::set_structure(const Structure *structure)
{
    m_structure_owner.reset();
    m_structure = structure;
}

void DFTMesh::set_structure(const Structure &structure)
{
    m_structure_owner.reset();
    m_structure = &structure;
}

void DFTMesh::set_structure(std::shared_ptr<const Structure> structure)
{
    m_structure_owner = std::move(structure);
    m_structure = m_structure_owner.get();
}

mfem::Mesh &DFTMesh::mesh()
{
    ensure_mesh_();
    return *m_mesh;
}

const mfem::Mesh &DFTMesh::mesh() const
{
    ensure_mesh_();
    return *m_mesh;
}

void DFTMesh::save(const std::string &file) const
{
    ensure_mesh_();

    std::ofstream output_stream(file);
    if (!output_stream)
    {
        throw std::runtime_error("DFTMesh::save: cannot open " + file);
    }

    m_mesh->Print(output_stream);
}

void DFTMesh::save_vtu(const std::string &file) const
{
    ensure_mesh_();

    std::ofstream output_stream(file);
    if (!output_stream)
    {
        throw std::runtime_error("DFTMesh::save_vtu: cannot open " + file);
    }

    output_stream << "<?xml version=\"1.0\"?>\n";
    output_stream << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    output_stream << "<UnstructuredGrid>\n";
    m_mesh->PrintVTU(output_stream);
    output_stream << "</Piece>\n";
    output_stream << "</UnstructuredGrid>\n";
    output_stream << "</VTKFile>\n";
}

void DFTMesh::init_periodic_cell_from_lattice(const Structure::LatticeVectors &lattice_vectors, int nx, int ny, int nz,
                                              mfem::Element::Type element_type, bool use_sfc_ordering)
{
    if (nx < 3 || ny < 3 || nz < 3)
    {
        throw std::runtime_error("Periodic mesh needs >= 2 interior vertices per direction. Use nx, ny, nz >= 3.");
    }

    mfem::Mesh base_mesh = mfem::Mesh::MakeCartesian3D(nx, ny, nz, element_type, 1.0, 1.0, 1.0, use_sfc_ordering);

    for (int vertex_index = 0; vertex_index < base_mesh.GetNV(); ++vertex_index)
    {
        double *vertex = base_mesh.GetVertex(vertex_index);
        const double u = vertex[0];
        const double v = vertex[1];
        const double w = vertex[2];

        vertex[0] = u * lattice_vectors[0, 0] + v * lattice_vectors[1, 0] + w * lattice_vectors[2, 0];
        vertex[1] = u * lattice_vectors[0, 1] + v * lattice_vectors[1, 1] + w * lattice_vectors[2, 1];
        vertex[2] = u * lattice_vectors[0, 2] + v * lattice_vectors[1, 2] + w * lattice_vectors[2, 2];
    }

    m_lattice_vectors = lattice_vectors;
    m_translations = make_translations_from_lattice_(lattice_vectors);
    m_base_mesh = std::make_unique<mfem::Mesh>(std::move(base_mesh));

    rebuild_periodic_from_base_();
}

void DFTMesh::init_periodic_cell_from_structure(int nx, int ny, int nz, mfem::Element::Type element_type,
                                                bool use_sfc_ordering)
{
    ensure_structure_();
    init_periodic_cell_from_lattice(m_structure->lattice(), nx, ny, nz, element_type, use_sfc_ordering);
}

void DFTMesh::refine_near_atoms(double radius, int levels, double shrink)
{
    ensure_base_mesh_();
    ensure_structure_();
    if (radius <= 0.0 || levels <= 0)
    {
        return;
    }

    m_base_mesh->EnsureNCMesh();
    const std::vector<mfem::Vector> atom_positions = atom_positions_cartesian_();

    for (int level = 0; level < levels; ++level)
    {
        const double radius_squared = radius * radius;
        mfem::Array<int> refine_list;
        refine_list.Reserve(m_base_mesh->GetNE() / 8 + 1);

        for (int element_index = 0; element_index < m_base_mesh->GetNE(); ++element_index)
        {
            const mfem::Vector center = element_center_from_transform_(*m_base_mesh, element_index);

            bool is_near_atom = false;
            for (const auto &atom_position : atom_positions)
            {
                if (periodic_dist2_(center, atom_position) < radius_squared)
                {
                    is_near_atom = true;
                    break;
                }
            }

            if (is_near_atom)
            {
                refine_list.Append(element_index);
            }
        }

        if (refine_list.Size() == 0)
        {
            break;
        }

        m_base_mesh->GeneralRefinement(refine_list);

        radius *= shrink;
        if (radius <= 0.0)
        {
            break;
        }
    }

    rebuild_periodic_from_base_();
}

void DFTMesh::ensure_mesh_() const
{
    if (!m_mesh)
    {
        throw std::runtime_error("DFTMesh: mesh not initialized");
    }
}

void DFTMesh::ensure_base_mesh_() const
{
    if (!m_base_mesh)
    {
        throw std::runtime_error("DFTMesh: base mesh not initialized");
    }
}

void DFTMesh::ensure_structure_() const
{
    if (!m_structure)
    {
        throw std::runtime_error("DFTMesh: structure not set");
    }
}

std::vector<mfem::Vector> DFTMesh::make_translations_from_lattice_(const Structure::LatticeVectors &lattice_vectors)
{
    std::vector<mfem::Vector> translations;
    translations.reserve(kSpatialDimension);

    for (int row = 0; row < kSpatialDimension; ++row)
    {
        mfem::Vector translation(kSpatialDimension);
        for (int dim = 0; dim < kSpatialDimension; ++dim)
        {
            translation[dim] = lattice_vectors[row, dim];
        }
        translations.push_back(translation);
    }

    return translations;
}

void DFTMesh::rebuild_periodic_from_base_()
{
    ensure_base_mesh_();
    if (m_translations.empty())
    {
        throw std::runtime_error("DFTMesh: translations not set");
    }

    const std::vector<int> vertex_mapping = m_base_mesh->CreatePeriodicVertexMapping(m_translations);
    mfem::Mesh periodic_mesh = mfem::Mesh::MakePeriodic(*m_base_mesh, vertex_mapping);
    m_mesh = std::make_unique<mfem::Mesh>(std::move(periodic_mesh));
}

mfem::Vector DFTMesh::fractional_to_cartesian_(const Structure::LatticeVectors &lattice_vectors,
                                               const Structure::AtomicPosition &fractional_position)
{
    mfem::Vector cartesian_position(kSpatialDimension);
    cartesian_position[0] = fractional_position[0] * lattice_vectors[0, 0] +
                            fractional_position[1] * lattice_vectors[1, 0] +
                            fractional_position[2] * lattice_vectors[2, 0];
    cartesian_position[1] = fractional_position[0] * lattice_vectors[0, 1] +
                            fractional_position[1] * lattice_vectors[1, 1] +
                            fractional_position[2] * lattice_vectors[2, 1];
    cartesian_position[2] = fractional_position[0] * lattice_vectors[0, 2] +
                            fractional_position[1] * lattice_vectors[1, 2] +
                            fractional_position[2] * lattice_vectors[2, 2];
    return cartesian_position;
}

std::vector<mfem::Vector> DFTMesh::atom_positions_cartesian_() const
{
    std::vector<mfem::Vector> atom_positions;
    atom_positions.reserve(m_structure->num_atoms());

    const auto &lattice_vectors = m_structure->lattice();
    for (const auto &atom : m_structure->atoms())
    {
        for (std::size_t position_index = 0; position_index < atom.num_positions(); ++position_index)
        {
            const auto atomic_position = atom.position(position_index);
            if (m_structure->coord_type() == Structure::CoordType::Fractional)
            {
                atom_positions.push_back(fractional_to_cartesian_(lattice_vectors, atomic_position));
                continue;
            }

            mfem::Vector cartesian_position(kSpatialDimension);
            for (int dim = 0; dim < kSpatialDimension; ++dim)
            {
                cartesian_position[dim] = atomic_position[dim];
            }
            atom_positions.push_back(cartesian_position);
        }
    }

    return atom_positions;
}

double DFTMesh::dist2_(const mfem::Vector &lhs, const mfem::Vector &rhs)
{
    const double dx = lhs[0] - rhs[0];
    const double dy = lhs[1] - rhs[1];
    const double dz = lhs[2] - rhs[2];
    return dx * dx + dy * dy + dz * dz;
}

double DFTMesh::periodic_dist2_(const mfem::Vector &point, const mfem::Vector &atom_position) const
{
    const auto &lattice_vectors = m_structure->lattice();

    mfem::Vector a1(kSpatialDimension);
    mfem::Vector a2(kSpatialDimension);
    mfem::Vector a3(kSpatialDimension);
    for (int dim = 0; dim < kSpatialDimension; ++dim)
    {
        a1[dim] = lattice_vectors[0, dim];
        a2[dim] = lattice_vectors[1, dim];
        a3[dim] = lattice_vectors[2, dim];
    }

    double min_distance_squared = std::numeric_limits<double>::max();
    for (int i = -1; i <= 1; ++i)
    {
        for (int j = -1; j <= 1; ++j)
        {
            for (int k = -1; k <= 1; ++k)
            {
                mfem::Vector translated_atom_position(atom_position);
                translated_atom_position[0] += i * a1[0] + j * a2[0] + k * a3[0];
                translated_atom_position[1] += i * a1[1] + j * a2[1] + k * a3[1];
                translated_atom_position[2] += i * a1[2] + j * a2[2] + k * a3[2];

                min_distance_squared = std::min(min_distance_squared, dist2_(point, translated_atom_position));
            }
        }
    }

    return min_distance_squared;
}

mfem::Vector DFTMesh::element_center_from_transform_(mfem::Mesh &mesh, int element_index)
{
    mfem::ElementTransformation *transformation = mesh.GetElementTransformation(element_index);
    mfem::IntegrationPoint integration_point;
    integration_point.Set3(0.5, 0.5, 0.5);

    mfem::Vector center(kSpatialDimension);
    transformation->Transform(integration_point, center);
    return center;
}

} // namespace dft

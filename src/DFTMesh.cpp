#include "DFTMesh.h"

#include <limits>
#include <utility>

namespace dft
{

void DFTMesh::set_structure(const Structure *s)
{
    structure_owner_.reset();
    structure_ = s;
}

void DFTMesh::set_structure(const Structure &s)
{
    structure_owner_.reset();
    structure_ = &s;
}

void DFTMesh::set_structure(std::shared_ptr<const Structure> s)
{
    structure_owner_ = std::move(s);
    structure_ = structure_owner_.get();
}

mfem::Mesh &DFTMesh::mesh()
{
    ensure_mesh();
    return *mesh_;
}

const mfem::Mesh &DFTMesh::mesh() const
{
    ensure_mesh();
    return *mesh_;
}

void DFTMesh::save(const std::string &file) const
{
    ensure_mesh();
    std::ofstream ofs(file);
    if (!ofs)
    {
        throw std::runtime_error("DFTMesh::save: cannot open " + file);
    }
    mesh_->Print(ofs);
}

void DFTMesh::saveVTU(const std::string &file) const
{
    ensure_mesh();

    std::ofstream ofs(file);
    if (!ofs)
    {
        throw std::runtime_error("DFTMesh::saveGmsh: cannot open " + file);
    }

    ofs << "<?xml version=\"1.0\"?>\n";
    ofs << "<VTKFile type=\"UnstructuredGrid\" "
           "version=\"0.1\" byte_order=\"LittleEndian\">\n";
    ofs << "<UnstructuredGrid>\n";
    mesh_->PrintVTU(ofs);
    ofs << "</Piece>\n";
    ofs << "</UnstructuredGrid>\n";
    ofs << "</VTKFile>\n";
}

void DFTMesh::init_periodic_cell_from_lattice(const Structure::Mat3 &lat, int nx, int ny, int nz,
                                              mfem::Element::Type type, bool sfc_ordering)
{
    if (nx < 3 || ny < 3 || nz < 3)
    {
        throw std::runtime_error("Periodic mesh needs >=2 interior vertices per direction. "
                                 "Use nx,ny,nz >= 3 (elements).");
    }

    mfem::Mesh m = mfem::Mesh::MakeCartesian3D(nx, ny, nz, type, 1.0, 1.0, 1.0, sfc_ordering);

    for (int i = 0; i < m.GetNV(); ++i)
    {
        double *v = m.GetVertex(i);
        const double u = v[0];
        const double vv = v[1];
        const double w = v[2];

        const double x = u * lat[0][0] + vv * lat[1][0] + w * lat[2][0];
        const double y = u * lat[0][1] + vv * lat[1][1] + w * lat[2][1];
        const double z = u * lat[0][2] + vv * lat[1][2] + w * lat[2][2];

        v[0] = x;
        v[1] = y;
        v[2] = z;
    }

    lattice_ = lat;
    translations_ = make_translations_from_lattice_(lat);
    base_mesh_ = std::make_unique<mfem::Mesh>(std::move(m));

    rebuild_periodic_from_base_();
}

void DFTMesh::init_periodic_cell_from_structure(int nx, int ny, int nz, mfem::Element::Type type, bool sfc_ordering)
{
    ensure_structure();
    init_periodic_cell_from_lattice(structure_->lattice(), nx, ny, nz, type, sfc_ordering);
}

void DFTMesh::refine_near_atoms(double radius, int levels, double shrink)
{
    ensure_base_mesh();
    ensure_structure();
    if (radius <= 0.0 || levels <= 0)
    {
        return;
    }

    base_mesh_->EnsureNCMesh();
    std::vector<mfem::Vector> atom_pos = atom_positions_cartesian_();

    for (int lvl = 0; lvl < levels; ++lvl)
    {
        const double r2 = radius * radius;

        mfem::Array<int> refine_list;
        refine_list.Reserve(base_mesh_->GetNE() / 8 + 1);

        for (int el = 0; el < base_mesh_->GetNE(); ++el)
        {
            mfem::Vector c = element_center_via_transform_(*base_mesh_, el);

            bool hit = false;
            for (const auto &rp : atom_pos)
            {
                if (periodic_dist2_(c, rp) < r2)
                {
                    hit = true;
                    break;
                }
            }
            if (hit)
            {
                refine_list.Append(el);
            }
        }

        if (refine_list.Size() == 0)
        {
            break;
        }

        base_mesh_->GeneralRefinement(refine_list);

        radius *= shrink;
        if (radius <= 0.0)
        {
            break;
        }
    }

    rebuild_periodic_from_base_();
}

void DFTMesh::ensure_mesh() const
{
    if (!mesh_)
    {
        throw std::runtime_error("DFTMesh: mesh not initialized");
    }
}

void DFTMesh::ensure_base_mesh() const
{
    if (!base_mesh_)
    {
        throw std::runtime_error("DFTMesh: base mesh not initialized");
    }
}

void DFTMesh::ensure_structure() const
{
    if (!structure_)
    {
        throw std::runtime_error("DFTMesh: structure not set");
    }
}

std::vector<mfem::Vector> DFTMesh::make_translations_from_lattice_(const Structure::Mat3 &lat)
{
    std::vector<mfem::Vector> translations;
    translations.reserve(3);

    for (int i = 0; i < 3; ++i)
    {
        mfem::Vector t(3);
        t[0] = lat[i][0];
        t[1] = lat[i][1];
        t[2] = lat[i][2];
        translations.push_back(t);
    }
    return translations;
}

void DFTMesh::rebuild_periodic_from_base_()
{
    ensure_base_mesh();
    if (translations_.empty())
    {
        throw std::runtime_error("DFTMesh: translations not set (init first)");
    }

    std::vector<int> v2v = base_mesh_->CreatePeriodicVertexMapping(translations_);
    mfem::Mesh pm = mfem::Mesh::MakePeriodic(*base_mesh_, v2v);
    mesh_ = std::make_unique<mfem::Mesh>(std::move(pm));
}

mfem::Vector DFTMesh::frac_to_cart_(const Structure::Mat3 &lat, const Vec3 &f)
{
    mfem::Vector r(3);
    r[0] = f[0] * lat[0][0] + f[1] * lat[1][0] + f[2] * lat[2][0];
    r[1] = f[0] * lat[0][1] + f[1] * lat[1][1] + f[2] * lat[2][1];
    r[2] = f[0] * lat[0][2] + f[1] * lat[1][2] + f[2] * lat[2][2];
    return r;
}

std::vector<mfem::Vector> DFTMesh::atom_positions_cartesian_() const
{
    std::vector<mfem::Vector> pos;
    pos.reserve(structure_->num_atoms());
    const auto &lat = structure_->lattice();

    for (const auto &a : structure_->atoms())
    {
        const auto &p = a.position();
        if (structure_->coord_type() == Structure::CoordType::Fractional)
        {
            pos.push_back(frac_to_cart_(lat, p));
        }
        else
        {
            mfem::Vector r(3);
            r[0] = p[0];
            r[1] = p[1];
            r[2] = p[2];
            pos.push_back(r);
        }
    }
    return pos;
}

double DFTMesh::dist2_(const mfem::Vector &a, const mfem::Vector &b)
{
    const double dx = a[0] - b[0];
    const double dy = a[1] - b[1];
    const double dz = a[2] - b[2];
    return dx * dx + dy * dy + dz * dz;
}

double DFTMesh::periodic_dist2_(const mfem::Vector &x, const mfem::Vector &a) const
{
    const auto &lat = structure_->lattice();

    mfem::Vector a1(3), a2(3), a3(3);
    for (int d = 0; d < 3; ++d)
    {
        a1[d] = lat[0][d];
        a2[d] = lat[1][d];
        a3[d] = lat[2][d];
    }

    double min_d2 = std::numeric_limits<double>::max();

    for (int i = -1; i <= 1; ++i)
    {
        for (int j = -1; j <= 1; ++j)
        {
            for (int k = -1; k <= 1; ++k)
            {
                const double dx = x[0] - (a[0] + i * a1[0] + j * a2[0] + k * a3[0]);
                const double dy = x[1] - (a[1] + i * a1[1] + j * a2[1] + k * a3[1]);
                const double dz = x[2] - (a[2] + i * a1[2] + j * a2[2] + k * a3[2]);

                const double d2 = dx * dx + dy * dy + dz * dz;
                if (d2 < min_d2)
                {
                    min_d2 = d2;
                }
            }
        }
    }
    return min_d2;
}

mfem::Vector DFTMesh::element_center_via_transform_(mfem::Mesh &mesh, int el_id)
{
    mfem::ElementTransformation *T = mesh.GetElementTransformation(el_id);
    mfem::IntegrationPoint ip;
    ip.Set3(0.5, 0.5, 0.5);

    mfem::Vector x(3);
    T->Transform(ip, x);
    return x;
}

} // namespace dft

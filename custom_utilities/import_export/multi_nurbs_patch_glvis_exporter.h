//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGLVisMETRIC_APPLICATION_MULTI_NURBS_PATCH_GLVIS_EXPORTER_H_INCLUDED)
#define  KRATOS_ISOGLVisMETRIC_APPLICATION_MULTI_NURBS_PATCH_GLVIS_EXPORTER_H_INCLUDED

// System includes
#include <vector>
#include <fstream>
#include <iomanip>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/nurbs/bsplines_fespace.h"
#include "custom_utilities/import_export/multipatch_exporter.h"

namespace Kratos
{

/**
Export NURBS patch/multipatch to GLVis to visualize with NURBS toolbox by M. Spink
 */
template<int TDim>
class MultiNURBSPatchGLVisExporterWriter : public MultiPatchExporter<TDim>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MultiNURBSPatchGLVisExporterWriter);

    /// Type definition
    typedef MultiPatchExporter<TDim> BaseType;
    typedef typename BaseType::knot_container_t knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;

    /// Default constructor
    MultiNURBSPatchGLVisExporterWriter() : BaseType() {}

    /// Destructor
    virtual ~MultiNURBSPatchGLVisExporterWriter() {}

    /// Export a single patch
    virtual void Export(typename Patch<TDim>::Pointer pPatch, const std::string& filename) const
    {
        BaseType::Export(pPatch, filename);
    }

    /// Export a multipatch
    virtual void Export(typename MultiPatch<TDim>::Pointer pMultiPatch, const std::string& filename) const
    {
        std::ofstream outfile;
        outfile.open(filename, std::ios::out);
        outfile << std::setprecision(BaseType::Accuracy());

        outfile << "MFEM NURBS mesh v1.0\n\n";

        outfile << "#\n";
        outfile << "# MFEM Geometry Types (see mesh/geom.hpp):\n";
        outfile << "#\n";
        outfile << "# SEGMENT     = 1\n";
        outfile << "# SQUARE      = 3\n";
        outfile << "# CUBE        = 5\n";
        outfile << "#\n\n";

        outfile << "dimension\n" << TDim << "\n\n";

        std::size_t nvertices;
        std::vector<std::vector<std::size_t> > elements;
        std::vector<std::vector<std::size_t> > boundary;
        std::vector<std::tuple<std::size_t, std::size_t, std::size_t, int> > edges;
        std::map<std::size_t, std::vector<double> > knotvec;
        // KRATOS_WATCH(__LINE__) // error with 1D here
        if (TDim == 1)
        {
            KRATOS_THROW_ERROR(std::logic_error, "Here there will be an error, no more proceeding", "")
        }
        this->GenerateCornerTopology(*pMultiPatch, nvertices, elements, boundary, edges, knotvec);
        // KRATOS_WATCH(__LINE__)

        outfile << "elements\n" << elements.size() << "\n";
        for (std::size_t i = 0; i < elements.size(); ++i)
        {
            outfile << "1";
            if (elements[i].size() == 4)
                outfile << " 3";
            else if (elements[i].size() == 8)
                outfile << " 5";
            else
                KRATOS_THROW_ERROR(std::logic_error, "Invalid number of nodes for an element:", elements[i].size())

            for (std::size_t j = 0; j < elements[i].size(); ++j)
                outfile << " " << elements[i][j];
            outfile << "\n";
        }
        outfile << "\n";

        std::size_t nboundary = 0;
        for (std::size_t i = 0; i < edges.size(); ++i)
            if (std::get<3>(edges[i]) != 0)
                ++nboundary;
        outfile << "boundary\n" << nboundary << "\n";
        for (std::size_t i = 0; i < edges.size(); ++i)
        {
            if (std::get<3>(edges[i]) != 0)
                outfile << "1 1 " << std::get<0>(edges[i]) << " " << std::get<1>(edges[i]) << "\n";
        }
        outfile << "\n\n";

        outfile << "edges\n" << edges.size() << "\n";
        for (std::size_t i = 0; i < edges.size(); ++i)
        {
            outfile << std::get<2>(edges[i]) << " "
                    << std::get<0>(edges[i]) << " " << std::get<1>(edges[i]) << "\n";
        }
        outfile << "\n\n";

        outfile << "vertices\n" << nvertices << "\n\n";

        outfile << "patches\n\n";
        for (typename MultiPatch<TDim>::PatchContainerType::iterator it = pMultiPatch->Patches().begin();
                it != pMultiPatch->Patches().end(); ++it)
        {
            outfile << "# patch " << it->Id() << "\n\n";
            outfile << "knotvectors\n" << TDim << "\n";

            if (it->pFESpace()->Type() != BSplinesFESpace<TDim>::StaticType())
                KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "does not support non-NURBS patch")

            typename BSplinesFESpace<TDim>::Pointer pFESpace = boost::dynamic_pointer_cast<BSplinesFESpace<TDim> >(it->pFESpace());
            for (std::size_t dim = 0; dim < TDim; ++dim)
            {
                outfile << pFESpace->Order(dim) << " " << pFESpace->Number(dim);
                for (std::size_t i = 0; i < pFESpace->KnotVector(dim).size(); ++i)
                    outfile << " " << pFESpace->KnotVector(dim)[i];
                outfile << "\n";
            }
            outfile << "\n";

            outfile << "dimension\n" << TDim << "\n\n";

            typename ControlGrid<ControlPoint<double> >::Pointer pControlGrid = it->pControlPointGridFunction()->pControlGrid();
            outfile << "controlpoints\n";
            for (std::size_t i = 0; i < pControlGrid->size(); ++i)
            {
                for (std::size_t dim = 0; dim < TDim; ++dim)
                    outfile << " " << (*pControlGrid)[i][dim];
                outfile << " " << (*pControlGrid)[i][3] << "\n";
            }
            outfile << "\n";
        }

        outfile << std::endl;

        outfile.close();

        std::cout <<" Multipatch is exported to " << filename << " successfully" << std::endl;
    }

private:

    /// Generate the multipatch topology that can be read in GLVis
    /// It is noted that only the BSplines patch is supported. If the patch is hierarchical BSplines or TSplines, this will generate error
    void GenerateCornerTopology(const MultiPatch<TDim>& r_multipatch,
        std::size_t& nvertices,
        std::vector<std::vector<std::size_t> >& elements,
        std::vector<std::vector<std::size_t> >& boundary,
        std::vector<std::tuple<std::size_t, std::size_t, std::size_t, int> >& edges,
        std::map<std::size_t, std::vector<double> >& knotvecs) const
    {
        typedef typename Patch<TDim>::vertex_t vertex_t;
        typedef typename Patch<TDim>::edge_t edge_t;
        typedef typename Patch<TDim>::face_t face_t;
        typedef typename Patch<TDim>::volume_t volume_t;
        typedef typename MultiPatch<TDim>::PatchContainerType PatchContainerType;

        std::map<std::size_t, std::vector<vertex_t> > patch_vertices;
        std::map<std::size_t, std::vector<edge_t> > patch_edges;
        std::map<std::size_t, std::vector<face_t> > patch_faces;
        std::map<std::size_t, std::vector<volume_t> > patch_volumes;
        std::map<std::size_t, std::vector<std::size_t> > patch_knotv;

        // generate vertices, edges, faces for all the patches
        std::size_t start_vertex_id = 0;
        std::size_t start_knotv_id = 0;
        std::map<std::size_t, std::vector<double> > all_knotvec;
        for (typename PatchContainerType::const_iterator it = r_multipatch.begin(); it != r_multipatch.end(); ++it)
        {
            const std::size_t& id = it->Id();
            it->GenerateTopolgyData(start_vertex_id, patch_vertices[id], patch_edges[id], patch_faces[id], patch_volumes[id], start_knotv_id, patch_knotv[id]);

            typename BSplinesFESpace<TDim>::Pointer pFESpace = boost::dynamic_pointer_cast<BSplinesFESpace<TDim> >(it->pFESpace());

            // collect the knot vector in respective dimension
            for (std::size_t i = 0; i < TDim; ++i)
            {
                for (std::size_t j = 0; j < pFESpace->KnotVector(i).size(); ++j)
                    all_knotvec[patch_knotv[id][i]].push_back(pFESpace->KnotVector(i)[j]);
            }
        }

        #ifdef DEBUG_GLVIS_EXPORT
        for (typename PatchContainerType::const_iterator it = r_multipatch.begin(); it != r_multipatch.end(); ++it)
        {
            const std::size_t& id = it->Id();

            std::cout << "edge (p1) for patch " << id << std::endl;
            for (std::size_t i = 0; i < patch_edges[id].size(); ++i)
            {
                std::cout << std::get<2>(patch_edges[id][i])
                          << " " << std::get<0>(patch_edges[id][i])
                          << " " << std::get<1>(patch_edges[id][i])
                          << " " << std::get<3>(patch_edges[id][i])
                          << std::endl;
            }
        }
        #endif

        // for all patch, account for the corners and then renumbering the vertex and edge
        for (typename PatchContainerType::const_iterator it = r_multipatch.begin(); it != r_multipatch.end(); ++it)
        {
            const std::size_t& id = it->Id();

            if (it->pNeighbor(_LEFT_) != NULL)
            {
                const std::size_t& other_id = it->pNeighbor(_LEFT_)->Id();
                // TODO accommodate for different type of synchronization. Because the left neighbor does not always tie by the right boundary.
                SynchronizeVertices(TDim, _LEFT_, patch_vertices[id], _RIGHT_, patch_vertices[other_id], patch_edges[other_id], patch_faces[other_id], patch_volumes[other_id]);
                if (TDim == 2)
                {
                    patch_knotv[other_id][1] = patch_knotv[id][1];
                    std::get<2>(patch_edges[other_id][0]) = patch_knotv[id][1]; // change the knot index vector
                    std::get<2>(patch_edges[other_id][1]) = patch_knotv[id][1];
                    std::get<3>(patch_edges[other_id][1]) = 0; // disable the boundary flag
                    std::get<3>(patch_edges[id][0]) = 0; // disable the boundary flag
                }
                else if (TDim == 3)
                {
                    patch_knotv[other_id][1] = patch_knotv[id][1];
                    patch_knotv[other_id][2] = patch_knotv[id][2];
                    // TODO change knot vector index for edge
                    // TODO disable the boundary flag for edges and faces
                }
            }

            if (it->pNeighbor(_RIGHT_) != NULL)
            {
                const std::size_t& other_id = it->pNeighbor(_RIGHT_)->Id();
                // TODO accommodate for different type of synchronization. Because the right neighbor does not always tie by the left boundary.
                SynchronizeVertices(TDim, _RIGHT_, patch_vertices[id], _LEFT_, patch_vertices[other_id], patch_edges[other_id], patch_faces[other_id], patch_volumes[other_id]);
                if (TDim == 2)
                {
                    patch_knotv[other_id][1] = patch_knotv[id][1];
                    std::get<2>(patch_edges[other_id][0]) = patch_knotv[id][1];
                    std::get<2>(patch_edges[other_id][1]) = patch_knotv[id][1];
                    std::get<3>(patch_edges[other_id][0]) = 0; // disable the boundary flag
                    std::get<3>(patch_edges[id][1]) = 0; // disable the boundary flag
                }
                else if (TDim == 3)
                {
                    patch_knotv[other_id][1] = patch_knotv[id][1];
                    patch_knotv[other_id][2] = patch_knotv[id][2];
                    // TODO change knot vector index for edge
                    // TODO disable the boundary flag for edges and faces
                }
            }

            if (it->pNeighbor(_TOP_) != NULL)
            {
                const std::size_t& other_id = it->pNeighbor(_TOP_)->Id();
                SynchronizeVertices(TDim, _TOP_, patch_vertices[id], _BOTTOM_, patch_vertices[other_id], patch_edges[other_id], patch_faces[other_id], patch_volumes[other_id]);
                if (TDim == 2)
                {
                    patch_knotv[other_id][0] = patch_knotv[id][0];
                    std::get<2>(patch_edges[other_id][2]) = patch_knotv[id][0];
                    std::get<2>(patch_edges[other_id][3]) = patch_knotv[id][0];
                    std::get<3>(patch_edges[other_id][2]) = 0; // disable the boundary flag
                    std::get<3>(patch_edges[id][3]) = 0; // disable the boundary flag
                }
                else if (TDim == 3)
                {
                    patch_knotv[other_id][0] = patch_knotv[id][0];
                    patch_knotv[other_id][1] = patch_knotv[id][1];
                    // TODO change knot vector index for edge
                    // TODO disable the boundary flag for edges and faces
                }
            }

            if (it->pNeighbor(_BOTTOM_) != NULL)
            {
                const std::size_t& other_id = it->pNeighbor(_BOTTOM_)->Id();
                SynchronizeVertices(TDim, _BOTTOM_, patch_vertices[id], _TOP_, patch_vertices[other_id], patch_edges[other_id], patch_faces[other_id], patch_volumes[other_id]);
                if (TDim == 2)
                {
                    patch_knotv[other_id][0] = patch_knotv[id][0];
                    std::get<2>(patch_edges[other_id][2]) = patch_knotv[id][0];
                    std::get<2>(patch_edges[other_id][3]) = patch_knotv[id][0];
                    std::get<3>(patch_edges[other_id][3]) = 0; // disable the boundary flag
                    std::get<3>(patch_edges[id][2]) = 0; // disable the boundary flag
                }
                else if (TDim == 3)
                {
                    patch_knotv[other_id][0] = patch_knotv[id][0];
                    patch_knotv[other_id][1] = patch_knotv[id][1];
                    // TODO change knot vector index for edge
                    // TODO disable the boundary flag for edges and faces
                }
            }

            if (it->pNeighbor(_FRONT_) != NULL)
            {
                const std::size_t& other_id = it->pNeighbor(_FRONT_)->Id();
                SynchronizeVertices(TDim, _FRONT_, patch_vertices[id], _BACK_, patch_vertices[other_id], patch_edges[other_id], patch_faces[other_id], patch_volumes[other_id]);
                patch_knotv[other_id][0] = patch_knotv[id][0];
                patch_knotv[other_id][2] = patch_knotv[id][2];
                // TODO change knot vector index for edge
                // TODO disable the boundary flag for edges and faces
            }

            if (it->pNeighbor(_BACK_) != NULL)
            {
                const std::size_t& other_id = it->pNeighbor(_BACK_)->Id();
                SynchronizeVertices(TDim, _BACK_, patch_vertices[id], _FRONT_, patch_vertices[other_id], patch_edges[other_id], patch_faces[other_id], patch_volumes[other_id]);
                patch_knotv[other_id][0] = patch_knotv[id][0];
                patch_knotv[other_id][2] = patch_knotv[id][2];
                // TODO change knot vector index for edge
                // TODO disable the boundary flag for edges and faces
            }
        }

        #ifdef DEBUG_GLVIS_EXPORT
        for (typename PatchContainerType::const_iterator it = r_multipatch.begin(); it != r_multipatch.end(); ++it)
        {
            const std::size_t& id = it->Id();

            std::cout << "edge (p2) for patch " << id << std::endl;
            for (std::size_t i = 0; i < patch_edges[id].size(); ++i)
            {
                std::cout << std::get<2>(patch_edges[id][i])
                          << " " << std::get<0>(patch_edges[id][i])
                          << " " << std::get<1>(patch_edges[id][i])
                          << " " << std::get<3>(patch_edges[id][i])
                          << std::endl;
            }
        }
        #endif

        // collect all the vertices in the multipatch
        std::set<vertex_t> all_vertices;
        for (typename PatchContainerType::const_iterator it = r_multipatch.begin(); it != r_multipatch.end(); ++it)
        {
            const std::size_t& id = it->Id();
            all_vertices.insert(patch_vertices[id].begin(), patch_vertices[id].end());
        }

        // reassign each id a new index
        std::map<vertex_t, vertex_t> old_to_new;
        start_vertex_id = 0;
        for (typename std::set<vertex_t>::iterator it = all_vertices.begin(); it != all_vertices.end(); ++it)
            old_to_new[*it] = start_vertex_id++;

        // finally reassign the new id for all patches
        for (typename PatchContainerType::const_iterator it = r_multipatch.begin(); it != r_multipatch.end(); ++it)
        {
            const std::size_t& id = it->Id();

            for (std::size_t i = 0; i < patch_vertices[id].size(); ++i)
            {
                patch_vertices[id][i] = old_to_new[patch_vertices[id][i]];
            }

            for (std::size_t i = 0; i < patch_edges[id].size(); ++i)
            {
                std::get<0>(patch_edges[id][i]) = old_to_new[std::get<0>(patch_edges[id][i])];
                std::get<1>(patch_edges[id][i]) = old_to_new[std::get<1>(patch_edges[id][i])];
            }

            for (std::size_t i = 0; i < patch_faces[id].size(); ++i)
            {
                std::get<0>(patch_faces[id][i]) = old_to_new[std::get<0>(patch_faces[id][i])];
                std::get<1>(patch_faces[id][i]) = old_to_new[std::get<1>(patch_faces[id][i])];
                std::get<2>(patch_faces[id][i]) = old_to_new[std::get<2>(patch_faces[id][i])];
                std::get<3>(patch_faces[id][i]) = old_to_new[std::get<3>(patch_faces[id][i])];
            }

            for (std::size_t i = 0; i < patch_volumes[id].size(); ++i)
            {
                std::get<0>(patch_volumes[id][i]) = old_to_new[std::get<0>(patch_volumes[id][i])];
                std::get<1>(patch_volumes[id][i]) = old_to_new[std::get<1>(patch_volumes[id][i])];
                std::get<2>(patch_volumes[id][i]) = old_to_new[std::get<2>(patch_volumes[id][i])];
                std::get<3>(patch_volumes[id][i]) = old_to_new[std::get<3>(patch_volumes[id][i])];
                std::get<4>(patch_volumes[id][i]) = old_to_new[std::get<4>(patch_volumes[id][i])];
                std::get<5>(patch_volumes[id][i]) = old_to_new[std::get<5>(patch_volumes[id][i])];
                std::get<6>(patch_volumes[id][i]) = old_to_new[std::get<6>(patch_volumes[id][i])];
                std::get<7>(patch_volumes[id][i]) = old_to_new[std::get<7>(patch_volumes[id][i])];
            }
        }

        // collect all the knot vector index in all the patches
        std::set<std::size_t> all_knotvs;
        for (typename PatchContainerType::const_iterator it = r_multipatch.begin(); it != r_multipatch.end(); ++it)
        {
            const std::size_t& id = it->Id();
            for (std::size_t i = 0; i < TDim; ++i)
                all_knotvs.insert(patch_knotv[id][i]);
        }

        #ifdef DEBUG_GLVIS_EXPORT
        std::cout << "all_knotvs:";
        for (std::set<std::size_t>::iterator it = all_knotvs.begin(); it != all_knotvs.end(); ++it)
            std::cout << " " << *it;
        std::cout << std::endl;
        #endif

        // reassign each knot vector a new index
        std::map<std::size_t, std::size_t> old_to_new_knotv;
        start_knotv_id = 0;
        for (std::set<std::size_t>::iterator it = all_knotvs.begin(); it != all_knotvs.end(); ++it)
            old_to_new_knotv[*it] = start_knotv_id++;

        // finally reassign the new id for all patches
        for (typename PatchContainerType::const_iterator it = r_multipatch.begin(); it != r_multipatch.end(); ++it)
        {
            const std::size_t& id = it->Id();
            for (std::size_t i = 0; i < TDim; ++i)
                patch_knotv[id][i] = old_to_new_knotv[patch_knotv[id][i]];
        }

        // reassign the knot vector index in each edge
        for (typename PatchContainerType::const_iterator it = r_multipatch.begin(); it != r_multipatch.end(); ++it)
        {
            const std::size_t& id = it->Id();
            for (std::size_t i = 0; i < patch_edges[id].size(); ++i)
            {
                std::get<2>(patch_edges[id][i]) = old_to_new_knotv[std::get<2>(patch_edges[id][i])];
            }
        }

        // assign the knot vector accordingly
        for (typename PatchContainerType::const_iterator it = r_multipatch.begin(); it != r_multipatch.end(); ++it)
        {
            const std::size_t& id = it->Id();
            for (std::size_t i = 0; i < TDim; ++i)
                knotvecs[ old_to_new_knotv[patch_knotv[id][i]] ] = all_knotvec[ patch_knotv[id][i] ];
        }

        // collect all the edges in all the patches
        std::set<edge_t> all_edges;
        for (typename PatchContainerType::const_iterator it = r_multipatch.begin(); it != r_multipatch.end(); ++it)
        {
            const std::size_t& id = it->Id();
            for (std::size_t i = 0; i < patch_edges[id].size(); ++i)
            {
                all_edges.insert(patch_edges[id][i]);
            }
        }
        edges.assign(all_edges.begin(), all_edges.end());

        ///////////////////////////////////////////////////

        // export the information
        nvertices = all_vertices.size();

        if (TDim == 2)
        {
            for (typename PatchContainerType::const_iterator it = r_multipatch.begin(); it != r_multipatch.end(); ++it)
            {
                const std::size_t& id = it->Id();

                std::vector<std::size_t> elem(4);
                elem[0] = std::get<0>(patch_faces[id][0]);
                elem[1] = std::get<1>(patch_faces[id][0]);
                elem[2] = std::get<3>(patch_faces[id][0]); // here we switch the role of the vertex because GLVis only accepts the sequence 0-1-3-2 for quadrilateral
                elem[3] = std::get<2>(patch_faces[id][0]);

                elements.push_back(elem);
            }
        }
        else if (TDim == 3)
        {
            // TODO
            KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented for 3D", "")
        }
    }

    /// Synchronize from 1->2
    template<typename vertex_t, typename edge_t, typename face_t, typename volume_t>
    void SynchronizeVertices( const int& Dim,
        const BoundarySide& side1,
        std::vector<vertex_t>& vertices1,
        const BoundarySide& side2,
        std::vector<vertex_t>& vertices2,
        std::vector<edge_t>& edges2,
        std::vector<face_t>& faces2,
        std::vector<volume_t>& volumes2 ) const
    {
        if (vertices1.size() != vertices2.size())
            KRATOS_THROW_ERROR(std::logic_error, "The number of vertices is not compatible", "")

        std::vector<int> map = GetJointMapping(Dim, side1, side2);

        std::map<std::size_t, std::size_t> old_to_new;
        for (std::size_t i = 0; i < vertices2.size(); ++i)
            old_to_new[vertices2[i]] = vertices2[i];
        for (std::size_t i = 0; i < map.size()/2; ++i)
        {
            old_to_new[vertices2[map[i*2+1]]] = vertices1[map[i*2]];
            vertices2[map[i*2+1]] = vertices1[map[i*2]];
        }

        for (std::size_t i = 0; i < edges2.size(); ++i)
        {
            std::get<0>(edges2[i]) = old_to_new[std::get<0>(edges2[i])];
            std::get<1>(edges2[i]) = old_to_new[std::get<1>(edges2[i])];
        }

        for (std::size_t i = 0; i < faces2.size(); ++i)
        {
            std::get<0>(faces2[i]) = old_to_new[std::get<0>(faces2[i])];
            std::get<1>(faces2[i]) = old_to_new[std::get<1>(faces2[i])];
            std::get<2>(faces2[i]) = old_to_new[std::get<2>(faces2[i])];
            std::get<3>(faces2[i]) = old_to_new[std::get<3>(faces2[i])];
        }

        for (std::size_t i = 0; i < volumes2.size(); ++i)
        {
            std::get<0>(volumes2[i]) = old_to_new[std::get<0>(volumes2[i])];
            std::get<1>(volumes2[i]) = old_to_new[std::get<1>(volumes2[i])];
            std::get<2>(volumes2[i]) = old_to_new[std::get<2>(volumes2[i])];
            std::get<3>(volumes2[i]) = old_to_new[std::get<3>(volumes2[i])];
            std::get<4>(volumes2[i]) = old_to_new[std::get<4>(volumes2[i])];
            std::get<5>(volumes2[i]) = old_to_new[std::get<5>(volumes2[i])];
            std::get<6>(volumes2[i]) = old_to_new[std::get<6>(volumes2[i])];
            std::get<7>(volumes2[i]) = old_to_new[std::get<7>(volumes2[i])];
        }
    }

    /// Support function to generate an array to map node from left to right, top to bottom, etc...
    /// The node indexing convention follows Burstedde paper of P4est
    std::vector<int> GetJointMapping(const int& dim, const BoundarySide& side1, const BoundarySide& side2) const
    {
        if (dim == 2)
        {
            if (side1 == _LEFT_)
            {
                if (side2 == _RIGHT_)
                    return std::vector<int>{ 0, 1, /**/ 2, 3 };
                else
                    //TODO
                    KRATOS_THROW_ERROR(std::logic_error, "Mapping for other side is not yet implemented.", "")
            }
            else if (side1 == _RIGHT_)
            {
                if (side2 == _LEFT_)
                    return std::vector<int>{ 1, 0, /**/ 3, 2 };
                else
                    //TODO
                    KRATOS_THROW_ERROR(std::logic_error, "Mapping for other side is not yet implemented.", "")
            }
            else if (side1 == _TOP_)
            {
                if (side2 == _BOTTOM_)
                    return std::vector<int>{ 2, 0, /**/ 3, 1 };
                else
                    // TODO
                    KRATOS_THROW_ERROR(std::logic_error, "Mapping for other side is not yet implemented.", "")
            }
            else if (side1 == _BOTTOM_)
            {
                if (side2 == _TOP_)
                    return std::vector<int>{ 0, 2, /**/ 1, 3 };
                else
                    // TODO
                    KRATOS_THROW_ERROR(std::logic_error, "Mapping for other side is not yet implemented.", "")
            }
        }
        else if (dim == 3)
        {
            // TODO
            KRATOS_THROW_ERROR(std::logic_error, "Mapping for 3D is not implemented yet.", "")
        }
    }

}; // end class MultiNURBSPatchGLVisExporterWriter



class MultiNURBSPatchGLVisExporter
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(MultiNURBSPatchGLVisExporter);

    template<int TDim>
    static void Export(typename Patch<TDim>::Pointer pPatch, const std::string& filename)
    {
        MultiNURBSPatchGLVisExporterWriter<TDim> dummy;
        dummy.Export(pPatch, filename);
    }

    template<int TDim>
    static void Export(typename MultiPatch<TDim>::Pointer pMultiPatch, const std::string& filename)
    {
        MultiNURBSPatchGLVisExporterWriter<TDim> dummy;
        dummy.Export(pMultiPatch, filename);
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MultiNURBSPatchGLVisExporter";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};


/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const MultiNURBSPatchGLVisExporter& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGLVisMETRIC_APPLICATION_MULTI_NURBS_PATCH_GLVIS_EXPORTER_H_INCLUDED defined


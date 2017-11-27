//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 3 June 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_DOMAIN_MANAGER_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_DOMAIN_MANAGER_H_INCLUDED

// System includes
#include <vector>
#include <set>
#include <map>
#include <iterator>
#include <iostream>
#include <fstream>

// Project includes
#include "includes/define.h"

namespace Kratos
{

/**
This class represents a union of domains (rectangle/cuboid) in the parameter space of NURBS in 2D. THis is useful to manage the refinement domain in the hierarchical NURBS mesh
 */
class DomainManager
{

public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(DomainManager);

    /// Type definition
    typedef std::set<double> coords_container_t;

    /// Default constructor
    DomainManager(const std::size_t& Id) : mId(Id) {}

    /// Destructor
    virtual ~DomainManager() {}

    /// Get the Id of this domain manager
    const std::size_t& Id() const {return mId;}

    /// Set the id for this domain manager
    void SetId(const std::size_t& Id) {mId = Id;}

    /// Fill the internal array of X & Y-coordinates. It must be done before cells are added to this container.
    virtual void AddXcoord(const double& X) {mXcoords.insert(X);}
    virtual void AddYcoord(const double& Y) {mYcoords.insert(Y);}
    virtual void AddZcoord(const double& Z) {mZcoords.insert(Z);}

    /// Add the rectangle to the set
    virtual void AddCell(const double& Xmin,
                         const double& Xmax,
                         const double& Ymin,
                         const double& Ymax)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Call virtual method", __FUNCTION__)
    }

    /// Add the cuboid to the set
    virtual void AddCell(const double& Xmin,
                         const double& Xmax,
                         const double& Ymin,
                         const double& Ymax,
                         const double& Zmin,
                         const double& Zmax)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Call virtual method", __FUNCTION__)
    }

    /// Check if a rectangle if it is contained in the set.
    virtual bool IsInside(const double& Xmin,
                          const double& Xmax,
                          const double& Ymin,
                          const double& Ymax) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Call virtual method", __FUNCTION__)
    }

    /// Check if a rectangle if it is contained in the set.
    virtual bool IsInside(const double& Xmin,
                          const double& Xmax,
                          const double& Ymin,
                          const double& Ymax,
                          const double& Zmin,
                          const double& Zmax) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Call virtual method", __FUNCTION__)
    }

    /// Export the domain to matlab for visualization
    virtual void ExportDomain(std::string fn, std::string color, double distance) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Call virtual method", __FUNCTION__)
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const {}
    virtual void PrintData(std::ostream& rOStream) const {}

protected:

    coords_container_t mXcoords;
    coords_container_t mYcoords;
    coords_container_t mZcoords;

private:
    std::size_t mId;
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const DomainManager& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_DOMAIN_MANAGER_H_INCLUDED defined


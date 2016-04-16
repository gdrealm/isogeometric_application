//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 3 June 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_HN_CUBOID_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_HN_CUBOID_H_INCLUDED

// System includes
#include <vector>
#include <iostream>

// Project includes
#include "includes/define.h"

namespace Kratos
{

/**
    Represent a set of cuboid in the same level in hierarchical NURBS hierarchy
 */
class HnCuboid
{
private:
    class cuboid
    {
        public:
            KRATOS_CLASS_POINTER_DEFINITION(cuboid);
            cuboid(double Xmin, double Xmax, double Ymin, double Ymax, double Zmin, double Zmax) : mXmin(Xmin), mXmax(Xmax), mYmin(Ymin), mYmax(Ymax), mZmin(Zmin), mZmax(Zmax) {}
            virtual ~cuboid() {}

            double Xmin() const {return mXmin;}
            double Xmax() const {return mXmax;}
            double Ymin() const {return mYmin;}
            double Ymax() const {return mYmax;}
            double Zmin() const {return mZmin;}
            double Zmax() const {return mZmax;}

            bool IsInside(double X, double Y, double Z) const
            {
                return (X >= mXmin) && (X <= mXmax) && (Y >= mYmin) && (Y <= mYmax) && (Z >= mZmin) && (Z <= mZmax);
            }

        private:
            double mXmin, mXmax, mYmin, mYmax, mZmin, mZmax;
    };

public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(HnCuboid);

    /// Type definitions
    typedef cuboid::Pointer cuboid_t;
    typedef std::map<unsigned int, std::vector<cuboid_t> > cuboid_container_t;

    /// Default constructor
    HnCuboid() : mNumberOfSet(0) {}

    /// Destructor
    virtual ~HnCuboid() {}

    /// Add the cuboid to the cuboid set
    void AddCuboid(const double& Xmin,
                   const double& Xmax,
                   const double& Ymin,
                   const double& Ymax,
                   const double& Zmin,
                   const double& Zmax)
    {
        for(cuboid_container_t::iterator it = mp_cuboid_set.begin(); it != mp_cuboid_set.end(); ++it)
        {
            for(std::size_t i = 0; i < it->second.size(); ++i)
            {
                if( it->second[i]->IsInside(Xmin, Ymin, Zmin)
                 || it->second[i]->IsInside(Xmin, Ymin, Zmax)
                 || it->second[i]->IsInside(Xmin, Ymax, Zmin)
                 || it->second[i]->IsInside(Xmin, Ymax, Zmax)
                 || it->second[i]->IsInside(Xmax, Ymin, Zmin)
                 || it->second[i]->IsInside(Xmax, Ymin, Zmax)
                 || it->second[i]->IsInside(Xmax, Ymax, Zmin)
                 || it->second[i]->IsInside(Xmax, Ymax, Zmax))
                 {
                    it->second.push_back(cuboid_t(new cuboid(Xmin, Xmax, Ymin, Ymax, Zmin, Zmax)));
                    return;
                 }
            }
        }
        ++mNumberOfSet;
        mp_cuboid_set[mNumberOfSet].push_back(cuboid_t(new cuboid(Xmin, Xmax, Ymin, Ymax, Zmin, Zmax)));
    }

    /// Check if a cuboid if it is contained in the cuboid set
    bool IsInside(const double& Xmin,
                  const double& Xmax,
                  const double& Ymin,
                  const double& Ymax,
                  const double& Zmin,
                  const double& Zmax) const
    {
        bool point_is_inside;
        for(cuboid_container_t::const_iterator it = mp_cuboid_set.begin(); it != mp_cuboid_set.end(); ++it)
        {
            point_is_inside = false;
            for(std::size_t i = 0; i < it->second.size(); ++i)
            {
                point_is_inside = point_is_inside || it->second[i]->IsInside(Xmin, Ymin, Zmin);
                if(point_is_inside == true) break;
            }
            if(point_is_inside == false) continue;

            point_is_inside = false;
            for(std::size_t i = 0; i < it->second.size(); ++i)
            {
                point_is_inside = point_is_inside || it->second[i]->IsInside(Xmin, Ymin, Zmax);
                if(point_is_inside == true) break;
            }
            if(point_is_inside == false) continue;

            point_is_inside = false;
            for(std::size_t i = 0; i < it->second.size(); ++i)
            {
                point_is_inside = point_is_inside || it->second[i]->IsInside(Xmin, Ymax, Zmin);
                if(point_is_inside == true) break;
            }
            if(point_is_inside == false) continue;

            point_is_inside = false;
            for(std::size_t i = 0; i < it->second.size(); ++i)
            {
                point_is_inside = point_is_inside || it->second[i]->IsInside(Xmin, Ymax, Zmax);
                if(point_is_inside == true) break;
            }
            if(point_is_inside == false) continue;

            point_is_inside = false;
            for(std::size_t i = 0; i < it->second.size(); ++i)
            {
                point_is_inside = point_is_inside || it->second[i]->IsInside(Xmax, Ymin, Zmin);
                if(point_is_inside == true) break;
            }
            if(point_is_inside == false) continue;

            point_is_inside = false;
            for(std::size_t i = 0; i < it->second.size(); ++i)
            {
                point_is_inside = point_is_inside || it->second[i]->IsInside(Xmax, Ymin, Zmax);
                if(point_is_inside == true) break;
            }
            if(point_is_inside == false) continue;

            point_is_inside = false;
            for(std::size_t i = 0; i < it->second.size(); ++i)
            {
                point_is_inside = point_is_inside || it->second[i]->IsInside(Xmax, Ymax, Zmin);
                if(point_is_inside == true) break;
            }
            if(point_is_inside == false) continue;

            point_is_inside = false;
            for(std::size_t i = 0; i < it->second.size(); ++i)
            {
                point_is_inside = point_is_inside || it->second[i]->IsInside(Xmax, Ymax, Zmax);
                if(point_is_inside == true) break;
            }
            if(point_is_inside == false) continue;

            return true;
        }
        return false;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "Number of cuboid set: " << mNumberOfSet << ",";
        for(cuboid_container_t::const_iterator it = mp_cuboid_set.begin(); it != mp_cuboid_set.end(); ++it)
            rOStream << " " << it->second.size();
    }

private:
    unsigned int mNumberOfSet;
    cuboid_container_t mp_cuboid_set;
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const HnCuboid& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_HN_CUBOID_H_INCLUDED


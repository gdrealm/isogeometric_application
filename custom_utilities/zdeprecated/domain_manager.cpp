#include "domain_manager.h"

namespace Kratos
{

    void DomainManager::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DomainManager in 2D";
    }

    void DomainManager::PrintData(std::ostream& rOStream) const
    {
        rOStream << "X-coordinates:";
        for(std::set<double>::const_iterator it = mXcoords.begin(); it != mXcoords.end(); ++it)
            rOStream << " " << *it;
        rOStream << std::endl;

        rOStream << "Y-coordinates:";
        for(std::set<double>::const_iterator it = mYcoords.begin(); it != mYcoords.end(); ++it)
            rOStream << " " << *it;
        rOStream << std::endl;

        rOStream << "Cells:" << std::endl;
        for(std::map<std::size_t, std::set<std::size_t> >::const_iterator it1 = mActiveCells.begin(); it1 != mActiveCells.end(); ++it1)
        {
            rOStream << " column " << it1->first << ":";
            for(std::set<std::size_t>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
                rOStream << " " << *it2;
            rOStream << std::endl;
        }
    }

    void DomainManager::AddCell(const double& Xmin,
                                const double& Xmax,
                                const double& Ymin,
                                const double& Ymax)
    {
        std::set<double>::iterator it_x1 = mXcoords.find(Xmin);
        std::set<double>::iterator it_x2 = mXcoords.find(Xmax);
        std::set<double>::iterator it_y1 = mYcoords.find(Ymin);
        std::set<double>::iterator it_y2 = mYcoords.find(Ymax);

        if(it_x1 == mXcoords.end() || it_x2 == mXcoords.end())
            KRATOS_THROW_ERROR(std::runtime_error, "Cell does not align with x-coordinates", "")

        if(it_y1 == mYcoords.end() || it_y2 == mYcoords.end())
            KRATOS_THROW_ERROR(std::runtime_error, "Cell does not align with y-coordinates", "")

        std::size_t i1 = std::distance(mXcoords.begin(), it_x1);
        std::size_t i2 = std::distance(mXcoords.begin(), it_x2);
        std::size_t j1 = std::distance(mYcoords.begin(), it_y1);
        std::size_t j2 = std::distance(mYcoords.begin(), it_y2);

        for(std::size_t i = i1; i < i2; ++i)
            for(std::size_t j = j1; j < j2; ++j)
                mActiveCells[i].insert(j);
    }

    bool DomainManager::IsInside(const double& Xmin,
                                 const double& Xmax,
                                 const double& Ymin,
                                 const double& Ymax) const
    {
        double tol = 1.0e-10;

        // find the lower bound for the Xmin and upper bound for Xmax
        std::size_t i1 = 0, i2 = 0;
        for(std::set<double>::iterator it = mXcoords.begin(); it != mXcoords.end(); ++it)
            if(Xmin > *it - tol)
                ++i1;
        if(i1 == 0 || i1 == mXcoords.size())
            return false;
        else
            --i1;

        //
        for(std::set<double>::iterator it = mXcoords.begin(); it != mXcoords.end(); ++it)
            if(Xmax > *it + tol)
                ++i2;
        if(i2 == 0 || i2 == mXcoords.size())
            return false;

        // find the lower bound for the Ymin and upper bound for Ymax
        std::size_t j1 = 0, j2 = 0;
        for(std::set<double>::iterator it = mYcoords.begin(); it != mYcoords.end(); ++it)
            if(Ymin > *it - tol)
                ++j1;
        if(j1 == 0 || j1 == mYcoords.size())
            return false;
        else
            --j1;

        //
        for(std::set<double>::iterator it = mYcoords.begin(); it != mYcoords.end(); ++it)
            if(Ymax > *it + tol)
                ++j2;
        if(j2 == 0 || j2 == mYcoords.size())
            return false;

//        KRATOS_WATCH(i1)
//        KRATOS_WATCH(i2)
//        KRATOS_WATCH(j1)
//        KRATOS_WATCH(j2)

        // check if in the bound if all the cells are active
        for(std::size_t i = i1; i < i2; ++i)
        {
            std::map<std::size_t, std::set<std::size_t> >::const_iterator it1 = mActiveCells.find(i);
            if(it1 == mActiveCells.end())
                return false;

            for(std::size_t j = j1; j < j2; ++j)
            {
                std::set<std::size_t>::iterator it2 = it1->second.find(j);
                if(it2 == it1->second.end())
                    return false;
            }
        }

        return true;
    }

    void DomainManager::ExportDomain(std::string fn, std::string color, double z_coordinate)
    {
    
    }

}


//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 4 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_CONTROL_POINT_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_CONTROL_POINT_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>

// External includes 

// Project includes
#include "includes/define.h"

namespace Kratos
{

/**
    Represent a control point in isogeometric mesh topology.
 */
template<typename TDataType>
class ControlPoint
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ControlPoint);

    /// Default constructor
    ControlPoint() : mX(0.0), mY(0.0), mZ(0.0), mW(1.0) {}

    /// Destructor
    virtual ~ControlPoint() {}

    /// X-coordinate
    const TDataType& X0() const {return mX/mW;}

    /// Y-coordinate
    const TDataType& Y0() const {return mY/mW;}

    /// Z-coordinate
    const TDataType& Z0() const {return mZ/mW;}

    /// Weight
    const TDataType& W() const {return mW;}

    /// Set the coordinate. The input is the physical coordinates in 3D space.
    void SetCoordinates(const TDataType& X, const TDataType& Y, const TDataType& Z, const TDataType& W)
    {
        mX = W*X;
        mY = W*Y;
        mZ = W*Z;
        mW = W;
    }

    /// Add to the coordinate. The input is the increment of physical coordinates in 3D space.
    void AddCoordinates(const TDataType& X, const TDataType& Y, const TDataType& Z, const TDataType& W)
    {
        mX += W*X;
        mY += W*Y;
        mZ += W*Z;
        mW += W;
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "(" << X0() << ", " << Y0() << ", " << Z0() << ", " << W() << ")";
    }

private:
    TDataType mX, mY, mZ, mW;
};

/// output stream function
template<class TDataType>
inline std::ostream& operator <<(std::ostream& rOStream, const ControlPoint<TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_CONTROL_POINT_H_INCLUDED


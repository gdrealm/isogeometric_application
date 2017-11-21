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

    /// Constant constructor
    ControlPoint(const double& v) : mX(v), mY(v), mZ(v), mW(v) {}

    /// Destructor
    virtual ~ControlPoint() {}

    /// homogeneous X-coordinate
    TDataType& X() {return mX;}
    const TDataType& X() const {return mX;}

    /// X-coordinate
    const TDataType X0() const {return mX/mW;}

    /// homogeneous Y-coordinate
    TDataType& Y() {return mY;}
    const TDataType& Y() const {return mY;}

    /// Y-coordinate
    const TDataType Y0() const {return mY/mW;}

    /// homogeneous Z-coordinate
    TDataType& Z() {return mZ;}
    const TDataType& Z() const {return mZ;}

    /// Z-coordinate
    const TDataType Z0() const {return mZ/mW;}

    /// Weight
    TDataType& W() {return mW;}
    const TDataType& W() const {return mW;}

    /// Set the coordinate. The input is the physical coordinates in 3D space.
    void SetCoordinates(const TDataType& _X, const TDataType& _Y, const TDataType& _Z, const TDataType& _W)
    {
        mX = _W*_X;
        mY = _W*_Y;
        mZ = _W*_Z;
        mW = _W;
    }

    /// Add to the coordinate. The input is the increment of physical coordinates in 3D space.
    void AddCoordinates(const TDataType& _X, const TDataType& _Y, const TDataType& _Z, const TDataType& _W)
    {
        mX += _W*_X;
        mY += _W*_Y;
        mZ += _W*_Z;
        mW += _W;
    }

    // overload operator []
    TDataType& operator[] (const int& i)
    {
        if (i == 0) return X();
        else if (i == 1) return Y();
        else if (i == 2) return Z();
        else if (i == 3) return W();
    }

    const TDataType& operator[] (const int& i) const
    {
        if (i == 0) return X();
        else if (i == 1) return Y();
        else if (i == 2) return Z();
        else if (i == 3) return W();
    }

    // overload operator ()
    TDataType operator() (const int& i) const
    {
        if (i == 0) return X0();
        else if (i == 1) return Y0();
        else if (i == 2) return Z0();
        else if (i == 3) return W();
    }

    /// Assignment operator
    ControlPoint& operator=(const ControlPoint& rOther)
    {
        this->mX = rOther.mX;
        this->mY = rOther.mY;
        this->mZ = rOther.mZ;
        this->mW = rOther.mW;
        return *this;
    }

    /// Addition operator
    ControlPoint& operator+=(const ControlPoint& rOther)
    {
        this->mX += rOther.mX;
        this->mY += rOther.mY;
        this->mZ += rOther.mZ;
        this->mW += rOther.mW;
        return *this;
    }

    /// Addition operator
    friend ControlPoint operator+(ControlPoint lhs, const ControlPoint& rhs)
    {
        lhs += rhs;
        return lhs;
    }

    /// Multiplication operator
    ControlPoint& operator*=(const TDataType& alpha)
    {
        this->mX *= alpha;
        this->mY *= alpha;
        this->mZ *= alpha;
        this->mW *= alpha;
        return *this;
    }

    /// Multiplication operator
    friend ControlPoint operator*(const double& alpha, ControlPoint c)
    {
        c *= alpha;
        return c;
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Control Point";
    }

    void PrintData(std::ostream& rOStream) const
    {
        // print the control point in homogeneous coordinates
        rOStream << "(" << X() << ", " << Y() << ", " << Z() << ", " << W() << ")";
    }

private:
    TDataType mX, mY, mZ, mW;
};

/// output stream function
template<class TDataType>
inline std::ostream& operator <<(std::ostream& rOStream, const ControlPoint<TDataType>& rThis)
{
    // rThis.PrintInfo(rOStream);
    // rOStream << ": ";
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_CONTROL_POINT_H_INCLUDED


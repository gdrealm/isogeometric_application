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
#include "custom_utilities/trans/transformation.h"

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
    ControlPoint() : mWX(0.0), mWY(0.0), mWZ(0.0), mW(1.0) {}

    /// Constant constructor
    ControlPoint(const double& v) : mWX(v), mWY(v), mWZ(v), mW(v) {}

    /// Destructor
    virtual ~ControlPoint() {}

    /// homogeneous X-coordinate
    TDataType& WX() {return mWX;}
    const TDataType& WX() const {return mWX;}

    /// X-coordinate
    const TDataType X() const {return mWX/mW;}

    /// homogeneous Y-coordinate
    TDataType& WY() {return mWY;}
    const TDataType& WY() const {return mWY;}

    /// Y-coordinate
    const TDataType Y() const {return mWY/mW;}

    /// homogeneous Z-coordinate
    TDataType& WZ() {return mWZ;}
    const TDataType& WZ() const {return mWZ;}

    /// Z-coordinate
    const TDataType Z() const {return mWZ/mW;}

    /// Weight
    TDataType& W() {return mW;}
    const TDataType& W() const {return mW;}

    /// Set the coordinate. The input is the physical coordinates in 3D space.
    void SetCoordinates(const TDataType& _X, const TDataType& _Y, const TDataType& _Z, const TDataType& _W)
    {
        mWX = _W*_X;
        mWY = _W*_Y;
        mWZ = _W*_Z;
        mW  = _W;
    }

    /// Add to the coordinate. The input is the increment of physical coordinates in 3D space.
    void AddCoordinates(const TDataType& _X, const TDataType& _Y, const TDataType& _Z, const TDataType& _W)
    {
        mWX += _W*_X;
        mWY += _W*_Y;
        mWZ += _W*_Z;
        mW  += _W;
    }

    // overload operator []
    TDataType& operator[] (const int& i)
    {
        if (i == 0) return WX();
        else if (i == 1) return WY();
        else if (i == 2) return WZ();
        else if (i == 3) return W();
    }

    const TDataType& operator[] (const int& i) const
    {
        if (i == 0) return WX();
        else if (i == 1) return WY();
        else if (i == 2) return WZ();
        else if (i == 3) return W();
    }

    // overload operator ()
    TDataType operator() (const int& i) const
    {
        if (i == 0) return X();
        else if (i == 1) return Y();
        else if (i == 2) return Z();
        else if (i == 3) return W();
    }

    /// Assignment operator
    ControlPoint& operator=(const ControlPoint& rOther)
    {
        this->mWX = rOther.mWX;
        this->mWY = rOther.mWY;
        this->mWZ = rOther.mWZ;
        this->mW = rOther.mW;
        return *this;
    }

    /// Addition operator
    ControlPoint& operator+=(const ControlPoint& rOther)
    {
        this->mWX += rOther.mWX;
        this->mWY += rOther.mWY;
        this->mWZ += rOther.mWZ;
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
        this->mWX *= alpha;
        this->mWY *= alpha;
        this->mWZ *= alpha;
        this->mW *= alpha;
        return *this;
    }

    /// Apply the homogeneous transformation to the control point
    void ApplyTransformation(const Transformation<TDataType>& trans)
    {
        TDataType res[4];
        for (std::size_t i = 0; i < 4; ++i)
            res[i] = trans(i, 0)*this->mWX + trans(i, 1)*this->mWY + trans(i, 2)*this->mWZ + trans(i, 3)*this->mW;
        this->mWX = res[0]; this->mWY = res[1]; this->mWZ = res[2]; this->mW = res[3];
    }

    /// Multiplication operator
    friend ControlPoint operator*(const double& alpha, ControlPoint c)
    {
        c *= alpha;
        return c;
    }

    /// Multiplication operator
    friend ControlPoint operator*(const Transformation<TDataType>& trans, ControlPoint c)
    {
        c.ApplyTransformation(trans);
        return c;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Control Point";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        // print the control point in homogeneous coordinates
        rOStream << "(" << WX() << ", " << WY() << ", " << WZ() << ", " << W() << ")";
    }

private:
    TDataType mWX, mWY, mWZ, mW;
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


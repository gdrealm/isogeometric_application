#pragma once

#include "../geomodcore.h"

namespace geomodcore {
    
    static unsigned int NEXT_VERTEX_ID = 0;
    static unsigned int NEXT_EDGE_ID = 0;
    static unsigned int NEXT_FACE_ID = 0;

    typedef boost::tuple<PPoint2d, PPoint2d, GPoint>  ControlPoint;

}

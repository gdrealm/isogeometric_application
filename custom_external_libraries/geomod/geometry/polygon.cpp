
#include "polygon.h"

namespace geomodcore {

Polygon::Polygon(void) {
}

Polygon::Polygon(const std::vector<GPoint>& points) {
    this->assign(points.begin(), points.end());
}

Polygon::Polygon(const Polygon& curve) {
    this->assign(curve.begin(), curve.end());
}

Polygon& Polygon::operator= (const Polygon& curve) {
    this->assign(curve.begin(), curve.end());
    return *this;
}

Polygon::~Polygon(void) {
    clear();
}

double Polygon::getPolygonLength() const {
    double pLength = 0.0;
    for (unsigned int i = 1; i < size(); i++) {
        const GPoint& p0 = this->at(i-1);
        const GPoint& p1 = this->at(i);
        pLength += sqrt(pow(p1[0]-p0[0], 2) + pow(p1[1]-p0[1], 2) + pow(p1[2]-p0[2], 2));
    }
    return pLength;
}

}


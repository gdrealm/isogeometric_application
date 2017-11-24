#include "includes/define.h"
#include "custom_utilities/nurbs/knot_array_1d.h"

using namespace Kratos;

int main(int argc, char** argv)
{
    typedef KnotArray1D<double> knot_container_t;
    typedef knot_container_t::knot_t knot_t;

    knot_container_t knot_vector;

    knot_vector.pCreateKnot(0.0);
    knot_vector.pCreateKnot(0.0);
    knot_vector.pCreateKnot(0.0);
    knot_vector.pCreateKnot(0.0);
    knot_vector.pCreateKnot(1.0);
    knot_vector.pCreateKnot(1.0);
    knot_vector.pCreateKnot(1.0);
    knot_vector.pCreateKnot(1.0);

    KRATOS_WATCH(knot_vector)
    KRATOS_WATCH(knot_vector.nspans())
    std::tuple<knot_t, knot_t> span = knot_vector.span(1);
    KRATOS_WATCH(std::get<0>(span)->Index())
    KRATOS_WATCH(std::get<0>(span)->Value())
    KRATOS_WATCH(std::get<1>(span)->Index())
    KRATOS_WATCH(std::get<1>(span)->Value())

    return 0;
}


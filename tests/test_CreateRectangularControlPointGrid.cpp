#include "includes/define.h"
#include "custom_utilities/control_grid_library.h"

using namespace Kratos;

int main(int argc, char** argv)
{
    std::vector<double> start = {0.0, 0.0};

    std::vector<std::size_t> ngrid = {5, 5};

    std::vector<double> end = {1.0, 1.0};

    ControlGrid<ControlPoint<double> >::Pointer pGrid = ControlGridLibrary::CreateStructuredControlPointGrid<2>(start, ngrid, end);
    KRATOS_WATCH(*pGrid)

    return 0;
}


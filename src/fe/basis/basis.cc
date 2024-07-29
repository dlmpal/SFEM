#include "basis.h"
#include "../../common/error.h"

namespace sfem::fe::basis
{
    //=============================================================================
    Basis *CreateBasis(const mesh::Cell &cell)
    {
        Basis *basis = nullptr;

        switch (cell.type())
        {
        case mesh::CellType::point:
            basis = new basis::PointBasis();
            break;

        case mesh::CellType::line:
            switch (cell.order())
            {
            case 1:
                basis = new basis::LinearLineBasis();
                break;
            case 2:
                basis = new basis::QuadraticLineBasis();
                break;
            case 3:
                basis = new basis::CubicLineBasis();
                break;
            }
            break;

        case mesh::CellType::triangle:
            switch (cell.order())
            {
            case 1:
                basis = new basis::LinearTriangleBasis();
                break;
            case 2:
                basis = new basis::QuadraticTriangleBasis();
                break;
            case 3:
                basis = new basis::CubicTriangleBasis();
                break;
            }
            break;

        case mesh::CellType::quad:
            switch (cell.order())
            {
            case 1:
                basis = new basis::LinearQuadBasis();
                break;
            case 2:
                basis = new basis::QuadraticQuadBasis();
                break;
            case 3:
                basis = new basis::CubicQuadBasis();
                break;
            }
            break;

        case mesh::CellType::tet:
            switch (cell.order())
            {
            case 1:
                basis = new basis::LinearTetrahedralBasis();
                break;
            case 2:
                basis = new basis::QuadraticTehtrahedralBasis();
                break;
            case 3:
                basis = new basis::CubicTetrahedralBasis();
                break;
            }
            break;

        case mesh::CellType::hex:
            switch (cell.order())
            {
            case 1:
                basis = new basis::LinearHexahedralBasis();
                break;
            }
            break;

        default:
            break;
        }

        if (basis == nullptr)
        {
            error::invalid_cell_error(cell.idx(), static_cast<int>(cell.type()), cell.order(), __FILE__, __LINE__);
        }

        return basis;
    }
}
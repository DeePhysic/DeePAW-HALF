#include "ScaLapack.hpp"

#ifdef VASPML_USE_SCALAPACK
#include <thread>
#include <chrono>
#endif

using namespace vaspml;

namespace vaspml::ScaLapack
{

void vaspml_pdgemm([[maybe_unused]] transpose             transA,
                   [[maybe_unused]] transpose             transB,
                   [[maybe_unused]] ScaLapackArray<Real>& arrayA,
                   [[maybe_unused]] ScaLapackArray<Real>& arrayB,
                   [[maybe_unused]] ScaLapackArray<Real>& arrayC,
                   [[maybe_unused]] Real                  alpha,
                   [[maybe_unused]] Real                  beta)
{
#ifdef VASPML_USE_SCALAPACK
    Int nRowsTotA = arrayA.get_nRowsTot();
    Int nColsTotA = arrayA.get_nColsTot();
    Int nRowsTotB = arrayB.get_nRowsTot();
    Int nColsTotB = arrayB.get_nColsTot();

    Int startRowA = arrayA.get_startRow();
    Int startColA = arrayA.get_startCol();

    Int startRowB = arrayB.get_startRow();
    Int startColB = arrayB.get_startCol();

    Int startRowC = arrayC.get_startRow();
    Int startColC = arrayC.get_startCol();

    Int  ctxtA = arrayA.getContext();
    char transALoc;
    char transBLoc;
    if (transA == transpose::NoTrans) transALoc = 'N';
    else transALoc = 'T';
    if (transB == transpose::NoTrans) transBLoc = 'N';
    else transBLoc = 'T';

    std::cout << alpha << "  " << beta << std::endl;

    // op(A) is MxK = nRowsTotA x nColsTotA
    // op(B) is KxN = nColsTotA x nColsTotB
    // op(C) is MxN = nRowsTotA x nColsTotB
    //std::cout << "desc_A      " << arrayA.get_rank() << "   ";
    //for ( Int i = 0; i < 9; i++ ){
    //    std::cout << arrayA.get_descriptor()[i] << "   ";
    //}
    //std::cout << std::endl;
    //std::cout << "desc_B      " << arrayB.get_rank() << "   ";
    //for ( Int i = 0; i < 9; i++ ){
    //    std::cout << arrayB.get_descriptor()[i] << "   ";
    //}
    //std::cout << std::endl;
    std::cout << "entry desc_A desc_B << desc_C" << std::endl;
    for (Int i = 0; i < 9; i++)
    {
        std::cout << i << " " << arrayA.get_descriptor()[i] << " " << arrayB.get_descriptor()[i]
                  << " " << arrayC.get_descriptor()[i] << std::endl;
    }

    // this runs trough but gives wrong results because of column major
    //   transALoc = 'N';
    //   transBLoc = 'N';
    //
    //   Int rA = nRowsTotA;
    //   Int rB = nRowsTotB;
    //   Int cB = nColsTotB;
    //
    //std::cout << rA << "  " << rB << "    " << cB << std::endl;
    //pdgemm_( &transALoc, &transBLoc,
    //         &rA,
    //         &cB,
    //         &rB,
    //         &alpha,

    //         arrayA.get_data().data(),
    //         &startRowA, &startColA,
    //         arrayA.get_descriptor(),

    //         arrayB.get_data().data(),
    //         &startRowB, &startColB,
    //         arrayB.get_descriptor(),

    //         &beta,
    //         arrayC.get_data().data(),
    //         &startRowC, &startColC,
    //         arrayC.get_descriptor() ); // describe arrayC

    transALoc = 'T';
    transBLoc = 'T';

    Int n1 = nRowsTotA;
    Int n2 = nRowsTotB;
    Int n3 = nColsTotB;
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    std::cout << "matrix A " << arrayA.get_rank() << "   " << arrayA.get_nRowsLoc() << "  "
              << arrayA.get_nColsLoc() << "   " << arrayA.get_nRowsTot() << "  "
              << arrayA.get_nColsTot() << "   " << std::endl;
    std::cout << "matrix B " << arrayB.get_rank() << "   " << arrayB.get_nRowsLoc() << "  "
              << arrayB.get_nColsLoc() << "   " << arrayB.get_nRowsTot() << "  "
              << arrayB.get_nColsTot() << "   " << std::endl;
    std::cout << "matrix C " << arrayC.get_rank() << "   " << arrayC.get_nRowsLoc() << "  "
              << arrayC.get_nColsLoc() << "   " << arrayC.get_nRowsTot() << "  "
              << arrayC.get_nColsTot() << "   " << std::endl;

    pdgemm_(&transALoc,
            &transBLoc,
            &n3,
            &n2,
            &n1,
            &alpha,

            arrayA.get_data().data(),
            &startRowA,
            &startColA,
            arrayA.get_descriptor(),

            arrayB.get_data().data(),
            &startRowB,
            &startColB,
            arrayB.get_descriptor(),

            &beta,
            arrayC.get_data().data(),
            &startRowC,
            &startColC,
            arrayC.get_descriptor()); // describe arrayC

#else
    // put your cblas call here
    // NOTE (aS): If you add the implementation here, please remove the [[maybe_unused]] attributes
    // at the top in the function signature. I put them there just to temporarily silence compiler
    // warnings.
#endif

    return;
}

} //namespace vaspml::ScaLapack

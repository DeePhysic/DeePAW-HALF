#ifndef SCALAPACK_HPP
#define SCALAPACK_HPP

#ifdef VASPML_USE_SCALAPACK
#include "MlMPI.hpp"
#endif

#include "Tutor.hpp"
#include "debug.hpp"
#include "scalapackInterface.hpp"
#include "types.hpp"

#ifdef VASPML_USE_SCALAPACK
#include "math.hpp"

#include <algorithm> // for max
#endif

#include <iostream>

namespace vaspml
{

#ifndef VASPML_USE_SCALAPACK
class MlMPI;
#endif

using namespace scalapack_interface;

/*******************************************************************************************
 * Array class giving the possibility to distribute the memory with scalapack.
 *
 * The ScaLapackArray can be compiled with the or without scalapack. In case it is compiled
 * without scalapack the data will be stored in an ordinary std::vector
 *
 *
 * The class can be used in the following way:
 *
 * @code
 * ScaLapackArray array;
 * array.Init( nRowsTot, nColsTot, mpiIn, nProcessesRow, nProcessesCol );
 * @endcode
 *
 *******************************************************************************************/
template<class T>
class ScaLapackArray
{
  public:
    /*******************************************************************************************
     * default constructor, initializes all parameters to zero. Supply template
     * parameter when calling
     *
     * after calling class is not in a usable state. Only in a created state
     *******************************************************************************************/
    ScaLapackArray();
    /*******************************************************************************************
     * creating a ScaLapackArray class.
     *
     * after calling the local arrays are not allocated yet. The user has to call the
     * allocateArray function to allocate the local arrays. This gives the freedom
     * to set up the grid and context and only allocate array when needed.
     *
     * @param[in] nRowsTot total number of rows in non distributed matrix
     * @param[in] nColsTot total number of columns in non distributed matrix
     * @param[in] mpiIn mpi communicator from which the scalapack context will be created
     * @param[in] nProcessesRow number of processors for rows
     * @param[in] nProcessesCol number of processors for columns
     *
     * When mpi is not supplied a new scalapack context will be created.
     * If nProcessesRow and nProcessesCol are not supplied or are supplied in a faulty
     * way a Fermat factorization is used to determine these values.
     *******************************************************************************************/
    ScaLapackArray(const Int    nRowsTot,
                   const Int    nColsTot,
                   const MlMPI* mpiIn = nullptr,
                   const Int    nProcessesRow = 0,
                   const Int    nProcessesCol = 0);
    /*******************************************************************************************
     * Initializing class. Can be called after default constructor
     *
     * @param[in] nRowsTot total number of rows in non distributed matrix
     * @param[in] nColsTot total number of columns in non distributed matrix
     * @param[in] mpiIn mpi communicator from which the scalapack context will be created
     * @param[in] nProcessesRow number of processors for rows
     * @param[in] nProcessesCol number of processors for columns
     *
     * When mpi is not supplied a new scalapack context will be created.
     * If nProcessesRow and nProcessesCol are not supplied or are supplied in a faulty
     * way a Fermat factorization is used to determine these values.
     *******************************************************************************************/
    void init(const Int    nRowsTot,
              const Int    nColsTot,
              const MlMPI* mpiIn = nullptr,
              const Int    nProcessesRow = 0,
              const Int    nProcessesCol = 0);
    /*******************************************************************************************
     * function allocates the data array. After this function was called the class can be fully
     * used
     *******************************************************************************************/
    void allocateArray();
    /*******************************************************************************************
     * returns the index in the distributed scalapack array
     *
     * @param[in] rowIndex supply row index of matrix. First index in math notation
     * @param[in] colIndex supply row index of matrix. First second in math notation
     *
     * second index is fast index
     *******************************************************************************************/
    Int getLocalIndex(const Int rowIndex, const Int colIndex) const;
    /*******************************************************************************************
     * returns value stored in distributed matrix
     *
     * @param[in] rowIndex supply row index of matrix. First index in math notation
     * @param[in] colIndex supply row index of matrix. First second in math notation
     *
     * second index is fast index
     *******************************************************************************************/
    T getValue(const Int rowIndex, const Int colIndex) const;
    /*******************************************************************************************
     * returns the index in the distributed scalapack array
     *
     * @param[in] rowIndex supply row index of matrix. First index in math notation
     * @param[in] colIndex supply row index of matrix. First second in math notation
     * @param[in] value defines the value to which the array entry will be set
     *
     * second index is fast index
     *******************************************************************************************/
    void setValue(const Int rowIndex, const Int colIndex, const T value);
    /*******************************************************************************************
     * return the number of rows in the distributed array for the calling processor
     *******************************************************************************************/
    Int get_nRowsLoc() const;
    /*******************************************************************************************
     * return the number of columns in the distributed array for the calling processor
     *******************************************************************************************/
    Int get_nColsLoc() const;
    /*******************************************************************************************
     * return the total number of rows in the non distributed total array;
     *
     * this quantity is independent of the processor calling the getter
     *******************************************************************************************/
    Int get_nRowsTot() const;
    /*******************************************************************************************
     * return the total number of columns in the non distributed total array;
     *
     * this quantity is independent of the processor calling the getter
     *******************************************************************************************/
    Int get_nColsTot() const;
    /*******************************************************************************************
     * get the total number of processors over which the array is distributed
     *******************************************************************************************/
    Int get_nprocs() const;
    /*******************************************************************************************
     *  get rank of processor
     *******************************************************************************************/
    Int get_rank() const;
    /*******************************************************************************************
     * get row index of calling processor
     *******************************************************************************************/
    Int get_procRowIndex() const;
    /*******************************************************************************************
     * get column index of calling processor
     *******************************************************************************************/
    Int get_procColIndex() const;
    /*******************************************************************************************
     * get number of processor per row; independent of calling processors
     *
     * gives the number of row grid points
     *******************************************************************************************/
    Int get_nProcessesRow() const;
    /*******************************************************************************************
     * get number of processor per column; independent of calling processors
     *
     * gives the number of column grid points
     *******************************************************************************************/
    Int get_nProcessesCol() const;
    /*******************************************************************************************
     * get std::vector<T> in which the scalapack distributed data is stored
     *******************************************************************************************/
    std::vector<T>& get_data();
    /*******************************************************************************************
     * get descriptor on how the array is distributed in scalapack
     *******************************************************************************************/
    Int* get_descriptor();
    /*******************************************************************************************
     * get the context of the given scalapack array
     *******************************************************************************************/
    Int getContext() const;
    /*******************************************************************************************
     * get start index of row in total non distributed matrix
     *******************************************************************************************/
    Int get_startRow() const;
    /*******************************************************************************************
     * get start index of column in total non distributed matrix
     *******************************************************************************************/
    Int get_startCol() const;
    /*******************************************************************************************
     * write distributed matrix for given rank to screen
     *
     * @param[in] rank for which output is written
     *******************************************************************************************/
    void writeToScreen(const Int rank, const String header = "") const;

  private:
    /*******************************************************************************************
     *  do initializtion when vaspml is compiled without scalapack
     *******************************************************************************************/
    void NoScalapackInit(const Int nRowsTot,
                         const Int nColsTot,
                         const MlMPI* = nullptr,
                         const Int = 0,
                         const Int = 0);
    /// total number of processors
    Int nprocs;
    /// private rank in communicator
    Int rank;
    /// row index of current processor
    Int procRowIndex;
    /// column index of current processor
    Int procColIndex;
    /// context ( communicator ) of scalapack
    Int ctxt;
    /// number of processes acting along a row
    Int nProcessesRow;
    /// number of processes acting along a column
    Int nProcessesCol;
    /// number of rows in the total matrix
    Int nRowsTot;
    /// number of columns in the total matrix
    Int nColsTot;
    /// number of rows for given rank
    Int nRowsLoc;
    /// number of coulmns for given rank
    Int nColsLoc;
    /// leading dimension of the local matrix
    Int lld;
    /// data array storing the local data for every processor
    std::vector<T> data;
    /*******************************************************************************************
     * blocking factor for distribution of matrices
     * P4 optimal, larger values were even slightly better (160) but still
     *
     * @note can not be defined as const because of the scalapack interface
     *
     *******************************************************************************************/
    //Int blockingFactor = 16;
    Int blockingFactor = 2;
    Int zero = 0;
    Int one = 1;
    /*******************************************************************************************
     * scalapack descriptor
     *
     * for more information check the documentation for the function descinit_
     *******************************************************************************************/
    Int descriptor[9];
    /// control variable to check if data array was already allocated
    bool isAllocated;
    Int  startRow;
    Int  startCol;
};

template<class T>
ScaLapackArray<T>::ScaLapackArray()
{
    nprocs = 0;
    rank = 0;
    ctxt = -1;
    nProcessesRow = 0;
    nProcessesCol = 0;
    nRowsTot = 0;
    nColsTot = 0;
    nRowsLoc = 0;
    nColsLoc = 0;
    procRowIndex = 0;
    procColIndex = 0;
    isAllocated = false;
    startRow = 1;
    startCol = 1;
}

template<typename T>
ScaLapackArray<T>::ScaLapackArray(const Int    nRowsTot,
                                  const Int    nColsTot,
                                  const MlMPI* mpiIn,
                                  const Int    nProcessesRow,
                                  const Int    nProcessesCol)
{
    init(nRowsTot, nColsTot, mpiIn, nProcessesRow, nProcessesCol);
}

template<typename T>
void ScaLapackArray<T>::NoScalapackInit(const Int nRowsTot,
                                        const Int nColsTot,
                                        const MlMPI*,
                                        const Int,
                                        const Int)
{
    this->nRowsTot = nRowsTot;
    this->nColsTot = nColsTot;
    nRowsLoc = nRowsTot;
    nColsLoc = nColsTot;
    procRowIndex = 1;
    procColIndex = 1;
    isAllocated = false;

    return;
}

template<class T>
void ScaLapackArray<T>::init(const Int                     nRowsTot,
                             const Int                     nColsTot,
                             [[maybe_unused]] const MlMPI* mpiIn,
                             [[maybe_unused]] const Int    nProcessesRow,
                             [[maybe_unused]] const Int    nProcessesCol)
{
    String funcName = __func__;
    if (nRowsTot <= 0 or nColsTot <= 0)
    {
        String errorMessage = "ERROR: " + funcName + " invalid values of nRowsTot "
                            + std::to_string(nRowsTot) + " or nColsTot " + std::to_string(nColsTot)
                            + "\n";
        global::tutor.error(errorMessage);
    }
    this->nRowsTot = nRowsTot;
    this->nColsTot = nColsTot;
#ifdef VASPML_USE_SCALAPACK
    if (mpiIn->get_communicator() == MPI_COMM_NULL)
    {
        Cblacs_pinfo(&rank, &nprocs);
        Cblacs_get(0, 0, &ctxt);
    }
    else
    {
        // couple to MPI communicator
        rank = mpiIn->get_rank();
        nprocs = mpiIn->get_numberRanks();
        ctxt = Csys2blacs_handle(mpiIn->get_communicator());
    }
    // determine the number of processor rows and columns
    // when input is not appropriate
    if (nProcessesRow <= 0 or nProcessesCol <= 0)
    {
        math::fermatRazor(nprocs, this->nProcessesRow, this->nProcessesCol);
    }
    else if (nProcessesRow * nProcessesCol != nprocs)
    {
        String warningMessage = "WARNING: " + funcName + " invalid values of nProcessesRow "
                              + std::to_string(nProcessesRow) + " nProcessesCol "
                              + std::to_string(nProcessesRow)
                              + +" determining new values with fermatRazor \n";
        global::tutor.warning(warningMessage);
        math::fermatRazor<Int>(nprocs, this->nProcessesRow, this->nProcessesCol);
    }
    else
    {
        this->nProcessesRow = nProcessesRow;
        this->nProcessesCol = nProcessesCol;
    }
    Cblacs_gridinit(&ctxt, "Row-major", this->nProcessesRow, this->nProcessesCol);
    Cblacs_gridinfo(ctxt, &this->nProcessesRow, &this->nProcessesCol, &procRowIndex, &procColIndex);

    // determine number of elements per rank
    // nRowsLoc = numroc_( &this->nRowsTot, &blockingFactor, &procRowIndex, &zero,
    // &this->nProcessesRow ); nColsLoc = numroc_( &this->nColsTot, &blockingFactor, &procColIndex,
    // &zero, &this->nProcessesCol );

    nColsLoc =
        numroc_(&this->nRowsTot, &blockingFactor, &procRowIndex, &zero, &this->nProcessesRow);
    nRowsLoc =
        numroc_(&this->nColsTot, &blockingFactor, &procColIndex, &zero, &this->nProcessesCol);

    // determine number of elements per rank
    // determine leading dimension of local ( distributed ) matrix
    //lld = std::max( 1, numroc_( &this->nRowsTot, &blockingFactor, &procRowIndex, &zero,
    //&this->nProcessesRow ) );
    lld = std::max(
        1,
        numroc_(&this->nRowsTot, &blockingFactor, &procRowIndex, &zero, &this->nProcessesRow));
    Int info = 0;
    // couple mpi communicator to scalapack context
    if (mpiIn->get_communicator() != MPI_COMM_NULL)
        ctxt = Csys2blacs_handle(mpiIn->get_communicator());
    //descinit_( descriptor, &nRowsTot, &nColsTot, &blockingFactor, &blockingFactor, &zero, &zero,
    //&ctxt, &lld, &info );
    descinit_(descriptor,
              &nRowsTot,
              &nColsTot,
              &blockingFactor,
              &blockingFactor,
              &zero,
              &zero,
              &ctxt,
              &lld,
              &info);

    if (info != 0)
    {
        String warningMessage = "WARNING: " + funcName + " descinit_ returned error code "
                              + std::to_string(info) + "\n";
        global::tutor.warning(warningMessage);
    }
    isAllocated = false;
    startRow = 1;
    startCol = 1;
#else
    // initialization without scalapack memory distribution
    NoScalapackInit(nRowsTot, nColsTot);
#endif

    return;
}

template<class T>
void ScaLapackArray<T>::allocateArray()
{
    data.resize(nRowsLoc * nColsLoc);
    isAllocated = true;

    return;
}

template<class T>
Int ScaLapackArray<T>::getLocalIndex(const Int rowIndex, const Int colIndex) const
{
    return rowIndex * nColsLoc + colIndex;
}

template<class T>
void ScaLapackArray<T>::setValue(const Int rowIndex, const Int colIndex, const T value)
{
    VASPML_DEBUG_L2(
        if (!isAllocated)
        {
            String function = __func__;
            String errorMessage = "ERROR: " + function
                                + " trying to set value but data array not allocated yet.\n"
                                + "Call allocateArray( void ) first.";
            global::tutor.bug(errorMessage);
        }
    );
    data[rowIndex * nColsLoc + colIndex] = value;

    return;
}

template<class T>
T ScaLapackArray<T>::getValue(const Int rowIndex, const Int colIndex) const
{
    VASPML_DEBUG_L2(
        if (!isAllocated)
        {
            String function = __func__;
            String errorMessage = "ERROR: " + function
                                + " trying to set value but data array not allocated yet.\n"
                                + "Call allocateArray( void ) first. And then fill matrix";
            global::tutor.bug(errorMessage);
        }
    );

    return data[rowIndex * nColsLoc + colIndex];
}

template<class T>
Int ScaLapackArray<T>::get_nRowsLoc() const
{
    return nRowsLoc;
}

template<class T>
Int ScaLapackArray<T>::get_nColsLoc() const
{
    return nColsLoc;
}

template<class T>
Int ScaLapackArray<T>::get_nRowsTot() const
{
    return nRowsTot;
}

template<class T>
Int ScaLapackArray<T>::get_nColsTot() const
{
    return nColsTot;
}

template<typename T>
Int ScaLapackArray<T>::get_rank() const
{
    return rank;
}

template<typename T>
Int ScaLapackArray<T>::get_procRowIndex() const
{
    return procRowIndex;
}

template<typename T>
Int ScaLapackArray<T>::get_procColIndex() const
{
    return procColIndex;
}

template<typename T>
Int ScaLapackArray<T>::get_nProcessesRow() const
{
    return nProcessesRow;
}

template<typename T>
Int ScaLapackArray<T>::get_nProcessesCol() const
{
    return nProcessesCol;
}

template<typename T>
std::vector<T>& ScaLapackArray<T>::get_data()
{
    return data;
}

template<typename T>
Int* ScaLapackArray<T>::get_descriptor()
{
    return descriptor;
}

template<typename T>
Int ScaLapackArray<T>::getContext() const
{
    return ctxt;
}

template<typename T>
Int ScaLapackArray<T>::get_startRow() const
{
    return startRow;
}

template<typename T>
Int ScaLapackArray<T>::get_startCol() const
{
    return startCol;
}

template<typename T>
Int ScaLapackArray<T>::get_nprocs() const
{
    return nprocs;
}

template<typename T>
void ScaLapackArray<T>::writeToScreen(const Int rank, const String header) const
{
    if (this->rank == rank)
    {
        if (!header.empty()) { std::cout << header << std::endl; }
        for (Int i = 0; i < nRowsLoc; i++)
        {
            for (Int j = 0; j < nColsLoc; j++) { std::cout << getValue(i, j) << "   "; }
            std::cout << std::endl;
        }
    }

    return;
}

/*******************************************************************************************
 * namespace contains a collection of functionality which can be applied to the
 * ScaLapackArray class
 *******************************************************************************************/
namespace ScaLapack
{

enum class transpose
{
    Trans,
    NoTrans
};

template<typename T>
void vaspml_pdgemr2d(ScaLapackArray<T>& arrayA, ScaLapackArray<T>& arrayB)
{
    VASPML_DEBUG_L1(
        if (arrayA.get_nRowsTot() != arrayB.get_nRowsTot()
            or arrayA.get_nColsTot() != arrayB.get_nColsTot())
        {
            String functionName = __func__;
            global::tutor.bug("ERROR: " + functionName
                              + " dimensions of scalapack arrays do not agree");
        }
    );

    Int nRowsTot = arrayA.get_nRowsTot();
    Int nColsTot = arrayA.get_nColsTot();
    Int startRowA = arrayA.get_startRow();
    Int startColA = arrayA.get_startCol();
    Int startRowB = arrayB.get_startRow();
    Int startColB = arrayB.get_startCol();
    Int ctxtA = arrayA.getContext();

    pdgemr2d_(&nRowsTot,
              &nColsTot,
              arrayA.get_data().data(),
              &startRowA,
              &startColA,
              arrayA.get_descriptor(),
              arrayB.get_data().data(),
              &startRowB,
              &startColB,
              arrayB.get_descriptor(),
              &ctxtA);

    return;
}

void vaspml_pdgemm(transpose             transA,
                   transpose             transB,
                   ScaLapackArray<Real>& arrayA,
                   ScaLapackArray<Real>& arrayB,
                   ScaLapackArray<Real>& arrayC,
                   Real                  alpha,
                   Real                  beta);

} //namespace ScaLapack

} //namespace vaspml

#endif // SCALAPACK_HPP

#ifndef SCALAPACKINTERFACE_HPP
#define SCALAPACKINTERFACE_HPP

#ifdef VASPML_USE_SCALAPACK
#include "mpi.hpp"

#include "types.hpp"
#endif

namespace vaspml::scalapack_interface
{

#ifdef VASPML_USE_SCALAPACK
extern "C"
{
/*******************************************************************************************
 * process info routine sets up the scalapack process number and number of processes
 *
 * @param rank scalapack rank of processor as in MPI
 * @param nprocs scalapack number of processors as in MPI
 *******************************************************************************************/
void Cblacs_pinfo(Int* rank, Int* nprocs);
/*******************************************************************************************
 * Gets values that BLACS use for internal defaults.
 *
 * This routine gets the values that the BLACS are using for internal defaults. Some values are
 * tied to a BLACS context, and some are more general. The most common use is in retrieving a
 *default system context for input into blacs_gridinit or blacs_gridmap. Some systems, such as MPI*,
 *supply their own version of context. For those users who mix system code with BLACS code, a BLACS
 *context should be formed in reference to a system context. Thus, the grid creation routines take a
 *system context as input. If you wish to have strictly portable code, you may use blacs_get to
 *retrieve a default system context that will include all available processes. This value is not
 *tied to a BLACS context, so the parameter icontxt is unused. blacs_get returns information on
 *three quantities that are tied to an individual BLACS context, which is passed in as icontxt. The
 *information that may be retrieved is:
 *
 * @param icontxt On values of what that are tied to a particular context, this parameter is the
 * integer handle indicating the context. Otherwise, ignored.
 * @param what Indicates what BLACS internal(s) should be returned in val. Present options are:
 * <table>
 * <caption id="multi_row">values of what variable</caption>
 * <tr><th> what = 0   <th> Handle indicating default system context. <th>
 * <tr><th> what = 1   <th> The BLACS message ID range. <th>
 * <tr><th> what = 2   <th> The BLACS debug level the library was compiled with. <th>
 * <tr><th> what = 10  <th>  Handle indicating the system context used to define the BLACS context
 *                           whose handle is icontxt. <th>
 * <tr><th> what = 11  <th>  Number of rings multiring broadcast topology is presently using. <th>
 * <tr><th> what = 12  <th>  Number of branches general tree broadcast topology is presently using.
 *<th> <tr><th> what = 13  <th>  Number of rings multiring combine topology is presently using. <th>
 * <tr><th> what = 14  <th>  Number of branches general tree combine topology is presently using.
 *<th> <tr><th> what = 15  <th>  Whether topologies are forced to be repeatable or not. A non-zero
 *return value indicates that topologies are being forced to be repeatable. See Repeatability and
 *Coherence for more information about repeatability.<th> <tr><th> what = 16  <th>  Whether
 *topologies are forced to be heterogenous coherent or not. A non-zero return value indicates that
 *topologies are being forced to be heterogenous coherent. See Repeatability and Coherence for more
 *information about coherence.<th>
 * </table>
 * @param ctxt can be considered as a MPI communicator. It defines the contect on which
 * the remaining blacs and scalapack calls will operate. Usually every distributed array
 * should get it's own context.
 *******************************************************************************************/
void Cblacs_get(Int icontxt, Int what, Int* val);
/*******************************************************************************************
 * Assigns available processes into BLACS process grid.
 *
 * All BLACS codes must call this routine, or its sister routine blacs_gridmap. These routines take
 * the available processes, and assign, or map, them into a BLACS process grid.
 * In other words, they establish how the BLACS coordinate system maps into the native machine
 *process numbering system. Each BLACS grid is contained in a context, so that it does not interfere
 *with distributed operations that occur within other grids/contexts. These grid creation routines
 *may be called repeatedly to define additional contexts/grids. The creation of a grid requires
 *input from all processes that are defined to be in this grid. Processes belonging to more than one
 *grid have to agree on which grid formation will be serviced first, much like the globally blocking
 *sum or broadcast. These grid creation routines set up various internals for the BLACS, and one of
 *them must be called before any calls are made to the non-initialization BLACS. Note that these
 *routines map already existing processes to a grid: the processes are not created dynamically. On
 *most parallel machines, the processes are "created" when you run your executable. When using the
 *PVM BLACS, if the virtual machine has not been set up yet, the routine blacs_setup should be used
 *to create the virtual machine. This routine creates a simple nprow x npcol process grid. This
 *process grid uses the first nprow * npcol processes, and assigns them to the grid in a row- or
 *column-major natural ordering. If these process-to-grid mappings are unacceptable, call
 *blacs_gridmap.
 *
 * @param ictxt handle indicating the system context to be used in creating the BLACS context.
 * Call blacs_get to obtain a default system context. Can be thought of as the communicator created
 *by Cblacs_get
 * @param order Indicates how to map processes to BLACS grid. Options are 'R' for Use row-major
 *natural ordering 'C' for column-major natural ordering
 * @param nprow Indicates how many process rows the process grid should contain.
 * @param npcol Indicates how many process columns the process grid should contain.
 *******************************************************************************************/
void Cblacs_gridinit(Int* ictxt, const char* order, Int nprow, Int npcol);
/*******************************************************************************************
 * Returns information on the current grid.
 *
 * This routine returns information on the current grid. If the context handle does not point at a
 *valid context, all quantities are returned as -1
 *
 * @param ictxt Integer handle that indicates the context.
 * @param nprow Number of process rows in the current process grid.
 * @param Number of process columns in the current process grid.
 * @param myrow Row coordinate of the calling process in the process grid.
 * @param mycol Column coordinate of the calling process in the process grid.
 *******************************************************************************************/
void Cblacs_gridinfo(Int ictxt, Int* nprow, Int* npcol, Int* myrow, Int* mycol);
/*******************************************************************************************
 * Returns the row and column coordinates in the process grid.
 *
 * Given the system process number, this function returns the row and column coordinates in
 * the BLACS process grid.
 *
 * @param ictxt handle that indicates the context.
 * @param rank Process number the coordinates of which are to be determined. This parameter
 * stand for the process number of the underlying machine, that is, it is a tid for PVM.
 * @param myrow Row coordinates of the pnum process in the BLACS grid.
 * @param mycol Column coordinates of the pnum process in the BLACS grid.
 *******************************************************************************************/
void Cblacs_pcoord(Int ictxt, Int rank, Int* myrow, Int* mycol);
/*******************************************************************************************
 * This function computes the local number of rows or columns of a block-cyclically distributed
 *matrix contained in a process row or process column, respectively, indicated by the calling
 *sequence argument iproc.
 *
 * @param n is the number of rows M_ or columns N_ in a global matrix that has been block-cyclically
 *distributed
 * @param nb is the row block size MB_ or the column block size NB_
 * @param is process row index myrow or the process column index mycol.
 * @param is the process row RSRC_ or the process column CSRC_ over which the first row or column,
 *respectively, of the global matrix is distributed.
 * @param is the number of rows nprow or the number of columns npcol in the process grid.
 *
 * @return is the local number of rows or columns of a block-cyclically distributed matrix
 * contained in a process row or process column, respectively, indicated by the calling sequence
 *argument iproc
 *******************************************************************************************/
Int numroc_(Int* n, Int* nb, Int* iproc, Int* isrcproc, Int* nprocs);
/*******************************************************************************************
 * This subroutine initializes a type-1 array descriptor with error checking.
 *
 * @param desc
 * @param m is the number of rows M_ in a global matrix that has been block-cyclically distributed.
 * @param is the number of columns N_ in a global matrix that has been block-cyclically distributed.
 * @param mb is the row block size MB_.
 * @param nb is the column block size NB_.
 * @param irsrc is the process row RSRC_ over which the first row of the global matrix is
 *distributed.
 * @param icsrc is the process column CSRC_ over which the first column of the global matrix is
 *distributed.
 * @param ictxt is the BLACS context.
 * @param lld is the leading dimension of the local array.
 *******************************************************************************************/
void descinit_(Int*       desc,
               const Int* m,
               const Int* n,
               const Int* mb,
               const Int* nb,
               const Int* irsrc,
               const Int* icsrc,
               const Int* ictxt,
               const Int* lld,
               Int*       info);
/*******************************************************************************************
 * Holds up execution of all processes within the indicated scope.
 *
 * This routine holds up execution of all processes within the indicated scope until they have all
 *called the routine.
 *
 * @param ictxt handle that indicates the context.
 * @param type Parameter that indicates whether a process row (scope='R'), column ('C'), or
 * entire grid ('A') will participate in the barrier.
 *******************************************************************************************/
void Cblacs_barrier(Int ictxt, char* type);
/*******************************************************************************************
 * Copies a submatrix from one general rectangular matrix to another.
 *
 * The p?gemr2dfunction copies the indicated matrix or submatrix of A to the indicated matrix or
 *submatrix of B. It provides a truly general copy from any block cyclicly-distributed matrix or
 *submatrix to any other block cyclicly-distributed matrix or submatrix. With p?trmr2d, these
 *functions are the only ones in the ScaLAPACK library which provide inter-context operations: they
 *can take a matrix or submatrix A in context A (distributed over process grid A) and copy it to a
 *matrix or submatrix B in context B (distributed over process grid B).
 *
 * There does not need to be a relationship between the two operand matrices or submatrices other
 * than their global size and the fact that they are both legal block cyclicly-distributed matrices
 * or submatrices. This means that they can, for example, be distributed across different
 * process grids, have varying block sizes and differing matrix starting points, or be
 * contained in different sized distributed matrices.
 *
 * @param m global number of rows of matrix A
 * @param n global number of columns of matrix A
 * @param A Pointer to array of size lld_a* LOCc(ja+n-1) containing the source matrix A.
 * @param ia,ja The row and column index in the array A indicating the
 * first row and the first column, respectively, of the submatrix of A) to copy.
 * @param desca (global and local) array of size dlen_. The array descriptor for the distributed
 *matrix A. Only dtype_a = 1 is supported, so dlen_ = 9. If the calling process is not part of the
 *context of A, ctxt_a must be equal to -1.
 * @param ib, jb (global) The row and column indices in the array B indicating the first row and the
 * first column, respectively, of the submatrix B to which to copy the matrix.
 * 1 ≤ib≤total_rows_in_b - m +1, 1 ≤jb≤total_columns_in_b - n +1
 * @param descb (global and local) array of size dlen_. The array descriptor for the distributed
 *matrix B. Only dtype_b = 1 is supported, so dlen_ = 9. If the calling process is not part of the
 *context of B, ctxt_b must be equal to -1.
 * @param ictxt The context encompassing at least the union of all processes in context A and
 *context B. All processes in the context ictxt must call this function, even if they do not own a
 *piece of either matrix.
 * @param B Pointer into the local memory to array of size lld_b*LOCc(jb+n-1).
 * Overwritten by the submatrix from A.
 *******************************************************************************************/
void pdgemr2d_(const Int*  m,
               const Int*  n,
               const Real* A,
               const Int*  ia,
               const Int*  ja,
               const Int*  desca,
               Real*       B,
               const Int*  ib,
               const Int*  jb,
               const Int*  descb,
               const Int*  ictxt);
/*******************************************************************************************
 * Frees a BLACS context
 *
 * Release the resources when contexts are no longer needed. After freeing a context, the context no
 *longer exists, and its handle may be re-used if new contexts are defined.
 *
 * @param ictxt handle that indicates the BLACS context to be freed.
 *******************************************************************************************/
void Cblacs_gridexit(Int ictxt);
/*******************************************************************************************
 * Frees all BLACS contexts and releases all allocated memory.
 *
 * This routine should be called when a process has finished all use of the BLACS. The
 * continue parameter indicates whether the user will be using the underlying communication
 * platform after the BLACS are finished. This information is most important for the PVM BLACS.
 * If continue is set to 0, then pvm_exit is called; otherwise, it is not called.
 * Setting continue not equal to 0 indicates that explicit PVM send/recvs will
 * be called after the BLACS routines are used. Make sure your code calls pvm_exit.
 * PVM users should either call blacs_exit or explicitly call pvm_exit to avoid PVM problems.
 *
 * @param Flag indicating whether message passing continues after the BLACS are done. If continue
 * is non-zero, the user is assumed to continue using the machine after completing the BLACS.
 * Otherwise, no message passing is assumed after calling this routine.
 *
 * @note not compatible with MPI_Finalize()
 *******************************************************************************************/
void Cblacs_exit(Int* error_code);
/*******************************************************************************************
 * convert MPI communicator to an Int.
 *
 * This function can be used to be the blacs/scalapack context to a certain MPI communicator
 *
 * @param[in] communicator MPI communicator to use in combination with scalapack
 *******************************************************************************************/
Int Csys2blacs_handle(MPI_Comm communicator);

/*******************************************************************************************
 * Computes a scalar-matrix-matrix product and adds the result to a scalar-matrix product
 * for distributed matrices.
 *
 * where:
 * op(x) is one of op(x) = x, or op(x) = x',
 * alpha and beta are scalars,
 * sub(A)=A(ia:ia+m-1, ja:ja+k-1), sub(B)=B(ib:ib+k-1, jb:jb+n-1), and sub(C)=C(ic:ic+m-1,
 *jc:jc+n-1), are distributed matrices.
 *
 * @param[in] transA Specifies the form of op(sub(A)) used in the matrix multiplication:
 * if transa = 'N' or 'n', then op(sub(A)) = sub(A);
 * if transa = 'T' or 't', then op(sub(A)) = sub(A)';
 * if transa = 'C' or 'c', then op(sub(A)) = sub(A)'.
 * @param[in] transB Specifies the form of op(sub(B)) used in the matrix multiplication:
 * if transb = 'N' or 'n', then op(sub(B)) = sub(B);
 * if transb = 'T' or 't', then op(sub(B)) = sub(B)';
 * if transb = 'C' or 'c', then op(sub(B)) = sub(B)'.
 * @param[in] N1 Specifies the number of rows of the distributed matrices op(sub(A)) and sub(C), m≥
 *0.
 * @param[in] N2 Specifies the number of columns of the distributed matrices op(sub(B)) and sub(C),
 *n≥ 0.
 * @param[in] N3 Specifies the number of columns of the distributed matrix op(sub(A)) and the number
 *of rows of the distributed matrix op(sub(B)).
 * @param[in] alpha When alpha is equal to zero, then the local entries of the arrays matrixA and
 *matrixB corresponding to the entries of the submatrices sub(A) and sub(B) respectively need not be
 *set on input.
 * @param[in] matrixA Array, size lld_a by kla, where kla is LOCc(ja+k-1) when transa = 'N' or 'n',
 *and is LOCq(ja+m-1) otherwise. Before entry this array must contain the local pieces of the
 *distributed matrix sub(A).
 * @param[in] ia row index in the distributed matrix A indicating the first row of the submatrix
 *sub(A)
 * @param[in] ja column index in the distributed matrix A indicating the first column of the
 *submatrix sub(A)
 * @param[in] descA array of dimension 9. The array descriptor of the distributed matrix A.
 * @param[in] matrixB Array, size lld_b by klb, where klb is LOCc(jb+n-1) when transb = 'N' or 'n',
 *and is LOCq(jb+k-1) otherwise. Before entry this array must contain the local pieces of the
 *distributed matrix sub(B).
 * @param[in] ib the row index in the distributed matrix B indicating the first row and the first
 *column of the submatrix sub(B)
 * @param[in] jb the column index in the distributed matrix B indicating the first row and the first
 *column of the submatrix sub(B)
 * @param[in] descB array of dimension 9. The array descriptor of the distributed matrix B.
 * @param[in] beta When beta is equal to zero, then sub(C) need not be set on input.
 * @param[in] ic the row index in the distributed matrix C indicating the first row of the submatrix
 *sub(C)
 * @param[in] jc the column index in the distributed matrix C indicating the first column of the
 *submatrix sub(C) param[in] descC array of dimension 9. The array descriptor of the distributed
 *matrix C.
 * @param[inout] matrixC array, size (lld_a, LOCq(jc+n-1)). Before entry this array must contain the
 *local pieces of the distributed matrix sub(C). Overwritten by the m-by-n distributed matrix
 *alpha*op(sub(A))*op(sub(B)) + beta*sub(C).
 *******************************************************************************************/
void pdgemm_(char*       transA,
             char*       transB,
             const Int*  N1,
             const Int*  N2,
             const Int*  N3,
             const Real* alpha,
             const Real* matrixA,
             const Int*  ia,
             const Int*  ja,
             const Int*  descA,
             const Real* matrixB,
             const Int*  ib,
             const Int*  jb,
             const Int*  descx,
             const Real* beta,
             Real*       matrixC,
             const Int*  ic,
             const Int*  jc,
             const Int*  descy);

void pdgemv_(char*       transA,
             const Int*  n1,
             const Int*  n2,
             const Real* alpha,
             const Real* A,
             const Int*  ia,
             const Int*  ja,
             const Int*  descA,
             const Real* x,
             const Int*  ix,
             const Int*  jx,
             const Int*  descX,
             const Int*  incx,
             const Real* beta,
             Real*       y,
             const Int*  iy,
             const Int*  jy,
             const Int*  descY,
             const Int*  incy);
}
#endif

} //namespace vaspml::scalapack_interface

#endif

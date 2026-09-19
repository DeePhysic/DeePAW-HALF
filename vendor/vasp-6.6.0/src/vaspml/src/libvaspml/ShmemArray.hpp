#ifndef SHMEMARRAY_HPP
#define SHMEMARRAY_HPP

#include "MlMPI.hpp"

#include "Tutor.hpp"
#include "debug.hpp"
#include "math.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <algorithm>  // for all_of
#include <fstream>    // for basic_ofstream
#include <numeric>    // for accumulate
#include <stdexcept>  // for out_of_range, runtime_error
#include <tuple>      // for tie, tuple

// debug
#include <iostream>

namespace vaspml
{

/*******************************************************************************************
 * @class ShmemArray
 * @brief Can be used to store arrays in a shared memory either with sysv or MPI.
 *
 * The class can be used to allocate an array in shared memory. To allocate the
 * array in shared memory the user has to supply a shared_ptr to an MlMPI class
 * from which a shared memory communicator can be created.
 * The class can for example be created with \n
 * @code
 * ShmemArray<Real> array( n, mpiIn ); \n
 * @endcode
 * And then the shared memory array can be filled with \n
 * @code
 * for ( Size i = 0; i < n; n++ ) \n
 * if ( array[ "mpiShmem" ].get_rank() == 0 ) array.set_value( i, value ); \n
 * array[ "mpiShmem" ].barrier(); \n
 * @endcode
 * @note when the functionality should be used with system V the code has to
 * be compiled with -Dsysv
 *******************************************************************************************/
template<class T>
class ShmemArray
{
  public:
    /*******************************************************************************************
     * default constructor
     *
     * initializes variables either to zero or nullptr
     *******************************************************************************************/
    ShmemArray();
    /*******************************************************************************************
     * constructor
     * @param n number of elements of data array
     * @param mpiShared supply MPI class with shared memory support if shared memory should be used
     *******************************************************************************************/
    ShmemArray(const Int n, const std::shared_ptr<MlMPI>& mpiIn = nullptr);
    /*******************************************************************************************
     * deconstruct the class.
     *
     * deallocates and frees the occupied memory
     *******************************************************************************************/
    ~ShmemArray();
    // if one of these are needed check folder
    // /fsc/home/jona/DATA/Documents/work/Redesign/Graveyard
    // there are implementations for the serial part
    /// deleted copy constructor
    ShmemArray(const ShmemArray& other) = delete;
    /// deleted copy assignement constructor
    ShmemArray& operator=(const ShmemArray& other) = delete;
    /// deleted move constructor
    ShmemArray(ShmemArray&& other) noexcept = delete;
    /// deleted move assignement operator
    ShmemArray& operator=(ShmemArray&& other) noexcept = delete;
    /*******************************************************************************************
     * write contents of the Shmem array to the screen.
     *
     * array will be written in three column format. First column is rank and
     * second is index in array and the last column is the value stored
     *******************************************************************************************/
    void print() const;
    /*******************************************************************************************
     * get the total size of the shmem array in number of elements
     *
     * @note since the class is parent class for ShmemArray2D and ShememArray2DVariableLen
     * the  value printed will then be the total number of entries.
     *******************************************************************************************/
    Size get_size() const;
    /*******************************************************************************************
     * set value in array
     * @param n index which sould be changed
     * @param value determines to value entry will be set
     * @param rank determines which rank does the writting
     *
     * @note after function call the user has to call the mpi barrier function.
     * This can be done with arrayName[ "mpiShmem" ].barrier(); arrayName denotes the
     * instance of the created ShmemArray
     *******************************************************************************************/
    void set_value(const Size n, T value, Int rank = 0);
    /*******************************************************************************************
     * return value stored in array
     *
     * @param indx denotes the index of the element which should be returned
     *
     * @note value will be copied
     *******************************************************************************************/
    T get_value(Size indx) const;
    /*******************************************************************************************
     * directly access mpi comunicator used in the class from outside the class.
     *
     * @note these are meant for convienience if a mpi barrier is used
     *******************************************************************************************/
    const MlMPI* operator->() const;
    /*******************************************************************************************
     * directly access mpi comunicator used in the class from outside the class.
     *
     * @note these are meant for convienience if a mpi barrier is used.
     * @note MlMPI is returned by reference without const. The user can not alter the
     * MlMPI class from outside. The const is neglected since
     * the mlMPI class has to set it's status info variable reporting if an error had
     * happened
     *******************************************************************************************/
    const MlMPI& operator[](const String& key) const;
    /*******************************************************************************************
     * distribute the shared memory segement internode. This function has to be called after
     * shared memory was filled
     *******************************************************************************************/
    void distributeInternode();
    /*******************************************************************************************
     * return a reference to the stored memory segment.
     *
     * @warning user can change the memory segment from the outside.
     * With great power comes great responsibility
     * @note this function is meant to be used in combination with cblas routines
     *******************************************************************************************/
    T*& get_dataPointer();
    /*******************************************************************************************
     * return a pointer to the stored memory segment.
     *
     * can be used as any ordinary c++ const raw pointer
     *
     * @note a copy of the pointer pointing to the data will be made
     *******************************************************************************************/
    const T* get_dataPointer() const;
    bool     get_useMPI() const;
    /*******************************************************************************************
     * Write contents of array into a file.
     *
     * @param fileName name of file to which the array is written.
     *******************************************************************************************/
    void writeToFile(const String& fileName) const;

  protected:
    /*******************************************************************************************
     * allocating the data pointer.
     *
     * The data pointer will be allocated as a standard c++ raw pointer when no MlMPI
     * class was supplied to the constructor. When MlMPI is supplied the data
     * will be either allocated with the system V approach or MPI to shared memory
     *******************************************************************************************/
    void allocate();
    /*******************************************************************************************
     * data pointer in which data is stored in either shared memory or c++ raw pointer format
     *******************************************************************************************/
    T* data;
    T* dataOffset;
    /*******************************************************************************************
     * total number of elements of type T stored in shared memory array
     *******************************************************************************************/
    Size size;
    /*******************************************************************************************
     * storing mpi communicator
     *
     * This map stores a shared memory communicator and a inter node communicator.
     * The shared mmemory communicator can be accessed with mpiShmem and the internode
     * communicator can be accessed with key mpiInter.
     *******************************************************************************************/
    std::map<String, MlMPI> MPIcomms;
    /*******************************************************************************************
     * determines if an MlMPI class was supplied so that data can be stored to shared memory
     *
     * @note variable will also be set to false when MlMPI was supplied but the code
     * was compiled in a way that no shared memory is used.
     *******************************************************************************************/
    // TODO useMPI is in principle useShmem
    bool useMPI;
};

template<class T>
ShmemArray<T>::ShmemArray()
{
    size = 0;
    data = nullptr;
    useMPI = false;
}

template<class T>
ShmemArray<T>::ShmemArray(const Int n, const std::shared_ptr<MlMPI>& mpiIn)
{
    if (n < 0)
    {
        throw std::runtime_error(
            "ERROR:ShmemArray<T>::ShmemArray( const Int n, const std::shared_ptr<MlMPI>& mpiIn )\n"
            "array size smaller zero. This does not make sense\n");
    }
    size = n;
    if (mpiIn == nullptr) { useMPI = false; }
    else
    {
        std::tie(MPIcomms["mpiShmem"], MPIcomms["mpiInter"]) = mpiIn->make_splitSharedInternode();
        useMPI = true;
    }
#ifndef use_shmem
    useMPI = false;
#endif
    allocate();
}

template<class T>
ShmemArray<T>::~ShmemArray()
{
    if (useMPI)
    {
        if (MPIcomms["mpiShmem"].get_sharedAllocCheck()) { MPIcomms["mpiShmem"].freeWindow(); }
    }
    else { delete[] data; }
}

template<class T>
Size ShmemArray<T>::get_size() const
{
    return size;
}

template<class T>
void ShmemArray<T>::set_value(const Size indx, T value, Int rank)
{
    if (!useMPI) { data[indx] = value; }
    else
    {
        if (rank == MPIcomms["mpiShmem"].get_rank()) { data[indx] = value; }
    }

    return;
}

template<class T>
T ShmemArray<T>::get_value(Size indx) const
{
    return data[indx];
}

template<class T>
void ShmemArray<T>::print() const
{
    for (Size i = 0; i < size; i++)
    {
        std::cout << " proc = " << MPIcomms.at("mpiShmem").get_rank() << "   " << i << "  "
                  << data[i] << std::endl;
    }

    return;
}

template<class T>
void ShmemArray<T>::allocate()
{
#ifdef use_shmem
    if (MPIcomms["mpiShmem"].get_sharedCommunicatorCheck()
        || MPIcomms["mpiShmem"].get_numberRanks() > 1)
    {
        // allocate array with shared memory
        MPIcomms["mpiShmem"].allocateShmemSpace(data, size);
    }
    else
    {
        if (MPIcomms["mpiShmem"].get_sharedCommunicatorCheck())
        {
            std::cout << "Note you are allocating a ShmemArray<T> without " << std::endl;
            std::cout << "supplying a shared memory communicator" << std::endl;
            std::cout << "Memory will be allocated on all cores" << std::endl;
        }
        // allocate array with non shared memory or single core
        data = new T[size];
    }
#else
    data = new T[size];
#endif

    return;
}

template<class T>
const MlMPI* ShmemArray<T>::operator->() const
{
    VASPML_DEBUG_L1(
        if (!useMPI)
        {
            global::tutor.bug("ERROR: " + flf(VASPML_FLF)
                              + ".Trying to access MPI communicator in ShmemArray but array was "
                                "set up without MPI.");
        }
    );

    return &MPIcomms.at("mpiShmem");
}

template<class T>
const MlMPI& ShmemArray<T>::operator[](const String& key) const
{
    VASPML_DEBUG_L1(
        if (!useMPI)
        {
            global::tutor.bug("ERROR: " + flf(VASPML_FLF)
                              + ".Trying to access MPI communicator in ShmemArray but array was "
                                "set up without MPI.");
        }
        auto it = MPIcomms.find(key);
        if (it == MPIcomms.end())
        {
            global::tutor.bug("ERROR: " + flf(VASPML_FLF) + ".\n"
                              + "key not found in MPIcomms dictionary. Possible values are \n"
                              + "mpiShmem or mpiInter");
        }
    );

    return MPIcomms.at(key);
}

template<class T>
void ShmemArray<T>::distributeInternode()
{
    if (useMPI)
    {
        // use inter node communicator to send from root to other
        if (MPIcomms["mpiShmem"].get_rank() == 0) MPIcomms["mpiInter"].bcast(data, size, 0);
        MPIcomms["mpiShmem"].barrier();
    }

    return;
}

template<class T>
T*& ShmemArray<T>::get_dataPointer()
{
    return data;
}

template<class T>
const T* ShmemArray<T>::get_dataPointer() const
{
    return data;
}

template<class T>
bool ShmemArray<T>::get_useMPI() const
{
    return useMPI;
}

template<class T>
void ShmemArray<T>::writeToFile(const String& fileName) const
{
    bool write = false;
    if (useMPI)
    {
        if (MPIcomms.at("mpiShmem").get_rank() == 0) write = true;
    }
    else write = true;

    if (write)
    {
        auto file = file_io::openFileO(fileName);
        for (Size i = 0; i < size; i++) { file << str("%24.16E ", data[i]) << std::endl; }
        file.close();
    }

    return;
}

/*******************************************************************************************
 * @class ShmemArray2D
 * @brief Can be used to store arrays in a shared memory either with sysv or MPI.
 *
 * The class can be used to allocate an array in shared memory. To allocate the
 * array in shared memory the user has to supply a shared_ptr to an MlMPI class
 * from which a shared memory communicator can be created.
 * The class can for example be created with \n
 * @code
 * ShmemArray<Real> array( n0, n1, mpiIn ); \n
 * @endcode
 * And then the shared memory array can be filled with \n
 * @code
 * for ( Size i = 0; i < n0; n++ )
 * for ( Size j = 0; j < n1; n++ )
 * if ( array[ "mpiShmem" ].get_rank() == 0 ) array.set_value( i,j, value );\n
 * array[ "mpiShmem" ].barrier();\n
 * @endcode
 * @note the faster index of the ShmemArray2D will always be the second index
 * @note when the functionality should be used with system V the code has to
 * be compiled with -Dsysv
 *******************************************************************************************/
template<class T>
class ShmemArray2D : public ShmemArray<T>
{
  public:
    /*******************************************************************************************
     * When using in default constructed way data will not be stored in shared memory
     *******************************************************************************************/
    ShmemArray2D();
    /*******************************************************************************************
     * construct the ShmemArray2D
     *
     * @param n0 first dimension of the shared memory array
     * @param n1 second dimension of the shared memory array
     * @mpiIn mpi communicator used create the shared memory communicators
     *******************************************************************************************/
    ShmemArray2D(const Int n0, const Int n1, const std::shared_ptr<MlMPI>& mpiIn = nullptr);
    /*******************************************************************************************
     * set value in array
     * @param n index which sould be changed
     * @param value determines to value entry will be set
     * @param rank determines which rank does the writting
     *
     * @note always call an mpi barrier with the proper communicator present in your class
     * after assigning a value
     *******************************************************************************************/
    void set_value(const Size indx0, const Size indx1, T value, Int rank = 0);
    /*******************************************************************************************
     * @brief Adds a value to a specific element in a 2D shared memory array.
     *
     * This function updates the value at the specified indices in the shared memory array.
     * If MPI is enabled, the operation is performed only on the rank specified.
     * If MPI is not enabled, the operation is performed directly on the local data.
     *
     * @tparam T The type of the elements in the shared memory array.
     * @param indx0 The row index of the element to be updated.
     * @param indx1 The column index of the element to be updated.
     * @param value The value to be added to the specified element.
     * @param rank (Optional) The MPI rank that performs the operation. Defaults to 0.
     *
     * @throws std::out_of_range If `indx0` is greater than or equal to `dimension0` or
     *                           if `indx1` is greater than or equal to `dimension1`.
     *
     * @details
     * - If `useMPI` is `false`, the value is directly added to the local data array.
     * - If `useMPI` is `true`, the value is added only if the current rank matches the specified
     *`rank`.
     * - The function performs bounds checking on `indx0` and `indx1` to ensure they are within
     *valid ranges.
     *
     * Example usage:
     * @code
     * ShmemArray2D<double> array;
     * array.add_value(2, 3, 5.0); // Adds 5.0 to the element at (2, 3)
     * @endcode
     *******************************************************************************************/
    void add_value(const Size indx0, const Size indx1, T value, Int rank = 0);
    /*******************************************************************************************
     * return value stored in array
     *
     * @param indx index of element which should be returned from array
     * @note a copy of the value will be made
     *******************************************************************************************/
    T        get_value(const Size indx0, const Size indx1) const;
    const T* get_slice(const Size indx0) const;
    /*******************************************************************************************
     * getting the first dimension of the shared memory array
     *******************************************************************************************/
    Size get_size0() const;
    /*******************************************************************************************
     * getting the second dimension of the shared memory array
     *******************************************************************************************/
    Size get_size1() const;

  private:
    /*******************************************************************************************
     * variable to store the first dimension of the shared memory array
     *
     * total size of the shared memory array will be dimension0 * dimension1
     *******************************************************************************************/
    Size dimension0;
    /*******************************************************************************************
     * variable to store the second dimension of the shared memory array
     *
     * total size of the shared memory array will be dimension0 * dimension1
     *******************************************************************************************/
    Size dimension1;
};

template<class T>
ShmemArray2D<T>::ShmemArray2D() : ShmemArray<T>()
{
    dimension0 = 0;
    dimension1 = 0;
}

template<class T>
ShmemArray2D<T>::ShmemArray2D(const Int n0, const Int n1, const std::shared_ptr<MlMPI>& mpiIn) :
    ShmemArray<T>(n0 * n1, mpiIn)
{
    if (n0 < 0)
    {
        throw std::runtime_error(
            "ERROR:ShmemArray2D<T>::ShmemArray2D( const Int n0, const Int n1 \n"
            "                                     const std::shared_ptr<MlMPI>& mpiIn )\n"
            "array size smaller zero. This does not make sense\n");
    }
    if (n1 < 0)
    {
        throw std::runtime_error(
            "ERROR:ShmemArray2D<T>::ShmemArray2D( const Int n0, const Int n1 \n"
            "                                     const std::shared_ptr<MlMPI>& mpiIn )\n"
            "array size smaller zero. This does not make sense\n");
    }
    dimension0 = n0;
    dimension1 = n1;
}

template<class T>
void ShmemArray2D<T>::set_value(const Size indx0, const Size indx1, T value, Int rank)
{
    VASPML_DEBUG_L1(
        if (indx0 >= dimension0)
        {
            throw std::out_of_range("ERROR: ShmemArray2D<T>::set_value( const Size indx0, "
                                    "const Size indx1, T value, Int rank )\n"
                                    "indx0 is out of bounds");
        }
        if (indx1 >= dimension1)
        {
            throw std::out_of_range("ERROR: ShmemArray2D<T>::set_value( const Size indx0, "
                                    "const Size indx1, T value, Int rank )\n"
                                    "indx1 is out of bounds");
        }
    );
    if (!this->useMPI)
    {
        Size indx = indx0 * dimension1 + indx1;
        this->data[indx] = value;
    }
    else
    {
        if (rank == this->MPIcomms["mpiShmem"].get_rank())
        {
            Size indx = indx0 * dimension1 + indx1;
            this->data[indx] = value;
        }
        //this -> MPIcomms[ "mpiShmem" ].barrier();
    }
}

template<class T>
void ShmemArray2D<T>::add_value(const Size indx0, const Size indx1, T value, Int rank)
{
    VASPML_DEBUG_L1(
        if (indx0 >= dimension0)
        {
            throw std::out_of_range("ERROR: ShmemArray2D<T>::add_value( const Size indx0, "
                                    "const Size indx1, T value, Int rank )\n"
                                    "indx0 is out of bounds");
        }
        if (indx1 >= dimension1)
        {
            throw std::out_of_range("ERROR: ShmemArray2D<T>::add_value( const Size indx0, "
                                    "const Size indx1, T value, Int rank )\n"
                                    "indx1 is out of bounds");
        }
    );
    if (!this->useMPI)
    {
        Size indx = indx0 * dimension1 + indx1;
        this->data[indx] += value;
    }
    else
    {
        if (rank == this->MPIcomms["mpiShmem"].get_rank())
        {
            Size indx = indx0 * dimension1 + indx1;
            this->data[indx] += value;
        }
        //this -> MPIcomms[ "mpiShmem" ].barrier();
    }

    return;
}

template<class T>
T ShmemArray2D<T>::get_value(Size indx0, Size indx1) const
{
    VASPML_DEBUG_L1(
        if (indx0 >= dimension0)
        {
            throw std::out_of_range("ERROR: ShmemArray2D<T>::get_value( const Size indx0, "
                                    "const Size indx1 )\n"
                                    "indx0 is out of bounds");
        }
        if (indx1 >= dimension1)
        {
            throw std::out_of_range("ERROR: ShmemArray2D<T>::get_value( const Size indx0, "
                                    "const Size indx1 )\n"
                                    "indx1 is out of bounds");
        }
    );
    Size indx = indx0 * dimension1 + indx1;

    return this->data[indx];
}

template<class T>
const T* ShmemArray2D<T>::get_slice(Size indx0) const
{
    VASPML_DEBUG_L1(
        if (indx0 >= dimension0)
        {
            throw std::out_of_range("ERROR: ShmemArray2D<T>::get_value( const Size indx0, "
                                    "const Size indx0 )\n"
                                    "indx0 is out of bounds");
        }
    );
    Size indx = indx0 * dimension1;

    return &this->data[indx];
}

template<class T>
Size ShmemArray2D<T>::get_size0() const
{
    return dimension0;
}

template<class T>
Size ShmemArray2D<T>::get_size1() const
{
    return dimension1;
}

/*******************************************************************************************
 * @class ShmemArray2DVariableLen
 * @brief Can be used to store arrays in a shared memory either with sysv or MPI.
 *
 * The class can be used to allocate an array in shared memory. To allocate the
 * array in shared memory the user has to supply a shared_ptr to an MlMPI class
 * from which a shared memory communicator can be created.
 * The class can for example be created with
 * @code
 * ShmemArray<Real> array( n0, n1, mpiIn ); \n
 * @endcode
 * And then the shared memory array can be filled with \n
 * @code
 * for ( Size i = 0; i < n0; n++ )
 * for ( Size j = 0; j < n1[ i ]; n++ )
 * if ( array[ "mpiShmem" ].get_rank() == 0 ) array.set_value( i,j, value ); \n
 * array[ "mpiShmem" ].barrier(); \n
 * @endcode
 * @note the faster index of the ShmemArray2DVariableLen will always be the second index
 * @note when the functionality should be used with system V the code has to
 * be compiled with -Dsysv
 *******************************************************************************************/
template<class T>
class ShmemArray2DVariableLen : public ShmemArray<T>
{
  public:
    /*******************************************************************************************
     * When using in default constructed way data will not be stored in shared memory
     *******************************************************************************************/
    ShmemArray2DVariableLen();
    /*******************************************************************************************
     * constructing a ShmemArray2DVariableLen
     *
     * @param lengths the lengths of the second dimensions of the ShmemArray. The size of lengths
     *determines the first dimension
     * @param mpiIn mpi class to make shared memory access and internode communication possible
     *
     *******************************************************************************************/
    ShmemArray2DVariableLen(const Vec1Int& lengths, const std::shared_ptr<MlMPI>& mpiIn = nullptr);
    /*******************************************************************************************
     * set single value in array
     *
     * @param indx0 index which should be changed ( slow index )
     * @param indx1 index which should be changed ( fast index )
     * @param value determines to value entry will be set
     * @param shmemRank determines which mpiShmem rank does the writing
     * @param interRank determines which mpiInter rank does the writing
     *
     * @note only the cores for which both shmemRank and interRank are satisfied
     * will do the writing
     *******************************************************************************************/
    void set_value(const Size indx0, const Size indx1, const T value, const Int rank = 0);
    /*******************************************************************************************
     * copy data array to the shared memory array
     *
     * @param inputData data array which will be copied to the shared memory array
     * @param shmemRank determines which mpiShmem rank does the writing
     * @param interRank determines which mpiInter rank does the writing
     *
     * @note only the cores for which both shmemRank and interRank are satisfied
     * will do the writing
     *******************************************************************************************/
    void set_value(const std::vector<std::vector<T>>& inputData,
                   const Int                          shemRank = 0,
                   const Int                          interRank = 0);
    /*******************************************************************************************
     * return value stored in array
     *
     * @param indx0 entry which should be returned from array ( slow index )
     * @param indx1 entry which should be returned from array ( fast index )
     *******************************************************************************************/
    T get_value(const Size indx0, const Size indx1) const;
    /*******************************************************************************************
     * get a const pointer to a slice in the second dimension of ShmemArray
     *
     * @param indx0 index of the first dimension for which the slice should be obtained
     *
     * returns a pointer to a slice chosen by indx0
     * as data[ indx0,: ],
     * @note that a pointer has to be created
     *******************************************************************************************/
    const T* get_slice(const Size indx0) const;
    /**
     * get a pointer to a slice in the second dimension of ShmemArray
     *
     * @param indx0 index of the first dimension for which the slice should be obtained
     *
     * returns a pointer to a slice chosen by indx0
     * as data[ indx0,: ], note that a pointer has to be created
     *******************************************************************************************/
    T*  get_slice(const Size indx0);
    T*& get_sliceReference(const Size indx0);
    /*******************************************************************************************
     * returns the number of elements stored in a certain slice
     *
     * @param indx0 index for the first entry of the array
     *******************************************************************************************/
    Size get_size1(const Size indx0) const;
    /*******************************************************************************************
     * get size of zeroth dimension of the shared memory array
     *******************************************************************************************/
    Size get_size0() const;

  private:
    /*******************************************************************************************
     * length of first dimension of ShmemArray
     *******************************************************************************************/
    Size dimension0;
    /*******************************************************************************************
     * lengths of the second dimensions.
     *
     * lengths.size() matches dimension0
     *******************************************************************************************/
    Vec1Size lengths;
    /*******************************************************************************************
     * offsets which tell at which entries of the 1D contiguous shared memory array the slow index
     * is increased
     *******************************************************************************************/
    Vec1Size offsets;
};

template<class T>
Size ShmemArray2DVariableLen<T>::get_size0() const
{
    return dimension0;
}

template<class T>
ShmemArray2DVariableLen<T>::ShmemArray2DVariableLen() : ShmemArray<T>()
{
    dimension0 = 0;
    lengths.resize(0);
    offsets.resize(0);
}

template<class T>
ShmemArray2DVariableLen<T>::ShmemArray2DVariableLen(const Vec1Int&                lengths,
                                                    const std::shared_ptr<MlMPI>& mpiIn) :
    ShmemArray<T>(math::sumVector(lengths), mpiIn)
{
    VASPML_DEBUG_L1(
        if (!std::all_of(lengths.begin(), lengths.end(), [](Int i) { return i >= 0; }))
        {
            throw std::runtime_error("ERROR:ShmemArray2DVariableLen<T>::ShmemArray2DVariableLen: "
                                     "Some of your vector dimensions are below zero");
        }
    );

    dimension0 = lengths.size();
    offsets.resize(dimension0);
    this->lengths.resize(dimension0);
    this->lengths[0] = lengths[0];
    for (Size i = 1; i < lengths.size(); i++)
    {
        offsets[i] = offsets[i - 1] + lengths[i - 1];
        this->lengths[i] = lengths[i];
    }
}

template<class T>
void ShmemArray2DVariableLen<T>::set_value(const Size indx0,
                                           const Size indx1,
                                           const T    value,
                                           const Int  rank)
{
    VASPML_DEBUG_L1(
        if (indx0 >= dimension0)
        {
            throw std::out_of_range(
                "ERROR: ShmemArray2DVariableLen<T>::set_value( const Size indx0, const "
                "Size indx1, const T value, const Int rank )\n"
                "indx0 is out of bounds");
        }
        if (indx1 >= lengths[indx0])
        {
            throw std::out_of_range(
                "ERROR: ShmemArray2DVariableLen<T>::set_value( const Size indx0, const "
                "Size indx1, const T value, const Int rank )\n"
                "indx1 is out of bounds");
        }
    );
    Size indx = offsets[indx0] + indx1;
    if (!this->useMPI) { this->data[indx] = value; }
    else
    {
        if (rank == this->MPIcomms["mpiShmem"].get_rank()) { this->data[indx] = value; }
        //this->MPIcomms[ "mpiShmem" ].barrier();
    }

    return;
}

template<class T>
void ShmemArray2DVariableLen<T>::set_value(const std::vector<std::vector<T>>& inputData,
                                           [[maybe_unused]] const Int         shmemRank,
                                           [[maybe_unused]] const Int         interRank)
{
    VASPML_DEBUG_L1(
        if (this->MPIcomms["mpiShmem"].get_rank() == shmemRank
            && this->MPIcomms["mpiInter"].get_rank() == interRank)
        {

            Size inputSize = std::accumulate(inputData.begin(),
                                             inputData.end(),
                                             0,
                                             [](Size sum, const std::vector<T>& inner_vec)
                                             { return sum + inner_vec.size(); });

            if (inputSize != this->size)
            {
                throw std::out_of_range("ERROR: ShmemArray2DVariableLen<T>::set_value( const "
                                        "std::vector<std::vector<T>>& inputData, const Int "
                                        "shmemRank, const Int interRank )\n"
                                        "Size of inputData and number of elements in "
                                        "ShmemArray2DVariableLen does not match\n");
            }
        }
    );

    if (!this->useMPI)
    {
        Size counter = 0;
        for (Size indx0 = 0; indx0 < dimension0; indx0++)
        {
            for (Size indx1 = 0; indx1 < lengths[indx0]; indx1++)
            {
                this->data[counter] = inputData[indx0][indx1];
                counter++;
            }
        }
    }
    else
    {
#ifdef use_shmem
        if (this->MPIcomms["mpiShmem"].get_rank() == shmemRank
            and this->MPIcomms["mpiInter"].get_rank() == interRank)
        {
#endif
            Size counter = 0;
            for (Size indx0 = 0; indx0 < dimension0; indx0++)
            {
                for (Size indx1 = 0; indx1 < lengths[indx0]; indx1++)
                {
                    this->data[counter] = inputData[indx0][indx1];
                    counter++;
                }
            }
#ifdef use_shmem
        }
        this->MPIcomms["mpiShmem"].barrier();
        this->distributeInternode();
#endif
    }

    return;
}

template<class T>
T ShmemArray2DVariableLen<T>::get_value(Size indx0, Size indx1) const
{
    VASPML_DEBUG_L1(
        if (indx0 >= dimension0)
        {
            throw std::out_of_range(
                "ERROR: ShmemArray2DVariableLen<T>::get_value( const Size indx0, const "
                "Size indx1, T value, Int rank )\n"
                "indx0 is out of bounds");
        }
        if (indx1 > lengths[indx0])
        {
            throw std::out_of_range(
                "ERROR: ShmemArray2DVariableLen<T>::get_value( const Size indx0, const "
                "Size indx1, T value, Int rank )\n"
                "indx1 is out of bounds");
        }
    );
    Size indx = offsets[indx0] + indx1;

    return this->data[indx];
}

template<class T>
const T* ShmemArray2DVariableLen<T>::get_slice(const Size indx0) const
{
    VASPML_DEBUG_L1(
        if (indx0 >= dimension0)
        {
            throw std::out_of_range("ERROR: const T*& ShmemArray2DVariableLen<T>::get_slice( const "
                                    "Size indx0 ) const\n"
                                    "indx0 is out of bounds");
        }
    );

    return this->data + offsets[indx0];
}

template<class T>
T* ShmemArray2DVariableLen<T>::get_slice(const Size indx0)
{
    VASPML_DEBUG_L1(
        if (indx0 >= dimension0)
        {
            throw std::out_of_range(
                "ERROR: T* ShmemArray2DVariableLen<T>::get_slice( const Size indx0 )\n"
                "indx0 is out of bounds");
        }
    );

    return this->data + offsets[indx0];
}

template<class T>
T*& ShmemArray2DVariableLen<T>::get_sliceReference(const Size indx0)
{
    VASPML_DEBUG_L1(
        if (indx0 >= dimension0)
        {
            throw std::out_of_range(
                "ERROR: T* ShmemArray2DVariableLen<T>::get_slice( const Size indx0 )\n"
                "indx0 is out of bounds");
        }
    );
    this->dataOffset = this->data + offsets[indx0];

    return this->dataOffset;
}

template<class T>
Size ShmemArray2DVariableLen<T>::get_size1(const Size indx0) const
{
    VASPML_DEBUG_L1(
        if (indx0 >= dimension0)
        {
            throw std::out_of_range("ERROR: Size ShmemArray2DVariableLen<T>::get_size1( const Size "
                                    "indx0 ) const \n"
                                    "indx0 is out of bounds");
        }
    );

    return lengths[indx0];
}

} //namespace vaspml

#endif

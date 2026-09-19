!> pm:  implementation of the Lanczos algorithm for computing the absorption spectrum
!> Implementation is based on the reformulation of the problem of computing the dielectric function, \epsilon_{ij}(\omega) as a continued fraction
!> For the diagonal elements, {ii},  we have
!>                                                    |u_0^{ii}|^2  
!> \epsilon_{ii}(\omega) = 1 - ----------------------------------------------------------        (1)
!>                                                                    b_1
!>                             \omega  - a_1 + i\eta + ----------------------------------
!>                                                                                b_2
!>                                                     \omega - a_2 + i\eta + -----------
!>                                                                                ...
!
!> where:
!> i)   |u_0^{ii}|^2 is the amplitude of the initial guess vector, before being normalized to 1. We take the dipole moments as the initial guess, 
!> from OPTMAT, so {ii} is the direction of the dipole.
!> ii)  a_i and b_i are the Lanczos coefficients, and are real numbers.
!> iii) \omega is the frequency, and \eta is the width, controlled by CSHIFT inside the INCAR.
!> 
!> For the off-diagonal elements, this has to be modified. The assumption at the start is that \epsilon has an operational form like
!>                                 1
!> \epsilon(\omega) = 1 - <P|------------|P> (2),
!>                            \omega - H
!
!> where H is the BSE Hamiltonian, and <P|P> = |u_0|^2 is real. But if the left and right vectors are not the same, i.e. they are dipole moments 
!> along different directions, we cannot follow this and there is no guarantee that <P|P> is real. We solve this by rotating the matrix 45 degrees 
!> and using
!>                                         1
!> \bar\epsilon_{m}(\omega) = - <R_m|------------|R_m> (3),
!>                                    H - \omega
!> where |R_m> is the rotated dipole vector. We will have
!> |R_1> = 1/\sqrt{2}(|P_x> + |P_y>)   (4a)
!> |R_2> = 1/\sqrt{2}(|P_y> + |P_z>)   (4b)
!> |R_3> = 1/\sqrt{2}(|P_z> + |P_x>)   (4c)
!
!> This way the vectors on the left and right side of (3) are the same, and <R_m|R_m> is a real number. At the end the three \bar\epsilon auxiliary 
!> functions are rotated back so that we can recover the {xy}, {xz}, and {yz} components
!> 
! 
#include "symbol.inc"
#include "bse.inc"
module bse_lanczos

  use prec     
  use constant 
  use wave_high     
  use bse_base
  use bse_struct
  use vaspxml
  use chi_glb
  use iso_c_binding
#ifdef VASP_HDF5
  use vhdf5
#endif

  implicit none

! Some further (local) variables, kept here until the BSE refactoring branch is merged 
  integer,private  :: matdim !< Matrix dimension used for the calculation (-> rank of bse matrix from bse::calculate_bse)
  integer,private, allocatable  :: matdim_node(:) !< Number of rows at each node (bse matrix distributed for MPI jobs ... !)
  integer,private, allocatable  :: matstart(:)    !< Corresponding "global index" of first local matrix element at each node

  complex(mq),private, parameter  :: imun     = (0.0_mq,1.0_mq) !< Imaginary unit
  complex(mq),private, parameter  :: czero    = (0.0_mq,0.0_mq) !< complex zero

  integer,private :: node_lead=0    !< number of the leading process.
  integer,private :: node_me        !< number of my process (set by vasp).
  integer,private :: nranks         !< total number of processes (set by vasp).

#ifdef scaLAPACK

contains


!****************** subroutine calculate_bse_lanczos    ******************************
!
!> main driver to setup the Haydock-lanczos solver
!> it should check whether we are using the TD approximation or not and send the 
!> necessary parameters to the subroutine that performs the iterative loop 
!
!*************************************************************************************

  subroutine calculate_bse_lanczos(comm,shift,omega,ispin,maxk,ncv,kptweight_in,encut_sp_in,lbeyondtd,nedos_in,io,nbseeig,qpoint,bsematrix,optmat)

    use base
    use constant
    use reader_tags
    use tutor, only: vtutor
    use chi_glb, only: lgwlf, bseprec, ismearbse, lanczoseigen
#ifdef VASP_HDF5
    use vhdf5
#endif
    implicit none 

    type(communic) :: comm              !< communicator to be used (set by vasp).
    real(q) :: shift                    !< complex shift (the i*eta) used in epsilon
    real(q) :: omega                    !< cell volume
    integer :: ispin,maxk,ncv           !< number of spin channels, maximum number of k-points, number of electron-hole pairs in the BSE
    real(q) :: kptweight_in(maxk)       !< weight of each k-point
    real(q) :: encut_sp_in              !< frequency cutoff 
    logical :: lbeyondtd                !< beyond TDA or not
    integer :: nedos_in                 !< number of frequencies
    type (in_struct) :: io              !< input-output descriptor
    integer :: nbseeig                  !< number of BSE eigenvalues to be computed
    real(q) :: qpoint(3)                !< q-point for finite-momentum BSE
#ifdef single_prec_bse
    GDEFS, intent(in) :: bsematrix(:)
    GDEFS, intent(in) :: optmat(:,:)
#endif
#ifdef double_prec_bse
    GDEF, intent(in) :: bsematrix(:)
    GDEF, intent(in) :: optmat(:,:)
#endif
     
    integer :: madtdin                  !< rank of bsematrix, 
    real(q) :: prec_thr                 !< precision in computing epsilon
    logical :: lopen                    

    integer :: it, ierr, id1, n_iter
    real(q), allocatable :: r_and_intensity(:,:)
    character(40):: sp_type
    
    !> Lanczos coefficients and results
    real(q), allocatable :: a(:,:,:), b(:,:,:), tmp_a(:), tmp_b(:), rwork(:)
    complex(q), allocatable :: spectra1(:,:), spectra2(:,:), eps(:,:,:), tmp_m(:,:)
    real(mq) :: norm(1:6), pref
    real(q) :: error(1:6)
#ifdef USE_ZHEEV
    complex(q), allocatable ::  work(:)
#else
    real(q), allocatable :: work(:)
#endif
    character(len=10), dimension(6) :: edir
    !> HDF5 output related variables
#ifdef VASP_HDF5
    character(len=:), allocatable :: prefix
#endif 

    ! precision control section
    ! set default threshold for epsilon
    prec_thr = 1.d-3  !> this is the equivalent of BSEPREC=(M)EDIUM/(m)edium

    if ((bseprec(1:1)=='L') .or. (bseprec(1:1)=='l') .or. (bseprec(1:1)=='F') .or. (bseprec(1:1)=='f')) then
      prec_thr = 1.d-2
    elseif ((bseprec(1:1)=='H') .or. (bseprec(1:1)=='h')) then
      prec_thr = 1.d-4
    elseif ((bseprec(1:1)=='A') .or. (bseprec(1:1)=='a')) then
      prec_thr = 1.d-5
    elseif ((bseprec(1:1)=='M') .or. (bseprec(1:1)=='m') .or. (bseprec(1:1)=='N') .or. (bseprec(1:1)=='n')) then
      prec_thr = 1.d-3
    else 
      call vtutor%alert("Undefined value for tag BSEPREC in INCAR. Default value ('Normal') will be used!")
    endif

    sp_type = 'l'        ! use Lorentzian functions to approximate the Dirac delta-function
    if (ismearbse == 0) then
      sp_type = 'g'      ! use Gaussian functions to approximate the Dirac delta-function
    endif 

    ! go beyond tamm-dancoff approximation (antires==2) or not (antires/=2)?
    if (lbeyondtd) then
      if (io%iu0>=0) write(io%iu0,*) 'ANTIRES was set to 2 at some point during execution.'
      lbeyondtd = .false.
    else
      if (io%iu0>=0) write(io%iu0,*) 'Calling Haydock-Lanczos solver to compute optical spectra.'
    endif

    matdim=ncv
    if (lbeyondtd) matdim = matdim*2

    ! because we use VASP native units, the prefactor is the same as in the diagonalisation solver, i.e. 4\pi e^2/volume
    pref = edeps/omega*real(kptweight_in(1),kind=mq)
    if (ispin==1) pref=2._mq*pref

    ! a and b coefficients
    ! the a and b coefficients have to be stored in order to reproduce the full spectra using the continued fraction
    ! but the maximum rank of their array is the rank of the bse matrix, as the number of iterations cannot exceed the
    ! matrix's rank
    allocate(a(1:matdim,1:6,1:6),b(1:matdim,1:6,1:6))

    PROFILING_START('bse_lanczos')
    if (io%iu0>=0) write(io%iu0,*) 'Starting the solver.' 
    if (.not. lbeyondtd) call lanczos_solver_tda(comm,matdim,encut_sp_in,nedos_in,shift,io,nbseeig,prec_thr,a,b,norm,n_iter,error,bsematrix,optmat)

    if (io%iu0>=0) write(io%iu0,*) 'Solver iterative cycle ended.'

    edir = [character(len=10) :: 'x', 'y', 'z', 'xy', 'xz', 'yz']
    ! write information on convergence to the OUTCAR file

    if (io%iu6>=0) then
      write(io%iu6,*) ""
      write(io%iu6,*) "======================================================================="
      write(io%iu6,*) ""
      write(io%iu6,*) "Parameters set for Lanczos algorithm"
      write(io%iu6,'(a,e14.7)') "   Threshold: ", prec_thr 
      write(io%iu6,'(a,f8.3)')  "   Maximum energy: ", encut_sp_in
      write(io%iu6,*) "   Number of energy steps: ", nedos_in
      write(io%iu6,*) "   Number of iterations needed to reach set accuracy:", n_iter
      write(io%iu6,*) "   Final accuracy along all directions:"
      do id1 = 1, 6
        write(io%iu6,'(a,e14.7)') edir(id1), error(id1)
      enddo
      write(io%iu6,*) ""
      write(io%iu6,*) "======================================================================="
    endif

    PROFILING_START('final_lanczos_spectrum')

    allocate(spectra1(1:nedos_in,1:6),spectra2(1:nedos_in,1:6))
    spectra1 = czero; spectra2 = czero

    PROFILING_START('lanczos_eigen_decomposition')
    ! hidden tag lanczoseigen, which is always set to true, to force the routine to compute eigenvalues and amplitudes
    ! if set to false it will skip this step and compute only the spectrum
    if ( .not. lanczoseigen ) then
      call lanczos_spectra(a, b, n_iter, shift, 0, encut_sp_in, nedos_in, spectra1, norm, 0)
    elseif ( lanczoseigen ) then
      if (io%iu0 >=0 ) then
        write(io%iu0,*) "Computing ", n_iter, "eigenvalues of the BSE"
      endif

      allocate(tmp_a(1:n_iter),tmp_b(1:n_iter-1),tmp_m(1:n_iter,1:n_iter))
      tmp_a = 0._q; tmp_b=0._q; tmp_m=0._q

#ifdef USE_ZHEEV
      allocate(rwork(1:3*n_iter-2),work(1:2*n_iter-1))
      rwork = 0._q; work = czero
#else
      allocate(work(1:2*n_iter-2))
      work = 0._q
#endif

      allocate(r_and_intensity(1:2,1:n_iter))
      r_and_intensity = 0._q

      do id1 = 1, 6
#ifndef USE_ZHEEV
        tmp_a(1:n_iter) = a(1:n_iter,id1,id1)
        forall (it=1:n_iter-1) tmp_b(it) = sqrt(b(it,id1,id1))
        call zsteqr('i',n_iter, tmp_a(1), tmp_b(1), tmp_m(1,1), n_iter, work, ierr)
        if (ierr > 0) then
          call vtutor%warning("STEQR failed at converging all eigenvalues. Plotting the dielectric function using eigenvalues is not &
          guaranteed to give accurate results. Try compiling VASP with -DUSE_ZHEEV to use another eigenvalue solver.")
        endif
#else
        forall (it=1:n_iter) tmp_m(it,it) = cmplx(a(it,id1,id1),0._q)
        forall (it=1:n_iter-1) tmp_m(it+1,it) = cmplx(sqrt(b(it,id1,id1)),0._q)
        call zheev('v','l', n_iter, tmp_m, n_iter, tmp_a, work, 2*n_iter-1, rwork, ierr)
        if (ierr > 0) then
          call vtutor%warning("HEEV failed at converging all eigenvalues. Plotting the dielectric function using eigenvalues is not &
          guaranteed to give accurate results.")
        endif
#endif

        if (.not. lbeyondtd ) call lanczos_spectra_eigen(tmp_a, n_iter, tmp_m(1,:), nedos_in, encut_sp_in, norm(id1), shift, spectra1(:,id1), sp_type, 'f')
        if ( lbeyondtd )      call lanczos_spectra_eigen(tmp_a, n_iter, tmp_m(1,:), nedos_in, encut_sp_in, norm(id1), shift, spectra1(:,id1), sp_type, 't')
        ! store the eigenvalues from xx direction into r_and_intensity, which will be writen to vaspout.h5
        if (id1 == 1) then
          r_and_intensity(1,1:n_iter) = tmp_a(1:n_iter)
        endif
        !> sum the diagonal terms to take the trace
        if (id1 <= 3) then
          r_and_intensity(2,1:n_iter) = r_and_intensity(2,1:n_iter) + real(tmp_m(1,1:n_iter)*conjg(tmp_m(1,1:n_iter))*norm(id1))**2
        endif
      enddo
      
      call xml_vecarray("opticaltransitions")
      call xml_array_real(r_and_intensity,"(F10.4,' ')")
      call xml_close_tag("varray")

      deallocate(tmp_a,tmp_b,tmp_m)
      deallocate(a,b)
#ifdef USE_ZHEEV
      deallocate(rwork)
#endif
      deallocate(work)
    endif
    
    PROFILING_STOP('lanczos_eigen_decomposition')
    ! multiply by the prefactor
    spectra1 = pref*spectra1

    ! diagonal components are dealt with here
    allocate(eps(1:nedos_in,1:3,1:3))
    eps= czero
    
    do id1 = 1,3
      if (any(abs(qpoint)>1E-6_q)) then
        spectra1(1:nedos_in,id1) = spectra1(1:nedos_in,id1)/(tpi**2*sum(qpoint**2))
        eps(1:nedos_in,id1,id1) = 1/(1 - spectra1(1:nedos_in,id1))
      else
        eps(1:nedos_in,id1,id1) = 1._q + spectra1(1:nedos_in,id1)
      endif
    enddo ! id1

    ! index for spectra1(:,k) is 
    ! id1 = 1, id2 = 2 => k = 4
    ! id1 = 1, id2 = 3 => k = 6
    ! id1 = 2, id2 = 3 => k = 5
    ! I rotate the auxilar \bar\epsilon here back to the xy, xz, and yz axis according to
    ! \epsilon_{ij} = \bar\epsilon_{ij} - 0.5[\epsilon_{ii} + \epsilon_{jj} - 2], i\neq j
    if (.not. any(abs(qpoint)>1E-6_q)) then
      eps(:,1,2) = spectra1(:,4) - 0.5_q*(eps(:,1,1) + eps(:,2,2) - 2._q)
      eps(:,1,3) = spectra1(:,6) - 0.5_q*(eps(:,1,1) + eps(:,3,3) - 2._q)
      eps(:,2,3) = spectra1(:,5) - 0.5_q*(eps(:,2,2) + eps(:,3,3) - 2._q)
      eps(:,2,1) = eps(:,1,2); eps(:,3,2) = eps(:,2,3); eps(:,3,1) = eps(:,1,3)
    endif

    ! write the dielectric function in the xml file
    if (any(abs(qpoint)>1E-6_q)) then
      call xml_epsilon_scalar_w(encut_sp_in/(nedos_in-1),real(eps(:,1,1),q),aimag(eps(:,1,1)),nedos_in)
    else
      call xml_epsilon_w(encut_sp_in/(nedos_in-1),real(eps),aimag(eps),nedos_in)
    endif

#ifdef VASP_HDF5
    if (lgwlf) then
      prefix = "bse"
    else
      prefix = "tdhf"
    end if
    call vh5_write_optical_transitions(ih5outfileid, r_and_intensity)
    call vh5_write_dielectric_energies(ih5outfileid, encut_sp_in/(nedos_in-1), nedos_in, prefix)
    if (any(abs(qpoint)>1E-6_q)) then
      call vh5_write_dielectric_scalar_dynamic(ih5outfileid, eps(1:nedos_in,1,1), prefix)
    else
      call vh5_write_dielectric_dynamic(ih5outfileid, eps, prefix)
    endif
#endif

    deallocate(r_and_intensity)
    deallocate(eps,spectra1)

    PROFILING_STOP('final_lanczos_spectrum')
    PROFILING_STOP('bse_lanczos')

  end subroutine calculate_bse_lanczos

!****************** subroutine lanczos_bse_solver    ******************************
!
!> computes the haydock-lanczos coefficients and dielectric function 
!
!******************************************************************************

  subroutine lanczos_solver_tda(comm,fulldim,omegamax,nedos_in,shift,io,nbseeig,prec_eps,a,b,norm_v0,itercount,err,bsematrix,optmat)
    use base
    use string, only: str
    use tutor, only: vtutor
    USE c2f_interface, ONLY : vtime
#ifdef VASP_HDF5
    use vhdf5
#endif
    implicit none

    type(communic) :: comm                     !< communicator to be used (set by vasp).
    integer :: fulldim                         !< rank of bsematrix (i.e. ncv if Tamm-Dancoff, 2*ncv if beyond Tamm-Dancoff)
    real(q) :: omegamax                        !< maximum frequency of epsilon
    integer :: nedos_in                        !< number of frequencies in epsilon 
    real(q) :: shift                           !< complex shift used in epsilon
    type (in_struct) :: io                     !< input-output descriptor 
    integer :: nbseeig                         !< number of BSE eigenvalues
    real(q) :: prec_eps                        !< convergence threshold for epsilon
    real(q) :: a(:,:,:)                        !< a-Lanczos coefficients
    real(q) :: b(:,:,:)                        !< b-Lanczos coefficients
    real(mq) :: norm_v0(1:6)                   !< norm of the Lanczos vectors
    integer :: itercount                       !< iteration counter
    real(q) :: err(1:6)                       !< root-mean-square of epsilon
#ifdef single_prec_bse
    GDEFS, intent(in) :: bsematrix(:)
    GDEFS, intent(in) :: optmat(:,:)
#endif

#ifdef double_prec_bse
    GDEF, intent(in) :: bsematrix(:)
    GDEF, intent(in) :: optmat(:,:)
#endif
    ! like in the time-evolution code, every step computes a new vector from two previous iterations. These are stored in vector_0, vector_1, 
    ! and v_tmp_1
#ifdef gammareal
    real(mq), allocatable :: vector_0(:,:)               ! real part of initial and (n-1)th iteration vectors
    real(mq), allocatable :: vector_1(:,:)               ! real part of second and n-th iteration vectors 
    real(mq), allocatable :: v_tmp_1(:,:), v_tmp_2(:,:)  ! real part of the auxiliary vectors
    real(mq) :: ortho_test(1:6,1:6)                      ! used for all dot products between vectors
#else
    complex(mq), allocatable :: vector_0(:,:)                ! initial and (n-1)th iteration of the Lanczos vectors
    complex(mq), allocatable :: vector_1(:,:)                ! second and n-th iteration of the Lanczos vectors
    complex(mq), allocatable :: v_tmp_1(:,:), v_tmp_2(:,:)   ! auxiliary vectors 
    complex(mq) :: ortho_test(1:6,1:6)                       ! used for all dot products between vectors
#endif
    complex(q), pointer :: spectra1(:,:), spectra2(:,:)     ! spectra at different iterations 
    
    integer :: id1, id2                                      ! Cartesian directions
    integer :: it                                            ! iteration
    integer :: nomega                                        ! frequency index, number of frequencies  
    integer :: error(1:6), err_count                         ! error counter for when |<vector_0|vector_1>|==0 fails
    integer :: np, nq, cols, mycol, myrow, npcol, nprow      ! BLACS local descriptors for the grid
    real(mq) :: zz, z1                                       ! CPU time
    character(len=10), dimension(6) :: eps_dir               ! directions of epsilon 
    
    real(q) :: cput0,vput0,cput,vput
#ifdef single_prec_bse
    real(qs), parameter :: prec_dot_prod=1.d-3   ! maximum value for <vector_0|vector_1> 
#else
    real(q), parameter  :: prec_dot_prod=1.d-8   ! maximum value for <vector_0|vector_1>
#endif
    ! pointer to force hermicity in bsematrix if it has finite imaginary components in the diagonal elements
#ifndef gammareal
    complex(mq), pointer :: bsematrix_p(:) => NULL() ! pointer to force hermicity in bsematrix if it has finite imaginary components in the 
                                                     ! diagonal elements
#endif
    integer  :: matdim_node
    integer, allocatable :: matdim_rank(:)
    integer, external :: numroc
    ! we can determine the local number of rows from the descriptor of the matrix
    call blacs_gridinfo(bse_desc(CTXT_), nprow, npcol, myrow, mycol)
    matdim_node = numroc(bse_desc(M_),bse_desc(MB_),comm%node_me-1,0,nprow)
    ! used for the fancier method, where the number of elements per row done by each rank is computed from the 
    ! descriptor vector
#ifdef single_prec_bse
    allocate(matdim_rank(1:comm%ncpu+1))
    matdim_rank = 0
    matdim_rank(comm%node_me+1) = matdim_node
    call M_sum_i(comm,matdim_rank,comm%ncpu+1)
#endif

    eps_dir = [character(len=10) :: 'xx', 'yy', 'zz', 'xy', 'xz', 'yz']

    PROFILING_START("lanczos_iterative_cycle")

    allocate(vector_0(1:fulldim,1:6),vector_1(1:fulldim,1:6))
    allocate(v_tmp_1(1:matdim_node,1:6),v_tmp_2(1:fulldim,1:6))

    vector_0 = zerobse; vector_1 = zerobse
    v_tmp_1 = zerobse; ortho_test = zerobse; v_tmp_2 = zerobse

    a=0._mq; b=0._mq

    ! spectra1 and spectra2 are allocated on a sparser frequency grid that will be used for convergence
    ! frequency grid for convergence 
    nomega = int(sqrt(real(nedos_in)))
    allocate(spectra1(1:nomega,1:6),spectra2(1:nomega,1:6))
    ! the spectra is always a complex array
    spectra1 = czero; spectra2 = czero
    
    norm_v0 = 0._mq

    vector_0(1:fulldim,1:3) = optmat(1:fulldim,1:3)

#ifndef gammareal    
    ! enforce hermicity by removing spurious imaginary part in the diagonal elements
    call hermicity(fulldim, bsematrix, bse_desc)
#endif

     ! off-diagonal components xy, yz, zx, in that order
     vector_0(1:fulldim,4) = (optmat(1:fulldim,1) + optmat(1:fulldim,2))/sqrt(2._mq)  ! \epsilon_{xy}
     vector_0(1:fulldim,5) = (optmat(1:fulldim,2) + optmat(1:fulldim,3))/sqrt(2._mq)  ! \epsilon_{yz}
     vector_0(1:fulldim,6) = (optmat(1:fulldim,1) + optmat(1:fulldim,3))/sqrt(2._mq)  ! \epsilon_{zx}

    ! compute the norm along each direction and normalize the vector
    do id1 = 1, 6
#ifdef gammareal
      ! all quantities are real by default, no need to call real function
      norm_v0(id1) = sqrt(dot_product(vector_0(1:fulldim,id1),vector_0(1:fulldim,id1)))
#else
      norm_v0(id1) = sqrt(real(dot_product(vector_0(1:fulldim,id1), vector_0(1:fulldim,id1)),mq))
#endif
      vector_1(1:fulldim,id1) = vector_0(1:fulldim,id1)/norm_v0(id1)
    enddo

    !> call start time 
    if (io%iu0>=0) call vtime(vput0,cput0)

    ! starting the loop
    it = 1
    itercount = 1
    err_count = 0
    do while (it <= fulldim)

      ! temporary vectors, set to czero for precaution
      v_tmp_1 = zerobse; v_tmp_2 = zerobse
      ! |v_tmp_1> = bsematrix|vector_1>
      ! tests showed that this operation can be left as a gemm call in both single and double precision, without significantly
      ! changing the value of <0|1> and keep execution time as short as possible
      if (matdim_node>0) &
          call GGEMMbse('n', 'n', matdim_node, 6, fulldim, onebse, bsematrix, matdim_node, vector_1, fulldim, &
                        zerobse, v_tmp_1, matdim_node)
      ! for MPI v_tmp_1 (disordered) needs to be copied into vector_1 (ordered) (copy+order all to task 1, then broadcast from task 1 to all)
      call PGGEMR2Dbse(fulldim, 6, v_tmp_1, 1, 1, dvec_desc, v_tmp_2, 1, 1, vec_desc, dvec_desc(2))
      call M_bcastbse_g_from(comm, v_tmp_2, size(v_tmp_2), 1)

      ! a(it) = <vector_1|bsematrix|vector_1> = <vector_1|v_tmp_1>
      ! for single precision this operation must call Kahan
#ifndef single_prec_bse 
      call GGEMMbse(trans,'n', 6, 6, fulldim, onebse, vector_1, fulldim, v_tmp_2, fulldim, zerobse, ortho_test, 6)
#else
#ifdef gammareal
      call matmul_kahan(comm,transpose(vector_1), v_tmp_2, ortho_test, 6, fulldim, fulldim, 6, 6, 6, 't',matdim_rank)
#else
      call matmul_kahan(comm,conjg(transpose(vector_1)), v_tmp_2, ortho_test, 6, fulldim, fulldim, 6, 6, 6, 't',matdim_rank)
#endif
#endif
      a(it,1:6,1:6) = real(ortho_test(1:6,1:6),mq)

      ! v_tmp_1 <- bsematrix|vector_1> - a(it)|vector_1> - b(it-1)|v0> = |v_tmp_1> - a(it)|vector_1> - b(it-1)|vector_0>
      do id1 = 1,6
        ! while a(it) = real(<vector_1|bsematrix|vector_1>), and since bsematrix is hermitian this is not a problem at the
        ! conceptual level, numerically it is better to use <vector_1|bsematrix|vector_1> when orthogonalising the vectors
        v_tmp_2(:,id1) = v_tmp_2(:,id1) - ortho_test(id1,id1)*vector_1(:,id1)
        if (it > 1) then
          v_tmp_2(:,id1) = v_tmp_2(:,id1) - sqrt(b(it-1,id1,id1))*vector_0(:,id1)
        endif
      enddo

      ! b(it)^2 = <v_tmp_2|v_tmp_2>
#ifndef single_prec_bse
      call GGEMMbse(trans,'n', 6, 6, fulldim, onebse, v_tmp_2, fulldim, v_tmp_2, fulldim, zerobse, ortho_test, 6)
#else
#ifdef gammareal
      call matmul_kahan(comm,transpose(v_tmp_2), v_tmp_2, ortho_test, 6, fulldim, fulldim, 6, 6, 6, 't',matdim_rank)
#else
      call matmul_kahan(comm,conjg(transpose(v_tmp_2)), v_tmp_2, ortho_test, 6, fulldim, fulldim, 6, 6, 6, 't',matdim_rank)
#endif
#endif
      b(it,1:6,1:6) = real(ortho_test(1:6,1:6),mq)

      ! |vector_0> <- |vector_1>
      vector_0 = vector_1

      ! |vector_1> <- |v_tmp_1>/sqrt(b(it))
      ! old vector_1 is copied to vector_0, but now I need to reassign the distributed vector v_tmp_1 to vector_1
      ! seems easier to zero vector_1, reassign the column-stripes contained in v_tmp_1, and then sum and bcast them
      vector_1 = zerobse
      do id1=1,6
        vector_1(:,id1) = v_tmp_2(:,id1)/sqrt(b(it,id1,id1))
      enddo 

      ! orthogonality test between |vector_0> and |vector_1>, <vector_0|vector_1>
#ifndef single_prec_bse
      call GGEMMbse(trans,'n', 6, 6, fulldim, onebse, vector_0, fulldim, vector_1, fulldim, zerobse, ortho_test, 6)
#else
#ifdef gammareal
      call matmul_kahan(comm,transpose(vector_0), vector_1, ortho_test, 6, fulldim, fulldim, 6, 6, 6, 't',matdim_rank)
#else
      call matmul_kahan(comm,conjg(transpose(vector_0)), vector_1, ortho_test, 6, fulldim, fulldim, 6, 6, 6, 't',matdim_rank)
#endif
#endif

      error = 0

      ! check if there are errors and if so change error-array
      do id1 = 1,6
        ! Break in the <0|1> = 0 condition.
        if (abs(ortho_test(id1,id1)) > prec_dot_prod) then
          error(id1) = 1
          call vtutor%warning("At iteration " // str(it) // " the orthogonality condition between the two consequtive Lanczos vector, |0> and |1>, was &
               broken. It is not certain that this is recoverable, so the code will now stop! &
               Most likely this comes from having a non-hermitian BSE hamiltonian. Please try running the DFT steps with &
               LORBITALREAl=.TRUE. in the INCAR.")
        endif
      enddo
      
      ! force all ranks to sum error-array, serves as check point
      call M_sum_i(comm, error, 6) 

      ! check for errors and choose what to do next depending on IDIOT
      if (any(error(1:6) .gt. 0 )) then
        call M_sumbse_v(comm, ortho_test, 36)
        ! Standard procedure: the iterative procedure will stop now. No output is written.
        if (io%idiot > 0) then
          call vtutor%warning("At iteration " // str(it) //  " the value of check <0|1>=0 was larger than the allowed threshold. ")
          do id1 = 1, 6
            if (io%iu6>=0) write(io%iu6, "(a,a,e14.7)")  'error along ', eps_dir(id1), abs(ortho_test(id1,id1))
          enddo
          call vtutor%error("At iteration " // str(it) // " the orthogonality condition between the two consecutive Lanczos vectors, |0> and |1>, was &
                 broken. It is not certain that this is recoverable, so the code will now stop! You can try to force the code to continue runnig by & 
                 setting IDIOT=0 in the INCAR file. Alternatively, try re-running after recomputing the DFT steps & 
                 with LORBITALREAL=.TRUE. in the INCAR. If the problem persists, please contact the developers.")
          ! set it to fulldim just in case
          it = fulldim
        ! user wants the code to continue
        else
          ! 1st error occurrence: print message about the error and report values 
          if (err_count == 0) then
            call vtutor%alert("At iteration " // str(it) // " the orthogonality condition between the two consecutive Lanczos & 
                 vectors, |0> and |1>, was broken. Please try re-running after recomputing the DFT steps with LORBITALREAL=.TRUE. in the INCAR. & 
                 VASP will continue the calculation, but with no guarantees that the final result is physical. The values of |<0|1>| & 
                 will be reported in the OUTCAR.")
            if (io%iu6>=0) write(io%iu6, "(a,/,a,/,a)") ' Value of |<0|1>| along all directions when the |<0|1>|==0 check failed ',&
&             '  Iteration      X         Y         Z        XY        YZ        ZX', &
&             '  --------------------------------------------------------------------------------------------------------------'

          endif
          ! 2nd or further occurrence: print only values 
          if (io%iu6>=0) write(io%iu6, "(i6,6e14.7)") it, abs(ortho_test(1,1)),abs(ortho_test(2,2)),abs(ortho_test(3,3)),& 
&                                                         abs(ortho_test(1,2)),abs(ortho_test(1,3)),abs(ortho_test(2,3))
          err_count = err_count + 1
        endif
      endif

      if (mod(it,10) == 0) then
        ! lanczos spectra at current iteration
        call lanczos_spectra(a, b, it  , shift, 0, omegamax, nomega, spectra2, norm_v0, 0)
        ! lanczos spectra at previous iteration
        call lanczos_spectra(a, b, it-1, shift, 0, omegamax, nomega, spectra1, norm_v0, 0)
        ! rms between spectra(it)=spectra2 and spectra(it-1) = spectra1
        call spectra_convergence(spectra1, spectra2, nomega, err)        
        ! if error condition is satisfied set make the cycle break
        if (sum(err(1:3))/3._q < prec_eps) then
          it = fulldim
          ! force broadcast from node 1 to all nodes in case not all exit the loop due to some differences in precision
          call M_bcast_i_from(comm, it, 1, 1)
        endif
      endif

      ! write to OUTCAR and stdout only every 100 iterations
      if (mod(it,100) == 0) then
        ! get time
        if (io%iu0>=0) call vtime(vput,cput)
        !> write information about convergence to OUTCAR and stdout
        if (io%iu0>=0) write(io%iu0,'(a,i6,a,e14.8,a,f13.2,a)') ' Lanczos iteration ', it ,' done: RMS ', sum(err(1:3))/3._q, '; Time ', cput-cput0, ' sec'
      endif

      it = it + 1
      itercount = itercount + 1

    enddo ! it

    ! since the next iteration is not done once convergence is reached
    itercount = itercount - 1

    !> deallocate spectra and vector arrays
    deallocate(spectra1,spectra2,vector_0,vector_1,v_tmp_1,v_tmp_2)

    !> print to std out information on convergence
    if (io%iu0>=0) then
      call vtime(vput,cput)
      if (io%iu0>=0) write(io%iu0,'(a,i6,a,e14.8,a,f13.2,a)') ' Lanczos cycle converged after ', itercount ,' iterations. Trace error is ', &
&                    sum(err(1:3))/3._q, '; Time ', cput-cput0, ' sec'
    endif
    
    PROFILING_STOP("lanczos_iterative_cycle")

  end subroutine lanczos_solver_tda


!****************** subroutine lanczos_spectra    *****************************
!
!> computes the spectra from haydock-lanczos method using a continued fraction
!
!******************************************************************************

  subroutine lanczos_spectra(a, b, iter, shift, term, omegamax, nomega, spectra, u0, res)
    implicit none

    real(q), intent(in)  :: a(:,:,:)                   !< Lanczos a-coefficients
    real(q), intent(in)  :: b(:,:,:)                   !< Lanczos b-coefficients
    real(mq), intent(in) :: u0(:)                      !< norm of the initial Lanczos vector
    integer, intent(in)  :: iter                       !< maximum iteration
    integer, intent(in)  :: term                       !< terminator option
    integer, intent(in)  :: res                        !< resonant or anti-resonant spectra
    integer, intent(in)  :: nomega                     !< number of frequencies
    real(q),intent(in)   :: omegamax                   !< maximum frequency
    real(q),intent(in)   :: shift                      !< complex shift in epsilon 
    complex(q), intent(out) :: spectra(:,:)           !< epsilon
    ! workspace
    integer  :: iw, it, id1,id2                      ! frequency, iteration, and Cartesian direction indices
    real(q) :: deltaw                                ! frequency
    complex(mq) :: c1, c2, f, g                      ! terminator parameters
    complex(mq) :: p_amp, n_amp

    c1 = czero
    c2 = czero
    f = czero
    g = czero

    deltaw = omegamax/(nomega-1)
    spectra = czero
    do id1 = 1,6
      do iw = 0, nomega-1
        if (term > 0) then
          c1 = sum(a(1:iter,id1,id1))/iter
          c1 = sum(b(1:iter,id1,id1))/iter
          f = (iw*deltaw)**2 - c1**2 + c2**2
          g = 2.0*iw*deltaw*c2**2
          spectra(iw+1,id1) = (f**2 - 2.0*iw*deltaw*g)/g
        endif
        p_amp = czero; n_amp = czero
        do it = 0, iter-1, 1
          p_amp = 1.0/( iw*deltaw - a(iter-it,id1,id1) + imun*shift - b(iter-it,id1,id1)*p_amp)
          n_amp = 1.0/(-iw*deltaw - a(iter-it,id1,id1) + imun*shift - b(iter-it,id1,id1)*n_amp)
        enddo
        spectra(iw+1,id1) = p_amp + conjg(n_amp)
      enddo
      spectra(:,id1) = -u0(id1)**2*spectra(:,id1)
    enddo
  end subroutine lanczos_spectra

!****************** subroutine lanczos_spectra_eigen **************************
!
!> computes the spectra from Haydock-Lanczos method using eigenvalues
!> supports both Lorentzian and Gaussian functions
!
!******************************************************************************

  subroutine lanczos_spectra_eigen(a, itercount, lanc_vec, nomega, omegamax, u0, shift, spectra, sp_kind, ares)
    implicit none

    real(q), intent(in)  :: a(:)
    real(mq), intent(in) :: u0
    integer, intent(in)  :: itercount, nomega
    real(q),intent(in)   :: omegamax, shift
    complex(q), intent(in) :: lanc_vec(:)
    complex(q), intent(out) :: spectra(:)
    character(len=1) :: sp_kind, ares

    !> workspace
    integer  :: iw, it, iw1
    real(q)  :: deltaw, c_amp, eigen

    deltaw = omegamax/(nomega-1)
    spectra = czero

    do it = 1, itercount
      c_amp = real(lanc_vec(it)*conjg(lanc_vec(it)))
      eigen = a(it)
      if (ares == 't') then
        c_amp = c_amp/sqrt(a(it))
        eigen = sqrt(a(it))
      endif
      !> if Lorentizans are used we can build the real and imaginary parts together
      if (sp_kind == 'l') then
        do iw = 0, nomega-1
          spectra(iw+1) = spectra(iw+1) + c_amp*(1._q/( iw*deltaw - eigen + imun*shift) &
                                               - 1._q/( iw*deltaw + eigen + imun*shift))
        enddo
      !> for gaussians we need to build the imaginary part and the use Kramers-Kronig for the real part
      elseif (sp_kind == 'g') then
        do iw = 0, nomega - 1
          spectra(iw+1) = spectra(iw+1) - imun*c_amp*sqrt(pi)/shift*(exp(-(iw*deltaw - eigen)**2/(shift**2)) &
                                                                   - exp(-(iw*deltaw + eigen)**2/(shift**2)))
        enddo
      endif
    enddo

    !> if Gaussians are used, compute the real part at the end
    if (sp_kind == 'g') call kramers_kronig_bse(nomega,deltaw,spectra)
    
    !> renormalise the spectra with the optical transition matrix elements
    spectra(:) = -u0**2*spectra(:)

  end subroutine lanczos_spectra_eigen


!****************** subroutine kramers_kronig_bse    **************************
!
!> computes the real part of epsilon from it's imaginary part using the Kramers-Kronig relation
!
!******************************************************************************
  subroutine kramers_kronig_bse(nomega,deltaw,spectral)
    integer, intent(in) :: nomega
    real(q), intent(in)    :: deltaw
    complex(q), intent(inout) :: spectral(:)

    !> workspace
    integer :: iw1, iw2
    real(q) :: tmp

    do iw1 = 0, nomega - 1
      tmp = 0._q
      do iw2 = 0, nomega - 1
        if (iw2 == iw1) cycle
        tmp = tmp + 2._q/pi*deltaw*aimag(spectral(iw2+1))*iw2*deltaw/((iw2*deltaw)**2 - (iw1*deltaw)**2)
      enddo
      spectral(iw1+1) = spectral(iw1+1) + tmp
    enddo

  end subroutine kramers_kronig_bse

 
!****************** subroutine spectra_convergence    *************************
!
!> computes the RMS error between two different iterations of the spectra from haydock-lanczos
!
!******************************************************************************
  subroutine spectra_convergence(spectra1, spectra2, nomega, error)
    complex(q), intent(in) :: spectra1(:,:)                   !< spectra at previous iteration
    complex(q), intent(in) :: spectra2(:,:)                   !< spectra at current iteration
    integer, intent(in) :: nomega                              !< number of frequencies
    real(q) :: error(:)                                       !< root-mean-square between spectra1 and spectra2
    ! workspace           
    integer :: iw, id                                          ! frequency and Cartesian direction indices 

    error = 0._mq 
    do id = 1,6
      do iw = 1, nomega
        error(id) = error(id) + (aimag(spectra1(iw,id)) - aimag(spectra2(iw,id)))**2
      enddo
      error(id) = sqrt(error(id)/nomega)
    enddo

  end subroutine spectra_convergence


!****************** subroutine write_bse_lanczos    **********************************
!
!> writes the BSE fatbands from Lanczos to a text file and to the vaspout.h5 file
!
!*************************************************************************************

  subroutine write_bse_lanczos(nbseeig,itercount,eigenval,eigenvec)
    implicit none
    real(q), intent(in) :: eigenval(:), eigenvec(:,:)
    integer, intent(in) :: nbseeig, itercount
    ! local
    integer ::  i, j
#ifdef VASP_HDF5
    complex(q) :: fatbands(itercount,1)
#endif

    call vtutor%warning('The Lanczos solver can compute exciton eigenvalues, but it cannot decompose an exciton into {v,k} -> {c,k} &
    transitions. The BSEFATBANDS file and vaspout.h5 will contain the eigenvalues and their respective amplitude, but not the &
    discriminated list of transitions that form the excitonic state.')
    
    open(unit=78, file='BSEFATBAND', status='unknown')
    write(78,'(I18)') nbseeig
    do i = 1, nbseeig
      write(78,'(I6,A14,F14.8,A20,F14.6,A5,F14.6)') i, 'BSE eigenvalue:', eigenval(i), 'BSE eigenvector:', eigenvec(i,1), ' +i* ', eigenvec(i,2)
#ifdef VASP_HDF5
      fatbands(i,1) = cmplx(eigenvec(i,1),eigenvec(i,2))
#endif
    enddo      
    close(78)

#ifdef VASP_HDF5
    do i = nbseeig + 1, itercount
      fatbands(i,1) = cmplx(eigenvec(i,1),eigenvec(i,2))
    enddo
    call vh5_write_bsefatband(IH5OUTFILEID, fatbands)
#endif

  end subroutine write_bse_lanczos


!****************** subroutine lanczos_start_time   **********************************
!
!> obtains the start time of the iteration loop and prints it to standard output and OUTCAR
!
!*************************************************************************************
  subroutine lanczos_start_time(z,io)
    use base

    implicit none
    real(mq) :: z
    type (in_struct)        :: io

    if (node_me==0) call cpu_time(z)

  end subroutine lanczos_start_time  

#endif
end module bse_lanczos

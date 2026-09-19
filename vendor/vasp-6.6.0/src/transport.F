#include "symbol.inc"

module electron_transport
    use prec
    use constant
    use mymath
    use mathtools, ONLY: det_3x3,matrix_inverse_3x3
    use wave_struct_def, ONLY: wavespin, wavedes
    use mkpoints_struct_def, only: kpoints_struct
    use, intrinsic :: iso_c_binding
    use elphon_kgrid
    use poscar, ONLY: latt
    use elphon_base, ONLY: elph_settings, elph_selfen_approx_t
    use tutor, ONLY: vtutor
    use bandgap_struct, ONLY: bandgap_info
    use mpimy, only: communic
    implicit none

    real(q), parameter :: energy_0=1e-9_q
    !real(q), parameter :: energy_0=8.617333262e-5 ! 10^(-5) Kelvin
    real(q), parameter :: hbar = 6.582119569e-16_q ![eV.s]
    real(q), parameter :: transport_factor = EVTOJ / hbar / 1.e-10_q ![C^2/(J.s.m)]

    !> @brief Container for the electronic transport information
    !> it computes the transport coefficients for a fixed carrier concentration,
    !> doping and temperature
    type transport_type
        integer :: driver                                      !< driver 1 linear grids and simpson itnegration or driver 2 gauss-legendre integration
        character*3 :: boltzmann                               !< approximation to the boltzmann equation
        real(q) :: a(3,3)                                      !< lattice vectors in Ang, as columns
        real(q) :: vol                                         !< volume of unit cell in Ang^3
        integer :: nrotk                                       !< number of space group operations
        integer :: isymop(3,3,48)                              !< space group operations
        integer :: nbands                                      !< number of bands considered for transport
        integer :: nkpoints                                    !< number of points in IBZ
        integer :: nspin                                       !< number of spin index in energy bands
        integer :: rspin                                       !< spin degeneracy factor (2 for spin-unpolarized, 1 for collinear
                                                               !< and non-collinear spin)
        real(q),allocatable :: temps(:)                        !< list of temperatures at which to compute the transport function
        real(q),allocatable :: mu(:)                           !< chemical potentials for each temperature
        real(q),allocatable :: nelect(:)                                      !< number of electron per cell (including doping)
        real(q) :: nelect0                                     !< number of electron per cell (without doping)
        type(elph_selfen_approx_t) :: selfen_approx            !< approximation used to compute tau
        real(q),allocatable :: ne(:)                           !< number of electrons in the conduction bands
        real(q),allocatable :: nh(:)                           !< number of holes in the valence bands
        real(q) :: ncarrier                                    !< number of carriers per cell
        real(q) :: dfermi_tol                                  !< define a tolerance for the derivative of the fermi-dirac distribution
                                                               !  beyond which the state is considered to be unoccupied
                                                               !  this is used to determine the energy window where to compute the
                                                               !  transport function
        integer,allocatable :: id_idx(:)                       !< indexes for which this accumulator was generated
        character(len=15),allocatable :: id_name(:)             !< name of the fields for which this accumulator was generated

        type(bandgap_info),allocatable :: bandgap(:)           !< information about the bandgaps

        integer :: nedos                                       !< number of energy grid points
        real(q),allocatable :: emin(:)                         !< minimun energy for transport function
        real(q),allocatable :: emax(:)                         !< maximum energy for transport function
        logical :: do_plot                                     !< compute the transport function for plotting
        logical :: do_hall                                     !< compute hall coefficients
        integer :: nedos_plot                                  !< number of energy grid points
        real(q) :: emin_plot                                   !< minimun energy for plotting the transport function
        real(q) :: emax_plot                                   !< maximum energy for plotting the transport function
        real(q) :: emin_userinput                              !< minimun energy for integrating the transport function in driver 1
        real(q) :: emax_userinput                              !< maximum energy for integrating the transport function in driver 1
        real(q) :: alpha                                       !< exponent of energy grid
        real(q),allocatable :: energy(:,:)                     !< energy grid
        real(q),allocatable :: energy_plot(:)                  !< energy grid for plotting the transport function
        real(q),pointer :: velocity(:,:,:,:)                   !< electron group velocity (nbands,nkpoints,nspin,dir)
                                                               !  only a reference, the memory must be allocated somewhere else
        real(q),pointer :: effective_mass(:,:,:,:,:)           !< electron effective masses (nbands,nkpoints,nspin,dir)
                                                               !  only a reference, the memory must be allocated somewhere else
        real(q),pointer :: tau(:,:,:)                          !< electron lifetime (doping,temperature,nbands,nkpoints,nspin)
                                                               !  only a reference, the memory must be allocated somewhere else
        real(q),allocatable :: tau_average(:)                              !< compute an average tau around mu (to compare with CRTA calculations)
        real(q),allocatable :: mean_free_path(:,:,:,:)         !< (nbands,nkpoints,nspin,dir)
        real(q),allocatable :: mean_free_path_mag(:,:,:,:,:)   !< (nbands,nkpoints,nspin,dir,dir)
        real(q),allocatable :: transport_function(:,:,:,:,:)     !< (nedos,3,3,2,ntemps)
        real(q),allocatable :: transport_function_e(:,:,:,:,:)   !< (nedos,3,3,2,ntemps)
        real(q),allocatable :: transport_function_h(:,:,:,:,:)   !< (nedos,3,3,2,ntemps)
        real(q),allocatable :: transport_function_plot(:,:,:,:,:)  !< (nedos_plot,3,3,2,ntemps)
        real(q),allocatable :: transport_function_hall(:,:,:,:,:,:)     !< (nedos,3,3,3,2,ntemps)
        real(q),allocatable :: transport_function_hall_plot(:,:,:,:,:)  !< (nedos_plot,3,3,3,2,ntemps)
        real(q),allocatable :: lab(:,:,:,:,:)                                !< Onsager coefficients
        real(q),allocatable :: lab_e(:,:,:,:,:)                              !< Onsager coefficients for electrons
        real(q),allocatable :: lab_h(:,:,:,:,:)                              !< Onsager coefficients for holes
        real(q),allocatable :: lab_hall(:,:,:,:)
        real(q),allocatable :: e_conductivity(:,:,:)                         !< electronic conductivity (ohm m)^{-1}
        real(q),allocatable :: e_conductivity_e(:,:,:)                       !< electronic conductivity (ohm m)^{-1}
        real(q),allocatable :: e_conductivity_h(:,:,:)                       !< electronic conductivity (ohm m)^{-1}
        real(q),allocatable :: mobility(:,:,:)                               !< electron or hole mobility cm^2/(V.s)
        real(q),allocatable :: mobility_e(:,:,:)                             !< electron mobility cm^2/(V.s)
        real(q),allocatable :: mobility_h(:,:,:)                             !< hole mobility cm^2/(V.s)
        real(q),allocatable :: e_t_conductivity(:,:,:)                       !< electronic thermal conductivity W/(m.K)
        real(q),allocatable :: seebeck(:,:,:)                                !< Seebeck coefficient V/K
        real(q),allocatable :: peltier(:,:,:)                                !< Peltier coefficient V
        real(q),allocatable :: hall_conductivity(:,:,:,:)                    !< hall conductivity
    end type

    !> Container for chemical potential computed at different temperatures and doping levels
    type chemical_potential
        real(q),allocatable :: selfen_carrier_den(:) !< carrier density is used to specify the number of carriers per cell
        real(q),allocatable :: selfen_carrier_per_cell(:) !< number of carriers per cell
        real(q),allocatable :: selfen_mu(:) !< chemical potential is fixed so we determine the number of carriers per cell
        real(q) :: efermi0 !< fermi level at temperature indicated by the SIGMA parameter
        real(q),allocatable :: temperatures(:) !< list of temperatures
        real(q),allocatable :: selfen_muij(:,:) !< chemical potential in eV for each (temperature,carrier)
        real(q),allocatable :: selfen_nij(:,:) !< number of electrons per unit cell for each (temperature,carrier)
    end type chemical_potential

    contains

    !> Initialize the transport datastructure
    subroutine transport_init(self,a,nrotk,isymop,temps,nelect,nelect0,selfen_approx,rspin,mu,emin_plot,emax_plot,&
                              emin_userinput,emax_userinput,nedos,nedos_plot,dfermi_tol,bandgap,id_idx,id_name,tau,velocity)
        use string, only : str
        use bandgap_struct, only: bandgap_info
        type(transport_type) :: self
        integer,intent(in) :: nrotk !< number of space group operations
        integer,intent(in) :: isymop(:,:,:) !< space group operations
        real(q),intent(in) :: a(3,3) !< lattice vectors
        real(q),intent(in) :: temps(:) !< temperatures at which to compute the transport coefficients (K)
        real(q),intent(in) :: nelect(:)  !< number of electrons per unit cell (including doping)
        real(q),intent(in) :: nelect0 !< number of electrons per unit cell (without doping)
        type(elph_selfen_approx_t), intent(in) :: selfen_approx !< approximation used to compute the lifetimes
        real(q),intent(in) :: rspin !< spin degeneracy factor
        real(q),intent(in) :: mu(:) !< chemical potential
        integer,intent(in) :: nedos !< number of energy points at which to evaluate the transport function
        integer,intent(in) :: nedos_plot !< number of energy points at which to evaluate the transport function
        real(q),intent(in) :: emin_plot !< min energy points at which to evaluate the transport function for plotting
        real(q),intent(in) :: emax_plot !< max energy points at which to evaluate the transport function for plotting
        real(q),intent(in) :: emin_userinput !< min energy points at which to integrate the transport function in driver 1
        real(q),intent(in) :: emax_userinput !< max energy points at which to integrate the transport function in driver 1
        real(q),intent(in) :: dfermi_tol !< tolerance beyond which the derivative of the fermi-dirac function is considered zero
        type(bandgap_info),intent(in) :: bandgap(:)
        integer,intent(in) :: id_idx(:) !< indexes for which this accumulator was generated
        character(len=*),intent(in) :: id_name(:) !< name of the fields for which this accumulator was generated
        real(q),target,intent(in) :: tau(:,:,:)
        real(q),target,intent(in) :: velocity(:,:,:,:)
        !local variables
        integer :: ntemps

        self%a(:,:) = a(:,:)
        ! cell volume
        self%vol = abs(det_3x3(self%a))

        self%nrotk = nrotk
        ! in vasp the rotation matrices of rotations are defined as r'_j = \sum_i ISYMOP(i,j, ..) r_i
        self%isymop(:,:,:) = isymop(:,:,:)

        self%temps=temps
        self%nelect=nelect
        self%nelect0=nelect0
        self%mu=mu
        self%selfen_approx=selfen_approx

        self%rspin = rspin
        self%bandgap = bandgap

        ! generate energy grids for computing transport coefficients
        self%nedos=nedos
        self%dfermi_tol=dfermi_tol

        self%id_idx=id_idx
        self%id_name=id_name

        !self%emin !set in transport_compute_kgrid or transport_compute
        !self%emax !set in transport_compute_kgrid or transport_compute
        !self%energy !set in transport_compute_kgrid or transport_compute

        ! generate energy grids for for plotting transport coefficients
        self%do_plot = nedos_plot>0
        self%do_hall = .false.
        self%emin_plot=emin_plot
        self%emax_plot=emax_plot
        !set in transport_compute_kgrid or transport_compute
        !self%energy_plot=linspace(self%emin_plot,self%emax_plot,self%nedos_plot)
        self%nedos_plot=nedos_plot
        self%emin_userinput=emin_userinput
        self%emax_userinput=emax_userinput

        self%nspin=size(tau,3)
        self%nkpoints=size(tau,2)
        self%nbands=size(tau,1)
        allocate(self%mean_free_path(self%nbands,self%nkpoints,self%nspin,3))

        ntemps = size(self%temps)
        allocate(self%energy(self%nedos,ntemps))
        allocate(self%lab(2,2,3,3,ntemps))
        allocate(self%lab_e(2,2,3,3,ntemps))
        allocate(self%lab_h(2,2,3,3,ntemps))
        allocate(self%lab_hall(3,3,3,ntemps))
        allocate(self%e_conductivity(3,3,ntemps))
        allocate(self%e_conductivity_e(3,3,ntemps))
        allocate(self%e_conductivity_h(3,3,ntemps))
        allocate(self%mobility(3,3,ntemps))
        allocate(self%mobility_e(3,3,ntemps))
        allocate(self%mobility_h(3,3,ntemps))
        allocate(self%e_t_conductivity(3,3,ntemps))
        allocate(self%seebeck(3,3,ntemps))
        allocate(self%peltier(3,3,ntemps))

        allocate(self%ne(ntemps))
        allocate(self%nh(ntemps))

        allocate(self%emin(ntemps))
        allocate(self%emax(ntemps))

        allocate(self%transport_function(self%nedos,3,3,self%nspin,ntemps))
        allocate(self%transport_function_e(self%nedos,3,3,self%nspin,ntemps))
        allocate(self%transport_function_h(self%nedos,3,3,self%nspin,ntemps))
        allocate(self%transport_function_hall(self%nedos,3,3,3,self%nspin,ntemps))
        if (self%do_plot) allocate(self%transport_function_plot(self%nedos_plot,3,3,self%nspin,ntemps))

        if (self%nbands /= size(velocity,1)) then
           call vtutor%bug('Wrong dimensions in first dimension of the velocity array got '//str(size(velocity,1))//&
                           ' expected '//str(self%nbands),__FILE__,__LINE__)
        endif

        if (self%nkpoints /= size(velocity,2)) then
           call vtutor%bug('Wrong dimensions in second dimension of the velocity array got '//str(size(velocity,2))//&
                           ' expected '//str(self%nkpoints),__FILE__,__LINE__)
        endif

        if (self%nspin /= size(velocity,3)) then
           call vtutor%bug('Wrong dimenions in third dimension of the velocity array got '//str(size(velocity,3))//&
                           ' expected '//str(self%nspin),__FILE__,__LINE__)
        endif

        self%tau=>tau
        allocate(self%tau_average(ntemps))
        self%tau_average = huge(1.0_q)
        self%velocity=>velocity
    end subroutine transport_init

    !> Compute average tau weighed by the derivative of the fermi-dirac occupation factor
    !> This is only to display to OUTCAR and write to the hdf5 file
    subroutine transport_compute_tau_average(self,celtot,itemp)
        use elphon_base, only: dfermi_dirac
        type(transport_type) :: self
        complex(q),intent(in) :: celtot(:,:,:)
        integer,intent(in) :: itemp !< index of the temperature which we are computing
        ! local variables
        integer :: isp, ikpt, ib
        real(q) :: weight, total_weight, eig
        total_weight = 0.0_q
        self%tau_average(itemp) = 0.0_q
        ! TODO HM: find a better function to average tau: the states around the fermi level
        ! have a larger contribution according to the derivative of the fermi-dirac distribution
        ! so this should be accounted for. I think one should compute a transport relaxation time density of states
        ! with delta evaluated using the tetrahedron method:
        ! tau_average = sum_nk tau_nk dfermi(enk-mu,T) delta(enk-mu)
        ! I will not implement this for now so instead I evaluate:
        ! tau_average = ( sum_nk tau_nk dfermi(enk-mu,T) ) / ( sum_nk dfermi(enk-mu,T) )
        ! the point is that when CRTA is used tau_average = TRANSPORT_RELAXATION_TIME and its written to the hdf5 file
        do isp=1,self%nspin
            do ikpt=1,self%nkpoints
                do ib=1,self%nbands
                    eig = real(celtot(ib,ikpt,isp),q)
                    if (self%temps(itemp)*BOLKEV < 1e-12) then
                        weight = dfermi_dirac(eig-self%mu(itemp),1.0_q*BOLKEV)
                    else
                        weight = dfermi_dirac(eig-self%mu(itemp),self%temps(itemp)*BOLKEV)
                    endif
                    self%tau_average(itemp) = self%tau_average(itemp) + self%tau(ib,ikpt,isp)*weight
                    total_weight = total_weight + weight
                enddo
            enddo
        enddo
        if (total_weight>tiny(total_weight)) then
            self%tau_average(itemp) = self%tau_average(itemp)/total_weight
        endif
    end subroutine transport_compute_tau_average

    !> Now compute the transport quantities
    subroutine transport_compute(self,itemp,wdes_dense_ibz,w_dense_ibz,kpoints,comm)
        use wave_struct_def, only: wavespin, wavedes
        use base, only: symmetry
        use string, only: str
        use elphon_base, only: transport_emin_emax, linspace
        type(transport_type) :: self
        integer,intent(in) :: itemp
        type(wavedes), intent(in) :: wdes_dense_ibz
        type(wavespin), intent(in) :: w_dense_ibz
        type(kpoints_struct), intent(in) :: kpoints
        type(communic), intent(in) :: comm
        ! local variables
        integer :: i, iv
        integer :: ivbm, icbm
        real(q), allocatable :: vov(:,:,:,:,:)
        real(q) :: vbm, cbm, midgap
        real(q) :: emin,emax,df

        self%driver = 1
        PROFILING_START('transport_compute')

        ! lambda = v * tau
        call mean_free_path(self)

        ! symmetrized product of mean free path and velocity
        call sym_tensor_prod(self%mean_free_path,self%velocity,self%a,self%nrotk,self%isymop,vov,comm)

        ! for plotting
        if (self%do_plot) then
            call find_energy_grid_linear(w_dense_ibz, self%nedos_plot, self%emin_plot, self%emax_plot, self%energy_plot)
            call compute_transport_function_tetrahedron(self%energy_plot,self%transport_function_plot(:,:,:,:,itemp),kpoints,w_dense_ibz%celtot,wdes_dense_ibz,&
                                                        vov,self%vol,comm)
        endif

        ! prepare linear energy grid for computation of transport coefficients
        call transport_emin_emax(self%mu(itemp),self%temps(itemp)*BOLKEV,self%dfermi_tol,self%emin_userinput,self%emax_userinput,self%emin(itemp),self%emax(itemp))
        if (self%temps(itemp)>tiny(self%temps(itemp))) then
            self%energy(:,itemp) = linspace(self%emin(itemp),self%emax(itemp),self%nedos)
        else
            !self%nedos=1
            self%energy(:,itemp) = self%mu(itemp)
        endif

        ! for calculation of transport coefficients
        call compute_transport_function_tetrahedron(self%energy(:,itemp),self%transport_function(:,:,:,:,itemp),kpoints,w_dense_ibz%celtot,wdes_dense_ibz,&
                                                    vov,self%vol,comm)
        deallocate(vov)

        ! here instead of recomputing by limiting the bands we split the transport function into contributions from the conduction and valence
        cbm = self%bandgap(1)%fundamental_conduction%eigenvalue
        vbm = self%bandgap(1)%fundamental_valence%eigenvalue
        midgap = vbm+(cbm-vbm)/2
        call transport_function_cb_vb(self%energy(:,itemp),self%transport_function(:,:,:,:,itemp),midgap,self%transport_function_e(:,:,:,:,itemp),self%transport_function_h(:,:,:,:,itemp))

        ! compute number of electrons and holes
        ivbm = self%bandgap(1)%fundamental_valence%band
        icbm = self%bandgap(1)%fundamental_conduction%band
        self%nh(itemp) = number_of_electrons_lineargrid(kpoints,wdes_dense_ibz,w_dense_ibz,self%mu(itemp),self%energy(:,itemp),self%temps(itemp),comm,nbstart=1,nbstop=ivbm)
        self%nh(itemp) = self%nelect0-self%nh(itemp)
        self%ne(itemp) = number_of_electrons_lineargrid(kpoints,wdes_dense_ibz,w_dense_ibz,self%mu(itemp),self%energy(:,itemp),self%temps(itemp),comm,nbstart=icbm,nbstop=wdes_dense_ibz%nb_tot)

        ! get Onsager coefficients in (J/C)^(i+j-2)* C^2/(m s J) and remenber that (J/C = eV)
        call get_onsager_coef_uniform(self%temps(itemp),self%mu(itemp),self%transport_function(:,:,:,:,itemp),self%lab(:,:,:,:,itemp),self%energy(:,itemp))
        call get_onsager_coef_uniform(self%temps(itemp),self%mu(itemp),self%transport_function_h(:,:,:,:,itemp),self%lab_h(:,:,:,:,itemp),self%energy(:,itemp))
        call get_onsager_coef_uniform(self%temps(itemp),self%mu(itemp),self%transport_function_e(:,:,:,:,itemp),self%lab_e(:,:,:,:,itemp),self%energy(:,itemp))

        ! get transport coefficients in S.I to be compared to experiments
        call get_transport_coef(self,self%temps(itemp),self%lab(:,:,:,:,itemp),&
                                self%e_conductivity(:,:,itemp),self%seebeck(:,:,itemp),&
                                self%peltier(:,:,itemp),self%e_t_conductivity(:,:,itemp))
        call get_mobility(self,self%lab_e(:,:,:,:,itemp),self%ne(itemp),self%mobility_e(:,:,itemp))
        call get_mobility(self,self%lab_h(:,:,:,:,itemp),self%nh(itemp),self%mobility_h(:,:,itemp))
        call set_mobility(self,self%nelect(itemp),self%mobility_e(:,:,itemp),self%mobility_h(:,:,itemp),self%mobility(:,:,itemp))

       PROFILING_STOP('transport_compute')
    end subroutine transport_compute

    subroutine transport_compute_hall_gausslegendre(self,itemp,wdes_dense_ibz,w_dense_ibz,kpoints,comm)
        type(transport_type) :: self
        integer,intent(in) :: itemp
        type(wavedes), intent(in) :: wdes_dense_ibz
        type(wavespin), intent(in) :: w_dense_ibz
        type(kpoints_struct), intent(in) :: kpoints
        type(communic), intent(in) :: comm
        ! local variables
        real(q), allocatable :: vovcrosseffmass(:,:,:,:,:,:)
        real(q),allocatable :: energy(:),gauss_weights(:)

        ! lambda = v x effective_mass * tau**2
        call mean_free_path_mag(self)

        call sym_tensor3_prod(self%mean_free_path_mag,self%velocity,self%a,self%nrotk,self%isymop,vovcrosseffmass,comm)

         ! for plotting
        if (self%do_plot) then
            call compute_transport_function_hall_tetrahedron(self%energy_plot,self%transport_function_hall_plot,kpoints,w_dense_ibz,wdes_dense_ibz,&
                                                             vovcrosseffmass,self%vol,comm)
        endif

        ! prepare gauss-legendre energy grid for computation of transport coefficients
        call find_energy_grid_gausslegendre(self%temps(itemp),self%nedos,energy,gauss_weights)
        self%energy(:,itemp) = energy + self%mu
        self%emin = minval(self%energy)
        self%emax = maxval(self%energy)

        ! for calculation of transport coefficients
        call compute_transport_function_hall_tetrahedron(self%energy(:,itemp),self%transport_function_hall(:,:,:,:,:,itemp),kpoints,w_dense_ibz,wdes_dense_ibz,&
                                                         vovcrosseffmass,self%vol,comm)
        deallocate(vovcrosseffmass)

        ! compute hall conductivity
        call get_hall_conductivity_gausslegendre(self,itemp,energy,gauss_weights)

    end subroutine transport_compute_hall_gausslegendre

    !> Compute onsager coefficients using gauss legendre integration
    subroutine get_hall_conductivity_gausslegendre(self,itemp,energies,gauss_weights)
        type(transport_type) :: self
        integer,intent(in) :: itemp
        real(q),intent(in) :: energies(:)
        real(q),intent(in) :: gauss_weights(:)
        ! local variables
        integer :: i,j,k
        integer :: isp
        real(q) :: inv_e_conductivity(3,3)
        PROFILING_START('get_hall_conductivity_gausslegendre')
        ! compute onsager coefficients
        self%lab_hall = 0.0_q
        do i=1,3
          do j=1,3
            do k=1,3
              do isp=1,self%nspin
                  self%lab_hall(i,j,k,itemp) = &
                  self%lab_hall(i,j,k,itemp) +0.5_q*sum(gauss_weights(:)*self%transport_function_hall(:,i,j,k,isp,itemp))
              enddo
            enddo
          enddo
        enddo

        ! compute hall conductivity
        inv_e_conductivity = matrix_inverse_3x3(self%e_conductivity)
        do k=1,3
            self%hall_conductivity(:,:,k,itemp) = matmul(inv_e_conductivity(:,:),matmul(self%lab_hall(:,:,k,itemp),inv_e_conductivity(:,:)))
        enddo
        PROFILING_STOP('get_hall_conductivity_gausslegendre')
    end subroutine get_hall_conductivity_gausslegendre

    !> Allocate transport function and compute it using tetrahedron method
    subroutine compute_transport_function_tetrahedron(energies,transport_function,kpoints,celtot,wdes_dense_ibz,&
                                                      vov,vol,comm,nbstart,nbstop)
        real(q),intent(in) :: energies(:)
        real(q),target,intent(out) :: transport_function(:,:,:,:)
        type(kpoints_struct),intent(in) :: kpoints
        type(wavedes), intent(in) :: wdes_dense_ibz
        complex(q), intent(in) :: celtot(:,:,:)
        real(q), target, intent(in) :: vov(:,:,:,:,:)
        real(q), intent(in) :: vol !< volume of the unit cell
        type(communic), intent(in) :: comm
        integer,optional,intent(in) :: nbstart
        integer,optional,intent(in) :: nbstop
        ! local variables
        integer :: nspin,nkpoints,nbands
        real(q),pointer :: dos(:,:,:)
        real(q),pointer :: matrix_elements(:,:,:,:)

        nbands = size(vov,3)
        nkpoints = size(vov,4)
        nspin = size(vov,5)
        call c_f_pointer(c_loc(vov),matrix_elements,[9,nbands,nkpoints,nspin])
        ! get 1/(V0 N) \sum_{kn} [tau_{kn}/hbar] [hbar v_{kn}] [hbar v_{kn}] delta(e-e_{kn}) in 1/Ang.
        !allocate(transport_function(size(energies), 3, 3, nspin))
        call c_f_pointer(c_loc(transport_function),dos,[size(energies), 9, nspin])
        !call compute_transport_function_tetrahedron_fbz(kpoints,wdes_dense_ibz,celtot,self%energy_plot,vov,self%transport_function_plot)
        !call compute_transport_function_tetrahedron_ibz(kpoints,wdes_dense_ibz,celtot,latt_cur,self%energy_plot,vov,self%transport_function_plot)
        call tetrahedron_dos(kpoints,wdes_dense_ibz,celtot,energies,matrix_elements,dos,comm,nbstart,nbstop)
        ! get transport function in C^2/(m s J)=( ohm m )^-1 if v is in eV Ang
        ! need to multiply by e^2/hbar to get a conductivity
        transport_function=transport_function*transport_factor/vol
    end subroutine compute_transport_function_tetrahedron

    !> Allocate transport function and compute it using tetrahedron method
    subroutine compute_transport_function_hall_tetrahedron(energies,transport_function_hall,kpoints,w_dense_ibz,wdes_dense_ibz,&
                                                           vovcrosseffmass,vol,comm)
        real(q),intent(in) :: energies(:)
        real(q),target,intent(out) :: transport_function_hall(:,:,:,:,:)
        type(kpoints_struct),intent(in) :: kpoints
        type(wavedes), intent(in) :: wdes_dense_ibz
        type(wavespin), intent(in) :: w_dense_ibz
        real(q), target, intent(in) :: vovcrosseffmass(:,:,:,:,:,:)
        real(q), intent(in) :: vol !< volume of the unit cell
        type(communic) :: comm
        ! local variables
        integer :: nspin,nkpoints,nbands
        real(q),pointer :: dos(:,:,:)
        real(q),pointer :: matrix_elements(:,:,:,:)

        nbands = size(vovcrosseffmass,4)
        nkpoints = size(vovcrosseffmass,5)
        nspin = size(vovcrosseffmass,6)
        call c_f_pointer(c_loc(vovcrosseffmass),matrix_elements,[18,nbands,nkpoints,nspin])
        ! get 1/(V0 N) \sum_{kn} [tau_{kn}/hbar] [hbar v_{kn}] [hbar v_{kn}] delta(e-e_{kn}) in 1/Ang.
        !allocate(transport_function_hall(size(energies), 3, 3, 3, wdes_dense_ibz%ispin))
        call c_f_pointer(c_loc(transport_function_hall),dos,[size(energies), 18, wdes_dense_ibz%ispin])
        !call compute_transport_function_tetrahedron_fbz(kpoints,wdes_dense_ibz,w_dense_ibz,self%energy_plot,vov,self%transport_function_plot)
        !call compute_transport_function_tetrahedron_ibz(kpoints,wdes_dense_ibz,w_dense_ibz,latt_cur,self%energy_plot,vov,self%transport_function_plot)
        call tetrahedron_dos(kpoints,wdes_dense_ibz,w_dense_ibz%celtot,energies,matrix_elements,dos,comm)
        ! get transport function in C^2/(m s J)=( ohm m )^-1 if v is in eV Ang
        ! need to multiply by e^2/hbar to get a conductivity
        transport_function_hall=transport_function_hall*transport_factor/vol
    end subroutine compute_transport_function_hall_tetrahedron

    !> Compute the DOS weigthed by a matrix elemnt using the trahedron method
    subroutine tetrahedron_dos(kpoints,wdes_dense_ibz,celtot,energy,&
                               matrix_elements,dos,comm,nbstart,nbstop)
        USE mkpoints, ONLY: kpoints_struct
        USE tet, ONLY: deltas_single_energy
        USE tet_macdonald, ONLY: sort_etetra
        !use tet_macdonald, only: tetralite_theta
        type(kpoints_struct), intent(in) :: kpoints
        type(wavedes), intent(in) :: wdes_dense_ibz
        complex(q), intent(in) :: celtot(:,:,:)
        real(q), intent(in) :: energy(:) !< energy at which to compute the transport function
        real(q), intent(in) :: matrix_elements(:,:,:,:)
        real(q), intent(out) :: dos(:,:,:)
        type(communic), intent(in) :: comm
        integer,optional,intent(in) :: nbstart
        integer,optional,intent(in) :: nbstop
        ! local variables
        logical :: skip_hall !< skip the computation of the hall effect
        integer :: i, j, k, isp, ib, itet, iw
        integer :: kindexes(4), kindexes_sorted(4), perm(4)
        integer :: nbstart_, nbstop_
        real(q) :: emin, emax
        real(q) :: eig(4)
        real(q) :: mels(4)
        real(q) :: kweights(4,size(energy))

        PROFILING_START('tetrahedron_dos')
        nbstart_ = 1;                    if (present(nbstart)) nbstart_ = nbstart
        nbstop_ = wdes_dense_ibz%nb_tot; if (present(nbstop))  nbstop_ = nbstop
        emin = minval(energy)
        emax = maxval(energy)
        dos=0.0_q
        ! loop over tetrahedrons
        do itet=1,kpoints%ntet
#ifdef MPI
           if (mod(itet-1,comm%ncpu).ne.comm%node_me-1) cycle
#endif
           kindexes = kpoints%idtet(1:,itet)
           do isp=1,wdes_dense_ibz%ispin
              do ib=nbstart_,nbstop_
                 eig(1) = real(celtot(ib,kindexes(1),isp),q)
                 eig(2) = real(celtot(ib,kindexes(2),isp),q)
                 eig(3) = real(celtot(ib,kindexes(3),isp),q)
                 eig(4) = real(celtot(ib,kindexes(4),isp),q)

                 !early skip
                 if (all(eig>emax).or.all(eig<emin)) cycle

                 perm = [1,2,3,4]
                 call sort_etetra(eig,perm)
                 kindexes_sorted(1) = kindexes(perm(1))
                 kindexes_sorted(2) = kindexes(perm(2))
                 kindexes_sorted(3) = kindexes(perm(3))
                 kindexes_sorted(4) = kindexes(perm(4))

                 ! compute the weights for all the energies
                 ! rspin=2 for collinear non spin polarized, 1 otherwise
                 ! idtet(0,itet) is the degeneracy of the tetra
                 !call tetralite_delta(kpoints%volwgt*wdes_dense_ibz%rspin*kpoints%idtet(0,itet),eig,size(energy),energy,alpha,intdos(:,isp))
                 do iw=1,size(energy)
                    call deltas_single_energy(kpoints%volwgt*wdes_dense_ibz%rspin*kpoints%idtet(0,itet),energy(iw),eig,kweights(:,iw))
                 enddo

                 ! loop over matrix elements and accumulate
                 do i=1,size(matrix_elements,1)
                    mels(1) = matrix_elements(i, ib, kindexes_sorted(1), isp)
                    mels(2) = matrix_elements(i, ib, kindexes_sorted(2), isp)
                    mels(3) = matrix_elements(i, ib, kindexes_sorted(3), isp)
                    mels(4) = matrix_elements(i, ib, kindexes_sorted(4), isp)
                    ! accumulate contribution from this matrix element
                    do iw=1,size(energy)
                       dos(iw, i, isp) = &
                       dos(iw, i, isp) + sum(kweights(:,iw)*mels)
                    enddo
                 enddo
              enddo
           enddo
        enddo
        CALLMPI( M_sum_d(comm, dos, size(dos)) )
        PROFILING_STOP('tetrahedron_dos')
    end subroutine tetrahedron_dos

    !> Compute transport function by summing the contribution of the corner of each of the 24 tetrahedron
    !> by looping over all the k-points in the ibz
    subroutine compute_transport_function_tetrahedron_ibz(kpoints,wdes_dense_ibz,celtot,latt_cur,energy,vov,transport_function)
        use mkpoints, only: get_indexes, get_tetdisp_24tet, kpoint_integers
        use tet_macdonald, only: weight_24tetra_delta
        type(kpoints_struct), intent(in) :: kpoints
        type(latt), intent(in) :: latt_cur
        type(wavedes), intent(in) :: wdes_dense_ibz
        complex(q), intent(in) :: celtot(:,:,:)
        real(q), intent(in) :: energy(:)
        real(q), intent(in) :: vov(:,:,:,:,:)
        real(q), intent(inout) :: transport_function(:,:,:,:)
        ! local variables
        integer :: i, j
        integer :: ik, isp, itet, ib
        integer :: ic, ikibz
        integer :: kcut24(3,4,24) !< displacements of 24 tetrahedra that contribute to a k-point
        integer :: idx_k(4,24)
        integer :: kint(3), origin(24)
        real(q) :: eleig_tetra(4,24)

        PROFILING_START('transport_function_tetrahedron_ibz')

        transport_function = 0.0_q
        ! get tetrahedron displacements
        call get_tetdisp_24tet(kpoints%b,kcut24)
        ! transform kcut
        do itet=1,24
          do ic=1,4
            kcut24(:,ic,itet) = matmul(kpoints%ilsnf,kcut24(:,ic,itet))
          enddo
        enddo
        ! get [0,0,0] point of each kcut
        do itet=1,24
            do ic=1,4
                if( all(kcut24(:,ic,itet)==[0,0,0])) exit
            enddo
            origin(itet)=ic
        enddo

        do ik=1,kpoints%nkpts
            ! find the index of the k  point
            call kpoint_integers(kpoints,latt_cur,kpoints%vkpt(:,ik),kint)
            ! find indexes on the corners of the 24 tetrahedra
            call get_indexes(kpoints,kint,kcut24,idx_k)
            do isp=1,wdes_dense_ibz%ispin
                do ib=1,wdes_dense_ibz%nb_tot
                   do itet=1,24
                       do ic=1,4
                           ikibz = idx_k(ic,itet)
                           eleig_tetra(ic,itet) = real(celtot(ib,ikibz,isp),q)
                       enddo
                   enddo !itet
                   do j=1,3
                       do i=1,3
                           transport_function(:,i,j,isp) = &
                           transport_function(:,i,j,isp) + &
                           weight_24tetra_delta(kpoints%volwgt * wdes_dense_ibz%rspin,origin,size(energy),energy,eleig_tetra)*&
                           kpoints%wtkpt(ik)*kpoints%nkpx*kpoints%nkpy*kpoints%nkpz*vov(i, j, ib, ik, isp)
                       enddo
                   enddo
                enddo !ib
            enddo !isp
        enddo !kpoints
        PROFILING_STOP('transport_function_tetrahedron_ibz')
    end subroutine compute_transport_function_tetrahedron_ibz

    !> Compute transport function by summing the contribution of the corner of each of the 24 tetrahedron
    !> by looping over all the k-points in the fbz
    subroutine compute_transport_function_tetrahedron_fbz(kpoints,wdes_dense_ibz,celtot,energy,vov,transport_function)
        use mkpoints, only: get_indexes, get_tetdisp_24tet, kint_mod, kint_kindex
        use tet_macdonald, only: weight_24tetra_delta
        type(kpoints_struct), intent(in) :: kpoints
        type(wavedes), intent(in) :: wdes_dense_ibz
        complex(q), intent(in) :: celtot(:,:,:)
        real(q), intent(in) :: energy(:)
        real(q), intent(in) :: vov(:,:,:,:,:)
        real(q), intent(inout) :: transport_function(:,:,:,:)
        ! local variables
        integer :: ik,jk,kk,itet
        integer :: ib, isp
        integer :: ic, ikibz, ikfbz, i, j
        integer :: kcut24(3,4,24) !< displacements of 24 tetrahedra that contribute to a k-point
        integer :: idx_k(4,24)
        integer :: kint(3), origin(24)
        real(q) :: eleig_tetra(4,24)

        PROFILING_START('transport_function_tetrahedron_fbz')
        transport_function = 0.0_q
        ! get tetrahedron displacements
        call get_tetdisp_24tet(kpoints%b,kcut24)
        ! transform kcut
        do itet=1,24
          do ic=1,4
            kcut24(:,ic,itet) = matmul(kpoints%ilsnf,kcut24(:,ic,itet))
          enddo
        enddo
        ! get [0,0,0] point of each kcut
        do itet=1,24
            do ic=1,4
                if( all(kcut24(:,ic,itet)==[0,0,0])) exit
            enddo
            origin(itet)=ic
        enddo

        ! Compute the transport function by looping over all the k-points in the FBZ
        do ik=1,kpoints%nkpx
        do jk=1,kpoints%nkpy
        do kk=1,kpoints%nkpz
            ! generate the index of the k  point
            kint = [ik,jk,kk]
            ikfbz = kint_kindex(kpoints,kint)
            ikibz = kpoints%fbz2ibz(1,ikfbz)
            ! find indexes on the corners of the 24 tetrahedra
            call get_indexes(kpoints,kint,kcut24,idx_k)
            do isp=1,wdes_dense_ibz%ispin
                do ib=1,wdes_dense_ibz%nb_tot
                   do itet=1,24
                       do ic=1,4
                           eleig_tetra(ic,itet) = real(celtot(ib,idx_k(ic,itet),isp),q)
                       enddo
                   enddo !ITET
                   do j=1,3
                       do i=1,3
                           transport_function(:,i,j,isp) = &
                           transport_function(:,i,j,isp) + &
                           weight_24tetra_delta(kpoints%volwgt * wdes_dense_ibz%rspin,origin,size(energy),energy,eleig_tetra)*&
                           vov(i, j, ib, ikibz, isp)
                       enddo
                   enddo
                enddo !ib
            enddo !isp
        enddo !kpoints
        enddo !kpoints
        enddo !kpoints
        PROFILING_STOP('transport_function_tetrahedron_fbz')

    end subroutine compute_transport_function_tetrahedron_fbz

    !> Find the chemical potential mu_ij for a list of i temperatures and j carrier densities
    subroutine find_chemical_potential_gausslegendre(wdes_dense_ibz,w_dense_ibz,latt_cur,kpoints, &
                                       nelect,nedos,temperatures,nb_electrons_input, &
                                       nij,muij,comm)

        use constant
        use elphon_base, only: celtot_max
        use elphon_kgrid
        type (wavedes)  :: wdes_dense_ibz
        type (wavespin) :: w_dense_ibz
        type (latt) :: latt_cur
        type (kpoints_struct) :: kpoints
        real(q), intent(in) :: nelect !< number of electrons without doping
        integer, intent(in) :: nedos !< number of energy points to use in the evaluation of the integrated density of states
        real(q), intent(in) :: temperatures(:)  !< list of temperatures at which to compute the chemical potential
        real(q), intent(in) :: nb_electrons_input(:) !< extra number of electrons per unit cell
        real(q), allocatable, intent(out) :: nij(:,:)    !< carrier concentration per cell
        real(q), allocatable, intent(out) :: muij(:,:)   !< chemical potential
        type(communic) :: comm !< communicator to paralelize the work
        ! local variables
        real(q), allocatable :: energy(:), gauss_weights(:)
        integer :: itemp,in,imu
        integer :: iv
        integer :: nbands,nkpoints,nspin,ntemps,nelec
        real(q) :: n0,e0
        real(q) :: nh, ne
        real(q),parameter :: number_electrons_tolerance=1e-10_q ! to be consistent with converged_num_electron value in tet.f
        real(q) :: emin,emax
        real(q) :: emin_save,emax_save

        ! size array
        ntemps=size(temperatures)
        nelec=size(nb_electrons_input)
        nbands=size(w_dense_ibz%celtot,1)
        nkpoints=size(w_dense_ibz%celtot,2)
        nspin=size(w_dense_ibz%celtot,3)

        ! carrier concentration matrix
        allocate(nij(ntemps,nelec))
        do itemp=1,ntemps
           do imu=1,nelec
              nij(itemp,imu)=nb_electrons_input(imu)+nelect
            enddo
        enddo

        ! save energy range
        emin_save=minval(real(w_dense_ibz%celtot))
        emax_save=celtot_max(w_dense_ibz%celtot,w_dense_ibz%wdes%nb_totk)

        ! chemical potential matrix
        allocate(muij(ntemps,nelec))

        do itemp=1,ntemps
           ! allocate energy grid
           call find_energy_grid_gausslegendre(temperatures(itemp),nedos,energy,gauss_weights)

           do in=1,nelec
              !write(*,*) 'computing', itemp, '/', ntemps, '-', in, '/', nelec
              ! find chemical potential from bisection method
              emin=emin_save
              emax=emax_save
              !write(*,'(a,2f18.10)') 'start', e0, nij(itemp,in)
              do
                 ! iter = iter + 1
                 e0=0.5_q*(emin+emax)
                 n0 = number_of_electrons_gausslegendre(kpoints,wdes_dense_ibz,w_dense_ibz%celtot,e0,energy,gauss_weights,comm)
                 !write(*,'(a,4f18.10)') 'result', e0,emin,emax,n0
                 !write(*,'(2i3,4f12.6)') iter, comm%node_me, emax, emin, n0, e0
                 if ( n0 < nij(itemp,in)) emin=e0
                 if ( n0 > nij(itemp,in)) emax=e0
                 if ( (emax-emin) < 1e-6 ) exit ! 0.1 mev precision
                 if (abs(n0-nij(itemp,in)) < number_electrons_tolerance) exit
              enddo
              muij(itemp,in)=e0
           enddo ! loop over carrier densities

           !free energy grid
           deallocate(energy)
           deallocate(gauss_weights)
        enddo ! loop over temperatures

    end subroutine find_chemical_potential_gausslegendre

    !> Compute the number of electrons for a given temperature and chemical potential using Gauss-Legendre integration
    !> The temperature and  chemical potential are encoded in the energy grid (energy array)
    function number_of_electrons_gausslegendre(kpoints,wdes_dense_ibz,celtot,mu,energy,gauss_weights,comm,nbstart,nbstop) result(nij)
        type(kpoints_struct), intent(in) :: kpoints
        type(wavedes)  :: wdes_dense_ibz
        complex(q), intent(in) :: celtot(:,:,:)
        real(q), intent(in) :: mu !< chemical potential
        real(q), intent(in) :: energy(:) !< gauss-legendre energy grid at which to evaluate the integrated density of states
        real(q), intent(in) :: gauss_weights(:) !< weights of the gauss-legendre integration
        real(q) :: nij
        type(communic) :: comm
        integer,optional :: nbstart
        integer,optional :: nbstop
        ! local variables
        integer :: isp, i
        real(q), allocatable :: energy_mu(:)
        real(q), allocatable :: dos(:,:)

        ! shift energy grid
        allocate(energy_mu(size(energy)))
        energy_mu = energy+mu

        ! compute dos
        allocate(dos(size(energy),wdes_dense_ibz%ispin))
        call compute_intdos_tetrahedron(kpoints,wdes_dense_ibz,celtot,energy_mu,dos,comm,nbstart,nbstop)

        ! compute number of electrons
        nij = 0
        do isp=1,wdes_dense_ibz%ispin
            nij = nij + 0.5_q*sum(dos(:,isp)*gauss_weights(:))
        enddo

        ! deallocate
        deallocate(energy_mu)
        deallocate(dos)
    end function number_of_electrons_gausslegendre

    !> Compute the number of electrons for a given temperature and chemical potential using Gauss-Legendre integration
    !> The temperature and  chemical potential are encoded in the energy grid (energy array)
    function number_of_electrons_lineargrid(kpoints,wdes_dense_ibz,w_dense_ibz,mu,energy,temp,comm,nbstart,nbstop) result(nij)
        use elphon_base, only: dfermi_dirac
        type(kpoints_struct), intent(in) :: kpoints
        type(wavedes)  :: wdes_dense_ibz
        type(wavespin) :: w_dense_ibz
        real(q), intent(in) :: mu !< chemical potential
        real(q), intent(in) :: energy(:) !< energy grid at which to evaluate the integrated density of states
        real(q), intent(in) :: temp !< temperature to compute the number of electrons
        real(q) :: nij
        type(communic) :: comm
        integer,optional :: nbstart
        integer,optional :: nbstop
        ! local variables
        integer :: isp, i
        real(q) :: h
        real(q), allocatable :: dos(:,:)

        ! special case of T=0
        if (temp==0) then
            allocate(dos(1,wdes_dense_ibz%ispin))
            call compute_intdos_tetrahedron(kpoints,wdes_dense_ibz,w_dense_ibz%celtot,[mu],dos,comm,nbstart,nbstop)

            nij = 0.0_q
            do isp=1,wdes_dense_ibz%ispin
                nij = nij + dos(1,isp)
            enddo

            ! deallocate
            deallocate(dos)
            return
        endif

        ! compute dos
        allocate(dos(size(energy),wdes_dense_ibz%ispin))
        call compute_intdos_tetrahedron(kpoints,wdes_dense_ibz,w_dense_ibz%celtot,energy,dos,comm,nbstart,nbstop)

        ! compute number of electrons
        nij = 0.0_q
        h = (energy(size(energy))-energy(1))/(size(energy)-1)
        do isp=1,wdes_dense_ibz%ispin
            do i=1,size(energy)-2,2
                nij = nij + dos(i+0,isp)*dfermi_dirac(energy(i+0)-mu,temp*BOLKEV)*h/3
                nij = nij + dos(i+1,isp)*dfermi_dirac(energy(i+1)-mu,temp*BOLKEV)*h*4/3
                nij = nij + dos(i+2,isp)*dfermi_dirac(energy(i+2)-mu,temp*BOLKEV)*h/3
            enddo
        enddo

        ! deallocate
        deallocate(dos)
    end function number_of_electrons_lineargrid

    !> Compute the integral of the density of states using the tetrahedrom method
    subroutine compute_intdos_tetrahedron(kpoints,wdes_dense_ibz,celtot,energy,intdos,comm,nbstart,nbstop)
        use tutor, only: vtutor
        use mpimy, only: communic
        use mkpoints, only: get_indexes, get_tetdisp_24tet, kpoint_integers
        use tet_macdonald, only: weight_24tetra_delta, tetralite_theta, sort_etetra
        type(kpoints_struct), intent(in) :: kpoints
        type(wavedes), intent(in) :: wdes_dense_ibz
        complex(q), intent(in) :: celtot(:,:,:)
        real(q), intent(in) :: energy(:)
        real(q), intent(inout) :: intdos(:,:)
        type(communic) :: comm
        integer,optional :: nbstart
        integer,optional :: nbstop
        ! local variables
        integer :: ik, isp, itet, ib
        integer :: ic, ikibz, iw
        integer :: my_nbstart, my_nbstop
        integer :: idx_k(4,24)
        integer :: kcut24(3,4,24) !< displacements of 24 tetrahedra that contribute to a k-point
        integer :: kint(3), origin(24)
        integer :: kindexes(4), kindexes_sorted(4), perm(4)
        real(q) :: emin, emax
        real(q) :: eleig_tetra(4,24), eig(4), alpha(4)

        PROFILING_START('intdos_tetrahedron')
        emin = minval(energy)
        emax = maxval(energy)
        intdos=0.0_q
        alpha=1.0_q
        my_nbstart = 1
        if (present(nbstart)) my_nbstart = nbstart
        my_nbstop = wdes_dense_ibz%nb_tot
        if (present(nbstop)) my_nbstop = nbstop
        if (kpoints%ntet==0) then
            call vtutor%bug('kpoint mapping was not initialized',__FILE__,__LINE__)
        endif
        ! loop over tetrahedrons
        do itet=1,kpoints%ntet
#ifdef MPI
           if (mod(itet-1,comm%ncpu).ne.comm%node_me-1) cycle
#endif
           kindexes = kpoints%idtet(1:,itet)
           do isp=1,wdes_dense_ibz%ispin
              do ib=my_nbstart,my_nbstop
                 eig(1) = real(celtot(ib,kindexes(1),isp),q)
                 eig(2) = real(celtot(ib,kindexes(2),isp),q)
                 eig(3) = real(celtot(ib,kindexes(3),isp),q)
                 eig(4) = real(celtot(ib,kindexes(4),isp),q)

                 !early skip
                 if (all(eig>emax)) cycle
                 if (all(eig<emin)) then
                    intdos(:,isp) = intdos(:,isp) + kpoints%volwgt*wdes_dense_ibz%rspin*kpoints%idtet(0,itet)
                    cycle
                 endif

                 ! sort eigenvalues
                 perm = [1, 2, 3, 4]
                 call sort_etetra(eig,perm)

                 ! compute the weights for all the energies
                 call tetralite_theta(kpoints%volwgt*wdes_dense_ibz%rspin*kpoints%idtet(0,itet),eig,size(energy),energy,alpha,intdos(:,isp))
              enddo
           enddo
        enddo
        CALLMPI( M_sum_d(comm, intdos, size(intdos)) )
        PROFILING_STOP('intdos_tetrahedron')
    end subroutine compute_intdos_tetrahedron

    !> Compute the dos and integral of the density of states using the tetrahedrom method
    subroutine compute_dosintdos_tetrahedron(kpoints,wdes_dense_ibz,celtot,latt_cur,energy,dos,intdos,comm)
        use mpimy, only: communic
        use mkpoints, only: get_indexes, get_tetdisp_24tet, kpoint_integers
        use tet_macdonald, only: weight_24tetra_delta, tetralite_theta, tetralite_delta, sort_etetra
        type(kpoints_struct), intent(in) :: kpoints
        type(latt), intent(in) :: latt_cur
        type(wavedes), intent(in) :: wdes_dense_ibz
        complex(q), intent(in) :: celtot(:,:,:) !< array with the eigenvalues in the IBZ
        real(q), intent(in) :: energy(:) !< energy grid at which to compute the DOS
        real(q), intent(inout) :: dos(:,:) !< density of states
        real(q), intent(inout) :: intdos(:,:) !< integral of the density of states
        type(communic) :: comm
        ! local variables
        integer :: ik, isp, itet, ib
        integer :: ic, ikibz, iw
        integer :: idx_k(4,24)
        integer :: kcut24(3,4,24) !< displacements of 24 tetrahedra that contribute to a k-point
        integer :: kint(3), origin(24)
        integer :: kindexes(4), kindexes_sorted(4), perm(4)
        real(q) :: emin, emax
        real(q) :: eleig_tetra(4,24), eig(4), alpha(4)

        PROFILING_START('dosintdos_tetrahedron')
        emax = maxval(energy)
        dos=0.0_q
        intdos=0.0_q
        alpha=1.0_q
        ! loop over tetrahedrons
        do itet=1,kpoints%ntet
#ifdef MPI
           if (mod(itet-1,comm%ncpu).ne.comm%node_me-1) cycle
#endif
           kindexes = kpoints%idtet(1:,itet)
           do isp=1,wdes_dense_ibz%ispin
              do ib=1,wdes_dense_ibz%nb_tot
                 eig(1) = real(celtot(ib,kindexes(1),isp),q)
                 eig(2) = real(celtot(ib,kindexes(2),isp),q)
                 eig(3) = real(celtot(ib,kindexes(3),isp),q)
                 eig(4) = real(celtot(ib,kindexes(4),isp),q)

                 !early skip
                 if (all(eig>emax)) cycle

                 ! sort eigenvalues
                 perm = [1, 2, 3, 4]
                 call sort_etetra(eig,perm)

                 ! compute the weights for all the energies
                 call tetralite_delta(kpoints%volwgt*wdes_dense_ibz%rspin*kpoints%idtet(0,itet),eig,size(energy),energy,alpha,dos(:,isp))
                 call tetralite_theta(kpoints%volwgt*wdes_dense_ibz%rspin*kpoints%idtet(0,itet),eig,size(energy),energy,alpha,intdos(:,isp))
              enddo
           enddo
        enddo
        CALLMPI( M_sum_d(comm, dos, size(intdos)) )
        CALLMPI( M_sum_d(comm, intdos, size(intdos)) )
        PROFILING_STOP('dosintdos_tetrahedron')
    end subroutine compute_dosintdos_tetrahedron

    !> Compute the transport functiom, onsager coefficients and transport quantities
    !> Uses vasp implementation of the tetrahedron method
    !> Use Gauss-legendre integration to integrate the transport function and get the onsager coefficients
    subroutine transport_compute_gausslegendre(self,itemp,wdes_dense_ibz,w_dense_ibz,kpoints,comm)
        use base, only: symmetry
        use wave_struct_def, only: wavespin, wavedes
        use string, only: str
        use mpimy, only: communic
        type(transport_type) :: self
        integer,intent(in) :: itemp
        type(wavedes), intent(in) :: wdes_dense_ibz
        type(wavespin), intent(in) :: w_dense_ibz
        type(kpoints_struct), intent(in) :: kpoints
        type(communic), intent(in) :: comm
        ! local variables
        integer :: iv
        integer :: ivbm, icbm
        integer :: i,j,alpha,beta,isp
        real(q) :: cbm, vbm, midgap
        real(q),allocatable :: energy(:),gauss_weights(:)
        real(q),allocatable :: vov(:,:,:,:,:)
        real(q),allocatable :: vovcrosseffmass(:,:,:,:,:,:)

        self%driver = 2
        PROFILING_START('transport_compute_gausslegendre')

        ! compute mean free path in pm
        call mean_free_path(self)

        ! symmetrized product of mean free path and velocity
        call sym_tensor_prod(self%mean_free_path,self%velocity,self%a,self%nrotk,self%isymop,vov,comm)

        ! for plotting
        if (self%do_plot) then
            call find_energy_grid_linear(w_dense_ibz, self%nedos_plot, self%emin_plot, self%emax_plot, self%energy_plot)
            call compute_transport_function_tetrahedron(self%energy_plot,self%transport_function_plot(:,:,:,:,itemp),kpoints,w_dense_ibz%celtot,wdes_dense_ibz,&
                                                        vov,self%vol,comm)
        endif

        ! prepare gauss-legendre energy grid for computation of transport coefficients
        call find_energy_grid_gausslegendre(self%temps(itemp),self%nedos,energy,gauss_weights)
        self%energy(:,itemp) = energy(:) + self%mu(itemp)
        self%emin(itemp) = minval(self%energy(:,itemp))
        self%emax(itemp) = maxval(self%energy(:,itemp))

        ! get transport function in C^2/(m s J)=( ohm m )^-1
        call compute_transport_function_tetrahedron(self%energy(:,itemp),self%transport_function(:,:,:,:,itemp),kpoints,w_dense_ibz%celtot,wdes_dense_ibz,&
                                                    vov,self%vol,comm)
        ! for holes
        !call compute_transport_function_tetrahedron(self%energy,self%transport_function_h,kpoints,w_dense_ibz,wdes_dense_ibz,&
        !                                            vov,self%vol,.false.,comm,nbstart=1,nbstop=self%index_vbm)
        ! for electrons
        !call compute_transport_function_tetrahedron(self%energy,self%transport_function_e,kpoints,w_dense_ibz,wdes_dense_ibz,&
        !                                            vov,self%vol,.false.,comm,nbstart=self%index_cbm,nbstop=wdes_dense_ibz%nb_tot)
        deallocate(vov)

        ! here instead of recomputing by limiting the bands we split the transport function into contributions from the conduction and valence
        vbm = self%bandgap(1)%fundamental_valence%eigenvalue
        cbm = self%bandgap(1)%fundamental_conduction%eigenvalue
        midgap = vbm+(cbm-vbm)/2
        call transport_function_cb_vb(self%energy(:,itemp),self%transport_function(:,:,:,:,itemp),midgap,self%transport_function_e(:,:,:,:,itemp),self%transport_function_h(:,:,:,:,itemp))

        ! compute number of electrons and holes
        ivbm = self%bandgap(1)%fundamental_valence%band
        icbm = self%bandgap(1)%fundamental_conduction%band
        self%nh(itemp) = number_of_electrons_gausslegendre(kpoints,wdes_dense_ibz,w_dense_ibz%celtot,self%mu(itemp),energy,gauss_weights,comm,nbstart=1,nbstop=ivbm)
        self%nh(itemp) = self%nelect0-self%nh(itemp)
        self%ne(itemp) = number_of_electrons_gausslegendre(kpoints,wdes_dense_ibz,w_dense_ibz%celtot,self%mu(itemp),energy,gauss_weights,comm,nbstart=icbm,nbstop=wdes_dense_ibz%nb_tot)

        ! get onsager coefficients
        call get_onsager_coef_gausslegendre(self%transport_function(:,:,:,:,itemp),self%lab(:,:,:,:,itemp),energy,gauss_weights)
        call get_onsager_coef_gausslegendre(self%transport_function_h(:,:,:,:,itemp),self%lab_h(:,:,:,:,itemp),energy,gauss_weights)
        call get_onsager_coef_gausslegendre(self%transport_function_e(:,:,:,:,itemp),self%lab_e(:,:,:,:,itemp),energy,gauss_weights)
        deallocate(gauss_weights)
        deallocate(energy)

        ! get transport coefficients in S.I to be compared to experiments
        call get_transport_coef(self,self%temps(itemp),self%lab(:,:,:,:,itemp),self%e_conductivity(:,:,itemp),self%seebeck(:,:,itemp),&
                                self%peltier(:,:,itemp),self%e_t_conductivity(:,:,itemp))
        call get_mobility(self,self%lab_e(:,:,:,:,itemp),self%ne(itemp),self%mobility_e(:,:,itemp))
        call get_mobility(self,self%lab_h(:,:,:,:,itemp),self%nh(itemp),self%mobility_h(:,:,itemp))
        call set_mobility(self,self%nelect(itemp),self%mobility_e(:,:,itemp),self%mobility_h(:,:,itemp),self%mobility(:,:,itemp))

        PROFILING_STOP('transport_compute_gausslegendre')
    end subroutine transport_compute_gausslegendre

    !> Split the transport function into contributions from the valence and conduction bands
    !> this is done by copying data
    subroutine transport_function_cb_vb(energies,transport_function,midgap,transport_function_e,transport_function_h)
        real(q),intent(in) :: energies(:)
        real(q),intent(in) :: transport_function(:,:,:,:)
        real(q),intent(in) :: midgap
        real(q),intent(out) :: transport_function_e(:,:,:,:)
        real(q),intent(out) :: transport_function_h(:,:,:,:)
        !local variables
        integer :: nspin, i
        real(q) :: e
        nspin = size(transport_function,4)
        !allocate(transport_function_h(size(energies), 3, 3, nspin))
        !allocate(transport_function_e(size(energies), 3, 3, nspin))
        transport_function_h = 0.0_q
        transport_function_e = 0.0_q
        do i=1,size(energies)
            e = energies(i)
            if (e<midgap) then
                transport_function_h(i,:,:,:) = transport_function(i,:,:,:)
            else
                transport_function_e(i,:,:,:) = transport_function(i,:,:,:)
            endif
        enddo
    end subroutine transport_function_cb_vb

    !> Write the header of the transport results
    subroutine transport_print_header(iu,naccumulator)
        integer,intent(in) :: iu
        integer,optional,intent(in) :: naccumulator
        if (iu<0) return
        write(iu,*)
        if (present(naccumulator)) write(iu,'(A,I6)') 'Transport for self-energy accumulator N=', naccumulator
        write(iu,'(7A20)') 'T K','mu eV','sigma S/m','mob cm^2/(V.s)','seebeck μV/K','peltier μV','kappa_e W/(m.K)'
    end subroutine transport_print_header

    !> Write information about transport calculation
    subroutine transport_print(self,idx,iu)
        use tutor, only: vtutor
        use elphon_base, only: elph_selfen_approx_str, elph_scattering_approx_crta
        type(transport_type) :: self
        integer :: idx !< index of the transport calculator
        integer :: iu !< index of the file
        integer :: itemp
        !local variables
        real(q) :: de
        if (iu<0) return
        write(iu,'(A,I3)') 'Transport calculator N =',idx
        write(iu,'(A)')    '----------------------------------'
        select case (self%driver)
            case(1)
                write(iu,'(A)') ' transport driver: 1 ! linear grid and Simpson integration'
            case(2)
                write(iu,'(A)') ' transport driver: 2 ! Gauss-Legendre integration'
            case default
                call vtutor%bug('unknown transport driver selected',__FILE__,__LINE__)
        end select
        write(iu,'(A,A)') ' Scattering approximation: ', elph_selfen_approx_str(self%selfen_approx,verbose=.true.)
        write(iu,'(A,L)') ' Static self-energy: ', self%selfen_approx%static
        write(iu,'(A,I6)') ' Transport number of points:', self%nedos
        do itemp=1,size(self%temps)
            write(iu,'(A,F8.3,A)') ' temperature: ', self%temps(itemp), ' K'
            write(iu,'(A,F8.3,A,F8.3,A,F8.3,A)') '   Transport energy range:   [',self%emin(itemp),':',&
                                                                      self%emax(itemp),']  wich corresponds to ',&
                                                                      (self%emax(itemp)-self%emin(itemp)), ' eV'
            write(iu,'(A,ES13.4,A)') '   Average relaxation time: ', self%tau_average(itemp)*HBAR, ' s'
            write(iu,'(A,ES13.4)') '   Number of electrons: ', self%ne(itemp)
            write(iu,'(A,ES13.4)') '   Number of holes:     ', self%nh(itemp)
            if (self%driver==1.and.self%nedos>1) then
                de = (self%emax(itemp)-self%emin(itemp))/(self%nedos-1)
                write(iu,'(A,F8.3,A,F8.3,A)')            '   Transport energy sampling:', de*1000, ' meV or', de/BOLKEV,' K'
            endif
        enddo
        if (self%do_plot) then
            write(iu,'(A,F8.3,A,F8.3,A,F8.3,A)') ' Transport energy range for plotting:   [',self%emin_plot,':',&
                                                                      self%emax_plot,']  wich corresponds to ',&
                                                                      (self%emax_plot-self%emin_plot), ' eV'
            if (self%driver==1) then
            de = (self%emax_plot-self%emin_plot)/(self%nedos_plot-1)
            write(iu,'(A,F8.3,A,F8.3,A)')            ' Transport energy sampling for plottting:', de*1000, ' meV or', de/BOLKEV,' K'
            endif
        endif
        write(iu,'(A)')
    end subroutine transport_print

    !> Write the results of a transport computation
    subroutine transport_print_results(self,idir1,idir2,itemp,iu)
        type(transport_type) :: self
        integer,intent(in) :: idir1
        integer,intent(in) :: idir2
        integer,intent(in) :: itemp
        integer,intent(in) :: iu
        if (iu<0) return
        write(iu,'(7F20.8)') self%temps(itemp), &
                             self%mu(itemp), &
                             self%e_conductivity(idir1,idir2,itemp), &
                             self%mobility(idir1,idir2,itemp), &
                             self%seebeck(idir1,idir2,itemp)*1e6, &
                             self%peltier(idir1,idir2,itemp)*1e6, &
                             self%e_t_conductivity(idir1,idir2,itemp)
    end subroutine transport_print_results

    subroutine transport_set_results_array_line(self,idir1,idir2,itemp,transport_results)
        type(transport_type) :: self
        integer,intent(in) :: idir1
        integer,intent(in) :: idir2
        integer,intent(in) :: itemp
        real(q),intent(inout) :: transport_results(:)
        transport_results(1) = self%temps(itemp)
        transport_results(2) = self%mu(itemp)
        transport_results(3) = self%e_conductivity(idir1,idir2,itemp)
        transport_results(4) = self%mobility(idir1,idir2,itemp)
        transport_results(5) = self%seebeck(idir1,idir2,itemp)*1e6
        transport_results(6) = self%peltier(idir1,idir2,itemp)*1e6
        transport_results(7) = self%e_t_conductivity(idir1,idir2,itemp)
    end subroutine transport_set_results_array_line

    subroutine print_transport_results_outcar(transport_results,tag,iu)
        real(q),intent(in) :: transport_results(:,:)
        character(len=*),intent(in) :: tag
        integer,intent(in) :: iu
        !local variables
        integer :: iline
        if (iu<0) return
        do iline=1,size(transport_results,2)
            write(iu,'(7F20.8,A,A)') transport_results(:,iline), ' ',tag
        enddo
    end subroutine print_transport_results_outcar

    !> @brief Deallocate all members of the transport structure allocated in transport_init
    subroutine transport_free(self)
        type(transport_type) :: self
        deallocate(self%mean_free_path)
        deallocate(self%energy)
        deallocate(self%emin)
        deallocate(self%emax)

        deallocate(self%transport_function)
        deallocate(self%transport_function_e)
        deallocate(self%transport_function_h)
        deallocate(self%transport_function_hall)

        deallocate(self%lab)
        deallocate(self%lab_e)
        deallocate(self%lab_h)
        deallocate(self%lab_hall)

        deallocate(self%e_conductivity)
        deallocate(self%e_conductivity_e)
        deallocate(self%e_conductivity_h)
        deallocate(self%mobility)
        deallocate(self%mobility_e)
        deallocate(self%mobility_h)
        deallocate(self%e_t_conductivity)
        deallocate(self%seebeck)
        deallocate(self%peltier)

        deallocate(self%ne)
        deallocate(self%nh)

        deallocate(self%tau_average)

        if (self%do_plot) then
            deallocate(self%energy_plot)
            deallocate(self%transport_function_plot)
        endif
    end subroutine transport_free

#ifdef VASP_HDF5
    !> @brief write transport to hdf5 file
    subroutine transport_hdf5write(self,fileid,idx,group,subgroup)
        use vhdf5
        use string, only: str
        use elphon_base, only: elph_selfen_approx_hdf5write
        type(transport_type) :: self
        integer,intent(in) :: idx
        integer(HID_T),intent(in) :: fileid
        character(len=*),optional :: group
        character(len=*),optional :: subgroup
        ! local variables
        integer(HID_T) :: groupid, subgroupid
        character(len=MAX_LEN_GROUP) :: my_group
        character(len=MAX_LEN_GROUP) :: my_subgroup

        my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
        my_subgroup = GRP_ELPHON;   if (present(subgroup)) my_subgroup = trim(subgroup)

        ! open group
        VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group)//'/'//trim(my_subgroup)//'/electrons', groupid))
        VH5_CHECK(vh5_group_open_or_create(groupid, 'transport_'//str(idx), subgroupid))

        ! write metadata
        VH5_CHECK(vh5_write(subgroupid,'mu',self%mu))
        VH5_CHECK(vh5_write(subgroupid,'n',self%nelect))
        VH5_CHECK(vh5_write(subgroupid,'n0',self%nelect0))
        VH5_CHECK(vh5_write(subgroupid,'ne',self%ne))
        VH5_CHECK(vh5_write(subgroupid,'nh',self%nh))
        VH5_CHECK(vh5_write(subgroupid,'cell_volume',self%vol))

        VH5_CHECK(vh5_write(subgroupid,'id_idx',self%id_idx))
        VH5_CHECK(vh5_write(subgroupid,'id_name',self%id_name))

        call elph_selfen_approx_hdf5write(self%selfen_approx,subgroupid)
        VH5_CHECK(vh5_write(subgroupid,'tau_average',self%tau_average*HBAR))

        VH5_CHECK(vh5_write(subgroupid,'temps',self%temps))
        VH5_CHECK(vh5_write(subgroupid,'nedos',self%nedos))
        VH5_CHECK(vh5_write(subgroupid,'dfermi_tol',self%dfermi_tol))

        ! transport function data
        VH5_CHECK(vh5_write(subgroupid,'energy',self%energy))
        VH5_CHECK(vh5_write(subgroupid,'emin',self%emin))
        VH5_CHECK(vh5_write(subgroupid,'emax',self%emax))
        VH5_CHECK(vh5_write(subgroupid,'transport_function',self%transport_function))

        if (self%do_plot) then
            VH5_CHECK(vh5_write(subgroupid,'energy_plot',self%energy_plot))
            VH5_CHECK(vh5_write(subgroupid,'transport_function_plot',self%transport_function_plot))
        endif

        ! final quantities
        VH5_CHECK(vh5_write(subgroupid,'lab',self%lab))
        VH5_CHECK(vh5_write(subgroupid,'e_conductivity',self%e_conductivity))
        VH5_CHECK(vh5_write(subgroupid,'mobility',self%mobility))
        VH5_CHECK(vh5_write(subgroupid,'e_t_conductivity',self%e_t_conductivity))
        VH5_CHECK(vh5_write(subgroupid,'seebeck',self%seebeck))
        VH5_CHECK(vh5_write(subgroupid,'peltier',self%peltier))

        VH5_CHECK(vh5_group_close_writing(subgroupid))
        VH5_CHECK(vh5_group_close_writing(groupid))
    end subroutine transport_hdf5write

    !> @brief write transport to hdf5 file
    subroutine transports_hdf5write(fileid,id_name,id_size,group)
        use vhdf5
        use string, only: str
        use elphon_base, only: elph_selfen_approx_hdf5write
        integer(HID_T),intent(in) :: fileid
        character(len=15) :: id_name(:)
        integer,intent(in) :: id_size(:)
        character(len=*),optional :: group
        ! local variables
        integer(HID_T) :: groupid
        character(len=MAX_LEN_GROUP) :: my_group

        my_group    = GRP_RESULTS//'/'//GRP_ELPHON//'/electrons/transport_meta'
        if (present(group))    my_group    = trim(group)

        ! open group
        VH5_CHECK(vh5_group_open_or_create(fileid, my_group, groupid))

        ! write metadata
        VH5_CHECK(vh5_write(groupid,'ncalculators',product(id_size)))
        VH5_CHECK(vh5_write(groupid,'id_size',id_size))
        VH5_CHECK(vh5_write(groupid,'id_name',id_name))

        VH5_CHECK(vh5_group_close_writing(groupid))
    end subroutine transports_hdf5write
#endif

    !> Compute energy grid for Gauss-Legendre integration
    subroutine find_energy_grid_gausslegendre(temperature,n_energy_points,energy,weights)
        use constant, only: bolkev
        use gauss_quad, only: gauss_legendre
        real(q), intent(in) :: temperature
        integer, intent(in) :: n_energy_points !< number of energy points to use for gauss legendre quadrature
        real(q), allocatable, intent(out) :: energy(:) !< energies for the gauss-legendre integration
        real(q), allocatable, intent(out) :: weights(:) !< weights for the gauss-legendre integration
        !local variables
        integer :: n
        real(q) :: smearing
        real(q), allocatable :: x(:)

        PROFILING_START('energy_grid_gausslegendre')
        n=n_energy_points            ! (exact if n_energy_points is odd but add 1 point if even)
        if (mod(n,2)==0) n=n+1       ! make sure N is odd

        ! compute temperature kbT
        smearing = bolkev*temperature

        ! generate energy grid
        if (smearing > 0.0_q) then
            allocate(x(n))
            allocate(energy(n))
            allocate(weights(n))
            call gauss_legendre(-1.0_q,1.0_q,x,weights,n,eps=1e-12_q)
            energy = smearing * log((1+x)/(1-x))
            deallocate(x)
        else
            allocate(energy(n))
            allocate(weights(n))
            energy(:) = 0.0_q
            weights(:) = 0.0_q
            weights(1) = 2.0_q
        endif
        PROFILING_STOP('energy_grid_gausslegendre')
    end subroutine find_energy_grid_gausslegendre

    !> Get a linear grid of energies between emin and emax
    subroutine find_energy_grid_linear(w, nedos, emin, emax, energy)
        use wave_struct_def, only: wavespin
        use elphon_base, only: linspace, celtot_max
        type(wavespin),intent(in) :: w
        integer,intent(in) :: nedos
        real(q),intent(inout) :: emin
        real(q),intent(inout) :: emax
        real(q),allocatable,intent(out) :: energy(:)
        if (emin==huge(emin)) then
            emin = minval(real(w%celtot,q))
        endif
        if (emax==huge(emax)) then
            emax = celtot_max(w%celtot,w%wdes%nb_totk)
        endif
        energy = linspace(emin,emax,nedos)
    end subroutine find_energy_grid_linear

    !> @brief Compute mean free path
    !> Currently we are using the linearized Boltzman equation so the mean free path amounts to v*tau
    !> when implementing an iterative solver of the Bolzman equation we expand this routine.
    subroutine mean_free_path(self)
        type(transport_type) :: self
        ! local variables
        integer :: i

        do i=1,3
           self%mean_free_path(:,:,:,i)=self%velocity(:,:,:,i)*self%tau(:,:,:)
        enddo
    end subroutine mean_free_path

    subroutine mean_free_path_mag(self)
        type(transport_type) :: self
        !local variables
        integer :: i

        allocate(self%mean_free_path_mag(self%nbands,self%nkpoints,self%nspin,3,3))

        do i=1,3
           self%mean_free_path_mag(:,:,:,i,1)=(self%velocity(:,:,:,2)*self%effective_mass(:,:,:,i,3)&
                                              -self%velocity(:,:,:,3)*self%effective_mass(:,:,:,i,2))*self%tau(:,:,:)**2
           self%mean_free_path_mag(:,:,:,i,2)=(self%velocity(:,:,:,3)*self%effective_mass(:,:,:,i,1)&
                                              -self%velocity(:,:,:,1)*self%effective_mass(:,:,:,i,3))*self%tau(:,:,:)**2
           self%mean_free_path_mag(:,:,:,i,3)=(self%velocity(:,:,:,1)*self%effective_mass(:,:,:,i,2)&
                                              -self%velocity(:,:,:,2)*self%effective_mass(:,:,:,i,1))*self%tau(:,:,:)**2
        enddo
    end subroutine mean_free_path_mag

    !> Return the (symetrized (over SPG operations)) tensor product of 2 quantities which transform as group velocity
    subroutine sym_tensor_prod(velocity_1,velocity_2,a,nrotk,isymop,vov,comm)
        use mpimy, only: communic
        real(q), intent(in) :: velocity_1(:,:,:,:) !< (nbands,nkpoints,nspin,dir)
        real(q), intent(in) :: velocity_2(:,:,:,:) !< (nbands,nkpoints,nspin,dir)
        real(q), intent(in) :: a(3,3) !< lattice vectors
        integer, intent(in) :: nrotk !< number of symmetry operations
        integer, intent(in) :: isymop(3,3,48) !< integer representation of the symmetry operations
        real(q), allocatable, intent(out) :: vov(:,:,:,:,:)     !< (nbands,nkpoints,nspin,dir1,dir2)
        type(communic), intent(in) :: comm
        !local variables
        integer :: nbands,nkpoints,nspin,dir1,dir2
        integer :: ik,ib,isp

        PROFILING_START('sym_tensor_prod')
        nbands=size(velocity_1,1)
        nkpoints=size(velocity_1,2)
        nspin=size(velocity_1,3)

        if (nbands .ne. size(velocity_2,1)) then
           call vtutor%bug('Error in sym_tensor_prod(1)',__FILE__,__LINE__)
        endif

        if (nkpoints .ne. size(velocity_2,2)) then
           call vtutor%bug('Error in sym_tensor_prod(2)',__FILE__,__LINE__)
        endif

        if (nspin .ne. size(velocity_2,3)) then
           call vtutor%bug('Error in sym_tensor_prod(3)',__FILE__,__LINE__)
        endif

        if ((size(velocity_1,4) .ne. 3) .or. (size(velocity_2,4) .ne. 3) ) then
            call vtutor%bug('Error in sym_tensor_prod(4)',__FILE__,__LINE__)
        endif

        allocate(vov(3,3,nbands,nkpoints,nspin))
        vov=0.0_q

        do isp=1,nspin
           do ik=1,nkpoints
#ifdef MPI
              if (mod(ik-1,comm%ncpu).ne.comm%node_me-1) cycle
#endif
              do ib=1,nbands

                 do dir2=1,3
                    do dir1=1,3
                       vov(dir1,dir2,ib,ik,isp)=velocity_1(ib,ik,isp,dir1)*velocity_2(ib,ik,isp,dir2)
                    enddo
                 enddo
                 CALL TSYM(vov(:,:,ib,ik,isp),isymop,nrotk,a) ! return 1/#R sum Rv otime Rv
              enddo
           enddo
        enddo
        CALLMPI( M_sum_d(comm, vov, size(vov)) )
        PROFILING_STOP('sym_tensor_prod')
    end subroutine sym_tensor_prod

    !> Return the symmetrized third order tensor product of a group velocity and a second order tensor
    subroutine sym_tensor3_prod(tensor2,velocity,a,nrotk,isymop,vovcrosseff_mass,comm)
        use mpimy, only: communic
        real(q), intent(in) :: tensor2(:,:,:,:,:) !< (nbands,nkpoints,nspin,dir,dir)
        real(q), intent(in) :: velocity(:,:,:,:) !< (nbands,nkpoints,nspin,dir)
        real(q), intent(in) :: a(3,3) !< lattice vectors
        integer, intent(in) :: nrotk !< number of symmetry operations
        integer, intent(in) :: isymop(3,3,48) !< integer representation of the symmetry operations
        real(q), allocatable, intent(out) :: vovcrosseff_mass(:,:,:,:,:,:)     !< (nbands,nkpoints,nspin,dir1,dir2,dir3)
        type(communic), intent(in) :: comm
        !local variables
        integer :: nbands,nkpoints,nspin,dir1,dir2,dir3
        integer :: ik,ib,isp

        PROFILING_START('sym_tensor3_prod')
        nbands=size(velocity,1)
        nkpoints=size(velocity,2)
        nspin=size(velocity,3)

        if (nbands /= size(tensor2,1)) then
           call vtutor%bug('Wrong number of bands in tensor2',__FILE__,__LINE__)
        endif

        if (nkpoints /= size(tensor2,2)) then
           call vtutor%bug('Wrong number of kpoints in tensor2',__FILE__,__LINE__)
        endif

        if (nspin /= size(tensor2,3)) then
           call vtutor%bug('Wrong number of spin components in tensor2',__FILE__,__LINE__)
        endif

        if ((size(tensor2,4) /= 3) .or. (size(tensor2,5) /= 3) ) then
            call vtutor%bug('Wrong dimensions of tensor2',__FILE__,__LINE__)
        endif

        allocate(vovcrosseff_mass(3,3,3,nbands,nkpoints,nspin))
        vovcrosseff_mass=0.0_q

        do isp=1,nspin
           do ik=1,nkpoints
#ifdef MPI
              if (mod(ik-1,comm%ncpu).ne.comm%node_me-1) cycle
#endif
              do ib=1,nbands

                 do dir3=1,3
                    do dir2=1,3
                       do dir1=1,3
                          vovcrosseff_mass(dir1,dir2,dir3,ib,ik,isp)=velocity(ib,ik,isp,dir1)*tensor2(ib,ik,isp,dir2,dir3)
                       enddo
                    enddo
                 enddo
                 CALL TSYM3(vovcrosseff_mass(:,:,:,ib,ik,isp),isymop,nrotk,a) ! return 1/#R sum Rv otime Rv
              enddo
           enddo
        enddo
        CALLMPI( M_sum_d(comm, vovcrosseff_mass, size(vovcrosseff_mass)) )
        PROFILING_STOP('sym_tensor3_prod')
    end subroutine sym_tensor3_prod

    !> Get derivatives of the fermi occupation functions with respect to a change in the electric field
    !> in the linear response regime
    !> \f(
    !> \partial_{E_\beta} f_{n\mathbf{k}} =
    !> e \frac{\partial f_{n\mathbf{k}}}{\partial \epsilon_{n\mathbf{k}}}
    !> \tau_{n\mathbf{k}}
    !> \frac{\partial \epsilon_{n\mathbf{k}}}{\partial \mathbf{k}_\beta}
    !> \f)
    subroutine get_linear_response_dfnk_de
    end subroutine get_linear_response_dfnk_de

    !> Get derivatives of the occupation functions with respect to a change in the electric and magnetic fields
    !> in the linear response regime
    !> \f(
    !>\partial_{E_\beta} \partial_{B_\gamma} f_{n\mathbf{k}} =
    !>-\frac{e}{\hbar} \tau_{n\mathbf{k}}
    !>\left(
    !>    \frac{\partial \epsilon_{n\mathbf{k}}}{\partial \mathbf{k}} \times\frac{\partial}{\partial \mathbf{k}}
    !> \right)_\gamma
    !> \partial_{E_\beta} f_{n\mathbf{k}}
    !> \f)
    subroutine get_linear_response_d2fnk_dedb
    end subroutine get_linear_response_d2fnk_dedb

    !> take a self-energy, in eV, computed in betwwen band_start,band_stop among nbands,
    !> and return the lifetime/relaxation tau in 1/eV
    !> I would like to read a file outside and add here
    subroutine relaxation_time(band_start,band_stop,nbands,nkpoints,nspin,self_energy,tau)
        use constant
        integer, intent(in) :: band_start
        integer, intent(in) :: band_stop
        integer, intent(in) :: nbands
        integer, intent(in) :: nkpoints
        integer, intent(in) :: nspin
        complex(q), allocatable, intent(in) :: self_energy(:,:,:)  !< (nbands,nkpoints,nspin)
        real(q), intent(out) :: tau(nbands,nkpoints,nspin)
        ! local variables
        real(q), parameter :: TINY = 1e-16    ! roughly correspond to a second for tau, probably safe upper bound
        integer ib,ik,isp
        real(q) :: im_s

        tau=0.0_q

        do isp=1,nspin
           do ik=1,nkpoints                ! irreducible points
              do ib=band_start,band_stop
                 im_s=abs(aimag(self_energy(ib,ik,isp)))
                 if (im_s > TINY) then
                    tau(ib,ik,isp)=1.0_q/(2.0_q*im_s)
                 endif
              enddo
           enddo
        enddo
    end subroutine relaxation_time

    !> Compute the Onsager coefficients by integrating the transport function on a linear grid using the simpson rule
    subroutine get_onsager_coef_uniform(temperature,mu,transport_function,lab,energy)
    real(q),intent(in) :: temperature
    real(q),intent(in) :: mu
    real(q),intent(in) :: transport_function(:,:,:,:)
    real(q),intent(out) :: lab(:,:,:,:)
    real(q),intent(in) :: energy(:)
    !local variables
    integer :: i,j,ie,alpha,beta,isp
    integer :: nspin,nedos
    real(q) :: df,h

    nedos = size(transport_function,1)
    nspin = size(transport_function,4)
    lab=0.0_q

    h=(energy(nedos)-energy(1))/(nedos-1)
    do isp=1,nspin
       do i=1,2
          do j=1,2
             do alpha=1,3
                do beta=1,3
                   ! for the T=0 case we use a single point in the integration grid because the derivative of the fermi dirac
                   ! distribution becomes a delta function
                   if (nedos == 1) then
                       lab(i,j,alpha,beta) = transport_function(1,alpha,beta,isp)
                   else
                   do ie=1,nedos-2,2

                       if (dfermi(energy(ie)-mu,bolkev*temperature,energy_0,df)) then
                          lab(i,j,alpha,beta)=lab(i,j,alpha,beta)+h/3*transport_function(ie,alpha,beta,isp)* & 
                                                    df*(mu-energy(ie))**(i+j-2)
                       endif
                       if (dfermi(energy(ie+1)-mu,bolkev*temperature,energy_0,df)) then
                          lab(i,j,alpha,beta)=lab(i,j,alpha,beta)+h*4/3*transport_function(ie+1,alpha,beta,isp)* &
                                                    df*(mu-energy(ie+1))**(i+j-2)
                       endif
                       if (dfermi(energy(ie+2)-mu,bolkev*temperature,energy_0,df)) then
                          lab(i,j,alpha,beta)=lab(i,j,alpha,beta)+h/3*transport_function(ie+2,alpha,beta,isp)* &
                                                    df*(mu-energy(ie+2))**(i+j-2)
                       endif

                   enddo
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo

    end subroutine get_onsager_coef_uniform

    !> Compute onsager coefficients using gauss legendre integration
    subroutine get_onsager_coef_gausslegendre(transport_function,lab,energies,gauss_weights)
        real(q),intent(in) :: transport_function(:,:,:,:)
        real(q),intent(out) :: lab(:,:,:,:)
        real(q),intent(in) :: energies(:)
        real(q),intent(in) :: gauss_weights(:)
        ! local variables
        integer :: alpha,beta
        integer :: i,j
        integer :: isp, nspin
        PROFILING_START('get_onsager_coef_gausslegendre')
        nspin = size(transport_function,4)
        ! compute onsager coefficients
        lab = 0.0_q
        do beta=1,3
          do alpha=1,3
            do isp=1,nspin
            do j=1,2
              do i=1,2
                lab(i,j,alpha,beta) = &
                lab(i,j,alpha,beta) +0.5_q*sum(gauss_weights(:)*transport_function(:,alpha,beta,isp)*(-energies(:))**(i+j-2))
              enddo
            enddo
            enddo
          enddo
        enddo
        PROFILING_STOP('get_onsager_coef_gausslegendre')
    end subroutine get_onsager_coef_gausslegendre

    !> Compute the transport coefficients from the Onsager coefficients
    subroutine get_transport_coef(self,temperature,lab,e_conductivity,seebeck,peltier,e_t_conductivity)
        type(transport_type) :: self
        real(q),intent(in) :: temperature
        real(q),intent(in) :: lab(:,:,:,:)
        real(q),intent(out) :: e_conductivity(:,:)
        real(q),intent(out) :: seebeck(:,:)
        real(q),intent(out) :: peltier(:,:)
        real(q),intent(out) :: e_t_conductivity(:,:)

        ! electronic conductivity in  C^2/(m s J)=(ohm m)^-1
        e_conductivity(1:3,1:3)=lab(1,1,1:3,1:3)

        ! Seebeck coefficient in J/C/K = V/K
        if (det_3x3(e_conductivity) == 0.0_q) then
           seebeck(:,:) = 0.0_q
        else
           seebeck=matrix_inverse_3x3(e_conductivity)
        endif
        seebeck=matmul(seebeck,lab(1,2,1:3,1:3))
        if (temperature > tiny(temperature) ) then
           seebeck=seebeck/temperature
        else
           seebeck=0.0_q
        endif

        ! Peltier coefficient in V --> change the computation in a way which also to check PI = T * S^t
        peltier=transpose(seebeck)*temperature

        ! electronic thermal conductivity in J /(m s)/K = W/(m K)
        if (det_3x3(e_conductivity) == 0.0_q) then
           e_t_conductivity(:,:) = 0.0_q
        else
           e_t_conductivity=matrix_inverse_3x3(e_conductivity)
        endif
        e_t_conductivity=matmul(matmul(lab(2,1,1:3,1:3),e_t_conductivity),lab(1,2,1:3,1:3))
        e_t_conductivity=lab(2,2,1:3,1:3)-e_t_conductivity
        if (temperature > tiny(temperature)) then
           e_t_conductivity=e_t_conductivity/temperature
        else
           e_t_conductivity=0.0_q
        endif

    end subroutine get_transport_coef

    !> Compute the mobility from the electronic conductivity
    subroutine get_mobility(self,lab,n,mobility)
        type(transport_type) :: self
        real(q),intent(in) :: lab(:,:,:,:)
        real(q),intent(in) :: n !< number of carriers (electrons or holes) per unit cell
        real(q),intent(out) :: mobility(:,:)
        !local variables
        real(q) :: carrier_density,volume

        volume = 1e-30*self%vol ! from Angstroem^3 to m^3
        carrier_density = EVTOJ * n / volume !C/m^3

        ! electronic mobility [cm^2/(V.s)] = electronic conductivity [ohm m]^-1
        ! / carrier_density [C/m^3] * 100^2 [cm^2/m^2]
        if (carrier_density > tiny(carrier_density)) then
            mobility(1:3,1:3) = lab(1,1,1:3,1:3) / carrier_density * 100**2 ! m^2/(V.s) to cm^2/(V.s)
        else
           mobility(1:3,1:3) = 0.0_q
        endif
    end subroutine get_mobility

    ! set mobility depending on whether we are considering electron or hole doping
    subroutine set_mobility(self,nelect,mobility_e,mobility_h,mobility)
        type(transport_type) :: self
        real(q),intent(in) :: nelect
        real(q),intent(in) :: mobility_e(3,3)
        real(q),intent(in) :: mobility_h(3,3)
        real(q),intent(out) :: mobility(3,3)
        if (nelect-self%nelect0>0) then
            mobility = mobility_e
        else
            mobility = mobility_h
        endif
    end subroutine set_mobility

    !> Compute the carrier density and chemical potential for different dopings and temperatures.
    !> This routine does something similar to GET_Nij_MUij except that this one uses DENSTA
    subroutine get_nij_muij_densta(elph_set, kpoints, latt_cur, wdes_dense_ibz, w_dense_ibz, info, selfen_nij, selfen_muij, efermi0, comm)
        use constant, only: bolkev
        use mpimy, only: communic
        use base, only: info_struct, in_struct
        use elphon_base, only: fermi_dirac, celtot_max
        use tet, only: tetrahedron_calculator, new_calculator, fermi_weights_state
        use mkpoints, only: tetirr, kpoints_struct, generatinglattice
        use fermi_energy, only : delstp
        use tutor, only: vtutor
        use string, only: str
        type (elph_settings)  :: elph_set
        type (kpoints_struct) :: kpoints
        type (latt)           :: latt_cur
        type (wavedes)  :: wdes_dense_ibz
        type (wavespin) :: w_dense_ibz
        type (info_struct) :: info
        real(q), allocatable, intent(out) :: selfen_muij(:,:) !< chemical potential matrix
        real(q), allocatable, intent(out) :: selfen_nij(:,:)  !< number of electron per cell matrix
        real(q) :: efermi0 !< fermi energy
        type (communic) :: comm
        ! Local variables
        integer,parameter :: nedos=301 ! not really used because densta only computed the fermi energy
        integer :: i,j,ib,ik
        integer :: ispin, ntemps, ncarriers
        real(q) :: sigma, nelect
        real(q) :: efermi
        real(q) :: eps, dfun, sfun
        real(q) :: entropy, par(1,1,1,1,wdes_dense_ibz%ncdij)
        real(q) :: dospar(1,1,1,wdes_dense_ibz%ncdij)
        real(q) :: dos(nedos,wdes_dense_ibz%ispin),dosi(nedos,wdes_dense_ibz%ispin)
        real(q),parameter :: SIGMA_MIN = 1*bolkev
        real(q),allocatable :: energy(:), gauss_weights(:) !only used when gauss-legendre integration is activated here
        real(q),allocatable :: nelect_delta(:)
        real(q),allocatable :: nij_tmp(:,:), muij_tmp(:,:)
        type(tetrahedron_calculator) calculator
        type(fermi_weights_state) state

        ! consistency check
        if (kpoints%ismear/=elph_set%ismear) then
            call vtutor%bug('The value of ismear is inconsistent in get_nij_muij_densta',__FILE__,__LINE__)
        endif

        ! to be safe set emin and emax of kpoints here
        kpoints%emin = minval(real(w_dense_ibz%celtot,q))
        kpoints%emax = celtot_max(w_dense_ibz%celtot,w_dense_ibz%wdes%nb_totk)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! First we determine the intrinsic doping (at zero temperature)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        select case(elph_set%ismear)
           case( 0:,-1,-4,-5,-14,-15)
             call densta( -1, -1, wdes_dense_ibz, w_dense_ibz, kpoints, &
                  info%nelect, info%nup_down, entropy, efermi0, kpoints%sigma, &
                  .false., nedos, 0, 0, dos, dosi, par, dospar)
           case(-24)
             ! call this routine with the temperature that corresponds to sigma in the INCAR and intrinsic number of electrons
             call find_chemical_potential_gausslegendre(wdes_dense_ibz,w_dense_ibz,latt_cur,kpoints,&
                                          info%nelect,elph_set%fermi_nedos,[kpoints%sigma/BOLKEV],[0.0_q],nij_tmp,muij_tmp,comm)
             efermi0 = muij_tmp(1,1)
             deallocate(nij_tmp)
             deallocate(muij_tmp)
           case default
             call invalid_ismear
        end select

        ! compute number of temperatures
        ntemps = size(elph_set%selfen_temps)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Use chemical potential to determine the number of carriers
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (elph_set%selfen_nmu > 0) then

            allocate(selfen_muij(ntemps,elph_set%selfen_nmu))
            allocate(selfen_nij(ntemps,elph_set%selfen_nmu))

            ! loop over temperatures
            do i=1,ntemps
              sigma = bolkev*elph_set%selfen_temps(i)
              if (elph_set%ismear>=-1) sigma = max(sigma,SIGMA_MIN)
              ! loop over the chemical potentials
              do j=1,elph_set%selfen_nmu
                 efermi = efermi0 + elph_set%selfen_mu(j)
                 ! Here we have the different cases
                 select case(elph_set%ismear)
                   ! Methfessel-Paxon
                   case (0:)
                      do ispin=1,wdes_dense_ibz%ispin
                        do ik=1,wdes_dense_ibz%nkpts
                          do ib=1,wdes_dense_ibz%nbands
                             eps = real(w_dense_ibz%celtot(ib,ik,ispin),q)
                             call delstp(kpoints%ismear,(efermi-eps)/sigma,dfun,sfun)
                             w_dense_ibz%ferwe(ib,ik,ispin) = wdes_dense_ibz%rspin*sfun*wdes_dense_ibz%wtkpt(ik)
                          enddo
                        enddo
                      enddo
                      ! Compute the number of electrons
                      nelect = sum(w_dense_ibz%ferwe)

                   ! Fermi smearing
                   case(-1)
                      do ispin=1,wdes_dense_ibz%ispin
                        do ik=1,wdes_dense_ibz%nkpts
                          do ib=1,wdes_dense_ibz%nbands
                             eps = real(w_dense_ibz%celtot(ib,ik,ispin),q)
                             w_dense_ibz%ferwe(ib,ik,ispin) = wdes_dense_ibz%rspin*fermi_dirac(eps-efermi,sigma)*wdes_dense_ibz%wtkpt(ik)
                          enddo
                        enddo
                      enddo
                      ! Compute the number of electrons
                      nelect = sum(w_dense_ibz%ferwe)

                   ! tetrahedron method
                   case(-4,-14)
                      calculator = new_calculator(real(w_dense_ibz%celtot), wdes_dense_ibz%rspin, kpoints%volwgt * kpoints%idtet(0,:),&
                                                  kpoints%idtet(1:4,:), smearing=sigma, use_deltas=.true.)
                      state%fermi_energy = efermi
                      state%use_bloechl = .false.
                      state%unit_out = -1 !io%iu6
                      call calculator%fermi_weights(state, w_dense_ibz%fertot)
                      ! Compute the number of electrons
                      nelect = sum(w_dense_ibz%fertot)

                   ! tetrahedron method with Bloechl corrections
                   case(-5,-15)
                      calculator = new_calculator(real(w_dense_ibz%celtot), wdes_dense_ibz%rspin, kpoints%volwgt * kpoints%idtet(0,:),&
                                                  kpoints%idtet(1:4,:), smearing=sigma, use_deltas=.true.)
                      state%fermi_energy = efermi
                      state%use_bloechl = .true.
                      state%unit_out = -1 !io%iu6
                      call calculator%fermi_weights(state, w_dense_ibz%fertot)
                      ! Compute the number of electrons
                      nelect = sum(w_dense_ibz%fertot)

                   ! Gauss-Legendre
                   case (-24)
                      ! allocate energy grid
                      call find_energy_grid_gausslegendre(elph_set%selfen_temps(i),elph_set%fermi_nedos,energy,gauss_weights)
                      nelect = number_of_electrons_gausslegendre(kpoints,wdes_dense_ibz,&
                                                      w_dense_ibz%celtot,efermi,energy,gauss_weights,comm)
                      deallocate(energy)
                      deallocate(gauss_weights)
                   case default
                      call invalid_ismear
                 end select

                 ! save the results
                 selfen_nij(i,j)=nelect
                 selfen_muij(i,j)=efermi

              enddo ! loop over the chemical potentials
            enddo ! loop over the temperatures

        endif !(elph_set%selfen_nmu > 0)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! use number of electrons to determine the chemical potential
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (elph_set%selfen_ncarrier_den > 0) then
            allocate(nelect_delta(elph_set%selfen_ncarrier_den))
            nelect_delta = 1.0e-24*elph_set%selfen_carrier_den*latt_cur%omega
        elseif (elph_set%selfen_ncarrier_per_cell > 0) then
            allocate(nelect_delta(elph_set%selfen_ncarrier_per_cell))
            nelect_delta = elph_set%selfen_carrier_per_cell
        endif

        if (elph_set%selfen_ncarrier_den>0 .or. elph_set%selfen_ncarrier_per_cell>0) then
            ncarriers = size(nelect_delta,1)
            allocate(selfen_muij(ntemps,ncarriers))
            allocate(selfen_nij(ntemps,ncarriers))

            select case(elph_set%ismear)
               case(0:,-1,-4,-5,-14,-15)
                  ! loop over temperatures
                  do i=1,ntemps
                     sigma  = bolkev*elph_set%selfen_temps(i)
                     if (elph_set%ismear>=-1) sigma = max(sigma,SIGMA_MIN)
                     ! loop over the number of carriers
                     do j=1,ncarriers
                        nelect = info%nelect + nelect_delta(j)
                        ! compute density of states on the dense grid
                        call densta( -1, -1, wdes_dense_ibz, w_dense_ibz, kpoints, &
                             nelect, info%nup_down, entropy, efermi, sigma, &
                             .true., nedos, 0, 0, dos, dosi, par, dospar)
                        ! store carrier density and fermi level
                        selfen_nij(i,j)=nelect
                        selfen_muij(i,j)=efermi
                     enddo !icarriers
                  enddo !itemps
               case(-24)
                  call find_chemical_potential_gausslegendre(wdes_dense_ibz,w_dense_ibz,latt_cur,kpoints,&
                                               info%nelect,elph_set%fermi_nedos,elph_set%selfen_temps,nelect_delta,selfen_nij,selfen_muij,comm)
               case default
                  call invalid_ismear
            end select

        endif !elph_set%selfen_ncarrier_den>0 .or. elph_set%selfen_ncarrier_per_cell>0

        contains
        subroutine invalid_ismear
          call vtutor%error('Invalid option for ELPH_ISMEAR='//str(elph_set%ismear)//' found in the INCAR.\n '//&
                            'The valid options are: \n'//&
                            '  >1: methfessel-paxton\n'//&
                            '   0: gaussian smearing\n'//&
                            '  -1: fermi-dirac smearing\n'//&
                            '  -4: tetrahedron T=0\n'//&
                            '  -5: tetrahedron T=0 with Bloechl corrections\n'//&
                            ' -14: tetrahedron finite T \n'//&
                            ' -15: tetrahedron finite T with Bloechl corrections\n'//&
                            ' -24: same as –14 but using a faster and memory saving algorithm')
        end subroutine

    end subroutine get_nij_muij_densta

    !> Write header of the chemical potential calculation
    subroutine write_nij_muij_header(elph_set,iu)
        type(elph_settings) :: elph_set
        integer :: iu
        if (iu>0) then
            write(iu,'(a)') 'Chemical potential calculation:'
            write(iu,'(a)') '==============================='
            write(iu,'(a)') ''
            write(iu,'(a,i3)') 'elph_ismear=', elph_set%ismear
            if (elph_set%ismear==-24) &
            write(iu,'(a,i8)') 'elph_fermi_nedos=', elph_set%fermi_nedos
            if (elph_set%selfen_ncarrier_per_cell>0) then
                write(iu,'(a)') 'elph_selfen_carrier_per_cell='
                call write_array(elph_set%selfen_carrier_per_cell,'(f12.3)',iu)
            endif
            if (elph_set%selfen_ncarrier_den>0) then
                write(iu,'(a)') 'elph_selfen_carrier_den='
                call write_array(elph_set%selfen_carrier_den,'(e12.3)',iu)
            endif
            if (elph_set%selfen_nmu>0) then
                write(iu,'(a)') 'elph_selfen_mu='
                call write_array(elph_set%selfen_mu,'(f12.3)',iu)
            endif
            if (size(elph_set%selfen_temps)>0) then
                write(iu,'(a)') 'elph_selfen_temps='
                call write_array(elph_set%selfen_temps,'(f12.3)',iu)
            endif
            write(iu,'(a)') ''
        endif
        contains
        subroutine write_array(a,fmt,iu)
            real(q) :: a(:)
            character(len=*) :: fmt
            integer :: iu
            !local
            integer :: i
            do i=1,size(a)
                write(iu,fmt) a(i)
            enddo
        end subroutine write_array
    end subroutine write_nij_muij_header

    !> Write results of the chemical potential calculation
    subroutine write_nij_muij_final(temperatures, selfen_nij, selfen_muij, sigma, efermi, iu)
        real(q) :: temperatures(:) !< temperatures at which the chemical potential and number of electrons are computed
        real(q) :: selfen_muij(:,:) !< chemical potential matrix (itemp,inelect)
        real(q) :: selfen_nij(:,:)  !< number of electron per cell matrix (itemp,inelect)
        real(q) :: sigma !< temperature at which to compute the fermi enegy
        real(q) :: efermi !< fermi energy on the dense grid
        integer :: iu !< unit of the file to write to
        ! local variables
        integer :: i, j

        if (iu<0) return
        write(iu,'(a,f20.10,a)') ' Fermi energy: ',efermi,' eV (dense k-point grid)'
        write(iu,'(a,f12.6,a,f12.6,a)') 'computed at a temperature of ',sigma/BOLKEV,' K or sigma=',sigma,' eV'
        write(iu,'(a)') ''
        write(iu,'(a)') '                   Number of electrons per cell'
        write(iu,'(a)') '                   ----------------------------'
        do i=1,size(selfen_nij(:,1))
           write(iu,'(a2,1x,f15.8,100f15.8)') 'T=',temperatures(i), (selfen_nij(i,j),j=1,size(selfen_nij(1,:)))
        enddo
        write(iu,'(a)') '                   ----------------------------'
        write(iu,'(a)') '                       Chemical potential'
        write(iu,'(a)') '                   ----------------------------'
        do i=1,size(selfen_muij(:,1))
           write(iu,'(a2,1x,f15.8,100f15.8)') 'T=',temperatures(i),(selfen_muij(i,j),j=1,size(selfen_muij(1,:)))
        enddo
        write(iu,'(a)') '                   ----------------------------'
        write(iu,'(a)') ''
    end subroutine write_nij_muij_final

#ifdef VASP_HDF5
    !> Write information about the electron-phonon calculation to the hdf5 file
    subroutine chemical_potential_hdf5write(temperatures,selfen_carrier_den,selfen_carrier_per_cell,selfen_mu,efermi0,selfen_nij,selfen_muij,groupid)
        use vhdf5
        use string, only: str
        use elphon_base, only: elph_selfen_approx_hdf5write
        real(q),intent(in) :: temperatures(:) !< temperatures at which to compute number of carriers and chemical potential
        real(q),allocatable,intent(in) :: selfen_carrier_den(:) !< list of carrier densities
        real(q),allocatable,intent(in) :: selfen_carrier_per_cell(:) !< list of carriers per unit cell
        real(q),allocatable,intent(in) :: selfen_mu(:) !< list of chemical potentials
        real(q),intent(in) :: efermi0 !< fermi energy at zero temperature
        real(q),intent(in) :: selfen_nij(:,:) !< number of carriers for every (temperature, number of carriers)
        real(q),intent(in) :: selfen_muij(:,:) !< chemical potential for every (temperature, number of carriers)
        integer(HID_T),intent(in) :: groupid

        if (allocated(selfen_carrier_den)) then
            VH5_CHECK(vh5_write(groupid,'selfen_carrier_den',selfen_carrier_den))
        endif
        if (allocated(selfen_carrier_per_cell)) then
            VH5_CHECK(vh5_write(groupid,'selfen_carrier_per_cell',selfen_carrier_per_cell))
        endif
        if (allocated(selfen_mu)) then
            VH5_CHECK(vh5_write(groupid,'selfen_mu',selfen_mu))
        endif
        VH5_CHECK(vh5_write(groupid,'temperatures',temperatures))
        VH5_CHECK(vh5_write(groupid,'efermi',efermi0))
        VH5_CHECK(vh5_write(groupid,'muij',selfen_muij))
        VH5_CHECK(vh5_write(groupid,'nij',selfen_nij))
    end subroutine chemical_potential_hdf5write
#endif

    !> Derivative of the Fermi-Dirac distribution
    function dfermi(e,kbt,de,df) result(succeeded)
        real(q), intent(in) :: e
        real(q), intent(in) :: kbt
        real(q), intent(in) :: de
        real(q), intent(out) :: df
        logical :: succeeded
        real(q) x

        if (kbt < de) then
           succeeded = .false.
        else
           x=e/kbt
           df=1.0d0/(2.0d0+exp(x)+exp(-x))
           df=df/kbt
           if ( df < huge(1.0_q)) then
              succeeded = .true.
           else
              succeeded = .false.
           endif
        endif
    end function dfermi

end module electron_transport

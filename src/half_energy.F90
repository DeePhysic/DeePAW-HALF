module half_energy
  use half_kinds,only:dp
  use half_constants,only:pi,felect
  use half_types,only:crystal_t
  use half_kpoints,only:kpoint_set_t
  implicit none
  private
  public::compute_occupations,ewald_energy,read_vasp_eigenval
contains
  subroutine read_vasp_eigenval(path,set,eigenvalues,nelect)
    character(len=*),intent(in)::path
    type(kpoint_set_t),intent(out)::set
    real(dp),allocatable,intent(out)::eigenvalues(:,:)
    real(dp),intent(out)::nelect
    character(len=2048)::line
    integer::unit,ios,i,ik,ib,nelectron,nk,nb,index
    open(newunit=unit,file=path,status='old',action='read',iostat=ios)
    if(ios/=0)error stop 'HALF: cannot open EIGENVAL'
    do i=1,5;read(unit,'(A)',iostat=ios)line;if(ios/=0)error stop 'HALF: incomplete EIGENVAL header';end do
    read(unit,*,iostat=ios)nelectron,nk,nb
    if(ios/=0.or.nk<1.or.nb<1)error stop 'HALF: invalid EIGENVAL dimensions'
    set%nk=nk;set%full_count=nk;allocate(set%points(nk,3),set%weights(nk),set%multiplicities(nk),eigenvalues(nk,nb))
    do ik=1,nk
      do
        read(unit,'(A)',iostat=ios)line;if(ios/=0)error stop 'HALF: missing EIGENVAL k point'
        if(len_trim(line)>0)exit
      end do
      read(line,*,iostat=ios)set%points(ik,:),set%weights(ik)
      if(ios/=0)error stop 'HALF: invalid EIGENVAL k point'
      do ib=1,nb
        read(unit,'(A)',iostat=ios)line;if(ios/=0)error stop 'HALF: missing EIGENVAL band'
        read(line,*,iostat=ios)index,eigenvalues(ik,ib)
        if(ios/=0.or.index/=ib)error stop 'HALF: invalid EIGENVAL band row'
      end do
    end do
    close(unit)
    if(sum(set%weights)<=0)error stop 'HALF: invalid EIGENVAL weights'
    set%weights=set%weights/sum(set%weights);set%multiplicities=1;nelect=real(nelectron,dp)
  end subroutine

  subroutine compute_occupations(eigenvalues,weights,nelect,sigma,occupation,chemical_potential,band_energy,entropy_term)
    real(dp),intent(in)::eigenvalues(:,:),weights(:),nelect,sigma
    real(dp),allocatable,intent(out)::occupation(:,:)
    real(dp),intent(out)::chemical_potential,band_energy,entropy_term
    integer::nk,nb,nstates,s,i,j,position,stop,ik,ib
    integer,allocatable::order(:),state_k(:),state_b(:)
    real(dp),allocatable::energy(:)
    real(dp)::remaining,level,tolerance,capacity,fraction,filled,lo,hi,mid,count_e,x,p
    nk=size(eigenvalues,1);nb=size(eigenvalues,2)
    if(size(weights)/=nk)error stop 'HALF: occupation k-point weight mismatch'
    if(abs(sum(weights)-1.0_dp)>1e-10_dp)error stop 'HALF: k-point weights must sum to one'
    if(nelect< -1e-12_dp.or.nelect>2.0_dp*nb+1e-12_dp)error stop 'HALF: electrons do not fit in requested bands'
    if(sigma<0.or.sigma/=sigma)error stop 'HALF: sigma must be finite and non-negative'
    allocate(occupation(nk,nb));occupation=0.0_dp
    if(sigma==0.0_dp)then
      nstates=nk*nb;allocate(order(nstates),state_k(nstates),state_b(nstates),energy(nstates));s=0
      do ik=1,nk;do ib=1,nb;s=s+1;order(s)=s;state_k(s)=ik;state_b(s)=ib;energy(s)=eigenvalues(ik,ib);end do;end do
      do i=2,nstates
        s=order(i);j=i-1
        do while(j>=1)
          if(energy(order(j))<=energy(s))exit
          order(j+1)=order(j);j=j-1
        end do
        order(j+1)=s
      end do
      remaining=nelect;chemical_potential=minval(eigenvalues);position=1
      do while(position<=nstates.and.remaining>1e-13_dp)
        level=energy(order(position));stop=position+1;tolerance=1e-10_dp*max(1.0_dp,abs(level))
        do while(stop<=nstates)
          if(abs(energy(order(stop))-level)>tolerance)exit;stop=stop+1
        end do
        capacity=0.0_dp
        do i=position,stop-1;capacity=capacity+2.0_dp*weights(state_k(order(i)));end do
        fraction=min(1.0_dp,max(0.0_dp,remaining/capacity))
        do i=position,stop-1
          occupation(state_k(order(i)),state_b(order(i)))=2.0_dp*fraction
        end do
        filled=fraction*capacity;if(filled>0)chemical_potential=level;remaining=remaining-filled;position=stop
      end do
      if(remaining>1e-9_dp)error stop 'HALF: not enough bands for electron count'
      band_energy=sum(spread(weights,2,nb)*occupation*eigenvalues);entropy_term=0.0_dp;return
    end if
    lo=minval(eigenvalues)-max(10.0_dp,100.0_dp*sigma);hi=maxval(eigenvalues)+max(10.0_dp,100.0_dp*sigma)
    do i=1,256
      mid=0.5_dp*(lo+hi);count_e=electron_count(mid)
      if(count_e<nelect)then;lo=mid;else;hi=mid;end if
      if(hi-lo<1e-13_dp*max(1.0_dp,abs(mid)))exit
    end do
    chemical_potential=0.5_dp*(lo+hi);entropy_term=0.0_dp
    do ik=1,nk;do ib=1,nb
      x=max(-700.0_dp,min(700.0_dp,(eigenvalues(ik,ib)-chemical_potential)/sigma))
      p=1.0_dp/(1.0_dp+exp(x));occupation(ik,ib)=2.0_dp*p
      p=max(1e-300_dp,min(1.0_dp-1e-16_dp,p))
      entropy_term=entropy_term+2.0_dp*sigma*weights(ik)*(p*log(p)+(1-p)*log(1-p))
    end do;end do
    band_energy=sum(spread(weights,2,nb)*occupation*eigenvalues)
  contains
    real(dp) function electron_count(mu)result(value)
      real(dp),intent(in)::mu;integer::jk,jb;real(dp)::xx
      value=0
      do jk=1,nk;do jb=1,nb
        xx=max(-700.0_dp,min(700.0_dp,(eigenvalues(jk,jb)-mu)/sigma))
        value=value+weights(jk)*2.0_dp/(1.0_dp+exp(xx))
      end do;end do
    end function
  end subroutine

  subroutine ewald_energy(crystal,charges,total,real_energy,reciprocal_energy,self_energy,background_energy,eta,tolerance)
    type(crystal_t),intent(in)::crystal
    real(dp),intent(in)::charges(:)
    real(dp),intent(out)::total
    real(dp),intent(out),optional::real_energy,reciprocal_energy,self_energy,background_energy
    real(dp),intent(in),optional::eta,tolerance
    real(dp)::alpha,tol,scale,rcut,gcut,smin_r,smin_g,positions(crystal%nions,3),rv(3),gv(3),d(3),dist,gn
    real(dp)::re,ge,se,be,sr,si,phase,total_charge
    integer::rb,gb,n1,n2,n3,i,j
    if(size(charges)/=crystal%nions)error stop 'HALF: Ewald charge count mismatch'
    tol=1e-11_dp;if(present(tolerance))tol=tolerance
    if(tol<=0.or.tol>=1)error stop 'HALF: Ewald tolerance must lie between zero and one'
    alpha=sqrt(pi)/crystal%volume**(1.0_dp/3.0_dp);if(present(eta))alpha=eta
    if(alpha<=0.or.alpha/=alpha)error stop 'HALF: Ewald eta must be positive and finite'
    scale=sqrt(-log(tol));rcut=scale/alpha;gcut=2*alpha*scale
    smin_r=smallest_singular(crystal%lattice);smin_g=smallest_singular(crystal%reciprocal)
    rb=ceiling(rcut/smin_r)+1;gb=ceiling(gcut/smin_g)+1
    positions=matmul(crystal%positions,crystal%lattice);re=0.0_dp
    do i=1,crystal%nions;do j=1,crystal%nions
      do n1=-rb,rb;do n2=-rb,rb;do n3=-rb,rb
        rv=matmul([real(n1,dp),real(n2,dp),real(n3,dp)],crystal%lattice)
        d=positions(i,:)-positions(j,:)+rv;dist=sqrt(sum(d*d))
        if(dist>1e-13_dp.and.dist<=rcut)re=re+charges(i)*charges(j)*erfc(alpha*dist)/dist
      end do;end do;end do
    end do;end do;re=0.5_dp*felect*re;ge=0.0_dp
    do n1=-gb,gb;do n2=-gb,gb;do n3=-gb,gb
      gv=matmul([real(n1,dp),real(n2,dp),real(n3,dp)],crystal%reciprocal);gn=sqrt(sum(gv*gv))
      if(gn>1e-13_dp.and.gn<=gcut)then
        sr=0;si=0
        do i=1,crystal%nions
          phase=dot_product(gv,positions(i,:));sr=sr+charges(i)*cos(phase);si=si-charges(i)*sin(phase)
        end do
        ge=ge+exp(-gn*gn/(4*alpha*alpha))*(sr*sr+si*si)/(gn*gn)
      end if
    end do;end do;end do
    ge=2*pi*felect/crystal%volume*ge;se=-felect*alpha/sqrt(pi)*sum(charges*charges)
    total_charge=sum(charges);be=-felect*pi*total_charge**2/(2*alpha*alpha*crystal%volume)
    total=re+ge+se+be
    if(present(real_energy))real_energy=re;if(present(reciprocal_energy))reciprocal_energy=ge
    if(present(self_energy))self_energy=se;if(present(background_energy))background_energy=be
  end subroutine

  real(dp) function smallest_singular(matrix)result(value)
    real(dp),intent(in)::matrix(3,3)
    real(dp)::a(3,3),c,s,t,tau,temp,off
    integer::sweep,p,q,r
    a=matmul(matrix,transpose(matrix))
    do sweep=1,32
      off=abs(a(1,2))+abs(a(1,3))+abs(a(2,3));if(off<1e-15_dp*max(1.0_dp,maxval(abs(a))))exit
      do p=1,2;do q=p+1,3
        if(abs(a(p,q))>0)then
          tau=(a(q,q)-a(p,p))/(2*a(p,q));t=sign(1.0_dp,tau)/(abs(tau)+sqrt(1+tau*tau));c=1/sqrt(1+t*t);s=t*c
          do r=1,3
            if(r/=p.and.r/=q)then
              temp=a(r,p);a(r,p)=c*temp-s*a(r,q);a(p,r)=a(r,p);a(r,q)=s*temp+c*a(r,q);a(q,r)=a(r,q)
            end if
          end do
          temp=a(p,p);a(p,p)=c*c*temp-2*s*c*a(p,q)+s*s*a(q,q)
          a(q,q)=s*s*temp+2*s*c*a(p,q)+c*c*a(q,q);a(p,q)=0;a(q,p)=0
        end if
      end do;end do
    end do
    value=sqrt(max(0.0_dp,min(a(1,1),a(2,2),a(3,3))))
  end function
end module half_energy

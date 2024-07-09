program gradpsi
implicit none
!include "fftw_f77.i"
!	integer::i,n,k1,ireal,iimag,t
!	double precision::x,dx,L,kx,scale,pi,nu,dt,rk2
!	double precision,allocatable,dimension(:)::phi
!	double precision,allocatable,dimension(:)::gradphi!,gradsqrphi
!	integer,dimension(1) :: dim
!	integer*8 :: pfor,pinv
!	dim(1)=n
!	n=64
!	pi=4.0d0*datan(1.0d0)
!	scale=1/(sqrt(dfloat(n)))
!	allocate(phi(n+2))
	!!allocate(gradphi(n+2),gradsqrphi(n+2))
!	nu=0.00467
!	L=2*pi
!	dt=0.005
!	dx=L/(dfloat(n))
	!call rfftwnd_f77_create_plan(pfor,1,dim,FFTW_REAL_TO_COMPLEX,FFTW_MEASURE + FFTW_IN_PLACE)
	!call rfftwnd_f77_create_plan(pinv,1,dim,FFTW_COMPLEX_TO_REAL,FFTW_MEASURE + FFTW_IN_PLACE)
!	call rfftwnd_f77_create_plan(pfor,1,dim,FFTW_REAL_TO_COMPLEX,FFTW_MEASURE + FFTW_IN_PLACE)
!        call rfftwnd_f77_create_plan(pinv,1,dim,FFTW_COMPLEX_TO_REAL,FFTW_MEASURE + FFTW_IN_PLACE)
!	open(unit=12,file='input1',status='unknown')
	!do i=1,n
	!	x=(i-n/2)*dx
	!	phi(i)=cos(x)
	!	write(12,*) x,phi(i)
!	!enddo
!        do i=1,n
!                x=(i-n/2)*dx
!                phi(i)=exp(-x**2/(0.1))+sin(x)
!                write(12,*)x,phi(i)
!        enddo

!	close(12)
        
        include "fftw_f77.i"
        integer::i,n,k1,ireal,iimag,t
        double precision::length,kx,dx,dt,nu,pi,rk2,scale,x
        double precision,allocatable,dimension(:)::phi,gradphi,gradsqrphi
        integer,dimension(1) :: dim
        integer*8 :: pfor,pinv
        character*80::fname
        n=64
        allocate(phi(n+2),gradphi(n+2),gradsqrphi(n+2))
        phi=0.0
        gradphi=0.0
        gradsqrphi=0.0
        nu=0.00467
        dt=0.005
        dim(1)=n
        pi=4.0d0*datan(1.0d0)
        length=2*pi
        dx=length/(dfloat(n))
        scale=1.0d0/(sqrt(dfloat(n)))
        call rfftwnd_f77_create_plan(pfor,1,dim,FFTW_REAL_TO_COMPLEX,FFTW_MEASURE + FFTW_IN_PLACE)
        call rfftwnd_f77_create_plan(pinv,1,dim,FFTW_COMPLEX_TO_REAL,FFTW_MEASURE + FFTW_IN_PLACE)
        open(unit=32,file='1dinput',status='unknown')
        do i=1,n
                x=(i-n/2)*dx
                phi(i)=exp(-x**2/(0.1))!+sin(x)
                write(32,*)x,phi(i)
        enddo
        close(32)
        !call rfftwnd_f77_one_real_to_complex(pfor,phi,0)
        !phi=phi*scale

	do t=1,100
	call rfftwnd_f77_one_real_to_complex(pfor,phi,0)
	phi=phi*scale
	open(unit=14,file='output1',status='unknown')
	do i=1,n/2+1
		k1=i-1
		ireal=2*i-1
		iimag=2*i
		kx=(2*pi)/(length)*((k1))
		gradphi(ireal)=-kx*phi(iimag)
		gradphi(iimag)=kx*phi(ireal)
	enddo
	call rfftwnd_f77_one_complex_to_real(pinv,gradphi,0)
	gradphi=gradphi*scale
!	do i=1,n
!		x=(i-n/2)*dx
!		write(14,*) x,gradphi(i)
!	enddo
!	close(14)
	

	do i=1,n/2+1
		k1=i-1
		kx=dfloat(k1)
		rk2=dfloat(k1**2)
		ireal=2*i-1
		iimag=2*i
		gradsqrphi(ireal)=-rk2*phi(ireal)
		gradsqrphi(iimag)=-rk2*phi(iimag)
	enddo
	call rfftwnd_f77_one_complex_to_real(pinv,gradsqrphi,0)
	gradsqrphi=gradsqrphi*scale
!	
!	
	!do i=1,n
	!	phi(i)=phi(i)-phi(i)*gradphi(i)*dt+nu*gradsqrphi(i)*dt
	!enddo
!		
	call rfftwnd_f77_one_complex_to_real(pinv,phi,0)
	phi=phi*scale

        do i=1,n
                phi(i)=phi(i)-phi(i)*gradphi(i)*dt+nu*gradsqrphi(i)*dt
        enddo
        print*,mod(t,10)
        if(mod(t,10).eq.0)then
        write(fname,'(i4)')(t)
        open(unit=141, file='output_dir/file'//trim(adjustl(fname))//'.dat', status='unknown')
        do i=1,n
                write(141,*)i,phi(i)
        enddo
        close(141)
        endif
	enddo
        do i=1,n
                write(444,*)i,phi(i)
        enddo
        
		
end program gradpsi

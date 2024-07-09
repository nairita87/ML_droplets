program gradsqrin1d
	implicit none
	include "fftw_f77.i"
	integer::i,n,k1,ireal,iimag
	double precision::length,kx,dx,pi,rk2,scale,x
	double precision,allocatable,dimension(:)::phi,gradsqrphi
	integer,dimension(1) :: dim
	integer*8 :: pfor,pinv
	n=256
	allocate(phi(n+2),gradsqrphi(n+2))
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
		phi(i)=exp(-x**2/(0.1))+sin(x)
		write(32,*)x,phi(i)
	enddo
	close(32)
	call rfftwnd_f77_one_real_to_complex(pfor,phi,0)
	phi=phi*scale

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
	open(unit=35,file='1doutput',status='unknown')
	do i=1,n
		x=(i-n/2)*dx
		write(35,*)x,gradsqrphi(i)
	enddo
	close(35)
end program

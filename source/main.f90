!==========================================================================
program main
  implicit none
  real*8, parameter :: au2cm=219474d0
  real*8, parameter :: ang2rad=0.0174533
  character(3) :: sym
  real*8 :: anums
  real*8, dimension(:), allocatable :: e
  real*8, dimension(:,:), allocatable :: h,gref
  real*8, dimension(:,:,:), allocatable :: geom
  real*8, dimension(:,:,:), allocatable :: cg, dcg
  real*8 :: r,phi
  real*8 :: eshift
  integer :: i,j,k,total,ios,idx,natoms,nstates
  call initPotential()
  call getinfo(natoms, nstates)
  total=500
  allocate(e(nstates))
  allocate(h(nstates,nstates))
  allocate(cg(3*natoms,nstates,nstates))
  allocate(dcg(3*natoms,nstates,nstates))

  allocate(geom(3,natoms,total))
  allocate(gref(3,natoms))

  gref(:,1)=(/0.00000000,0.00000000,0.00000000/)
  gref(:,2)=(/1.91568930,-0.00000000,0.00000000/)
  gref(:,3)=(/-0.90158290,-1.62123135,0.00000000/)
  gref(:,4)=(/-0.90158290,1.62123135,0.00000000/)

  open(200,file='output')
  do i=1,10
    do j=0,30,5
      geom(:,:,i)=gref
      r=gref(1,2)+0.4*i
      phi=j*ang2rad
      geom(1,2,i)=r*dcos(phi)
      geom(3,2,i)=r*dsin(phi)
      call EvaluateSurfgen(geom(:,:,i),e,cg,h,dcg)
      write(200,"(2f8.3,f18.8)") r,dble(j),(e(1)+e(2))/2d0*au2cm
    end do
  end do

  stop
end
!==========================================================================

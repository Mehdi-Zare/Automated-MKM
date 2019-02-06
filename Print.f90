program Print
implicit none
integer  ,  parameter                   ::      dp = selected_real_kind(15, 307)
integer                                 ::      ierror                        !Deallocation status for allocatable arrays
integer                                 ::      i, numRxn                                       !The number of reactants, products, or transition states        
!character(20)                           ::      filename                                !The variable for the name of input file
real(dp) ,  allocatable, dimension (:)  ::      Grxn
!,Gact,Keq,kf,kr             !equilibrium constant, forward and reverse rate constant


numRxn=3

!allocate (RxnName(numRxn))
allocate (Grxn(numRxn))
!allocate (Gact(numRxn))
!allocate (Keq(numRxn))
!allocate (kf(numRxn))
!allocate (kr(numRxn))

open (unit=100, file = 'Grxn', status='old', action = 'read', iostat = ierror)
do i=1,numRxn
read(100,'(ES20.13)',iostat = ierror) Grxn(i)
end do
write(*,*) "Grxn=" , Grxn
close(100)

!open (unit=200, file = 'Gact', status='old', action = 'read', iostat = ierror)
!do i=1,numRxn
!read(200,*,iostat = ierror) Gact(i)
!end do
!close(200)


!open (unit=300, file = 'Keq', status='old', action = 'read', iostat = ierror)
!do i=1,numRxn
!read(300,*,iostat = ierror) Keq(i)
!end do
!close(300)

!open (unit=400, file = 'kf', status='old', action = 'read', iostat = ierror)
!do i=1,numRxn
!read(400,*,iostat = ierror) kf(i)
!end do
!close(400)

!open (unit=500, file = 'kr', status='old', action = 'read', iostat = ierror)
!do i=1,numRxn
!read(500,*,iostat = ierror) kr(i)
!end do
!close(500)





!write(*,*) "RxnName" , RxnName
!write(*,*) Grxn
!write(*,*) Gact
!write(*,*) Keq
!write(*,*) kr
!write(*,*) kf















end program Print

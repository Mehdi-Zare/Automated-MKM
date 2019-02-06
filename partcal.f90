! This is a FORTRAN code to calculate vibrational, Rotational and translational 
!partition function from the optimized CONTCAR file
! Written by Osman Mamun on 7.15.15
Program partcalc
implicit none
integer, parameter:: dp = selected_real_kind(15, 307)
!real,parameter::hh=4.135e-15,kb=8.617e-5
Character(len=1)::name(3)
Real(dp),allocatable::T(:),x(:),y(:),z(:),qrot(:),qvib(:),qtrans(:),lnqrot(:),lnqtrans(:),f(:),lnqvib(:),qv(:,:),lnq(:)
Real(dp)::multiplier,Nav,Pi,kB,h,Imatrix(3,3),IXX,IYY,IZZ,IXY,IYZ,IZX,A(3,3),E(3),xmass,ymass,qvi,ZPE
Real(dp)::zmass,msum,mxsum,mysum,mzsum,n(3),m(3),mu,sigma,d,P
Integer(dp)::LT,HT,IT,nTemp,ntotal,i,niter,nnn,Linearity,ndiff,nvib,j,ierror
Nav=6.02214076e23
h=6.62607015e-34
kB=1.380649e-23
Pi=3.14159265359
P=101325

!call intro()

open (unit=888, file = 'Gas-Data', status='old', action = 'read', iostat = ierror)
read(888,*,iostat = ierror)
read(888,*,iostat = ierror)
read(888,*,iostat = ierror) ndiff
read(888,*,iostat = ierror) Linearity
read(888,*,iostat = ierror) sigma

write(*,*) "sigma", sigma, "ndiff", ndiff, "linearity", Linearity
close(888)

open (unit=8881, file = 'Temp-File', status='old', action = 'read', iostat = ierror)
read(8881,*,iostat = ierror) LT
read(8881,*,iostat = ierror) HT
read(8881,*,iostat = ierror) IT
close(8881) 
!Write(*,*) 'Please give the temperature range you are interested in?'
!Write(*,*) 'Lower value of Temp.  Higher value of Temp.  Temp. interval'
!Read(*,*) LT, HT, IT
!LT=413
!HT=473
!IT=60
!Write(*,*) 'How many different types of atoms do you have in your system? (Maximum is 3)'
!Read (*,*) ndiff  
!Write(*,*) 'What is the linearity of the molecule? 1=linear 2=non-linear'
!Read (*,*) Linearity
nTemp=((HT-LT)/IT)+1
allocate (T(nTemp))
allocate (qrot(nTemp))
allocate (lnqrot(nTemp))
allocate (qvib(nTemp))
allocate (lnqvib(nTemp))
allocate (qtrans(nTemp))
allocate (lnqtrans(nTemp))
allocate (lnq(nTemp))
do i=1,nTemp
  T(i)=LT+(i-1)*IT
end do
!Write(*,*) 'What is the symmetry number of the molecule involved?'
!read(*,*) sigma

Open(unit=10,file='CONTCAR',status='old')
    read(10,*)
    read(10,*)
    read(10,*) multiplier
    read(10,*)
    read(10,*)
    if (ndiff==1) then
          read(10,*) name(1)
          read(10,*) n(1)
          n(2)=0
          n(3)=0
          m(2)=0
          m(3)=0
          ntotal=n(1)
    else if (ndiff==2) then
          read(10,*) name(1),name(2)
          read(10,*) n(1),n(2)
          ntotal=n(1)+n(2)
          n(3)=0
          m(3)=0
    else
          read(10,*) name(1),name(2),name(3)
          read(10,*) n(1),n(2),n(3)
          ntotal=n(1)+n(2)+n(3)
    end if
    allocate (x(ntotal))
    allocate (y(ntotal))
    allocate (z(ntotal))
    read(10,*)
    read(10,*)
    do i=1,ntotal
      read(10,*) x(i),y(i),z(i)
      x(i)=x(i)*multiplier
      y(i)=y(i)*multiplier
      z(i)=z(i)*multiplier
    end do
close(10)
if (Linearity==1) then
    nvib=3*ntotal-5
	else 
	nvib=3*ntotal-6
end if
allocate (f(nvib))
allocate (qv(ntotal,nvib))
ZPE=0
Open(unit=20,file='freq.dat',status='old')
   do i=1,nvib
      read(20,*) f(i)
      ZPE=ZPE+0.5*4.135e-15*3e10*f(i)
   end do
close(20)

do i=1,nTemp
   qvi=1
   do j=1,nvib
      qvi=qvi*(1/(1-exp(-4.135e-15*f(j)*3e10/8.617e-5/T(i))))
   end do
   qvib(i)=qvi
   lnqvib(i)=log(qvib(i))
 end do


do i=1,ndiff
  if (name(i)=='C') then
    m(i)=12.0107
  else if (name(i)=='O') then
    m(i)=15.9994
  else
    m(i)=1.00794
  end if
end do
mxsum=0
mysum=0
mzsum=0
do i=1,ntotal
  if (i<=n(1)) then
    mxsum=mxsum+m(1)*x(i)
    mysum=mysum+m(1)*y(i)
    mzsum=mzsum+m(1)*z(i)
  else if (i>n(1) .and. i<=n(1)+n(2)) then
    mxsum=mxsum+m(2)*x(i)
    mysum=mysum+m(2)*y(i)
    mzsum=mzsum+m(2)*z(i)
  else
    mxsum=mxsum+m(3)*x(i)
    mysum=mysum+m(3)*y(i)
    mzsum=mzsum+m(3)*z(i)
  end if
end do
msum=m(1)*n(1)+m(2)*n(2)+m(3)*n(3)
xmass=mxsum/msum
ymass=mysum/msum
zmass=mzsum/msum
if (Linearity==2) then
    IXX=0
    IYY=0
    IZZ=0
    IXY=0
    IYZ=0
    IZX=0
    do i=1,ntotal
       if (i<=n(1)) then
          IXX=IXX+m(1)*((y(i)-ymass)**2+(z(i)-zmass)**2)
          IYY=IYY+m(1)*((z(i)-zmass)**2+(x(i)-xmass)**2)
          IZZ=IZZ+m(1)*((x(i)-xmass)**2+(y(i)-ymass)**2)
          IXY=IXY+m(1)*(x(i)-xmass)*(y(i)-ymass)
          IYZ=IYZ+m(1)*(y(i)-ymass)*(z(i)-zmass)
          IZX=IZX+m(1)*(z(i)-zmass)*(x(i)-xmass)
       else if (i>n(1) .and. i<=n(1)+n(2)) then
          IXX=IXX+m(2)*((y(i)-ymass)**2+(z(i)-zmass)**2)
          IYY=IYY+m(2)*((z(i)-zmass)**2+(x(i)-xmass)**2)
          IZZ=IZZ+m(2)*((x(i)-xmass)**2+(y(i)-ymass)**2)
          IXY=IXY+m(2)*(x(i)-xmass)*(y(i)-ymass)
          IYZ=IYZ+m(2)*(y(i)-ymass)*(z(i)-zmass)
          IZX=IZX+m(2)*(z(i)-zmass)*(x(i)-xmass)
       else
          IXX=IXX+m(3)*((y(i)-ymass)**2+(z(i)-zmass)**2)
          IYY=IYY+m(3)*((z(i)-zmass)**2+(x(i)-xmass)**2)
          IZZ=IZZ+m(3)*((x(i)-xmass)**2+(y(i)-ymass)**2)
          IXY=IXY+m(3)*(x(i)-xmass)*(y(i)-ymass)
          IYZ=IYZ+m(3)*(y(i)-ymass)*(z(i)-zmass)
          IZX=IZX+m(3)*(z(i)-zmass)*(x(i)-xmass)
       end if
    end do
    Imatrix(1,1)=IXX
    Imatrix(2,2)=IYY
    Imatrix(3,3)=IZZ
    Imatrix(1,2)=IXY
    Imatrix(1,3)=IZX
    Imatrix(2,3)=IYZ
    Imatrix(2,1)=IXY
    Imatrix(3,2)=IYZ
    Imatrix(3,1)=IZX  
    A=Imatrix
    nnn=3
    call Eigenval(A,nnn,E,niter)
    do i=1,ntemp
      qrot(i)=(pi*pi*pi*pi**0.5/sigma)*((8*kB*T(i)*1e-23/Nav)**1.5)*((E(1)*E(2)*E(3))**0.5)/h/h/h
    end do
else 
  if (ndiff==1) then
    do i=1,ntemp
         d=((x(2)-x(1))**2+(y(2)-y(1))**2+(z(2)-z(1))**2)
         qrot(i)=(1/sigma)*(8*pi*pi*T(i)*kB/Nav/h/h)*(m(1)/2)*d*1e-23
    end do
  else 
    mu=m(1)*m(2)/(m(1)+m(2))
    do i=1,ntemp
         d=((x(2)-x(1))**2+(y(2)-y(1))**2+(z(2)-z(1))**2)
         qrot(i)=(1/sigma)*(8*pi*pi*T(i)*kB/Nav/h/h)*(mu/4)*d*1e-23
    end do 
  end if
end if
lnqrot=log(qrot)
qtrans=((2*pi*msum*kB*T/h/h/1000/Nav)**1.5)*kB*T/P
lnqtrans=log(qtrans)
lnq=lnqtrans+lnqrot+lnqvib
!Write(*,*) lnqrot,lnqtrans
     

Open(unit=500,file='rot-trans-vib.txt',status='new')
Write(500,501)
501 format(7X,"T",15X,"qrot",16X,"qtrans",16X,"qvib",16X,"ln(qrot)",16X,"ln(qtrans)",16X,"ln(qvib)",16X,"lnq")
Write(500,502)
502 format (3X,'===========================================================',  &
               '===========================================================',  &
			   '===========================================================')
do i=1,nTemp
  Write(500,503) T(i),qrot(i),qtrans(i),qvib(i),lnqrot(i),lnqtrans(i),lnqvib(i),lnq(i)
  503 format (3X,F8.2,3X,7(ES20.13,3X))
end do
Write(500,504)
504 format (3X,'===========================================================',  &
               '===========================================================',  &
			   '===========================================================')
Write(500,*) 'ZPE=',ZPE,'eV'

close(500) 
 

Open(unit=8882,file='lnq',status='new')
do i=1,nTemp
  Write(8882,8883) lnq(i)
  8883 format (ES20.13)
end do

close(8882)




end program partcalc
!********************************************************************************!


!************************ Give a intro of the code ******************************!
   Subroutine intro()
   integer:: status=0
   character(10)::answer
   write(*,100)
         100 format(1X,"This code requires a POSCAR formatted file that contains coordinate of" ,&
                       "carbon, oxygen and hydrogen. " ,& 
                       "If you want to use it for atoms other than the" ,&
                       "three specified, please talk to the appropriate person to make arrangement or " ,&
                       "you can also implement the required change by editing this code")
600  write (*,200)
           200 format(1X,"Do you still want to use this code? type yes or no")
           read(*,*) answer
           if (answer == 'no') then
                write (*,300)
                      300 format(1X,"Aborting code-partfunc")
                call exit(status)
           elseif (answer == 'yes') then
               write(*,400)
                     400 format(1X,"Starting code-partfunc")
           else
               write(*,500)
                    500 format(1X,"WTF, I couldn't recognize the command")
               go to 600
           end if
   end subroutine intro
!*********************************************************************************!

!************************* Eigenvalue subroutine  ********************************!
Subroutine eigenval(A,nnn,E,niter)
implicit none
!This program will find eigenvalue of a n by n real symmetric matrix using jacobi's algorithm
!The main idea is to rotate the matrix until all the off diagonal elements are zero
integer, parameter:: dp = selected_real_kind(15, 307)
Integer(dp),intent(in):: nnn
real(dp),intent(inout)::A(nnn,nnn)
Real(dp),intent(out)::E(nnn)
real(dp)::ADM(nnn,nnn),Theta,t,c,s,S1(nnn,nnn),S1T(nnn,nnn),D1(nnn,nnn),D(nnn,nnn),Amax
integer(dp)::m(2),i,j
integer(dp),intent(out)::niter

niter=0

do
  niter=niter+1
!************************* Find maximum Aij ******************************!
     do i=1,nnn
        do j=1,nnn
           if (i ==j ) then
              ADM(i,j)=0
           else
              ADM(i,j)=abs(A(i,j))
           end if
        end do
     end do

    
     Amax=maxval(ADM)
     do i=1,nnn
       do j=1,nnn
         if(abs(ADM(i,j))==Amax) then
           m(1)=j
           m(2)=i
         end if
       end do
     end do 
     if (Amax<0.0001) exit
     Theta=(A(m(1),m(1))-A(m(2),m(2)))/2/A(m(2),m(1))
     if (Theta <0) then
        t=-1/(abs(Theta)+(Theta*Theta+1)**0.5)
     else
        t=1/(abs(Theta)+(Theta*Theta+1)**0.5)
     end if
     c=1/(t**2+1)**0.5
     s=c*t

     do i =1,nnn
       do j=1,nnn
         if (i==j) then
           S1(i,j)=1
         else
           S1(i,j)=0
         end if
       end do
     end do
     S1(m(1),m(1))=c
     S1(m(2),m(2))=c
     S1(m(2),m(1))=s
     S1(m(1),m(2))=-s
     S1T=transpose(S1)
     D1=matmul(S1T,A)
     D=matmul(D1,S1)
     do i=1,nnn
       do j=1,nnn
         if (abs(D(i,j))<0.0001) then
           D(i,j)=0
         else
           D(i,j)=D(i,j)
         end if
       end do
     end do
    A=D
end do
do i=1,nnn
  do j=1,nnn
    if (i==j) then
       E(i)=A(i,j)
    end if
  end do
end do

End subroutine

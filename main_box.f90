!============================================================================================
!    SDBOX - FORTRAN90 implementation of a Derivative-Free algorithm for bound 
!    constrained optimization problems 
!    Copyright (C) 2011  G.Liuzzi, S. Lucidi
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!    S. Lucidi, M. Sciandrone. A Derivative-Free Algorithm for Bound Constrained Optimization, 
!    Computational Optimization and Applications, 21(2): 119-142 (2002)
!    DOI: 10.1023/A:1013735414984
!
!============================================================================================

program main_box
      implicit none

	  integer   :: n
!-----------------------------------------------------------------------------	  
	  integer	:: i, istop, icheck
	  real		:: tbegin, tend

      real*8, allocatable :: x(:),bl(:),bu(:)


!-----------------------------------------------------------------------------

	  integer ::            num_funct,num_iter
	  real*8             :: f,alfamax,delta 
	  real*8			   :: fob,fapp
	  real*8			   :: violiniz, violint, finiz, fex
      real*8             :: alfa_stop
      integer            :: nf_max,iprint

!------------------------------------------------------------------------------

!-----------------------------------------------------------------------
!      Starting point and bound calculation
!-----------------------------------------------------------------------

      write(1,*) ' '

	  call setdim(n)

	  allocate(x(n), bl(n), bu(n))

	  call setbounds(n,bl,bu)	
	  call startp(n,x)
 
	  do i=1,n
		write(1,*) ' x(',i,')=',x(i)
	  enddo


	  do i=1,n

         if((x(i).lt.bl(i)).or.(x(i).gt.bu(i))) then
		   write(*,*) ' punto iniziale non in box'
		   stop
		 endif
     
	  enddo

2002  format(2d20.10)

!-----------------------------------------------------------------------
!     calculate starting point violation 
!-----------------------------------------------------------------------

    violiniz=0.0d0
    
	do i=1,n
       violiniz=max(violiniz,x(i)-bu(i),bl(i)-x(i)) 
	end do
    

	  	call funct(n,x,fob)


        finiz=fob
   
        write(*,*) ' ------------------------------------------------- '
        write(1,*) ' ------------------------------------------------- '

        write(*,*) ' objective function at xo = ',fob
        write(1,*) ' objective function at xo = ',fob

!       ---- objective function value ----

        write(*,*) ' ------------------------------------------------- '
        write(1,*) ' ------------------------------------------------- '


!      ---- x(i) =  i-th variable value----

        write(*,*) ' ------------------------------------------------- '
        write(1,*) ' ------------------------------------------------- ' 		      
	    do i=1,n
		   write(*,*) 'xo(',i,') =',x(i)
		   write(1,*) 'xo(',i,') =',x(i)
        enddo








!------------------------------------

      
      call funct(n,x,fob)
      

	  write(*,*) 'fob = ',fob

	   num_funct   = 1 
	   num_iter    = 1



!      set minimum stepsize to 0.001
1	   alfa_stop=1.d-6
!      set maximum number of function evaluations
  	   nf_max=20000
!      set output verbosity
       iprint=0

	   open(78,file='solution.dat',status='replace')

 	   call cpu_time(tbegin)
	   
       call sd_box(n,x,f,bl,bu,alfa_stop,nf_max,num_iter,num_funct,iprint,istop)

	   call cpu_time(tend)

	   call funct(n,x,fob) 


      !-----------------------------------------------------------------------

       write(2,987) n,violiniz,finiz,fob,num_funct,num_iter
 
 987 format(' & ', i3,' & ',es9.2,' & ',es14.7,' & ',es14.7,' & ',i5,' & ',i5,'\\')

!	   close(77)
      
	   write(*,*) '------------------------------------------------------------------------------'     
	   if(istop.eq.1) then
         write(*,*)  ' END - stopping criterion satisfied '
	   endif
       if(istop.eq.2) then
         write(*,*)  ' END - maximum number of function calculation reached  =',nf_max
	   endif

       if(istop.eq.3) then
         write(*,*)  ' END -  maximum number of iterations reached =',nf_max
	   endif

       write(78,*) 'objective function=', fob

	   write(78,*) ''

	   write(78,*) 'variables'
	   do i=1,n
         write(78,*)' x(',i,')=',x(i)
       enddo
	   do i=1,n
         if((x(i)-bl(i)).le.1.d-24) write(78,*)' la variabile x(',i,') al suo limite inferiore'
		 if((bu(i)-x(i)).le.1.d-24) write(78,*)' la variabile x(',i,') al suo limite superiore'
       enddo
	   write(78,*) ''
       write(78,*) 'CPU time=', tend-tbegin

	   close(78)

	   write(*,*) ' total time:',tend-tbegin
	   write(*,*) '------------------------------------------------------------------------------'  

        write(*,*) ' ------------------------------------------------- '
        write(1,*) ' ------------------------------------------------- '

!       ---- fo = objective function value ----

        write(*,*) ' ------------------------------------------------- '
        write(1,*) ' ------------------------------------------------- '

        write(*,*) ' objective function = ',fob
        write(1,*) ' objective function = ',fob


!      ---- x(i) = i-th variable ----

        write(*,*) ' ------------------------------------------------- '
        write(1,*) ' ------------------------------------------------- ' 		      
	    do i=1,n
		   write(*,*) 'x(',i,') =',x(i)
		   write(1,*) 'x(',i,') =',x(i)
        enddo

!      ---- nftot = number of function evaluations ----

        write(*,*) ' ------------------------------------------------- '
        write(1,*) ' ------------------------------------------------- '

        write(*,*) ' number of function evaluations = ',num_funct 
        write(1,*) ' number of function evaluations = ',num_funct     		    
	   
	   do i=1,n
         if((x(i)-bl(i)).le.1.d-24) write(*,*)' variable x(',i,') is at lower bound'
		 if((bu(i)-x(i)).le.1.d-24) write(*,*)' variable x(',i,') is at upper bound'
         if((x(i)-bl(i)).le.1.d-24) write(1,*)' variable x(',i,') is at lower bound'
		 if((bu(i)-x(i)).le.1.d-24) write(1,*)' variable x(',i,') is at upper bound'
       enddo

	   deallocate(x,bl,bu)

end program main_box

      subroutine sd_box(n,x,f,bl,bu,alfa_stop,nf_max,ni,nf,iprint,istop)
      implicit none
	  logical :: cambio_eps
      integer :: n,i,j,i_corr,nf,ni,nf_max
      integer :: num_fal,istop
      integer :: iprint,i_corr_fall
	  integer :: flag_fail(n)

      real*8 :: x(n),z(n),d(n)
      real*8 :: alfa_d(n),alfa,alfa_max, alfa_d_old
      real*8 :: f,fz , eta
	  real*8 :: bl(n),bu(n),alfa_stop,maxeps

!     values of f calculated on a n+1 simplex

      real*8 :: fstop(n+1)


!     num_fal number of failures

!     i_corr is the index of the current direction


!     initialization

	  eta = 1.d-6

      flag_fail=0

	  num_fal=0

      istop = 0

      fstop=0.d0

!     ---- choice of the starting stepsizes along the directions --------

      do i=1,n

           alfa_d(i)=dmax1(1.d-3,dmin1(1.d0,dabs(x(i))))
      
           if(iprint.ge.1) then
              write(*,*) ' alfainiz(',i,')=',alfa_d(i)
              write(1,*) ' alfainiz(',i,')=',alfa_d(i)
           endif

      end do

      do i=1,n      
        d(i)=1.d0 
      end do
!     ---------------------------------------------------------------------  
     
      
      call funct(n,x,f)
      

	  nf=nf+1

	  i_corr=1

      fstop(i_corr)=f

      do i=1,n
	    z(i)=x(i)
      end do

      if(iprint.ge.2) then
        write(*,*) ' ----------------------------------'
        write(1,*) ' ----------------------------------'
        write(*,*) ' finiz =',f
        write(1,*) ' finiz =',f
        do i=1,n
          write(*,*) ' xiniz(',i,')=',x(i)
          write(1,*) ' xiniz(',i,')=',x(i)
        enddo
      endif

!---------------------------   
!     main loop
!---------------------------

      do 

         if(iprint.ge.0) then
         
           write(*,100) ni,nf,f,alfa_max
           write(1,100) ni,nf,f,alfa_max
100        format(' ni=',i4,'  nf=',i5,'   f=',d12.5,'   alfamax=',d12.5)
         endif
         if(iprint.ge.2) then
	       do i=1,n
                write(*,*) ' x(',i,')=',x(i)
                write(1,*) ' x(',i,')=',x(i)
            enddo
         endif
!-------------------------------------
!    sampling along coordinate i_corr
!-------------------------------------
                call linesearchbox_cont(n,x,f,d,alfa,alfa_d,z,fz,i_corr,num_fal,&
                           alfa_max,i_corr_fall,iprint,bl,bu,ni,nf)

         if(dabs(alfa).ge.1.d-12) then
		    
			flag_fail(i_corr)=0
		               
            x(i_corr) = x(i_corr)+alfa*d(i_corr)
            f=fz
 	        fstop(i_corr)=f
			     
            num_fal=0
            ni=ni+1
      
         else

			flag_fail(i_corr)=1

	        if(i_corr_fall.lt.2) then 

		      fstop(i_corr)=fz         

              num_fal=num_fal+1
              ni=ni+1

	        endif

	     end if

		 z(i_corr) = x(i_corr)

         if(i_corr.lt.n) then
            i_corr=i_corr+1
         else
            i_corr=1
         end if 

         call stop(n,alfa_d,istop,alfa_max,nf,ni,fstop,f,alfa_stop,nf_max,flag_fail)

         if (istop.ge.1) then
			write(*,100) ni,nf,f,alfa_max
			write(1,100) ni,nf,f,alfa_max
			exit
		 endif

      enddo
      return
    


      end
        


!     #######################################################

      subroutine stop(n,alfa_d,istop,alfa_max,nf,ni,fstop,f,alfa_stop,nf_max, flag_fail)
      implicit none
      
      integer :: n,istop,i,nf,ni,nf_max
	  integer :: flag_fail(n)

      real*8 :: alfa_d(n),alfa_max,fstop(n+1),ffstop,ffm,f,alfa_stop

	  logical :: test

      istop=0

      alfa_max=0.0d0


      do i=1,n
          if(alfa_d(i).gt.alfa_max) then
            alfa_max=alfa_d(i)
          end if
      end do
     

      if(ni.ge.(n+1)) then
        ffm=f
        do i=1,n
          ffm=ffm+fstop(i)
        enddo
        ffm=ffm/dfloat((n+1))

        ffstop=(f-ffm)*(f-ffm)
        do i=1,n
           ffstop=ffstop+(fstop(i)-ffm)*(fstop(i)-ffm)
        enddo
 
        ffstop=dsqrt(ffstop/dfloat(n+1))



	  endif


        
      if(alfa_max.le.alfa_stop) then
	    test=.true.
        if (test.eqv..true.) then
		   istop = 1
		end if
        
	  end if
      


      if(nf.gt.nf_max) then
        istop = 2
      end if

      if(ni.gt.nf_max) then
        istop = 3
      end if

      return

      end




!     *********************************************************
!     *         
!     *                 Continuous Linesearch
!     *
!     *********************************************************
           
 
      subroutine linesearchbox_cont(n,x,f,d,alfa,alfa_d,z,fz,i_corr,num_fal,&
                                 alfa_max,i_corr_fall,iprint,bl,bu,ni,nf)
      
      implicit none

      integer :: n,i_corr,nf
      integer :: i,j
      integer :: ni,num_fal
      integer :: iprint,i_corr_fall
	  integer :: ifront,ielle
      real*8 :: x(n),d(n),alfa_d(n),z(n),bl(n),bu(n)
      real*8 :: f,alfa,alfa_max,alfaex, fz,gamma, gamma_int
      real*8 :: delta,delta1,fpar,fzdelta

	  
	  gamma=1.d-6

      delta =0.5d0
      delta1 =0.5d0

      i_corr_fall=0

	  ifront=0

!     index of current direction

      j=i_corr

	  if(iprint.ge.1) then
			write(*,*) 'variabile continua  j =',j,'    d(j) =',d(j),' alfa=',alfa_d(j)
			write(1,*) 'variabile continua  j =',j,'    d(j) =',d(j),' alfa=',alfa_d(j)
	  endif


	  if(dabs(alfa_d(j)).le.1.d-3*dmin1(1.d0,alfa_max)) then
			alfa=0.d0
			if(iprint.ge.1) then
				 write(*,*) '  alfa piccolo'
				 write(1,*) '  alfa piccolo'
				 write(*,*) ' alfa_d(j)=',alfa_d(j),'    alfamax=',alfa_max
				 write(1,*) ' alfa_d(j)=',alfa_d(j),'    alfamax=',alfa_max
			endif
			return
	  endif
      
!     choice of the direction

	  do ielle=1,2

		 if(d(j).gt.0.d0) then

		     if((alfa_d(j)-(bu(j)-x(j))).lt.(-1.d-6)) then                 
   			    alfa=dmax1(1.d-24,alfa_d(j))
			 else
			    alfa=bu(j)-x(j)
				ifront=1
				if(iprint.ge.1) then
					   write(*,*) ' point on the boundary. *'
					   write(1,*) ' point on the boundary. *'
				endif
			 endif

		  else

			 if((alfa_d(j)-(x(j)-bl(j))).lt.(-1.d-6)) then
			    alfa=dmax1(1.d-24,alfa_d(j))
			 else
				alfa=x(j)-bl(j)
				ifront=1
				if(iprint.ge.1) then
					   write(*,*) ' punto espan. sulla front. *'
					   write(1,*) ' punto espan. sulla front. *'
				endif
			 endif

		  endif

		  if(dabs(alfa).le.1.d-3*dmin1(1.d0,alfa_max)) then
  
			 d(j)=-d(j)
			 i_corr_fall=i_corr_fall+1
			 ifront=0

			 if(iprint.ge.1) then
				   write(*,*) ' direzione opposta per alfa piccolo'
				   write(1,*) ' direzione opposta per alfa piccolo'
				   write(*,*) ' j =',j,'    d(j) =',d(j)
				   write(1,*) ' j =',j,'    d(j) =',d(j)
				   write(*,*) ' alfa=',alfa,'    alfamax=',alfa_max
				   write(1,*) ' alfa=',alfa,'    alfamax=',alfa_max
			  endif
			  alfa=0.d0
			  cycle

		  endif

		  alfaex=alfa

		  z(j) = x(j)+alfa*d(j)
    
		 
	      call funct(n,z,fz)
		  

		  nf=nf+1

		  if(iprint.ge.1) then
				write(*,*) ' fz =',fz,'   alfa =',alfa
				write(1,*) ' fz =',fz,'   alfa =',alfa
		  endif
		  if(iprint.ge.2) then
			  do i=1,n
				  write(*,*) ' z(',i,')=',z(i)
				  write(1,*) ' z(',i,')=',z(i)
			  enddo
		  endif

		  fpar= f-gamma*alfa*alfa

!         test on the direction

		  if(fz.lt.fpar) then

!         expansion step

			 do

		   	   if((ifront.eq.1)) then

			         if(iprint.ge.1) then
				         write(*,*) ' accetta punto sulla frontiera fz =',fz,'   alfa =',alfa
				         write(1,*) ' accetta punto sulla frontiera fz =',fz,'   alfa =',alfa
			         endif
				     alfa_d(j)=delta*alfa

				     return

				 end if

				 if(d(j).gt.0.d0) then
							
					 if((alfa/delta1-(bu(j)-x(j))).lt.(-1.d-6)) then
						 alfaex=alfa/delta1
					 else
						 alfaex=bu(j)-x(j)
						 ifront=1
						 if(iprint.ge.1) then
							write(*,*) ' punto espan. sulla front.'
							write(1,*) ' punto espan. sulla front.'
						 endif
					 end if

				 else

					 if((alfa/delta1-(x(j)-bl(j))).lt.(-1.d-6)) then
						 alfaex=alfa/delta1
					 else
						 alfaex=x(j)-bl(j)
						 ifront=1
						 if(iprint.ge.1) then
							write(*,*) ' punto espan. sulla front.'
							write(1,*) ' punto espan. sulla front.'
						 endif
					 end if

				 endif
						 
				 z(j) = x(j)+alfaex*d(j) 
				   
     
				
			     call funct(n,z,fzdelta)
							      
				
				 nf=nf+1

				 if(iprint.ge.1) then
					  write(*,*) ' fzex=',fzdelta,'  alfaex=',alfaex  
					  write(1,*) ' fzex=',fzdelta,'  alfaex=',alfaex
				 endif
				 if(iprint.ge.2) then
					  do i=1,n
						 write(*,*) ' z(',i,')=',z(i)
						 write(1,*) ' z(',i,')=',z(i)
					  enddo
				 endif

				 fpar= f-gamma*alfaex*alfaex

				 if(fzdelta.lt.fpar) then

					 fz=fzdelta
					 alfa=alfaex

				 else               
					 alfa_d(j)=delta*alfa
			         if(iprint.ge.1) then
				         write(*,*) ' accetta punto fz =',fz,'   alfa =',alfa
				         write(1,*) ' accetta punto fz =',fz,'   alfa =',alfa
			         endif
					 return
				 end if

		     enddo 

		  else   !opposite direction    

			 d(j)=-d(j)
			 ifront=0

			 if(iprint.ge.1) then
				   write(*,*) ' direzione opposta'
				   write(1,*) ' direzione opposta'
				   write(*,*) ' j =',j,'    d(j) =',d(j)
				   write(1,*) ' j =',j,'    d(j) =',d(j)
			 endif

		  endif       ! test on the direction
			  
	  enddo       

	  if(i_corr_fall.eq.2) then
			 alfa_d(j)=alfa_d(j)
	  else
			 alfa_d(j)=delta*alfa_d(j)
	  end if

	  alfa=0.d0

	  if(iprint.ge.1) then
			write(*,*) ' failure along the direction'
			write(1,*) ' failure along the direction'
	  endif

	  return      
	  
      end

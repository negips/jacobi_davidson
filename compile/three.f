c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,ieg)

      implicit none  
  
      include 'SIZE'
      include 'PARALLEL'
!      include 'TOTAL'
      include 'INPUT'
      include 'TSTEP'
      include 'NEKUSE'

      integer e,ix,iy,iz,ieg


      e = gllel(ieg)
      utrans = 1.0   
      udiff = param(2)

      if (ifield .eq. 2) then
        udiff = param(8)
      endif
 
      return
      end
c-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,ieg)

      implicit none        
  
      include 'SIZE'
      include 'NEKUSE'

      integer ix,iy,iz,ieg

      ffx = 0.00
      ffy = 0.0
      ffz = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      qvol =  0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userchk

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'SOLN'
      include 'MVGEOM'
      include 'INPUT'
      include 'TSTEP'
      include 'MASS'
      include 'PARALLEL'      ! nelgv

      include 'F3D'
      include 'FS_ALE'
      include 'TEST'
      include 'DOMAIN'
      include 'WZ'

      integer lt
      parameter (lt=lx1*ly1*lz1*lelt)

      integer i,j,k,e,n,n2

      n  = lx1*ly1*lz1*nelv
      n2 = lx2*ly2*lz2*nelv

      if (istep.eq.0) then


        call frame_start

12    format(A4,2x,16(E12.5,2x))

      endif 

      call frame_monitor

      call chkpt_main


      if (istep.eq.nsteps.or.lastep.eq.1) then
        call frame_end
      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)

      implicit none

      include 'SIZE'
      include 'NEKUSE'
      include 'GEOM'

      integer ix,iy,iz,iside,ieg
      real pi

      integer jp
      common /ppointr/ jp

      real rmid

      real glmin,glmax
      integer n

      real rad1,rad2,omega1,omega2
      real a1,a2
      common /cylindrical/ rad1,rad2,omega1,omega2,a1,a2

      ux   = 0.0
      uy   = 0.0
      uz   = 0.0

      rmid = (rad1+rad2)/2
      if (y.lt.(rmid)) uz = omega1*rad1

      return
      end
c-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)

      implicit none

      include 'SIZE'
      include 'INPUT'         ! if3d
      include 'PARALLEL'
      include 'NEKUSE'
      include 'GEOM'

      include 'F3D'

      integer ix,iy,iz,ieg
      real pi

      integer jp
      common /ppointr/ jp

      real fcoeff(3)
      real xl(3)
      real mth_ran_dst

      logical ifcouette
      logical ifpoiseuille
      logical iftaylor
      logical iftestmvb

      real rad1,rad2,omega1,omega2
      real a1,a2
      common /cylindrical/ rad1,rad2,omega1,omega2,a1,a2

      real rmid,y0,z0
      real rad

      ifcouette         = .false.
      ifpoiseuille      = .false.
      iftaylor          = .true.
      iftestmvb         = .false.

      pi = 4.0*atan(1.0)


      if (ifpoiseuille) then
        ux = 1.0 - y**2
        uy = 0.
        uz = 0.0 + 0.0
      elseif (ifcouette) then
        ux = 0.0 + 1.0*y
        uy = 0.
        uz = 0.0 + 0.0
      elseif (iftaylor) then
        ux = 0.0 ! sin(pi*(y-1.0)/0.8)
        uy = 0.
        uz = a1*y + a2/y
      elseif (iftestmvb) then
        rmid = (rad1+rad2)/2.0
        ux   = -uparam(1)*exp(-((y-rmid)/0.25)**2)
        uy   = -0.1*uparam(1)*exp(-((y-rmid)/0.25)**2)
        uz   = 0.1*(a1*y + a2/y)
      endif  


      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat   ! This routine to modify element vertices

      implicit none  
  
      include 'SIZE'      ! _before_ mesh is generated, which 
      include 'INPUT'
      include 'GEOM'
      include 'TSTEP'
      include 'PARALLEL'

!      ifaxis = .true.   ! just for initialization
      param(42)=0       ! 0: GMRES (nonsymmetric), 1: PCG w/o weights
      param(43)=0       ! 0: Additive multilevel (param 42=0), 1: Original 2 level
      param(44)=0       ! 0: E based Schwartz (FEM), 1: A based Schwartz



      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2   ! This routine to modify mesh coordinates

      implicit none

      include 'SIZE'
      include 'INPUT'      ! cbc
      include 'PARALLEL'
      include 'GEOM'

      integer iel,ifc
      integer n

      if (if3d.and..not.ifcyclic) then
        n = lx1*ly1*lz1*nelv
!        call cmult(zm1,0.20,n)
      endif  

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat3

      implicit none        

      include 'SIZE'
      include 'SOLN'    ! tmult
      include 'INPUT'
      include 'GEOM'


      return
      end
c-----------------------------------------------------------------------








c automatically added by makenek
      subroutine usrdat0() 

      return
      end

c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end

c automatically added by makenek
      subroutine userqtl

      call userqtl_scig

      return
      end

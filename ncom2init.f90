!==================================================================!
! <> HYCOM dataset interpolation to NCOM grid (ver.2.0)            !
!------------------------------------------------------------------!
! Description :                                                    !
! 1) This program is used for interpolating HYCOM coarse data to   !
!    NCOM grid.                                                    !
!------------------------------------------------------------------!
!                               Reference : NCOMv4.0 source code   !
!                              Programmar : S.H.Liu                !
!                                 Update  : 2015-11-19             !
!==================================================================!
        program ncom2init
        use netcdf
        implicit none
        integer i,j,k
        integer nrecl,iunit,irec,ios

        !! NCOM Horizontal & Vertical Grid file
        character(len=9) Hgrid_file
        character(len=9) Vgrid_file
        integer ndx,ndy,ndz
        parameter (ndx=409,ndy=409,ndz=50)
        real elon(ndx,ndy),alat(ndx,ndy) 
        real dx(ndx,ndy),dy(ndx,ndy)
        real ncom_depth(ndx,ndy),ang(ndx,ndy)

        !! NCOM I/O Array
        character(len=20) oinit_file
        real nlon(ndx),nlat(ndy),ndep(ndz)
        real ne(ndx,ndy),nu(ndx,ndy,ndz),nv(ndx,ndy,ndz)
        real nt(ndx,ndy,ndz),ns(ndx,ndy,ndz)

        !! HYCOM I/O Array
        character(len=35) hycomfile
        integer hdx,hdy,hdz
        parameter (hdx=626,hdy=876,hdz=40)
        real hlon(hdx),hlat(hdy),hdep(hdz)
        real he(hdx,hdy),hu(hdx,hdy,hdz),hv(hdx,hdy,hdz)
        real ht(hdx,hdy,hdz),hs(hdx,hdy,hdz)

        !! Input Segement
        integer icount
        real ework(ndx*ndy)
        real uwork(ndx*ndy)
        real vwork(ndx*ndy)
        real twork(ndx*ndy)
        real swork(ndx*ndy)

!! (1) Read the NCOM horizontal grid data (Get lon,lat,depth)
!!-----------------------------------------------------------------
        Hgrid_file='ohgrd_1.A'
        Vgrid_file='ovgrd_1.D'
        call read_hgrid(trim(Hgrid_file),ndx,ndy,elon,alat,&
                             dx,dy,ncom_depth,ang)
        call read_vgrid(trim(Vgrid_file),50,35,ndep)
        ndep=-ndep  !! negetive to positive

        !!! Convert to 1-d array
        do i=1,ndx
        nlon(i)=elon(i,1)
        end do

        do j=1,ndy
        nlat(j)=alat(1,j)
        end do

!! (2) Read HYCOM dataset and Regrid to NCOM grid
!!-----------------------------------------------------------------
        hycomfile='hycom_glb_regp05_2015111400_t000.nc'
        call read_hycom(trim(hycomfile),hlon,hlat,hdep, &
                        he,hu,hv,ht,hs)

!! (3) Interpolation HYCOM grid to NCOM grid
!!-----------------------------------------------------------------
        call IntpHYCOM(nlon,nlat,ndep,ne,nu,nv,nt,ns,& 
                       ncom_depth,&
                       hlon,hlat,hdep,he,hu,hv,ht,hs)

!! (4) Put vars into NCOM Binary Format file (oinit_1.A.XXXX)
!!-----------------------------------------------------------------
        iunit = 99
        irec  = 01
        nrecl = ndx * ndy
        oinit_file='oinit_1.A'
        open(iunit,file=trim(oinit_file),form='UNFORMATTED',&
                   status='unknown',access='DIRECT',recl=4,&
                   action='WRITE')

        write(*,*) "Test the speed"
!        do j = 1, ndy
!        do i = 1, ndx
!          ework(icount) = ne(i,j)
          call wtNCOM(ne, nrecl, iunit, irec, ios)
!        end do
!        end do
        stop

        !! Total 197 records in oinit file
        do irec = 1, 197
        icount = 0

          if (  irec .eq.  1 ) &
                print *, "Write elevation to NCOM oinit"
          if (( irec .ge.  2 ) .and. ( irec .le. 50 )) &
                print *, "Write U component current to NCOM oinit"
          if (( irec .ge. 51 ) .and. ( irec .le. 99 )) &
                print *, "Write V component current to NCOM oinit"
          if (( irec .ge.100 ) .and. ( irec .le.148 )) &
                print *, "Write Temperature to NCOM oinit"
          if (( irec .ge.149 ) .and. ( irec .le.197 )) &
                print *, "Write Salinity to NCOM oinit"

          do j = 1, ndy
          do i = 1, ndx
          icount = icount + 1
          
          !! Write elevation to NCOM oinit
          if (  irec .eq. 1 ) then
          ework(icount) = ne(i,j)
          call wtNCOM(ework, nrecl, iunit, irec, ios)
          end if

          !! Write U component current to NCOM oinit
          if (( irec .ge. 2 ).and. ( irec .le. 50 )) then
          uwork(icount) = nu(i,j,irec)
          call wtNCOM(uwork, nrecl, iunit, irec, ios)
          end if

          !! Write V component current to NCOM oinit
          if (( irec .ge. 51 ).and. ( irec .le. 99 )) then
          vwork(icount) = nv(i,j,irec)
          call wtNCOM(vwork, nrecl, iunit, irec, ios)
          end if

          !! Write Temperature to NCOM oinit
          if (( irec .ge.100 ).and. ( irec .le.148 )) then
          twork(icount) = nt(i,j,irec)
          call wtNCOM(twork, nrecl, iunit, irec, ios)
          end if

          !! Write Salinity to NCOM oinit
          if (( irec .ge.149 ).and. ( irec .le.197 )) then
          swork(icount) = ns(i,j,irec)
          call wtNCOM(swork, nrecl, iunit, irec, ios)
          end if

          end do
          end do

        end do

!! (5) Output to netCDF for test (If needs,turn ncout = 1)
!!-----------------------------------------------------------------

        print *,"NCOM Initial Field Finished"
        end program



!-----------------------------------------------------------------!
! (1) Read NCOM Horizontal & Vertical Grid                        !
!-----------------------------------------------------------------!
! Description:                                                    !
! (Horizontal)                                                    !
! 1.infile      : input filename                                  !
! 2.n           : Horizontal X Grid numbers                       !
! 3.m           : Horizontal Y Grid numbers                       !
! 4.elon(n,m)   : Longitude                                       !
! 5.alat(n,m)   : Latitude                                        !
! 6.dx(n,m)     : Distance Between Two X points                   !
! 7.dy(n,m)     : Distance Between Two Y points                   !
! 8.h(n,m)      : Smooth Depth (Not real Depth)                   !
! 9.ang(n,m)    : Rotation Angle at the location.                 !
!                 (For Huge scalar only)                          !
!-----------------------------------------------------------------!
! (Vertical)                                                      !
! 1.infile      : input filename                                  !
! 2.l           : Total layers + 1                                !
! 3.ls          : sigma layers                                    !
! 4.zw          : Sigma-Z depth                                   !
!-----------------------------------------------------------------!
        subroutine read_hgrid(infile, n,m, elon,alat,dx,dy,h,ang)
        implicit none
        character*(*) infile
        integer n,m
        real    elon(n,m),alat(n,m),dx(n,m),dy(n,m),h(n,m),ang(n,m)
        integer iunit,nrecl

        !! Open file
        iunit=199
        nrecl=4*n*m
        open(iunit,file=infile,form='unformatted', &
             access='direct',recl=nrecl,status='old')
  
        !! Read field
        read(iunit,rec=1) elon
        read(iunit,rec=2) alat
        read(iunit,rec=3) dx
        read(iunit,rec=4) dy
        read(iunit,rec=5) h
        read(iunit,rec=6) ang

        !! Close file
        close(iunit)
        return
        end

        subroutine read_vgrid(infile,l,ls,zw)
        implicit none
        character*(*) infile
        integer l,ls
        real    zw(l)
        integer iunit,ierr,k
        real*4  al,als,azw(l)

        iunit=199
        open(iunit,file=infile,form='unformatted',status='old')
        read(iunit) al,als,azw
        close(iunit)

        ierr=0
        if (nint(al) .ne. l) then
        write(6,'(/a)') 'Error in read_vgrid: value of l does not agree'
        ierr=1
        endif
        if (nint(als) .ne. ls) then
        write(6,'(/a)')'Error in read_vgrid: value of ls does not agree'
        ierr=1
        endif
        if (ierr .ne. 0) then
        stop 'Error in read_vgrid:  value of l and/or ls does not agree'
        endif

        do k=1,l
        zw(k)=azw(k)
        enddo

        return
        end

!-----------------------------------------------------------------!
! (2) Read HYCOM dataset                                          !
!-----------------------------------------------------------------!
        subroutine read_hycom(infile,lon,lat,dep,e,u,v,t,s) 
        use netcdf
        implicit none
        character (len=*) infile
        integer i,j,k
        integer status,ncid
        integer err
        parameter (err=0)

        !! Dimension ID
        integer LonDimId, LatDimId, DepDimId
        integer  E_DimId,  U_DimId, V_DimId
        integer  T_DimId,  S_DimId

        !! Dimension Var ID
        integer LonVarId, LatVarId, DepVarId
        integer  E_VarId,  U_VarId, V_VarId
        integer  T_VarId,  S_VarId

        ! HYCOM Dimensions 
        integer hdx,hdy,hdz 
        parameter (hdx=626,hdy=876,hdz=40)

        !! Get HYCOM Variables
        real lon(hdx),lat(hdy),dep(hdz)
        real e(hdx,hdy)
        real u(hdx,hdy,hdz),v(hdx,hdy,hdz)
        real t(hdx,hdy,hdz),s(hdx,hdy,hdz)

        !! Test open file
        STATUS=NF90_OPEN(trim(infile),NF90_NOWRITE,NCID)
        IF(STATUS.ne.0) THEN
        write(*,*) "HYCOM data not found !!",ncid,trim(infile)
        write(*,*) NF90_STRERROR(STATUS)
        STOP
        EndIF

        !! Get Dimension ID
        STATUS=NF90_INQ_DIMID(NCID,'lon',LonDimId)
        STATUS=NF90_INQ_DIMID(NCID,'lat',LatDimId)
        STATUS=NF90_INQ_DIMID(NCID,'depth',DepDimId)
        STATUS=NF90_INQ_DIMID(NCID,'surf_el',E_DimId)
        STATUS=NF90_INQ_DIMID(NCID,'water_u',U_DimId)
        STATUS=NF90_INQ_DIMID(NCID,'water_v',V_DimId)
        STATUS=NF90_INQ_DIMID(NCID,'water_temp',T_DimId)
        STATUS=NF90_INQ_DIMID(NCID,'salinity',S_DimId)

        !! Get Var ID
        STATUS=NF90_INQ_VARID(NCID,'lon',LonVarId)
        STATUS=NF90_INQ_VARID(NCID,'lat',LatVarId)
        STATUS=NF90_INQ_VARID(NCID,'depth',DepVarId)
        STATUS=NF90_INQ_VARID(NCID,'surf_el',E_VarId)
        STATUS=NF90_INQ_VARID(NCID,'water_u',U_VarId)
        STATUS=NF90_INQ_VARID(NCID,'water_v',V_VarId)
        STATUS=NF90_INQ_VARID(NCID,'water_temp',T_VarId)
        STATUS=NF90_INQ_VARID(NCID,'salinity',S_VarId)

        !! Get HYCOM Variables
        STATUS=NF90_GET_VAR(NCID,LonVarId,lon)
        STATUS=NF90_GET_VAR(NCID,LatVarId,lat)
        STATUS=NF90_GET_VAR(NCID,DepVarId,dep)
        STATUS=NF90_GET_VAR(NCID,E_VarId,e)
        STATUS=NF90_GET_VAR(NCID,U_VarId,u)
        STATUS=NF90_GET_VAR(NCID,V_VarId,v)
        STATUS=NF90_GET_VAR(NCID,T_VarId,t)
        STATUS=NF90_GET_VAR(NCID,S_VarId,s)
        STATUS=NF90_CLOSE(NCID)

        !! Shift the scalar in HYCOM netcdf
        do k=1,hdz
        do j=1,hdy
        do i=1,hdx

        if(k.eq.1) then 
        if(e(i,j).ne.-30000) e(i,j)=e(i,j)*0.001
        end if
        
        if(u(i,j,k).ne.-30000) u(i,j,k)=u(i,j,k)*0.001
        if(v(i,j,k).ne.-30000) v(i,j,k)=v(i,j,k)*0.001
        if(t(i,j,k).ne.-30000) t(i,j,k)=t(i,j,k)*0.001+20
        if(s(i,j,k).ne.-30000) s(i,j,k)=s(i,j,k)*0.001+20

        end do
        end do
        end do

        return
        
        end subroutine

!-----------------------------------------------------------------!
! (3) Interpolate HYCOM coarse grid to NCOM                       !
!-----------------------------------------------------------------!
        subroutine IntpHYCOM(nlon,nlat,ndep,ne,nu,nv,nt,ns,&
                             ncom_depth,&
                             hlon,hlat,hdep,he,hu,hv,ht,hs)
        implicit none
        integer i,j,k
        integer ndx,ndy,ndz
        integer hdx,hdy,hdz
        parameter (ndx=409,ndy=409,ndz=49)
        parameter (hdx=626,hdy=876,hdz=40)

        real ncom_depth(ndx,ndy)
        real nlon(ndx),nlat(ndy)
        real ndep(ndz),ne(ndx,ndy)
        real nu(ndx,ndy,ndz),nv(ndx,ndy,ndz)
        real nt(ndx,ndy,ndz),ns(ndx,ndy,ndz)
        real hlon(hdx),hlat(hdy)
        real hdep(hdz),he(hdx,hdy)
        real hu(hdx,hdy,hdz),hv(hdx,hdy,hdz)
        real ht(hdx,hdy,hdz),hs(hdx,hdy,hdz)

        !! Interpolation Arrays
        real ilon(138),ilat(138),ilayer(40)
        real iu(138,138,40),iv(138,138,40)
        real itemp(138,138,40),isal(138,138,40)
        real iel(138,138)

        !! Horizontal Interpolation (Pass Arrays)
        real ou(ndx,ndy,40),ov(ndx,ndy,40)
        real otemp(ndx,ndy,40),osal(ndx,ndy,40)
        real oel(ndx,ndy)

        !! Resize to approach the NCOM grid
        do k=1,40
        do j=1,138
        do i=1,138
        ilon(i)=hlon(213+i)
        ilat(j)=hlat(325+j)
        itemp(i,j,k)=ht(213+i,325+j,k)
        isal(i,j,k)=hs(213+i,325+j,k)
        iu(i,j,k)=hu(213+i,325+j,k)
        iv(i,j,k)=hv(213+i,325+j,k)
        iel(i,j)=he(213+i,325+j)

        if(ilon(i).eq.-30000) ilon(i)=0
        if(ilat(j).eq.-30000) ilat(j)=0
        if(itemp(i,j,k).eq.-30000) itemp(i,j,k)=0
        if(isal(i,j,k).eq.-30000) isal(i,j,k)=0
        if(iu(i,j,k).eq.-30000) iu(i,j,k)=0
        if(iv(i,j,k).eq.-30000) iv(i,j,k)=0
        if(iel(i,j).eq.-30000) iel(i,j)=0
        end do
        end do
        end do
        
        !! Filled the blank data
        write(*,*) "Start NCOM Data Interpolation:"
        write(*,*) "----------------------------------------"
 
        do k=1,hdz
          write(*,*) "Horizontal Grid Data Interpolation: ",k
          call fillm(138,138,itemp(:,:,k),0.0)
          call fillm(138,138,isal(:,:,k),0.0)
          call fillm(138,138,iu(:,:,k),0.0)
          call fillm(138,138,iv(:,:,k),0.0)
          if(k.eq.1) call fillm(138,138,iel(:,:),0.0)

          !! Horizontal Interpolation
          call SpAk2D(itemp(:,:,k),ilon,ilat,138,138,&
                      otemp(:,:,k),nlon,nlat,409,409)
          call SpAk2D(isal(:,:,k),ilon,ilat,138,138,&
                      osal(:,:,k),nlon,nlat,409,409)
          call SpAk2D(iu(:,:,k),ilon,ilat,138,138,&
                      ou(:,:,k),nlon,nlat,409,409)
          call SpAk2D(iv(:,:,k),ilon,ilat,138,138,&
                      ov(:,:,k),nlon,nlat,409,409)
          if(k.eq.1)  call SpAk2D(iel(:,:),ilon,ilat,138,138,&
                      oel(:,:),nlon,nlat,409,409)
        end do

        !! Vertical Interpolation
        write(*,*)
        write(*,*) "Vertical Grid Data Interpolation : "
        do j=1,ndy
        do i=1,ndx
        call SPAK1D(hdep(:),otemp(i,j,:),40,&
                    ndep(:),nt(i,j,:),49)
        call SPAK1D(hdep(:),osal(i,j,:),40,&
                    ndep(:),ns(i,j,:),49)
        call SPAK1D(hdep(:),ou(i,j,:),40,&
                    ndep(:),nu(i,j,:),49)
        call SPAK1D(hdep(:),ov(i,j,:),40,&
                    ndep(:),nv(i,j,:),49)
        end do
        end do

        ne=oel

        !! Mask the Land & Sea grid (By NCOM depth)
        do k=1,49
        write(*,*) "Now masking the vertical grid: ",k
        do j=1,409
        do i=1,409
        if((ncom_depth(i,j).gt.ndep(k)).or.&
           (ncom_depth(i,j).eq.0)) then
           nu(i,j,k)=0
           nv(i,j,k)=0
           nt(i,j,k)=0
           ns(i,j,k)=0

        if(k.eq.1)then
          if((ncom_depth(i,j).gt.ndep(k)).or.&
             (ncom_depth(i,j).eq.0)) then
             ne(i,j)=0
          end if
        end if

        end if

        end do
        end do
        end do

        end subroutine IntpHYCOM




!-----------------------------------------------------------------!
!  Write dimensions to oinit_1.B                                  !
!-----------------------------------------------------------------!
        subroutine wrt_init_B(nt,mt,l,ls,ncomdate)
        implicit none
        integer ncomdate
        integer nt,mt,l,ls
        character (len=20) filename

        ! Write filename
        write(filename,1) "oinit_1.B.",ncomdate,"00"

        ! Open the file and write the parameter
        open(10,file=filename) 
        write(10,2) nt,mt,l,ls
        
        ! Close the file
        close(10)

1       FORMAT(A10,I8,A2)
2       FORMAT(4(I6))
        end subroutine
       
!====================================================================!
! <> 填補網格空洞用Sapiro Filter副程式 <>                            !
!====================================================================!
        Subroutine fillm(nx,ny,F,ERR)
        implicit none
        integer nx,ny,i,j,i1,il,j1,jl,ii,jj,iw,jw
        integer nbadb,count_bad
        real ERR
        real F(nx,ny), tmp(nx,ny),wt,hst

        real W(-1:1,-1:1)
        data W / 1.,2.,1.,&
                 2.,4.,2.,&
                 1.,2.,1. /
        nbadb = nx*ny
1       count_bad = 0
        do j = 1,ny
          do i = 1,nx
            tmp(i,j) = F(i,j)
          end do
        end do

        do j = 1,ny
          j1 = max0(j-1,1)
          jl = min0(j+1,ny)
        do i = 1, nx
          i1 = max0(i-1,1)
          il = min0(i+1,nx)

          if ( tmp(i,j) .eq. ERR ) then
            count_bad = count_bad + 1
            WT = 0
            HST = 0

            do jj = j1,jl
               jw = jj-j
               do ii = i1,il
                  iw = ii-i

                 if ( tmp(ii,jj) .ne. ERR ) then
                   HST = HST + W(iw,jw)*tmp(ii,jj)
                   WT = WT + W(iw,jw)
                 end if

               end do
            end do

            if (WT.gt.0) F(i,j) = HST/WT
            end if

        end do
        end do
        if ( count_bad .eq. 0 .or. count_bad .eq. nbadb ) return
        nbadb = count_bad
        goto 1

        end
!===================================================================!
! <> 內插用副程式 Akima Spline Interpolations
!--------------------------------------------------------------------
!    Description: 
!    不同於Cubic Spline Interpolation 會有邊界值問題,在缺少資料或轉
!    折處會有較大的誤差,Akima Spline Interpolation則解決了在此種問題
!    ,並將誤差降至最低.
!--------------------------------------------------------------------
!    References:
!    H. Akima, "A New Method of Interpolation and Smooth Curve
!              Fitting Based on local procedures", J. Ass. for
!              Computing Machinery, Vol. 17, No. 4, Oct. 1970,
!              pp. 589-602.
!--------------------------------------------------------------------
!    Copied from NCOM source code (ncom_setup_spln.F)
!                            Oringin Programmar: D.S.Ko   1997-11-05 
!                        Modified to f90 format: S.H.Liu  2015-01-15
!===================================================================!
        Subroutine  SpAk2D(F,X,Y,NX,NY, FI,XI,YI,NXI,NYI)
        !## input
        real      F(NX,NY),X(NX),Y(NY)
        real      XI(NXI),YI(NYI)
        !## Output
        real      FI(NXI,NYI)
        !## Temp
        parameter (Nmax=10000)    ! max dimension in x/y
        real      FU(Nmax)
        real      FIU(Nmax)
        real      FX(Nmax,Nmax), FY(Nmax,Nmax)

        !## Interpolate F(X,Y) to F(Xi,Y)
        do j = 1, NY
        do i = 1, NX
        FU(i) = F(i,j)
        end do
        call SpAk1d (X,FU,NX,XI,FIU,NXI)
        do i = 1, NXI
        FX(i,j) = FIU(i)
        end do
        end do

        !## Interpolate F(Xi,Y) to F(Xi,Yi)
        do i = 1, NXI
        do j = 1, NY
        FU(j) = FX(i,j)
        end do
        call spak1d (Y,FU,NY,YI,FIU,NYI)
        do j = 1, NYI
        FX(i,j) = FIU(j)
        end do
        end do

        !## Interpolate F(X,Y) to F(X,Yi)
        do i = 1, NX
        do j = 1, NY
        FU(j) = F(i,j)
        end do
        call spak1d (Y,FU,NY,YI,FIU,NYI)
        do j = 1, NYI
        FY(i,j) = FIU(j)
        end do
        end do

        !## Interpolate F(X,Yi) to F(Xi,Yi)
        do j = 1, NYI
        do i = 1, NX
        FU(i) = FY(i,j)
        end do
        call spak1d (X,FU,NX,XI,FIU,NXI)
        do i = 1, NXI
        FI(i,j) = .5*(FX(i,j)+FIU(i))
        end do
        end do

        return
        end

! -------------------------------------------------------------------
        Subroutine SpAk1D(X,Y,N,XI,YI,NI)
        real X(*),Y(*),XI(*),YI(*)
        parameter (Nmax=10000) ! modify the no. as needed
        real Coef(4,Nmax)
        real Temp1(Nmax+4), Temp2(Nmax+4)
        call SplAkm (X,Y,N,Coef,Temp1,Temp2)
        call Splder (XI,YI,NI, N,X,Coef)
        return
        end

!--------------------------------------------------------------------
        Subroutine  SplAkm (x,y,n,coef,s,ss)
        real    x(n), y(n)
        real    coef(4,n)
        real    s(-1:n+1), ss(n)
        real*8  a,b,ab, dx
        !## Compute slope for each interval of points
        do i = 1, n-1
        s(i) = (y(i+1)-y(i))/(x(i+1)-x(i))
        end do

        !## Compute 2 extra points at each end
        s(-1 ) = 3.*s(1  ) - 2.*s(2  )
        s( 0 ) = 2.*s(1  ) -    s(2  )
        s(n  ) = 2.*s(n-1) -    s(n-2)
        s(n+1) = 3.*s(n-1) - 2.*s(n-2)

        !## Compute Akima slope for every point (each slope needs 5
        !points)
        do i = 1, n
        a = abs(s(i-1) - s(i-2))
        b = abs(s(i+1) - s(i  ))
        ab = a+b
        if (ab .eq. 0.) then
        ss(i) = (s(i-1) + s(i)) * .5
        else
        ss(i) = (b*s(i-1) + a*s(i)) / ab
        endif
        end do

        !## Calculate coefficients: const 1,1,2,6 have been factored out
        do i = 1, n-1
        dx = x(i+1) - x(i)
        coef(1,i) = y(i)
        coef(2,i) = ss(i)
        coef(3,i) = (3.*s(i)-2.*ss(i)-ss(i+1))/dx
        coef(4,i) = (ss(i)+ss(i+1)-2.*s(i))/(dx*dx)
        end do

        !## fill in with linear slope at end (but they should not be
        !used)
!       coef(1,n) = y(n)
!       coef(2,n) = ss(n)
!       coef(3,n) = 0.
!       coef(4,n) = 0.

        return
        end

! ----------------------------------------------------------------------
        Subroutine Splder (xi,yi,ni, n,x,coef)
        real  xi(*),yi(*)
        real  x(*),coef(4,*)
        i1 = 1
        il = n-1

        do ii = 1, ni

        if (xi(ii).lt.x(i1+1)) then
          i = i1
        else if (xi(ii).ge.x(il)) then
          i = il
        else
          i = i1+1
        !## Bracket xi with breakpoints x(i)
          do while (xi(ii).ge.x(i))
            i = i + 1
          end do
          i = i - 1
        end if

        dx = xi(ii) - x(i)

        yi(ii) = ((coef(4,i)*dx+coef(3,i))*dx+coef(2,i))*dx+coef(1,i)

        end do

        return
        end

!----------------------------------------------------------------!
! Write single record to NCOM Binary                             !
!----------------------------------------------------------------!
! Description:                                                   !
! 1)  direct access write a single record.                       !
! 2)  expressed as a subroutine because i/o with                 !
!     implied do loops can be slow on some machines.             !
!----------------------------------------------------------------!
!                              AUTHOR : Alan J. Wallcraft, NRL   !
!                       Re-Programmar : S.H.Liu                  !
!                                Date : 2015-11-17               !
!----------------------------------------------------------------!
      subroutine wtNCOM(a,n, iunit,irec,ios)
      implicit none
      integer n,iunit,irec,ios
      real*4  a(n)
      write(unit=iunit, rec=irec, iostat=ios) a
      return
      end


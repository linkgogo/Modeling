        program CTD_Statistics
        implicit none
!#==============================================================#!
! <> 國科會 CTD 統計分析 及 數據擷取程式 (ver.2.1)               !
!----------------------------------------------------------------!
! 說明:                                                          !
! 1.本程式將依照(經緯度+月份)擷取CTD資料,並計算CTD歷年統計月平均 !
!   資料,用來提供驗證用的分析資料.                               !
! 2.注意ship-CTD施放時,表層水會因為上下擺盪或氣泡導致量測失誤,請 !
!   斟酌使用,本人建議使用 5m 以下的資料.                         !
! 3.本資料處理垂直分層依序為:                                    !
!   0,   2,   4,   6,   8,  10,  12,  15,  20,  25,   30,        !
!  35,  40,  45,  50,  60,  70,  80,  90, 100, 125,  150,        !
! 200, 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000         !
!   若要自行定義,請自行調整layer層的各層深度,並重新編譯程式.     !
!----------------------------------------------------------------!
! 檔案格式:                                                      !
! 經度  緯度  年  月  日  時  分  深度  溫度  鹽度  密度(dbar)   !
!  (degree)                       (m)  (degC) (PSU) (kg/m3)      !
!----------------------------------------------------------------!
! 更新歷程:                                                      !
! v1.0 (2014-08-04)                                              !
! copied original program from CTD_SORT_ALL.f90.                 !
!                                                                !
! v1.1 (2014-09-05)                                              !
! set mode options:                                              !
! mode=1,                                                        !
! select the monthly data of the CTD raw dataset, then average   !
! the total data into 1 file format.                             !
! mode=2,                                                        !
! select the monthly data of the CTD raw dataset, then list all  !
! data into different files. filename rule: f00000-f99999        !
!                                                                !
! v1.2 (2014-10-01)                                              !
! set more mode options:                                         !
! mode=3,                                                        !
! adds the format for yearly output. filename count: f000-f999   !
! mode=4,                                                        !
! adds the format for YYYYMMDD output. filename ctd_out_YYYYMMDD !
!                                                                !
! v1.3 (2014-11-05)                                              !
! adds "mode" option check, if you type number not between 1~4,  !
! the program will return to the question and ask your input     !
! again.                                                         !
!                                                                !
! v1.4 (2014-12-05)                                              !
! 1.adds sound velocity analysis format in mode 1.               !
!                                                                !
! v1.5 (2015-03-30)                                              !
! 1.Test gfortran compiler                                       !
! 2.Test PGI Fortran compiler                                    !
! 3.Adds instrinsic "SYSTEM" to creat ./allout directory.        !
! 4.Adds "Death Time" in this program.                           !
!                                                                !
! v2.0 (2015-09-13)                                              !
! 1.Clean up the program. (Refresh)                              !
! 2.Adds 0.5 degree cell NETCDF output                           !
! 3.Adds standard-deviation for data analysis.                   !
! 4.Adds/Fixed new Sound Speed solver.                           !
! 5.Remove unneccessary mode(only remained mode 1).              !
! 6.Output to the program direction.                             !
!   Regard the allout direction.                                 !
! 7.Renew Filename to:                                           !
!   (1) ASCII   -> CTD-Mean.txt                                  !
!   (2) netCDF  -> CTDStd.nc                                     !
! 8.To keep the data reliable, I won't interpolate CTD data.     !
!   However, it may makes the data lack of data in somewhere.    !
!                                                                !
! v2.1 (2015-11-06)                                              !
! 1.Adds the center algorithm in point calculation.              !
!----------------------------------------------------------------!
!                                            Create:  S.H.Liu    !
!                                            Update:  2015-11-08 !
!#==============================================================#!
        !! Iteral/Count numbers
        integer i,j,k,l,m,n
        integer, parameter :: layer=33  !! depth layers
        integer icount,tcount(layer)
        integer imax,idepth(layer)
        integer itemp(layer),isal(layer),idensity(layer)

        !! Dataset Array
        integer, parameter :: ln=10422479
        real lon(ln),lat(ln)
        integer year(ln),  month(ln), day(ln)
        integer hour(ln), minute(ln), sec(ln)
        integer depth(ln)
        real temp(ln),sal(ln),density(ln)

        !! ASCII Files
        integer mode
        integer :: ios=0

        !! Dimension Defined
        real ilon1,ilon2
        real ilat1,ilat2
        integer iyear1,iyear2
        integer imonth1,imonth2

        !! Analysis Array
        real ftemp(layer), fsal(layer)
        real fdensity(layer),fsv(layer)
        real,allocatable :: t2lon(:),t2lat(:),t2year(:),t2month(:)
        real,allocatable :: t2day(:),t2hour(:),t2minute(:),t2sec(:)
        real,allocatable :: t2depth(:),t2temp(:),t2sal(:),t2density(:)
        real,allocatable :: ttemp(:,:),tsal(:,:),tdensity(:,:)

        !! Others
        real tlon,tlat  !! Lat-Lon in the Center of Grid 
        real proc       !! Processing --% complete
        integer depth_index(layer) / 0,2,4,6,8,10,12,15,20,25,30,35,&
                40,45,50,60,70,80,90,100,125,150,200,250,300,350,&
                400,500,600,700,800,900,1000/
        real T,S,Dep !! <- SoundSpeed Calculate

        !! Looking for TWN-CTD 1985-2010 dataset
        open(10,file='1985_2010_ctd.txt',status='old',&
                action='read',iostat=ios)
        if(ios.ne.0)then
        print *, "Error: Can't find 1985_2010_ctd.txt"
        print *, "Please make sure the file in the same direction. "
        print *, "And re-execute the program. "
        stop
        endif

        call system("clear")
        print *, "#===================================================#"
        print *, "# <> TWN-CTD Analysis program :                     #"
        print *, "#---------------------------------------------------#"
        print *, "#                  ---------------                  #"
        print *, "#                  >>> Warning <<<                  #"
        print *, "#                  ---------------                  #"
        print *, "# The old data may not include the GPS accurate     #"
        print *, "# position. Be avoid if the data range select is    #"
        print *, "# in these years.                                   #"
        print *, "#---------------------------------------------------#"
        print *, "#                 Program Manager : S.H.Liu         #"
        print *, "#                 Program Updated : 2015-11-08      #"
        print *, "#===================================================#"
        print *, "Press <Enter> to continue."
        read(*,*)
       
        call system("clear")
        print *, "====================================================="
        print *, " <>>>> Select the processing : ( Enter 1 ~ 3 ) <<<<> "
        print *, "-----------------------------------------------------"
        print *, "  mode 1 :                                           "
        print *, "  Average in User defined range. Give 4 points       "
        print *, "  longitude and latitude to get the Std. and Mean    "
        print *, "  analyst in ASCII FORMAT file.                      "
        print *, "-----------------------------------------------------"
        print *, "  mode 2 :                                           "
        print *, "  The same as mode 1. Only need input one point,     "
        print *, "  then calculating the square 1 degree average.      "
        print *, "-----------------------------------------------------"
        print *, "  mode 3 :                                           "
        print *, "  netCDF file output. It cost expensive computing.   "
        print *, "  All analysis data will be calculated into netcdf,  "
        print *, "  so you only need to run once.                      "
        print *, "====================================================="
9001    read(*,*) mode
        if((mode.lt.1).or.(mode.gt.3))then
        call system("clear")
        print *,"Error: No such processing. ( mode: 1 ~ 3 )           "
        print *,"Please <ENTER> 1 ~ 3 to continue the process         "
        goto 9001
        endif

!-----------------------------------------------------------------!
! <> mode 1: Average + Standard Deviation                         !
!-----------------------------------------------------------------!
        select case(mode)
        case(1)
        call system("clear")

9101    print *,'------------------------------------------------'
        print *,' (1)       (2)                                  '
        print *,'  +--------+                                    '
        print *,'  | o     o|      1. Input Area (1) - (4)       '
        print *,'  |   + o  |  --> 2. Average CTD observations   '
        print *,'  |o     o |      3. Put data in intersection   '
        print *,'  +--------+                                    '
        print *,' (3)       (4)                                  '
        print *,'------------------------------------------------'
        print *,''
        print *,' <> Enter the longitude Range (Start, End)      '
        print *,'    Range: 100E ~ 150E'
        print *,"------------------------------------------------"
        read(*,*) ilon1,ilon2
        if (((ilon1.lt.100).or.(ilon1.gt.150)) .or. &
            ((ilon2.lt.100).or.(ilon2.gt.150)) .or. &
            ((ilon1.gt.ilon2))) then
        call system("clear")
        print *, " Error: Out of the range. ( 100 ~ 150 ).          "
        print *, " Please Enter longitude(Start, End) again!        "
        goto 9101
        stop
        end if

        call system("clear")
        print *,' <> Enter the latitude Range (Start, End)     '
        print *,'    Range: 0N ~ 40N'
        print *,"------------------------------------------------"
9102    read(*,*) ilat1,ilat2
        if (((ilat1.lt.0).or.(ilat1.gt.40)) .or. &
            ((ilat2.lt.0).or.(ilat2.gt.40)) .or. &
            ((ilat1.gt.ilat2))) then
        call system("clear")
        print *, " Error: Out of the range. ( 0 ~ 40 ).             "
        print *, " Please Enter latitude(Start, End) again!         "
        goto 9102
        stop
        end if

        call system("clear")
        print *, ' <> Enter the Year Range (Start, End)         '
        print *, '    Range: 1985 - 2010'
        print *, "------------------------------------------------"
        print *, ' If you want only 1 year, you can enter the same'
        print *, ' Year to start the process. (ex. 2000 2002 )    '
9103    read(*,*) iyear1,iyear2
        if (((iyear1.lt.1985).or.(iyear1.gt.2010)) .or. &
            ((iyear2.lt.1985).or.(iyear2.gt.2010)) .or. &
            ((iyear1.gt.iyear2))) then
        call system("clear")
        print *, " Error: Out of the range. ( 1985 ~ 2010 ).      "
        print *, " Please Enter Year(Start, End) again!           "
        goto 9103
        end if

        call system("clear")
        print *, ' <> Enter the Month Range (Start, End)        '
        print *, "    Range: 1 - 12"
        print *, "-------------------------------------------------"
        print *, ' If you want only 1 month, you can enter the same'
        print *, ' Month to start the process. (ex. 7 7)           '
9104    read(*,*) imonth1,imonth2
        if (((imonth1.lt.1).or.(imonth1.gt.12)) .or. &
            ((imonth2.lt.1).or.(imonth2.gt.12)) .or. &
            ((imonth1.gt.imonth2))) then
        call system("clear")
        print *, " Error: Out of the range. ( 1 ~ 12 ).           "
        print *, " Please Enter Month(Start, End) again!          "
        goto 9104
        end if

        call system("clear")
        print *,"+---------------+","=================================="
        print *,"|               |","        < User Define >           "
        print *,"|               |","----------------------------------"
        print *,"|               |"," Year  : ",iyear1,"~",iyear2
        print *,"|   (Average)   |"," Month : ",imonth1,"~",imonth2
        print *,"|               |"," Lat   : ",ilat1,"~",ilat2 
        print *,"|               |"," Lon   : ",ilon1,"~",ilon2
        print *,"|               |"
        print *,"+---------------+","=================================="
        print *,"Target Longitude : ",(ilon1+ilon2)/2
        print *,"Target Latitude  : ",(ilat1+ilat2)/2
        print *,"+---------------+","=================================="
        write(*,*) "Start Loading dataset...."

        tlon=(ilon1+ilon2)/2
        tlat=(ilat1+ilat2)/2


        !! Initialize the accumulate matrix
        idepth(:)=0

        !! Read the CTD dataset 
        do i=1,ln
        read(10,*) lon(i),lat(i),year(i),month(i),day(i),hour(i),&
                   minute(i),sec(i),depth(i),temp(i),sal(i),&
                   density(i)

        !! Data Select in the User define range.
        if(((lon(i).ge.ilon1).and.(lat(i).ge.ilat1)).and.&
           ((lon(i).le.ilon2).and.(lat(i).le.ilat2)).and.&
           ((month(i).ge.imonth1).and.(month(i).le.imonth2)).and.&
           ((year(i).ge.iyear1).and.(year(i).le.iyear2))) then 

          !! Count times of data in each depth layer
          do l=1,33
            if (depth(i).eq.depth_index(l))then
                idepth(l) = idepth(l) + 1     
            endif
          end do

        endif

        !!! Processing % complete
          do m=1,10
            if(i.eq.999999*m)then
              proc = (real(i)/real(ln)) * 100
              write(*,1001) "Data Processing Complete: ",proc,"%"
1001          FORMAT(A26,F9.4,A3)
            end if
          end do

        end do

        !! Close the TWN-CTD dataset file
        close(10)

        !! Determine the new matrix size
        imax=max(idepth(1),idepth(2),idepth(3),idepth(4),idepth(5),&
               idepth(6),idepth(7),idepth(8),idepth(9),idepth(10),&
               idepth(11),idepth(12),idepth(13),idepth(14),idepth(15),&
               idepth(16),idepth(17),idepth(18),idepth(19),idepth(20),&
               idepth(21),idepth(22),idepth(23),idepth(24),idepth(25),&
               idepth(26),idepth(27),idepth(28),idepth(29),idepth(30),&
               idepth(31),idepth(32),idepth(33))
        allocate(ttemp(imax,layer))
        allocate(tsal(imax,layer))
        allocate(tdensity(imax,layer))

        call system("clear")
        print *, "Start Analysis .........."

        !! Initialize the calculator
        itemp(:)=0       !---------------!
        isal(:)=0        !  Data counts  !
        idensity(:)=0    !---------------!

        ttemp(:,:)=0     !---------------!
        tsal(:,:)=0      !  Data Summary !
        tdensity(:,:)=0  !---------------!

        do i=1,ln
        
        if(((lon(i).ge.ilon1).and.(lat(i).ge.ilat1)).and.&
           ((lon(i).le.ilon2).and.(lat(i).le.ilat2)).and.&
           ((month(i).ge.imonth1).and.(month(i).le.imonth2)).and.&
           ((year(i).ge.iyear1).and.(year(i).le.iyear2))) then
          
          do l=1,33
            if (depth(i).eq.depth_index(l))then
              itemp(l)=itemp(l)+1
              isal(l)=isal(l)+1
              idensity(l)=idensity(l)+1
              ttemp(itemp(l),l)=temp(i)
              tsal(isal(l),l)=sal(i)
              tdensity(idensity(l),l)=density(i)
            end if
          end do

        end if

        end do

        !! Display on the Screen
        write(*,*) "<<<-------- Temperature Average --------->>>"
        write(*,*) "Depth(m)     count      Temp.AVG(Degree)"
        do j=1,layer
        tcount(j)=0
        do i=1,imax
        if (ttemp(i,j).ne.0) tcount(j)=tcount(j)+1
        end do
        ftemp(j)=sum(ttemp(:,j))/real(tcount(j))
        if (tcount(j).eq.0) ftemp(j)=0
        write(*,1002) depth_index(j),'m',tcount(j),ftemp(j)
        end do

        write(*,*) "<<<<<-------- Salinity Average --------->>>>>"
        write(*,*) "Depth(m)     count     Sal.AVG(PSU)"
        do j=1,layer
        tcount(j)=0
        do i=1,imax
        if (tsal(i,j).ne.0) tcount(j)=tcount(j)+1
        end do
        fsal(j)=sum(tsal(:,j))/real(tcount(j))
        if (tcount(j).eq.0) fsal(j)=0
        write(*,1002) depth_index(j),'m',tcount(j),fsal(j)
        end do
1002    FORMAT(I4,A1,5X,I6,5X,F10.4)

        write(*,*) "<<<<<------ Pressure Average ------->>>>>"
        write(*,*) "Depth(m)     count      dbar.AVG(kg/m^2)"
        do j=1,layer
        icount=0
        do i=1,imax
        if (tdensity(i,j).ne.0) icount=icount+1
        end do
        fdensity(j)=(sum(tdensity(:,j))/real(icount))
        if (icount.eq.0) fdensity(j)=0
        write(*,1002) depth_index(j),'m',icount,fdensity(j)
        end do

        !!! *Lat-Lon is stand for the center of the data cell.
        tlon=(ilon1+ilon2)/2
        tlat=(ilat1+ilat2)/2

        !! Open the output filename and write the result in.
        open(21,file='TWN-CTD_STATISTICS.txt',status='unknown',&
                action='write',iostat=ios)
        
        do j=1,layer
        T=ftemp(j)
        S=fsal(j)
        Dep=depth_index(j)
        fsv(j)=1449.2+(4.6*T)-&
               (0.055*(T**2))+&
               (0.0003*(T**3))+&
               (1.34-0.01*T)*(S-35)+&
               (0.016*Dep)

        ! If the temperature or salinity has no data, turn velocity to 0
        if (( ftemp(j) .eq. 0 ) .or. ( fsal(j) .eq. 0 )) fsv(j)=0
        if ( fdensity(j) .eq. 0 ) fdensity(j) = -1000

        write(21,1003) tlon,tlat,depth_index(j),tcount(j),ftemp(j),&
                      fsal(j),fdensity(j)+1000,fsv(j)
        end do
        close(21)
1003    FORMAT(2(F9.4,1X),2(I4,3X),4(F10.4,1X))

        print *, " * Process Finished. "
        print *, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
        print *, ">>>>   Output to TWN-CTD_STATISTICS.txt    >>>>"  
        print *, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

!---------------------------------------------------------------------!
! <> mode = 2 Output to netCDF ( 0.5 degree / grid )                  !
!---------------------------------------------------------------------!
! < Data Range >                                                      !
! Year  : 1985 - 2010                                                 !
! Month :    1 -   12                                                 !
! Lon.  : 115E -  135E                                                !
! Lat.  :  15N -   35N                                                !
! Depth : Standard Depth for Physical Oceanographic ( ~1000m )        !
!---------------------------------------------------------------------!
!  Variables:                                                         !
!     Dimension(1) : lon                                              !
!     Dimension(2) : lat                                              !
!     Dimension(3) : Year                                             !
!     Dimension(4) : Month                                            !
!     Vars(01)  : AVG. Temperature          (25 years / 12 months)    !
!     Vars(02)  : AVG. Salinity             (25 years / 12 months)    !
!     Vars(03)  : AVG. Pressure(dbar)       (25 years / 12 months)    !
!     Vars(04)  : AVG. Sound Velocity       (25 years / 12 months)    !
!     Vars(05)  : STD. Temperature          (25 years / 12 months)    !
!     Vars(06)  : STD. Salinity             (25 years / 12 months)    !
!     Vars(07)  : STD. Pressure(dbar)       (25 years / 12 months)    !
!     Vars(08)  : STD. Sound Velocity       (25 years / 12 months)    !
!     Vars(09)  : Data Counts               (25 years / 12 months)    !
!     Vars(10)  : Annal AVG. Temperature    ( 12 months )             !
!     Vars(11)  : Annal AVG. Salinity       ( 12 months )             !
!     Vars(12)  : Annal AVG. Pressure       ( 12 months )             !
!     Vars(13)  : Annal AVG. Sound Velocity ( 12 months )             !
!     Vars(14)  : Annal STD. Temperature    ( 12 months )             !
!     Vars(15)  : Annal STD. Salinity       ( 12 months )             !
!     Vars(16)  : Annal STD. Pressure       ( 12 months )             !
!     Vars(17)  : Annal STD. Sound Velocity ( 12 months )             !
!     Vars(18)  : Annal Data Counts         ( 12 months )             !
!     Vars(19)  : Total Data Counts         ( 1 layer )               !
!---------------------------------------------------------------------!

!-----------------------------------------------------------------!
! <> mode 2: Square Average                                       !
!-----------------------------------------------------------------!
        case(2)
        call system("clear")
9201    print *,'------------------------------------------------'
        print *,' (1)       (2)                                  '
        print *,'  +--------+                                    '
        print *,"  | o     o|      1. Input '+' location         "
        print *,'  |   + o  |  --> 2. Average CTD observations   '
        print *,'  |o     o |      3. Output to the file         '
        print *,'  +--------+                                    '
        print *,' (3)       (4)                                  '
        print *,'------------------------------------------------'
        print *,' <> Enter the center point (longitude, latitude)'
        print *,'    Range: 100E ~ 150E                          '
        print *,'             0N ~  40N                          '
        print *,"------------------------------------------------"
        read(*,*) tlon,tlat
        if (((tlon.lt.100).or.(tlon.gt.150)) .or. &
            ((tlat.lt.0).or.(tlat.gt.40))) then
        call system("clear")
        print *, " Error: Out of the range. ( 100 ~ 150 ).          "
        print *, " Please Enter longitude(Start, End) again!        "
        goto 9201
        stop
        end if

        call system("clear")
        print *, ' <> Enter the Year Range (Start, End)         '
        print *, '    Range: 1985 - 2010'
        print *, "------------------------------------------------"
        print *, ' If you want only 1 year, you can enter the same'
        print *, ' Year to start the process. (ex. 2000 2002 )    '
9202    read(*,*) iyear1,iyear2
        if (((iyear1.lt.1985).or.(iyear1.gt.2010)) .or. &
            ((iyear2.lt.1985).or.(iyear2.gt.2010)) .or. &
            ((iyear1.gt.iyear2))) then
        call system("clear")
        print *, " Error: Out of the range. ( 1985 ~ 2010 ).      "
        print *, " Please Enter Year(Start, End) again!           "
        goto 9202
        end if

        call system("clear")
        print *, ' <> Enter the Month Range (Start, End)        '
        print *, "    Range: 1 - 12"
        print *, "-------------------------------------------------"
        print *, ' If you want only 1 month, you can enter the same'
        print *, ' Month to start the process. (ex. 7 7)           '
9203    read(*,*) imonth1,imonth2
        if (((imonth1.lt.1).or.(imonth1.gt.12)) .or. &
            ((imonth2.lt.1).or.(imonth2.gt.12)) .or. &
            ((imonth1.gt.imonth2))) then
        call system("clear")
        print *, " Error: Out of the range. ( 1 ~ 12 ).           "
        print *, " Please Enter Month(Start, End) again!          "
        goto 9203
        end if

        ilat1=tlat - 0.5
        ilat2=tlat + 0.5
        ilon1=tlon - 0.5
        ilon2=tlon + 0.5
        call system("clear")
        print *,"+---------------+","=================================="
        print *,"|               |","        < User Define >           "
        print *,"|               |","----------------------------------"
        print *,"|               |"," Year  : ",iyear1,"~",iyear2
        print *,"|   (Average)   |"," Month : ",imonth1,"~",imonth2
        print *,"|               |"," Lat   : ",ilat1,"~",ilat2 
        print *,"|               |"," Lon   : ",ilon1,"~",ilon2
        print *,"|               |"
        print *,"+---------------+","=================================="
        print *,"Target Longitude : ",(ilon1+ilon2)/2
        print *,"Target Latitude  : ",(ilat1+ilat2)/2
        print *,"+---------------+","=================================="
        write(*,*) "Start Loading dataset...."

        !! Initialize the accumulate matrix
        idepth(:)=0

        !! Read the CTD dataset 
        do i=1,ln
        read(10,*) lon(i),lat(i),year(i),month(i),day(i),hour(i),&
                   minute(i),sec(i),depth(i),temp(i),sal(i),&
                   density(i)

        !! Data Select in the User define range.
        if(((lon(i).ge.ilon1).and.(lat(i).ge.ilat1)).and.&
           ((lon(i).le.ilon2).and.(lat(i).le.ilat2)).and.&
           ((month(i).ge.imonth1).and.(month(i).le.imonth2)).and.&
           ((year(i).ge.iyear1).and.(year(i).le.iyear2))) then 

          !! Count times of data in each depth layer
          do l=1,33
            if (depth(i).eq.depth_index(l))then
                idepth(l) = idepth(l) + 1     
            endif
          end do

        endif

        !!! Processing % complete
          do m=1,10
            if(i.eq.999999*m)then
              proc = (real(i)/real(ln)) * 100
              write(*,1201) "Data Processing Complete: ",proc,"%"
1201          FORMAT(A26,F9.4,A3)
            end if
          end do

        end do

        !! Close the TWN-CTD dataset file
        close(10)

        !! Determine the new matrix size
        imax=max(idepth(1),idepth(2),idepth(3),idepth(4),idepth(5),&
               idepth(6),idepth(7),idepth(8),idepth(9),idepth(10),&
               idepth(11),idepth(12),idepth(13),idepth(14),idepth(15),&
               idepth(16),idepth(17),idepth(18),idepth(19),idepth(20),&
               idepth(21),idepth(22),idepth(23),idepth(24),idepth(25),&
               idepth(26),idepth(27),idepth(28),idepth(29),idepth(30),&
               idepth(31),idepth(32),idepth(33))
        allocate(ttemp(imax,layer))
        allocate(tsal(imax,layer))
        allocate(tdensity(imax,layer))

        call system("clear")
        print *, "Start Analysis .........."

        !! Initialize the calculator
        itemp(:)=0       !---------------!
        isal(:)=0        !  Data counts  !
        idensity(:)=0    !---------------!

        ttemp(:,:)=0     !---------------!
        tsal(:,:)=0      !  Data Summary !
        tdensity(:,:)=0  !---------------!

        do i=1,ln
        
        if(((lon(i).ge.ilon1).and.(lat(i).ge.ilat1)).and.&
           ((lon(i).le.ilon2).and.(lat(i).le.ilat2)).and.&
           ((month(i).ge.imonth1).and.(month(i).le.imonth2)).and.&
           ((year(i).ge.iyear1).and.(year(i).le.iyear2))) then
          
          do l=1,33
            if (depth(i).eq.depth_index(l))then
              itemp(l)=itemp(l)+1
              isal(l)=isal(l)+1
              idensity(l)=idensity(l)+1
              ttemp(itemp(l),l)=temp(i)
              tsal(isal(l),l)=sal(i)
              tdensity(idensity(l),l)=density(i)
            end if
          end do

        end if

        end do

        !! Display on the Screen
        write(*,*) "<<<-------- Temperature Average --------->>>"
        write(*,*) "Depth(m)     count      Temp.AVG(Degree)"
        do j=1,layer
        tcount(j)=0
        do i=1,imax
        if (ttemp(i,j).ne.0) tcount(j)=tcount(j)+1
        end do
        ftemp(j)=sum(ttemp(:,j))/real(tcount(j))
        if (tcount(j).eq.0) ftemp(j)=0
        write(*,1002) depth_index(j),'m',tcount(j),ftemp(j)
        end do

        write(*,*) "<<<<<-------- Salinity Average --------->>>>>"
        write(*,*) "Depth(m)     count     Sal.AVG(PSU)"
        do j=1,layer
        tcount(j)=0
        do i=1,imax
        if (tsal(i,j).ne.0) tcount(j)=tcount(j)+1
        end do
        fsal(j)=sum(tsal(:,j))/real(tcount(j))
        if (tcount(j).eq.0) fsal(j)=0
        write(*,1202) depth_index(j),'m',tcount(j),fsal(j)
        end do
1202    FORMAT(I4,A1,5X,I6,5X,F10.4)

        write(*,*) "<<<<<------ Pressure Average ------->>>>>"
        write(*,*) "Depth(m)     count      dbar.AVG(kg/m^2)"
        do j=1,layer
        icount=0
        do i=1,imax
        if (tdensity(i,j).ne.0) icount=icount+1
        end do
        fdensity(j)=(sum(tdensity(:,j))/real(icount))
        if (icount.eq.0) fdensity(j)=0
        write(*,1202) depth_index(j),'m',icount,fdensity(j)
        end do

        !! Open the output filename and write the result in.
        open(22,file='TWN-CTD_STATISTICS.txt',status='unknown',&
                action='write',iostat=ios)

        do j=1,layer
        T=ftemp(j)
        S=fsal(j)
        Dep=depth_index(j)
        fsv(j)=1449.2+(4.6*T)-&
               (0.055*(T**2))+&
               (0.0003*(T**3))+&
               (1.34-0.01*T)*(S-35)+&
               (0.016*Dep)

        ! If the temperature or salinity has no data, turn velocity to 0
        if (( ftemp(j) .eq. 0 ) .or. ( fsal(j) .eq. 0 )) fsv(j)=0
        if ( fdensity(j) .eq. 0 ) fdensity(j) = -1000

        write(22,1203) tlon,tlat,depth_index(j),tcount(j),ftemp(j),&
                      fsal(j),fdensity(j)+1000,fsv(j)
        end do
        close(22)
1203    FORMAT(2(F9.4,1X),2(I4,3X),4(F10.4,1X))

        print *, " * Process Finished. "
        print *, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
        print *, ">>>>   Output to TWN-CTD_STATISTICS.txt    >>>>"  
        print *, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
!!-------------------------------------------------------------------!!
!! <> End Case selection                                             !!
!!-------------------------------------------------------------------!!
        case default
        print *,"Under-Construction !!"
        end select 

        end program


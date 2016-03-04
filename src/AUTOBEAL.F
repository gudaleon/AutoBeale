cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  File name:    DBeale.f
C
C ***************************************************************
C *                                                             *
C * R. Peter Richards                                           *
C * Water Quality Lab                                           *
C * Heidelberg College                                          *
C * Tiffin, Ohio 44883                                          *
C * 419 448-2198                                                *
C *                                                             *
C *                                                             *
C *                                                             *
C * Macintosh implementation of IJC Beale ratio estimator       *
C * Program.                                                    *
C *                                                             *
C ***************************************************************
C
C  Purpose: A subroutine to calculate substance loadings using the stratified
C           Beale ratio estimator
C
C  Compile with:    (makefile beale.make)
C
C
C  Edit History:
C
C  ????      IJC Fortran program written             Kevin McGunagle         IJC
C  Dec-1993  Converted to MacFortran II, code "modernized"                   RPR
C            Code added in to allow time strata as well as flow              RPR
C            Code modified to allow discontinuous time strata, more useful   RPR
C            summary file.  Sort based on stratum number, assigned
C            in READAT. Other procedures reorganized...
C            Code added to calculate the user-chosen percentile of flow      RPR
C            for each year, and use it as the stratification point.
C  Nov 1994  Code added to use finite population correction, third order     RPR
C            variance estimate, and frational d.f. in calculating c.i.
C  Sept 1998 Conversion to non-graphics IBM version begun                    RPR
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine beale

      implicit none
      include 'common.fi'

      integer firstpass
            real*4         t05(30),tval1,tval2,frac
      real*8 fpc
            integer*4  kk
c  nstratblocks: number of stratum "areas" defined by user.  More than one area can be
c      assigned to the same stratum.  E.g. "Winter" may be Jan-April+Oct-Dec.
c  nstrata: number of actual strata, ²nstratblocks
      integer*4    i,l
      character*4  concunits,flowunits,loadunits
      character*1  ff

      data t05 /12.71,4.3,3.18,2.78,2.57,  2.45,2.37,2.31,2.26,2.23,
     & 2.2,2.18,2.16,2.14,2.13,   2.12,2.11,2.10,2.09,2.09,
     & 2.08,2.07,2.07,2.06,2.06,  2.06,2.05,2.05,2.05,2.04/
      ff=char(12)

1      format(a1,'TRIBUTARY NAME: ',a15/' YEAR: ',a12/' PARAMETER: ',a60/)
2      format(' Stratum',i3,/
     + ' This stratum includes dates from ',i8,' to ',i8)
22      format('      It also includes dates from ',i8,' to ',i8)
14      format(//'              Loadings',10x,'Flows     Concentrations'/
     + '    Date     ',a2,'/day     m3/sec       cfs',8x,a4/)
3      format(3x,i8,3(1x,f10.3),1x,f12.5)
4      format(//' Number of days in the stratum:',i4/
     + ' Mean stratum flow:',f10.3,' m3/sec or',f10.0,' cfs'//
     + ' Mean sample loading:',f15.1,1x,a2,'/day'/
     + ' Mean sample  flow:',f10.3,' m3/sec or',f10.0,' cfs'//
     + ' Ratio of mean stratum flow to mean sample flow:',f6.2//
     + ' Biased estimate:',f15.1,1x,a2,'/day'/' Unbiased estimate:',
     + f15.1,1x,a2,'/day'/' Bias correction:',f15.1,1x,a2,'/day'//
     + ' Mean square error:',f15.3,' (',a2,'/day)**2'/
     + ' based on',f6.0,' degrees of freedom'/
     + ' Sum of squares error:',f15.3,' (',a2,'/day)**2')
5      format(' No flow data available. Load calculation impossible')
6      format(a16,' No load can be calculated.  n < 2 for stratum',i3)
7      format(a1)

8      format('   River               Load, ',a2,'/yr         MSE       df         ±95% c.i.'/)
9      format(' Unstratified'/' The date range is from',i10,' to',i10//
     +       '            Loadings',12x,'Flows     Concentrations'/
     +       '    Date     ',a2,'/day     m3/sec       cfs',7x,a4,/)
10      format(//' Summary over',i4,' strata:'/
     +       ' Mean daily loading:',f15.3,' ',a2,'/day'/
     +       20x,f15.0,' ',a2,'/year'/20x,f15.3,1x,a2,'/year'//
     +       ' Mean square error:',f16.3,' (',a2,'/day)**2'/
     +       20x,f16.3,' (',a2,'/year)**2'/20x,f16.3,
     +       ' (',a2,'/year)**2'//
     +       ' based on',f10.3,' degrees of freedom')

11      format (a16,2f16.4,f10.3,f16.4)
12      format (//' Stratification applied:'//,
     +      '    block  stratum   start      end   '/)
13      format (2i8,2i9)
15      format (' The 95% confidence interval half-width is ',f16.4,
     &   1x,a2,'/year')

      select case (concfac)
      case (.000001d0);concunits='ng/L'
      case (.001d0);concunits='µg/L'
      case (1.0d0); concunits='mg/L'
      case (1000.0d0);concunits='g/L'
      case default
        write (*,*) 'Unexpected concentration units multiplier:', concfac
      end select

      if (loadfac.eq.1000.0d0) then
      loadunits='kg g'      !per year,per day
      else
      loadunits='mtkg'      !per year,per day
      end if

      if (flowfac.eq.1.0d0) then
      flowunits='cms'
      else
      flowunits='cfs'
      end if

      if (firstyear) write (3,8) loadunits(1:2)

      cumfl=0.0d0
      cumse=0.0d0
      df=0.0d0
      do 200 i=1,nstrata            !calculate stratum loads
      r(i)=0.0d0
      suml(i)=0.0d0
      suml2(i)=0.0d0
      sumf(i)=0.0d0
      sumf2(i)=0.0d0
      sumfl(i)=0.0d0
      sumx3(i)=0.0d0
      sumx2y(i)=0.0d0
      sumxy2(i)=0.0d0
      nf(i)=0
      flowmu(i)=0.0d0

      if (year2.ne.year(3:4)) year=year//'-'//year2
      if (bigio) then
        if (calculation_number.eq.1.and.i.eq.1) then
          write (2,1) char(32),oldtribname,year,param
        else
          write (2,1) ff,oldtribname,year,param
        end if
      end if

      if (bigio) then
        if(nstrata.gt.1) then
          firstpass=0
          do (kk=1,nstratblocks)
          if (strat(kk,5).eq.i) then
          if (firstpass.eq.0) then
            write (2,2) i,strat(kk,1),strat(kk,2)
            firstpass=1
          else
            write (2,22) strat(kk,1),strat(kk,2)
          end if
          end if
          repeat
          write (2,14) loadunits(3:4),concunits
        else
          write (2,9) strat(i,1),strat(i,2),loadunits(3:4),concunits
        end if
      end if

      do (l=1,ndays)
        if (stratum(l).ne.i) cycle
        nf(i)=nf(i)+1
        flowmu(i)=flowmu(i)+flow(l)
          if (conc(l) .ne. -1.0d0) then
          templ=flow(l)*conc(l)*86.4d0
          r(i)=r(i)+1.0d0
          suml(i)=suml(i)+templ
          suml2(i)=suml2(i)+templ*templ
          sumf(i)=sumf(i)+flow(l)
          sumf2(i)=sumf2(i)+flow(l)*flow(l)
          sumfl(i)=sumfl(i)+flow(l)*templ
          fl1=flow(l)/0.028316849d0                  !put it also in cfs for listing
          if (bigio) write (2,3) dates(l),templ,flow(l),fl1,conc(l)/(concfac*loadfac)
          end if
      repeat
      if (r(i).le.1.0d0) then      !insufficient data to calculate a load for stratum
        if (bigio) write (2,6)  oldtribname,strat(i,5)
        write (3,6) oldtribname,strat(i,5)
        needs_work=.true.
        go to 300
      end if
      if (nf(i).le.0) then            !no flow data - can't calculate a load
        if (bigio) write (2,5)
        write (3,5)
        call showalert(3, 'No flow data in one or more strata!'//
     &          'Please adjust the strata...')
        go to 300
      end if
      fpc=1.0d0-(r(i)/dble(nf(i)))
      fpc=1.0d0
      flowmu(i)=flowmu(i)/dble(nf(i))
      avflow=sumf(i)/r(i)
      avload=suml(i)/r(i)

      do (l=1,ndays)        !now calculate the third order terms
        if (stratum(l).ne.i) cycle
        if (conc(l).ne.-1.0d0) then
          sumx3(i)=sumx3(i) + (flow(l)-avflow)**3
          sumx2y(i)=sumx2y(i) + (flow(l)-avflow)**2 * (flow(l)*conc(l)*86.4d0-avload)
          sumxy2(i)=sumxy2(i) + (flow(l)-avflow) * (flow(l)*conc(l)*86.4d0-avload)**2
        end if
      repeat

      if(r(i) .gt. 1.0d0) then
        sxy=(sumfl(i)-(r(i)*avflow*avload))/(r(i)-1.0d0)
        sx2=(sumf2(i)-(r(i)*avflow*avflow))/(r(i)-1.0d0)
        sy2=(suml2(i)-(r(i)*avload*avload))/(r(i)-1.0d0)
        if (avload.gt.0.0d0) then
          sx2y=(sumx2y(i)/(r(i)-1.0d0))/(avflow**2*avload)
          sxy2=(sumxy2(i)/(r(i)-1.0d0))/(avflow*avload**2)
          sx3= (sumx3(i)/(r(i)-1.0d0))/(avflow**3)
        else
          sx3=0.0d0
          sx2y=0.0d0
          sxy2=0.0d0
        end if
      else
        sxy=0.0d0
        sx2=0.0d0
        sy2=0.0d0
        sx2y=0.0d0
        sxy2=0.0d0
        sx3=0.0d0
      end if

      if (avload.gt.0.0d0) then            !whew! some conc data
        f1=1.0d0+((fpc/r(i))*(sxy/(avload*avflow)))
        f2=1.0d0+((fpc/r(i))*(sx2/(avflow*avflow)))
        flowes(i)=avload*(flowmu(i)/avflow)*(f1/f2)
        s1=sx2/(avflow*avflow)
        s2=sy2/(avload*avload)
        s3=sxy/(avflow*avload)
        f1=(fpc/r(i))*(s1+s2-s3-s3)
        f2=(fpc**2/r(i)**2)*(2.0d0*s1*s1-4.0d0*s1*s3+s3*s3+s1*s2)
        f3=((2.0d0*fpc/r(i))/nf(i)) * (sx3-2.0d0*sx2y+sxy2)
        rmse(i)=(avload*flowmu(i)/avflow)**2*(f1+f2+f3)
        fl1=flowmu(i)/flowfac
        fl2=avflow/flowfac
        tmp1=flowmu(i)/avflow
        tmp2=avload*tmp1
        tmp3=flowes(i)-tmp2
        tmp4=rmse(i)*r(i)
      else                        !gads! no conc data in this stratum
        flowes(i)=0.0d0
        rmse(i)=0.0d0
        fl1=flowmu(i)/flowfac
        fl2=avflow/flowfac
        tmp1=flowmu(i)/avflow
        tmp2=0.0d0
        tmp3=0.0d0
        tmp4=0.0d0
      end if


      if (bigio) write (2,4) nf(i),flowmu(i),fl1,avload,loadunits(3:4),avflow,fl2,tmp1,tmp2,loadunits(3:4),flowes(i),
     +        loadunits(3:4),tmp3,loadunits(3:4),rmse(i),loadunits(3:4),r(i)-1.0d0,tmp4,loadunits(3:4),
200     continue

      do (i=1,nstrata)
      tmp1=r(i)-1.0d0
      rf=nf(i)        !convert nf to real to avoid integer overflow
      cumfl=cumfl+flowes(i)*(rf/rn)
      cumse=cumse+(rmse(i)*rf*rf/(rn*rn))
      df=df+(rf**4*rmse(i)**2/tmp1)
      repeat

        if(nstrata .gt. 1) then
         if (df.eq.0) df=-1
          df=(cumse*rn*rn)*(cumse*rn*rn)/df
          ce1=cumfl*365.0d0
          ce2=cumfl*0.365d0
          mse1=cumse*133225.0d0
          mse2=cumse*0.133225d0
          if (bigio) write (2,1) ff,oldtribname,year,param
          if (bigio) write (2,12)
          do (i=1,nstratblocks)            ! write a stratification summary
           if (bigio) write (2,13) i,strat(i,5),strat(i,1),strat(i,2)
          repeat

          if (bigio) write (2,10) nstrata,cumfl,loadunits(3:4),ce1,loadunits(3:4),ce2,loadunits(1:2),
     +        cumse,loadunits(3:4),mse1,loadunits(3:4),mse2,loadunits(1:2),df
        end if

! calculate the 95% confidence interval half-width (i.e. the ±number)
      if (df.ge.30) then
      tval=1.96d0+2.4d0/df
      else
      tval1=t05(idint(dmax1(1.0d0,df)))
      tval2=t05(idint(dmax1(1.0d0,df)+1))
      frac=df-idint(df)
      tval=tval1+frac*(tval2-tval1)
      end if
      if (nstrata.gt.1) then
      ci=tval*sqrt(mse2)
      else
      ci=tval*sqrt(rmse(1)*0.133225d0)
      end if
      if (bigio) write (2,15) ci,loadunits(1:2)

      if (nstrata.gt.1) then
      write (3,11) oldtribname,ce2,mse2,df,ci
      else
      write (3,11) oldtribname,flowes(1)*0.365d0,rmse(1)*0.133225d0,r(1)-1.,ci
      end if

300      return
      end

***************************************************************************


      subroutine read_data

      implicit none
      include 'common.fi'

      real*4 concin,concsum,cnt
      integer*4 i,junk,daynum,datein,wh,itest1,itest2,itest3
      character*8 localdate
      character*132 mess
      logical*1 nomess,yy


  1         format(f12.0)
  2         format(a48)
  7         format(i8)
  8      format (i2)
  122      format (a10,6x,i8,2(1x,f12.4))
  123      format (a10,6x,i8,1x,f12.4,13x)
      yy=.false.;nomess=.false.
      if (maxflow.gt.0) goto 1000      !flow data has already been read
      rn=0.0;ndays=0;maxconc=-1.0;maxflow=-1.0;firstday=0

c read the flow data for a given run, which determines the span of the run (first time only)

      do 100 i=1,367
      read (1,2,end=200,err=200) inputline
      inputline=inputline//repeat(' ',48)        !force line to be full
      ndays=ndays+1

      if (inputline(1:2).ne."19".and.inputline(1:2).ne.'20') then     !assure YYYYMMDD format
        if (.not.yy) then
          call showalert(1,'Flow file appears to use two-digit years, 20th century assumed!')
          yy=.true.
        end if
        inputline='19'//inputline
      end if
      if (i.eq.1) then
        oldfirstdate=inputline(1:8)
        read (oldfirstdate,7) itest1
        year=inputline(1:4)
      end if

      localdate=inputline(1:8)
          read (localdate,7) itest2
        itest3=itest2-itest1
        if (itest3.ge.10000) then
          ndays=ndays-1
!          write (*,*) itest1,'    ',localdate
          goto 200                  !stop when a full year is read....think about it
        end if
        if (i.gt.1) then
          itest3=itest2-dates(i-1)
        else
          itest3=1
        end if
        if (itest3.gt.1.and.itest3.lt.65) then
          write (mess,*) "No flow for ",dates(i-1)+1,'. There MUST be a flow for every day of the year!'//
     &          ' Aborting this run!'
          call showalert(3,mess)
          stop
        end if
        if (itest2.le.dates(i)) then
          write (mess,*) itest2,'  ',dates(i),'    Dates out of order (or duplicate dates) in flow file!'
          call showalert(3,mess)
          stop
        end if
      read (localdate,7) dates(i)
      read (inputline,*) junk,flow(i)        !will this work w/ spaces and/or tabs?!
      select case (flow(i))
        case (-1,-9,-3,-99,-999,-9999)
          write (mess,*) 'Missing value code for flow on ',dates(i),'. There MUST be a flow for every day of the year!'//
     &          ' Aborting this run!'
               call showalert(3,mess)
          stop
        case default
          flow(i)=flow(i)*flowfac
          rn=rn+1
      end select
!      if (flow(i).gt.0) then
!        flow(i)=flow(i)*flowfac
!        rn=rn+1
!      else
!        flow(i)=-1
!      end if
      if (flow(i).gt.maxflow) maxflow=flow(i)
100      continue

            write (*,*) 'WARNING - just read 367 flow values - too many for one year'
200         write (year2,8) dates(ndays)/10000


c read the concentrations for a given load calculation, converting the units as you go
c determine the maximum concentration for scaling purposes in the GUI

1000      daynum=1;yy=.false.
      do (i=1,366)
      conc(i)=-1.0d0
      repeat
      do 1100 i=1,1000
      read (4,2,end=1200) inputline
      inputline=adjustl(inputline)
      inputline=inputline//repeat(' ',48)        !force line to be full
      wh=index(inputline," ")
      if (wh.eq.0) then
        call showalert(3,'Separator missing in concentration file!')
        stop
      end if

      tribname=inputline(1:wh-1)//repeat(' ',40)
      if (i.eq.1) then
        oldtribname=tribname
      end if

      if (tribname.ne.oldtribname) then
        backspace 4
        goto 1300
      end if

      inputline=inputline(wh+1:48)//repeat(' ',48)
      inputline=adjustl(inputline)
      if (inputline(1:2).ne."19".and.inputline(1:2).ne.'20') then     !assure YYYYMMDD format
        if (.not.yy) then
          call showalert(1,'Conc file appears to use two-digit years, 20th century assumed!')
          yy=.true.
        end if
      end if
      if (yy) inputline='19'//inputline
      wh=index(inputline," ")
      if (wh.gt.9) inputline=inputline(1:8)//inputline(wh:50)//repeat(' ',16)    !dump time info
      read (inputline,*) datein, concin
      select case (concin)
        case (-1,-3,-9,-99,-999,-9999)
          goto 1100
        case default
          continue
      end select

      if (datein.lt.dates(daynum)) then
        write (mess,*) datein,'  ',dates(daynum),'   Dates out of order in concentration file!'
        call showalert(3,mess)
        stop
      end if

      if (datein.eq.dates(daynum).and.i.gt.1) then            !more conc sample(s) on the same date; average them
        if (cnt.eq.1.and.nomess) then
          nomess=.false.
          write (mess,*) datein
          mess='There is more than one chemical sample on '//trim(mess)
     &        //'.  Concentrations for this and similar days will being averaged.'
          call showalert(1,mess)
        end if
        concsum=concsum+concin
        cnt=cnt+1
        goto 1100
      else
        if (cnt.gt.1) then
          concsum=concsum/cnt
          conc(daynum)=concsum*concfac*loadfac
          lastday=daynum
        end if
        concsum=0;cnt=0
      end if

      do                        !find the day on which the conc belongs and put it there
        if (daynum.gt.ndays) exit            !this can happen if the conc file dates exceed the flow file...
        if (dates(daynum).eq.datein) then
          concsum=concin
          cnt=1
          conc(daynum)=concin*concfac*loadfac
          if (firstday.eq.0) firstday=daynum
          lastday=daynum
          exit
        end if
        daynum=daynum+1
      repeat
1100      continue

      write (*,*) 'WARNING - reached end of read-concentration loop without interruption...'
            go to 1300
1200         inputline(1:15)='end of file    '


1300      maxconc=-1
      do (i=1,366)
      if (conc(i).gt.maxconc) maxconc=conc(i)
      repeat

      return
      end


******************

      subroutine assign_stratmem

      implicit none
      include 'common.fi'

      integer*4 i,j
      character*256 warning


      do (i=1,ndays)            !assign stratum membership to each day's observations
      if (flow(i).eq.-1.0d0) then      !flows missing
        call showalert(2, 'missing flow in assign_stratmem...')
        write (warning,*) 'No flow for ',dates(i),'. There MUST be a flow for every day of the year!'//
     &          ' Aborting this run!'
             call showalert(3,warning)
!        stratum(i)=0
!        goto 405
      end if
      do (j=1,nstratblocks)
        if (dates(i).ge.strat(j,1).and.dates(i).le.strat(j,2)) then
          stratum(i)=strat(j,5)
          goto 405
        end if
      repeat
      write (warning,10) int(dates(i)),flow(i)/flowfac
      if (bigio) write (2,10) int(dates(i)),flow(i)/flowfac
      write (3,10) int(dates(i)),flow(i)/flowfac
10        format ('Assign_stratmem is unable to assign a stratum to ',i8,' '
     &          f9.3,'.')
           call showalert(3,warning)
           stop
405         repeat
            return
            end

***************************************************************************
! Calculate with no output; for interim report to GUI and other temporary uses
      subroutine bealecalc

      implicit none
      include 'common.fi'

            real*4         t05(30),tval1,tval2,frac
      real*8 fpc
      integer*4    i,l

      data    t05 /12.71,4.3,3.18,2.78,2.57,  2.45,2.37,2.31,2.26,2.23,
     &           2.2,2.18,2.16,2.14,2.13,   2.12,2.11,2.10,2.09,2.09,
     &             2.08,2.07,2.07,2.06,2.06,  2.06,2.05,2.05,2.05,2.04/


      cumfl=0.0d0
      cumse=0.0d0
      df=0.0d0
      do 200 i=1,nstrata            !calculate stratum loads
      r(i)=0.0d0
      suml(i)=0.0d0
      suml2(i)=0.0d0
      sumf(i)=0.0d0
      sumf2(i)=0.0d0
      sumfl(i)=0.0d0
      sumx3(i)=0.0d0
      sumx2y(i)=0.0d0
      sumxy2(i)=0.0d0
      nf(i)=0
      flowmu(i)=0.0d0

      do (l=1,ndays)
        if (stratum(l).ne.i) cycle
        nf(i)=nf(i)+1
        flowmu(i)=flowmu(i)+flow(l)
          if (conc(l) .ne. -1.0d0) then
          templ=flow(l)*conc(l)*86.4d0
          r(i)=r(i)+1.0d0
          suml(i)=suml(i)+templ
          suml2(i)=suml2(i)+templ*templ
          sumf(i)=sumf(i)+flow(l)
          sumf2(i)=sumf2(i)+flow(l)*flow(l)
          sumfl(i)=sumfl(i)+flow(l)*templ
          sumx3(i)=sumx3(i)+flow(l)**3
          sumx2y(i)=sumx2y(i)+flow(l)**2*templ
          sumxy2(i)=sumxy2(i)+flow(l)*templ**2
          end if
      repeat

      if (r(i).le.1.0d0) then      !insufficient data to calculate a load for stratum
        cint=1.0e20
        llod=-999
        if (i.eq.nstrata) llod=-888            !flag for iterations to say "you can stop now"
        go to 300
      end if
      if (nf(i).le.0) then            !no flow data - can't calculate a load
        cint=1.0e20
        llod=-999
        go to 300
      end if
      if (sumf(i).eq.0) then            !no flow on any sampling day - can't calculate a BRE load
        cint=1.0e20
        llod=-999
        go to 300
      end if
      fpc=1.0d0-(r(i)/dble(nf(i)))
      fpc=1.0d0
      flowmu(i)=flowmu(i)/dble(nf(i))
      avflow=sumf(i)/r(i)
      avload=suml(i)/r(i)

      if(r(i) .gt. 1.0d0) then
        sxy=(sumfl(i)-(r(i)*avflow*avload))/(r(i)-1.0d0)
        sx2=(sumf2(i)-(r(i)*avflow*avflow))/(r(i)-1.0d0)
        sy2=(suml2(i)-(r(i)*avload*avload))/(r(i)-1.0d0)

        if (avload.eq.0) then
          sx3=0.0d0
          sx2y=0.0d0
          sxy2=0.0d0
        else
          sx2y=sumx2y(i) - 2.0d0*sumf(i)*sumfl(i)/r(i) - sumf2(i)*suml(i)/r(i)
     &           + 2.0d0*sumf(i)**2*suml(i)/r(i)**2
          sx2y=(sx2y/(r(i)-1.0d0))/(avflow**2*avload)
          sxy2=sumxy2(i) - 2.0d0*suml(i)*sumfl(i)/r(i) - suml2(i)*sumf(i)/r(i)
     &           + 2.0d0*suml(i)**2*sumf(i)/r(i)**2
          sxy2=(sxy2/(r(i)-1.0d0))/(avflow*avload**2)
          sx3=sumx3(i) - 3.0d0*sumf2(i)*sumf(i)/r(i) + 2.0d0*sumf(i)**3/r(i)**2
          sx3=(sx3/(r(i)-1.0d0))/(avflow**3)
        end if
      else
        sxy=0.0d0
        sx2=0.0d0
        sy2=0.0d0
        sx2y=0.0d0
        sxy2=0.0d0
        sx3=0.0d0
      end if

      if (avload.gt.1.0d-7) then            !whew! some conc data
        f1=1.0d0+((fpc/r(i))*(sxy/(avload*avflow)))
        f2=1.0d0+((fpc/r(i))*(sx2/(avflow*avflow)))
        flowes(i)=avload*(flowmu(i)/avflow)*(f1/f2)
        s1=sx2/(avflow*avflow)
        s2=sy2/(avload*avload)
        s3=sxy/(avflow*avload)
        f1=(fpc/r(i))*(s1+s2-s3-s3)
        f2=(fpc**2/r(i)**2)*(2.0d0*s1*s1-4.0d0*s1*s3+s3*s3+s1*s2)
        f3=((2.0d0*fpc/r(i))/dble(nf(i))) * (sx3-2.0d0*sx2y+sxy2)
        rmse(i)=(avload*flowmu(i)/avflow)**2*(f1+f2+f3)
!        if (rmse(i).lt.-1.0d-8) then
!          call SysBeep(val2(0))
!          call SysBeep(val2(0))
!          write (*,*) 'in bealecalc, RMSE less than zero in stratum',i,'  RMSE',rmse(i)
!        end if
      else                        !gads! no conc data in this stratum, or all 0
        flowes(i)=0.0d0
        rmse(i)=0.0d0
      end if
200     continue

      do (i=1,nstrata)
      tmp1=r(i)-1.0d0
      rf=nf(i)        !convert nf to real to avoid integer overflow
      cumfl=cumfl+flowes(i)*(rf/rn)
      cumse=cumse+(rmse(i)*rf*rf/(rn*rn))
      df=df+(rf**4*rmse(i)**2/tmp1)
      repeat

        if(nstrata .gt. 1) then
      if (df.eq.0) df=-1
          df=(cumse*rn*rn)*(cumse*rn*rn)/df
          ce1=cumfl*365.0d0
          ce2=cumfl*0.365d0
          mse1=cumse*133225.0d0
          mse2=cumse*0.133225d0
        end if

! calculate the 95% confidence interval half-width (i.e. the ±number)
      if (df.ge.30) then
      tval=1.96d0+2.4d0/df
      else
      tval1=t05(idint(dmax1(1.0d0,df)))
      tval2=t05(idint(dmax1(1.0d0,df)+1))
      frac=df-idint(df)
      tval=tval1+frac*(tval2-tval1)
      end if
      if (nstrata.gt.1) then
      ci=tval*sqrt(mse2)
      else
      ci=tval*sqrt(rmse(1)*0.133225d0)
      end if

      if (nstrata.gt.1) then            !send back ci and either ce2 or flowes*.365d0 and df or r-1
      llod=ce2;degf=df
      else
      llod=flowes(1)*0.365d0;degf=r(1)-1.0
      end if
      cint=ci

300      return
      end
******************************
! Calculate with no output, only shifting the current observation from before stratum to after stratum
      subroutine bc2(dday)

      implicit none
      include 'common.fi'

      real*8     op,tval1,tval2,frac,ff1,ff2,fpc
            real*4         t05(30)
      integer*4  i,dday

      data    t05 /12.71,4.3,3.18,2.78,2.57,  2.45,2.37,2.31,2.26,2.23,
     &           2.2,2.18,2.16,2.14,2.13,   2.12,2.11,2.10,2.09,2.09,
     &             2.08,2.07,2.07,2.06,2.06,  2.06,2.05,2.05,2.05,2.04/



      cint=1.0e20
      cumfl=0.0d0
      cumse=0.0d0
      df=0.0d0
      do 200 i=mybound,mybound+1      !diferentially calculate stratum loads for current stratum and next stratum
      op=1.0d0;if (i.eq.mybound+1) op=-1.0d0
      if (i.le.0) write (*,*)  'i²0, i=',i
      if (i.ge.64) write (*,*) 'i³64,i=',i
      flowmu(i)=dble(nf(i))*flowmu(i)        !Get it back to a sum, rather than an average
      nf(i)=nf(i)+op
      flowmu(i)=flowmu(i)+op*flow(dday)
      if (conc(dday).ne.-1.0d0) then
        templ=flow(dday)*conc(dday)*86.4d0
        r(i)=r(i)+op*1.0d0
        suml(i)=suml(i)+op*templ
        suml2(i)=suml2(i)+op*templ*templ
        sumf(i)=sumf(i)+op*flow(dday)
        sumf2(i)=sumf2(i)+op*flow(dday)*flow(dday)
        sumfl(i)=sumfl(i)+op*flow(dday)*templ
        sumx3(i)=sumx3(i)+op*flow(dday)**3
        sumx2y(i)=sumx2y(i)+op*flow(dday)**2*templ
        sumxy2(i)=sumxy2(i)+op*flow(dday)*templ**2
      end if
      if (r(i).le.1.0d0) then      !insufficient data to calculate a load for stratum
        flowmu(i)=flowmu(i)/dble(nf(i))        !get it back like it was
        cint=1.0e20
        llod=-999
        if (i.eq.nstrata) llod=-888            !flag for iterations to say "you can stop now"
        go to 300
      end if
      if (nf(i).le.0) then            !no flow data - can't calculate a load
        cint=1.0e20
        llod=-999
        go to 300
      end if
      if (sumf(i).eq.0) then            !no flow on any sampling day - can't calculate a BRE load
        cint=1.0e20
        llod=-999
        go to 300
      end if
      fpc=1.0d0-(r(i)/dble(nf(i)))
      fpc=1.0d0
      flowmu(i)=flowmu(i)/dble(nf(i))
      avflow=sumf(i)/r(i)
      avload=suml(i)/r(i)

      if(r(i) .gt. 1.0d0) then
        sxy=(sumfl(i)-(r(i)*avflow*avload))/(r(i)-1.0d0)
        sx2=(sumf2(i)-(r(i)*avflow*avflow))/(r(i)-1.0d0)
        sy2=(suml2(i)-(r(i)*avload*avload))/(r(i)-1.0d0)
        if (avload.eq.0.0d0) then
          sx3=0.0d0
          sx2y=0.0d0
          sxy2=0.0d0
        else
          sx2y=sumx2y(i) - 2.0d0*sumf(i)*sumfl(i)/r(i) - sumf2(i)*suml(i)/r(i)
     &           + 2.0d0*sumf(i)**2*suml(i)/r(i)**2
          sx2y=(sx2y/(r(i)-1.0d0))/(avflow**2*avload)
          sxy2=sumxy2(i) - 2.0d0*suml(i)*sumfl(i)/r(i) - suml2(i)*sumf(i)/r(i)
     &           + 2.0d0*suml(i)**2*sumf(i)/r(i)**2
          sxy2=(sxy2/(r(i)-1.0d0))/(avflow*avload**2)
          sx3=sumx3(i) - 3.0d0*sumf2(i)*sumf(i)/r(i) + 2.0d0*sumf(i)**3/r(i)**2
          sx3=(sx3/(r(i)-1.0d0))/(avflow**3)
        end if
      else
        sxy=0.0d0
        sx2=0.0d0
        sy2=0.0d0
        sx2y=0.0d0
        sxy2=0.0d0
        sx3=0.0d0
      end if

      if (avload.gt.1.0d-7) then            !whew! some conc data (use 10e-7 as definition of 0 to avoid
        ff1=1.0d0+((fpc/r(i))*(sxy/(avload*avflow)))    !rounding error problems...  min. real load~.001*.05*86.4*.0283=10e-4
        ff2=1.0d0+((fpc/r(i))*(sx2/(avflow*avflow)))    !                 max. observed rounding error ~10e-9
        flowes(i)=avload*(flowmu(i)/avflow)*(ff1/ff2)
        s1=sx2/(avflow*avflow)
        s2=sy2/(avload*avload)
        s3=sxy/(avflow*avload)
        f1=(fpc/r(i))*(s1+s2-s3-s3)
        f2=(fpc**2/(r(i)*r(i)))*(2.0d0*s1*s1-4.0d0*s1*s3+s3*s3+s1*s2)
        f3=((2.0d0*fpc/r(i))/dble(nf(i))) * (sx3-2.0d0*sx2y+sxy2)
         rmse(i)=(avload*flowmu(i)/avflow)**2*(f1+f2+f3)
!        if (rmse(i).lt.-1.0d-8) then
!          write (*,*) 'In bc2, RMSE less than zero on day ',dday,', stratum',i,'  RMSE',rmse(i)
!        call SysBeep(val2(0))
!        end if
      else                        !gads! no conc data in this stratum, or all 0
        flowes(i)=0.0d0
        rmse(i)=0.0d0
      end if

200     continue

      do (i=1,nstrata)
      tmp1=r(i)-1.0d0
      rf=nf(i)        !convert nf to real to avoid integer overflow
      cumfl=cumfl+flowes(i)*(rf/rn)
      cumse=cumse+(rmse(i)*rf*rf/(rn*rn))
      df=df+(rf**4*rmse(i)**2/tmp1)
      repeat

        if(nstrata .gt. 1) then
      if (df.eq.0) df=-1
          df=(cumse*rn*rn)*(cumse*rn*rn)/df
          ce1=cumfl*365.0d0
          ce2=cumfl*0.365d0
          mse1=cumse*133225.0d0
          mse2=cumse*0.133225d0
        end if

! calculate the 95% confidence interval half-width (i.e. the ±number)

      if (df.ge.30) then
      tval=1.96d0+2.4d0/df
      else                              !interpolate for fractional degrees of freedom....
      tval1=t05(idint(dmax1(1.0d0,df)))
      tval2=t05(idint(dmax1(1.0d0,df)+1))
      frac=df-idint(df)
      tval=tval1+frac*(tval2-tval1)
      end if
      if (nstrata.gt.1) then
      ci=tval*sqrt(mse2)
      else
      ci=tval*sqrt(rmse(1)*0.133225d0)
      end if

      if (nstrata.gt.1) then            !send back ci and either ce2 or flowes*.365d0 and df or r-1
      llod=ce2;degf=df
      else
      llod=flowes(1)*0.365d0;degf=r(1)-1.0
      end if
      cint=ci

300      return
      end

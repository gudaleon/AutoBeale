
c  THIS SET OF ROUTINES PROBABLY NEED COMPLETE REVISION, OR AT LEAST EXTENSIVE REVISION...


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  File name:    io.f
c
c		   ***************************************************************
c		   *								 *
c		   *	R. Peter Richards					 *
c		   *	154 Morgan Street					 *
c		   *	Oberlin, Ohio 44883					 *
c		   *	216 775-1020						 *
c		   *							 	 *
c		   ***************************************************************
c
c  Purpose: 	
c
c  Compile with:	Buildprogram beale 
c
c  Edit History:
c    11/30/93		modified from SHAPE:shapeio.f			RPR
c	1994		MANY unrecorded changes and developments	RPR
c    10/26/94		flowfac changed from .0283 to 0.028316849	RPR
c    9/14/1998		Conversion to non-graphics IBM version begun	RPR
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     	SUBROUTINE basics_dialog

c This is needed but bring it back to a simple text dialog....

	
	implicit none
	include "common.fi"
	
	integer num
	character*1 lett
	


10	write (*,*)
	write (*,*) 'Choose units for concentration by entering the number:'
	write (*,*) '   1. g/L'
	write (*,*) '   2. mg/L'
	write (*,*) '   3. ug/L'
	read (*,*)  num
	select case (num)
	  case (1);concfac=1000.0d0
	  case (2);concfac=1.0d0
	  case (3);concfac=.001d0
	  case default
	    write (*,*) 'Inappropriate choice; please enter again!'
	    goto 10
	end select

20	write (*,*)
	write (*,*) 'Choose units for flow by entering the number:'
	write (*,*) '   1. cubic feet per second'
	write (*,*) '   2. cubic meters per second'
	read (*,*)  num
	select case (num)
	  case (1);flowfac=0.028316849d0
	  case (2);flowfac=1.0d0
	  case default
	    write (*,*) 'Inappropriate choice; please enter again!'
	    goto 20
	end select

30	write (*,*)
	write (*,*) 'Choose units for loads by entering the number:'
	write (*,*) '   1. metric tons per year'
	write (*,*) '   2. kilograms per year'
	read (*,*)  num
	select case (num)
	  case (1);loadfac=1.0d0
	  case (2);loadfac=1000.0d0
	  case default
	    write (*,*) 'Inappropriate choice; please enter again!'
	    goto 30
	end select

	write (*,*)
	write (*,*) 'Do you want a long output file as well as the summary results? (Y/N)'
	read (*,*) lett
	if (lett.eq.'y'.or.lett.eq.'Y') bigio=.true.

	param='Loads for the data in file "'//trim(concfilename)//'"'

	return
	end
	
*******************************

	subroutine showalert(num,text)
	
	implicit none
	include "common.fi"
	
	integer*4 num
	character*(*) text
	
	write (*,*) text
	
	return
	end	
	
	
**************************


	subroutine nextdate(yyyymmdd)
	
	implicit none
	
	integer*4 yyyymmdd,year,month,day,daysinmonth(12)
	
	data daysinmonth/31,28,31,30,31,30,31,31,30,31,30,31/

	if (yyyymmdd.lt.18000000) then
	  write (*,*) 'Illegal date:',yyyymmdd
	  write (*,*) 'Choose Exit from File menu to quit'
	  stop
	end if
	
	if (yyyymmdd.ge.24000000) then
	  write (*,*) 'Illegal date:',yyyymmdd
	  write (*,*) 'Choose Exit from File menu to quit'
	  stop
	end if
	
	year=yyyymmdd/10000
	month=yyyymmdd/100 - 100*(year)
	day=yyyymmdd - 100*month - 10000*(year)
	daysinmonth(2)=28
	if (mod(year,4).eq.0) daysinmonth(2)=29
	
	day=day+1
	
	if (day.gt.daysinmonth(month)) then
	  day=1
	  month=month+1
	  if (month.eq.13) then
	    month=1
	    year=year+1
	  end if
	end if
	    
	yyyymmdd=(year)*10000 + month*100 + day
 
	return
	end

**************************

	subroutine datedays(yyyymmdd,ndays)

c returns the number of days since 1900

	implicit none
	
	integer*4 j,ndays,yyyymmdd,year,month,day,daysinmonth(12)
	
	data daysinmonth/31,28,31,30,31,30,31,31,30,31,30,31/

	year=yyyymmdd/10000
	month=yyyymmdd/100 - 100*(year)
	day=yyyymmdd - 100*month - 10000*(year)
	daysinmonth(2)=28
	if (mod(year,4).eq.0) daysinmonth(2)=29

	ndays=365*(year-1900) + (year-1901)/4
	ndays=ndays+day
	do (j=1,month-1)
	  ndays=ndays+daysinmonth(j)
	repeat
	
	return
	end

**************************

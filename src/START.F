

C	WINDOW-RELATED STUFF PROBABLY NEEDS TO BE REVISED COMPLETELY
C	SOME ROUTINES CAN PROBABLY BE DUMPED ENTIRELY


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  File name:    Start.f
c
c		   ***************************************************************
c		   *								 *
c		   *    MORPHOGENESIS						 *
c		   *	R. Peter Richards					 *
c		   *	154 Morgan Street					 *
c		   *	Oberlin, Ohio 44883					 *
c		   *	216 775-1020						 *
c		   *							 	 *
c		   ***************************************************************
c
c  Purpose: 	Get us going with an initial window...
c
c  Compile with:	buildprogram beale
c
c  Edit History:
c    30-Nov-1993	Modified from SHAPE:start.f			RPR
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	Program start
	
c	Opens the window, calls beale subroutine(s)

	implicit none
	include 'common.fi'
	
	integer*4 ios,wh
	
	character*260 outfilename
      character*260 path
	character*260 name

c_________________________ End of declarations ___________________________

	bigio=.false.
	firstyear=.true.;firsttime=.true.;is_auto=.true.
	
! get the data file name

	path = "concfile.txt"

      open (4,err=4,iostat=ios,file=path,status='OLD')
      concfilename=path

! get the flow file name

	path = "flowfile.txt"
	open (1,err=5,iostat=ios,file=path,status='OLD')
	
	wh=index(concfilename,".")
	if (wh.gt.0) then
	  concfilename=concfilename(1:wh-1)
	else
	  concfilename=concfilename(1:len(trim(concfilename))-1)
	end if
	
	outfilename=trim(concfilename)//'_long.txt'
	if (bigio) open (2,file=trim(outfilename),filetype="TEXT",creator="????")

	outfilename=trim(concfilename)//'_summary.txt'
	open (3,file=trim(outfilename),filetype="TEXT",creator="????")


!!!! Now that we have an input file, we need to read it's contents, etc etc....	
10	format (a48)

	do						!for each year's data
	  doneflag=.false.
	  calculation_number=calculation_number+1	!needed by output files
	  call read_data				!read the year's data
	  if (firstday.eq.0) then
	    if (bigio) write (2,*) oldtribname,
     &	' No load was calculated because no concentration data were present'
	    write (3,*) oldtribname,
     &	' No load was calculated because no concentration data were present'
	    goto 120
	  end if
	  strat(1,1)=dates(1)
	  strat(1,2)=dates(ndays)
	  strat(1,5)=1
	  nstratblocks=1;nstrata=1
	  call assign_stratmem
	  call bealecalc
	  kci(1)=cint
!	  call SetCursor(wait)
!	  call ShowCursor()
	  needs_adjusted=.false.
	  currentstratum=1
!	  call PmBackColor(val2(backgroundColor))
!	  call EraseRect(new_my_rect)
!	  call drawtitle
	  do			!THE MAIN ENGINE TO FIND N-STRATUM MINIMUM ci
	    currentstratum=currentstratum+1;if (currentstratum.gt.64) exit
	    call pick_stratum_bound			!set up strata using sequential minimum approach
	    if (doneflag) then
	      if (.not.needs_adjusted) exit
	      call adjust_bounds
	      if (improved) then	!if there were improvements, let it try some more....
	        doneflag=.false.
		needs_adjusted=.false.
	      else
	        exit
	      end if
	    end if
	  end do
	  call beale					!calculate the loads and record them
	  oldtribname=tribname
	  firstyear=.false.
	  if (inputline(1:11).eq.'end of file') exit
120	end do

130	continue
	
	write (*,*) 'DONE!!!! Choose Exit from File menu to finish this run'
	end
	
**************************
	
	subroutine pick_stratum_bound
	
	implicit none
	include "common.fi"
	
	integer*4 k,day,strt
	integer*2 bestday,chemobs,curbound,didit
	real minci
	
	
	chemobs=0;didit=0
	do (k=1,ndays)			!Initialize dayci array to a large number
	  dayci(k)=1.0e20
	end do

	chemobs=1;curbound=1;strt=101+firstday
	
	mybound=curbound;day=firstday
	
10	call add_boundary(strt)		!add a new stratum boundary as soon as possible after the first observation

	if (mybound.eq.-1) then		!means add_boundary hit a boundary on the head and refused to act
	  strt=strt+2;day=day+2
	  curbound=curbound+1
	  goto 10
	end if
c Use a custom-counted loop because of complexities for which a do loop is not well suited
	do 
	  if (day.ge.lastday-1) goto 110	!Can't put a new boundary on lastday; also it provokes an error...
	  day=day+1				!move to the next day

c if the current day is within one of the end of the stratum, remove current stratum and
c shift into the next stratum block
	  if (dates(day+1).eq.strat(curbound+1,2)) then
	    chemobs=0;didit=0
!	    copyrect.left=day+99;copyrect.right=day+100
	    mybound=curbound
	    call remove_boundary
27	    curbound=curbound+1
	    day=day+3
	    if (day.ge.lastday) goto 115
28	    if (dates(day+1).eq.strat(curbound,2)) goto 27	!watch out you don't jump over a short stratum
29	    if (dates(day).eq.strat(curbound,2)) then		!watch out you don't jump on a short stratum boundary
	      curbound=curbound+1
	      day=day+2
	      goto 28
	    end if
	    
	    call add_boundary(day+100)
	    if (conc(day-1).ne.-1) chemobs=chemobs+1
	  end if
c else move the current stratum boundary to this day and calculate the estimator
c if we have enough chemistry
	  if (conc(day).ne.-1) chemobs=chemobs+1	!count up the chemical observations
	  if (flow(day).ne.-1) stratum(day)=curbound
	  if (chemobs.lt.2) goto 100
	  strat(curbound,2)=dates(day)
	  strat(curbound+1,1)=dates(day+1)
	  
	  if (didit.eq.0) then
	    call bealecalc
!	    write (8,*) 'Bealecalc ',day,curbound,mybound
	    if (llod.ge.0) didit=1
!	    if (llod.gt.0) didit=1
	  else
!	    write (8,*) 'BC2 ',day,dates(day),curbound,mybound
	    call bc2(day)
!	    call bealecalc 
	  end if
!	  if (llod.eq.-888.or.day.ge.ndays) then
!	    exit
!	  end if
	  if (llod.ge.0) then
	    dayci(day)=cint
	  end if
100	end do
	
110	mybound=curbound		!Remove last boundary used
	call remove_boundary
115	minci=1.0e20
	do (k=1,ndays)
	  if (dayci(k).lt.minci) then
	    minci=dayci(k)
	    bestday=k
	  end if
	end do
	
! place the stratum boundary at the best day
!	if (minci.lt.kci(currentstratum-1)) then
	if (minci.lt.0.995*kci(currentstratum-1)) then
	  call add_boundary(bestday+100)
	  kci(currentstratum)=minci
	  needs_adjusted=.true.
	else 	
	  doneflag=.true.
	  currentstratum=currentstratum-1
	  return
	end if

	return
	end

*****************************************

	subroutine adjust_bounds
	
	implicit none
	include "common.fi"
	
	integer*4 k,odaynum,d1,d2,j,day1,day2
	integer*4 bestday,chemobs,bestbound,bound(366),didit
	real minci,bestci
	
!	write (8,*) 'In adjust_bounds'
!	do (k=1,nstrata)
!	  write (8,*) k,int(strat(k,1)),int(strat(k,2))
!	end do
	bestci=kci(currentstratum);improved=.false.
	do				!until it no longer does any good
!	call drawdata
!	call drawstrata	
	chemobs=0
	do (k=1,ndays)			!Initialize dayci array to a large number
	  dayci(k)=1.0e20
	end do
	
	do (k=1,nstrata-1)
	  call datedays(strat(k,2),d2)
	  call datedays(strat(1,1),d1)
	  odaynum=1+(d2-d1)		!the original day number of the stratum which is being shifted through the range
	  mybound=k
	  call remove_boundary
	  
	  call datedays(strat(k,1),d2)
	  day1=d2-d1 + 1
	  call datedays(strat(k,2),d2)
	  day2=d2-d1
	  if (day1.le.firstday) day1=firstday+1
	  if (day2.ge.lastday)  day2=lastday-1

	  chemobs=0
	  if (conc(day1-1).ne.-1) chemobs=chemobs+1
	  call add_boundary(101+day1)		!Includes bealecalc
	  didit=0
	  do (j=day1+1,day2)
	    if (conc(j).ne.-1) chemobs=chemobs+1
	    if (flow(j).ne.-1) stratum(j)=k
	    if (chemobs.ge.2) then
	      strat(k,2)=dates(j)
	      strat(k+1,1)=dates(j+1)
	      
	      if (didit.eq.0) then
	        call bealecalc
		if (llod.ge.0) didit=1
!		if (llod.gt.0) didit=1
	      else
	         call bc2(j)
!	         call bealecalc
	      end if
	      
	      if (llod.ge.0) then
		if (dayci(j).gt.cint-bestci) then
		  dayci(j)=cint-bestci
		  bound(j)=mybound
		end if
	      end if
	    end if
	  end do
	  
	  mybound=k
	  call remove_boundary
	  call add_boundary(100+odaynum)		!put original boundary back for now
	end do
	minci=1.0e20
	do (k=1,ndays)
	  if (dayci(k).lt.minci) then
	    minci=dayci(k)
	    bestday=k
	    bestbound=bound(k)
	  end if
	end do
	
	
! move the stratum boundary to the best day, if any adjustment leads to improvement
!	if (minci.lt.0) then
	if (minci+bestci.lt.0.995*bestci) then
	  mybound=bestbound
	  call remove_boundary
	  call add_boundary(bestday+100)
	  bestci=minci+bestci
	  improved=.true.
	else 	
	  kci(currentstratum)=bestci		!Update best kci value so further comparisons will use right standard
	  return
	end if

	end do
	
	return
	end
	

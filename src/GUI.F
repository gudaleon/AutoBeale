
C	PROBABLY NEED TO KEEP LAST TWO SUBROUTINES

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  File name:    GUI.f
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
c  Purpose: 	Graphic user interface (GUI) for BEALE, now doing only a few of its former things...
c
c  Compile with:	buildprogram beale
c
c  Edit History:
c    5-Jan-1994		Split off from SHAPE:start.f			RPR
c  14-Sept-1998		Conversion to non-graphics IBM version begun	RPR
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


	subroutine add_boundary(loc)
	
c Splits a stratum into two and reassigns stratum numbers, preserving joined strata
c If you split half of a joined stratum, the joined pieces stay with the split piece nearest to them...UGH!

	implicit none
	include 'common.fi'
	
	integer*4 i,j,k,n,loc
	integer*4 z,w,old,curstrat
	integer*4 temp
	
	w=loc-100
	do (k=1,nstratblocks)
	  if (dates(w).gt.strat(k,1).and.dates(w).lt.strat(k,2)) then		!This one gets split
	    z=nstratblocks+1			!tack the new one on the end
	    strat(z,1)=dates(w+1)
	    strat(z,2)=strat(k,2)
	    strat(z,5)=z
	    strat(k,2)=dates(w)
	    old=k
	    goto 123				!Only need to do this for this stratum
	  end if
	repeat
!	call SysBeep();call SysBeep();call SysBeep()
	mybound=-1
	return
123	nstrata=nstrata+1
	nstratblocks=nstratblocks+1

	do (k=old+1,nstratblocks)			!Update stratum membership for later joined blocks
	  if (strat(k,5).eq.strat(old,5)) strat(k,5)=strat(z,5)
	repeat

	do (j=1,nstratblocks-1)		! sort by increasing date
	  do (k=j+1,nstratblocks)
	    if (strat(k,1).lt.strat(j,1)) then
	     do (i=1,5)
	      temp=strat(k,i)
	      strat(k,i)=strat(j,i)
	      strat(j,i)=temp
	     repeat
	    end if
	  repeat
	repeat
	
	do (k=1,nstratblocks)			!reassign stratum numbers
	  strat(k,5)=strat(k,5)+100	!temporarily increase them all
	repeat
	curstrat=1
	do (k=1,nstratblocks)
	  if (strat(k,5).gt.curstrat) then
	    temp=strat(k,5)
	    strat(k,5)=curstrat
	    do (n=k+1,nstratblocks)
	      if (strat(n,5).eq.temp) strat(n,5)=curstrat
	    repeat
	    curstrat=curstrat+1
	  end if
	repeat
		  
	call assign_stratmem
	call bealecalc

! redefine the stratum rectangles
	  curstrat=1
	  do (j=1,ndays)
	    if (dates(j).eq.strat(curstrat,1)) strrct(curstrat,2)=100+j	!begin date
	    if (dates(j).eq.strat(curstrat,2)) then
	      strrct(curstrat,4)=100+j					!end date
	      strrct(curstrat,1)=47
	      strrct(curstrat,3)=249
	      curstrat=curstrat+1
	    end if
	  repeat
	  strrct(curstrat,4)=101+ndays		!have to set the last stratum's end (and top and bot) separately
	  strrct(curstrat,1)=47			!since strat(last,2)=date(ndays+1)...
	  strrct(curstrat,3)=249

	approach=-1
	mybound=old
	blo=strrct(old,2)-100
	bhi=strrct(old+1,4)-100
	
	return
	end
	
****************************

c Joins two adjacent strata and reassigns stratum numbers, preserving other joined stratum memberships
c Requires that the boundary to remove be known by stratblock number, input through common as "mybound"
c Preceeded usually by a call to subroutine nearest, which finds mybound.

	subroutine remove_boundary
	
	implicit none
	include 'common.fi'
	
	integer*4 k,loc,temp,n
	integer*2 oldstratumnum,curstrat	

!	copyrect.top=50;copyrect.bottom=249	
!	copyrect.left=		must be set in calling routine
!	copyrect.right=		must be set in calling routine
!	call PmForeColor(val2(blackColor))
!	call PmBackColor(val2(whiteColor))
!	call CopyBits(val4(osPort.p_portPixMap),val4(my_window.port.p_portPixMap),
!     &		copyrect,copyrect,
!     &		val2(srcCopy + ditherCopy),val4(0))


	strat(mybound,2)=strat(mybound+1,2)
	oldstratumnum=strat(mybound+1,5)
		
	do (k=mybound+1,nstratblocks)		!move all other stratum blocks up one
	  strat(k,1)=strat(k+1,1)
	  strat(k,2)=strat(k+1,2)
	  strat(k,5)=strat(k+1,5)
	repeat
	nstratblocks=nstratblocks-1
	nstrata=nstrata-1

	do (k=mybound+1,nstratblocks)		!possibly reassign merged stratum numbers for this stratum
	  if (strat(k,5).eq.oldstratumnum) strat(k,5)=strat(mybound,5)
	repeat

	do (k=1,nstratblocks)			!reassign stratum numbers
	  strat(k,5)=strat(k,5)+100	!temporarily increase them all
	repeat
	curstrat=1
	do (k=1,nstratblocks)
	  if (strat(k,5).gt.curstrat) then
	    temp=strat(k,5)
	    strat(k,5)=curstrat
	    do (n=k+1,nstratblocks)
	      if (strat(n,5).eq.temp) strat(n,5)=curstrat
	    repeat
	    curstrat=curstrat+1
	  end if
	repeat


	call assign_stratmem	
	call bealecalc
	
	return
	end

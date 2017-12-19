        PROGRAM write_archive
c read   hycom archive file [ab]
c read   binary file after LETKF
c write  hycom archive file [ab]
        USE MOD_ZA
        IMPLICIT NONE
        EXTERNAL MINMAX
        integer check
        real v2d(1500,1100,3),tmp(1500,1100)
        integer layer, nn,nv,LLL
        character*240 file_in
        integer jversn,iexpt,yrflag
        REAL montg1(1500,1100)
        REAL srfhgt(1500,1100),surflx(1500,1100)
        REAL salflx(1500,1100),bl_dpth(1500,1100)
        REAL mix_dpth(1500,1100)
        REAL covice(1500,1100)
        REAL thkice(1500,1100)
        REAL temice(1500,1100)
        REAL u_btrop(1500,1100)
        REAL v_btrop(1500,1100)
        REAL u(1500,1100,32),v(1500,1100,32)
        REAL dp(1500,1100,32),temp(1500,1100,32)
        REAL saln(1500,1100,32)
        INTEGER MSK(1500,1100)
        real HMINA,HMAXA,dist,XMIN,XMAX
        integer i,j,k,ix,jx
        integer nstep
        double precision time(3)
        real thbase,thet,hminb,hmaxb
        real th(32)
        integer LL(32)
        character cline*80
        character ctitle(4)*80
        call XCSPMD
        CALL ZAIOST
       open(100,file="data.bin",form='unformatted',
     1 access='sequential')
c***************************************************
c    v2d(i,j,1) = SSH
c    v2d(i,j,2) = x_barotropic velocity  (not eastward) 
c    v2d(i,j,3) = y_barotropic velocity (not northward)
c    u(i,j,k) = x_baroclinic velocity (k=1 is top layer)
c    v(i,j,k) = y_baroclinic velocity  
c    dp(i,j,k) = layer thickness
c    temp(i,j,k) = temperature
c    saln(i,j,k) = salinity
c***************************************************
      do nn=1,3
      read(100) tmp
      do i=1,1500
      do j=1,1100
      v2d(i,j,nn)=tmp(i,j)
      enddo
      enddo
      enddo
      do LLL=1,32
        read(100) tmp
        do i=1,1500
        do j=1,1100
        u(i,j,LLL)=tmp(i,j)
        enddo
        enddo
      enddo
      do LLL=1,32
        read(100) tmp
        do i=1,1500
        do j=1,1100
        v(i,j,LLL)=tmp(i,j)
        enddo
        enddo
      enddo
      do LLL=1,32
        read(100) tmp
        do i=1,1500
        do j=1,1100
        dp(i,j,LLL)=tmp(i,j)*9806.0
        enddo
        enddo
      enddo
      do LLL=1,32
        read(100) tmp
        do i=1,1500
        do j=1,1100
        temp(i,j,LLL)=tmp(i,j)
        enddo
        enddo
      enddo
      do LLL=1,32
        read(100) tmp
        do i=1,1500
        do j=1,1100
        saln(i,j,LLL)=tmp(i,j)
        enddo
        enddo
      enddo
c***************************************************
c open old archive file
        CALL ZAIOPF('input.a','OLD',21)
        open(110,file='input.b',form='formatted',
     &  status='old',action='read')
c open new archive file
        CALL ZAIOPF('output.a','NEW',51)
        open(120,file='output.b',form='formatted',
     &  status='new',action='write')
        read(110,116) ctitle,jversn,iexpt,yrflag,idm,jdm
        write(120,116) ctitle,jversn,iexpt,yrflag,idm,jdm
 116    format (a80/a80/a80/a80/
     &   i5,4x,'''iversn'' = hycom version number x10'/
     &   i5,4x,'''iexpt '' = experiment number x10'/
     &   i5,4x,'''yrflag'' = days in year flag'/
     &   i5,4x,'''idm   '' = longitudinal array size'/
     &   i5,4x,'''jdm   '' = latitudinal  array size'/
     &   'field       time step  model day',
     &   '  k  dens        min              max')

 117  format (a8,' =',i11,f11.3,i3,f7.3,1p2e16.7)
        CALL ZAIORD(montg1,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        CALL ZAIOWR(montg1,MSK,.false.,xmin,xmax,51,.true.)
        write(120,117) 'montg1  ',nstep,time(1),layer,thbase,xmin,xmax

        CALL ZAIORD(srfhgt,MSK,.FALSE.,HMINA,HMAXA,21)
        write(*,*) 'SSH',srfhgt(500,1050),v2d(500,1050,1)*9.806
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        do i=1,1500
        do j=1,1100
        srfhgt(i,j)=v2d(i,j,1)*9.806
        enddo
        enddo
        CALL ZAIOWR(srfhgt,MSK,.false.,xmin,xmax,51,.true.)
        write(120,117) 'srfhgt  ',nstep,time(1),layer,thbase,xmin,xmax

        CALL ZAIORD(surflx,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        CALL ZAIOWR(surflx,MSK,.false.,xmin,xmax,51,.true.)
        write(120,117) 'surflx  ',nstep,time(1),layer,thbase,xmin,xmax

        CALL ZAIORD(salflx,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        CALL ZAIOWR(salflx,MSK,.false.,xmin,xmax,51,.true.)
        write(120,117) 'salflx  ',nstep,time(1),layer,thbase,xmin,xmax

        CALL ZAIORD(bl_dpth,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        CALL ZAIOWR(bl_dpth,MSK,.false.,xmin,xmax,51,.true.)
        write(120,117) 'bl_dpth ',nstep,time(1),layer,thbase,xmin,xmax

        CALL ZAIORD(mix_dpth,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        CALL ZAIOWR(mix_dpth,MSK,.false.,xmin,xmax,51,.true.)
        write(120,117) 'mix_dpth',nstep,time(1),layer,thbase,xmin,xmax

        CALL ZAIORD(covice,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        CALL ZAIOWR(covice,MSK,.false.,xmin,xmax,51,.true.)
        write(120,117) 'covice  ',nstep,time(1),layer,thbase,xmin,xmax

        CALL ZAIORD(thkice,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        CALL ZAIOWR(thkice,MSK,.false.,xmin,xmax,51,.true.)
        write(120,117) 'thkice  ',nstep,time(1),layer,thbase,xmin,xmax

        CALL ZAIORD(temice,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        CALL ZAIOWR(temice,MSK,.false.,xmin,xmax,51,.true.)
        write(120,117) 'temice  ',nstep,time(1),layer,thbase,xmin,xmax

        CALL ZAIORD(u_btrop,MSK,.FALSE.,HMINA,HMAXA,21)
        write(*,*) 'UB',u_btrop(500,1050),v2d(500,1050,2) 
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        do i=1,1500
        do j=1,1100
        u_btrop(i,j)=v2d(i,j,2)
        enddo
        enddo
        CALL ZAIOWR(u_btrop,MSK,.false.,xmin,xmax,51,.true.)
        write(120,117) 'u_btrop  ',nstep,time(1),layer,thbase,xmin,xmax

        CALL ZAIORD(v_btrop,MSK,.FALSE.,HMINA,HMAXA,21)
        write(*,*) 'VB',v_btrop(500,1050),v2d(500,1050,3) 
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        do i=1,1500
        do j=1,1100
        v_btrop(i,j)=v2d(i,j,3)
        enddo
        enddo
        CALL ZAIOWR(v_btrop,MSK,.false.,xmin,xmax,51,.true.)
        write(120,117) 'v_btrop  ',nstep,time(1),layer,thbase,xmin,xmax

        do k=1,32
        read(110,'(a)') cline
        read(110,'(a)') cline
        read(110,'(a)') cline
        read(110,'(a)') cline
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),LL(k),th(k),hminb,hmaxb
        write(*,*) LL(k),th(k)
        enddo

        do k=1,32
        CALL ZAIOWR(u(1,1,k),MSK,.false.,xmin,xmax,51,.true.)
        write(120,117) 'u-vel.  ',nstep,time(1),k,th(k),xmin,xmax
        CALL ZAIOWR(v(1,1,k),MSK,.false.,xmin,xmax,51,.true.)
        write(120,117) 'v-vel.  ',nstep,time(1),k,th(k),xmin,xmax
        CALL ZAIOWR(dp(1,1,k),MSK,.false.,xmin,xmax,51,.true.)
        write(120,117) 'thknss  ',nstep,time(1),k,th(k),xmin,xmax
        CALL ZAIOWR(temp(1,1,k),MSK,.false.,xmin,xmax,51,.true.)
        write(120,117) 'temp    ',nstep,time(1),k,th(k),xmin,xmax
        CALL ZAIOWR(saln(1,1,k),MSK,.false.,xmin,xmax,51,.true.)
        write(120,117) 'salin   ',nstep,time(1),k,th(k),xmin,xmax
        enddo
        check=0
        if (check.eq.1) then
      write(*,*) 'u'
      do layer=1,32
      write(*,*) layer,u(500,1050,layer) 
      enddo
      write(*,*) 'v'
      do layer=1,32
      write(*,*) layer,v(500,1050,layer) 
      enddo
      write(*,*) 'dp'
      do layer=1,32
      write(*,*) layer,dp(500,1050,layer) 
      enddo
      write(*,*) 't'
      do layer=1,32
      write(*,*) layer,temp(500,1050,layer) 
      enddo
      write(*,*) 's'
      do layer=1,32
      write(*,*) layer,saln(500,1050,layer) 
      enddo
        endif
        CALL ZAIOCL(21)
        CALL ZAIOCL(51)
        close(21)
        close(51)
        close(110)
        close(120)
        stop
        end

        PROGRAM read_archive
c read  hycom archive file [ab]
c write binary file 
        USE MOD_ZA
        IMPLICIT NONE
        EXTERNAL MINMAX
        real v2d(1500,1100,3),tmp(1500,1100)
        real v3d(1500,1100,32,5)
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
        integer LL(32)
        real th(32)
        character cline*80
        character ctitle(4)*80
        call XCSPMD  
        CALL ZAIOST
c***************************************************
c open old archive file
        CALL ZAIOPF('input.a','OLD',21)
        open(110,file='input.b',form='formatted',
     &  status='old',action='read')
c open new archive file
        read(110,116) ctitle,jversn,iexpt,yrflag,idm,jdm
c        write(120,116) ctitle,jversn,iexpt,yrflag,idm,jdm
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

        CALL ZAIORD(srfhgt,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        do i=1,1500
        do j=1,1100
        v2d(i,j,1)=srfhgt(i,j)/9.806
        enddo
        enddo

        CALL ZAIORD(surflx,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        CALL ZAIORD(salflx,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        CALL ZAIORD(bl_dpth,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        CALL ZAIORD(mix_dpth,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        CALL ZAIORD(covice,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        CALL ZAIORD(thkice,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        CALL ZAIORD(temice,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        CALL ZAIORD(u_btrop,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        do i=1,1500
        do j=1,1100
        v2d(i,j,2)=u_btrop(i,j)
        enddo
        enddo

        CALL ZAIORD(v_btrop,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        do i=1,1500
        do j=1,1100
        v2d(i,j,3)=v_btrop(i,j)
        enddo
        enddo
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
        do k=1,32
        CALL ZAIORD(u(1,1,k),MSK,.FALSE.,HMINA,HMAXA,21)
        CALL ZAIORD(v(1,1,k),MSK,.FALSE.,HMINA,HMAXA,21)
        CALL ZAIORD(dp(1,1,k),MSK,.FALSE.,HMINA,HMAXA,21)
        CALL ZAIORD(temp(1,1,k),MSK,.FALSE.,HMINA,HMAXA,21)
        CALL ZAIORD(saln(1,1,k),MSK,.FALSE.,HMINA,HMAXA,21)
        enddo

c write binary file
       open(100,file="data.bin",form='unformatted',
     1 access='sequential')
        do k=1,3
        do i=1,1500
        do j=1,1100
        tmp(i,j)=v2d(i,j,k)
        enddo
        enddo
        write(100) tmp
        enddo

        do k=1,32
        do i=1,1500
        do j=1,1100
        tmp(i,j)=u(i,j,k)
        enddo
        enddo
        write(100) tmp
        enddo
        do k=1,32
        do i=1,1500
        do j=1,1100
        tmp(i,j)=v(i,j,k)
        enddo
        enddo
        write(100) tmp
        enddo
        do k=1,32
        do i=1,1500
        do j=1,1100
        tmp(i,j)=dp(i,j,k)/9806.
        enddo
        enddo
        write(100) tmp
        enddo
        do k=1,32
        do i=1,1500
        do j=1,1100
        tmp(i,j)=temp(i,j,k)
        enddo
        enddo
        write(100) tmp
        enddo
        do k=1,32
        do i=1,1500
        do j=1,1100
        tmp(i,j)=saln(i,j,k)
        enddo
        enddo
        write(100) tmp
        enddo
        close(100)
        stop
        end

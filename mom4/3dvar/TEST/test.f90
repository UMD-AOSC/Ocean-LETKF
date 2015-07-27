program test
INTEGER :: i,j,k
INTEGER :: isd=0,ied=721,jsd=0,jed=411
INTEGER :: ni=720,nj=410,nk=1

do k=1,nk
  do j=jsd,jed
    do i=isd,ied
      if (i .ne. MODULO(i-1,ni)+1 .OR. j .ne. MODULO(j-1,nj)+1) then
        print *, i,j,k," :: ", MODULO(i-1,ni)+1,MODULO(j-1,nj)+1,k
      endif
     !T_prog(1)%field(i,j,k,taup1) = temp(i,j,k)
     !T_prog(2)%field(i,j,k,taup1) = salt(i,j,k)
    enddo
  enddo
enddo


end program test


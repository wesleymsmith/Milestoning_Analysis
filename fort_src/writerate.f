program writerate

c====================================================
c       a program to bin the trj with respect to the 
c       cells defined by the string images 
c====================================================

        implicit none
        integer nimage,nxx,jmin,jmin0,imother,tmother
        !nimage: number of cells 
        parameter(nimage=20)
        parameter(nxx=10000000)
        integer icell(nimage),k,kk,ii,j,i,ij
        double precision mc1,cd1,cdist(nimage)
        double precision mindist,mindistr
        double precision c1av0(nimage),nmatr(nimage,nimage)

        open(10,status='old',file='imother.txt')
          read(10,*) imother
        close(10)

        do i=1,nimage
         do j=1,nimage
          nmatr(i,j)=0.d0
         enddo
        enddo


        jmin=999
        jmin0=jmin
        tmother=0

        open(10,status='old',file='fort.37')

        do i=1,nxx

          read(10,*,END=20) jmin

           if((jmin.ne.imother).and.(jmin0.eq.imother)) then
            nmatr(imother,jmin)=nmatr(imother,jmin)+1   ! escape matrix (off-diag)
           endif

           if(jmin.eq.imother) then
            tmother=tmother+1               ! time in cell
           endif

         jmin0=jmin

        enddo
        
20      continue

          write(38,100) (dble(nmatr(imother,j))/dble(tmother), 
     c     j=1,nimage)

100     format(20(f8.6,2x))

        end


          write(38,*) (dble(nmatr(imother,j))/dble(tmother), j=1,nimage)

        end
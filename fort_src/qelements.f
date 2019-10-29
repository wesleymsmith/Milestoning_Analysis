        program qelements

c====================================================
c       a program to bin the trj
c       with respect to the walls btwn cells
c====================================================

        implicit none
        integer nimage,nxx,jmin,jmin0,n12,n21,n11,n22,r1,r2
        integer nimagem,nimagep
        integer nm,nm0
        parameter(nxx=200)
        integer j,i


c=======================
c set counters

        jmin=999
        nm  = 99999
        n12 = 0
        n21 = 0
        n11 = 0
        n22 = 0
        r1  = 0
        r2  = 0
c=======================

c=======================
c define neighbor cells

c must know which cell is being analyzed, can be put in first line of
cell index traj
        open(5,status='old',file='fort.5')
        read(5,*) nimage
c        write(6,*) nimage

        nimagep=nimage+1
        nimagem=nimage-1
c=======================
        write(*,*) "starting up"

        do i=1,nxx
          jmin0=jmin
          nm0=nm

          read(5,*,END=20) jmin

             if((jmin0.eq.nimagem).and.(jmin.eq.nimage)) nm=1
             if((jmin0.eq.nimage).and.(jmin.eq.nimagem)) nm=1
             if((jmin0.eq.nimage).and.(jmin.eq.nimagep)) nm=2
             if((jmin0.eq.nimagep).and.(jmin.eq.nimage)) nm=2

           if((nm0.eq.1).and.(nm.eq.1).and.(jmin.eq.nimage)) r1=r1+1
           if((nm0.eq.2).and.(nm.eq.2).and.(jmin.eq.nimage)) r2=r2+1

           if((nm0.eq.1).and.(nm.eq.2)) n12=n12+1
           if((nm0.eq.2).and.(nm.eq.1)) n21=n21+1


        enddo

20      continue

c  fort.39 will be tmpallN
c  fort.40 will be tmpallR

          open(39,file="tmpallN")
          open(40,file="tmpallR")
           write(39,*) n12,n21
           write(40,*) r1,r2
          close(39)
          close(40)

        end

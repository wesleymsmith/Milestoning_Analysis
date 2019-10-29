      program writeindex

c=================================
c       a program to bin the free traj
c       with respect to the cells defined by
c       images loaded in ca1v0 
c=================================

        implicit none
        integer nimage,nxx,jmin

c       nimage is number of cells
        parameter(nimage=17)
        parameter(nxx=200)
        integer icell(nimage),k,kk,ii,j,i,ij
        double precision mc1,mc2,cd1,cd2,cdist(nimage)
        double precision mindist,mindistr
        double precision c1av0(nimage),c2av0(nimage)


         ! open trajectory
         open(10,status='old',file='tmp.traj')

        jmin=999

        do j=1,nimage
        icell(j)=0
        enddo

        mindistr=0.d0
 
        ! centers (AKA images) values here
     
        c1av0(17)=  4.0d0
        c1av0(16)=  3.5d0
        c1av0(15)=  3.0d0
        c1av0(14)=  2.5d0
        c1av0(13)=  2.0d0
        c1av0(12)=  1.5d0
        c1av0(11)=  1.0d0
        c1av0(10)=  0.5d0
        c1av0( 9)=  0.0d0
        c1av0( 8)= -0.5d0
        c1av0( 7)= -1.0d0
        c1av0( 6)= -1.5d0
        c1av0( 5)= -2.0d0
        c1av0( 4)= -2.5d0
        c1av0( 3)= -3.0d0
        c1av0( 2)= -3.5d0
        c1av0( 1)= -4.0d0
        
        
        do i=1,nxx
          ! read trajectory
          read(10,*) ij,mc1

          mindist=999999.d0

           do j=1,nimage         !  for every frame of trj, 
             cd1=mc1-c1av0(j)    ! loop on images to find the closest
             cdist(j) = dsqrt(cd1**2)
             mindist=min(mindist,cdist(j))
             if(mindist.ne.mindistr) jmin=j   ! find the cell where the sampler is
             mindistr=mindist
           enddo

          write(37,*) jmin

        enddo

        end

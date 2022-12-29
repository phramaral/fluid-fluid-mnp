PROGRAM grafted

    IMPLICIT NONE
    REAL :: R,dtheta,dphi,pi,x,y,z,atotal,amol,dist,theta,phi,dij,radius
    INTEGER :: N_beads,Mphi,Mtheta,N_count,NIONS
    INTEGER :: i,j,k,l
    REAL, ALLOCATABLE :: pos(:,:)
    REAL, ALLOCATABLE :: coord(:,:)
    character(4), allocatable :: symbols(:)
    REAL, ALLOCATABLE :: tmp_angle(:,:)
    REAL, ALLOCATABLE :: angle(:,:)
    REAL :: v1(3),v2(3)
    character(256)    :: xyzfg="etilenoglicol-O.xyz"         ! file .xyz of functional group
    !character(256)    :: xyzfg="polietilenoglicol-O.xyz"         ! file .xyz of functional group
    !character(256)    :: xyzfg="stearic_acid-O.xyz"         ! file .xyz of functional group
    !character(256)    :: xyznp="np_Fe3O4.xyz"                      ! file .xyz of nanoparticle
    character(256)    :: xyznp="np-anchor_points.xyz"                      ! file .xyz of nanoparticle  
    character(4) :: string
    character(4), allocatable :: element_fg(:),element_np(:)
    integer :: natoms_group,iatom,nancor,natoms_np,natoms
    real :: rx,ry,rz
    real, allocatable :: r_group(:,:),r_np(:,:),r_ancor(:,:),r_f(:,:)
    real, allocatable :: r_rotate(:,:)
    
    open(8,file=xyzfg)
    read(8,*)natoms_group
    read(8,*)
    allocate(r_group(natoms_group,3))
    allocate(element_fg(natoms_group))
    allocate(r_rotate(natoms_group,3))
    
    do iatom=1,natoms_group
       read(8,*)string,rx,ry,rz
       element_fg(iatom)=string
       r_group(iatom,1)=rx
       r_group(iatom,2)=ry
       r_group(iatom,3)=rz
    end do
    close(8)

    nancor=154
    
    open(8,file=xyznp)
    read(8,*)natoms_np
    read(8,*)
    allocate(r_np(natoms_np,3))
    allocate(element_np(natoms_np))
    allocate(r_ancor(nancor,3))
    i=0
    do iatom=1,natoms_np
       read(8,*)string,rx,ry,rz
       element_np(iatom)=string
       r_np(iatom,1)=rx
       r_np(iatom,2)=ry
       r_np(iatom,3)=rz
       if(string == "ap") then
         i=i+1
         r_ancor(i,1)=rx
         r_ancor(i,2)=ry
         r_ancor(i,3)=rz 
         write(6,*)string,rx,ry,rz
       end if
    end do
    close(8)
    
    allocate(r_f(natoms_group,3))
    allocate(coord(nancor*natoms_group,3))
    allocate(symbols(nancor*natoms_group))
    
    k=0      
    do j=1,nancor
        
       r_f=0
       R=sqrt(r_ancor(j,1)*r_ancor(j,1)+r_ancor(j,2)*r_ancor(j,2)+r_ancor(j,3)*r_ancor(j,3))   
       radius=R+2.08
       theta=acos(r_ancor(j,3)/R)
       phi=atan2(r_ancor(j,2),r_ancor(j,1))
       
       ! transalate the group in the z direction on sphere surface
       do iatom=1,natoms_group
          r_f(iatom,1)=r_group(iatom,1)
          r_f(iatom,2)=r_group(iatom,2)
          r_f(iatom,3)=r_group(iatom,3)+radius
       end do
     
       call rotation(r_f,r_rotate,natoms_group,theta,phi)
       do iatom=1,natoms_group
          write(6,*)element_fg(iatom),r_rotate(iatom,1),r_rotate(iatom,2),r_rotate(iatom,3)
          k=k+1
          symbols(k)=element_fg(iatom)
          coord(k,1)=r_rotate(iatom,1)
          coord(k,2)=r_rotate(iatom,2)
          coord(k,3)=r_rotate(iatom,3)
       end do
    end do
    
    natoms=natoms_group*nancor+natoms_np
    
    1000 format (I5)
    2000 format (A11)
    3000 format (A4,2X,F17.12,2X,F17.12,2X,F17.12)
    
    !open(unit=9,file="magnetite_np_grafted_O.xyz",action="write")
    open(unit=9,file="np_grafted-etilenoglicol-O.xyz",action="write")
    
    
    write(9,1000)natoms
    write(9,2000)"magnetite"
    do iatom=1,natoms_np-nancor
       string=element_np(iatom)
       rx=r_np(iatom,1)
       ry=r_np(iatom,2)
       rz=r_np(iatom,3)
       !if (string /= "Fe5") then
         write(9,3000)string,rx,ry,rz
       !end if
    end do
    
    l=0
    do i=1,nancor
       j=natoms_np-nancor+i
       string=element_np(j)
       rx=r_np(j,1)
       ry=r_np(j,2)
       rz=r_np(j,3)
       write(9,3000)string,rx,ry,rz 
       do k=1,natoms_group 
          l=l+1
          string=symbols(l)
          rx=coord(l,1)
          ry=coord(l,2)
          rz=coord(l,3)
          write(9,3000)string,rx,ry,rz
       end do
    end do
    
    close(9)
    
    contains
    
    subroutine rotation(positions,rotate,natoms,theta,phi)
    implicit none
    integer, intent(in)    :: natoms
    real,    intent(in)    :: positions(natoms,3)
    real,    intent(in)    :: theta
    real,    intent(in)    :: phi
    real,    intent(out)   :: rotate(natoms,3)
    real :: coord(natoms,3)
    integer :: iatom,i,j,k
    
    coord=0.0
    
    ! rotation 'y' theta  angle
    do iatom=1,natoms
       coord(iatom,1)=cos(theta)*positions(iatom,1)+sin(theta)*positions(iatom,3)
       coord(iatom,2)=positions(iatom,2)
       coord(iatom,3)=-1.0*sin(theta)*positions(iatom,1)+cos(theta)*positions(iatom,3)
    end do
    
    ! rotation 'z' phi  angle
    do iatom=1,natoms
       rotate(iatom,1)=cos(phi)*coord(iatom,1)-sin(phi)*coord(iatom,2)
       rotate(iatom,2)=sin(phi)*coord(iatom,1)+cos(phi)*coord(iatom,2)
       rotate(iatom,3)=coord(iatom,3)
    end do
    
    end subroutine rotation
    
    END PROGRAM grafted
    
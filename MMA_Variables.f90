module MMA_variables
    implicit none
    integer                               :: m                        
    integer                               :: n                        
    integer                               :: kkttol                   
    integer                               :: maxoutit                 
    integer                               :: outeriter                
    integer                               :: outit                   
    real                                  :: a0                       
    real                                  :: epsimin                   
    real                                  :: kktnorm                  
    real                                  :: residumax                
    real, dimension(:), allocatable       :: residu                   
    real, dimension(:), allocatable       :: outvector1               
    real, dimension(:), allocatable       :: outvector2               
    real, dimension(:,:), allocatable     :: a                         
    real, dimension(:,:), allocatable     :: c                         
    real, dimension(:,:), allocatable     :: d                         
    real, dimension(:,:), allocatable     :: xmma                      
    real, dimension(:,:), allocatable     :: ymma                      
    real                                  :: zmma                      
    real, dimension(:,:), allocatable     :: lam                       
    real, dimension(:,:), allocatable     :: xsi                       
    real, dimension(:,:), allocatable     :: eta                       
    real, dimension(:,:), allocatable     :: mu                        
    real                                  :: zet                       
    real, dimension(:,:), allocatable     :: s                         
    real                                  :: f0val                     
    real, dimension(:,:), allocatable     :: fval                      
    real, dimension(:,:), allocatable     :: xval                      
    real, dimension(:,:), allocatable     :: xmax                      
    real, dimension(:,:), allocatable     :: xmin                      
    real, dimension(:,:), allocatable     :: xold1                      
    real, dimension(:,:), allocatable     :: xold2                     
    real, dimension(:,:), allocatable     :: df0dx                     
    real, dimension(:,:), allocatable     :: dfdx                      
    real, dimension(:,:), allocatable     :: low                       
    real, dimension(:,:), allocatable     :: upp                       

contains

    subroutine ObjectiveFunction(x,f0val,df0dx,fval,dfdx)
        implicit none
        real, dimension(:,:), allocatable, intent(inout)    :: x
        real, intent(inout)                                 :: f0val
        real, dimension(:,:), allocatable, intent(inout)    :: df0dx
        real, dimension(:,:), allocatable, intent(inout)    :: fval
        real, dimension(:,:), allocatable, intent(inout)    :: dfdx
        ! ----------------------------------------------------------- !
        !  This file calculates function values and gradients         !
        !  for the following "toy problem":                           !
        !                                                             !
        !    minimize x(1)^2 + x(2)^2 + x(3)^2                        !
        !  subject to (x(1)-5)^2 + (x(2)-2)^2 + (x(3)-1)^2 =< 9       !
        !             (x(1)-3)^2 + (x(2)-4)^2 + (x(3)-3)^2 =< 9       !
        !              0 =< x(j) =< 5, for j=1,2,3.                   !
        ! ----------------------------------------------------------- !
        ! f0val = x(1)^2 + x(2)^2 + x(3)^2;
        f0val = sum(x**2)
        
        ! df0dx = [2*x(1)
	    !          2*x(2)
	    !          2*x(3)];
        df0dx = 2*x

        !fval  = [(x(1)-5)^2+(x(2)-2)^2+(x(3)-1)^2-9
        !         (x(1)-3)^2+(x(2)-4)^2+(x(3)-3)^2-9];
        fval(1,1) = (x(1,1)-5)**2 + (x(2,1)-2)**2 + (x(3,1)-1)**2 - 9.0
        fval(2,1) = (x(1,1)-3)**2 + (x(2,1)-4)**2 + (x(3,1)-3)**2 - 9.0
        
        !dfdx  = 2*[x(1)-5  x(2)-2  x(3)-1
        !           x(1)-3  x(2)-4  x(3)-3];
        dfdx(1,:) = 2.0*[x(1,1) - 5.0,x(2,1) - 2.0,x(3,1) - 1.0]
        dfdx(2,:) = 2.0*[x(1,1) - 3.0,x(2,1) - 4.0,x(3,1) - 3.0]    
    end subroutine ObjectiveFunction
end module MMA_variables

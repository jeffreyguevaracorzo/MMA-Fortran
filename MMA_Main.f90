program Example
    use MMA_variables
    use MMA_Routines
    implicit none
    ! --------------------------------------------------------------- !
    ! these subroutines are essentially an adaptation/translation of  !
    ! the code presented by Professor Svanberg originally programmed  !
    ! in Matlab. The definition of the model is presented below and   !
    ! in the subroutine MMA_variables.f90 (avoids modifying other     !
    ! elements of the code).                                          !
    !                                                                 !
    ! Note: for the execution of the code it is necessary to use the  !
    !       Lapack library to execute the SGESV and DGESV routines    !
    !       and additionally a makefile is presented to facilitate    !
    !       the compilation and execution.                            !
    ! --------------------------------------------------------------- !
    m = 2
    n = 3
    allocate(xval(n,1),xmin(n,1),xmax(n,1))
    allocate(a(m,1),c(m,1),d(m,1))
    allocate(df0dx(n,1),fval(m,1),dfdx(m,n))
    epsimin = 0.0000001
    outeriter = 0
    maxoutit = 1
    kkttol = 0
    xval(:,1) = [4.0,3.0,2.0]
    xold1 = xval
    xold2 = xval
    xmin(:,1) = [0.0,0.0,0.0]
    xmax(:,1) = [5.0,5.0,5.0]
    low = xmin
    upp = xmax
    a0 = 1.0
    a(:,1) = [0.0,0.0]
    d(:,1) = [1.0,1.0]
    c(:,1) = [1000.0,1000.0]

    ! --------------------------------------------------------------- !
    !               Here star the MMA optimization process            !
    ! --------------------------------------------------------------- !
    if (outeriter.lt.0.5) then
        call ObjectiveFunction(xval,f0val,df0dx,fval,dfdx)
        outvector1 = [real(outeriter), xval]
        outvector2 = [f0val, fval]
    end if
    kktnorm = kkttol + 10
    outit = 0        

    ! Optimization loop
    do while ((kktnorm.gt.kkttol).and.(outit.lt.maxoutit))
        ! counter
        outit = outit + 1
        outeriter = outeriter + 1
        
        ! MMA subproblem is solved
        call MMA_Sub(xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp,m,n,outeriter,xval, &                                
                     xmin,xmax,xold1,xold2,f0val,df0dx,fval,dfdx,a0,a,c,d)
        
        ! Updating of some vectors
        xold2 = xold1
        xold1 = xval
        xval  = xmma

        ! Getting the f0val, df0dx, fval and dfdx
        call ObjectiveFunction(xval,f0val,df0dx,fval,dfdx)

        ! Residual vector of KKT conditions
        call kktcheck(residu,kktnorm,residumax,xmma,ymma,zmma,lam,xsi,eta, &
                      mu,zet,s,xmin,xmax,df0dx,fval,dfdx,a0,a,c,d)
        outvector1 = [real(outeriter), xval(:,1)]
        outvector2 = [f0val, fval(:,1)]
    end do

    ! Print results
    write(unit=*, fmt=*) 'outvector1'
    write(unit=*, fmt=*) outvector1
    write(unit=*, fmt=*) 'outvector2'
    write(unit=*, fmt=*) outvector2
    write(unit=*, fmt=*) 'xval'
    write(unit=*, fmt=*) xval
end program Example

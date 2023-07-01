! ----------------------------------------------------------------------------
SUBROUTINE klhsc (Kmat, nobs, y, nlam, ulam, &
  & tol, maxit, gam, anlam, npass, jerr, alpmat)
! ----------------------------------------------------------------------------
  IMPLICIT NONE
  ! - - - arg types - - -
  INTEGER :: nobs, nlam, anlam, jerr, maxit, npass (nlam)
  DOUBLE PRECISION :: tol, Kmat (nobs, nobs), y (nobs), ulam (nlam)
  DOUBLE PRECISION :: gam, alpmat (nobs+1, nlam)
  ! - - - local declarations - - -
  INTEGER :: j, l, info, lwork, nlmax
  DOUBLE PRECISION, ALLOCATABLE :: works (:)
  DOUBLE PRECISION :: Umat (nobs, nobs), eigens (nobs)
  DOUBLE PRECISION :: lpinv (nobs), Usum (nobs), lpUsum (nobs)
  DOUBLE PRECISION :: rnobs, gval, hval, gamvec (nobs)
  DOUBLE PRECISION :: zvec (nobs), svec (nobs), vvec (nobs)
  DOUBLE PRECISION :: r (nobs), dif (nobs+1), alpvec (nobs+1)
  ! - - - begin - - -
  nlmax = 8 * nobs
  rnobs = Real(nobs)
  
  Umat = Kmat
  ALLOCATE(works(nlmax))
  lwork = -1
  CALL DSYEV('vectors', 'upper', nobs, Umat, nobs, eigens, works, lwork, info)
  lwork = MIN(nlmax, INT(works(1)))
  CALL DSYEV('vectors', 'upper', nobs, Umat, nobs, eigens, works, lwork, info)
  DEALLOCATE(works)
  eigens = eigens + gam
   
  npass = 0
  r = 0.0D0
  alpmat = 0.0D0
  alpvec = 0.0D0
  zvec = 0.0D0
  svec = 0.0D0
  vvec = 0.0D0
  
  Do j = 1, nobs
     Usum(j) = Sum(Umat(:,j))
  ENDDO
  
  lambda_loop: DO l = 1, nlam
    dif = 0.0D0
    lpinv = 1.0D0 / (eigens + 2.0D0 * rnobs * ulam(l))
    lpUsum = lpinv * Usum
    vvec = Matmul(Umat, eigens * lpUsum) !: g_K
    svec = Matmul(Umat, lpUsum) !: v_K
    gval = 1.0D0 / (rnobs - Sum(vvec))
    ! - - - update alpha - - - 
    update_alpha: DO
        DO j = 1, nobs
          IF (r(j) > 1.0D0) THEN
            zvec(j) = - y(j) * r(j) ** (-1)
          ELSE
            zvec(j) = - y(j)
          END IF
        ENDDO
      gamvec = zvec + 2.0D0 * rnobs * ulam(l) * alpvec(2:(nobs + 1))
      hval = sum(zvec) - Dot_product(vvec, gamvec)
      dif(1) = - gval * hval
      dif(2:(nobs+1)) = -dif(1) * svec - &
        & Matmul(Umat, Matmul(gamvec, Umat) * lpinv)
      alpvec = alpvec + dif
      r = r + y * (dif(1) + Matmul(Kmat, dif(2:(nobs + 1))))
      npass(l) = npass(l) + 1 
      IF (Sum(dif * dif) < tol) EXIT
      IF (Sum(npass) > maxit) EXIT
    ENDDO update_alpha
    alpmat(:, l) = alpvec
    IF (Sum(npass) > maxit) THEN
      jerr = -l
      EXIT
    ENDIF
    anlam = l
  ENDDO lambda_loop
END SUBROUTINE klhsc


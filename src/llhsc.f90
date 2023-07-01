! ----------------------------------------------------------------------------
SUBROUTINE llhsc (Xmat, nobs, np, y, nlam, ulam, &
  & tol, maxit, gam, anlam, npass, jerr, btmat)
! ----------------------------------------------------------------------------
  IMPLICIT NONE
  ! - - - arg types - - -
  INTEGER :: nobs, np, nlam, anlam, jerr, maxit, npass (nlam)
  DOUBLE PRECISION :: tol, Xmat (nobs, np), y (nobs), ulam (nlam)
  DOUBLE PRECISION :: gam, btmat (np+1, nlam)
  ! - - - local declarations - - -
  INTEGER :: j, l, info, npone, lwork, nlmax
  DOUBLE PRECISION, ALLOCATABLE :: works (:)
  DOUBLE PRECISION :: cval, gval, eigens (np+1), pinv (np+1), XX (np, np)
  DOUBLE PRECISION :: Umat (np+1, np+1), Aione (np+1), Xsum (np), gamvec (np+1)
  DOUBLE PRECISION :: rnobs, Gmat (np+1, np+1)
  DOUBLE PRECISION :: r (nobs), zvec (nobs), dif (np+1), btvec (np+1)
  ! - - - begin - - -
  npone = np + 1   !: intercept + number of variables
  nlmax = 8 * npone
  npass = 0
  r = 0.0D0
  btmat = 0.0D0
  btvec = 0.0D0

  rnobs = Real(nobs)
  
  !: column sum
  !Xsum = Sum(Xmat, dim=1) 
  Do j = 1, np
     Xsum(j) = Sum(Xmat(:,j))
  ENDDO

  !: Gmat = P0
  Gmat(1, 1) = rnobs
  Gmat(1, 2:(np + 1)) = Xsum
  Gmat(2:(np + 1), 1) = Xsum
  CALL DGEMM('T', 'N', np, np, nobs, 1.0D0, &
    & Xmat, nobs, Xmat, nobs, 0.0D0, XX, np) !: XX = X^T X
  Gmat(2:(np + 1), 2:(np + 1)) = XX

  !: eigens
  Umat = Gmat
  ALLOCATE(works(nlmax))
  lwork = -1
  CALL DSYEV('vectors', 'upper', npone, Umat, npone, eigens, works, lwork, info)
  lwork = MIN(nlmax, INT(works(1)))
  CALL DSYEV('vectors', 'upper', npone, Umat, npone, eigens, works, lwork, info)
  DEALLOCATE(works)
  eigens = eigens + gam
  
  lambda_loop: DO l = 1, nlam
    dif = 0.0D0
    cval = 2.0D0 * rnobs * ulam(l)
    pinv = (eigens + cval) ** (-1)
    Aione = Matmul(Umat, pinv * Umat(1,:)) !: Aione = v = Q_lamb_inv
    gval = cval / (1.0D0 - cval * Aione(1)) !: gval = g
    ! - - - update beta - - - 
    update_beta: DO
      DO j = 1, nobs
        IF (r(j) > 1.0D0) THEN
          zvec(j) = - y(j) * r(j) ** (-1) 
        ELSE
          zvec(j) = - y(j)
        END IF
      ENDDO
      !: gamvec = n * gam
      gamvec(1) = sum(zvec)
      gamvec(2:npone) = Matmul(zvec, Xmat) + 2 * rnobs * ulam(l) * btvec(2:npone)
      
      gamvec(1) = gamvec(1) + gval * Dot_product(Aione, gamvec)
      dif = - Matmul(Umat, pinv * Matmul(gamvec, Umat))
      btvec = btvec + dif
      r = r + y * (dif(1) + Matmul(Xmat, dif(2:npone)))
      npass(l) = npass(l) + 1
      IF (Sum(dif * dif) < tol) EXIT
      IF (Sum(npass) > maxit) EXIT
    ENDDO update_beta
    btmat(:, l) = btvec
    IF (Sum(npass) > maxit) THEN
      jerr = -l
      EXIT
    ENDIF
    anlam = l
  ENDDO lambda_loop
END SUBROUTINE llhsc

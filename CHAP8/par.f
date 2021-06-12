C  File par.f: parameters for the CCC program  
C  The size of the runfile is determined by the following declarations in the
C  main program. 
C      real vmat(kmax*nchan,kmax*(nchan+1))
C      dimension chil(maxr,kmax,nchan)
C  The following choice results in a file of 20M on a SUN
C     The MAXR parameter determines the maximum number of radial points
      PARAMETER (maxr=1000)
C  The nchan = ncmax parameters determine the maximum number of states
      parameter (nchan=20,ncmax=nchan,nnmax=nchan)
C  The KMAX parameter determines the max number of k-grid points      
      parameter (kmax=51)

      
C  The following parameters should be left alone
      parameter (lamax=0,nchane2e=1,nchanop=1)
      parameter (lnabmax=lamax,ltmax=1)
      parameter (lmax=1)

C  File main.f
C  The following program is the much reduced version of Igor Bray's CCC
C  program for the calculation of electron/positron scattering on atoms/ions.
C  This version is specific to the e-H Temkin-Poet model, where only
C  states with zero orbital angular momentum are treated.
      program CCC_TP
      include 'par.f'
      real vmat(kmax*nchan,kmax*(nchan+1))
      complex wk(kmax*nchan)
      real ovlp(ncmax,0:lamax)
      common/meshrr/ meshr,rmesh(maxr,3)
      common/cont/ ce(ncmax),cint(ncmax),ncstates,energy
      common /pspace/ nabot,labot,natop,latop,ntype,ipar,nze,ninc,linc
      dimension nabot(0:lamax), natop(0:lamax)
      dimension npk(nchan+1)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      common/smallr/ formcut,regcut,expcut,fast
      logical fast
      common /double/id,jdouble(20)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
      dimension chil(maxr,kmax,nchan),minchil(kmax,nchan)
      dimension nk(5,0:lmax), sk(5,0:lmax), 
     >   alpha(0:lamax), nps(0:lamax), nbnd(0:lmax)
      complex phasel(kmax,nchan),sigma(nchan)
      dimension ui(maxr), weightk(kmax)
      dimension gk(kmax,nchan),noprint(nchan),temp(maxr)
      character chan(100)*3, target*6, projectile*8
      common /charchan/ chan
      common /chanens/ enchan(100)
      logical noprint,hlike,uba(nchan)
      data uba/nchan*.false./
C
C  the following units correspond to files used by the program
C  3 ccc.in: the input file
C  42 totalcs: file of total cross sections
C  42 singlet.xx, triplet.xx: files storing half-on-shell K and V matrix
C                             elements
      ry= 13.6058      
      pi= acos(-1.)
      anorm= 8./pi
      ntype = 0
      
C  read the input. Many of the following parameters are not used in this
C  code.
      call readin(labot,latop,nabot,natop,lnabtop,nnbtop,energy,
     >   ntstart,ntstop,lstart,lstop,lttop,nbnd,npsbnd,albnd,alpha,
     >   ncstates,npot,lpot,nptop,lptop,formcut,regcut,expcut,theta,
     >   ifirst,isecond,nold,nqm,qcut,rmax,fast,ldw,nk,sk,nznuc,nze,
     >   itail,gamma,r0,ninc,linc,ipar,nent,zasym,nunit,ndbl,nps,
     >   ne2e,lslow,lfast,slowe,enion,enlevel,target,projectile)

C Define the projectile energy ery in Rydbergs
      ery = energy / ry
      
C  Setup the FACLOG array
      call arrange
      
C  Setup the R integration grids
      call grids(qcut, ndbl, rmax, rmesh, maxr, meshr, jdouble, 20, id)

C  Define arrays containing the powers of R, as well as 
C  the CNTFUG array which is L * (L + 1) / R**2.
      call setpow(0,0)

      hlike = .true.
      
C  Define the target states, first the eigenstates then the pseudostates
      nchanmax = 0
      nstmax = 0
      do la = 0, latop
         l = la
         print'('' Defining eigenstates to nnbtop:'',i3)', nnbtop
         call makeeigenpsi(nznuc,zasym,nnbtop,l,
     >      ery,ry,enion,enlevel)
         if (nold.gt.0) then
            al = alpha(l)
            if (nps(l).gt.0) then
               npstat = nps(l) - l
            else
               npstat = natop(l) - l
            endif 
            print'('' Defining pseudostates with NPS:'',i3,
     >         '' and ALPHA = LAMBDA / 2:'',f6.3)', npstat, al
             call makepspsi(nznuc,npstat,nnbtop,l,
     >         al,ovlp,ery,ry,enion,enlevel,slowe)
         endif 
         print*
         if (ipar.eq.0) then
            nchanmax = nchanmax + (natop(l) - nabot(l) + 1) * (l + 1)
            nstmax = nstmax + natop(l) - nabot(l) + 1
         else
            nchanmax = nchanmax + (natop(l) - nabot(l) + 1) * l
            if (l.gt.0) nstmax = nstmax + natop(l) - nabot(l) + 1
         endif 
         if (nchanmax.gt.nchan.and.l.le.lstop) then
            print*,'NCHAN should be at least NCHANMAX',nchan,nchanmax
            stop 'NCHAN should be at least NCHANMAX'
         endif
      end do
      print'('' Number of states and max number of channels:'',2i3)',
     >   nstmax,nchanmax
      
C  The valence electron is assumed to be initially in the state
C  NINC, LINC.
      etot = ery + enpsinb(ninc,linc)
      
      print'('' Total energy of the collision system:'',
     >   f9.4,'' eV ('',f8.4,'' Ryd )'')', ry * etot, etot
      n = 1
      ipar = 0
      nsmax = 0
      do while(n.le.nent.and.n.ne.0)
         call getchinfo (n,nchp,0,temp,maxpsi,enpsi,la,na,lp)
         if (n.ne.0) then
            print'('' '',a8,'' scattering on '',a6,a3,'' at  '',f9.4,
     >         '' eV ('',f8.4,'' Ryd )'')', projectile,target,
     >         chan(n),ry * (etot - enpsi), etot - enpsi
            if (chan(n)(1:1).eq.' '.or.chan(n)(1:1).eq.'t') nsmax = 1
            n = n + 1
         endif 
      enddo
C  For positrons the spin index will range from 0 to 0 only
      if (nze.eq.1) nsmax = 0
      call pwrite(nopen,energy,nznuc,zasym,ry,noprint,ovlp,lstart,
     >   lstop,projectile,target)

      do i = 1, maxr
         ui(i) = 0.0
      end do
      nqmold = nqm
C  Start the partial wave LG (=J total orbital angular momentum) loop
      do 770 lg = lstart, lstop
C  Define the k-grids
         call kgrid(nk,sk,etot,gk,wk,weightk,nbnd,nqm,lg,nchtop,nchopt,
     >      npk,nold,nptop,lptop,ifirst,ntstop,nnbtop,hlike,uba)

         if ((npk(nchtop+1)-1) * npk(nchtop+1) .gt.
     >      (nchan * kmax) * (nchan * kmax + 1))
     >      stop 'Increase either KMAX or NCHAN in par.f'
         if (meshr * (npk(nchtop+1)-1).gt. maxr * kmax * nchan)
     >      stop 'Increase either MAXR, KMAX, or NCHAN in par.f'

C  Define the plane waves
         call makechil(lg,gk,npk,minchil,chil,npk(nchtop+1)-1,phasel,
     >      nchtop,sigma)
         
C  Initialise the first order V matrix.
         call initv(vmat,npk(nchtop+1)-1)

C  Define the V-matrix elements
         call first(ifirst,nold,etot,lg,gk,npk,chil,minchil,ui,ldw,
     >      phasel,itail,nznuc,nchtop,qcut,vdondum,
     >      vmat,npk(nchtop+1)-1,theta,lfast,
     >      lslow,slowery,td,te1,te2)
         
C  Form and solve the linear equations
         if (ifirst.ge.0)
     >      call solvet(ifirst,vmat,gk,wk,weightk,nchtop,nqm,noprint,
     >      nopen,etot,lg,vdon,phasel,nent,isecond,nunit,sigma,ovlp,
     >      lfast,lslow,slowery,nznuc,zasym,npk,
     >      projectile,target,uba,nsmax)
C  End of LG loop
 770  continue
         
      stop 'job completed'
      end

C  Define hydrogen eigenstates
      subroutine makeeigenpsi(nznuc,zasym,nnbtop,l,
     >   ery,ry,enion,enlevel)
      include 'par.f'
      common/meshrr/ meshr,rmesh(maxr,3)
      common /pspace/ nabot,labot,natop,latop,ntype,ipar,nze,ninc,linc
      dimension nabot(0:lamax), natop(0:lamax)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      common/cont/ ce(ncmax),cint(ncmax),ncstates,energy
      common/smallr/ formcut,regcut,expcut,fast
      logical fast
      character opcl*10

      print'(''  l  n   rn                    '',
     >   ''e(eV)    |<phi|phi>|**2'')'
      rmax = rmesh(meshr,1)

      do n = max(nabot(l),l+1), min(nnbtop, nnmax)
         if (nznuc-nint(zasym).eq.1) then
C  We are here for hydrogen like ions as well as hydrogen
            call rnl(nznuc,n,l,psinb(1,n,l),enpsinb(n,l),
     >         istoppsinb(n,l))
         end if
         sum = 0.0
         do i=1,istoppsinb(n,l)
            sum = sum + psinb(i,n,l)*psinb(i,n,l) * rmesh(i,3)
         end do
         elevel = (enpsinb(n,l) * ry + enion) * 
     >      enlevel / enion
         if (ery+enpsinb(ninc,linc)-enpsinb(n,l).gt.0.0)
     >      then
            opcl = 'open'
         else
            opcl = 'closed'
         endif 
         print'(2i3,i6,f12.1,2f14.8,2x,a6)',
     >      l,n,istoppsinb(n,l),elevel,enpsinb(n,l)*ry,
     >      sum,opcl
         if (abs(sum-1.0).gt.1e-4.and.n.eq.nabot(l))
     >      stop 'Problem with wavefunction'
      end do
      return
      end
      
C  Define pseudostates
      subroutine makepspsi(nznuc,nps,nnbtop,l,
     >   al,ovlp,ery,ry,enion,enlevel,slowe)
      include 'par.f'
      common/meshrr/ meshr,rmesh(maxr,3)
      common /pspace/ nabot,labot,natop,latop,ntype,ipar,nze,ninc,linc
      dimension nabot(0:lamax), natop(0:lamax)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      common/cont/ ce(ncmax),cint(ncmax),ncstates,energy
      common/smallr/ formcut,regcut,expcut,fast
      logical fast
      common /worksp/
     >   ps2(maxr,ncmax),psen2(ncmax),minps2(ncmax),maxps2(ncmax)
      dimension als(ncmax),
     >   vnucl(maxr),ovlp(ncmax,0:lamax)
      character opcl*10,chin*1
      data vnucl,maxpot/maxr*0.0,0/
      save etot
      
      z0 = float(nznuc) 
      alorig = al
      al = abs(al)
      do n = 1, nps
         als(n) = al * 2.0
      enddo
      call makeps(z0, al, l, expcut, nps, ps2,
     >   psen2, minps2, maxps2, rmesh, vnucl, maxpot,
     >   meshr)
      if (l.eq.linc) then
         etot = ery + psen2(ninc-l)
         if (slowe.lt.0.0.and.etot.gt.0.0) slowe = etot / 2.0 * ry
         slowery = slowe / ry
      endif
 
C  Some iteration of alpha to get one of the states have energy slowery
      if (slowery.ne.0.0.and.etot.gt.0.0.and.psen2(nps).gt.etot) then
         almin = al
         almax = al
         small = 1e10
         diffmin = -1e10
         diffmax = 1e10
         it = 0
         alstep = al/50.0
         do while (small.gt.1e-5.and.it.lt.30)
            it = it + 1
            call makeps(z0, al, l, expcut, nps, ps2,
     >         psen2, minps2, maxps2, rmesh, vnucl,
     >         maxpot, meshr)
            small = 1e10               
            do n = 1, nps
               diff = (slowery-psen2(n)) / slowery
               if (abs(diff).lt.small) then
                  small = abs(diff)
                  nsmall = n
c$$$  print*, n, diff, almin, al, almax, psen2(n) * ry
               endif
            enddo
            diff = (slowery-psen2(nsmall)) / slowery
            if (diff.lt.0.0) then
               almax = al
               diffmax = diff
               if (diffmin.gt.0.0) then
c$$$  al = (almax + almin) / 2.0
                  al = (almax * diffmin - almin * diffmax) /
     >               (diffmin - diffmax)
               else
                  al = al - alstep
                  almin = al
               endif 
            else
               almin = al
               diffmin = diff
               if (diffmax.lt.0.0) then
                  al = (almax * diffmin - almin * diffmax) /
     >               (diffmin - diffmax)
c$$$  al = (almax + almin) / 2.0
               else
                  al = al + alstep
                  almax = al
               endif 
            endif
         enddo                   
c$$$            print*, it, nsmall, almin, al, almax, psen2(nsmall)*ry,
c$$$     >         slowery * ry
         print'('' ALPHA redefined to be:'',f10.6,'' after'',i3,
     >      '' iterations'')', al, it
      endif 
      print'(''  l  n   rn                    '',
     >   ''e(eV)    |<phi|phi>|**2  projection'')'

      ntop = min(nps + l,nnmax)
C  Define the projections of the pseudostate onto the discrete eigenstates
      do n = nabot(l), ntop
         ovlp(n,l) = 0.0
         do ne = nabot(l), nnbtop
            if (psen2(n-l).lt.0.0) then
               sum = 0.0
               do i = 1, min(istoppsinb(ne,l),maxps2(n-l))
                  sum = sum + psinb(i,ne,l) * ps2(i,n-l) * rmesh(i,3)
               enddo
               ovlp(n,l) = ovlp(n,l) + sum * sum
            endif 
         enddo
         do ne = nnbtop + 1, nnbtop + ncstates
               sum = 0.0
               do i = 1, min(istoppsinb(ne,l),maxps2(n-l))
                  sum = sum + psinb(i,ne,l) * ps2(i,n-l) * rmesh(i,3)
               enddo
               ovlp(n,l) = ovlp(n,l) + sum * sum * cint(ne-nnbtop)
          enddo 
      enddo

C  Redefine the eigenstate array with the pseudostates
      do n = nabot(l), ntop
         enpsinb(n,l) = psen2(n-l)
         istoppsinb(n,l) = maxps2(n-l)
         sum = 0.0
         do i=1,istoppsinb(n,l)
            psinb(i,n,l)=ps2(i,n-l)
            sum = sum + psinb(i,n,l)*psinb(i,n,l) * rmesh(i,3)
         end do
         do i = istoppsinb(n,l) + 1, maxr
            psinb(i,n,l) = 0.0
         enddo 
         en = enpsinb(n,l)
C  sqrt(2/pi)
         const = 0.79788456
         elevel = (enpsinb(n,l) * ry + enion) * 
     >      enlevel / enion
         if (ery+enpsinb(ninc,linc)-enpsinb(n,l).gt.0.0) then
            opcl = 'open'
         else
            opcl = 'closed'
         endif
C  If natop(l) is -N, then open channels plus following N - 1 closed ones
C  will be included.
         if (natop(l).lt.0.and.opcl.eq.'closed') 
     >      natop(l) = min(ntop, n - natop(l) - 2)
         if (n.le.natop(l).or.natop(l).le.0) then
            chin = '+'
         else
            chin = '-'
         endif 
         print'(2i3,i6,f12.1,3f14.8,2x,a6,1x,a)',
     >      l,n,istoppsinb(n,l),elevel,enpsinb(n,l)*ry,
     >      sum,ovlp(n,l),opcl,chin
      enddo 
C  If natop(l) is 0 then all open and closed channels are included. The <=
C  is used in case all channels are open, but input latop < 0.
      if (natop(l).le.0) natop(l) = ntop
      return
      end

      subroutine initv(vmat,n)
      real vmat(n,n+1)
      do ki = 1, n
         do kf = ki, n
            vmat(kf,ki) = 0.0
            vmat(ki,kf+1) = 0.0
         end do
      end do
      return
      end
      
      
C  The following program diagonalises a one-electron Hamiltonian with a given 
C  local potential using a Laguerre basis. It was written by Dmitry Konovalov
C  and first used in the P.R.A 43, 1301 (1991) paper. Igor Bray's modifications
C  for non-local potential have been removed.
c======================================================================
c  The following routine diagonalizes the Hamiltonian
C             2
C            d         l (l + 1)      Zinf
C     H = - -----   +  ---------  -  ------  + pot0(r)
c               2          2            r
c            d r          r
C
C  in a Laguerre basis of size NPS and exponential fall off factor ALF.
C
C                          2l+2           l+1
C   Laguerre basis:       L   (2 alf r)  r    exp (- alf r)
C                          n-1
C
C  To get the first bound state for each l exactly set alf to be
C  Z of atom divided by (l+1) and take NPS sufficiently large.
C  This program has been tried with NPS up to 100.
C===================================================================
C                         i    i   i    i      i    o    o       o
      subroutine makeps(Zinf, alf, l, expcut, nps, ps, energy, jmin, 
     >   jmax, gridr, pot0, maxpot, nr)
C         o      i     i      i    i    i
C  EXPCUT is used to avoid underflows at the start and end of the wave
C         functions PS.
C  ENERGY is the eigenvalues array in Rydbergs
C  GRIDR  is the grid containing NR R points and integration weights.
C  POT0(i) is defined for i = 1, maxpot is the input potential in Rydbergs
      include 'par.f'
      parameter (ndim=ncmax,nld=maxr)
      implicit double precision (a-h, o-z)
c
      real ps(maxr,ndim), gridr(maxr,3), pot0(maxr),energy(ndim),
     >   expcut, Zinf, alf
c
      dimension c(ndim, ndim), h(ndim, ndim),ch(ndim,ndim)
      dimension enrgr(ndim), dnorm(ndim), work(ndim), 
     >   dpot(nld)
      dimension grid(nld,3), fac(0:1000)
      dimension f(nld, ndim), f1(nld, ndim), jmin(ndim), jmax(ndim)
c

      dimension potmtrx(ndim, ndim), pl(ndim), weight(nld)
      dimension dnorm2(ndim), work2(8*ndim)

C
C     Parameter checking
C
      if (nr.gt.nld) then
         print*,'Increase NLD to at least:',nr
         stop 'Increase NLD'
      end if 
      if (nps.gt.ndim) then
         print*,'Increase NDIM to at least:',nps
         stop 'Increase NDIM'
      end if 
c
c     Define double precision input arrays
c
      r0 = 0d0
      i = 1
      m = 0
      dr = gridr(i,1)
      do while (i.le.nr)
         grid(i,1) = r0 + dfloat(i-m) * dr
         grid(i,3) = dfloat(2 * mod(i,2) + 2) * dr / 3d0
         if (abs((grid(i,1)-gridr(i,1))/grid(i,1)).gt.1e-5) then
            grid(i,1) = grid(i,1) + dr
            r0 = grid(i,1)
            dr = 2d0 * dr
            grid(i-1,3) = grid(i-1,3) * 1.5d0
            grid(i,3) = dfloat(2 * mod(i,2) + 2) * dr / 3d0
            m = i
         endif
         i = i + 1
      enddo
      grid(nr,3) = grid(nr,3) / 2d0
            
      test = 0d0
      do i=1,nr
C  DPOT is in a.u.
         dpot(i) = dble(pot0(i)) / 2.0
         weight(i) = grid(i,3)
      end do
      dalf = dble(alf)
      dz = dble(Zinf)
      dexpcut = dble(expcut)
c
c     Define Matrix of the potential
c
      call potl (dpot, maxpot, L, dalf, potmtrx, nps, ndim, 
     >   f, pl, grid, weight, nr, nld)
c
      nfac = 2 * nps + 2 * L
C  Diagonalize the Hamiltonian
      call basisl (potmtrx, dz, L, dalf, nps, ndim, enrgr, c, h,
     >   dnorm, dnorm2, work2, fac, nfac)
C  Make the PS wave functions on GRID
      call psdlgr (L, dalf, nps, ndim, c, f, f1, dnorm, enrgr, work,
     >   grid, nr, nld, dexpcut, jmin, jmax)
C  Check that the the wave functions are orthonormal. Useful only
C  if the PS states are sufficiently small at last R.
      call ortogc (f, nld, nr, nps, ch, ndim, grid(1,3), jmin, jmax)

C  Define the output arrays
      do n = 1, nps
         do i = 1, maxr
            ps(i,n) = 0.0
         enddo 
         if (f(jmin(n),n).gt.0) then
            do i = jmin(n), jmax(n)
               ps(i,n) = real(f(i,n))
            end do
         else
            do i = jmin(n), jmax(n)
               ps(i,n) = - real(f(i,n))
            end do
         end if 
C  Energy is in Rydbergs
         energy(n) = 2.0 * real(enrgr(n))
      end do 
      end

C=====================================================================
C  Define pseudostates in Laguerre basis using expansion coefficients C.
C=====================================================================
      subroutine psdlgr (L, alf, nmax, nload, c, f, f1, dnorm,enrgr,pl,
     >   grid, nr, nld, epscut, jmin, jmax)
      include 'par.f'
      implicit double precision (a-h, o-z)
      dimension c(nload, nmax), f(nld, nmax), dnorm(nmax), pl(nmax),
     >   grid(nr), jmin(nmax), jmax(nmax), f1(nld, nmax), enrgr(nmax)
      common /flogs/ faclog(1000)
      dimension powr(ncmax,ncmax,0:lamax), en(ncmax,0:lamax),
     >   alpha(0:lamax), nps(0:lamax)
      common /pspace/ nabot,labot,natop,latop,ntype,ipar,nze,ninc,linc
      common /worksp/
     >   ps2(maxr,ncmax),psen2(ncmax),minps2(ncmax),maxps2(ncmax)
      dimension nabot(0:lamax), natop(0:lamax)
      real ps2, psen2
C
C INPUT:
C -----      
C  L - orbital momentum.
C  alf - parameter of basis.
C  nmax - number of pseudostates.
C  nload, nld -  dimensions of arrays
C  C    - matrix of eigenvectors.
C  grid - "R" - grid of "NR" points
C
C OUTPUT:
C ------
C  f(i,n) - "n"-radial wave functions.
C  f1(i,n) = f(i,n) * r**(-L-1) * exp(alf*r) for laguerre integration.
C     jmin(n) - the first point of "n"-state.
C     jmax(n) - the last  point of "n"-state.
C  Pl - work space for the program.
C
      do n = 1, nmax
         en(n,l) = enrgr(n)
      enddo 
      nps(l) = nmax
      alpha(l) = alf
      L2 = 2 * L
      alfl = alf 
      a2 = 2.0d0 * alfl
C
C     Do some checks
C
      if (nr .gt. nld) then
         print*, 'number of points is more than dimension of GRID'
         print*, 'NR=', nr, ' NLD=', nld
         stop    'number of points is more than dimension of GRID'
      end if 
C
C     Define JMIN(N)
C
      xfirst = exp( log(epscut/dnorm(1)) / dble(l+1))
      rfirst = xfirst / a2
      i = 1
      do while ((grid(i) .le. rfirst) .and. (i .lt. nr))
         i = i + 1
      end do 
      ifirst = i
      do n = 1, nmax
         jmin(n) = ifirst
         jmax(n) = nr
      end do 
C  Get the Slater type orbital coefficients      
      do n = 1, nmax
         do m = 1, nmax
            powr(m,n,l) = 0d0
         enddo
      enddo

      do n = 1, nmax
         do m = 1, nmax
            const = c(m,n) * dnorm(m)
            do k = 0, m - 1
               powr(k+1,n,l) = powr(k+1,n,l) + const * (-a2)**k *
     >            exp(faclog(m+2+l2) - faclog(k+1) - faclog(m-k) -
     >            faclog(3+l2+k))
            enddo
         enddo
      enddo 
C
C     Loop by  r-grid
C
      cc1 = a2**(L+1)
      do i = 1, nr
         r = grid(i)
         x = a2 * r
         c1 = exp(-0.5d0 * x  +  dble(L+1) * log(x))
C
C        Define Laguerre polynomials
C
         pl1   = 1.0d0
         pl(1) = pl1
         pl2   = dble(L2 + 3) - x
         if (nmax.gt.1) pl(2) = pl2
         do n = 3, nmax
            pl3 = ((dble(2*n-1+L2)-x)*pl2 - dble(n+L2)*pl1) /
     >         dble(n-1)
            pl(n) = pl3
            pl1 = pl2
            pl2 = pl3
         end do
C
C        Loop by number of states.
C
         do n = 1, nmax
            f(i,n) = 0.0d0
            if (i .le. jmax(n)) then
               sum = 0.0d0
               do m = 1, nmax
                  sum = sum + c(m, n) * dnorm(m) * pl(m)
               end do
               f(i,  n) = c1  * sum
               f1(i, n) = cc1 * sum
C  The commented code below evaluates the wave function using a Slater 
C  representation and compares with the Laguerre form
c$$$               sum = 0d0
c$$$               do m = 1, nmax
c$$$                  sum = sum + powr(m,n,l) * r ** (m+l)
c$$$               enddo
c$$$               test = sum * exp(-alf*r) * a2 ** (l+1)
c$$$               if (abs((test-f(i,n))/(test+f(i,n))).gt.1e-3) print*,
c$$$     >            l,n,i,test,f(i,n)
               if (i .gt. jmin(n)+3+100) then
                  if (     (abs(f(i-3, n))  .lt.  epscut)
     >               .and. (abs(f(i-2, n))  .lt.  epscut)
     >               .and. (abs(f(i-1, n))  .lt.  epscut)
     >               .and. (abs(f(i,   n))  .lt.  epscut))
     >               jmax(n) = i - 3
               end if 
            end if 
         end do 
      end do 
      do n = 1, nmax
         if (f(jmin(n),n).lt.0d0) then
            do m = 1, nmax
               powr(m,n,l) = - powr(m,n,l)
            enddo 
         endif
      enddo
      return
      end

      subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)
c
      integer n,nm,ierr,matz
      double precision a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real symmetric matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a.
c
c        a  contains the real symmetric matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        z  contains the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tqlrat
c           and tql2.  the normal completion code is zero.
c
c        fv1  and  fv2  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  tred1(nm,n,a,w,fv1,fv2)
*  tqlrat encounters catastrophic underflow on the Vax
*     call  tqlrat(n,w,fv2,ierr)
      call  tql1(n,w,fv1,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  tred2(nm,n,a,w,fv1,z)
      call  tql2(nm,n,w,fv1,z,ierr)
   50 return
      end
      subroutine tql1(n,d,e,ierr)
c
      integer i,j,l,m,n,ii,l1,l2,mml,ierr
      double precision d(n),e(n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c
c     this subroutine is a translation of the algol procedure tql1,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues of a symmetric
c     tridiagonal matrix by the ql method.
c
c     on input
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct and
c          ordered for indices 1,2,...ierr-1, but may not be
c          the smallest eigenvalues.
c
c        e has been destroyed.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  sqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
c
      do 290 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 210
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  210    p = d(l) + f
c     .......... order eigenvalues ..........
         if (l .eq. 1) go to 250
c     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue
c
  250    i = 1
  270    d(i) = p
  290 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
      subroutine tql2(nm,n,d,e,z,ierr)
c
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c
c     this subroutine is a translation of the algol procedure tql2,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c
c        e has been destroyed.
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  sqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
c
      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
c
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
c
  300 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
      subroutine tred1(nm,n,a,d,e,e2)
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),e2(n)
      double precision f,g,h,scale
c
c     this subroutine is a translation of the algol procedure tred1,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix
c     to a symmetric tridiagonal matrix using
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        a contains information about the orthogonal trans-
c          formations used in the reduction in its strict lower
c          triangle.  the full upper triangle of a is unaltered.
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
         d(i) = a(n,i)
         a(n,i) = a(i,i)
  100 continue
c     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))
c
         if (scale .ne. 0.0d0) go to 140
c
         do 125 j = 1, l
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = 0.0d0
  125    continue
c
  130    e(i) = 0.0d0
         e2(i) = 0.0d0
         go to 300
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         e2(i) = scale * scale * h
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
         if (l .eq. 1) go to 285
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0
c
         do 240 j = 1, l
            f = d(j)
            g = e(j) + a(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + a(k,j) * d(k)
               e(k) = e(k) + a(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0d0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         h = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - h * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       a(k,j) = a(k,j) - f * e(k) - g * d(k)
c
  280    continue
c
  285    do 290 j = 1, l
            f = d(j)
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = f * scale
  290    continue
c
  300 continue
c
      return
      end
      subroutine tred2(nm,n,a,d,e,z)
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),z(nm,n)
      double precision f,g,h,hh,scale
c
c     this subroutine is a translation of the algol procedure tred2,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix to a
c     symmetric tridiagonal matrix using and accumulating
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        z contains the orthogonal transformation matrix
c          produced in the reduction.
c
c        a and z may coincide.  if distinct, a is unaltered.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
c
         do 80 j = i, n
   80    z(j,i) = a(j,i)
c
         d(i) = a(n,i)
  100 continue
c
      if (n .eq. 1) go to 510
c     .......... for i=n step -1 until 2 do -- ..........
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 2) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))
c
         if (scale .ne. 0.0d0) go to 140
  130    e(i) = d(l)
c
         do 135 j = 1, l
            d(j) = z(l,j)
            z(i,j) = 0.0d0
            z(j,i) = 0.0d0
  135    continue
c
         go to 290
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0
c
         do 240 j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + z(k,j) * d(k)
               e(k) = e(k) + z(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0d0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         hh = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)
c
            d(j) = z(l,j)
            z(i,j) = 0.0d0
  280    continue
c
  290    d(i) = h
  300 continue
c     .......... accumulation of transformation matrices ..........
      do 500 i = 2, n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = 1.0d0
         h = d(i)
         if (h .eq. 0.0d0) go to 380
c
         do 330 k = 1, l
  330    d(k) = z(k,i) / h
c
         do 360 j = 1, l
            g = 0.0d0
c
            do 340 k = 1, l
  340       g = g + z(k,i) * z(k,j)
c
            do 360 k = 1, l
               z(k,j) = z(k,j) - g * d(k)
  360    continue
c
  380    do 400 k = 1, l
  400    z(k,i) = 0.0d0
c
  500 continue
c
  510 do 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = 0.0d0
  520 continue
c
      z(n,n) = 1.0d0
      e(1) = 0.0d0
      return
      end
      double precision function pythag(a,b)
      double precision a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
C============================================================      
C  ch(i,j) = < f(i) I f(j) > - delta(i,j)
C============================================================      
      subroutine ortogc (f, nload, nr, nmax, ch, nch, weight,
     >   jmin, jmax)
      implicit double precision (a-h, o-z)
      dimension f(nload,nmax), ch(nch,nmax), weight(nr),
     >   jmin(nmax), jmax(nmax)
      ncount = 0
      err = 0d0
      do m = 1, nmax 
         do n = 1, nmax
            sum = 0d0
            r = 0d0
            func = 0d0
            do i = max(jmin(m), jmin(n)), min(jmax(m), jmax(n))
               w = weight(i)
               func1 = f(i,m)
               func2 = f(i,n)
               sum = sum + func1 * func2 * w
            end do 
            diag = 0.0d0
            if (n .eq. m)  diag = 1.0d0
            ch(n,m) = abs(sum - diag)
c$$$            print*,n,m,ch(n,m)
            if (ch(n,m).gt.err) err = ch(n,m)
c$$$            err = err + ch(n,m)
         end do 
         do n = 1, nmax
            if (abs(ch(n,m)).gt.1e-3) ncount = ncount + 1
c$$$     >         print*,'Warning: two PS states are not orthonormal'//
c$$$     >         ' or not zero at last R',
c$$$     >         n,m,ch(n,m)
         end do 
      end do
      if (ncount.gt.0) then
         print*,'Number of non orthonormal (to 1e-3) states:', ncount
         print*,'Max orthogonality error:',err
      endif 
      return
      end

C=====================================================================
C        Define matrix of some potential in a Laguerre basis.
C        Normalization done later.
C=====================================================================
      subroutine potl (pot0, maxpot, L, alf, potmtrx, nmax, 
     >    nload, f, pl, grid, weight, nr, nld)
      implicit double precision (a-h, o-z)
      dimension f(nld, nload), grid(nld), pot0(nld)
      dimension potmtrx(nload, nload), pl(nload), weight(nld)
C
C INPUT:
C -----      
C  L     - orbital momentum.
C  alf   - parameter of basis.                     alfl = alf 
C  nmax  - number of basis functions.              ------------------
C  nld   - dimension of r-arrays.
C  nload - dimension of matrices.
C  f, pl - work arrays
C OUTPUT:
C ------
C  potmtrx - matrix of the potential. 
C
      L2 = 2 * L
      L22 = L2 + 2
      alfl = alf 
      a2 = 2.0d0 * alfl
C
C     Do some checks
C
      if (nr .gt. nld) then
         print*, 'number of points is more than dimension of GRID'
         print*, 'NR=', nr, ' NLD=', nld
         print*, 'Stop in potl'
         stop    'number of points is more than dimension of GRID'
      end if 
c
c     Define Laguerre polynomials on r-grid
c
      call lagpol8 (L22, a2, f, nld, nmax, grid, nr, pl)
      do n = 1, nmax
         do i = 1, nr
            r = grid(i)
            f(i,n) = exp(- alfl * r + dble(l + 1) * log(r)) * f(i,n)
         end do
      end do 
c
c     Define Matrix of the potential
c
      do n = 1, nmax
         do m = 1, n
            call int3_8 (f(1,n), 1, nr, f(1,m), 1, nr, pot0, 1, maxpot,
     >         res, weight, nr)
            potmtrx(n,m) = res * a2**L22
            potmtrx(m,n) = potmtrx(n,m)
         end do
      end do 
      return
      end


C=====================================================================
C          Define basis functions in local potential -Z0/r + potmtrx(m,n)
C=====================================================================
      subroutine basisl (potmtrx, Z0, L, alf, nmax, nload, enrg, c, h,
     >   dnorm, dnorm2, work, fac, nfac)
      implicit double precision (a-h, o-z)
      dimension c(nload, nmax), h(nload, nmax), fac(0:nfac),
     >   potmtrx(nload, nload)
      dimension enrg(nmax), dnorm(nmax), dnorm2(nmax), work(8*nmax)
C
C      f(n,l,r) = dnorm(n,l) * (2*alfl*r)**(l+1) * exp(-alfl*r) *
C                 * Laguerre(2*l+2;n-1;2*alfl*r)
C
C INPUT:
C -----      
C  potmtrx   - matrix of some potential. 
C  L     - orbital momentum.
C  alf   - parameter of basis.              alfl = alf 
C  nmax  - number of pseudostates.          ------------------
C  nload - dimension of arrays
c  Z0    - The charge of the coulomb potential.
C
C OUTPUT:
C ------
C  enrg - eigenvalues
C  H    - Hamiltonian matrix.
C  C    - matrix of eigenvectors.
C  Work, dnorm2 - work space for the program
C
      L2 = 2 * L
      alfl = alf 
      a2 = 2.0d0 * alfl
c
c     Do some checks
c
      if (Z0 .lt. 0d0) then
         print*, 'Z0 is out of range,  Z0=', Z0
         stop    'Stop in BASISL'
      end if 
C
C     Define factorials as n!=exp(fac(n))
C
      fac(0) = 0.0d0
      fact = 0.0d0
      tmp = 1.0d0
      ntmp = 2 * nmax + L2
      if (nfac .lt. ntmp) then
         print*, 'nfac is not enough to store array of factorials,'
         print*, 'nfac has to be more than   2 * NMAX + 2 * L =', ntmp,
     >      ' but  nfac =', nfac
         stop 'nfac has to be more than   2 * NMAX + 2 * L'
      end if 
      do n = 1, ntmp
         fact   = fact + log(dble(n))
         fac(n) = fact
      end do
C
C     Define normalization coeff. Dnorm 
C
      c1 = sqrt(a2)
      c2 = 1.0d0 / a2
      do i = 1, nmax
         dnorm(i)  = c1 * exp(0.5d0 * (fac(i - 1)  -  fac(L2 + 1 + i)))
         dnorm2(i) = c2 * exp(fac(L2 + i)  -  fac(i - 1))
      end do 
C
C     Define Hamiltonian matrix 
C
      c2 = -a2 * a2 * 0.5d0
      do i = 1, nmax
         do j = 1, i
            sm2 = 0.0d0
            do jj = 1, j - 1
               do jjj = 1, min(i, jj)
                  sm2 = sm2 + dnorm2(jjj)
               end do
            end do
            sm1 = 0.0d0
            do jj = 1, min(i, j)
               sm1 = sm1 + dnorm2(jj)
            end do 
            diag = 0.d0
            if (i .eq. j)  diag = 0.25d0
            potlz = c2 * Z0 / alfl * sm1
            potl2  =  potmtrx(i,j)
            res = (dnorm(i) * dnorm(j) * (-dble(l+j) * sm1
     >         + sm2) + diag) * c2 
     >         + dnorm(i) * dnorm(j) * (potl2 + potlz)
            h(i, j) = res
            h(j, i) = res
         end do
      end do 
C         
C  if matc = 0, then  only eigenvalues.
      matc = 1
      call rs(nload, nmax, h, enrg, matc, c, dnorm2, work, ierr)
      lwork = 8 * nmax
      if (ierr .ne. 0)  then
         print*, 'Program "RS" finished abnormaly, ierr =', ierr
         stop    'Program "RS" finished abnormaly'
      end if 
      return
      end
C===================================================================
C                                     m
C   Laguerre's polinomials  f(i,n) = L   (dlambda * grid(i))
C                                     n-1
C===================================================================
      subroutine lagpol8 (m, dlambda, f, nload, nmax, grid, nr, pl)
      implicit double precision (a-h, o-z)
      dimension f(nload, nmax), grid(nr),  pl(nmax)
C
C INPUT:
C -----
C  m - parameter of Laguerre polinomial.
C  nmax  - the max number of n.
C  nload - dimension of f.
C  grid  - "R" - grid of "NR" points.
C  pl(nmax) - work space for this program.
C
C OUTPUT:
C ------
C  f(i,n) - Laguerre polinomials.
C
      L2 = m - 2
C
C     Loop by  r-grid
C
      do i = 1, nr
         r = grid(i)
         x = dlambda * r
C
C        Define Laguerre's polinomials
C        and store them in array 'f'
C
         pl1    = 1.0d0
         pl(1)  = pl1
         f(i,1) = pl1
c
         pl2    = dble(L2 + 3) - x
         if (nmax.gt.1) then
            pl(2)  = pl2
            f(i,2) = pl2
         endif 
c
         do n = 3, nmax
            pl3 = ((dble(2*n-1+L2)-x)*pl2 - dble(n+L2)*pl1) /
     >         dble(n-1)
            pl(n) = pl3
            f(i,n) = pl3
            pl1 = pl2
            pl2 = pl3
         end do
      end do 
      return
      end
C
C==================================================================
C
C   int3_8    Program for integration of three double precision functions.
C
C      call int3_8 (f1, nf1, nl1, f2, nf2, nl2, f3, nf3, nl3,
C     >   res, weight, nr)
C
C==================================================================
C   double precision  Programe for integration of three functions
C==================================================================
      subroutine int3_8 (f1, nf1, nl1, f2, nf2, nl2, f3, nf3, nl3,
     >   res, weight, nr)
      implicit double precision (a-h, o-z)
      dimension f1(nr), f2(nr), f3(nr), weight(nr)
      if ((nl1 .gt. nr)  .or.  (nl2 .gt. nr) .or.  (nl3 .gt. nr)) then
         print*, 'Number of points are not enough in r-mesh'
         print*, 'the last point of "F1"=', nl1
         print*, 'the last point of "F2"=', nl2
         print*, 'the last point of "F3"=', nl3
         print*, 'the last point of "WEIGHT"=', nr
         stop 'stop  in  "INT3_8"'
      end if 
      res = 0.0d0
      do i = max(nf1, nf2, nf3), min(nl1, nl2, nl3)
         res = res + f1(i) * f2(i) * f3(i) * weight(i)
      end do 
      return
      end
C  File vmat.f, define V-matrix elements
      subroutine first(ifirst,nold,etot,lg,gk,npk,chil,minchil,u,ldw,
     >   phasel,itail,nznuc,nchtop,qcut,vdon,vmat,ichi,
     >   theta,lfast,lslow,slowery,td,te1,te2)
      include 'par.f'
      integer npk(nchtop+1)
      dimension chil(meshr,ichi),minchil(ichi),
     >   u(maxr),gk(kmax,nchan),vdon(nchan,nchan,0:1),ui(maxr),uf(maxr)
      complex phasel(kmax,nchan)
      real vmat(ichi,ichi+1)
      common/meshrr/ meshr,rmesh(maxr,3)
      common /pspace/ nabot,labot,natop,latop,ntype,ipar,nze,ninc,linc
      dimension nabot(0:lamax), natop(0:lamax)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),maxpsinb(nnmax,0:lnabmax)
      common /double/id,jdouble(20)
      common /radpot/ ucentr(maxr)
      common/smallr/ formcut,regcut,expcut,fast
      logical fast
      dimension psii(maxr), psif(maxr), slowery(ncmax)
      real ud(maxr)
      data pi/3.1415927/
      data ud,uf,ui/maxr*0.0,maxr*0.0,maxr*0.0/

      td = 0.0
      te1 = 0.0
      te2 = 0.0
      rnorm = 2.0/pi
      do nchi = 1, nchtop
         call getchinfo (nchi, nt, lg, psii, maxpsii, ei, lia, nia, li)
         nqmi = npk(nchi+1) - npk(nchi)
C  Start calculation of V matrix elements
         do nchf = nchi, nchtop
            call getchinfo (nchf,nt,lg, psif, maxpsif, ef, lfa, nfa, lf)
            nqmf = npk(nchf+1) - npk(nchf)
            call clock(s1)
            nchtopf = nchtop
            nchtopi = nchtop
C  Define direct elements            
            call makev31d(nqmi,psii,maxpsii,lia,li,
     >         chil(1,npk(nchi)),minchil(npk(nchi)),gk(1,nchi),
     >         phasel(1,nchi),npk(nchtop+1)-1,nqmf,psif,maxpsif,
     >         lfa,lf,chil(1,npk(nchf)),minchil(npk(nchf)),
     >         gk(1,nchf),phasel(1,nchf),npk(nchtop+1)-1,lg,rnorm,
     >         ud,itail,nznuc,vdon,nchf,nchi,.true.,npk,0,vmat)
            call clock(s2)
            td = td + s2 - s1
C  Define exchange terms if IFIRST = 1
            if (ifirst.eq.1) then
C  Define two-electron exchange elements
               call makev3e(psii,maxpsii,lia,nchi,psif,maxpsif,lfa,
     >            nchf,li,lf,chil(1,npk(nchi)),minchil(npk(nchi)),nqmi,
     >            chil(1,npk(nchf)),minchil(npk(nchf)),nqmf,lg,rnorm,
     >            npk,vmat,npk(nchtop+1)-1,nchtop)
               call clock(s3)
               te1 = te1 + s3 - s2
C  Define energy dependent exchange terms
               call makev1e(nqmi,psii,maxpsii,ei,lia,li,
     >            chil(1,npk(nchi)),minchil(npk(nchi)),gk(1,nchi),
     >            npk(nchtop+1)-1,etot,theta,0,nqmf,psif,maxpsif,
     >            ef,lfa,lf,chil(1,npk(nchf)),minchil(npk(nchf)),
     >            gk(1,nchf),npk(nchtop+1)-1,lg,rnorm,
     >            uf,ui,nchf,nchi,nold,nznuc,npk,vmat)
               call clock(s4)
               te2 = te2 + s4 - s3
C  End of exchange
            end if 
C  End of NCHF loop
         end do
C  End of NCHI loop
      end do 
      return
      end

      subroutine getchnl(ch,n,l)
      character ch*3

      n = ichar(ch(2:2)) - ichar('0')
      if (ch(3:3).eq.'S') then
         l = 0
      elseif (ch(3:3).eq.'P') then
         l = 1
      elseif (ch(3:3).eq.'D') then
         l = 2
      elseif (ch(3:3).eq.'F') then
         l = 3
      elseif (ch(3:3).eq.'G') then
         l = 4
      elseif (ch(3:3).eq.'H') then
         l = 5
      elseif (ch(3:3).eq.'I') then
         l = 6
      endif
      return
      end
      
      function getchchar(n,l)
      character ch*1, getchchar*2
      ch(n) = char(n + ichar('0'))

      if (l.eq.0) then
         getchchar = ch(n)//'S'
      elseif (l.eq.1) then
         getchchar = ch(n)//'P'
      elseif (l.eq.2) then
         getchchar = ch(n)//'D'
      elseif (l.eq.3) then
         getchchar = ch(n)//'F'
      elseif (l.eq.4) then
         getchchar = ch(n)//'G'
      elseif (l.eq.5) then
         getchchar = ch(n)//'H'
      elseif (l.eq.6) then
         getchchar = ch(n)//'I'
      endif
      return
      end

C  Defines and returns channel information
      subroutine getchinfo (nch,nchp,lg,psi,maxpsi,enpsi,la,na,lp)
      include 'par.f'
      dimension psi(maxr)
      common /pspace/ nabot,labot,natop,latop,ntype,ipar,nze,ninc,linc
      dimension nabot(0:lamax), natop(0:lamax)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),maxpsinb(nnmax,0:lnabmax)
      character chan(100)*3, getchchar*2, chspin(3)*1
      common /charchan/ chan
      common /chanens/ enchan(100)

C  Do not change the lower case t below as it is used for checking
C  the spin of the incident state in the SOLVET routine
      data chspin/'s',' ','t'/
      
      iparity = (-1)**(lg + ipar)
      ncht = 0
      nchp = 1
      la = linc
      na = ninc
      enchan(nchp) = enpsinb(na,la)
      chan(nchp) = chspin(2)//getchchar(na,la)
      do 10 lp = abs(lg-la), lg+la
         if ((-1)**(la+lp).ne.iparity) go to 10
c$$$      do lp = abs(lg-la), lg+la, 2
         ncht = ncht + 1
         if (nch.eq.ncht) then
            enpsi = enpsinb(na,la)
            maxpsi = maxpsinb(na,la) 
            do i = 1, maxpsi
               psi(i) = psinb(i,na,la)
            end do
            return
         end if 
 10   continue 
c$$$      do la = labot, latop
c$$$         do na = nabot(la), natop(la)
C  The above looping is preferable, but the CROSS program expects it in
C  the order below. 
      nastart = nabot(labot)
      nastop  = natop(labot)
      do la = labot, latop
         if (nabot(la).le.la) then
            print*,'NABOT(LA) must be > LA, here NABOT(LA), LA:',
     >         nabot(la),la
         endif
         if (nabot(la).lt.nastart) nastart = nabot(la)
         if (natop(la).gt.nastop)  nastop  = natop(la)
      enddo
      
      do natmp = nastart, nastop
         do la = labot, min(natmp - 1, latop)
C  The following three lines are there for potassium to stop 3d being 
C  the first channel
            na = natmp
            if (la.le.1.and.nastart.lt.nabot(la).and.la+1.lt.nabot(la))
     >         na = natmp + nabot(la) - nastart
            if (na.ge.nabot(la).and.na.le.natop(la).and.
     >         (la.ne.linc.or.na.ne.ninc)) then
               nchp = nchp + 1
               chan(nchp) = chspin(2)//getchchar(na,la)
               enchan(nchp) = enpsinb(na,la)
               do 20 lp = abs(lg-la), lg+la
                  if ((-1)**(la+lp).ne.iparity) go to 20
c$$$               do lp = abs(lg-la), lg+la, 2
                  ncht = ncht + 1
                  if (nch.eq.ncht) then
                     enpsi = enpsinb(na,la)
                     maxpsi = maxpsinb(na,la) 
                     do i = 1, maxpsi
                        psi(i) = psinb(i,na,la)
                     end do
                     return
                  end if 
 20            continue 
            end if 
         end do
      end do
      nch = 0
      return
      end

C  One-electron exchange terms. Theta dependence is here.
      subroutine makev1e(nqmi,psii,maxpsii,ei,lia,li,
     >   chii,minchii,gki,ni,etot,theta,ne2e,
     >   nqmf,psif,maxpsif,ef,lfa,lf,chif,minchif,
     >   gkf,nf,lg,const,uf,ui,nchf,nchi,nold,nznuc,npk,vmat)
      include 'par.f'
      common/meshrr/ meshr,rmesh(maxr,3)
      dimension chii(meshr,nqmi),minchii(nqmi),gki(kmax)
      dimension chif(meshr,nqmf),minchif(nqmf),gkf(kmax)
      dimension psii(maxr), psif(maxr), uf(maxr), ui(maxr), npk(nchan+1)
      real vmat(nf,ni+1), ovlpf(kmax), ovlpkf(kmax)
      common /pspace/ nabot,labot,natop,latop,ntype,ipar,nze,ninc,linc
      dimension nabot(0:lamax), natop(0:lamax)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),maxpsinb(nnmax,0:lnabmax)

      if (li.eq.lf.and.lia.eq.lfa.and.li.le.latop) then

         ovlp = 0.0
         do i = 1, min(maxpsif,maxpsii)
           ovlp = ovlp + psii(i) * psif(i) * rmesh(i,3)
         enddo
         tmp = - ovlp * const * etot * theta / 2.0
         do n = nabot(li), natop(li)
            ep = enpsinb(n,li)
            do kf = 1, nqmf
               ovlpf(kf) = 0.0
               do i = 1, maxpsinb(n,li)
                  ovlpf(kf) = ovlpf(kf) + psinb(i,n,li) * chif(i,kf)
               enddo
            enddo 
            do ki = 1, nqmi
               kii = npk(nchi) + ki - 1
               ovlpi = 0.0
               do i = 1, maxpsinb(n,li)
                  ovlpi = ovlpi + psinb(i,n,li) * chii(i,ki)
               enddo
               do 5 kf = 1, nqmf
                 if (ne2e.eq.0) then
                     kff = npk(nchf) + kf - 1
                     if (kff.lt.kii) go to 5
                     vmat(kii,kff+1) = vmat(kii,kff+1) +
     >                  tmp * ovlpf(kf) * ovlpi
                  else
                     kff = nchf
                  endif 
                  vmat(kff,kii) = vmat(kff,kii) +
     >               tmp * ovlpf(kf) * ovlpi
 5             continue 
            enddo
         enddo 
      endif
      if (li.eq.lfa.and.lf.eq.lia) then
         sign = (-1) ** (li + lia - lg)         
         do kf = 1, nqmf
            ekf = abs(gkf(kf)) * gkf(kf)
            ovlpf(kf) = 0.0
            ovlpkf(kf) = 0.0
            do i = minchif(kf), maxpsii
               ovlpf(kf) = ovlpf(kf) + psii(i) * chif(i,kf)
               ovlpkf(kf) = ovlpkf(kf) + psii(i) * chif(i,kf) *
     >            (uf(i) / 2.0 + 1.0 / rmesh(i,1))
            enddo
         enddo 
         do ki = 1, nqmi
            kii = npk(nchi) + ki - 1
            ovlpi = 0.0
            ovlpki = 0.0
            eki = abs(gki(ki)) * gki(ki)
            do i = minchii(ki), maxpsif
               ovlpi = ovlpi + psif(i) * chii(i,ki)
               ovlpki = ovlpki + psif(i) * chii(i,ki) *
     >            (ui(i) / 2.0 + 1.0 / rmesh(i,1)) 
            enddo
            do 10 kf = 1, nqmf
               kff = npk(nchf) + kf - 1
               if (kff.lt.kii) go to 10
               ekf = abs(gkf(kf)) * gkf(kf)

               eterm = 0.5 * (eki + ekf - etot * (1.0 - theta))
               
               tmp = eterm * ovlpi * ovlpf(kf) - ovlpki * ovlpf(kf) -
     >            ovlpkf(kf) * ovlpi
  
               if (nold.eq.0)
     >            tmp = (ei + ef - etot * (1.0 - theta)) / 2.0 *
     >            ovlpi * ovlpf(kf)    
               
C  subtract exchange from the direct part for triplet scattering
               vmat(kii,kff+1) = vmat(kii,kff+1) - tmp * const * sign
C  add exchange to the direct part for singlet scattering
               vmat(kff,kii) = vmat(kff,kii) + tmp * const * sign
 10         continue 
         enddo 
      endif 
      end

C  Define the direct matrix elements
      subroutine makev31d(nqmi,psii,maxpsii,lia,li,chii,minchii,
     >   gki,phasei,ni,nqmf,psif,maxpsif,lfa,lf,chif,minchif,gkf,
     >   phasef,nf,lg,rnorm,u,itail,nznuc,vdon,nchf,nchi,torf,
     >   npk,ne2e,vmat)
      include 'par.f'
      common/meshrr/ meshr,rmesh(maxr,3)
      dimension chii(meshr,nqmi),minchii(nqmi),gki(kmax),
     >   phasei(kmax)
      dimension chif(meshr,nqmf),minchif(nqmf),gkf(kmax),
     >   phasef(kmax),npk(nchan+1)
      complex phasei, phasef
      dimension fun(maxr), temp(maxr), psii(maxr), psif(maxr),
     >   temp3(maxr), u(maxr)
      real vmat(nf,ni+1),vdon(nchan,nchan,0:1)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
      common /pspace/ nabot,labot,natop,latop,ntype,ipar,nze,ninc,linc
      dimension nabot(0:lamax), natop(0:lamax)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),maxpsinb(nnmax,0:lnabmax)
      common/matchph/rphase(kmax,nchan)
      real chitemp(maxr,kmax)
      logical torf
      common/smallr/ formcut,regcut,expcut,fast
      logical fast
      
      minfun = 1
      maxfun = min(maxpsii,maxpsif)
      do i = minfun, maxfun
         fun(i) = psii(i) * psif(i) * rmesh(i,3)
      end do 

      do i = 1, maxr
         temp3(i) = 0.0
      enddo
      mini = maxr
      maxi = 1
      ctemp = 0.0

      do 10 ilt = -lia, lia, 2
         lt = lfa + ilt
         const = 1.0
         call form(fun,minfun,maxfun,rpow1(1,lt),rpow2(1,lt),
     >      minrp(lt),maxrp(lt),meshr,temp,i1,i2)

C  Subtract 1/r, but only for scattering calculations as opposed to e2e
         if (lt.eq.0) then
            call nuclear(fun,torf,minfun,maxfun,i2,u,nznuc,temp)
         endif 
         do i = i1, i2
            temp3(i) = const * temp(i) + temp3(i)
         enddo
         mini = min(mini,i1)
         maxi = max(maxi,i2)
 10   continue
         
C  As both CHII and CHIF contain the integration weights, we divide TEMP by   
C  them.
      do i = mini, maxi
         temp3(i) = - float(nze) * rnorm * temp3(i) / rmesh(i,3)
      end do
      mini1 = mini

      do ki = 1, nqmi
         minki = max(mini1, minchii(ki))
         do i = minki, maxi
            chitemp(i,ki) = temp3(i) * chii(i,ki)
         enddo
      enddo 

c$par doall      
      do ki = 1, nqmi
         kii = npk(nchi) + ki - 1
         minki = max(mini1, minchii(ki))
         do 20 kf = 1, nqmf
            kff = npk(nchf) + kf - 1
            if (kff.lt.kii) go to 20
            mini = max(minki, minchif(kf))
            n = maxi - mini + 1
            tmp = 0.0
            if (n.gt.0) tmp = sdot(n,chif(mini,kf),1,chitemp(mini,ki),1)
c$$$            do i = mini, maxi
c$$$               tmp = tmp + chif(i,kf) * chitemp(i,ki)
c$$$            end do 
C  Store the direct matrix elements in both the singlet and triplet triangles
            vmat(kff,kii) = vmat(kff,kii) + tmp
            vmat(kii,kff+1) = vmat(kii,kff+1) + tmp
 20      continue 
      end do 
      end

C  Define two-electron exchange matrix elements
      subroutine makev3e(psii,maxpsii,lia,nchi,psif,maxpsif,lfa,nchf,
     >   li,lf,chii,minchii,nqmi,chif,minchif,nqmf,lg,rnorm,npk,
     >   vmat,ichi,nchtop)
      include 'par.f'
      dimension chii(meshr,nqmi),minchii(nqmi),chif(meshr,nqmf),
     >   minchif(nqmf)
      dimension fun(maxr), temp(maxr), psii(maxr), psif(maxr),
     >   temp2(maxr), npk(nchtop+1)
      real vmat(ichi,ichi+1)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
      common/meshrr/ meshr,rmesh(maxr,3)
      common /pspace/ nabot,labot,natop,latop,ntype,ipar,nze,ninc,linc
      dimension nabot(0:lamax), natop(0:lamax)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),maxpsinb(nnmax,0:lnabmax)

      do ki = 1, nqmi
         kii = npk(nchi) + ki - 1
         minfun = minchii(ki)
         maxfun = maxpsif
         do i = minfun, maxfun
            fun(i) = chii(i,ki) * psif(i)
         end do
         do i = 1, maxr
            temp2(i) = 0.0
         enddo
         mini = maxr
         maxi = maxpsii
         do 15 ilt = -lia, lia, 2
            lt = lf + ilt
            call form(fun,minfun,maxfun,rpow1(1,lt),
     >         rpow2(1,lt),minrp(lt),maxrp(lt),maxi,temp,i1,i2)
            const = rnorm
            do i = i1, i2
               temp2(i) = const * temp(i) + temp2(i)
            enddo
            mini = min(mini,i1)
 15      continue
         do i = mini, maxi
            temp2(i) = temp2(i) * psii(i)
         enddo
         mini1 = mini
         do 30 kf = 1, nqmf
            kff = npk(nchf) + kf - 1
            if (kff.lt.kii) go to 30
            tmp = 0.0
            mini = max(mini1, minchif(kf))
            n = maxi - mini + 1
            if (n.gt.0) tmp = sdot(n,chif(mini,kf),1,temp2(mini),1)
c$$$            do i = mini, maxi
c$$$               tmp = tmp + chif(i,kf) * temp2(i)
c$$$            end do 
C  subtract exchange from the direct part for triplet scattering
            vmat(kii,kff+1) = vmat(kii,kff+1) - tmp
C  add exchange to the direct part for singlet scattering
            vmat(kff,kii) = vmat(kff,kii) + tmp
 30      continue 
      end do
      end
C  File: solvet.f
C  The following routine solves the equation
C  T_fi(X_f,X_i) = V_fi(X_f,X_i) * phi(X_f) * phi(X_i) 
C   + Sum_l{ Sum_n w_nl*V_fn(X_f,X_n)*phi(X_f)*conjg(phi(X_n))*T_ni(X_n,X_i)
C   - i * pi * X_l * V_fl(X_f,X_l) * phi(X_f)*conjg(phi(X_l)) * T_fl(X_f,X_l)}.
C  Here sum over the l for the first term is the sum over all channels,
C  whereas the second term is summed over only the open channels.
C  We do this using real arithmetic by letting
C  Sum_i T_fi(X,X_i) * conjg(phi(X_f)) * conjg(phi(X_i)) *
C                 (I(i,i') + i * pi X_i * K_ii'(X_i,X_i')) = K_fi'(X,X_i'),
C  and then solving K(f,i) = V(f,i) + w(n) * V(f,n) * K(n,i),
C  where the sum over n is implied. We want to do this in a
C  symmetric manner so as w(f) = w(i) = 0 we solve for K(n,i) with n <> f
C  by writing (I(f',n)/w(n) - V(f',n)) (w(n) * K(n,i)) = V(f',i),
C  where f' ranges over all n <> f. This is solved for B(n,i) = w(n) * K(n,i),
C  and the required K(f,i) is then formed by K(f,i) = V(f,i) + V(f,n) * B(n,i).
C  This works fine, but gives enormous determinants, so instead we
C  solve (I(f',n) * sig(w(n)) - V'(f',n)) (sig(w(n)) * K'(n,i)) = V'(f',i),
C  where sig(w(n)) = w(n)/|w(n)|, V'(f',n) = V(f',n) * sqrt|w(f')*w(n)|,
C  K'(n,i) = K(n,i) * sqrt|w(n)| and V'(f',i) = V(f',i) * sqrt|w(f')|.
C  The required K(f,i) is recovered by the same equation as before.
      subroutine gettmat(lg,ns,nchtop,iex,npk,vmat,ichi,wk,gk,phasel,
     >   sigma,tdist,von,etot,kon,ton,err,vmatop,nchopt,nsmax,
     >   nds,v)
      include 'par.f'
      parameter (nmax=kmax*nchan, lwork = 70 * nmax)
      integer npk(nchtop+1)
      real vmat(ichi,ichi+1)
      complex vmatop(kmax,kmax,0:nchanop,nchanop),wk(kmax*nchan)
      real v(ichi,nchtop,2),work(lwork)
      complex ton(nchan,nchan),von(nchan,nchan), c,
     >   kon(nchan,nchan), ckern(nchan,nchan), cwork(nchan)
      complex phasel(kmax,nchan), tdist(nchan), sigma(nchan)
      dimension det(2), gk(kmax,nchan), err(nchan,nchan)
      character ud(0:1), ch
      common /pspace/ nabot,labot,natop,latop,ntype,ipar,nze,ninc,linc
      dimension nabot(0:lamax), natop(0:lamax)
      logical sprint, lprint
      data pi/3.1415927/
      data ud/'L','U'/
      data sprint,lprint/.false.,.false./
      ch(n) = char(n + ichar('0'))
      
      r = float((-1)**ns * min(1,iex))
      nd = npk(nchtop+1) - 1         
     
C  Form the driving vector and the first Born

      do nchi = 1, nchtop
         kf = 1
         call makev(vmat,kf,nchi,nchtop,npk(nchtop+1),npk,wk,ns,
     >      v(1,nchi,1))
         do n = 1, nd
            v(n,nchi,2) = v(n,nchi,1)
         enddo
         do nchf = 1, nchtop
            ton(nchf,nchi) = onshellk(vmat,nchtop,nchf,nchi,ns,npk,
     >         npk(nchtop+1))
         enddo
      enddo

      if (ns.eq.0) then
         do ki = 1, nd
            do kf = ki, nd
               vmat(kf,ki) = - vmat(kf,ki) *
     >            sqrt(abs(real(wk(kf))*real(wk(ki))))
C  Can comment out the above line and make a few more changes below
C  where sqrt(wk) occurs to get the simple documented way of solving 
C  the equations, but this gives more ill-conditioned matrices
            enddo
         enddo
      else
         do ki = 1, nd
            do kf = ki, nd         
               vmat(ki,kf+1) = - vmat(ki,kf+1)
     >            * sqrt(abs(real(wk(kf))*real(wk(ki))))
            enddo
         enddo 
      endif 
      
      epsil = 1e-4
C  Define the on-shell V, used for printing out only
      do nchi = 1, nchtop
         do nchf = 1, nchtop
            von(nchf,nchi) = ton(nchf,nchi) * phasel(1,nchf) *
     >         phasel(1,nchi) / gk(1,nchi) / gk(1,nchf)
         end do 
      end do 


C  Add the I matrix to -K
      do kf = 1, nd
         s = (real(wk(kf))+1e-30) / abs(real(wk(kf))+1e-30)
         vmat(kf,kf+ns) = vmat(kf,kf+ns) + s
      end do 

      call clock(s1)

      nds = npk(nchtop+1) - 1 - nchtop
      do nchi = 1, nchtop
         do nchf = 1, nchtop
            do kf = npk(nchf) + 1, npk(nchf+1) - 1
               v(kf-nchf,nchi,1) = v(kf,nchi,1)
               v(kf-nchf,nchi,2) = v(kf,nchi,1)
            enddo 
         enddo
      enddo      
      if (ns.eq.0) then
         do nchi = 1, nchtop
            do ki = npk(nchi) + 1, npk(nchi+1) - 1
               do kf = ki, nd
                  vmat(kf,ki-nchi) = vmat(kf,ki)
               enddo
            enddo
         enddo
         do nchf = 1, nchtop
            do kf = npk(nchf) + 1, npk(nchf+1) - 1
               do ki = 1, kf - nchf
                  vmat(kf-nchf,ki) = vmat(kf,ki)
               enddo
            enddo
         enddo
      else
         do nchf = 1, nchtop
            do kf = npk(nchf) + 1, npk(nchf+1) - 1
               do ki = kf - nchf + 1, nd + 1
                  vmat(kf-nchf,ki) = vmat(kf,ki)
               enddo
            enddo
         enddo
         do nchi = 1, nchtop
            do ki = npk(nchi) + 1, npk(nchi+1) - 1
               do kf = 1, ki - 1
                  vmat(kf,ki-nchi+1) = vmat(kf,ki+1)
               enddo
            enddo
         enddo
      endif
      
C  Solve the linear equations
      call matinv(vmat(1,1+ns),ud(ns),nd,nds,v,nchtop,work,lwork,
     >   det,rc,0)

      call clock(s2)
      do nchf = 1, nchtop
         do nchi = 1, nchtop
            if (nds.eq.0) then
               sum = 0.0
            else
               sum = sdot(nds,v(1,nchf,2),1,v(1,nchi,1),1)
            endif 
            ton(nchf,nchi) = ton(nchf,nchi) + sum
         enddo
      enddo 
C  Write out the half-off-shell K matrix
      do nchi = 1, 1  ! nchtop
         do nchf = 1, nchtop
         if (gk(1,nchi).gt.0.0.and.gk(1,nchf).gt.0.0) then
            if (ns.eq.0) then
               open(42,file='singlet.'//ch(nchf)//ch(nchi),
     >            status='unknown')
            else
               open(42,file='triplet.'//ch(nchf)//ch(nchi),
     >           status='unknown' )
            endif 
            write(42,'(''#       k'',8x,''K(k)'',8x,''V(k)     '',
     >         '' for transition:'', i2, '' <-'',i2)') nchf,nchi
            divk = gk(1,nchi) * gk(1,nchf)
               write(42,'(f12.6,1p,2e12.3,''   On-shell point'')') 
     >         gk(1,nchf),real(ton(nchf,nchi))/divk,real(von(nchf,nchi))
            do n = npk(nchf)+1, npk(nchf+1) - 1
               kf = n - npk(nchf) + 1
               divk = gk(1,nchi) * gk(kf,nchf)
               write(42,'(f12.6,1p,2e12.3)') gk(kf,nchf), 
     >            v(n-nchf,nchi,1)/divk/real(wk(n))*sqrt(abs(wk(n))),
     >            v(n-nchf,nchi,2)/sqrt(abs(wk(n)))/divk
            enddo
            close(42)
         endif 
         enddo
      enddo 
C  Get the T-matrix
      do nchi = 1, nchtop
C  The following lines are redundant
         do nchf = 1, nchtop
            if (gk(1,nchf).le.0.0.or.gk(1,nchi).le.0.0) then
               c = 0.0
            else
               c = 1.0
            endif 
            kon(nchf,nchi) = ton(nchf,nchi)/gk(1,nchi)/gk(1,nchf)
C  The variable KON is the K matrix K(f,i). To get the T matrix we use
C  Sum_i T_fi(X_f,X_i) * conjg(phi(X_f)) * conjg(phi(X_i)) *
C                 (I(i,i') + i * pi X_i * K_ii'(X_i,X_i')) = K_fi'(X_f,X_i').
C  As I, K and T are all symmetric we write this as (f <-> i')
C  Sum_i (I(f,i) + i * pi X_i * K_fi(X_f,X_i)) *
C   * T_ii'(X_i,X_i') * conjg(phi(X_i')) * conjg(phi(X_i)) = K_fi'(X_f,X_i'),
            ckern(nchf,nchi) = c * kon(nchf,nchi) *
     >         cmplx(0.0, pi*gk(1,nchi))
C  Maybe multiply sum by c in the next line also to get 0 for closed
C  channels when printing out NCHTOP T matrix elements
            ton(nchf,nchi) = c * kon(nchf,nchi)
         end do 
      end do
      do nch = 1, nchtop
         ckern(nch,nch) = ckern(nch,nch) + (1.0,0.0)
      end do 
C  Solve the linear equations
      nv = nchtop
      nd = nchtop
      call matinv2(ckern,nchan,nd,ton,nv,cwork,erfp,epsil)
      do nchi = 1, 1
      enddo 
      do nchi = 1, nchtop
         do nchf = 1, nchtop
            ton(nchf,nchi) = ton(nchf,nchi) * phasel(1,nchf) *
     >         phasel(1,nchi)
         enddo
C  In the case of Spin = 1.5 we musn't add TDIST to those TMATs which are 0
         if (abs(ton(nchi,nchi)).gt.0.0) then
            ton(nchi,nchi) = ton(nchi,nchi) + tdist(nchi)
            von(nchi,nchi) = von(nchi,nchi) + tdist(nchi)
         endif 
      enddo
C  The following inverts the T matrix to get the on-shell K matrix
C  It should always come out to be essentially real
      do nchi = 1, nchtop
         do nchf = 1, nchtop
            if (gk(1,nchf).le.0.0.or.gk(1,nchi).le.0.0) then
               c = 0.0
            else
               c = 1.0/sigma(nchf)/sigma(nchi)
            endif
            kon(nchf,nchi) = c * ton(nchf,nchi)
            ckern(nchf,nchi) = - c * ton(nchf,nchi) *
     >         cmplx(0.0, pi*gk(1,nchi))
         end do
         ckern(nchi,nchi) = ckern(nchi,nchi) + (1.0,0.0)
      enddo
      call matinv2(ckern,nchan,nd,kon,nv,cwork,erfp,epsil)
      return
      end
      
      subroutine makev(vmat,kf,nchf,nchtop,npktop,npk,wk,ns,vec)
      include 'par.f'
      integer npk(nchtop+1)
      real vmat(npktop-1,npktop),vec(*)
      complex wk(*)

      kff = npk(nchf) + kf - 1
      if (ns.eq.0) then
         do nchn = nchf, nchtop
            do knn = npk(nchn) + 1, npk(nchn+1) - 1
               vec(knn) = vmat(knn,kff) * sqrt(abs(wk(knn))) 
            end do
         end do
         do nchn = 1, nchf - 1
            do knn = npk(nchn) + 1, npk(nchn+1) - 1
               vec(knn) = vmat(kff,knn) * sqrt(abs(wk(knn))) 
            end do
         end do
      else
         do nchn = nchf, nchtop 
            do knn = npk(nchn) + 1, npk(nchn+1) - 1
               vec(knn) = vmat(kff,knn+1) * sqrt(abs(wk(knn))) 
            end do
         end do
         do nchn = 1, nchf - 1
            do knn = npk(nchn) + 1, npk(nchn+1) - 1
               vec(knn) = vmat(knn,kff+1) * sqrt(abs(wk(knn))) 
            end do
         end do
      endif 
      end
      
      function onshellk(vmat,nchtop,nchf,nchi,ns,npk,npktop)
      include 'par.f'
      integer npk(nchtop+1)
      real vmat(npktop-1,npktop)
      
      if (ns.eq.0) then
         if (nchf.ge.nchi) then
            onshellk = vmat(npk(nchf),npk(nchi))
         else 
            onshellk = vmat(npk(nchi),npk(nchf))
         endif 
      else
         if (nchf.ge.nchi) then
            onshellk = vmat(npk(nchi),npk(nchf)+1)
         else 
            onshellk = vmat(npk(nchf),npk(nchi)+1)
         endif 
      end if
      return
      end
      
C  The following routine solves the set of linear equations 
      subroutine matinv(kernel,ud,lda,n,v,m,work,lwork,determ,rcond,i)
      include 'par.f'
      parameter (nmax = kmax * nchan)
      real kernel(lda,lda)
      real v,work,determ(2)
      dimension v(lda,m,2), work(lwork), ipvt(nmax)
      character ud
      
      if (n.le.0) return
      ia=lda
      ib=lda
      ijob=i
      call ssysv(ud,n,m,kernel,lda,ipvt,v,lda,work,lwork,info)
      if (lwork.lt.work(1).or.nint(work(1)/n).ne.64)
     >   print*,'previous and optimal LWORK, NB:',lwork,work(1),
     >   work(1)/N
      if (info.ne.0) print*,'Warning: INFO, N:', info, N
      return
      end
C  File rest.f, where the linear equations are formed and solved
      subroutine solvet(iex,vmat,gk,wk,weightk,nchtop,nqm,noprint,
     >   nopen,etot,lg,vdon,phasel,nent,isecond,nunit,sigma,ovlp,
     >   lfast,lslow,slowery,nznuc,zasym,npk,
     >   projectile,target,uba,nsmax)
      include 'par.f'
      integer npk(*)
      real vmat(*),
     >   vdon(nchan,nchan,0:1),v(2*kmax*nchan*nchan)
      complex vmatop(kmax,kmax,0:nchanop,nchanop),wk(kmax*nchan),
     >   cwork(66 * nchan),eigv(nchan)
      complex sigma(nchan),
     >   ctmp,smat(nchan,nchan),wsmat(nchan,nchan)
      complex ton(nchan,nchan,0:1),von(nchan,nchan),t2nd(nchan,nchan)
      complex wton(nchan,nchan), wvon(nchan,nchan), wt2nd(nchan,nchan),
     >   wvdon(nchan,nchan), phasel(kmax,nchan), tdist(nchan)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),maxpsinb(nnmax,0:lnabmax)
      common /pspace/ nabot,labot,natop,latop,ntype,ipar,nze,ninc,linc
      dimension nabot(0:lamax), natop(0:lamax)
      common/meshrr/ meshr,rmesh(maxr,3)
      dimension jlm(nchan),non(nchan,2*lamax+1),partcs(nchan,nchan,0:1),
     >   weightk(kmax),gk(kmax,nchan), nychan(nchan), rwork(3*nchan),
     >   temp(maxr), err(nchan,nchan), slowery(ncmax)
      common/smallr/ formcut,regcut,expcut,fast
      common /double/id,jdouble(20)
      common /radpot/ ucentr(maxr)
      logical fast
      real sigtop(nchan,0:1)
      real ovlp(ncmax,0:lamax),esum(0:1)
      logical same, sames, noprint(nchan), uba(nchan)
      character spin(0:1), projectile*8,
     >   target*6, chan(100)*3
      common /charchan/ chan
      common /chanens/ enchan(100)
      save sames
      data pi,method,sames/3.1415927,1,.false./
      data spin/'+','-'/

      if (2*kmax*nchan*nchan.lt.(npk(nchtop+1)-1)*nchtop*2)
     >   stop 'Increase either KMAX or NCHAN in par.f'
      nch = 1
      n = 0
      print'('' JS f i  real(T)   imag(T)    real(S)   imag(S) '
     >   //'     V         K     eigenphase'')'
      
      same = .true.
      if (iex.eq.0) sames = .true.
      do ns = 0, nsmax
         if (sames.and.ns.eq.1) then
            do nchi = 1, nchtop
               do nchf = 1, nchtop
                  ton(nchf,nchi,1) = ton(nchf,nchi,0)
               enddo
            enddo
         else
C  The following routine forms and solves the linear equations
            call gettmat(lg,ns,nchtop,iex,npk,vmat,npk(nchtop+1)-1,
     >      wk,gk,phasel,sigma,tdist,von,etot,t2nd,ton(1,1,ns),err,
     >      vmatop,nchopt,nsmax,nds,v)
         endif 
         r = float((-1)**ns * min(1,iex))

C  The following initialization uses NCHAN instead of NCHTOP because
C  The partial cross sections are labeled by the state number which
C  may be greater than than NCHTOP when IPAR = 1.
         do nchi = 1, nchan
            sigtop(nchi,ns) = 0.0
            do nchf = 1, nchan
               partcs(nchf,nchi,ns) = 0.0
               smat(nchf,nchi) = (0.0,0.0)
            enddo
         enddo 
C  Define the S matrix from the on-shell T matrix
         do nchi = 1, nchtop
            smat(nchi,nchi) = (1.0,0.0)
            do nchf = 1, nchtop
               ctmp = sigma(nchf) * sigma(nchi) + (1e-10,0.0)
               gkfac = abs(gk(1,nchi) * gk(1,nchf))
               smat(nchf,nchi) = smat(nchf,nchi) - 2.0 * pi *
     >            (0.0,1.0) * sqrt(gkfac) * ton(nchf,nchi,ns) / ctmp
               wsmat(nchf,nchi) = smat(nchf,nchi)
            enddo
         enddo
         do nchi = 1, nchtop
            do nchf = 1, nchtop
               if (nchi.eq.nchf) then
                  ctmp = (-1.0,0.0)
               else
                  ctmp = (0.0,0.0)
               endif
               do nchn = 1, nchtop
                  ctmp = ctmp + conjg(smat(nchf,nchn)) * smat(nchn,nchi)
               enddo
               err(nchf,nchi) = abs(ctmp)
               if (abs(ctmp).gt.1e-4.and.nchopt.eq.0)
     >            print*,'S matrix not unitary', nchf,nchi,ctmp
            enddo
         enddo
         lwork = 66 * nchan
         ldvs = 1
 
C  Diagonalize the S matrix to obtain eigenphases
         call cgees('n','n',select,nchtop,smat,nchan,ndim,eigv,vs,ldvs,
     >      cwork,lwork,rwork,bwork,info)
         if (info.ne.0) print*,'INFO from CGEES:',info
         esum1 = 0.0
         esum2 = 0.0
         do nch = 1, nchtop
            esum1 = esum1 + real(log(eigv(nch))/2.0/(0.0,1.0))
            esum2 = esum2 +
     >         atan2(aimag(eigv(nch)),real(eigv(nch))) / 2.0
         enddo
         esum(ns) = esum1
         nchi = 0
         nchip = 0
         do while (nchip.le.nent.and.nchi.lt.nchtop)
            nchi = nchi + 1
c$$$         do nchi = 1, lent
            call getchinfo (nchi,nchip,lg,temp,maxpsi,ei,lia,nia,li)
            if (nchip.le.nent) then
               if (chan(nchip)(1:1).eq.'t') then
C  Incident on a triplet state of Helium-like target. Total spin S = ns + 0.5.
                  si = 1.0
                  Stot = float(ns) + 0.5
               elseif (chan(nchip)(1:1).eq.'s') then
C  Incident on a singlet state of Helium-like target. Total spin S = ns + 0.5.
                  si = 0.0
                  Stot = float(ns) + 0.5
               else
C  Incident on a hydrogenic target. Total spin S = ns.
                  si = 0.5
                  Stot = float(ns)
               endif

               if (nze.eq.1) then
                  spinw = 1.0
               else
                  spinw = (2.0 * Stot + 1.0) / 2.0 / (2.0 * si + 1.0)
               endif
               
               lent = nchi
               nymax = 0
               nfp = 0
               lfp = 0
               no = 0
               do nchf = 1, nchtop
                  call getchinfo (nchf,nchp,lg,temp,maxpsi,ef,lfa,nfa,
     >               lf)
                  const = spinw * (2.0 * pi) ** 4 / (4.0 * pi) *
     >               (2.0 * lg + 1.0) * gk(1,nchf) / gk(1,nchi) /
     >               (2.0 * lia + 1.0)
                  lt = lg + lfa
                  lb = abs(lg - lfa)
            
                  if (noprint(nchp)) then
                     nymax = nymax + 1
                     nychan(nymax) = nchf
                     if (nfa.ne.nfp.or.lfa.ne.lfp) no = no + 1
                     nfp = nfa
                     lfp = lfa
                     jlm(no) = lt - lb + 1
                     non(no,lf-lb+1) = nymax
                     gkfac = gk(1,nchi) * gk(1,nchf)
                     wvdon(nymax,nchi)= cmplx(vdon(nchf,nchi,ns)/gkfac)
                     wton(nymax,nchi) = ton(nchf,nchi,ns)
                     wt2nd(nymax,nchi) = t2nd(nchf,nchi)
                     wvon(nymax,nchi) = von(nchf,nchi)
                  end if
C  Define total cross sections for channel pairs
                  if (gk(1,nchf).gt.0.0) then
                     at2 = abs(ton(nchf,nchi,ns))**2
                     partcs(nchp,nchip,ns) = partcs(nchp,nchip,ns) +
     >                  at2 * const
C  On occasion there will be a closed channel for nchp < nchpmax
C  This is no big deal, as its partcs = 0.0
                     nchpmax = max(nchp,nent)
                  endif
C  Define the TOTAL cross section (TCS) from the imaginary part of the elastic
C  T matrix. 
                  if (nchp.eq.nchip.and.nchi.eq.nchf) then
                     sigtop(nchip,ns) = sigtop(nchip,ns) - const *
     >                  aimag(ton(nchi,nchi,ns)) / (pi * gk(1,nchf))
                  endif 
C  End of the "do nchf = 1, nchtop" loop
               end do
C  End of the "if (nchip.le.nent) then" statement
            endif 
            if (no.ne.nopen) then
               print*,'NO and NOPEN are not equal',no,nopen
               stop 'NO and NOPEN are not equal'
            endif
         end do
         sigt = 0.0
         sum = 0.0
         do nch = 1, nchpmax
            sigt = sigt + partcs(nch,1,ns)
         enddo
C  Print out the results
         j = mod(lg,100)
         do nchf = 1, nymax
            do nchi = 1, lent
               ctmp = log(eigv(nchf))/2.0/(0.0,1.0)
               rdel = real(ctmp)
               print '(i2,a,2i2,1p,3(2(e10.3),1x),0p,1f7.3)',
     >            j,spin(ns),nychan(nchf),nchi,
     >         wton(nchf,nchi), wsmat(nchf,nchi), real(wvon(nchf,nchi)),
     >         real(wt2nd(nchf,nchi)), rdel
               same = same.and.abs((wvon(nchf,nchi)-wton(nchf,nchi))/
     >            (wton(nchf,nchi)+1e-30)).lt.1e-3
            end do 
         end do
         if (abs((sigt-sigtop(1,ns))/(sigtop(1,ns)+1e-30)).gt.1e-4) 
     >      print*,'Warning: optical theorem is only satisfied to',
     >      abs((sigt-sigtop(1,ns))/(sigtop(1,ns)+1e-30))
         print '('' N='',i5,'' eigenphase sum:'',f7.3)', nds, esum(ns)
         call update(6)
C  End the spin loop
      end do 

C  Write out total cross sections
      if (nunit.gt.0) 
     >   call wrtcs(partcs,sigtop,nchpmax,nent,lg,etot,nsmax,
     >   ovlp,nunit,nznuc,zasym,projectile,target)
      end

      function nmpow(err)
      nmpow = 0
      do while (err * 10 ** nmpow .lt. 1.0.and.nmpow.lt.9)
         nmpow = nmpow + 1
      enddo
      return
      end

C  The following routine sets the K-grids      
      subroutine kgrid(nkor,skor,etot,gk,wk,weightk,nbnd,nqm,lg,nchtop,
     >   nchopt,npk,nold,nptop,lptop,ifirst,ntstop,nnbtop,hlike,uba)
      include 'par.f'
      complex wk(kmax*nchan)
      dimension nkor(5,0:lmax),skor(5,0:lmax),psin(maxr)
      dimension nk(10),sk(0:10),gridk(kmax),weightk(kmax),psi(maxr)
      double precision xx(kmax), ww(kmax),wf(2*kmax),dstart,dstop
      integer iwf(2*kmax),nbnd(0:lmax)
      dimension gk(kmax,nchan), gfixed(kmax), wfixed(kmax),npk(nchan+1)
      common /pspace/ nabot,labot,natop,latop,ntype,ipar,nze,ninc,linc
      common /worksp/
     >   ps2(maxr,ncmax),psen2(ncmax),minps2(ncmax),maxps2(ncmax)
      dimension nabot(0:lamax), natop(0:lamax), reg(maxr),ucentr(maxr)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      common/meshrr/ meshr,rmesh(maxr,3)
      common /double/njdouble,jdouble(20)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
      common/smallr/ formcut,regcut,expcut,fast
      logical fast,hlike,usetrapz,uba(nchan)
      common/cont/ ce(ncmax),cint(ncmax),ncstates,energy
      character chan(100)*3
      common /charchan/ chan
      common /chanens/ enchan(100)
      data pi,ucentr/3.1415927,maxr*0.0/

C  NBAD stores the number of channels |phi(n)> such that
C  Int dk <phi(n)|k><k|phi(n)> .ne. 1
      nbad = 0
      lprint = 0
      if (lprint.eq.0) lprint = max(latop,lg)
      
      if (lg.le.lprint) print '('' Kgrid quadrature set'')'
      nch = 1
      npk(nch) = 1
      call getchinfo (nch, ntmp, lg, psi, maxpsi, ea, la, na, li)
      if (nch.eq.0) stop 'NCH = 0 in KGRID'
C  The following is a loop of the form REPEAT ... UNTIL(condition)
 10   e = etot - ea
      if (e.ge.0.0) then
         rk = sqrt(e)
      else
         rk = - sqrt(-e)
      end if 
      nqk = 0
      lset = li
      nquad = 0
      do n = 1, 4
         nquad = nquad + nkor(n,lset)
      enddo 
      if (nqm.eq.1.or.uba(nch).or.nquad.eq.0) go to 20

            
      npoints = nkor(1,lset)
      endk = skor(1,lset)
      npoints2 = nkor(2,lset)
      endk2 = skor(2,lset)
      nendk = nkor(3,lset)
      endp = skor(3,lset)
      width = skor(4,lset)
      midnp = abs(nkor(4,lset))
      usetrapz = nkor(4,lset).lt.0
      if (width.lt.0.0) then
         dstart = 0d0
         nt = midnp
         nwf = 2 * nt
         niwf = 2 * nt
         dstop  = dble(-2.0 * width)
         call cgqf(nt,xx,ww,1,0d0,0d0,dstart,dstop,
     >      0,nwf,wf,niwf,iwf,ier)
         do i = 1, nt
            gfixed(i) = real(xx(i))
            wfixed(i) = real(ww(i))
         enddo
      end if
      
      if (lg.le.0)
     >   print '('' interval  points    start       stop          '//
     >   'test'')'
      enk = endk
 30   sumi=0.0
      sumip = 0.0
      nqk=0
C  Define the intervals and the number of points in each interval
      call makeints(mint,sk,nk,rk,npoints,width,midnp,enk,nendk,endp,
     >   npoints2,endk2)
C  Obtain Gaussian knots and weights within each interval
      dstop = 0d0
      do j=1,mint-1
         nqk=nqk+nk(j)
         nt=nk(j)
         nwf=2*nt
         niwf=2*nt
         if (nt.gt.0) then
C  Can use the following two lines here and in the i loop to make kgrid
C  symmetric about E rather than sqrt(E). A suitable change must be made
C  in the MAKEINTS routine.
c$$$            dstart = dstop**2
c$$$            dstop  = dble(sk(j))**2
            dstart = dstop
            dstop  = dble(sk(j))
            call cgqf(nt,xx,ww,1,0d0,0d0,dstart,dstop,
     >         0,nwf,wf,niwf,iwf,ier)
c$$$            dstart = sqrt(dstart)
c$$$            dstop  = sqrt(dstop)

C  Can try using the trapezoidal rule, but it is not as good.
            if (j.ge.2.and.rk.gt.dstart.and.rk.lt.dstop.and.usetrapz)
     >         then
               call trapez(dstart,dstop,nt,xx,ww)
            print*,'Will use trapezoidal rule around the on-shell point'
            endif 
            do i=nqk-nk(j)+1,nqk
               gridk(i)=xx(i-nqk+nk(j))
               weightk(i)=ww(i-nqk+nk(j))
               ecmn = gridk(i) ** 2
               eta = 0.0
               if (lg.eq.0) then
                  do k = 1, maxpsi
                     reg(k) = sin(gridk(i) * rmesh(k,1))
                  end do 
                  call getprod(psi,maxpsi,1,reg,la,rmesh,
     >               meshr,ntmp,tmp)
                  sumi = sumi + tmp * weightk(i) * 2.0 / pi
               endif 
            end do
            tmpold = tmp
            if (lg.eq.0) then
               if (dstart.lt.rk.and.rk.lt.dstop) then
                  print '(i5,i10,f12.6,'' *'',f10.6,f13.5)', j,nk(j),
     >               dstart, dstop, sumi - sumip
               else 
                  print '(i5,i10,2f12.6,f13.5)', j,nk(j),dstart,dstop,
     >               sumi - sumip
               end if
            endif 
            sumip = sumi
         end if 
      end do 
C  Here we have the last interval
      nt=nk(mint)
      if (nt.gt.0) then
         if (dstop.lt.rk) then
            print*,'The last interval must start > than RK'
            stop 'The last interval must start > than RK'
         end if  
         nqk=nqk+nt
         p=sk(mint)
         dstart=0d0
         dstopp = dstop
         dstop=dstop**(1.0-p)
         nwf=2*nt
         niwf=2*nt
         call cgqf(nt,xx,ww,1,0d0,0d0,dstart,dstop,0,nwf,wf,niwf,iwf,
     >      ier)
         if (ier.ne.0) print*,'KGRID IER:',ier
         do j=nt,1,-1
            jj=nqk-j+1
            gridk(jj)=xx(j)**(1.0/(1.0-p))
            weightk(jj)=ww(j)/(p-1d0)*gridk(jj)**p
            ecmn = gridk(jj) ** 2
            if (lg.eq.0) then
               eta = 0.0
               do i = 1, maxpsi
                  reg(i) = sin(gridk(jj) * rmesh(i,1))
               end do 
               call getprod(psi,maxpsi,1,reg,la,rmesh,
     >            meshr,ntmp,tmp)
               if (abs(tmp).gt.abs(tmpold)*1.5) then
                  endkt = gridk(jj)
                  tmpold = tmp
               endif 
               sumi = sumi + tmp * weightk(jj) * 2.0 / pi
            endif 
         end do 
         if (lg.eq.0) then
            print '(i5,i10,f12.6,''       oo'',f16.5)', mint,nk(mint),
     >         dstopp, sumi - sumip
            print'('' fall off power:'',f5.1,22x,''= '',f7.5)',
     >         sk(mint), sumi
         endif 
         if (abs(sumi - 1.0).gt.1e-2.and.lg.eq.0) nbad = nbad + 1
      end if
c$$$      close(60+nch)
      if (nqk+1.gt.kmax) then
         print*,'NQK + 1 > KMAX:',nqk+1,kmax
         stop 'Increase KMAX'
      end if 

 20   continue
C  Check that the integration rule will handle the principle value singularity
      j=1
      if (e.gt.0.0) then
         if (nqk.gt.0) then
            sum=0.0
            nt=0
            do while (sk(j).lt.rk-0.01.and.j.lt.mint)
               sum = 2.0 * atanh(sk(j)/rk)/rk
               nt=nt+nk(j)
               j=j+1
            end do
            j=mint-1
            tsum=0.0
            do while (sk(j).gt.rk+0.01.and.j.ge.1)
               tsum = - 2.0 * acoth(sk(j)/rk) / rk
               j=j-1
            end do
            sum=sum+tsum
            im=nt+1
            do while (gridk(im).le.(sk(j+1)).and.im.lt.kmax)
               sum=sum + 2.0*weightk(im)/(e-gridk(im)*gridk(im))
               im=im+1
            end do
         else
            sum = 0.0
         endif 
         if (lg.le.lprint)
     >      print '('' State, NCH, NA, LA, L, K, EA:'',a4,4i4,1p,
     >      2e15.5)',chan(ntmp),nch,na,la,li,rk,ea
      else
         if (nqm.gt.1) then
            sum = 0.0
            do i = 1, nk(1) + nk(2)
               sum = sum + weightk(i) * 2.0 / (e - gridk(i) ** 2)
            enddo
            print*,'Test of closed channel integral:',- 2.0 / sqrt(-e)*
     >         atan(sk(2)/sqrt(-e)) / sum, sum
         endif 
         if (lg.le.lprint)
     >      print'('' State, NCH, NA, LA, L, E:    '',a4,4i4,1p,
     >      e15.5,''    closed'')',chan(ntmp),nch,na,la,li,e
      end if 
      if (nqm.le.0.and.lg.eq.0.and.ifirst.ne.0.and.hlike) then
         nchp = 1
         call getchinfo (nchp, ntm, lg, psin, mpsn, eat, lat,nt,lit)
         do while (nchp.ne.0)
            if (lat.eq.la) then
               call getchnl(chan(nchp),n,l)
               print'('' '',a3,'' '',$)', chan(nchp)
            endif 
            nchp = nchp + 1
            call getchinfo (nchp, ntm, lg, psin, mpsn, eat, lat,nt,lit)
         enddo
         print*
c$$$         if (la.eq.0)print'(30(i3,''s'',1x))',(n,n=nabot(la),natop(la))
c$$$         if (la.eq.1)print'(30(i3,''p'',1x))',(n,n=nabot(la),natop(la))
c$$$         if (la.eq.2)print'(30(i3,''d'',1x))',(n,n=nabot(la),natop(la))
c$$$         if (la.eq.3)print'(30(i3,''f'',1x))',(n,n=nabot(la),natop(la))
c$$$         if (la.eq.4)print'(30(i3,''g'',1x))',(n,n=nabot(la),natop(la))
c$$$         if (la.eq.5)print'(30(i3,''h'',1x))',(n,n=nabot(la),natop(la))
c$$$         if (la.eq.6)print'(30(i3,''i'',1x))',(n,n=nabot(la),natop(la))
         nchp = 1
         call getchinfo (nchp, ntm, lg, psin, mpsn, eat, lat,nt,lit)
         do while (nchp.ne.0)
            if (lat.eq.la) then
               sumi = 0.0
               do i = 1, nqk
                  ecmn = gridk(i) ** 2
                  eta = 0.0
                  do j = 1, mpsn
                     reg(j) = sin(gridk(i) * rmesh(j,1))
                  end do 
                  tmp = 0.0
                  call getprod(psin,mpsn,
     >               1,reg,lat,rmesh,meshr,ntm,tmp)
                  sumi = sumi + tmp * weightk(i) * 2.0 / pi
               enddo
               print '('' '',f4.2,$)',sumi
               if (abs(sumi - 1.0).gt.1e-2) nbad = nbad + 1
            endif
            nchp = nchp + 1
            call getchinfo (nchp, ntm, lg, psin, mpsn, eat, lat,nt,lit)
         enddo 
         print*
      endif 

      do i = 1, nqk + 1
         kp = npk(nch) + i - 1
         if (i.eq.1) then
            gk(i,nch) = rk
            if (e.ge.0.0) then
C  The T(kf,ki) matrix has been divided by KF and KI
               wk(kp) = - sum - cmplx(0.0, pi / rk)
c$$$                        wk(kp) = - e * sum - cmplx(0.0, pi * rk)
            else
               wk(kp) = (0.0,0.0)
            end if 
            if (abs(real(wk(kp))).gt.1e-2) then
               print*,'Check K-grid. On-shell weight:',wk(kp)
               stop 'Real part of on shell weight must be < 1e-2'
            endif 
         else
            gk(i,nch) = gridk(i-1)
            ek = gridk(i-1) * gridk(i-1)
C  The T(kf,ki) matrix has been divided by KF and KI
            wk(kp) = 2.0 * weightk(i-1)/(e - ek)
c$$$                     wk(kp) = ek * 2.0 * weightk(i-1)/(e - ek)
         end if
      end do
      if (lg.le.lprint) print*
      nchtop = nch
      if (nqm.gt.0) nqk = min(nqk,nqm-1)
      nqk = nqk + nbnd(lset)
      npk(nchtop+1) = npk(nchtop) + nqk + 1 
      if (la.le.lptop.and.na.le.nptop) nchopt = nchtop
      nch = nch + 1
      call getchinfo (nch, ntmp, lg, psi, maxpsi, ea, la, na, li)
      if (nch.gt.nchan) then
         print*,'Recompile with NCHAN >=',nch
         stop 'Recompile with bigger NCHAN'
      end if 
      if (nch.ne.0) go to 10
      
      if (nqm.gt.kmax) stop 'Increase KMAX'
      if (nbad.gt.0.and.lg.le.lprint) print*,'Warning, NBAD:', nbad
      return
      end
      
      subroutine getprod(psi,maxpsi,minchi,chi,l,rmesh,meshr,
     >   nch,ovlp)
      include 'par.f'
      dimension psi(maxr), chi(maxr), rmesh(maxr,3)
      
      ovlp = 0.0
      ovlp1 = 0.0
      do i = minchi, maxpsi
         ovlp = ovlp + chi(i) * psi(i) * rmesh(i,3)
      enddo
      ovlp = ovlp * ovlp
      return
      end
                     
      
C  This routine reads the input file 'ccc.in'
      subroutine readin(labot,latop,nabot,natop,lnabtop,nnbtop,energy,
     >   ntstart,ntstop,lstart,lstop,lttop,nbnd,npsbnd,albnd,alpha,
     >   ncstates,npot,lpot,nptop,lptop,formcut,regcut,expcut,theta,
     >   ifirst,isecond,nold,nq,qcut,rmax,fast,ldw,nk,sk,nznuc,nze,
     >   itail,gamma,r0,ninc,linc,ipar,nent,zasym,nunit,ndbl,nps,
     >   ne2e,lslow,lfast,slowe,enion,enlevel,target,projectile)
      include 'par.f'
      logical fast
      dimension nk(5,0:lmax), sk(5,0:lmax), nbnd(0:lmax),
     >   nabot(0:lamax), natop(0:lamax), alpha(0:lamax), nps(0:lamax),
     >   slowe(ncmax), enions(20), enlevels(20)
      character targets(20)*6,target*(*),projectile*(*)
      data targets/'H I','He II','Li I','Be II','B III','C IV','N V',
     >   'O VI','F VII','Ne I','Na I','Mg II','Al III','Si IV',
     >   'P I','S I','Cl I','A VIII','K I','Ca II'/
      
C  The following energy ionization levels in Rydbergs are from
C  ATOMIC ENERGY LEVELS, Volume I, C. E. Moore, NSRDS (1971)
C  The data below corresponds to the Z of nearest hydrogen-like target.
C                H      HeII   Li    BeII   BIII   CIV    NV     OVI
      dataenions/13.595,54.400,5.390,18.206,37.920,64.476,97.863,138.08,
C        FVII    Ne     Na    MgII  AlIII SiIV  P    S      Cl    AVIII
     >   185.139,21.559,5.138,15.03,28.44,45.13,11.,10.357,13.01,143.46,
C        K     CaII
     >   4.339,11.87/
C   The data below corresponds to the Z of neutral target.
c$$$  data enions/13.595,24.580,5.390,9.320,8.296,11.264,14.54,13.614,
c$$$     >   17.42,21.559,5.138,7.644,5.984,8.149,11.0,10.357,13.01,15.755,
c$$$     >   4.339,6.111/
C  The energy levels are 0 for the ground state, and the given values below
C  corresponding to the ionization energy above.
C  The data below corresponds to the Z of nearest hydrogen-like target.
C                   H          HeII      Li       BeII     BIII
      data enlevels/109678.758,438889.04,43487.19,146881.7,305931.1,
C        CIV      NV       OVI       FVII    Ne       Na
     >   520177.8,789532.9,1113999.5,1493656,173931.7,41449.65,
C        MgII      AlIII     SiIV     P       S       Cl       AVIII
     >   121267.41,229453.99,364097.7,88560.0,83559.3,104991.0,1157400.,
C        K        CaII
     >   35009.78,95748.0/
C  The data below corresponds to the Z of neutral target.
c$$$      data enlevels/109678.758,198305.0,43487.19,75192.29,66930.0,
c$$$     >   90878.3,117345.0,109836.7,140553.5,173931.7,41449.65,
c$$$     >   61669.14,48279.16,65743.0,88560.0,83559.3,104991.0,127109.9,
c$$$     >   35009.78,49304.8/

      open(3,file='ccc.in',status='old')

      read(3,*) energy,nze,natop(0),nps(0),alpha(0)
      print '('' energy,nze,natop,nps,alpha:                '',
     >   f7.3,3i7,f7.3)',
     >   energy,nze,natop(0),nps(0),alpha(0)
      if (natop(0).gt.nchan) stop 'NATOP must be <= NCHAN'
      if (nps(0).gt.ncmax) stop 'NPS must be <= NCMAX'
      ninc = 1
      linc = 0
      zasym = 0.0
      nznuc = 1
      if (nze.eq.-1) then
         projectile = 'electron'
      else if (nze.eq.1) then
         projectile = 'positron'
      else 
         print*,'NZE should be -1 for electron scattering'
         print*,'NZE should be +1 for positron scattering'
         print*,'Here NZE:',nze
         stop 'Wrong value for NZE'
      end if
      target = targets(nznuc)
      enlevel = enlevels(nznuc)
      enion = enions(nznuc)
      n = nznuc - nint(zasym)

      labot = 0
      latop = 0
      nabot(0) = 1

      do l = labot, latop
         if (natop(l).gt.nnmax) then
            print*,'Must have NATOP <= NNMAX',natop(l),nnmax
            stop 'Must have NATOP <= NNMAX'
         endif          
      end do 
      read(3,*) nunit,nnbtop,nent,ifirst,nold,theta
      print '('' nunit,nnbtop,nent,ifirst,nold,theta:'',
     >   5i7,f7.3)',nunit,nnbtop,nent,ifirst,nold,theta

      if (nnbtop.gt.nnmax) stop 'NNTOP must be <= NNMAX'
      ntst = 0
      lnabtop = 0
      lttop = 0
      ncstates = 0
      lstart = 0
      lstop = 0
      ipar = 0
      
      nptop = 0
      lptop = 0
      if (alpha(latop).eq.0.0) stop 'ALPHA(LATOP) can''t be zero'
      
      npot = 0
      lpot = 0
      npsbnd = 0
      albnd = 0.0
      formcut = 0.0
      regcut = 0.0
      expcut = 1e-10
C  FORMCUT cuts the form factors in FORM
C  REGCUT determines the minimum value at which regular solutions start
C  EXPCUT determines the smallest value of functions containing EXP fall off
      
      r0 = 1.0
      
      isecond = -1
      itail = 0
      if (nze.eq.1) then
         if (ifirst.gt.0) ifirst = 0
         if (isecond.gt.0) isecond = 0
      end if 
      if (nptop.eq.0.and.isecond.ge.0) then
         print*,'ISECOND >= 0 is inconsistent with NPTOP = 0'
         stop 'ISECOND >= 0 is inconsistent with NPTOP = 0'
      endif 

      ne2e = 0
      lslow = 0
      lfast = 0

      read(3,*) nq,qcut,rmax,slowe(1)
      print '('' nq,qcut,rmax,slowe:                 '',
     >   i7,3f7.1)',nq,qcut,rmax,slowe(1)
      fast = .false.
      ndbl = 0

      lp = 0
      mint = 4
      l = 0
      do while (lp.le.lstop+latop)
         read(3,*) (nk(j,l),sk(j,l),j=1, mint)
         nbnd(l) = 0
         if (zasym.eq.0.0.and.l.gt.ldw) nbnd(l) = 0
         if (l.gt.lmax) stop 'L > LMAX'
         do while (lp.lt.l)
            nbnd(lp) = nbnd(l)
            do j = 1, mint
               nk(j,lp) = nk(j,l)
               sk(j,lp) = sk(j,l)
            enddo
            lp = lp + 1
         enddo
         lp = l + 1
         print '('' (nk(j),sk(j),j=1,4)  '',
     >      1x,4(i7,f7.2)///)', (nk(j,l),sk(j,l),j=1,mint)
         if (mod(nk(4,l),2).ne.0) stop 'NK(4,L) must be even'
      enddo 
      CLOSE(3)

      if (lttop.gt.ltmax) then
         print*,'LTTOP > LTMAX',lttop,ltmax
         print*,'Recompile with larger LTMAX'
         stop 'ABORTED'
      end if
      if (nq.gt.kmax) then
         print*,'Increase KMAX in par.f'
         stop 'Increase KMAX in par.f'
      endif 
      return
      end      
      
      subroutine arrange
c
c  this subroutine sets up the faclog array.
c
      include 'par.f'
      double precision faclog
      common /flogs/ faclog(1000)
      faclog(1)=0d0
      faclog(2)=0d0
      do 10 n=3,1000
   10 faclog(n)=faclog(n-1)+log(dfloat(n-1))

      return
      end
      
      subroutine setpow(lstop,lttop)
C  this subroutine sets up the arrays RPOW1=R**LNA and RPOW2=1/R**(LNA+1)
C  which are used in FORM, as well as CNTFUG=L(L+1)/XMESH(J)**2 which 
C  is used in Numerov integration.
      include 'par.f'
      common/meshrr/ meshr,rmesh(maxr,3)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   istartrp(0:ltmax),istoprp(0:ltmax),cntfug(maxr,0:lmax)
      common/smallr/ formcut,regcut,expcut,fast
      logical fast

      if (lmax.lt.lstop) then
         print*,'LSTOP, LMAX:',lstop, lmax
         stop 'Recompile with larger LMAX'
      endif 
      do j=1,meshr
         cntfug(j,0)=0.0
         x=rmesh(j,1)
         x2=x*x
         cntfug(j,1)=2.0/x2
      end do
      do l=2,lstop+lttop
         do j=1,meshr
            cntfug(j,l)=cntfug(j,l-1)*float(l+1)/float(l-1)
         end do
      end do

      istartrp(0)=1
      istoprp(0)=meshr
      do i=1,meshr
         rpow1(i,0)=1.0
         rpow2(i,0)=1.0/rmesh(i,1)
      end do
      do lna=1,lttop
         istartrp(lna)=istartrp(lna-1)
         istoprp(lna)=istoprp(lna-1)
         do i=istartrp(lna),istoprp(lna)
            rpow1(i,lna)=rpow1(i,lna-1)*rmesh(i,1)
            rpow2(i,lna)=rpow2(i,lna-1)/rmesh(i,1)
         end do
         do while (rpow1(istartrp(lna),lna).lt.regcut)
            istartrp(lna)=istartrp(lna)+1
         end do
         do while (rpow2(istoprp(lna),lna).lt.expcut*expcut)
            istoprp(lna)=istoprp(lna)-1
         end do
      end do 
      end
         
C
C======================================================================      
C  This subroutine sets up two grids: gridx and gridr. Gridx will be
C  used to solve the radial Schrodinger equation. Gridr, which is
C  typically every second point of gridx, is used to perform integrations
C  from zero to infinity using Simpson's rule. The structure is such
C  that a wave with momentum QMAX will have NPWAVE points per oscillation. This
C  determines HMAX, the largest step. At some point, XDBLE, the intervals,
C  dX, are progressively halved. This is done to accomodate irregular
C  solutions.
C  INPUT:
C   qmax  - the biggest momentum (a.u.) that can be calculated by this mesh.
C   rmax  - the largest "r=x" in the meshes.
C   nmaxx - array declaration parameter of first declaration of GRIDX
C   nmaxr - array declaration parameter of first declaration of GRIDR
C   nmaxj - array declaration parameter of first declaration of JDOUBLE
C  OUTPUT:
C   gridx - X grid
C   nx    - Number of X points
C   gridr - R grid
C   nr    - number of R points
C   jdouble - j points where dx doubles
C   njdouble - Number of doublings + 2
C======================================================================      
      subroutine grids(qmax,ndouble,rmax,gridr,nmaxr,nr,jdouble,
     >   nmaxj,njdouble)
C
C
      dimension gridr(nmaxr,3), jdouble(nmaxj)
C  jdouble stores the J index for GRIDR where DR doubles, with the first
C  point set to 1 and the last to NR. Used in the numerov integrations.

C  Set up GRIDR
C  The number of points per oscillation (half-wavelength): NPWAVE
      npwave = 6
C  The number of doubling of intervals is NDOUBLE
      njdouble = ndouble + 2
      if (njdouble.gt.nmaxj) then
         print*,'Increase NMAXJ in call to GRIDS to at least:',njdouble
         stop 'Increase NMAXJ in call to GRIDS'
      end if 
C  let NPDBL be the number of points with the same dx per interval
      npdbl = 40
      npdbl = 32

C  Make sure NPDBL is even
      npdbl=(npdbl/2) * 2
      if (npdbl.lt.4) then
         print*,'Need to have at least 4 equally spaced points in GRIDR'
         stop 'Need to have at least 4 equally spaced points in GRIDR'
      end if 
      hmax = 3.14/float(npwave)/qmax
      rdble = float(npdbl) * hmax * (2**ndouble - 1) / float(2**ndouble)
      rleft = rmax - rdble
      nleft = int(rleft / hmax) / 2 * 2
      ntot = nleft + npdbl * ndouble
c$$$      print*,'Estimated R max:',rdble + nleft * hmax,ntot,
c$$$     >   hmax * (float(npdbl) * (2**ndouble - 1) / float(2**ndouble)
c$$$     >   + nleft)
      if (ntot.le.nmaxr) then
c$$$         print*,'Old hmax:',hmax
         hmax = rmax / (float(npdbl) * (2**ndouble - 1) /
     >      float(2**ndouble) + nleft)
c$$$         print*,'New hmax:',hmax
      endif 
C  2**ndouble * hmin = hmax, therefore
      hmin = hmax/float(2**ndouble)
C  The value of the R point from which dR is constant, = HMAX, is RDBLE
C  RDBLE = NPDBL * hmin * (2**NDOUBLE-1)
      rdble = float(npdbl) * hmin * (2**ndouble-1)

c$$$      print*,'Grid R parameters:'
c$$$      print*,'NDOUBLE:',ndouble
c$$$      print*,'HMIN:',hmin
c$$$      print*,'HMAX:',hmax
c$$$      print*,'NPDBL:',npdbl
c$$$      print*,'RDBLE:',rdble
      
      dr = hmin
      jdouble(1) = 1
      do nd = 2, ndouble + 1
         jdouble(nd) = npdbl * (nd - 1)
      end do 

      r=0.0
      j = 0
      do nd = 1, ndouble
         do nj = 1, npdbl
            j = j + 1
            gridr(j,1) = r + float(nj) * dr
            gridr(j,2) = dr
C  Simpson's rule weights
            gridr(j,3) = float(mod(j,2) * 2 + 2) * dr / 3.0
         end do
         gridr(j,3) = dr
         r = gridr(j,1)
         dr = dr * 2.0
      end do

      if (abs((dr - hmax)/hmax).gt.1e-4) then
         print*,'DR and HMAX should be equal'
         stop 'DR and HMAX should be equal'
      end if 
      nj = nint((rmax-r)/hmax) / 2 * 2
c$$$      nj = nj / 32
c$$$      nj = nj * 32
      if (nj+j.gt.nmaxr) then
         nj = nmaxr - j
         print*,'Warning: NR had to be reduced to NMAXR'
      end if 
      do jj = 1, nj
         j = j + 1
         gridr(j,1) = r + float(jj) * dr
         gridr(j,2) = dr
C  Simpson's rule weights
         gridr(j,3) = float(mod(j,2) * 2 + 2) * dr / 3.0
      end do
C  Set the last weight so that the integrals ended exactly at rmax.
C  This is done so that the tail integrals could be done correctly.
      gridr(j,3) = dr / 3.0
      nr = j
      jdouble(ndouble+2) = j
      if (mod(j,2).ne.0)print*,'Warning expected J to be even in grids'

      do n = 0,0
         if (n.eq.0) then
            s = 2.0
         else
            s = 3.0
         endif
         w = 1.0
         test = 0.0
         do j = 2**n, nr, 2**n
            w = s - w
            test = test + gridr(j,1) ** 2 * gridr(j,3) * w
         enddo
         if (j-2**n.ne.nr)print*,'Warning last J should be NR',j-2**n,nr
         test = test * 2**n
c$$$         print*,2**n,test*3.0/gridr(nr,1)**3
      enddo 
      print*,'DR, RMAX and NR:', gridr(1,1), gridr(nr,1), nr

      return
      end

      subroutine improvek(e,sk1,nt,gridk,weightk,endf,endk)
      real gridk(nt),weightk(nt)

      endf = sk1
      if (e.lt.0.0) return
      rk = sqrt(e)
      if (sk1.lt.rk) then
         endf = (2.0 * rk - endk)
         if (endf.lt.0.0) stop 'Problem in IMPOVEK; ENDK too big'
         return
      endif 
      k = 1
      do while (gridk(k) .lt. rk)
         k = k + 1
      enddo
      call getint(e,sk1,sk1,gridk,weightk,nt,sum)
      if (sum.gt.0.0) then
         sk1max = sk1
         sk1min = sk1 * (9.0 * rk + gridk(k)) / 10.0 / gridk(k)
         call getint(e,sk1,sk1min,gridk,weightk,nt,sum)
      else 
         sk1min = sk1
         sk1max = sk1 * (9.0 * rk + gridk(k-1)) / 10.0 / gridk(k-1)
         call getint(e,sk1,sk1max,gridk,weightk,nt,sum)
      endif 

      ncount = 0
      do while(abs((sk1max-sk1min)/sk1min).gt.1e-8.and.ncount.lt.50)
         ncount = ncount + 1
         endf = (sk1max + sk1min) / 2.0
         call getint(e,sk1,endf,gridk,weightk,nt,sum)
         if (sum.gt.0.0) then
            sk1max = endf
         else
            sk1min = endf
         endif 
         if (ncount.ge.90) print*,sk1min,sk1max, sum
      enddo 
      end

      subroutine getint(e,sk1,sk1p,gridk,weightk,nt,sum)
      real gridk(nt),weightk(nt)
      r = sk1p / sk1
      sum = 0.0
      do i = 1, nt
         sum = sum + 2.0 * weightk(i) * r / (e - (gridk(i) * r)**2)
      end do 
      rk = sqrt(e)
      tsum = - 2.0 * acoth(sk1p/rk) / rk
      sum = sum + tsum
      end

C  Define the k-grid integration intervals and the number of points in
C  each interval
      subroutine makeints(mint,sk,nk,rk,np,width,midnp,endk,nendk,endp,
     >   np2,endk2)
      dimension sk(0:10), nk(10)
      nmin(npoints) = min(npoints/2,4)
      
      sk(0) = 0.0
      if (abs(width).lt.0.01) then
         print*,'Width must be greater than 0.01'
         stop 'Width must be greater than 0.01'
      endif
      if (endk+2.0*width.ge.endk2) then
         print*,'ENDK2 must be >= ENDK + 2 * WIDTH',endk2,endk+2.*width
         stop 'ENDK2 must be >= ENDK + 2 * WIDTH'
      endif 
      if (rk.lt.0.0.or.(width.lt.0.0.and.rk.lt.2.0*abs(width))) then
C  Here for closed channels
c$$$         nk(1) = midnp
c$$$         sk(1) = 2.0 * abs(width)
c$$$         nk(2) = np
         nk(1) = np
         sk(1) = endk
         nk(2) = np2
         sk(2) = endk2         
         mint = 3
      else if (rk.lt.width) then
C  Here if singularity is going to be in the first interval
         sk(1) = rk * 2.0
         nk(1) = midnp
         sk(2) = endk
         nk(2) = np
         sk(3) = endk2
         nk(3) = np2
         mint = 4
      else if (rk. lt. endk + 2.0 * width) then
C  Here if singularity is going to be in the second interval
         sk(1) = rk - width
         sk(2) = rk + width

C  Use the following to make the k grid symmetric about rk**2
c$$$            sk(1) = sqrt(rk**2 - 2.0 * rk * width)
c$$$            sk(2) = sqrt(rk**2 + 2.0 * rk * width)

         nk(2) = midnp
         if (rk .lt. endk - 2.0 * width ) then
            nk(1) = max(nmin(np),int(np * rk / (endk - 2.0 * width)))
            nk(1) = min(nk(1), np - nmin(np))
            npleft = np - nk(1)
            sk(3) = endk
            nk(3) = npleft
            sk(4) = endk2
            nk(4) = np2
            mint = 5
         else
            nk(1) = np
            sk(3) = endk2
            nk(3) = np2
            mint = 4
         endif 
      else if (rk. lt. endk2) then
         sk(1) = endk
         nk(1) = np
         sk(2) = rk - width
         sk(3) = rk + width
         nk(3) = midnp
C  The two factors of 3.0 below put extra points in the first interval
c$$$         nk(2) = max(4,int(np2 * (sk(2)-sk(1)) * 3.0 /
c$$$     >      (3.0 * (sk(2) - sk(1)) + (endk2-endk))))
         nk(2) = max(nmin(np2),int(np2 * (rk - endk) / (endk2-endk)))
         nk(2) = min(nk(2), np2 - nmin(np2))
         npleft = np2 - nk(2)
         if (rk .lt. endk2 - width) then
            sk(4) = endk2
            nk(4) = npleft
            mint = 5
         else
            mint = 4
         endif    
      else
         sk(1) = endk
         nk(1) = np
         sk(2) = rk - width
         nk(2) = np2
         sk(3) = rk + width
         nk(3) = midnp
         mint = 4
      endif
      nk(mint) = nendk
      sk(mint) = endp
      return
      end

C  Doesn't do much
      subroutine pwrite(nopen,energy,nznuc,zasym,ry,noprint,ovlp,
     >   jstart,jmax,projectile,target)
      include 'par.f'
      character projectile*8,target*6
      real energy,ry,ovlp(ncmax,0:lamax)
      dimension lea(nchan),exv(nchan),onshellk(nchan)
      common /pspace/ nabot,labot,natop,latop,ntype,ipar,nze,ninc,linc
      dimension nabot(0:lamax), natop(0:lamax), temp(maxr)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      logical noprint(nchan)
      
      npos = 0
      if (nze.eq.1) npos = 1
      zassymp = zasym

      if (ipar.eq.0) then
         jch = 0
      else
         jch = 1
      endif 
      nch = 1
      call getchinfo (nch, nchp, jch, temp, maxpsi, ea, la, na, l)
      nopen = 0
 10   continue 
      echan = (enpsinb(na,la)-enpsinb(ninc,linc)) * ry
      if (energy.ge.echan) then
         nopen = nopen + 1
         onshellk(nopen) = sqrt(energy/ry - echan/ry)
         lea(nopen) = 2 * la
C  The following line puts the overlap of a pseudostate with the I^-
C  projection operator in to EXV which will be used in the CROSS program
         exv(nopen) = ovlp(na,la)
C  The following change was done for the unnatural parity
c$$$         noprint(nch) = .true.
         if (nchp.ne.nch) print*,'Caution: NCH, NCHP in PWRITE',
     >      nch, nchp
         noprint(nchp) = .true.
      else
C  Need to check that the following line is correct
         onshellk(nch) = 0.0
C  The following change was done for the unnatural parity
         noprint(nchp) = .false.
c$$$         noprint(nch) = .false.
      end if
c$$$      exv(nch) = ovlp(na,la)
      nchtop = nch
      nch = nch + 1
      call getchinfo (nch, nchp, jch, temp, maxpsi, ea, la, na, l)
      if (nch.ne.0) go to 10

            if (npos.eq.0) then
         nsp = 2
      else
         nsp = 1
      endif
      end

C  Write out total cross sections to the totalcs file
      subroutine wrtcs(partcs,sigtop,nchpmax,nchipmax,lg,etot,nsmax,
     >   ovlp,nunit,nznuc,zasym,projectile,target)
      include 'par.f'
      dimension partcs(nchan,nchan,0:1), oldpj(nchan,nchan,0:1),
     >   oldp(nchan,nchan,0:1), sigionold(nchan,0:1), sigb(nchan,0:1), 
     >   sigion(nchan,0:1), sigt(nchan,0:1), sigtop(nchan,0:1),
     >   sigtopold(nchan,0:1), sigtopoldj(nchan,0:1),
     >   sigione(nchan,0:1), sigtope(nchan,0:1),
     >   sigionl(nchan,0:lamax,0:1),sigtl(nchan,0:lamax,0:1),
     >   sigbl(nchan,0:lamax,0:1)
      real ovlp(ncmax,0:lamax)
      logical canstop
      character*8 chunit(3)
      character*3 chan(100),chs(0:lamax)
      character projectile*(*), target*(*)
      common /charchan/ chan
      common /chanens/ enchan(100)
      asym(sing,trip,fac) = (sing - trip / fac) / (sing + trip + 1e-30)
      
      data pi,chunit/3.141593,' a0^2   ',' pi a0^2',' cm^2  '/

      nwf=2*ncmax
      niwf=2*ncmax

      fac = 2.0
      if (chan(1)(1:1).eq.' ') fac = 3.0
      small = 1e-3
      canstop = .true.
      if (lg.eq.0) canstop = .false.
C  Set default to a0**2
      unit = 1.0
      if (nunit.eq.2) then
         unit = 1.0 / pi
      elseif (nunit.eq.3) then
         unit = 5.2917 ** 2 * 1e-18
      else
         nunit = 1
      endif

C  The following Do loops are unnecessary, I think!
      do ns = 0, 1
         do nchip = 1, nchipmax
            sigtopold(nchip,ns) = 0.0
            do nchp = 1, nchpmax
               oldp(nchp,nchip,ns) = 0.0
            enddo
         enddo
      enddo

C  We redefine NCHPMAX and NCHIPMAX here because on occasion NCHPMAX for
C  the unnatural parity case comes out to be smaller than in the natural
C  parity case
      open(42,file='totalcs',status='unknown')
      write(42,'(a8,'' - '',a6,10(f7.2,''eV on '',a3))')
     >   projectile,target,(13.6058 * (etot-enchan(i)),chan(i),
     >   i=1,nchipmax)
      write(42,'(''Units:'',a8,'' spin-weights included'')')
     >   chunit(nunit)
      write(42,'(79a)') ('-',i=1,39)
      do ns = 0, nsmax
         do nchip = 1, nchipmax
            sigt(nchip,ns) = 0.0
            sigb(nchip,ns) = 0.0
            do l = 0, lamax
               sigionl(nchip,l,ns) = 0.0
               sigtl(nchip,l,ns) = 0.0
               sigbl(nchip,l,ns) = 0.0
            enddo
            ltop = 0
            do nchp = 1, nchpmax
C  SIGT will be the total cross section for all J, calculated by summing 
C  individual cross sections.
               sig = partcs(nchp,nchip,ns) * unit + oldp(nchp,nchip,ns)
               sigt(nchip,ns) = sigt(nchip,ns) + sig
               call getchnl(chan(nchp),n,l)
               if (l.gt.ltop) ltop = l
               sigtl(nchip,l,ns) = sigtl(nchip,l,ns) + sig

               if (enchan(nchp) .gt. 0.0) then
                  sigionl(nchip,l,ns) = sigionl(nchip,l,ns) + sig
               else
                  sigbl(nchip,l,ns) = sigbl(nchip,l,ns) + sig*ovlp(n,l)
                  sigb(nchip,ns) = sigb(nchip,ns) + sig * ovlp(n,l)
               endif 
            enddo

C  Redefine the ionization cross sections by SIGT - SIGB
            do l = 0, ltop
               sigionl(nchip,l,ns) = max(0.0,
     >            sigtl(nchip,l,ns)-sigbl(nchip,l,ns))
            enddo
            
            
            sum = 0.0
            sumold = 0.0
            sume = 0.0
            sumo = 0.0
            do nchp = 1, nchpmax
               sig = partcs(nchp,nchip,ns) * unit + oldp(nchp,nchip,ns)
               call getchnl(chan(nchp),n,l)
c$$$               if (enchan(nchp) .gt. 0.0) ovlp(n,l) = 0.0
               if (enchan(nchp) .lt. 1e-10) then
C  SUM will be the total non-break-up cross section for all J
                  sum = sum + ovlp(n,l) * sig
C  SUMOLD will be the total non-break-up cross section for all previous J
                  sumold = sumold + ovlp(n,l) * oldp(nchp,nchip,ns)
C  SUMO will be the total non-break-up cross section for the current J only
                  sumo = sumo + ovlp(n,l) * partcs(nchp,nchip,ns) * unit
C  SUME will be the total non-break-up cross section for the previous J only
                  sume = sume + ovlp(n,l) * oldpj(nchp,nchip,ns)
                  if (ovlp(n,l).gt.1.001) print*,
     >               'overlap too big, should be <= 1.0',n,l,ovlp(n,l)
               endif 
            enddo
C  Use the optical theorem to define the total and ionization cross sections
C  as then the code is suitable for both CCC and CCO. For CCC both 
C  forms should give much the same answer.
C  SIGION will be the total ionization cross section for all J
            sigion(nchip,ns) = max(sigtop(nchip,ns) * unit
     >         + sigtopold(nchip,ns) - sum,0.0)
C  SIGIONOLD will be the total ionization cross section for all previous J
            sigionold(nchip,ns) = max(sigtopold(nchip,ns) - sumold,0.0)
C  SIGIONE will be the extrapolated total ionization cross section
            sigione(nchip,ns) = sigionold(nchip,ns) +
     >         extrap(sigtop(nchip,ns) * unit - sumo,
     >         sigtopoldj(nchip,ns) - sume)
C  SIGTOPE will be the extrapolated total cross section
            sigtope(nchip,ns) = sigtopold(nchip,ns) +
     >         extrap(sigtop(nchip,ns) * unit,sigtopoldj(nchip,ns))
C  For pure Born approximation we get zero for the optical theorem. In this
C  case use the form below. Can't define sigionold in this case, need more work
C  Need to make SIGTOLD from the info above
            if (sigtop(nchip,ns).eq.0.0) then
               sigion(nchip,ns) = max(sigt(nchip,ns) - sum,0.0)
               sigionold(nchip,ns) = 0.0
            endif 
         enddo
      enddo
      
      do nchip = 1, nchipmax
         sumion = 0.0
         do l = 0, ltop
            chs(l) = ' + '
            if (l.eq.0) chs(l) = ' = '
            do ns = 0, nsmax
               sumion = sumion + sigionl(nchip,l,ns)
            enddo
         enddo 
         do ns = 0, nsmax
            write(42,'(f7.2,''eV on '',a3,'' TNBCS(S='',i1,''): '',
     >         1p,e8.2,9(a3,e8.2))') 13.6058 * (etot-enchan(nchip)),
     >         chan(nchip),ns,sigb(nchip,ns)
            write(42,'(f7.2,''eV on '',a3,''  TICS(S='',i1,''): '',
     >         1p,e8.2,9(a3,e8.2))') 13.6058 * (etot-enchan(nchip)),
     >         chan(nchip),ns,sigion(nchip,ns)
            write(42,'(f7.2,''eV on '',a3,''   TCS(S='',i1,''): '',
     >         1p,e8.2,9(a3,e8.2))') 13.6058 * (etot-enchan(nchip)),
     >         chan(nchip),ns,sigt(nchip,ns)
            write(42,'(''transition  cross section '',
     >         '' for S ='',i2)') ns
            do nchp = 1, nchpmax
               call getchnl(chan(nchp),n,l)
               write(42,'(a3,'' <-'',a3,1p,3e15.5,0p,f10.4)') 
     >            chan(nchp),chan(nchip),partcs(nchp,nchip,ns)*unit
            enddo
            write(42,'(79a)') ('-',i=1,39)
         enddo 
         write(42,'(f7.2,''eV on '',a3,'' TNBCS:      '',
     >      1p,e8.2,9(a3,e8.2))') 13.6058 * (etot-enchan(nchip)),
     >      chan(nchip),sigb(nchip,0)+sigb(nchip,1)
         write(42,'(f7.2,''eV on '',a3,''  TICS:      '',
     >      1p,e8.2,9(a3,e8.2))') 13.6058 * (etot-enchan(nchip)),
     >      chan(nchip),sigion(nchip,0)+sigion(nchip,1)
         write(42,'(f7.2,''eV on '',a3,''   TCS:      '',
     >      1p,e8.2,9(a3,e8.2))') 13.6058 * (etot-enchan(nchip)),
     >      chan(nchip),sigt(nchip,0)+sigt(nchip,1)
         write(42,'(79a)') ('-',i=1,39)
            
         write(42,'(f7.2,''eV on '',a3,
     >      '' TICS, spin asymmetry:'',
     >      1p,4e11.3)') 13.6058 * (etot-enchan(nchip)),
     >      chan(nchip),sigion(nchip,0) + sigion(nchip,1),
     >      asym(sigion(nchip,0), sigion(nchip,1),fac)
         
         sigtsum = sigt(nchip,0) + sigt(nchip,1)
         sigtopt = (sigtop(nchip,0) + sigtop(nchip,1)) * unit +
     >      sigtopold(nchip,0) + sigtopold(nchip,1)
         write(42,'(f7.2,''eV on '',a3,
     >      '' TCS (op. th.), asym: '',
     >      1p,3e11.3)') 13.6058 * (etot-enchan(nchip)),chan(nchip),
     >      sigtopt, 
     >      asym(sigtope(nchip,0),sigtope(nchip,1),fac)
         write(42,'(79a)') ('-',i=1,65)
         write(42,'(''transition  cross section '',
     >      ''   projection    spin asym      energy'')')
         do nchp = 1, nchpmax
            call getchnl(chan(nchp),n,l)
            summedcs = partcs(nchp,nchip,0)*unit + oldp(nchp,nchip,0) +
     >         partcs(nchp,nchip,1) * unit + oldp(nchp,nchip,1)
            asymcs = asym(partcs(nchp,nchip,0)*unit+oldp(nchp,nchip,0),
     >         partcs(nchp,nchip,1) * unit + oldp(nchp,nchip,1), fac)
            extrapcs = oldp(nchp,nchip,0) + oldp(nchp,nchip,1) +
     >         extrap((partcs(nchp,nchip,0)+partcs(nchp,nchip,1))*unit,
     >         oldpj(nchp,nchip,0) + oldpj(nchp,nchip,1))
            canstop = canstop.and.(extrapcs-summedcs)/(summedcs+1e-30).
     >         lt.small.and.extrapcs.ge.summedcs
            write(42,'(a3,'' <-'',a3,1p,3e15.5,0p,f10.4)') chan(nchp),
     >         chan(nchip),summedcs,  ovlp(n,l), asymcs,
     >         enchan(nchp)
         enddo
      enddo 
      close(42)
      if (canstop) then
         print*,'Stopping since have convergence for all states at J:',
     >      lg
         stop 'Have convergence for all states'
      endif 
      end
      
      function extrap(x,xp)
      extrap = x
      end
      
      subroutine trapez(dstart,dstop,nt,x,w)
      implicit double precision (a-h,o-z)
      dimension x(nt),w(nt)

      h = (dstop-dstart) / dfloat(nt - 1)
      do n = 1, nt
         x(n) = dstart + dfloat(n-1) * h
         w(n) = h
      enddo
      w(1) = h / 2d0
      w(nt) = h / 2d0
      if (abs(dstart).lt.1d-10) then
         x(1) = h / 2d0
         w(1) = h / 2d0
         w(2) = 0.75d0 * h
      endif
      sum = 0d0
      do n = 1, nt
         sum = sum + x(n) * w(n)
      enddo
      exact = (dstop**2 - dstart**2) / 2d0
      if (abs(exact/sum - 1d0).gt.1d-10) print*,'problems in TRAPEZ',
     >   exact/sum
      return
      end
      
      subroutine matinv2(api,n1,nmax,bhar,nb,wk,erfp,epsil)
      include 'par.f'
      complex api,bhar,wk
      real erfp,epsil
      dimension api(n1,n1),bhar(n1,nchan),wk(n1),kpvt(kmax*nchan)
      
      lda = n1
      n = nmax
      call cgesv(n,nb,api,lda,kpvt,bhar,lda,info)
      return
      end

      
C  This routine returns  FORM(R) = int dR' FUN(R') * F(R<) * G(R>), where
C  the range of integration is zero to infinity.
      subroutine form(fun,ifuns1,ifuns2,f,g,ifstart,igstop,irstop,
     >   formf,i1,i2)
      include 'par.f'
      dimension fun(maxr),formf(maxr),f(maxr),g(maxr),temp1(maxr),
     >   temp4(0:maxr),temp5(maxr+1),temp2(maxr),temp3(maxr)
      common/meshrr/ meshr,rmesh(maxr,3)
      common/smallr/ formcut,regcut,expcut,fast
      logical fast
      
      ifunstart=max(ifuns1,ifstart)
      ifunstop=min(ifuns2,igstop)
      if (ifunstop.le.ifunstart+2) then
         i1=2
         i2=1
         return
      end if 
      istop=min(ifunstop,irstop)
      if (fast) then
C  Find the integral for FUN(R') * F(R') from 0 to IRSTOP. The function FUN 
C  already contains the Simson's integration weights. 
         temp4(ifunstart-1)=0.0
         do i=ifunstart,istop
            temp4(i)=temp4(i-1)+fun(i)*f(i)
         end do 
C  Find the integral of FUN(R') * G(R') from infinity to IFUNSTART
         temp5(ifunstop+1)=0.0
         do i=ifunstop,ifunstart,-1
            temp5(i)=temp5(i+1)+fun(i)*g(i)
         end do 

C  Make the form factor
         do i=ifunstart,istop
            formf(i)=temp4(i)*g(i)+temp5(i+1)*f(i)
         end do
      else
C  This method is more accurate but takes longer.
         do i=ifunstart,istop
            temp1(i)=fun(i)*f(i)/rmesh(i,3)
         end do
         call maketemp(ifunstart,istop,temp1,temp2)
         do i=ifunstart,ifunstop
            temp1(i)=fun(i)*g(i)/rmesh(i,3)
         end do
         call maketempb(ifunstart,ifunstop,temp1,temp3)
         do i=ifunstart,istop
            formf(i)=temp2(i)*g(i)+temp3(i)*f(i)
         end do
         temp4(istop)=temp2(istop)
         temp5(ifunstart)=temp3(ifunstart)
      end if 
            
C  IRSTOP may be greater than IFUNSTOP
      i=istop
      do while (abs(formf(i)).gt.formcut.and.i.lt.min(irstop,igstop))
         i=i+1
         formf(i)=temp4(istop)*g(i)
      end do 
      i2=i
      do while (abs(formf(i2)).lt.formcut.and.i2.gt.ifunstart)
         i2=i2-1
      end do 

      i=ifunstart
      do while (i.gt.ifstart.and.abs(formf(i)).gt.formcut)
         i=i-1
         formf(i)=f(i)*temp5(ifunstart)
      end do
      i1=i
      do while (abs(formf(i1)).lt.formcut.and.i1.le.i2)
         i1=i1+1
      end do 
      end

      subroutine nuclear(fun,direct,ifuns1,ifuns2,irstop,u,nz,form)
      include 'par.f'
      dimension fun(maxr),form(maxr),u(maxr)
      common/meshrr/ meshr,rmesh(maxr,3)
      common /powers/ rpow1(maxr, 0:ltmax), rpow2(maxr, 0:ltmax),
     >   iminrp(0:ltmax), imaxrp(0:ltmax), cntfug(maxr, 0:lmax)
      common/smallr/ formcut,regcut,expcut,fast
      common /pspace/ nabot,labot,natop,latop,ntype,ipar,nze,ninc,linc
      dimension nabot(0:lamax), natop(0:lamax)

      logical direct
      z=float(nz)
      tmp=0.0
      if (direct) then
         do i=ifuns1,ifuns2
            tmp=tmp+fun(i)
         end do
         do i=1,irstop
            form(i) = form(i) - tmp * (u(i) / 2.0 + rpow2(i,0))
         end do
         do while (abs(form(irstop)).lt.formcut.and.irstop.gt.1)
            form(irstop) = 0.0
            irstop = irstop - 1
         end do 
      else
C  Exchange
         alpha = log(sqrt(expcut)) / rmesh(meshr,1) * 0.0
         do i=ifuns1,ifuns2
            tmp=tmp+fun(i)*(rpow2(i,0)+u(i)/2.0)
         end do
         do i=1,irstop
            form(i)=(form(i) - tmp) * exp(alpha * rmesh(i,1))
         end do
      end if
      end
      
C  The following routine returns the integral of FUN(I) from 
C  RMESH(I) to RMESH(ISTOP) for I=ISTART,ISTOP.
      subroutine maketempb(iistart,iistop,fun,temp)
      include 'par.f'
      common/meshrr/ meshr,rmesh(maxr,3)
      common /double/id,jdouble(20)
      dimension fun(maxr),temp(maxr)

      istop=iistop/2*2
      istart=(iistart+1)/2*2
      do i=1,meshr
         temp(i)=0.0
      end do
      if (istop-istart.lt.4) return
      
C  Define the integral at even points of RMESH using Simpson's rule.
C  Define the integral at odd points of RMESH using a combination of
C  Simpson's rule and Newton's 3/8 rule.
      temp(istop)=0.0
      do i=istop,istart+2,-2
         temp(i-2)=temp(i)+(fun(i-2)+4.*fun(i-1)+fun(i))*rmesh(i-1,2)/3.
         temp(i-3)=temp(i)+(fun(i-3)+3.0*fun(i-2)+3.0*fun(i-1)+fun(i))*
     >      0.375*rmesh(i,2)
      end do

C  The above Newton's 3/8 rule fails where dR doubles
C  I is the next point before the doubling. It should be odd.
      do n=2,id-1
         I=jdouble(n)-1
         temp(i)=temp(i-2)-(fun(i-2)+4.0*fun(i-1)+fun(i))*rmesh(i,2)/3.0
         if (i.eq.istart+1) temp(i)=temp(istart)
      end do

C  The following point had been left out
      temp(istop-1)=temp(istop-3)-(fun(istop-3)+4.0*fun(istop-2)+
     >   fun(istop-1))*rmesh(istop-1,2)/3.0

C  FUN(R) from R=RMESH(istop,1) to 'infinity' is assumed to be zero.
      do i=istop,meshr
         temp(i)=0.0
      end do

      if (iistart.gt.1) then
C  FUN(R) from R=0 to RMESH(istart-1,1) is assumed to be zero.
         do i=1,istart-1
            temp(i)=temp(istart)
         end do
      end if
      end

C  The following routine returns the integral of FUN(RP) from 0 to
C  RMESH(IRP) for IRP=1,IRPSTOP.
      subroutine maketemp(istart,istop,fun,temp)
      include 'par.f'
      common/meshrr/ meshr,rmesh(maxr,3)
      common /double/id,jdouble(20)
      dimension fun(maxr),temp(maxr)
      irpstart=(istart+1)/2*2
      irpstop=istop/2*2
      if (irpstart.eq.istart) fun(istart-1)=0.0
      
      do i=1,meshr
         temp(i)=0.0
      end do
      if (irpstop-irpstart.lt.4) return
         
      if (istart.eq.1) then
C  We assume that fun(0)=0. 
         temp(2)=(4.0*fun(1)+fun(2))*rmesh(1,2)/3.0
         temp(3)=(3.0*fun(1)+3.0*fun(2)+fun(3))*rmesh(1,2)*0.375
         temp(1)=temp(3)-(fun(1)+4.0*fun(2)+fun(3))*rmesh(1,2)/3.0
      else
         temp(irpstart)=(4.0*fun(irpstart-1)+fun(irpstart))*
     >      rmesh(irpstart-1,2)/3.0
         temp(irpstart+1)=(2.0*fun(irpstart-1)+
     >      4.0*fun(irpstart)+fun(irpstart+1))*rmesh(irpstart,2)/3.0
      end if
      
C  Define the integral at even points of RMESH using Simpson's rule.
C  Define the integral at odd points of RMESH using a combination of
C  Simpson's rule and Newton's 3/8 rule.
      do irp=irpstart+2,irpstop-2,2
         temp(irp)=temp(irp-2)+(fun(irp-2)+4.0*fun(irp-1)+fun(irp))
     >      *rmesh(irp-1,2)/3.0
         temp(irp+1)=temp(irp-2)+(fun(irp-2)+3.0*fun(irp-1)+
     >      3.0*fun(irp)+fun(irp+1))*0.375*rmesh(irp,2)
      end do

C  The above Newton's 3/8 rule fails where dR doubles.
C  I is the next point after the doubling. It should be odd. 
      do n=2,id-1
         i=jdouble(n)+1
         temp(i)=temp(i+2)-(fun(i)+4.0*fun(i+1)+fun(i+2))*rmesh(i,2)/3.0
         if (i.eq.irpstop-1) temp(i)=temp(i-1)
      end do

      if (istart.gt.1) then
C  FUN(R) from R=0 to RMESH(irpstart,1) is assumed to be zero.
         do irp=1,irpstart
            temp(irp)=0.0
         end do
      end if
C  FUN(R) from R=RMESH(irpstop,1) to 'infinity' is assumed to be zero.
      do i=irpstop,meshr
         temp(i)=temp(irpstop-1)
      end do
      end

      function fact(k)
      fact=1.0
      do i=2,k
         fact=fact*float(i)
      end do
      end

C  The following routine returns the hydrogenic radial functions.
      subroutine rnl(nznuc,n,l,chi,en,jstop)
      include 'par.f'
      common/meshrr/ meshr,rmesh(maxr,3)
      common/smallr/ formcut,regcut,expcut,fast
      logical fast
      dimension chi(maxr)
      Z = float(nznuc)
      do j=1,meshr
         chi(j)=0.0
      end do 
      const1 = sqrt(Z * fact(n-l-1) * fact(n+l)) / float(n)
      rho = 0.0
      exprho = exp(- rho / 2.0)
      j = 0
      do while (exprho.gt.expcut.and.j.lt.meshr)
         j = j + 1
         rlrho = 0.0
         rho = 2.0 * Z * rmesh(j,1) / float(n)
         exprho = exp(- rho / 2.0)
         do k = 0, n - l - 1
            const2=(-1)**k/fact(n-l-1-k)/fact(2*l+1+k)/fact(k)
            rlrho = rlrho + const2 * rho ** k
         enddo
         tmp = chi(j)
         chi(j) = const1 * rlrho * rho ** (l + 1) * exprho
      enddo
      en = - (Z / float(n)) ** 2
      jstop = j
      return
      end
      
      function acoth(x)
      if (abs(x).le.1.0) stop 'ACOTH defined for |ARG| > 1 only'
      acoth=log((1.0+x)/(x-1.0))/2.0
      end

      function atanh(x)
      if (abs(x).ge.1.0) stop 'ATANH defined for |ARG| < 1 only'
      atanh=log((1.0+x)/(1.0-x))/2.0
      end
      
      subroutine regular
      return
      end
      subroutine coulomb
      return
      end
      
C  Define L = 0 plane waves on a grid 
      subroutine makechil(lg,gk,npk,minchil,chil,ichi,phasel,nchtop,
     >   sigma)
      include 'par.f'
      common/meshrr/ meshr,rmesh(maxr,3)
      dimension gk(kmax,nchan), phasel(kmax,nchan),npk(nchtop+1),
     >   psi(maxr)
      dimension chil(meshr,ichi),minchil(ichi)
      complex phasel,sigma(nchan)

      do nch = 1, nchtop
         call getchinfo (nch, nt, lg, psi, maxpsi, ea, la, na, l)
         sigma(nch) = (1.0,0.0)
         do k = 1, npk(nch+1) - npk(nch)
            phasel(k,nch) = (1.0,0.0)
            kp = npk(nch) + k - 1
            if (gk(k,nch).gt.0.0) then
               do j = 1, meshr
                  chil(j,kp) = sin(gk(k,nch)*rmesh(j,1)) * rmesh(j,3)
               end do
               minchil(kp) = 1
            else
               minchil(kp) = max(meshr+1,maxr)
            end if 
         end do
      end do
      return
      end

C  This routine returns the CPU time in seconds since last call to CLOCK
C  IBM version
c$$$      subroutine clock( ctime )
c$$$      ctime = mclock() / 100.0
c$$$      end

      subroutine update(n)
c$$$      call flush(n)
      return
      end


C  This routine returns the CPU time in seconds since last call to CLOCK
C  SUN version
      subroutine clock( ctime )
      ctime = 0.0
c$$$      real tarray(2)
c$$$      ctime=etime(tarray)
      return
      end

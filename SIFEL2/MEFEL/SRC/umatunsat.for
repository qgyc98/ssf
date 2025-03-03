! Copyright (C)  2009  C. Tamagnini, E. Sellari, D. Masin, P.A. von Wolffersdorff
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
!  USA.

c------------------------------------------------------------------------------
      subroutine umatunsat(stress,statev,ddsdde,sse,spd,scd,
     &  rpl,ddsddt,drplde,drpldt,
     &  stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     &  ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     &  celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc,unsatvar)
c------------------------------------------------------------------------------
c user subroutine for Abaqus 6.3
c------------------------------------------------------------------------------
c
c Implemented constitutive law:
c -----------------------------
c MASIN HYPO - Masin hypoplastic model with intergranular strains
c D. Masin (2005) A hypoplastic constitutive model for clays. IJNAMG, 29:311-336
c
c Implementation based on:
c Fellin, W. and Ostermann, A. (2002):
c Consistent tangent operators for constitutive rate equations.
c International Journal for Numerical and Analytical Methods in Geomechanics
c
c ----------------------------------------------------------------------------
c The string for the material name may contain 9 characters.
c ----------------------------------------------------------------------------
c Material constants:
c       
c     ---------------------------------------------------------------------
c     props(j)      
c     ---------------------------------------------------------------------
c        1      phi 
c        2      p_t 
c        3      lam_star
c        4      kap_star
c        5      N_star
c        6      rr
c        7      n
c        8      l
c        9      m
c        10    s_e0 
c        11    e_0
c
c     ----------------------------------------------------------------------
c
c Solution dependent state variables (statev):
c definition via sdvini
c
c       1 ... void    void ratio
c       2 ... s    suction
c       3 ... S_r       degree of saturation
c       4 ... nfev    number of function evaluation
c       5 ... phi_mob phi_mob in degrees
c       6 ... dtsub   suggested substep size
c
c Authors: 
c          D. Masin (masin@natur.cuni.cz)
c     C. Tamagnini (tamag@unipg.it)
c
c----------------------------------------------------------------------------
c
      implicit none
c
      character*80 cmname
c
      integer ntens, ndi, nshr, nstatv, nprops, noel, npt,
     & layer, kspt, kstep, kinc
c
      double precision stress(ntens), statev(nstatv),
     &  ddsdde(ntens,ntens), ddsddt(ntens), drplde(ntens),
     &  stran(ntens), dstran(ntens), time(2), predef(1), dpred(1),
     &  props(nprops), coords(3), drot(3,3), dfgrd0(3,3), dfgrd1(3,3),
     &  unsatvar(5)
      double precision sse, spd, scd, rpl, drpldt, dtime, temp, 
     &  dtemp, pnewdt, celent, youngel, nuel
c
c ... 1. nasvdim    = maximum number of additional state variables
c     2. tolintT    = prescribed error tolerance for the adaptive 
c                     substepping scheme
c     3. maxnint    = maximum number of time substeps allowed.
c                     If the limit is exceeded abaqus is forced to reduce 
c                     the overall time step size (cut-back) 
c     4. DTmin      = minimum substeps size allowed.
c                     If the limit is exceeded abaqus is forced to reduce 
c                     the overall time step size (cut-back)
c     5. perturb    = perturbation parameter for numerical computation of Jacobian matrices
c     6. nfasv      = number of first additional state variable in statev field 
c     7. prsw       = switch for printing information
c
c ... declaration of local variables
c
      logical prsw,elprsw
c     
      integer i,error,maxnint,nfev,testnan,maxninttest,inittension
      integer nparms,nasvdim,nfasv,nydim,nasv,nyact,testing
c
      double precision dot_vect_hu
c       
      double precision parms(nprops),theta,tolintT,dtsub,DTmin,perturb
      double precision sig_n(6),sig_np1(6),DDtan(6,6),pore,dsuction
      double precision deps_np1(6),depsv_np1,norm_deps,tolintTtest
      double precision norm_deps2,pp,qq,cos3t,I1,I2,I3,norm_D2,norm_D
c
      parameter (nasvdim = 6)
      parameter (nydim = 6+nasvdim)
c       parameter (tolintT = 1.0d-3) ...orig value...
      parameter (tolintT = 1.0d-3) 
      parameter (tolintTtest = 1.0d-1) 
c
c       parameter (maxnint = 1000) ...orig value...
      parameter (maxnint = 10000)
      parameter (maxninttest = 1000)
      parameter (DTmin = 1.0d-17)
      parameter (perturb = 1.0d-5)
      parameter (nfasv = 1)
      parameter (prsw = .true.)

c
c ... additional state variables
c
      double precision  asv(nasvdim)
c
c ... solution vector (stresses, additional state variables)
c
      double precision  y(nydim),y_n(nydim),dy(nydim)
c
      common /z_nct_errcode/error
c
c ... Error Management:
c     ----------------
c     error =  0 ... no problem in time integration
c     error =  1 ... problems in evaluation of the time rate, (e.g. undefined 
c                    stress state), reduce time integration substeps
c     error =  3 ... problems in time integration, reduce abaqus load increment 
c                    (cut-back)
c     error = 10 ... severe error, terminate calculation
c
      error=0
c
c ... check problem dimensions
c
                
      if (ndi.ne.3) then
c
         write(1,*) 'ERROR: this UMAT can be used only for elm.'
         write(1,*) 'with 3 direct stress/strain components'
         write(1,*) 'noel = ',noel
         error=10
c     
      end if
c
c ... check material parameters and move them to array parms(nparms)
c
      call check_parms_hu(props,nprops,parms,nparms)
c
c ... print informations about time integration, useful when problems occur
c
      elprsw = .false.
      if (prsw) then
c
c ... print only in some defined elements
c
         if ((noel.eq.101).and.(npt.eq.1)) elprsw = .false.
      end if
c
c ... define number of additional state variables
c
      call define_hu(nasv)
      nyact = 6 + nasv
      if (nyact.gt.nydim) then
         write(1,*) 'ERROR: nasvdim too small in UMAT'
         error=10
      end if
c
c ... suggested time substep size, and initial excess pore pressure
c
      dtsub = statev(nasv+4)
      pore = 0

c
c ... vector of additional state variables
c
      do i=1,nasv
        asv(i) = statev(i-1+nfasv)
      end do
c
c ... compute volume strain increment and current net stress tensor
c
       do i=1,6        
          sig_n(i)=0
          deps_np1(i)=0
       end do
       call move_sig_hu(stress,ntens,pore,sig_n)
       call move_eps_hu(dstran,ntens,deps_np1,depsv_np1)
       dsuction=unsatvar(2)
          
       norm_D2=dot_vect_hu(2,deps_np1,deps_np1,6)
       norm_D=sqrt(norm_D2)

c ... check whether the strain rate from the ABAQUS is not NAN          

       testnan=0
       call umatisnan_hu(norm_D,testnan)
       if (testnan .eq. 1) then 
          call wrista_hu(3,y,nydim,deps_np1,
     &         dtime,coords,statev,nstatv,
     &         parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
          write(1,*) 'Error in integration, noel ',noel
          write(1,*) 'Try to decrease the global step size'
          call xit_hu
       end if

c
c --------------------
c ... Time integration
c --------------------
c
       call iniy_hu(y,nydim,nasv,ntens,sig_n,asv)
       call push_hu(y,y_n,nydim)
c
c ... check whether the initial state is not tensile
       inittension=0
       call check_RKF_hu(inittension,y,nyact,nasv,parms,nparms)
       
      if (elprsw) then
        write(1,*) '==================================================='
        write(1,*) 'Call of umat:'
        write(1,*) '==================================================='
        call wrista_hu(3,y,nydim,deps_np1,dtime,coords,statev,nstatv,
     &        parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
      end if
c
c ... parameter of the numerical differentiation: sqrt(macheps)*||deps||
c     double precision
c
      norm_deps2=dot_vect_hu(2,deps_np1,deps_np1,ntens)
      norm_deps=dsqrt(norm_deps2)
      theta=-perturb*max(norm_deps,1.0d-6)  ! negative sign for compression
c
c ... local integration using adaptive RKF-23 method, consistent Jacobian and error estimation
c
      if((dtsub.le.0.0d0).or.(dtsub.gt.dtime)) then
         dtsub = dtime
      end if
c
c     if testing==1 PLAXIS is testing for the initial strain increment.
      testing=0
c     For use in ABAQUS, comment out the following line
      if(kstep.eq.1 .AND. kinc.eq.1) testing=1
      if(norm_D.eq.0) testing=2
c     FEM asking for ddsdde only

      nfev = 0 ! initialisation

      if(inittension.eq.0) then
      if(testing.eq.1) then
         call rkf23_update_hu(y,nyact,nasv,dtsub,tolintTtest,
     &                        maxninttest,DTmin,
     &                        deps_np1,dsuction,parms,nparms,nfev,
     &                        elprsw,dtime,error)
c ... give original state if the model fails without substepping
         if(error.eq.3) then
            do i=1,nyact        
               y(i)=y_n(i)
            end do
            if(unsatvar(2).ne.0) y(6+2)=y_n(6+2)+unsatvar(2)
            error=0
         end if
      else if(testing.eq.2) then
         do i=1,nyact        
            y(i)=y_n(i)
         end do      
         if(unsatvar(2).ne.0) y(6+2)=y_n(6+2)+unsatvar(2)
c     ... Normal RKF23 integration
      else                      !inittension.eq.0 .and. testing.eq.0
         call rkf23_update_hu(y,nyact,nasv,dtsub,tolintT,
     &                        maxnint,DTmin,
     &                        deps_np1,dsuction,parms,nparms,nfev,
     &                        elprsw,dtime,error)
      end if

c
c ... error conditions (if any)
c
      if (error.eq.3) then
c
c ... reduce abaqus load increment
c     
         pnewdt= 0.25d0
c     
         write(1,*) 'subroutine UMAT: reduce step size in ABAQUS'
         call wrista_hu(1,y,nydim,deps_np1,dtime,
     &                  coords,statev,nstatv,
     &                  parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
c               call xit_hu
c               return
     
c ... do not do anything, we are the most likely close to the tensile region
         do i=1,nyact        
            y(i)=y_n(i)
         end do
c
      elseif (error.eq.10) then
c
         call wrista_hu(2,y,nydim,deps_np1,dtime,
     &                  coords,statev,nstatv,
     &                  parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
         call xit_hu
c
      end if
c
c ... update dtsub and nfev
c
      if(dtsub.le.0.0d0) then 
         dtsub = 0
      else if(dtsub.ge.dtime) then 
         dtsub = dtime
      end if
      statev(nasv+4)=dtsub
      statev(nasv+2)=dfloat(nfev)
c
c ... compute consistent tangent via numerical perturbation
c
        call perturbate_hu(y_n,y,nyact,nasv,dtsub,tolintT,maxnint,DTmin,
     &     deps_np1,parms,nparms,nfev,elprsw,theta,ntens,DDtan,dtime)

       else ! inittension.ne.0
c     we were initilly in the tensile stress, calc elastic
          youngel=100
          nuel=0.48
          call calc_elasti_hu(y,nyact,nasv,dtsub,tolintT,
     &                        maxnint,DTmin,
     &                        deps_np1,parms,nparms,nfev,elprsw,
     &                        dtime,DDtan,
     &                        youngel,nuel,error)
          if(unsatvar(2).ne.0) y(6+2)=y_n(6+2)+unsatvar(2)
       endif                    ! end inittension.eq.0

c
c ... convert solution (stress + cons. tangent) to abaqus format
c     update pore pressure and compute total stresses 
c
      call solout_hu(stress,ntens,asv,nasv,ddsdde,
     +            y,nydim,pore,depsv_np1,parms,nparms,DDtan)
     
c
c ... updated vector of additional state variables to abaqus statev vector
c
      do i=1,nasv
        statev(i-1+nfasv) = asv(i) 
      end do

      if(inittension.eq.0) then
      call calc_statev_hu(stress,statev,parms,nparms,nasv,
     &         nasvdim,deps_np1)
      end if
c
c -----------------------
c End of time integration
c -----------------------
c

c Temporary WRC for testing purposes

      unsatvar(1)=unsatvar(1)+unsatvar(2)
      unsatvar(3)=(10/unsatvar(1))**0.38
      unsatvar(4)=unsatvar(3)-(10/(unsatvar(1)-unsatvar(2)))**0.38
      unsatvar(5)=unsatvar(4)/unsatvar(2)

      return
      end
c------------------------------------------------------------------------------
c-----------------------------------------------------------------------------
      subroutine check_parms_hu(props,nprops,parms,nparms)
c-----------------------------------------------------------------------------
c checks input material parameters 
c
c written 10/2004 (Tamagnini & Sellari)
c-----------------------------------------------------------------------------
      implicit none
c
      integer nprops,nparms,i,error
c
      double precision props(nprops),parms(nprops)
      double precision zero,one,four,pi,pi_deg
      double precision phi_deg,phi,lam_star,kap_star,N_star,r_lc,p_ref
      double precision m_R,m_T,r_uc,beta_r,chi,bulk_w,p_t
      double precision lparam,mparam,nparam,separam
c
      parameter(zero=0.0d0,one=1.0d0,four=4.0d0,pi_deg=180.0d0)
c
      common /z_nct_errcode/error
c
      nparms=nprops
c     
      do i=1,nprops
         parms(i)=props(i)
      end do
c
c ... recover material parameters
c
      phi_deg=parms(1)
      lam_star=parms(3)
      kap_star=parms(4)
      N_star=parms(5)
      r_lc=parms(6)
      nparam=parms(7)
      lparam=parms(8)
      mparam=parms(9)
      p_ref=1.d0
      separam=parms(10) 
      p_t=parms(2)
        
c
      pi=four*datan(one)
      phi=phi_deg*pi/pi_deg
      parms(1)=phi
c
      if(phi.le.zero) then
c       
         write(1,*) 'ERROR: subroutine check_parms_hu:'
         write(1,*) 'phi = ',phi
         error = 10
         return 
c
      end if
c
      if(lam_star.le.zero) then
c       
         write(1,*) 'ERROR: subroutine check_parms_hu:'
         write(1,*) 'lam_star = ',lam_star
         error = 10 
         return 
c
      end if
c
      if(kap_star.le.zero) then
c       
         write(1,*) 'ERROR: subroutine check_parms_hu:'
         write(1,*) 'kap_star = ',kap_star
         error = 10 
         return 
c     
      end if
c
      if(N_star.le.zero) then
c     
         write(1,*) 'ERROR: subroutine check_parms_hu:'
         write(1,*) 'N_star = ',N_star
         error = 10 
         return 
c     
      end if
c
      if(r_lc.le.zero) then
c       
         write(1,*) 'ERROR: subroutine check_parms_hu:'
         write(1,*) 'r_lc = ',r_lc
         error = 10 
         return 
c     
      end if
c     
      if(p_ref.lt.zero) then
c     
         write(1,*) 'ERROR: subroutine check_parms_hu:'
         write(1,*) 'p_ref = ',p_ref
         error = 10 
         return 
c     
      end if
c
      if(separam.lt.zero) then
c       
         write(1,*) 'ERROR: subroutine check_parms_hu:'
         write(1,*) 'separam = ',separam
         error = 10 
         return 
c
      end if
c 
      if(p_t.lt.zero) then
c     
         write(1,*) 'ERROR: subroutine check_parms_hu:'
         write(1,*) 'p_t = ',p_t
         error = 10 
         return 
c     
      end if
c     
      return
      end
c-----------------------------------------------------------------------------
      subroutine define_hu(nasv)
c-----------------------------------------------------------------------------
      implicit none 
      integer nasv
c
c number of additional state variables 
c must be less than  18 (otherwise change nasvdim in umat)
c
c    nasv(1) ... void    void ratio
c    nasv(2) ... suction suction
c
c modified 6/2005 (Tamagnini, Sellari & Miriano)
c
      nasv = 2
      return
      end
c------------------------------------------------------------------------------
      double precision function dot_vect_hu(flag,a,b,n)
c------------------------------------------------------------------------------
c dot product of a 2nd order tensor, stored in Voigt notation
c created 10/2004 (Tamagnini & Sellari)
c
c flag = 1 -> vectors are stresses in Voigt notation
c flag = 2 -> vectors are strains in Voigt notation
c flag = 3 -> ordinary dot product between R^n vectors
c------------------------------------------------------------------------------
      implicit none
      integer i,n,flag
      double precision a(n),b(n)
      double precision zero,half,one,two,coeff
c
      parameter(zero=0.0d0,half=0.5d0,one=1.0d0,two=2.0d0)
c
      if(flag.eq.1) then
c
c ... stress tensor (or the like)
c
         coeff=two
c
      elseif(flag.eq.2) then
c
c ... strain tensor (or the like)
c     
         coeff=half
c
      else
c
c ... standard vectors
c
         coeff=one
c       
      end if
c     
      dot_vect_hu=zero
c
      do i=1,n
         if(i.le.3) then
            dot_vect_hu = dot_vect_hu+a(i)*b(i)
         else
            dot_vect_hu = dot_vect_hu+coeff*a(i)*b(i)
         end if
      end do
c     
      return
      end
c-----------------------------------------------------------------------------
      subroutine get_F_sig_q_hu(sig,q,nasv,parms,nparms,deps,
     &                          F_sig,F_q,dsuction,error)
c-----------------------------------------------------------------------------
c
c  finds vectors F_sigma and F_q in F(y)
c
c  written 6/2005 (Tamagnini, Sellari & Miriano)
c-----------------------------------------------------------------------------
      implicit none
      double precision dot_vect_hu
      
c 
      integer nparms,nasv,ii
c
      double precision sig(6),q(nasv),parms(nparms),deps(6)
      double precision MM(6,6),HH(nasv,6),F_sig(6),F_q(nasv)
      double precision LL(6,6),NN(6),norm_D,norm_D2,norm2
      double precision trD,depsh,depd,Aparam,edev(6),dsuction
      double precision Hunsat(6),separam,suction,zero,gamma
      double precision khalratefact,F_sig_new(6)
      integer j,ntens,error
      
      separam=parms(10)
      suction=q(2)
      zero=0
      ntens=6
      
c
c ... compute tangent operators
c
      call get_tan_hu(deps,sig,q,nasv,parms,nparms,MM,
     .                HH,LL,NN,Hunsat,error)
c
c ... compute F_sig=MM*deps
c
      call matmul_hu(LL,deps,F_sig,6,6,1)
      norm_D2=dot_vect_hu(2,deps,deps,6)
      norm_D=dsqrt(norm_D2)
c     
      do ii=1,6
         F_sig(ii)=F_sig(ii)+NN(ii)*norm_D
      end do
      
        khalratefact=-dsuction
        if (suction .gt. separam) then
           gamma=0.55D0
           khalratefact=-(1-gamma)*dsuction*
     &          (separam/suction)**gamma
        end if

        call move_sig_hu(F_sig,ntens,-khalratefact,F_sig_new)
        call move_sig_hu(F_sig_new,ntens,zero,F_sig)

c ...                   wetting collapse contribution                      
        if (suction .gt. separam) then
           if (dsuction .lt. zero) then
              do ii=1,6
                 F_sig(ii)=F_sig(ii)+Hunsat(ii)*dsuction
              end do
           end if
        end if
 
c ... compute F_q=HH*deps
c
        call matmul_hu(HH,deps,F_q,nasv,6,1)
c       F_q(1)=HH(1,1)*(deps(1)+deps(2)+deps(3))
        F_q(2)=dsuction

        return
        end
c-----------------------------------------------------------------------------
      subroutine get_tan_hu(deps,sig,q,nasv,parms,nparms,MM,HH,
     .                      LL,NN,Hunsat,error)
c-----------------------------------------------------------------------------
c  computes matrices M and H for Masin and Khalili (2008) model for unsaturated soils
c-----------------------------------------------------------------------------
      implicit none
c 
      integer nparms,nasv,i,j,error
c
      double precision dot_vect_hu
c     
      double precision sig(6),q(nasv),parms(nparms),deps(6)
      double precision eta(6),eta_dev(6),del(6),void,sig_star(6)
      double precision sensit,H_s(6),fdsbs
      double precision eta_del(6),eta_delta(6),eta_eps(6)
      double precision norm_del,norm_del2,norm_deps,norm_deps2,eta_dn2
      double precision pp,qq,cos3t,I1,I2,I3,tanpsi
      double precision a,a2,FF,alpha,fd,fdi,fs,c_1,c_2,Yi,YY
      double precision num,den,aF,Fa2,eta_n2,norm_m,norm_m2
      double precision II(6,6),IU(6,6),Hunsat(6),dpeds_divpe
      double precision MM(6,6),HH(nasv,6),LL(6,6),NN(6),AA(6,6),m(6)
      double precision NNpure(6),AAinv(6,6),AinvN(6),determ
      integer softmodel,ntens
      double precision m_dir(6),m_dir1(6),Leta(6),H_del(6,6),H_e(6)
      double precision load,rho,N_par,Stf,lparam,nparam,mparam,separam
      double precision kap_par,r_par,logsse,fufact,pehvor
      double precision zero,tiny,half,one,two,three,six,eight,nine
      double precision onethird,sqrt3,twosqrt2,sqrt2,oneeight,ln2m1
      double precision temp1,temp2,temp3,temp4,suction,lam_par
      double precision phi,lam_star,kap_star,N_star,r_lc,r_uc,p_ref
      double precision m_R,m_T,beta_r,chi,bulk_w,p_t,sinphi,sinphi2
      double precision sig_star_ef(6),ikron(6)
      real inv
c     
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0,six=6.0d0)
      parameter(tiny=1.0d-17,half=0.5d0,eight=8.0d0,nine=9.0d0)
      parameter(ntens=6)
c     
c ... initialize constants and vectors
c
      onethird=one/three
      sqrt3=dsqrt(three)
      twosqrt2=two*dsqrt(two)
      sqrt2=dsqrt(two)
      oneeight=one/eight
      onethird=one/three
      ln2m1=one/dlog(two)
c     
      do i=1,6
         do j=1,6
            MM(i,j)=zero
            LL(i,j)=zero
            II(i,j)=zero
            IU(i,j)=zero
            H_del(i,j)=zero
         end do
         eta_del(i)=zero
         eta_delta(i)=zero
         eta_eps(i)=zero
         ikron(i)=zero
      end do
c     
      do j=1,6
         HH(1,j)=zero
      end do
c     
c     ... fourth order identity tensors in Voigt notation
c     
      II(1,1)=one
      II(2,2)=one
      II(3,3)=one
      II(4,4)=half
      II(5,5)=half
      II(6,6)=half
c     
      IU(1,1)=one
      IU(2,2)=one
      IU(3,3)=one
      IU(4,4)=one
      IU(5,5)=one
      IU(6,6)=one
c     
      ikron(1)=one
      ikron(2)=one
      ikron(3)=one
c     
c     ... recover material parameters
c     
      phi=parms(1)
      lam_par=parms(3)
      kap_par=parms(4)
      N_par=parms(5)
      r_par=parms(6)
      nparam=parms(7)
      lparam=parms(8)
      mparam=parms(9)
      p_ref=1.d0
      separam=parms(10) 
      p_t=parms(2)
c     
      sinphi=dsin(phi)
      sinphi2=sinphi*sinphi
      
c     
c     ... recover internal state variables
c     
      void=q(1)
      suction=q(2)  
c     
c     ... axis translation due to cohesion (p_t>0)
c     
      sig_star(1)=sig(1)-p_t
      sig_star(2)=sig(2)-p_t
      sig_star(3)=sig(3)-p_t
      sig_star(4)=sig(4)
      sig_star(5)=sig(5)
      sig_star(6)=sig(6)
      
c     ... calculate Khalili stress
      call calc_khalili_stress(sig_star,sig_star_ef,
     &                         separam,suction,ntens)
      call move_sig_hu(sig_star_ef,ntens,zero,sig_star)
        
c     calculate N_star and lambda_star using current suction
      kap_star=kap_par
      r_lc=r_par
      
      if(suction>separam) then
         logsse=dlog(suction/separam)
      else
         logsse=0
      end if
      
      N_star=N_par+nparam*logsse
      lam_star=lam_par+lparam*logsse
c     
c     ... auxiliary stress tensors
c     
      call inv_sig_hu(sig_star,pp,qq,cos3t,I1,I2,I3)
c     
      
      eta(1)=sig_star(1)/I1
      eta(2)=sig_star(2)/I1
      eta(3)=sig_star(3)/I1
      eta(4)=sig_star(4)/I1
      eta(5)=sig_star(5)/I1
      eta(6)=sig_star(6)/I1
c     
      eta_dev(1)=eta(1)-onethird
      eta_dev(2)=eta(2)-onethird
      eta_dev(3)=eta(3)-onethird
      eta_dev(4)=eta(4)
      eta_dev(5)=eta(5)
      eta_dev(6)=eta(6)
c     
c     ... functions a and F
c     
      eta_dn2=dot_vect_hu(1,eta_dev,eta_dev,6)
      tanpsi=sqrt3*dsqrt(eta_dn2)
      temp1=oneeight*tanpsi*tanpsi+
     &      (two-tanpsi*tanpsi)/(two+sqrt2*tanpsi*cos3t)
      temp2=tanpsi/twosqrt2
c
      a=sqrt3*(three-sin(phi))/(twosqrt2*sin(phi))
      a2=a*a
      FF=dsqrt(temp1)-temp2
c     
c     ... barotropy and pyknotropy functions
c     
      temp1=(lam_star-kap_star)/(lam_star+kap_star)
      temp2=(three+a2)/(sqrt3*a)
      alpha=ln2m1*dlog(temp1*temp2)
c     
      temp1=two*I1/(three*p_ref)
      pehvor=dexp((N_star-dlog(one+void))/lam_star)
      fd=(-temp1/pehvor)**alpha
c     
      fdi=two**alpha
      temp1=-I1/lam_star
      temp2=three+a2-fdi*sqrt3*a
      fs=temp1/temp2
      
c     
c     ... tensor L
c     
      temp1=two/(nine*r_lc)
      c_1=temp1*temp2
      c_2=one+(one-c_1)*three/a2
c     
      do i = 1,6
         do j=1,6
            LL(i,j)=three*c_1*II(i,j)+
     &              three*c_2*a2*eta(i)*eta(j)
         end do
      end do
c     
c     ... function YY
c     
      
      Yi=sqrt3*a/(three+a2)
      num=(I1*I2+nine*I3)*(one-sinphi2)
      den=eight*I3*sinphi2
      YY=(Yi-one)*(num/den)+Yi
c     
c     ... tensor m and NN
c     
      aF=a/FF
      Fa2=FF*FF/a2
      eta_n2=dot_vect_hu(1,eta,eta,6)
      temp1=onethird*(six*eta_n2-one)/(Fa2+eta_n2)
c     
      do i=1,6
         m(i)=-aF*(eta(i)+eta_dev(i)-temp1*eta(i))
      end do
      norm_m2=dot_vect_hu(1,m,m,6)
      norm_m=sqrt(norm_m2)
c
      m_dir(1)=m(1)/norm_m
      m_dir(2)=m(2)/norm_m
      m_dir(3)=m(3)/norm_m
      m_dir(4)=m(4)/norm_m
      m_dir(5)=m(5)/norm_m
      m_dir(6)=m(6)/norm_m
c     
      m_dir1(1)=-YY*m_dir(1)
      m_dir1(2)=-YY*m_dir(2)
      m_dir1(3)=-YY*m_dir(3)
      m_dir1(4)=-YY*two*m_dir(4)
      m_dir1(5)=-YY*two*m_dir(5)
      m_dir1(6)=-YY*two*m_dir(6)
c     
      call matmul_hu(LL,m_dir1,NNpure,6,6,1)
      
      do i=1,6
         do j=1,6
            LL(i,j)=LL(i,j)*fs
         end do
         NN(i)=NNpure(i)*fs*fd
      end do        
      
      do j=1,6
         Hunsat(j)=0
      end do
      
      if(suction>separam) then
         dpeds_divpe=(nparam-log(pehvor)*lparam)/(lam_star*suction)
         do j=1,6
            Hunsat(j)=sig_star(j)*dpeds_divpe
         end do
c     
         do i=1,6
            do j=1,6
               AA(i,j)=LL(i,j)+sig_star(i)*ikron(j)/lam_star
            end do
         end do   
         
         call inverse(AA, ntens, ntens, AAinv)
         
         call matmul_hu(AAinv,NNpure,AinvN,6,6,1)
         do i=1,6
            AinvN(i)=AinvN(i)*fs
         end do        
         norm_m2=dot_vect_hu(1,AinvN,AinvN,6)
         norm_m=sqrt(norm_m2)
          
         fdSBS=1/norm_m        
         fufact=(fd/fdSBS)**(mparam/alpha)
         
          do j=1,6
             Hunsat(j)=fufact*Hunsat(j)
          end do
       end if

c
c ... void ratio evolution function (tension positive)
c        
        do i=1,6 
           if (i.le.3) then
              H_e(i)=one+void
           else
              H_e(i)=zero
           end if
        end do        
        
        do j=1,6
           HH(1,j)=H_e(j)
        end do
        
        return
        end
c-----------------------------------------------------------------------------
      subroutine iniy_hu(y,nydim,nasv,ntens,sig,qq)
c-----------------------------------------------------------------------------
c initializes the vector of state variables
c-----------------------------------------------------------------------------
      implicit none
c
      integer i,nydim,nasv,ntens
c
      double precision y(nydim),qq(nasv),sig(ntens)
c
      do i=1,nydim
        y(i) = 0
      end do
c
      do i=1,ntens
        y(i) = sig(i)
      end do
c
c additional state variables
c
      do i=1,nasv
        y(6+i) = qq(i)
      end do
c
      return
      end
c------------------------------------------------------------------------------
      subroutine inv_eps_hu(eps,eps_v,eps_s,sin3t)
c------------------------------------------------------------------------------
c calculate invariants of strain tensor
c------------------------------------------------------------------------------
c
      implicit none
c
      integer i
c
      double precision eps(6),edev(6),edev2(6),ev3
      double precision tredev3,eps_v,eps_s,sin3t
      double precision norm2,numer,denom
c
      double precision zero,one,two,three,six
      double precision onethird,twothirds,sqrt6
c
      data zero,one,two,three,six/0.0d0,1.0d0,2.0d0,3.0d0,6.0d0/
c
c ... some constants
c
      onethird=one/three
      twothirds=two/three
      sqrt6=dsqrt(six)
c     
c ... volumetric strain
c
      eps_v=eps(1)+eps(2)+eps(3)
c     
      ev3=onethird*eps_v
c     
c ... deviator strain
c     
      edev(1)=eps(1)-ev3
      edev(2)=eps(2)-ev3
      edev(3)=eps(3)-ev3
      edev(4)=eps(4)/two
      edev(5)=eps(5)/two
      edev(6)=eps(6)/two
c
c ... second invariant
c
      norm2=edev(1)*edev(1)+edev(2)*edev(2)+edev(3)*edev(3)+
     &      two*(edev(4)*edev(4)+edev(5)*edev(5)+edev(6)*edev(6))
c
      eps_s=dsqrt(twothirds*norm2)
c
c ... components of (edev_ij)(edev_jk)
c
      edev2(1)=edev(1)*edev(1)+edev(4)*edev(4)+edev(5)*edev(5)
      edev2(2)=edev(4)*edev(4)+edev(2)*edev(2)+edev(6)*edev(6)
      edev2(3)=edev(6)*edev(6)+edev(5)*edev(5)+edev(3)*edev(3)
      edev2(4)=two*(edev(1)*edev(4)+edev(4)*edev(2)+edev(6)*edev(5))
      edev2(5)=two*(edev(5)*edev(1)+edev(6)*edev(4)+edev(3)*edev(5))
      edev2(6)=two*(edev(4)*edev(5)+edev(2)*edev(6)+edev(6)*edev(3))
c     
c ... Lode angle
c
      if(eps_s.eq.zero) then 
c     
         sin3t=-one
c     
      else
c     
         tredev3=zero
         do i=1,6
            tredev3=tredev3+edev(i)*edev2(i)
         end do
c
         numer=sqrt6*tredev3
         denom=(dsqrt(norm2))**3
         sin3t=numer/denom
         if(dabs(sin3t).gt.one) then
            sin3t=sin3t/dabs(sin3t)
         end if
c     
      end if 
c
      return
      end
c------------------------------------------------------------------------------
      subroutine inv_sig_hu(sig,pp,qq,cos3t,I1,I2,I3)
c------------------------------------------------------------------------------
c calculate invariants of stress tensor
c
c NOTE: Voigt notation is used with the following index conversion
c
c       11 -> 1
c       22 -> 2
c    33 -> 3
c       12 -> 4
c       13 -> 5
c       23 -> 6
c
c------------------------------------------------------------------------------
c
      implicit none
c
      double precision sig(6),sdev(6)
      double precision eta(6),eta_d(6),eta_d2(6)
      double precision xmin1,xmin2,xmin3
      double precision tretadev3,pp,qq,cos3t,I1,I2,I3
      double precision norm2,norm2sig,norm2eta,numer,denom
c     
      double precision half,one,two,three,six
      double precision onethird,threehalves,sqrt6,tiny
c     
      double precision dot_vect_hu
c
      data half,one/0.5d0,1.0d0/
      data two,three,six/2.0d0,3.0d0,6.0d0/
      data tiny/1.0d-18/
c
c ... some constants
c
      onethird=one/three
      threehalves=three/two
      sqrt6=dsqrt(six)
c     
c ... trace and mean stress
c
      I1=sig(1)+sig(2)+sig(3)
      pp=onethird*I1
c
c ... deviator stress
c
      sdev(1)=sig(1)-pp
      sdev(2)=sig(2)-pp
      sdev(3)=sig(3)-pp
      sdev(4)=sig(4)
      sdev(5)=sig(5)
      sdev(6)=sig(6)
c     
c ... normalized stress and dev. normalized stress
c
      eta(1)=sig(1)/I1
      eta(2)=sig(2)/I1
      eta(3)=sig(3)/I1
      eta(4)=sig(4)/I1
      eta(5)=sig(5)/I1
      eta(6)=sig(6)/I1
c
      eta_d(1)=eta(1)-onethird
      eta_d(2)=eta(2)-onethird
      eta_d(3)=eta(3)-onethird
      eta_d(4)=eta(4)
      eta_d(5)=eta(5)
      eta_d(6)=eta(6)
c     
c ... second invariants
c     
      norm2=dot_vect_hu(1,sdev,sdev,6)
      norm2sig=dot_vect_hu(1,sig,sig,6)
      norm2eta=dot_vect_hu(1,eta_d,eta_d,6)
c     
      qq=dsqrt(threehalves*norm2)
      I2=half*(norm2sig-I1*I1)
c
c ... components of (eta_d_ij)(eta_d_jk)
c
      eta_d2(1)=eta_d(1)*eta_d(1)+eta_d(4)*eta_d(4)+eta_d(5)*eta_d(5)
      eta_d2(2)=eta_d(4)*eta_d(4)+eta_d(2)*eta_d(2)+eta_d(6)*eta_d(6)
      eta_d2(3)=eta_d(6)*eta_d(6)+eta_d(5)*eta_d(5)+eta_d(3)*eta_d(3)
      eta_d2(4)=eta_d(1)*eta_d(4)+eta_d(4)*eta_d(2)+eta_d(6)*eta_d(5)
      eta_d2(5)=eta_d(5)*eta_d(1)+eta_d(6)*eta_d(4)+eta_d(3)*eta_d(5)
      eta_d2(6)=eta_d(4)*eta_d(5)+eta_d(2)*eta_d(6)+eta_d(6)*eta_d(3)
c     
c ... Lode angle
c     
      if(norm2eta.lt.tiny) then 
c     
         cos3t=-one
c     
      else
c     
         tretadev3=dot_vect_hu(1,eta_d,eta_d2,6)
c     
         numer=-sqrt6*tretadev3
         denom=(dsqrt(norm2eta))**3
         cos3t=numer/denom
         if(dabs(cos3t).gt.one) then
            cos3t=cos3t/dabs(cos3t)
         end if
c     
      end if 
c
c ... determinant
c     
      xmin1=sig(2)*sig(3)-sig(6)*sig(6)
      xmin2=sig(4)*sig(3)-sig(6)*sig(5)
      xmin3=sig(4)*sig(6)-sig(5)*sig(2)
c     
      I3=sig(1)*xmin1-sig(4)*xmin2+sig(5)*xmin3
c
      return
      end
c------------------------------------------------------------------------------
      subroutine matmul_hu(a,b,c,l,m,n)
c------------------------------------------------------------------------------
c matrix multiplication
c------------------------------------------------------------------------------
      implicit none
c
      integer i,j,k,l,m,n
c
      double precision a(l,m),b(m,n),c(l,n)
c
      do i=1,l
        do j=1,n
          c(i,j) = 0.0d0
          do k=1,m
            c(i,j) = c(i,j) + a(i,k)*b(k,j)
          end do
        end do
      end do
c
      return
      end
c-----------------------------------------------------------------------------
      subroutine move_asv_hu(asv,nasv,qq_n)
c-----------------------------------------------------------------------------
c move internal variables in vector qq_n and changes intergranular strain 
c from continuum to soil mechanics convention
c
c NOTE: del has always 6 components
c
c written 6/2005 (Tamagnini, Sellari & Miriano)
c-----------------------------------------------------------------------------
      implicit none
      integer nasv,i
      double precision asv(nasv),qq_n(nasv),zero 
c
      parameter(zero=0.0d0)
c
      do i=1,nasv
         qq_n(i)=zero
      end do
c
c ... void ratio stored in qq_n(1)
c
      qq_n(1) = asv(1) 
c
c ... suction
c
      qq_n(2) = asv(2) 
c
      return
      end
c-----------------------------------------------------------------------------
      subroutine move_eps_hu(dstran,ntens,deps,depsv)
c-----------------------------------------------------------------------------
c Move strain increment dstran into deps and computes 
c volumetric strain increment
c
c NOTE: all strains negative in compression; deps has always 6 components
c
c written 7/2005 (Tamagnini, Sellari & Miriano)
c-----------------------------------------------------------------------------
      implicit none
      integer ntens,i
      double precision deps(6),dstran(ntens),depsv
c
      do i=1,ntens
         deps(i) = dstran(i)
      end do
c
      depsv=deps(1)+deps(2)+deps(3)
c
      return
      end
c-----------------------------------------------------------------------------
      subroutine move_sig_hu(stress,ntens,pore,sig)
c-----------------------------------------------------------------------------
c computes effective stress from total stress (stress) and pore pressure (pore)
c
c NOTE: stress = total stress tensor (tension positive)
c       pore   = exc. pore pressure (undrained conds., compression positive)
c       sig    = effective stress (tension positive)
c
c       sig has always 6 components
c
c written 7/2005 (Tamagnini, Sellari & Miriano)
c-----------------------------------------------------------------------------
      implicit none
      integer ntens,i
      double precision sig(6),stress(ntens),pore,zero 
c
      parameter(zero=0.0d0)
c
      do i=1,6
         sig(i)=zero
      end do
c
      do i=1,ntens
         if(i.le.3) then
            sig(i) = stress(i)+pore
         else
            sig(i) = stress(i)
         end if
      end do
c     
      return
      end
c-----------------------------------------------------------------------------
      subroutine norm_res_hu(y_til,y_hat,ny,nasv,norm_R)
c-----------------------------------------------------------------------------
c  evaluate norm of residual vector Res=||y_hat-y_til||
c
c  written 6/2005 (Tamagnini, Sellari & Miriano)
c-----------------------------------------------------------------------------
      implicit none
c 
      integer ny,nasv,ng,k,i,testnan
c
      double precision y_til(ny),y_hat(ny),void_til,void_hat,del_void
      double precision sensit_til,sensit_hat,del_sensit
      double precision err(ny),norm_R2,norm_R
      double precision norm_sig2,norm_q2,norm_sig,norm_q
      double precision sig_hat(6),sig_til(6),del_sig(6)
      double precision q_hat(nasv),q_til(nasv),del_q(nasv)
      double precision dot_vect_hu,zero
c
      parameter(zero=0.0d0)
c
      ng=6*nasv
      k=42+nasv
c
      do i=1,ny
         err(i)=zero
      end do
c
c ... recover stress tensor and internal variables
c
      do i=1,6
         sig_hat(i)=y_hat(i)
         sig_til(i)=y_til(i)
         del_sig(i)=dabs(sig_hat(i)-sig_til(i))
      end do
c     
c
      void_hat=y_hat(6+1)
      void_til=y_til(6+1)
      del_void=dabs(void_hat-void_til)
c
c ... relative error norms
c
      norm_sig2=dot_vect_hu(1,sig_hat,sig_hat,6)
      norm_sig=dsqrt(norm_sig2)
c
      if(norm_sig.gt.zero) then
         do i=1,6
            err(i)=del_sig(i)/norm_sig
         end do
      end if
c     
c     
      err(6+nasv-1)=del_void/void_hat
      
c
c ... global relative error norm
c
      norm_R2=dot_vect_hu(3,err,err,ny)
      norm_R=dsqrt(norm_R2)
c     
      testnan=0
      call umatisnan_hu(norm_sig,testnan)
c     T. Koudelka - norm_q is uninitialized
c     call umatisnan_hu(norm_q,testnan)
      call umatisnan_hu(void_hat,testnan)
c     T. Koudelka - sensit_hat is uninitialized
c     call umatisnan_hu(sensit_hat,testnan)
      if(testnan.eq.1) then
         norm_R=1.d20
      end if
      
      return
      end

c-----------------------------------------------------------------------------
      subroutine perturbate_hu(y_n,y_np1,n,nasv,dtsub,err_tol,maxnint,
     &    DTmin,deps_np1,parms,nparms,nfev,elprsw,theta,ntens,DD, dtime)
c-----------------------------------------------------------------------------
c
c  compute numerically consistent tangent stiffness
c
c  written 12/2005 (Tamagnini)
c-----------------------------------------------------------------------------
      implicit none
c 
      logical elprsw
c
      integer ntens,jj,kk,i
      integer n,nasv,nparms,nfev
      integer maxnint,error
c
      double precision y_n(n),y_np1(n),y_star(n),parms(nparms)
      double precision dtsub,err_tol,DTmin, dtime
      double precision theta,sig(6),q(nasv)
      double precision deps_np1(6),deps_star(6)
      double precision dsig(6),DD(6,6),HHtmp(nasv,6)
      double precision LL(6,6),NN(6),Hunsattmp(6)
      double precision zero
c
      parameter(zero=0.0d0)
c
      common /z_nct_errcode/error
c
c ... initialize DD and y_star
c 

        do kk=1,6
           do jj=1,6
              DD(kk,jj)=zero
           end do
        end do
        do i=1,6
           sig(i)=y_n(i)
        end do
        do i=1,nasv
           q(i)=y_n(6+i)
        end do
        
        call push_hu(y_n,y_star,n)
        
        if(error.ne.10) then
           call get_tan_hu(deps_np1,sig,q,nasv,parms,nparms,
     .                     DD,HHtmp,LL,NN,Hunsattmp,error)
        end if
        
        do kk=1,6
           do jj=1,6
              DD(kk,jj)=LL(kk,jj)
           end do
        end do

        return
        end        
        
c-----------------------------------------------------------------------------
      subroutine push_hu(a,b,n)
c-----------------------------------------------------------------------------
c push_hu vector a into vector b
c-----------------------------------------------------------------------------
      implicit none
      integer i,n
      double precision a(n),b(n) 
c
      do i=1,n
         b(i)=a(i)
      end do
c
      return
      end
c-----------------------------------------------------------------------------
      subroutine rhs_hu(y,ny,nasv,parms,nparms,deps,kRK,nfev,
     &                  dsuction,error)
c-----------------------------------------------------------------------------
c calculate coefficient kRK from current state y and strain increment deps
c Masin hypoplastic model for clays with intergranular strains
c
c written 12/2005 (Tamagnini & Sellari)
c-----------------------------------------------------------------------------
      implicit none
c
      integer error,ny,nparms,nasv,i,nfev
c
      double precision zero,one,two,four 
      double precision y(ny),kRK(ny),parms(nparms),deps(6)
      double precision sig(6),q(nasv),dsuction
      double precision F_sig(6),F_q(nasv)
c
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,four=4.0d0)
c
c ... update counter for the number of function f(y) evaluations
c
      nfev=nfev+1
c
c ... initialize kRK
c
      do i=1,ny
         kRK(i)=zero
      end do
c
c ... recover current state variables (sig,q)                   
c
      do i=1,6
         sig(i)=y(i)
      end do
c
      do i=1,nasv
         q(i)=y(6+i)
      end do
c     
c ... build F_sig(6) and F_q(nasv) vectors and move them into kRK
c
      call get_F_sig_q_hu(sig,q,nasv,parms,nparms,deps,
     &                    F_sig,F_q,dsuction,error)
      if(error.eq.10) return
c
      do i=1,6
c
         kRK(i)=F_sig(i)
c     
      end do                   
c       
      do i=1,nasv
c
         kRK(6+i)=F_q(i)
c
      end do                   
c
      return
      end

c-----------------------------------------------------------------------------
      subroutine rkf23_update_hu(y,n,nasv,dtsub,err_tol,maxnint,DTmin,
     &                           deps_np1,dsuction,parms,nparms,
     &                           nfev,elprsw,dtime,error)

c-----------------------------------------------------------------------------
c
c  numerical solution of y'=f(y)
c  explicit, adapive RKF23 scheme with local time step extrapolation
c
c  Tamagnini, Sellari & Miriano 6/2005
c
c-----------------------------------------------------------------------------
      implicit none
c
      logical elprsw
c
      integer n,nasv,nparms,i,ksubst,kreject,nfev
      integer maxnint,error,error_RKF
c
      double precision y(n),parms(nparms),dtsub,err_tol,DTmin
      double precision zero,half,one,two,three,four,six
      double precision ptnine,onesixth,onethird,twothirds,temp
c     
      double precision deps_np1(6),y_k(n),y_2(n),y_3(n),y_til(n)
      double precision y_hat(n)
      double precision T_k,DT_k,dtime
      double precision kRK_1(n),kRK_2(n),kRK_3(n)
      double precision norm_R,S_hull,dsuction
c
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0)
      parameter(four=4.0d0,six=6.0d0,half=0.5d0,ptnine=0.9d0)
c
c ... initialize y_k vector and other variables
c
        do i=1,n
           y_k(i)=zero
        end do
c
        onesixth=one/six
        onethird=one/three
        twothirds=two/three
c
c ... start of update process
c
                
        error_RKF=0
        T_k=zero      
        DT_k=dtsub/dtime
        ksubst=0
        kreject=0
        nfev=0
c
        do i=1,n
           y_k(i)=y(i)
        end do
c
c ... start substepping 
c
        do while(T_k.lt.one) 
c
           ksubst=ksubst+1
c
c ... write substepping info
c
c     write(*,1234) ksubst,T_k,DT_k
c     1234           format('Substep no.',i4,' -- T_k = ',d12.4,' -- DT_k = ',d12.4)
c
c ... check for maximum number of substeps
c
           if(ksubst.gt.maxnint) then
              write(1,*) 'number of substeps ',ksubst,
     &             ' is too big, step rejected'
              error=3
              return
           end if          
c
c ... build RK functions
c
           call check_RKF_hu(error_RKF,y_k,n,nasv,parms,nparms)
           if(error_RKF.eq.1) then 
              error=3
              return
           else
              call rhs_hu(y_k,n,nasv,parms,nparms,
     .                    deps_np1,kRK_1,nfev,dsuction,error)
           end if
           if(error.eq.10) return
c
c ... find y_2
c
           temp=half*DT_k
c     
           do i=1,n
              y_2(i)=y_k(i)+temp*kRK_1(i)
           end do

c               
           call check_RKF_hu(error_RKF,y_2,n,nasv,parms,nparms)
           if(error_RKF.eq.1) then 
              error=3
              return
           else
              call rhs_hu(y_2,n,nasv,parms,nparms,
     .                    deps_np1,kRK_2,nfev,dsuction,error)
           end if
           if(error.eq.10) return
c     
c ... find y_3
c

           do i=1,n
              y_3(i)=y_k(i)-DT_k*kRK_1(i)+two*DT_k*kRK_2(i)
           end do
c

           call check_RKF_hu(error_RKF,y_3,n,nasv,parms,nparms)
           if(error_RKF.eq.1) then 
              error=3
              return
           else
              call rhs_hu(y_3,n,nasv,parms,nparms,
     .                    deps_np1,kRK_3,nfev,dsuction,error)
           end if
           if(error.eq.10) return

c                               
c ... approx. solutions of 2nd (y_til) and 3rd (y_hat) order
c
           do i=1,n        
              y_til(i)=y_k(i)+DT_k*kRK_2(i)
              y_hat(i)=y_k(i)+DT_k*
     &        (onesixth*kRK_1(i)+twothirds*kRK_2(i)+onesixth*kRK_3(i))
           end do
c
c ... local error estimate
c

           call norm_res_hu(y_til,y_hat,n,nasv,norm_R)
c     check if output y_hat can be used as an input into the next step
           call check_RKF_hu(error_RKF,y_hat,n,nasv,parms,nparms)
           
           if (error_RKF.ne.0) then
c     error=1.d20
c     error_RKF=0
              error=3
              return
           end if
c
c ... time step size estimator according to Hull
c               
           if(norm_R .ne. 0) then
              S_hull=ptnine*DT_k*(err_tol/norm_R)**onethird
           else
              S_hull=1
           end if
c

           if (norm_R.lt.err_tol) then                             
c
c ... substep is accepted, update y_k and T_k and estimate new substep size DT_k
c
              do i=1,n        
                 y_k(i)=y_hat(i)
              end do
c     
              T_k=T_k+DT_k
              DT_k=min(four*DT_k,S_hull)
              dtsub=DT_k*dtime
              DT_k=min((one-T_k),DT_k)        
c     
           else
c
c ... substep is not accepted, recompute with new (smaller) substep size DT
c
              DT_k=max(DT_k/four,S_hull)
c
c ... check for minimum step size
c
              if(DT_k.lt.DTmin) then
                 write(1,*) 'substep size ',DT_k,
     &                ' is too small, step rejected'
                 error=3
                 return
              end if          
c                                       
           end if                                                  
c
c ... bottom of while loop
c
        end do
        
c
c ... recover final state
c
      do i=1,n
         y(i)=y_k(i)
      end do
c
      return
      end
c

c-----------------------------------------------------------------------------
      subroutine check_RKF_hu(error_RKF,y,ny,nasv,parms,nparms)
c-----------------------------------------------------------------------------
c Checks is RKF23 solout vector y is OK for hypoplasticity
c-----------------------------------------------------------------------------
      implicit none
c
      integer error_RKF,ny,nasv,i,nparms,testnan,iopt,ntens
c
      double precision y(ny),parms(nparms)
      double precision sig(6),pmean,sig_star(6),sig_star_ef(6)
      double precision xN1(3),xN2(3),xN3(3),S(3),P,Q,tmin
      double precision p_t,minstress,suction,zero,separam
c
      ntens=6
      zero=0.
c
      p_t    =parms(2)
      minstress=p_t/4.d0
      do i=1,6
         sig(i)=y(i)
      end do

      sig_star(1)=sig(1)-p_t
      sig_star(2)=sig(2)-p_t
      sig_star(3)=sig(3)-p_t
      sig_star(4)=sig(4)
      sig_star(5)=sig(5)
      sig_star(6)=sig(6)
        
      suction=y(8)
      separam=parms(10) 
      call calc_khalili_stress(sig_star,sig_star_ef,
     &                         separam,suction,ntens)
      call move_sig_hu(sig_star_ef,ntens,zero,sig_star)
                
      pmean=-(sig_star(1)+sig_star(2)+sig_star(3))/3
            
c     check for positive mean stress
      if(pmean .le. minstress) then
         error_RKF=1
      end if
c
c     calculate minimum principal stress
c     
      iopt=0
      Call PrnSig_hu(iopt, sig_star, xN1, xN2, xN3,
     &               S(1),S(2),S(3), P, Q)
      tmin     = 1.0d+20
      do i=1,3
         if(tmin .ge. -S(i)) then
            tmin=-S(i)
         endif         
      enddo 
c     
c     check for tension
c     
      if(tmin .le. minstress) then
         error_RKF=1
      end if
        
c     check for NAN
      testnan=0
      do i=1,ny
         call umatisnan_hu(y(i),testnan)
      end do
      if(testnan.eq.1) error_RKF=1
      
c
c     
      return
      end
c

c-----------------------------------------------------------------------------
      subroutine solout_hu(stress,ntens,asv,nasv,ddsdde,y,nydim,
     +                     pore,depsv_np1,parms,nparms,DD)
c-----------------------------------------------------------------------------
c copy the vector of state variables to umat output
c modified 7/2005 (Tamagnini, Sellari)
c
c NOTE: solid mechanics convention for stress and strain components
c       pore is always positive in compression
c-----------------------------------------------------------------------------
      implicit none
c
      integer nydim,nasv,nparms,ntens,i,j
c
      double precision y(nydim),asv(nasv),stress(ntens)
      double precision ddsdde(ntens,ntens),DD(6,6)
      double precision parms(nparms),bulk_w,pore,depsv_np1 

c
c updated total stresses (effective stresses stored in y(1:6))
c
      do i=1,ntens
         stress(i) = y(i)
      end do
c
c additional state variables (void ratio, suction)
c
      do i=1,nasv
         asv(i) = y(6+i)
      end do
c
c consistent tangent stiffness
c
      do j=1,ntens
        do i=1,ntens
          ddsdde(i,j) = DD(i,j)      
        end do
      end do

      return
      end

c-----------------------------------------------------------------------------
      subroutine wrista_hu(mode,y,nydim,deps_np1,dtime,coords,statev,
     &                 nstatv,parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
c-----------------------------------------------------------------------------
c ... subroutine for managing output messages
c
c     mode
c
c     all = writes:             kstep, kinc, noel, npt
c       2   = writes also:      error message,coords(3),parms(nparms),ndi,nshr,stress(nstress)
c                                               deps(nstress),dtime,statev(nstatv)
c     3   = writes also:        stress(nstress),deps(nstress),dtime,statev(nstatv)
c-----------------------------------------------------------------------------
      implicit none
c
      integer mode,nydim,nstatv,nparms,noel,npt,ndi,nshr,kstep,kinc,i
c
      double precision y(nydim),statev(nstatv),parms(nparms)
      double precision deps_np1(6),coords(3),dtime
c
c ... writes for mode = 2
c
      if (mode.eq.2) then
        write(1,*) '==================================================='
        write(1,*) 'ERROR: abaqus job failed during call of UMAT'
        write(1,*) '==================================================='
        write(1,*) 'state dump:'
        write(1,*) 
      end if
c
c ... writes for all mode values
c
      write(1,111) 'Step: ',kstep, 'increment: ',kinc,
     &     'element: ', noel, 'Integration point: ',npt
      write(1,*) 
c
c ... writes for mode = 2
c
      if (mode.eq.2) then
         write(1,*) 'Co-ordinates of material point:'
         write(1,104) 'x1 = ',coords(1),' x2 = ',coords(2),' x3 = ',
     &        coords(3)
         write(1,*) 
         write(1,*) 'Material parameters:'
         write(1,*) 
         do i=1,nparms
            write(1,105) 'prop(',i,') = ',parms(i)
         end do 
         write(1,*)
         write(1,102) 'No. of mean components:  ',ndi
         write(1,102) 'No. of shear components: ',nshr
         write(1,*)
      end if
c     
c ... writes for mode = 2 or 3
c     
      if ((mode.eq.2).or.(mode.eq.3)) then
        write(1,*) 'Stresses:'
        write(1,*) 
        write(1,101) 'sigma(1) = ',y(1)
        write(1,101) 'sigma(2) = ',y(2)
        write(1,101) 'sigma(3) = ',y(3)
        write(1,101) 'sigma(4) = ',y(4)
        write(1,101) 'sigma(5) = ',y(5)
        write(1,101) 'sigma(6) = ',y(6)
        write(1,*) 
        write(1,*) 'Strain increment:'
        write(1,*) 
        write(1,101) 'deps_np1(1) = ',deps_np1(1)
        write(1,101) 'deps_np1(2) = ',deps_np1(2)
        write(1,101) 'deps_np1(3) = ',deps_np1(3)
        write(1,101) 'deps_np1(4) = ',deps_np1(4)
        write(1,101) 'deps_np1(5) = ',deps_np1(5)
        write(1,101) 'deps_np1(6) = ',deps_np1(6)
        write(1,*) 
        write(1,*) 'Time increment:'
        write(1,*) 
        write(1,108) 'dtime = ',dtime
        write(1,*) 
        write(1,*) 'Internal variables:'
        write(1,*) 
        write(1,109) 'void = ',statev(1)
        write(1,*) 
        write(1,*) '==================================================='
c     
      end if
c     
101   format(1X,a15,e11.4)
102   format(1X,a25,i1)
103   format(1X,a7,i5)
104   format(1X,3(a5,f10.4,2X))
105   format(1X,a5,i2,a4,f20.3)
106   format(1X,3(a9,f12.4,2X))
107   format(1X,3(a10,f12.4,2X))
108   format(1X,a8,f12.4)
109   format(1X,a6,f10.4)
110   format(1X,a5,f10.4)
111   format(1X,a6,i4,2X,a11,i4,2X,a9,i10,2X,a19,i4)
c       
      return
      end

      
c-----------------------------------------------------------------------------
      subroutine calc_statev_hu(stress,statev,parms,nparms,nasv,
     & nasvdim,deps)
c-----------------------------------------------------------------------------
c
c  computes additional state variables for postprocessing
c
c-----------------------------------------------------------------------------
      implicit none
c 
      logical elprsw
c
      integer ntens,jj,kk,i
      integer n,nasv,nparms,nfev,nasvdim
      integer maxnint,error
c
      double precision parms(nparms),dot_vect_hu
      double precision stress(6),statev(nasvdim)
      double precision deps(6),tmax,tmin
      double precision MM(6,6),HHtmp(nasv,6)
      double precision LL(6,6),NN(6),sig_star_ef(6)
      double precision zero,two,four,iopt,three
      double precision I1,I2,I3,cos3t,pp,qq
      double precision sin2phi,sinphi,sig_star(6),p_t
      double precision norm_del,norm_del2,del(6)
c     
      parameter(zero=0.0d0,two=2.0d0,four=4.0d0,three=3.0d0,ntens=6)
c

c ... calc phimob (statev 5) from Matsuoka-Nakai YS

      p_t    =parms(2)
      do i=1,3
         sig_star(i)=stress(i)-p_t
      end do
      do i=4,6
         sig_star(i)=stress(i)
      end do
      call calc_khalili_stress(sig_star,sig_star_ef,
     &                         parms(10),statev(2),ntens)
      call move_sig_hu(sig_star_ef,ntens,zero,sig_star)

      call inv_sig_hu(sig_star,pp,qq,cos3t,I1,I2,I3)
      if(I3 .ne. 0) then
         sin2phi=(9.d0+I1*I2/I3)/(1.d0+I1*I2/I3)
      else 
         sin2phi=0
      end if
      if(sin2phi .lt. 0) then
         sin2phi=0
      end if 
      if(sin2phi .gt. 1) then
         sin2phi=1
      end if 
      sinphi=sqrt(sin2phi)
      
      statev(nasv+3)= asin(sinphi)*
     .         180.0d0/3.141592d0
      
      return
      end         
        
c-----------------------------------------------------------------------------
      subroutine umatisnan_hu(chcknum,testnan)
c-----------------------------------------------------------------------------
c
c  checks whether number is NaN
c
c-----------------------------------------------------------------------------
      double precision chcknum
      integer testnan
        
      if (.not.(chcknum .ge. 0. .OR. chcknum .lt. 0.)) testnan=1        
      if (chcknum .gt. 1.d30) testnan=1        
      if (chcknum .lt. -1.d30) testnan=1        
      if (chcknum .ne. chcknum) testnan=1        
      
      return
      end         
      
c-----------------------------------------------------------------------------
      subroutine xit_hu
c-----------------------------------------------------------------------------
      stop
c
      return
      end

C***********************************************************************
      Subroutine PrnSig_hu(IOpt,S,xN1,xN2,xN3,S1,S2,S3,P,Q)
      Implicit Double Precision (A-H,O-Z)
      Dimension S(*),xN1(*),xN2(*),xN3(*)

      If (iOpt.Eq.1) Then
        Call Eig_3_hu(0,S,xN1,xN2,xN3,S1,S2,S3,P,Q) ! with Eigenvectors
      Else
        Call Eig_3a_hu(0,S,S1,S2,S3,P,Q) ! no Eigenvectors
      End If
      Return
      End
C***********************************************************************
      Subroutine Eig_3_hu(iOpt,St,xN1,xN2,xN3,S1,S2,S3,P,Q)
      Implicit Double Precision (A-H,O-Z)
      Dimension St(6),A(3,3),V(3,3),
     *     xN1(3),xN2(3),xN3(3)
      !
      ! Get Eigenvalues/Eigenvectors for 3*3 matrix
      ! Wim Bomhof 15/11/'01
      ! PGB : adaption to Principal stress calculation
      !
      ! Applied on principal stresses, directions
      ! Stress vector St(): XX, YY, ZZ, XY, YZ, ZX
      !
      A(1,1) = St(1) ! xx
      A(1,2) = St(4) ! xy = yx
      A(1,3) = St(6) ! zx = xz

      A(2,1) = St(4) ! xy = yx
      A(2,2) = St(2) ! yy
      A(2,3) = St(5) ! zy = yz

      A(3,1) = St(6) ! zx = xz
      A(3,2) = St(5) ! zy = yz
      A(3,3) = St(3) ! zz

      ! Set V to unity matrix
      V(1,1) = 1
      V(2,1) = 0
      V(3,1) = 0

      V(1,2) = 0
      V(2,2) = 1
      V(3,2) = 0

      V(1,3) = 0
      V(2,3) = 0
      V(3,3) = 1


      abs_max_s=0.0
      Do i=1,3
         Do j=1,3
            if (abs(a(i,j)) .Gt. abs_max_s) abs_max_s=abs(a(i,j))
         End Do
      End Do
      Tol = 1d-20 * abs_max_s
      it = 0
      itmax = 50
      Do While ( it.Lt.itMax .And.
     *     abs(a(1,2))+abs(a(2,3))+abs(a(1,3)) .Gt. Tol )
         it=it+1
         Do k=1,3
            If (k .Eq. 1) Then
               ip=1
               iq=2
            Else If (k .Eq.2) Then
               ip=2
               iq=3
            Else
               ip=1
               iq=3
            End If
            If (abs(a(ip,iq)) .gt. Tol) Then
               tau=(a(iq,iq)-a(ip,ip))/(2.0*a(ip,iq))
               If (tau .Ge.0.0) Then
                  sign_tau=1.0
               Else
                  sign_tau=-1.0
               End If
               t=sign_tau/(abs(tau)+sqrt(1.0+tau*tau))
               c=1.0/sqrt(1.0+t*t)
               s=t*c
               a1p=c*a(1,ip)-s*a(1,iq)
               a2p=c*a(2,ip)-s*a(2,iq)
               a3p=c*a(3,ip)-s*a(3,iq)
               a(1,iq)=s*a(1,ip)+c*a(1,iq)
               a(2,iq)=s*a(2,ip)+c*a(2,iq)
               a(3,iq)=s*a(3,ip)+c*a(3,iq)
               a(1,ip)=a1p
               a(2,ip)=a2p
               a(3,ip)=a3p
               
               v1p=c*v(1,ip)-s*v(1,iq)
               v2p=c*v(2,ip)-s*v(2,iq)
               v3p=c*v(3,ip)-s*v(3,iq)
               v(1,iq)=s*v(1,ip)+c*v(1,iq)
               v(2,iq)=s*v(2,ip)+c*v(2,iq)
               v(3,iq)=s*v(3,ip)+c*v(3,iq)
               v(1,ip)=v1p
               v(2,ip)=v2p
               v(3,ip)=v3p
               
               ap1=c*a(ip,1)-s*a(iq,1)
               ap2=c*a(ip,2)-s*a(iq,2)
               ap3=c*a(ip,3)-s*a(iq,3)
               a(iq,1)=s*a(ip,1)+c*a(iq,1)
               a(iq,2)=s*a(ip,2)+c*a(iq,2)
               a(iq,3)=s*a(ip,3)+c*a(iq,3)
               a(ip,1)=ap1
               a(ip,2)=ap2
               a(ip,3)=ap3
            End If              ! a(ip,iq)<>0
         End Do                 ! k
      End Do                    ! While
      ! principal values on diagonal of a
      S1 = a(1,1)
      S2 = a(2,2)
      S3 = a(3,3)
      ! Derived invariants
      P = (S1+S2+S3)/3
      Q = Sqrt( ( (S1-S2)**2 + (S2-S3)**2 + (S3-S1)**2 ) / 2 )

      ! Sort eigenvalues S1 <= S2 <= S3
      is1 = 1
      is2 = 2
      is3 = 3
      if (s1.Gt.s2) Then
         t   = s2
         s2  = s1
         s1  = t
         it  = is2
         is2 = is1
         is1 = it
      End If
      if (s2.Gt.s3) Then
        t   = s3
        s3  = s2
        s2  = t
        it  = is3
        is3 = is2
        is2 = it
      End If
      if (s1.Gt.s2) Then
        t   = s2
        s2  = s1
        s1  = t
        it  = is2
        is2 = is1
        is1 = it
      End If
      Do i=1,3
         xN1(i) = v(i,is1)      ! first  column
         xN2(i) = v(i,is2)      ! second column
         xN3(i) = v(i,is3)      ! third  column
      End Do
      Return
      End                       ! Eig_3

      Subroutine Eig_3a_hu(iOpt,St,S1,S2,S3,P,Q) ! xN1,xN2,xN3,
      Implicit Double Precision (A-H,O-Z)
      Dimension St(6),A(3,3)    !  V(3,3),xN1(3),xN2(3),xN3(3)
      !
      ! Get Eigenvalues ( no Eigenvectors) for 3*3 matrix
      ! Wim Bomhof 15/11/'01
      !
      ! Applied on principal stresses, directions
      ! Stress vector XX, YY, ZZ, XY, YZ, ZX
      !
      A(1,1) = St(1)            ! xx
      A(1,2) = St(4)            ! xy = yx
      A(1,3) = St(6)            ! zx = xz
      
      A(2,1) = St(4)            ! xy = yx
      A(2,2) = St(2)            ! yy
      A(2,3) = St(5)            ! zy = yz
      
      A(3,1) = St(6)            ! zx = xz
      A(3,2) = St(5)            ! zy = yz
      A(3,3) = St(3)            ! zz
      
      abs_max_s=0.0
      Do i=1,3
         Do j=1,3
            if (abs(a(i,j)) .Gt. abs_max_s) abs_max_s=abs(a(i,j))
         End Do
      End Do
      Tol = 1d-20 * abs_max_s
      If (iOpt.Eq.1) Tol = 1d-50*abs_max_s
      it=0
      itmax = 50
      
      Do While ( it.lt.itmax .And.
     *          abs(a(1,2))+abs(a(2,3))+abs(a(1,3)) .Gt. Tol )
         
         it=it+1
         Do k=1,3
            If (k .Eq. 1) Then
               ip=1
               iq=2
            Else If (k .Eq.2) Then
               ip=2
               iq=3
            Else
               ip=1
               iq=3
            End If
            
            If (abs(a(ip,iq)) .gt. Tol) Then ! ongelijk nul ?
               tau=(a(iq,iq)-a(ip,ip))/(2.0*a(ip,iq))
               If (tau .Ge.0.0) Then
                  sign_tau=1.0
               Else
                  sign_tau=-1.0
               End If
               t=sign_tau/(abs(tau)+sqrt(1.0+tau*tau))
               c=1.0/sqrt(1.0+t*t)
               s=t*c
               a1p=c*a(1,ip)-s*a(1,iq)
               a2p=c*a(2,ip)-s*a(2,iq)
               a3p=c*a(3,ip)-s*a(3,iq)
               a(1,iq)=s*a(1,ip)+c*a(1,iq)
               a(2,iq)=s*a(2,ip)+c*a(2,iq)
               a(3,iq)=s*a(3,ip)+c*a(3,iq)
               a(1,ip)=a1p
               a(2,ip)=a2p
               a(3,ip)=a3p
               
               ap1=c*a(ip,1)-s*a(iq,1)
               ap2=c*a(ip,2)-s*a(iq,2)
               ap3=c*a(ip,3)-s*a(iq,3)
               a(iq,1)=s*a(ip,1)+c*a(iq,1)
               a(iq,2)=s*a(ip,2)+c*a(iq,2)
               a(iq,3)=s*a(ip,3)+c*a(iq,3)
               a(ip,1)=ap1
               a(ip,2)=ap2
               a(ip,3)=ap3
            End If              ! a(ip,iq)<>0
         End Do                 ! k
      End Do                    ! While
                                ! principal values on diagonal of a
      S1 = a(1,1)
      S2 = a(2,2)
      S3 = a(3,3)
                                ! Derived invariants
      P = (S1+S2+S3)/3
      Q = Sqrt( ( (S1-S2)**2 + (S2-S3)**2 + (S3-S1)**2 ) / 2 )
      
      if (s1.Gt.s2) Then
         t   = s2
         s2  = s1
         s1  = t
      End If
      if (s2.Gt.s3) Then
         t   = s3
         s3  = s2
         s2  = t
      End If
      if (s1.Gt.s2) Then
         t   = s2
         s2  = s1
         s1  = t
      End If
      Return
      End                       ! Eig_3a
      
c-----------------------------------------------------------------------------
      subroutine calc_khalili_stress(sig_net,sig_ef,
     &                               separam,suction,ntens)
c-----------------------------------------------------------------------------
c initializes the vector of state variables
c-----------------------------------------------------------------------------
      implicit none
c
      integer i,ntens
c
      double precision sig_net(6),sig_ef(6),separam,suction
      double precision khalfact
      double precision gamma
c
      gamma=0.55D0
      if(suction>separam) then
         khalfact=-suction*(separam/suction)**gamma
      else
         khalfact=-suction
      end if
      call move_sig_hu(sig_net,ntens,khalfact,sig_ef)
      
      return
      end

!****************************************************************************************
!********from http://www.algarcia.org/nummeth/Programs2E.html******************************
!****************************************************************************************/

      subroutine inverse(AAA, N, Nm, Ainv)
       
      integer N, Nm
      double precision AAA(Nm,Nm), Ainv(Nm,Nm)
! Compute inverse of matrix
! Input
!    AA   -    Matrix A (N by N)
!    N    -    Dimension of matrix A (used)
!    Nm   -    Dimension of matrix A (allocated memory)
! Outputs
!   Ainv  -    Inverse of matrix A (N by N)
!  determ -    Determinant of matrix A (return value)

      integer MAXN
      parameter( MAXN = 6 )
      integer i, j, k, index(MAXN), signDet, jPivot, indexJ
      double precision scale(MAXN), b(MAXN,MAXN)   ! Scale factor and work array
      double precision A(MAXN,MAXN) ! Working copy of input matrix
      double precision scalemax, ratio, ratiomax, coeff, determ, sum

      if( Nm .gt. MAXN ) then
         write(1,*) 'ERROR - Matrix is too large for inv routine'
         stop
      endif
      
      ! Copy matrix A so as not to modify original
      do i=1,N
         do j=1,N
            A(i,j) = AAA(i,j)
         enddo
      enddo
      
      !* Matrix b is initialized to the identity matrix
      do i=1,N
         do j=1,N
            if( i .eq. j ) then
               b(i,j) = 1.0
            else
               b(i,j) = 0.0
            endif
         enddo
      enddo
      
      !* Set scale factor, scale(i) = max( |a(i,j)| ), for each row
      do i=1,N
         index(i) = i           ! Initialize row index list
         scalemax = 0.0
         do j=1,N
            if( abs(A(i,j)) .gt. scalemax ) then
               scalemax = abs(A(i,j))
            endif
         enddo
         scale(i) = scalemax
      enddo
      
      !* Loop over rows k = 1, ..., (N-1)
      signDet = 1
      do k=1,N-1
         !* Select pivot row from max( |a(j,k)/s(j)| )
         ratiomax = 0.0
         jPivot = k
         do i=k,N
            ratio = abs(A(index(i),k))/scale(index(i))
            if( ratio .gt. ratiomax ) then
               jPivot=i
               ratiomax = ratio
            endif
         enddo
         !* Perform pivoting using row index list
         indexJ = index(k)
         if( jPivot .ne. k ) then ! Pivot
            indexJ = index(jPivot)
            index(jPivot) = index(k) ! Swap index jPivot and k
            index(k) = indexJ
            signDet = -1*signDet ! Flip sign of determinant
         endif
         !* Perform forward elimination
         do i=k+1,N
            coeff = A(index(i),k)/A(indexJ,k)
            do j=k+1,N
               A(index(i),j) = A(index(i),j) - coeff*A(indexJ,j)
            enddo
            A(index(i),k) = coeff
            do j=1,N
               b(index(i),j) = b(index(i),j) - A(index(i),k)*b(indexJ,j)
            enddo
         enddo
      enddo


      !* Compute determinant as product of diagonal elements
      determ = signDet          ! Sign of determinant
      do i=1,N
         determ = determ*A(index(i),i)
      enddo

      !* Perform backsubstitution
      do k=1,N
         Ainv(N,k) = b(index(N),k)/A(index(N),N)
         do i=N-1,1,-1
            sum = b(index(i),k)
            do j=i+1,N
               sum = sum - A(index(i),j)*Ainv(j,k)
            enddo
            Ainv(i,k) = sum/A(index(i),i)
         enddo
      enddo
c     inv = determ   ! Return the determinant
      return
      end
c-----------------------------------------------------------------------------
      subroutine calc_elasti_hu(y,n,nasv,dtsub,err_tol,maxnint,DTmin,
     &                          deps_np1,parms,nparms,nfev,elprsw,
     &                          dtime,DDtan,youngel,nuel,error)
c-----------------------------------------------------------------------------
c
c  numerical solution of y'=f(y)
c  explicit, adapive RKF23 scheme with local time step extrapolation
c
c  Tamagnini, Sellari & Miriano 6/2005
c
c-----------------------------------------------------------------------------
      implicit none
c
      logical elprsw
c
      integer n,nasv,nparms,i,ksubst,kreject,nfev
      integer maxnint,error,error_RKF,tension,j
c
      double precision y(n),parms(nparms),dtsub,err_tol,DTmin
      double precision zero,half,one,two,three,four,six
      double precision ptnine,onesixth,onethird,twothirds,temp
c     
      double precision deps_np1(6),y_k(n),y_2(n),y_3(n),y_til(n)
      double precision y_hat(n),DDtan(6,6)
      double precision T_k,DT_k,dtime,II(6,6),krondelta(6)
      double precision kRK_1(n),kRK_2(n),kRK_3(n)
      double precision norm_R,S_hull,youngel,nuel,F_sig(6)
c
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0)
      parameter(four=4.0d0,six=6.0d0,half=0.5d0,ptnine=0.9d0)
c
c ... initialize y_k vector and other variables
c
      do i=1,n
         y_k(i)=zero
      end do
c     
      onesixth=one/six
      onethird=one/three
      twothirds=two/three
      
c
c ... fourth order identity tensors in Voigt notation
c
      do i = 1,6
         do j=1,6
            II(i,j)=zero
         end do
      end do
      
      II(1,1)=one
      II(2,2)=one
      II(3,3)=one
      II(4,4)=half
      II(5,5)=half
      II(6,6)=half
c     
      krondelta(1)=one
      krondelta(2)=one
      krondelta(3)=one
      krondelta(4)=zero
      krondelta(5)=zero
      krondelta(6)=zero
c
c ... Elastic stiffness tensor 
c
      do i = 1,6
         do j=1,6
            DDtan(i,j)=(youngel/(1+nuel))*(II(i,j) + 
     &                 nuel/(1-2*nuel)*krondelta(i)*krondelta(j));
         end do
      end do
        
      call matmul_hu(DDtan,deps_np1,F_sig,6,6,1)
      do i=1,6
         y(i)=y(i)+F_sig(i)
      end do
      
      return
      end
c     

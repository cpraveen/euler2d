C Freestream values
      double precision mach_inf, aoa, aoa_deg, q_inf, u_inf, v_inf, 
     &                 r_inf, p_inf, T_inf, T_infd, ent_inf, 
     &                 prim_inf(nvar), a_inf, H_inf
      common/inf/mach_inf,aoa,aoa_deg,q_inf,u_inf,v_inf,r_inf,p_inf,
     &           T_inf,T_infd,ent_inf,prim_inf,a_inf,H_inf

C Vortex correction for farfield BC
      integer          vortex
      double precision xref, yref
      common/farfieldbc/xref, yref, vortex

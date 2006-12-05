C laminar paramters
C Rey = Reynolds number
C SCONST = Constant in Sutherland Law
      integer          iflow
      double precision Rey, prandtl, prandtl_turb, SCONST
      common/viscparam/Rey, prandtl, prandtl_turb, SCONST, iflow

      double precision Cb1, Cb2, sigma_sa, kolm, Cw1, Cw2, Cw3, Cv1, 
     &                 Cv2, Cv11, Cw31, Cw32, kolm2, Cb2Sig1, Cb2Sig2
      common/samodel/Cb1, Cb2, sigma_sa, kolm, Cw1, Cw2, Cw3, Cv1,
     &               Cv2, Cv11, Cw31, Cw32, kolm2, Cb2Sig1, Cb2Sig2

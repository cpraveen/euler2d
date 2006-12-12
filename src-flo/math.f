C------------------------------------------------------------------------------
C.....Error function, from Abromovitz and Stegun
C------------------------------------------------------------------------------
      double precision function ERRF(X)
      double precision X,ARG,E,VB,T,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

      ARG = X*X
      if(ARG .lt. 20.0d0)then
            E = dexp(-ARG)
      else
            E = 0.0d0
      endif
      VB = dabs(X)
      T = 1.0d0/(1.0d0 + 0.3275911d0*VB)
      tmp1 = 1.061405429d0*T
      tmp2 = (tmp1 - 1.453152027d0)*T
      tmp3 = (tmp2 + 1.421413741d0)*T
      tmp4 = (tmp3 - 0.284496736d0)*T
      tmp5 = (tmp4 + 0.254829592d0)*T
      tmp6 = 1.0d0 - tmp5*E
      if(X .lt. 0.0d0)then
            ERRF = -tmp6
      else
            ERRF =  tmp6
      endif
      return
      end

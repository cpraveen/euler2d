      integer no, yes
      parameter(no=0, yes=1)

c nvar = number of variables, fixed at 5, which includes turbulent
c viscosity
      integer nvar
      parameter(nvar=4)

C     GAMMA = ratio of specific heats
C     GAMMA1= GAMMA - 1
C     GAS_CONST = gas constant, this can be set to 1.0
C     M_PI = value of pi
      double precision GAMMA, GAMMA1, GAS_CONST, M_PI
      common/const/GAMMA, GAMMA1, GAS_CONST, M_PI


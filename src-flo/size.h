c npmax = number of points
c ntmax = number of triangles (elements in general)
c nemax = number of edges
      integer npmax, ntmax, nemax, nspmax, nfpmax, nopmax, nbpmax, 
     &        nbemax
      parameter(npmax = 24000,
     &          ntmax = 2*npmax,
     &          nemax = 3*npmax,
     &          nspmax= 2000,
     &          nfpmax= 2000,
     &          nopmax= 2000,
     &          nbpmax= 2000,
     &          nbemax= 2000)

C Actual number of points, elements, edges and boundary edges
C     np = number of points
C     nt = number of triangles/elements
C     ne = number of edges
C     nsp= number of solid boundary points
C     nfp= number of farfield boundary points
C     nop= number of outflow boundary points
C     nbp= number of boundary points, must equal nsp+nfp+nop
      integer np,nt,ne,nsp,nfp,nop,nbp,nbe,nip
      common/dims/np,nt,ne,nsp,nfp,nop,nbp,nbe,nip

      integer nsw, nin, nff, nif, nof,
     +        nsw1, nin1, nff1, nif1, nof1,
     +        nsw2, nin2, nff2, nif2, nof2
      common/nedge/nsw, nin, nff, nif, nof,
     +             nsw1, nin1, nff1, nif1, nof1,
     +             nsw2, nin2, nff2, nif2, nof2

c Maximum elements surrounding a point
      integer mesup
      parameter(mesup=10)

c Maximum elements surrounding a point
      integer mesubp
      parameter(mesubp=7)

c Maximum points surrounding a point
      integer mpsup
      parameter(mpsup=10)

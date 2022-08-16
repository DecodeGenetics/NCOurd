import numpy as np

class Grid:
  def __init__(self,type='linear',ext=False):
    self.type = type
    self.ext  = ext
    self.k    = self.init_kgrid()
    self.i    = self.init_igrid(self.k)
    
    if type == 'linear':
      self.w = self.init_linear()
    if type == 'cubic':
      self.w = self.init_cubic()

  def init_kgrid(self):
    # Set up the grid for tract lengths
    kgrid = np.arange(1,1000,1)
    kgrid = np.append(kgrid,np.arange(1000,10000,10))
    if self.ext:
      kgrid = np.append(kgrid,np.arange(10000,100000,100))
      kgrid = np.append(kgrid,np.arange(100000,300001,1000))
    else:
      kgrid = np.append(kgrid,np.arange(10000,100001,100))

    return kgrid

  def init_igrid(self,kgrid):
    # Index vector for extracting grid points from full range
    return kgrid - np.full(kgrid.shape,1)

  def init_linear(self):

    # Set up the coefficients for summation on the kgrid
    # These coefficients are accurate for up to second 
    # degree polynomials

    # First the basic coefficients
    wgrid = np.full(999,1.)
    wgrid = np.append(wgrid,np.full(900,10.))
    if self.ext:
      wgrid = np.append(wgrid,np.full(900,100.))
      wgrid = np.append(wgrid,np.full(201,1000.))
    else:
      wgrid = np.append(wgrid,np.full(901,100.))

    #----- Grid size 10 ----- 
    U = 9/2

    # Upper boundary coefficients at 1000
    idx = 999 
    wgrid[idx+0] = 1+U


    # Lower boundary coefficients at 10kb
    idx = 1899
    wgrid[idx-0] = 1+U

    #----- Grid size 100 -----
    U = 99/2

    # Upper boundary coefficients at 10kb
    idx = 1899
    wgrid[idx+0] += U

    # Lower boundary coefficients at 100kb
    idx = 2799
    wgrid[idx-0] = 1+U

    if self.ext:
      #----- Grid size 1000 -----
      U = 999/2

      # Upper boundary coefficients at 100kb
      idx = 2799
      wgrid[idx+0] += U

      # Lower boundary coefficients at 300kb
      idx = 2999
      wgrid[idx-0] = 1+U

    return wgrid

  def init_cubic(self):

    # Set up the coefficients for summation on the kgrid
    # These coefficients are accurate for up to second 
    # degree polynomials

    # First the basic coefficients
    wgrid = np.full(999,1.)
    wgrid = np.append(wgrid,np.full(900,10.))
    if self.ext:
      wgrid = np.append(wgrid,np.full(900,100.))
      wgrid = np.append(wgrid,np.full(201,1000.))
    else:
      wgrid = np.append(wgrid,np.full(901,100.))
    
    #----- Grid size 10 ----- 
    A = 33/80
    B = 393/80
    U = 261/80
    V = 591/80

    # Upper boundary coefficients at 1000
    idx = 999 
    wgrid[idx+0] = 1+U-A
    wgrid[idx+1] = 1+V+B-A
    wgrid[idx+2] = 1+2*B-6*A
    wgrid[idx+3] = 1+2*B-A

    # Lower boundary coefficients at 10kb
    idx = 1899
    wgrid[idx-3] = 1+2*B-A
    wgrid[idx-2] = 1+2*B-6*A
    wgrid[idx-1] = 1+V+B-A
    wgrid[idx-0] = 1+U-A

    #----- Grid size 100 -----
    A = 3333/800
    B = 42933/800
    U = 29601/800
    V = 62931/800

    # Upper boundary coefficients at 10kb
    idx = 1899
    wgrid[idx+0] += U-A
    wgrid[idx+1]  = 1+V+B-A
    wgrid[idx+2]  = 1+2*B-6*A
    wgrid[idx+3]  = 1+2*B-A

    # Lower boundary coefficients at 100kb
    idx = 2799
    wgrid[idx-3] = 1+2*B-A
    wgrid[idx-2] = 1+2*B-6*A
    wgrid[idx-1] = 1+V+B-A
    wgrid[idx-0] = 1+U-A

    if self.ext:
      #----- Grid size 1000 -----
      A = 333333/8000
      B = 4329333/8000
      U = 2996001/8000
      V = 6329331/8000

      # Upper boundary coefficients at 100kb
      idx = 2799
      wgrid[idx+0] += U-A
      wgrid[idx+1]  = 1+V+B-A
      wgrid[idx+2]  = 1+2*B-6*A
      wgrid[idx+3]  = 1+2*B-A

      # Lower boundary coefficients at 300kb
      idx = 2999
      wgrid[idx-3] = 1+2*B-A
      wgrid[idx-2] = 1+2*B-6*A
      wgrid[idx-1] = 1+V+B-A
      wgrid[idx-0] = 1+U-A

    return wgrid


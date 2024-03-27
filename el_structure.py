import xml.etree.ElementTree as ET
from numpy import exp
import numpy as np
from matplotlib import pyplot as plt
Ha_to_ev=13.605662285137*2
hartree2ry=2
def w0gauss(x,degauss=0.06):
  x2=x/degauss
  sqrtpm1= 1. / 1.77245385090551602729
  sqrt2=2.**0.5
  # cold smearing  (Marzari-Vanderbilt-DeVita-Payne)
#  arg = min (200./13.606, (x2 - 1.0 /  sqrt2 ) **2.)
#  return sqrtpm1 *exp ( - arg) * (2.0 - sqrt2 * x2)/degauss 
  #Methfessel-Paxton
  try: arg=min(200/13.606,x2**2)
  except: arg=np.array([min(200/13.606,m**2) for m in x2])
  return exp(-arg)*sqrtpm1/degauss 


class el_structure():
 def __init__(self,structure):
  self.ENE_fs=[]
  self.ENE=[]
  self.ef=0.
  self.bands_num=[]
  self.fermi_nbnd_el=0
  self.minband=0
  self.maxband=0
  self.prefix=structure.prefix
  self.tmp_dir=structure.tmp_dir
  self.soc=0

 def construct_tetra(self,structure):
  [nk1,nk2,nk3]= structure.no_of_kpoints
  self.tetra=[[0 for k in range(6*nk1*nk2*nk3)] for i in range(4)]
  for i in range(nk1):
   for j in range(nk2):
    for k in range(nk3):
           #  n1-n8 are the indices of k-point 1-8 forming a cube
      ip1 = (i+1)%nk1 #np.mod(i+1,nk1)
      jp1 = (j+1)%nk2 #np.mod(j+1,nk2)
      kp1 = (k+1)%nk3 #np.mod(k+1,nk3)
      n1 = k   + j  *nk3 + i   *nk2*nk3 
      n2 = k   + j  *nk3 + ip1 *nk2*nk3 
      n3 = k   + jp1*nk3 + i   *nk2*nk3 
      n4 = k   + jp1*nk3 + ip1 *nk2*nk3 
      n5 = kp1 + j  *nk3 + i   *nk2*nk3 
      n6 = kp1 + j  *nk3 + ip1 *nk2*nk3 
      n7 = kp1 + jp1*nk3 + i   *nk2*nk3 
      n8 = kp1 + jp1*nk3 + ip1 *nk2*nk3 
           #  there are 6 tetrahedra per cube (and nk1*nk2*nk3 cubes)
      n  = 6 * ( k + j*nk3 + i*nk3*nk2 )
      self.tetra [0][n] =n1
      self.tetra [1][n] = n2
      self.tetra [2][n] = n3
      self.tetra [3][n] = n6

      self.tetra [0][n+1] = n2
      self.tetra [1][n+1] = n3
      self.tetra [2][n+1] = n4
      self.tetra [3][n+1] = n6

      self.tetra [0][n+2] = n1
      self.tetra [1][n+2] = n3
      self.tetra [2][n+2] = n5
      self.tetra [3][n+2] = n6

      self.tetra [0][n+3] = n3
      self.tetra [1][n+3] = n4
      self.tetra [2][n+3] = n6
      self.tetra [3][n+3] = n8

      self.tetra [0][n+4] = n3
      self.tetra [1][n+4] = n6
      self.tetra [2][n+4] = n7
      self.tetra [3][n+4] = n8

      self.tetra [0][n+5] = n3
      self.tetra [1][n+5] = n5
      self.tetra [2][n+5] = n6
      self.tetra [3][n+5] = n7
 
 def calc_dos(tetra, gamma,et,ibnd, ef=0):
  ntetra=len(tetra[0])
  Tint = 0.0 
  o13 = 1.0 /3.0 
  eps  = 1.0e-14
  voll  = 1.0 /ntetra
  P1 = 0.0 
  P2 = 0.0 
  P3 = 0.0 
  P4 = 0.0 
  for nt in range(ntetra):
      #
      # etetra are the energies at the vertexes of the nt-th tetrahedron
      #
     etetra=[]
     for i in range(4):
        etetra.append( et[ibnd][tetra[i][nt]])

     itetra=np.argsort(etetra) #hpsort (4,etetra,itetra)
     etetra=[ etetra[i] for i in itetra]
      #
      # ...sort in ascending order: e1 < e2 < e3 < e4
      #
     [e1,e2,e3,e4] = [etetra[0],etetra[1],etetra[2],etetra[3]]

      #
      # kp1-kp4 are the irreducible k-points corresponding to e1-e4
      #
     ik1,ik2,ik3,ik4 = tetra[itetra[0]][nt],tetra[itetra[1]][nt],tetra[itetra[2]][nt],tetra[itetra[3]][nt]
     Y1,Y2,Y3,Y4  = gamma[ibnd][ik1],gamma[ibnd][ik2],gamma[ibnd][ik3],gamma[ibnd][ik4]

     if( e3 < ef and ef < e4): # THEN
        f14 = (ef-e4)/(e1-e4)
        f24 = (ef-e4)/(e2-e4)
        f34 = (ef-e4)/(e3-e4)

        G  =  3.0  * f14 * f24 * f34 / (e4-ef)
        P1 =  f14 * o13
        P2 =  f24 * o13
        P3 =  f34 * o13
        P4 =  (3.0  - f14 - f24 - f34 ) * o13

     elif ( e2 < ef and ef < e3 ):

        f13 = (ef-e3)/(e1-e3)
        f31 = 1.0  - f13
        f14 = (ef-e4)/(e1-e4)
        f41 = 1.0 -f14
        f23 = (ef-e3)/(e2-e3)
        f32 = 1.0  - f23
        f24 = (ef-e4)/(e2-e4)
        f42 = 1.0  - f24

        G   =  3.0  * (f23*f31 + f32*f24)
        P1  =  f14 * o13 + f13*f31*f23 / G
        P2  =  f23 * o13 + f24*f24*f32 / G
        P3  =  f32 * o13 + f31*f31*f23 / G
        P4  =  f41 * o13 + f42*f24*f32 / G
        G   =  G / (e4-e1)

     elif ( e1 < ef and ef < e2 ) :

        f12 = (ef-e2)/(e1-e2)
        f21 = 1.0  - f12
        f13 = (ef-e3)/(e1-e3)
        f31 = 1.0  - f13
        f14 = (ef-e4)/(e1-e4)
        f41 = 1.0  - f14

        G  =  3.0  * f21 * f31 * f41 / (ef-e1)
        P1 =  o13 * (f12 + f13 + f14)
        P2 =  o13 * f21
        P3 =  o13 * f31
        P4 =  o13 * f41

     else:

        G = 0.0 

     Tint = Tint + G * (Y1*P1 + Y2*P2 + Y3*P3 + Y4*P4) 

  dos_gam = Tint* voll

  return dos_gam  # #2 because DOS_ee is per 1 spin


 def find_ef(self):
  print('find EF for chosen degauss')
  minE=np.min(self.ENE[0])
  bottomE=0
  dE=0.01
  oldN=[999,999,999]
  emaxes=np.linspace(bottomE,.2,50)
  Es=np.arange(minE,bottomE,dE)
  DOS=[max(np.sum([np.array([w0gauss(E-j) for j in i])*self.kweights for i in self.ENE]),.1) for E in Es  ]
  Nstart=np.trapz(DOS,Es) 
  plt.plot(Es,DOS)
  plt.show()
  print(minE,Nstart,DOS[-1])
  for maxE in emaxes:
   Es=np.arange(bottomE,maxE,dE)
   DOS=[ np.sum([np.array([w0gauss(E-j) for j in i])*self.kweights for i in self.ENE]) for E in Es  ]
   N=Nstart+np.trapz(DOS,Es) 
   dif=abs(N-self.nel)
   if dif<oldN[1]: oldN=[maxE,dif,N]
  print('Ef,dif,N=',oldN)
  self.ENE=np.array(self.ENE)-0.01 #oldN[0]
  self.ENE_fs=np.array(self.ENE_fs)-0.01 #oldN[0]

 def read_el_structure(self):
  print(' read electronic structure')
  tree = ET.parse(self.tmp_dir+'/'+self.prefix+'.xml')
  root = tree.getroot()

  if_soc=root.find('input/spin/spinorbit').text
  if if_soc=='false': 
   self.soc=0
   print('no spinorbit detected')
  elif if_soc=='true':
   self.soc=1
   print('spinorbit detected')
  

  for i in root.findall('output/band_structure/fermi_energy'):
   self.ef=float(i.text.split()[0])
  print(self.ef)

  self.kweights=[]
  for i in root.findall('output/band_structure/ks_energies'):
    for actor in i.findall('eigenvalues'):
     self.ENE.append([(float(m)-self.ef)*hartree2ry  for m in  actor.text.split()])
    self.kweights.append(float(i.find('k_point').attrib['weight']))

  self.kweights=np.array(self.kweights)

  self.nel=float(root.find('output/band_structure/nelec').text.split()[0])

  maks=max([ len(i) for i in self.ENE])
  ENE2= [ [] for i in range(maks)]
  below_ENE=[ 0 for i in range(maks)]
  top_ENE=[ 0 for i in range(maks)]

#choose only complete bands
  for i in self.ENE:
   for (numj,j) in enumerate(i):
    ENE2[numj].append(j)
    if j<0: below_ENE[numj]=1
    elif j>0: top_ENE[numj]=1
#print('\nList of bands, which cross the Fermi level and its energy ranges:')
  for i in range(maks):
   if below_ENE[i] and top_ENE[i]: 
    self.bands_num.append(i)
#  self.bands_num=self.bands_num[::2]
  self.minband,self.maxband=self.bands_num[0],self.bands_num[-1]
  self.fermi_nbnd_el=len(self.bands_num)

  self.ENE=ENE2 
  self.ENE_fs=ENE2[self.minband:self.maxband+1]#np.transpose(np.array(ENE2)) #ENE=[] #ENE[i][j] , i - no of band, j-no of kpoint


 def calc_dos(self,structure):
  mini=np.min(self.ENE)
  maxi=np.max(self.ENE)
  ndos=1000
  dE=(maxi-mini)/ndos
  DOS=[]
  for m in range(ndos):
    E=mini+m*dE
    DOS.append([E,np.sum( [ w0gauss(E-np.array(band)*np.array(structure.WK)) for band in self.ENE])])
  h=open('eldos.dat','w')
  for i in DOS:

    
   h.write("{:.4} {:.6}\n".format(*i))


 def read_dense_ene(self):
  h=open(self.tmp_dir+'/'+self.prefix+'.a2Fsave')
  tmp=[i.split() for i in h.readlines()]
  h.close()
  self.ENE_dense=[]
  KPOINTS=[]
  WK=[]

  [nbnd_dense,nk_dense]=[int(m) for m in tmp[0]]
  for nl,line in enumerate(tmp[1:]):
   for m in line:
    self.ENE_dense.append(float(m))
   if len(self.ENE_dense)==nbnd_dense*nk_dense:
    break
  print(self.ef)
  self.ENE_dense=(np.array(self.ENE_dense).reshape((nk_dense,nbnd_dense)).transpose()-self.ef*hartree2ry) #*hartree2ry
  self.ENE_fs_dense=self.ENE_dense[self.minband:self.maxband+1]#np.transpose(np.array(ENE2)) #ENE=[] #ENE[i][j] , i - no of band, j-no of kpoint

  for nl2,line in enumerate(tmp[2+nl:]): 
    KPOINTS.append([float(m) for m in line]+[nl2])
    if len(KPOINTS)==nk_dense: break
  for nl3,line in enumerate(tmp[3+nl+nl2:]): 
    for m in line:
     WK.append(float(m))
    if len(WK)==nk_dense: break
  nk_dense=[int(m) for m in tmp[4+nl+nl2+nl3]]  
  return KPOINTS,WK,nk_dense
 


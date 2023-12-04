import xml.etree.ElementTree as ET
from numpy import exp
import numpy as np
from matplotlib import pyplot as plt
Ha_to_ev=13.605662285137*2
hartree2ry=2
def w0gauss(x):
  degauss=0.04
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
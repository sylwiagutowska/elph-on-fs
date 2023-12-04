Ha_to_ev=13.605662285137*2
PRECIS=6
import xml.etree.ElementTree as ET
import numpy as np
from operator import itemgetter
import glob

class structure():
 def __init__(self):
  self.allk=[]
  self.NONEQ=[]
  self.NONEQ_cryst=[]
  self.WK=[] #weights of kpoints
  self.pm=[0,1,-1]
  self.at=[]
  self.e=[]
  self.SYMM=[]
  self.no_of_kpoints=[]
  self.prefix=''
  self.prefixdyn=''
  self.at_pos=[]
  self.at_masses=[]
  self.at_names=[]
  self.irt=[] #we rotate atom by sym and check which atom he becomes
  self.rtau=[] #pos of atoms rotated by sym
  self.tmp_dir='../tmp_dir'
  self.scfout='../scfout' 

 def calc_noneq_cryst(self):
     einv=np.round(np.linalg.inv(np.transpose(self.e)),PRECIS)
     self.NONEQ_cryst=[list(np.round(np.matmul(einv,self.NONEQ[mmm][:3]),PRECIS))+[mmm] for mmm in range(len(self.NONEQ))]


 def read_structure(self):
  print(' read structure...')
  self.prefix=glob.glob(self.tmp_dir+'/*.xml')[0].replace('/','.').split('.')[-2]
  self.prefixdyn=glob.glob('../*dyn0')[0].split('.dyn')[0]
  print(self.prefixdyn)
  tree = ET.parse(self.tmp_dir+'/'+self.prefix+'.xml')
  root = tree.getroot()

  for i in root.findall('output/band_structure/starting_k_points/monkhorst_pack'):
   self.no_of_kpoints=[ int(i.get('nk1')),int(i.get('nk2')),int(i.get('nk3'))]
  print('Monkhorst grid of '+str( self.no_of_kpoints[0])+' '+str( self.no_of_kpoints[1])+' '+str( self.no_of_kpoints[2])+' points')



  #reciprocal lattice vectors
  for i in root.findall('output/basis_set/reciprocal_lattice'):
   self.e.append([round(float(m),PRECIS) for m in i.find('b1').text.split()])
   self.e.append([round(float(m),PRECIS) for m in i.find('b2').text.split()])
   self.e.append([round(float(m),PRECIS) for m in i.find('b3').text.split()])
  self.e=np.round(np.array(self.e),PRECIS)
  print (self.e)

  # lattice vectors
  for i in root.findall('output/atomic_structure/cell'):
   self.at.append([round(float(m),PRECIS) for m in i.find('a1').text.split()])
   self.at.append([round(float(m),PRECIS) for m in i.find('a2').text.split()])
   self.at.append([round(float(m),PRECIS) for m in i.find('a3').text.split()])
  self.at=np.array(self.at)
  print (self.at)

  mmm=0
  for i in root.findall('output/band_structure/ks_energies'):
    for actor in i.findall('k_point'):
     self.NONEQ.append([round(float(m),PRECIS) for m in  actor.text.split()])
     self.NONEQ[-1].append(mmm)
     mmm=mmm+1
  self.calc_noneq_cryst()
  print ('No of nonequivalent kpoints= '+str(len(self.NONEQ)))

  # atomic positions
  for i in root.findall('output/atomic_structure/atomic_positions/atom'):
   self.at_pos.append([round(float(m),PRECIS) for m in i.text.split()])
   self.at_names.append(i.get('name'))
  self.at_pos=np.array(self.at_pos)
  self.at_pos_crystal=[np.round(np.dot(np.linalg.inv(np.transpose(self.at)),i),PRECIS) for i in self.at_pos]
  self.at_pos_crystal_orig=[np.round(np.dot(np.linalg.inv(np.transpose(self.at)),i),PRECIS) for i in self.at_pos]
  for i in self.at_pos_crystal: 
    for m in range(3): 
     while i[m]>=1: i[m]=i[m]-1.
     while i[m]<0: i[m]=i[m]+1.
  print (self.at_pos) #,self.at_names)

  masses,names=[],[]
  #atomic masses
  for i in root.findall('input/atomic_species/species'):
   masses.append(round(float(i.find('mass').text),PRECIS)) #of spieces, not of positions
   names.append(i.get('name'))
  #print(masses,names)
  for i in self.at_names: # all positions
   for nj,j in enumerate(names): #all species
    if i==j: self.at_masses.append(masses[nj])

  #symmetry operations
  self.SYMM_crystal=[]
  self.SYMM=[]
  self.FRACTIONAL_TRANSLATION_crystal=[]
  self.FRACTIONAL_TRANSLATION_r=[]
  self.FRACTIONAL_TRANSLATION_k=[]
 # self.SYMM2=[]
  einv=np.linalg.inv(self.e)
 # atinv=np.linalg.inv(self.at)
  for neighbor in root.iter('rotation'):
     tmp=neighbor.text.split()
     tmp2=np.array([ [ float(m) for m in tmp[0:3]], [float(m) for m in tmp[3:6]], [float(m) for m in tmp[6:9]]])
     self.SYMM.append(np.transpose(np.dot(einv,(np.dot(tmp2,self.e)))))
     self.SYMM_crystal.append(np.round(tmp2,PRECIS))
  for neighbor in root.iter('fractional_translation'):
     tmp=neighbor.text.split()
     tmp2=np.round(np.array( [ float(m) for m in tmp[0:3]]),PRECIS)
     self.FRACTIONAL_TRANSLATION_crystal.append(tmp2)
     self.FRACTIONAL_TRANSLATION_r.append(np.round(np.matmul(tmp2,self.at),PRECIS) ) #in cart --- depends wether in reciprocal or real space
     self.FRACTIONAL_TRANSLATION_k.append(np.round(np.matmul(tmp2,self.e),PRECIS) ) #in cart --- depends wether in reciprocal or real space
  print ('No of symm. op.='+str(len(self.SYMM)))




  self.e=np.array(self.e)
  print (self.e)


 def calc_irt(self):
#subroutine sgam_lr (at, bg, nsym, s, irt, tau, rtau, nat)
#  !-----------------------------------------------------------------------
#  !
#  !     This routine computes the vector rtau which contains for each
#  !     atom and each rotation the vector S\tau_a - \tau_b, where
#  !     b is the rotated a atom, given by the array irt. These rtau are
#  !     non zero only if fractional translations are present.
#  !
# !   compute the atomic coordinates in crystal axis, xau
#  !
#  do na = 1, nat
#     do ipol = 1, 3
#        xau (ipol, na) = bg (1, ipol) * tau (1, na) + &
#                         bg (2, ipol) * tau (2, na) + &
#                         bg (3, ipol) * tau (3, na)
#     enddo
#  enddo
#  !
#  !    for each symmetry operation, compute the atomic coordinates
#  !    of the rotated atom, ft, and calculate rtau = Stau'-tau
#  !
#  rtau(:,:,:) = 0.0_dp
#  do isym = 1, nsym
#     do na = 1, nat
#        nb = irt (isym, na)
#        do ipol = 1, 3
#           ft (ipol) = s (1, ipol, isym) * xau (1, na) + &
#                       s (2, ipol, isym) * xau (2, na) + &
#                       s (3, ipol, isym) * xau (3, na) - xau (ipol, nb)
#        enddo
#	do ipol = 1, 3
#           rtau (ipol, isym, na) = at (ipol, 1) * ft (1) + &
#                                   at (ipol, 2) * ft (2) + &
#                                   at (ipol, 3) * ft (3)
  einv=np.round(np.linalg.inv(np.transpose(self.at)),PRECIS)
  self.irt=[[-1 for i in self.SYMM] for j in self.at_pos]
  self.rtau=[[[] for i in self.SYMM] for j in self.at_pos]
  self.rtau_cart=[[[] for i in self.SYMM] for j in self.at_pos]
  for ni,i in enumerate(self.at_pos_crystal_orig):
   for nj,j in enumerate(self.SYMM_crystal):
    pos22=np.round(np.dot(j,i)+self.FRACTIONAL_TRANSLATION_crystal[nj],PRECIS) #in cryst
    pos2=np.round(pos22,PRECIS-1) #in cryst
    self.rtau[ni][nj]=np.array([m for m in pos2])
    for m in range(3): 
     while pos2[m]>=1: pos2[m]-=1.
     while pos2[m]<0:  pos2[m]+=1.
    for nk,k in enumerate(self.at_pos_crystal):  
      pos2=np.round(pos2,PRECIS-1)
      k2=np.round(k,PRECIS-1)
      if pos2[0]==k2[0] and pos2[1]==k2[1] and pos2[2]==k2[2]:
       self.irt[ni][nj]=nk
       self.rtau[ni][nj]=np.round(self.rtau[ni][nj]-self.at_pos_crystal_orig[nk], PRECIS-1)
       break
    if self.irt[ni][nj]==-1: raise ValueError('WRONG SYMMETRY: IRT NOT FOUND',j,pos2,pos22,self.at_pos_crystal)
 # print(self.irt,self.rtau)

 def sorting(self,allk2):
  Xall=[]
  allk2=sorted(allk2, key=itemgetter(0))
  i=0
  while i<len(allk2): 
   X=[]
   x=allk2[i]
   while i<len(allk2) and x[0]==allk2[i][0]:
    X.append(allk2[i])
    i=i+1
   if len(X)>1: X=sorted(X, key=itemgetter(1))
   Xall.append(X)

  Yall=[]
  for X in Xall:
   x=X[0]
   i=0
   while i<len(X): 
    Y=[]
    x=X[i]
    while i<len(X) and x[1]==X[i][1]:
     Y.append(X[i])
     i=i+1
    Y=sorted(Y, key=itemgetter(2))
    Yall.append(Y)

  allk=[]
  for i in Yall:
   for j in i:
    allk.append(j)
  
#  print ' Sorting - Done!'
  return allk

 def remove_repeated_items(self,allk2):
  found=1
  while found:
   allk=[]
   for i in range(len(allk2)-1):
    x=allk2[i]
    y=allk2[i+1]
    if not((x[0]==y[0]) and (x[1]==y[1]) and (x[2]==y[2])):
     allk.append(x)
   if not((allk[-1][0]==allk2[-1][0]) and (allk[-1][1]==allk2[-1][1]) and (allk[-1][2]==allk2[-1][2])): allk.append(allk2[-1])
   if len(allk)==len(allk2): found=0
   else: allk2=list(allk)
  return allk
 

 def calc_weight_of_k(self):
  WK=[0 for i in range(len(self.NONEQ))]
  for i in range(len(self.NONEQ)):
   for j in self.allk:
    if j[3]==i:
     WK[i]+=1
  return WK #/len(self.allk)

 def make_kgrid(self,q=-1):
  print(' make whole kgrid')
  allk2=[]
  allk_in_crystal_coordinates=[]
  einv=np.round(np.linalg.inv(np.transpose(self.e)),PRECIS)
  no=0
  for nq in self.NONEQ_cryst:
  # print nq
#   nq2=[ round(sum([nq[m]*self.e[m][m2] for m in range(3)]),PRECIS)\
#                      for m2 in range(3)] #transform  to cartesian coordinates
   nq2=np.round(np.multiply(nq[:3],self.no_of_kpoints))
   for ns,sym in enumerate(self.SYMM_crystal):
    x=np.round(np.dot(np.transpose(sym),nq2[:3])) #[sum([sym[m1][m2]*nq[m2] for m2 in range(3)]) for m1 in range(3)]
    for mm in range(3):
      while x[mm]>=self.no_of_kpoints[mm]: x[mm]= x[mm]-self.no_of_kpoints[mm]
      while x[mm]<0: x[mm]= x[mm]+self.no_of_kpoints[mm]
    x=np.round(x/self.no_of_kpoints,PRECIS)
    k_point2=list(np.round(np.dot(x,self.e),PRECIS-1)) #[round(kk,PRECIS) for kk in 
                 #x[0]*self.e[0]+x[1]*self.e[1]+x[2]*self.e[2]] #from crystal to cartesian
  #  x=np.round(np.multiply(x,self.no_of_kpoints))
    allk2.append(k_point2)
    allk2[-1].append(nq[3])
    allk2[-1].append(no)
    allk_in_crystal_coordinates.append(list(x))
   # allk_in_crystal_coordinates[-1].append(nq[3])
    no+=1
 # print (len(allk2))
  allk2=self.sorting(allk2)
  self.allk=allk2
  #rearrange
  allk_in_crystal_coordinates=[ allk_in_crystal_coordinates[v[4]]+[v[3],nv] for nv,v in enumerate(self.allk)]
  allk_in_crystal_coordinates=self.sorting(allk_in_crystal_coordinates)
  self.allk_in_crystal_coordinates=self.remove_repeated_items( allk_in_crystal_coordinates )
  self.allk=[ self.allk[i[4]] for i in self.allk_in_crystal_coordinates]
  self.allk=[i[:4] for i in self.allk]
  self.allk_in_crystal_coordinates=[i[:4] for i in self.allk_in_crystal_coordinates]
  for ni,i in enumerate(self.allk):
   self.allk_in_crystal_coordinates[ni][3]=i[3]

  #for i in self.allk_in_crystal_coordinates: 
  # for m in range(3):
  #  i[m]=np.round(i[m]/self.no_of_kpoints[m],PRECIS-2)


  #calc weights of k
  self.WK=self.calc_weight_of_k()

  h=open('kpoints'+str(q)+'.dat','w')
  h.write(str(len(self.SYMM_crystal))+'\n')
  for ni,i in enumerate(self.NONEQ_cryst):
   for j in range(4):
    h.write(str(i[j])+' ')
   h.write('\n')
  h.close()

  if len(self.allk)!=self.no_of_kpoints[0]*self.no_of_kpoints[1]*self.no_of_kpoints[2]: 
  # for i in range(0,self.no_of_kpoints[0]):  
  #  for j in range(0,self.no_of_kpoints[1]):
  #   for k in range(0,self.no_of_kpoints[2]):
  #    found=0
  #    for kk in self.allk_in_crystal_coordinates:
  #     if i==np.round(kk[0]*self.no_of_kpoints[0],PRECIS-2) and j==np.round(kk[1]*self.no_of_kpoints[1],PRECIS-2) and k==np.round(kk[2]*self.no_of_kpoints[2],PRECIS-2): 
  #      found=1
  #      break
  #    if not found: print(q,': kpoint',i,j,k,'not found')
   raise ValueError(q,'wrong no of kpoints',len(self.allk),'!=',self.no_of_kpoints[0]*self.no_of_kpoints[1]*self.no_of_kpoints[2])



   

  '''
  h=open('kpoints.dat','w')
  for i in self.allk:
   for j in i:
    h.write(str(j)+' ')
   h.write('\n')
  h.close()



  h=open('kpoints_noneq.dat','w')
  for i in self.NONEQ:
   for j in i:
    h.write(str(j)+' ')
   h.write('\n')
  h.close()
  '''

 def check_symm(self,q_cryst,basic_kpoints,no_of_s_q,q_no):
 # SYMM2_crystal=[]
 # SYMM2=[]
  founded=[]
  einv=np.linalg.inv(np.transpose(self.e))
  q2=np.round(q_cryst,PRECIS) #transform from cartesian to crystal coordinates, then we have a cube of points 
  for i in range(3):
      while q2[i]>=1: q2[i]=q2[i]-1
      while q2[i]<0: q2[i]=q2[i]+1
  q1=np.round(q2,PRECIS-1)
  q2=np.round(q2,PRECIS-2)
  for si,sym in enumerate(self.SYMM_crystal):
     x=np.round(np.matmul(np.transpose(sym),q1),PRECIS) #[sum([sym[m1][m2]*q[m2] for m2 in range(3)]) for m1 in range(3)]   #in cryst
     for i in range(3):
      while x[i]>=1: x[i]=x[i]-1
      while x[i]<0: x[i]=x[i]+1
     x=np.round(x,PRECIS-2)
     if x[0]==q2[0] and  x[1]==q2[1] and  x[2]==q2[2]:
           founded.append(si)
#      SYMM2.append(sym)
#      SYMM2_crystal.append(self.SYMM_crystal[si])
  print(q_no,q_cryst,':',len(founded),'operations of symmetry')
  if len(founded)!=no_of_s_q: raise ValueError(q_no,' no of sym is wrong', len(founded), no_of_s_q)
  self.SYMM=[self.SYMM[si] for si in founded]
  self.SYMM_crystal=[self.SYMM_crystal[si] for si in founded]
  self.FRACTIONAL_TRANSLATION_crystal=[self.FRACTIONAL_TRANSLATION_crystal[si] for si in founded]
  self.FRACTIONAL_TRANSLATION_k=[self.FRACTIONAL_TRANSLATION_k[si] for si in founded]
  self.FRACTIONAL_TRANSLATION_r=[self.FRACTIONAL_TRANSLATION_r[si] for si in founded]

 def find_k_plus_q(self,k,allk,allk_cryst,q_cryst):
 # print(allk_cryst)
  kpq_no=-2
  #einv=np.linalg.inv(np.transpose(self.e))
  k2=k
  q2=q_cryst
  kpq0=np.round([k2[i]+q2[i] for i in range(3)],PRECIS-1) #in cryst coords
  for i in range(3): 
   while kpq0[i]>=1: kpq0[i]=kpq0[i]-1
   while kpq0[i]<0: kpq0[i]=kpq0[i]+1
  kpq0=np.round(kpq0,PRECIS-1) #in cryst coords
#  kpq=[round(sum([kpq0[m]*self.e[m][m2] for m in range(3)]),PRECIS-1) for m2 in range(3)] #in cart
  found=0
  for sym in self.SYMM_crystal:
   if found==1: break
   kpq2=np.round(np.matmul(np.transpose(sym),kpq0),PRECIS-2) #in c coord
   for i in range(3): 
    while kpq2[i]>=1: kpq2[i]=kpq2[i]-1
    while kpq2[i]<0: kpq2[i]=kpq2[i]+1
   kpq2=np.round(kpq2,PRECIS-2)
   for ki in allk_cryst:
      if found==1: break
      if kpq2[0]==ki[0] and kpq2[1]==ki[1] and kpq2[2]==ki[2]:
            kpq_no=ki[3]
            kpq_no_in_basic=ki[4]
            found=1
            break
  if kpq_no==-2: raise ValueError(q2,k2,kpq0,kpq2,' kpq not found',allk_cryst)
  return [kpq2,kpq_no,kpq_no_in_basic]

 def find_newkpoints_in_old_list(self,old_allk):
  for i in self.NONEQ: i.append(None)
  for i in self.NONEQ_cryst: i.append(None)
  for nk,k in enumerate(self.allk):
   k.append(old_allk[nk][3]) # k=kx,ky,kz,no_of_noneq_in_new_grid,no_of_noneq_in_basic_grid
   self.allk_in_crystal_coordinates[nk].append(old_allk[nk][3])
   self.NONEQ[k[3]][4]=old_allk[nk][3]
   self.NONEQ_cryst[k[3]][4]=old_allk[nk][3]

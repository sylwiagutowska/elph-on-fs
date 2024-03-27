import xml.etree.ElementTree as ET
import numpy as np
import structure,ph_structure,el_structure
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interpn
from multiprocessing import Process,Pool
import copy
from elph_structure_single_q import *
PRECIS=6
RY_TO_THZ=3289.8449

#def symmetrize(a):
"""
    Return a symmetrized version of NumPy array a.
    Values 0 are replaced by the array value at the symmetric
    position (with respect to the diagonal), i.e. if a_ij = 0,
    then the returned array a' is such that a'_ij = a_ji.
    Diagonal values are left untouched.
    a -- square NumPy array, such that a_ij = 0 or a_ji = 0, 
    for i != j.
"""
 #   return a + a.T - np.diag(a.diagonal())
#
class elph_structure():
 def __init__(self,phh_structure,lambda_or_a2f):
  self.ALL_COLORS=[[] for i in range(len(phh_structure.Q))] 
  self.KPOINTS_all_all_q=[]
  self.SUMMED_COLORS=[]
  self.lambda_or_a2f=lambda_or_a2f
 def parallel_job(self,structure,structure_dense,ell_structure,phh_structure):
  nq=len(phh_structure.Q)
  no_of_pool=nq
  with Pool(no_of_pool) as pol:
   results=pol.map(self.single_job,
              [[int(q+1),structure,structure_dense,phh_structure,ell_structure,self.lambda_or_a2f] for q in range(nq)])
   self.ALL_COLORS=[ i[0] for i in results]
   self.KPOINTS_all_all_q=[ i[1] for i in results]
   self.lambda_from_files_all=[i[2] for i in results]
   self.ALL_COLORS2=[ i[3] for i in results]
 def single_job(self,args):
  [q,structure,structure_dense,phh_structure,ell_structure,lambda_or_a2f]=args
  print('calculations for '+str(q)+'. of total '+str(len(phh_structure.Q))+' q points')
  elph_q=elph_structure_single_q(phh_structure,q )
  elph_q.make_kpoints_single_q(structure,phh_structure.Q[q-1],phh_structure.Q_crystal[q-1],phh_structure.qstar[q-1])
  elph_q.read_elph_single_q(phh_structure,ell_structure,structure)
 # elph_q.make_kpoints_dense_single_q(structure_dense,phh_structure.Q[q-1],phh_structure.Q_crystal[q-1],phh_structure.qstar[q-1])
 # elph_q.interpolate_elph_single_q(phh_structure,ell_structure,structure)
  #print(elph_q.ELPH.shape)
  elph_q.sum_elph_over_jband(phh_structure,ell_structure,structure)


  elph_q.calc_lambda_q_as_elphon_does(phh_structure,ell_structure,structure)
  elph_q.calc_lambda_summed_over_modes(structure,\
                phh_structure,ell_structure,lambda_or_a2f)  #'lambda' or 'elph')
  elph_q.elph_single_q_in_whole_kgrid(structure,\
                phh_structure,ell_structure,lambda_or_a2f)  #'lambda' or 'elph'
  return [elph_q.COLORS,elph_q.KPOINTS_all,elph_q.lambdas_from_file,elph_q.COLORS2]
  print(str(q)+'. point ended')
  


 def sum_over_q(self,phh_structure,structure,ell_structure):

  print('summing over q...')
  self.SUMMED_COLORS=np.array([[0 for k in range(len(structure.allk))] for jband in range(ell_structure.fermi_nbnd_el)])
  self.SUMMED_COLORS_nq=[[0 for k in range(len(structure.NONEQ))] for jband in range(ell_structure.fermi_nbnd_el)]

  suma=0
  for numq,col_q in enumerate(self.ALL_COLORS2):
   suma+=np.sum(col_q)*phh_structure.multiplicity_of_qs[numq]
  print('first summed over k',suma/sum(phh_structure.multiplicity_of_qs))

  


  '''
  for numq,col_q in enumerate(self.ALL_COLORS):
   for band in range(0,len(self.SUMMED_COLORS),2):
    col_q[band]=(col_q[band]+col_q[band+1])/2
    col_q[band+1]=col_q[band]
  
   for band in range(len(self.SUMMED_COLORS)):
     for numk,k in enumerate(structure.allk):
      col_q[band][k[3]]+=col_q[band][numk] /structure.WK[k[3]] 
   for band in range(len(col_q)):
     for numk,k in enumerate(structure.allk):
      col_q[band][numk]=col_q[band][k[3]]
  '''

  for band in range(len(self.SUMMED_COLORS)):
   for numk,k in enumerate(structure.allk):
    for numq,col_q in enumerate(self.ALL_COLORS):
     self.SUMMED_COLORS[band][numk]+=col_q[band][numk] *phh_structure.multiplicity_of_qs[numq] #*el_structure.w0gauss(-ell_structure.ENE_fs[band][k[3]]) 
 
  #imposing simmetry
  if ell_structure.soc==1:
   for band in range(0,len(self.SUMMED_COLORS),2):
    self.SUMMED_COLORS[band]=(self.SUMMED_COLORS[band]+self.SUMMED_COLORS[band+1])/2
    self.SUMMED_COLORS[band+1]=self.SUMMED_COLORS[band]
  
  for band in range(len(self.SUMMED_COLORS)):
   for numk,k in enumerate(structure.allk):
     self.SUMMED_COLORS_nq[band][k[3]]+=self.SUMMED_COLORS[band][numk] /structure.WK[k[3]] 
  for band in range(len(self.SUMMED_COLORS)):
   for numk,k in enumerate(structure.allk):
     self.SUMMED_COLORS[band][numk]=self.SUMMED_COLORS_nq[band][k[3]] #el_structure.w0gauss(-ell_structure.ENE_fs[band][k[3]]) 
 
 
  self.SUMMED_COLORS=self.SUMMED_COLORS/sum(phh_structure.multiplicity_of_qs) #/sum_wag #divide 2 because of spin
  ##verification
  sum_wag=0
  lamsum=0 
  lamsum2=0
  for band in range(len(self.SUMMED_COLORS)):
   for numk,k in enumerate(structure.allk):
    waga=el_structure.w0gauss(-ell_structure.ENE_fs[band][k[3]]) 
    sum_wag+=waga
    lamsum+=self.SUMMED_COLORS[band][numk] *waga  #*waga #--done for every q
    lamsum2+=self.ALL_COLORS[5][band][numk]*waga #/structure.WK[k[3]]  #*waga #--done for every q

  #if ell_structure.soc==1: sum_wag=sum_wag/2 #devide by 2 because we are suming not averaging over spin
  lamsum=lamsum/sum_wag 
  print('calc as elph',lamsum2/sum_wag)
  print ('summed lambda=',lamsum)
  print ('summed wagi=',sum_wag)
  lambda_from_file=np.sum(np.sum(np.array(self.lambda_from_files_all,dtype=float),axis=1)*phh_structure.multiplicity_of_qs/sum(phh_structure.multiplicity_of_qs))
  print('summed lambda from file=',lambda_from_file)
  dos=sum_wag/len(structure.allk)
  print ('summed dos=',dos,' 1/Ry = ',dos/13.606,' 1/eV')
 # print ('summed lambda/wagi=',lam/sum_wag/sum(phh_structure.multiplicity_of_qs))
  h=open('lambda_info','w')
  h.write('calc as elph'+str(lamsum2/sum_wag))
  h.write ('summed lambda='+str(lamsum))
  h.write ('summed wagi='+str(sum_wag))
  h.write('summed lambda from file='+str(lambda_from_file))
  h.write ('summed dos='+str(dos)+' 1/Ry = '+str(dos/13.606)+' 1/eV')
  h.close()

  self.SUMMED_COLORS_dense=[] #np.array([[0 for k in range(8*len(structure.allk))] for jband in range(ell_structure.fermi_nbnd_el)])
  xyz=[np.linspace(0,1,m) for m in structure.no_of_kpoints]
  for bnd in self.SUMMED_COLORS:
   bnd2=np.zeros(structure.no_of_kpoints) 
   m=0
   for i in range(structure.no_of_kpoints[0]):
    for j in range(structure.no_of_kpoints[1]):
     for k in range(structure.no_of_kpoints[2]):
      bnd2[i][j][k]=bnd[m]
      m+=1



   mult=2
#   bnd2=bnd.reshape(structure.no_of_kpoints)
   fn = RegularGridInterpolator((xyz[0],xyz[1],xyz[2]), bnd2)
   nk=np.array(structure.no_of_kpoints)*mult
   xyz2=[ [i/nk[0],j/nk[1],k/nk[2]] for k in range(nk[2]) for j in range(nk[1]) for i in range(nk[0])]
   bnd3_interp=fn(xyz2)
   self.SUMMED_COLORS_dense.append(bnd3_interp)

  
  ##save to file
  if self.lambda_or_elph=='a2f':
   h=open('a2f.frmsf','w')
  else:
   h=open('lambda.frmsf','w')
  h.write(str(structure.no_of_kpoints[0])+' '+str(structure.no_of_kpoints[1])+' ' +str(structure.no_of_kpoints[2])+'\n')
  h.write('1\n'+str(len(ell_structure.bands_num))+'\n')
  for i in structure.e:
   for j in i:
    h.write(str(j)+' ')
   h.write('\n')
  for bnd in ell_structure.ENE_fs:#_dense:
   for k in structure.allk:# _dense:
    h.write(str(bnd[k[3]])+'\n')
  for bnd in self.SUMMED_COLORS: #_dense:
   for k in bnd:
    h.write(str(np.abs(k))+'\n')
  h.close()


  ##save to file
  if self.lambda_or_elph=='a2f':
   h=open('elph.frmsf','w')
  else:
   h=open('lambda_dense.frmsf','w')
  h.write(str(nk[0])+' '+str(nk[1])+' ' +str(nk[2])+'\n')
  h.write('1\n'+str(len(ell_structure.bands_num))+'\n')
  for i in structure.e:
   for j in i:
    h.write(str(j)+' ')
   h.write('\n')
  for bnd in ell_structure.ENE_fs_dense:
   for k in structure.allk_dense:
    h.write(str(bnd[k[3]])+'\n')
  for bnd in self.SUMMED_COLORS_dense:
   for k in bnd:
    h.write(str(np.abs(k))+'\n')
  h.close()
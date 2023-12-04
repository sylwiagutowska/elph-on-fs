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
 def __init__(self,phh_structure,lambda_or_elph):
  self.ALL_COLORS=[[] for i in range(len(phh_structure.Q))] 
  self.KPOINTS_all_all_q=[]
  self.SUMMED_COLORS=[]
  self.lambda_or_elph=lambda_or_elph
 def parallel_job(self,structure,ell_structure,phh_structure):
  nq=len(phh_structure.Q)
  no_of_pool=nq
  with Pool(no_of_pool) as pol:
   results=pol.map(self.single_job,
              [[int(q+1),structure,phh_structure,ell_structure,self.lambda_or_elph] for q in range(nq)])
   self.ALL_COLORS=[ i[0] for i in results]
   self.KPOINTS_all_all_q=[ i[1] for i in results]
   self.lambda_from_files_all=[i[2] for i in results]
 def single_job(self,args):
  [q,structure,phh_structure,ell_structure,lambda_or_elph]=args
  print('calculations for '+str(q)+'. of total '+str(len(phh_structure.Q))+' q points')
  elph_q=elph_structure_single_q(phh_structure,q )
  elph_q.make_kpoints_single_q(structure,phh_structure.Q[q-1],phh_structure.Q_crystal[q-1],phh_structure.qstar[q-1])
  elph_q.read_elph_single_q(phh_structure,ell_structure,structure)
  elph_q.sum_elph_over_jband(phh_structure,ell_structure,structure)
  elph_q.calc_lambda_q_as_elphon_does(phh_structure,ell_structure,structure)
  elph_q.calc_lambda_summed_over_modes(structure,\
                phh_structure,ell_structure,lambda_or_elph)  #'lambda' or 'elph')
  elph_q.elph_single_q_in_whole_kgrid(structure,\
                phh_structure,ell_structure,lambda_or_elph)  #'lambda' or 'elph'
  return [elph_q.COLORS,elph_q.KPOINTS_all,elph_q.lambdas_from_file]
  print(str(q)+'. point ended')
  


 def sum_over_q(self,phh_structure,structure,ell_structure):

  print('summing over q...')
  self.SUMMED_COLORS=np.array([[0 for k in range(len(structure.allk))] for jband in range(ell_structure.fermi_nbnd_el)])
  self.SUMMED_COLORS_nq=[[0 for k in range(len(structure.NONEQ))] for jband in range(ell_structure.fermi_nbnd_el)]

  for band in range(len(self.SUMMED_COLORS)):
   for numk,k in enumerate(structure.allk):
    for numq,col_q in enumerate(self.ALL_COLORS):
     self.SUMMED_COLORS[band][numk]+=col_q[band][numk] *phh_structure.multiplicity_of_qs[numq]  #*waga #--done for every q


  #imposing simmetry
  if ell_structure.soc==1:
   for band in range(0,len(self.SUMMED_COLORS),2):
    self.SUMMED_COLORS[band]=(self.SUMMED_COLORS[band]+self.SUMMED_COLORS[band+1])/2
    self.SUMMED_COLORS[band+1]=self.SUMMED_COLORS[band]

  for band in range(len(self.SUMMED_COLORS)):
   for numk,k in enumerate(structure.allk):
     self.SUMMED_COLORS_nq[band][k[3]]+=self.SUMMED_COLORS[band][numk] /structure.WK[k[3]]
  #self.SUMMED_COLORS=[[0 for k in range(len(structure.allk))] for jband in range(ell_structure.fermi_nbnd_el)]
  for band in range(len(self.SUMMED_COLORS)):
   for numk,k in enumerate(structure.allk):
     self.SUMMED_COLORS[band][numk]=self.SUMMED_COLORS_nq[band][k[3]]
 
  self.SUMMED_COLORS=self.SUMMED_COLORS/2/sum(phh_structure.multiplicity_of_qs) #/sum_wag #divide 2 because of spin
  sum_wag=0
  lamsum=0
#  h=open('aa.frmsf','w')
#  h.write("12 12 12\n1\n4\n-1.0 -1.0 1.0\n1.0 1.0 1.0\n-1.0 1.0 -1.0\n")
  for band in range(len(self.SUMMED_COLORS)):
   for numk,k in enumerate(structure.allk):
    waga=el_structure.w0gauss(-ell_structure.ENE_fs[band][k[3]])
 #   h.write(str(ell_structure.ENE_fs[band][k[3]])+'\n')
    sum_wag+=waga
    lamsum+=self.SUMMED_COLORS[band][numk]*waga
  lamsum=lamsum/sum_wag
#  h.close()
#  lam=sum([ sum(i) for i in self.SUMMED_COLORS])
  print ('summed lambda=',lamsum)
  print ('summed wagi=',sum_wag)
  print('summed lambda from file=',np.sum(np.sum(np.array(self.lambda_from_files_all,dtype=float),axis=1)*phh_structure.multiplicity_of_qs/sum(phh_structure.multiplicity_of_qs)))
  dos=sum_wag/len(structure.allk)
  print ('summed dos=',dos,' 1/Ry = ',dos/13.606,' 1/eV')
 # print ('summed lambda/wagi=',lam/sum_wag/sum(phh_structure.multiplicity_of_qs))
  if self.lambda_or_elph=='elph':
   h=open('elph.frmsf','w')
  else:
   h=open('lambda.frmsf','w')
  h.write(str(structure.no_of_kpoints[0])+' '+str(structure.no_of_kpoints[1])+' ' +str(structure.no_of_kpoints[2])+'\n')
  h.write('1\n'+str(len(ell_structure.bands_num))+'\n')
  for i in structure.e:
   for j in i:
    h.write(str(j)+' ')
   h.write('\n')
  for bnd in ell_structure.ENE_fs:
   for k in structure.allk:
    h.write(str(bnd[k[3]])+'\n')
  for bnd in self.SUMMED_COLORS:
   for k in bnd:
    h.write(str(np.abs(k))+'\n')
  h.close()



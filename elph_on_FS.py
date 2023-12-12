
import xml.etree.ElementTree as ET
import numpy as np
from operator import itemgetter
import structure
import el_structure
import ph_structure
import elph_structure
import copy

 
ELECTRONMASS_SI  = 9.10938215e-31   # Kg
AMU_SI           = 1.660538782e-27  #Kg
AMU_AU           = AMU_SI / ELECTRONMASS_SI
AMU_RY           = AMU_AU / 2. #=911.44 
RY_TO_THZ=3289.8449


def round_complex(y):
 prec=9
 try:
  for i in range(len(y)):
   for j in range(len(y[i])):
    y[i][j]=round(y[i][j].real,prec)+1j*round(y[i][j].imag,prec)
 except:
  try:
   for i in range(len(y)):
    y[i]=round(y[i].real,prec)+1j*round(y[i].imag,prec)
  except:
    y=round(y.real,prec)+1j*round(y.imag,prec)
 return y

sstructure=structure.structure()
sstructure.read_structure()
sstructure.make_kgrid()
sstructure.calc_irt()


el_structure=el_structure.el_structure(sstructure)
el_structure.read_el_structure()
#el_structure.calc_dos(sstructure)
#el_structure.find_ef()
#exit()

KPOINTS,WK,nk_dense=el_structure.read_dense_ene()
#el_structure.find_ef()
structure_dense=copy.deepcopy(sstructure)
structure_dense.NONEQ=KPOINTS
structure_dense.calc_noneq_cryst()
structure_dense.no_of_kpoints=nk_dense
structure_dense.make_kgrid()
sstructure.allk_dense=structure_dense.allk
sstructure.allk_dense_in_crystal_coordinates=structure_dense.allk_in_crystal_coordinates

print(len(sstructure.allk_dense),el_structure.ENE_fs_dense.shape)
#print(structure_dense.allk_in_crystal_coordinates)
phh_structure=ph_structure.ph_structure(sstructure)
phh_structure.read_ph_structure()
#ph_structure.check_symm_of_q(structure)
print(phh_structure.multiplicity_of_qs)
print(np.array(el_structure.ENE_dense)[:,0])
print(np.array(el_structure.ENE)[:,0])
#print(phh_structure.PATT[1])
basic_structure=sstructure
#print(max([m[3] for m in sstructure.allk_dense]))


 #check if composing dynmat works - WORKS! GIVES PROPER omega^2


#pb_mass=207.2 #*(9.3975038) #**2
#pbbi_mass=207.2 #*6.6485187 #**2
amass=np.array(basic_structure.at_masses)*AMU_RY #[pbbi_mass,pbbi_mass]
print(amass)

#print(phh_structure.Q_crystal)
#sstructure.calc_irt()
#print(sstructure.irt)
#exit()
jedynki=[]
h=open('macierze_dyn.dat','w')
for qno in range(len(phh_structure.Q)):
 print(qno,phh_structure.Q_crystal[qno])
 structure_new=copy.deepcopy(sstructure)
 structure_new.check_symm(phh_structure.Q_crystal[qno],sstructure.NONEQ,phh_structure.no_of_s_q[qno],qno)
 structure_new.calc_irt()
#print(len(structure_new.SYMM))
# print(phh_structure.DYN2[qno],np.sum(phh_structure.DYN2[qno],axis=0))


 dyn=ph_structure.symmetrize(phh_structure.nat,np.array(phh_structure.PATT[qno]),
     (phh_structure.DYN2[qno]), sstructure.at,sstructure.e,
     sstructure.SYMM_crystal, structure_new.SYMM_crystal,
     structure_new.irt,structure_new.rtau,phh_structure.Q_crystal_orig[qno] )

 # print(round_complex(dyn)) ###now the dyn shoud be equal to this in  prefix.dyn$i file
 
 
 for i in range(3):
  for na in range(phh_structure.nat):
   mu=3*(na)+i
   for j in range(3):
    for nb in range(phh_structure.nat):
     nu=3*(nb)+j
     dyn [mu][nu] = dyn [mu][nu] / ( amass[na]*amass[nb])**0.5
 

 a=np.linalg.eig(dyn)
 dyn=a[1]
 prevfreq=phh_structure.FREQ[qno]
 phh_structure.FREQ[qno]=np.sqrt(np.abs(a[0])) 
 for ni,i in enumerate(a[0]):
  if i<0: phh_structure.FREQ[qno][ni]=0 #-phh_structure.FREQ[qno][ni]
 
 
 #exit()
 #now we are calculating the normal modes
 for nu in range(3*phh_structure.nat):
        for mu in range(3*phh_structure.nat):
         ma=int(mu/3)
  #       print(amass[na]**.5)
         dyn[mu][nu]= dyn[mu][nu] /   ((amass[na])**.5)
 
 
 phh_structure.DYN[qno]=[(dyn)]
 '''
 print(qno,dyn)
 for i in range(3):
  for na in range(phh_structure.nat):
   mu=3*(na)+i
   for j in range(3):
    for nb in range(phh_structure.nat):
     nu=3*(nb)+j
     phh_structure.DYN[qno][0] [mu][nu] = phh_structure.DYN[qno][0] [mu][nu]  / ( amass[na]*amass[nb])**0.5
 
 for nu in range(3*phh_structure.nat):
        for mu in range(3*phh_structure.nat):
         na=int(mu/3)
  #       print(amass[na]**.5)
         phh_structure.DYN[qno][0][nu][mu]= phh_structure.DYN[qno][0][nu][mu] /   ((amass[na])**.5)
 '''
 


 indx=np.argsort(phh_structure.FREQ[qno])
 phh_structure.ORDER_OF_IRR[qno]=indx
 if qno==0: 
  for m in range(3): phh_structure.FREQ[qno][indx[m]]=0
 h.write(str(qno)) 
 for i in phh_structure.DYN[qno][0]:
  for j in i:
   h.write(str(j)+' ')
  h.write('\n')
 jedynki.append([round(i/(prevfreq[ni]/RY_TO_THZ),4) for ni,i in enumerate(sorted(phh_structure.FREQ[qno]))])
 for i in jedynki[-1]:  h.write(str(i)+'\n')

 #if sum(jedynki[-1])!=len(jedynki[-1]): print(structure_new.SYMM_crystal,phh_structure.Q_crystal[qno],phh_structure.Q[qno])
h.close()
 #print(round_complex(a[0]*RY_TO_THZ*RY_TO_THZ))
 #print(round_complex(a[1]))


for ni,i in enumerate(jedynki):
 if sum(i)!=len(i): 
  print('Q',ni+1,phh_structure.Q_crystal[ni],'wrong! freq/freq=',i)
  if ni!=0: print(phh_structure.FREQ[qno])
 # print(phh_structure.FREQ[ni])
#print(sstructure.e)
#exit()


elphh_structure=elph_structure.elph_structure(phh_structure,'lambda')

'''
for q in range(1,len(ph_structure.Q)+1):
 print('calculations for '+str(q)+'. of total '+str(len(ph_structure.Q))+' q points')
 elph_structure.make_kpoints_single_q(q,structure,ph_structure.Q[q-1])
 elph_structure.read_elph_single_q(q,ph_structure,el_structure,structure)
 elph_structure.elph_single_q_in_whole_kgrid(q,structure,\
                ph_structure,el_structure,'lambda')  #'lambda' or 'elph'
elph_structure.sum_over_q(ph_structure,structure,el_structure)
'''
elphh_structure.parallel_job(sstructure,el_structure,phh_structure)
#elphh_structure.single_job([26,sstructure,phh_structure,el_structure,'lambda'])
elphh_structure.sum_over_q(phh_structure,sstructure,el_structure)
#elphh_structure.extrapolate_values(sstructure,el_structure)

'''
h=open('lambda.frmsf')
tmp=h.readlines()[6+len(el_structure.ENE_fs)*len(structure.allk):]
h.close()
m=0
elph_structure.SUMMED_COLORS=[[ 0 for i in range(len(structure.allk))] for j in range(len(el_structure.ENE_fs))]

for bnd in range(len(el_structure.ENE_fs)):
 for k in range(len(structure.allk)):
    elph_structure.SUMMED_COLORS[bnd][k]=float(tmp[m])
    m+=1
print(m)
elph_structure.extrapolate_values(structure,el_structure)
'''
###


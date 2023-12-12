import xml.etree.ElementTree as ET
import numpy as np
import structure,ph_structure,el_structure
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interpn
from multiprocessing import Process,Pool
import copy
PRECIS=6
ELECTRONMASS_SI  = 9.10938215e-31   # Kg
AMU_SI           = 1.660538782e-27  #Kg
AMU_AU           = AMU_SI / ELECTRONMASS_SI
AMU_RY           = AMU_AU / 2. #=911.44
RY_TO_THZ=3289.8449
RY_TO_GHZ=RY_TO_THZ*1000
RY_TO_CM_1=9.1132804961430295837e-6

class elph_structure_single_q():
 def __init__(self,phh_structure,q):
  self.ELPH_sum=[]
  self.elph_dir=phh_structure.elph_dir
  self.elph_dir2=phh_structure.elph_dir2
  self.prefix=phh_structure.prefix
  self.no_of_s_q=phh_structure.no_of_s_q[q-1]
  self.KPOINTS=[] #[0:2] list of noneq k at given q [3] its no in list of nonequivalent kpoints at given q  (self.KPOINTS)  and [4] its no in list of nonequivalent kpoints of electronic structure (ell_structure.NONEQ)
  self.WEIGHTS_OF_K=[]
  self.nbnd_el=0
  self.fermi_nbnd_el=0
  self.nkp=0
  self.KPOINTS_all=[] #[0:2] list of allk [3] its no in list of nonequivalent kpoints at given q  (self.KPOINTS)  
  self.KPQ=[] #[0] k+q, [1] its no in list of nonequivalent kpoints at given q  (self.KPOINTS)  
  self.lambda_or_elph='' #'lambda' or 'elph'
  self.COLORS=[]
  self.q_no=q

 def make_kpoints_single_q(self,basic_structure,q,q_cryst,qstar):
  print('make kgrid at given q')
  self.pm=[0,-1,1]
  tree = ET.parse(self.elph_dir+'elph.'+str(self.q_no)+'.1.xml')
  root = tree.getroot()
  self.nkp=int(root.find('PARTIAL_EL_PHON/NUMBER_OF_K').text)
  self.KPOINTS=[ [ round(float(m),PRECIS)\
         for m in root.find(
          'PARTIAL_EL_PHON/K_POINT.'+str(k)+'/COORDINATES_XK'
          ).text.split() ]\
           for k in range(1,self.nkp+1)]
  for numki,ki in enumerate(self.KPOINTS): 
   self.KPOINTS[numki].append(numki)
  structure_new=copy.deepcopy(basic_structure)
  structure_new.NONEQ=list(self.KPOINTS)
  structure_new.calc_noneq_cryst() #results in NONEQ_cryst
  structure_new.no_of_kpoints=basic_structure.no_of_kpoints
  structure_new.check_symm(q_cryst,basic_structure.NONEQ,self.no_of_s_q,self.q_no)
  structure_new.calc_irt()
  self.irt=structure_new.irt
  self.rtau=structure_new.rtau
  self.SYMM_q=structure_new.SYMM_crystal
  structure_new.make_kgrid(self.q_no)
  structure_new.find_newkpoints_in_old_list(basic_structure.allk) #aatach 4-th element to each k:  #it's nonequivalent in basic k-point list. Added to both allk AND NONEQ
  self.KPOINTS_all=structure_new.allk
  self.KPOINTS_all_cryst=structure_new.allk_in_crystal_coordinates
  self.KPOINTS_all_cryst_round2=[ list(np.round(i[:3],PRECIS-2))+i[3:] for i in structure_new.allk_in_crystal_coordinates] #np.round(structure_new.allk_in_crystal_coordinates,PRECIS-2)
  self.WEIGHTS_OF_K=structure_new.WK
  self.KPOINTS=structure_new.NONEQ
  self.KPOINTS_cryst=structure_new.NONEQ_cryst 
#  print(q_no,[ i for i in self.KPOINTS if [4]==None])
#  self.KPOINTS_all_all_q.append(structure_new.allk)

  #print(self.KPOINTS)
  '''
  print('star calculations')
  self.KPQ=[[] for qq in qstar]
  for nqq,qq in enumerate(qstar):
   for k in self.KPOINTS:
    self.KPQ[nqq].append(structure_new.find_k_plus_q(k, self.KPOINTS_all,qq)) #[kpq,no its nonequivalent in kpoints list for this q, no its nonequivalent in basic kpoint list]
  print('star ended')
  '''
  self.KPQ=[]
  for k in  self.KPOINTS_cryst:
   self.KPQ.append(structure_new.find_k_plus_q(k, self.KPOINTS_all,self.KPOINTS_all_cryst_round2,q_cryst))
  h=open('q_'+str(self.q_no)+'.info','w')
  h.write(str(len(self.SYMM_q))+'\n')
  #for i in self.KPQ:
  # for j in i:
  #  h.write(str(j)+' ')
  # h.write('\n')

  multik=np.zeros((len(basic_structure.NONEQ)))
  for k in self.KPQ: multik[k[2]]+=1
  for ni,i in  enumerate( self.KPQ):
    h.write(str(ni)+' '+str(i)+'\n')
  h.close()
  for i in multik: 
   if i==0: raise ValueError(str(self.q_no)+' len of multik = 0')

  '''
  #check if kpq are ok
  structure_newq=structure.structure()
  structure_newq.NONEQ_cryst=[list(m[0]) for nm,m in enumerate(self.KPQ)]
  for ni,i in enumerate(structure_newq.NONEQ_cryst): i.append(ni)
  structure_newq.NONEQ=[ [0,0,0,0] for i in range(len(  structure_newq.NONEQ_cryst))]
  structure_newq.SYMM=list(structure_new.SYMM)
  structure_newq.SYMM_crystal=list(structure_new.SYMM_crystal)
  structure_newq.e=basic_structure.e
  structure_newq.no_of_kpoints=basic_structure.no_of_kpoints
 # structure_newq.calc_noneq_cryst() #results in NONEQ_cryst
 # print("UWAGA!!",q_no,len(structure_newq.NONEQ_cryst))
  structure_newq.make_kgrid(q_no)
  structure_newq.find_newkpoints_in_old_list(basic_structure.allk) #aatach 4-th element to each k:  #it's nonequiv
  self.KPQ_all=structure_newq.allk
  '''
 def read_elph_file_before_qe7(self,q_point_no,phh_structure,ell_structure,structure):
  tree = ET.parse(self.elph_dir+'elph.'+str(q_point_no)+'.1.xml')
  root = tree.getroot()
  self.nbnd_el=int(root.find('PARTIAL_EL_PHON/NUMBER_OF_BANDS').text)
  #choose only bands which cross EF and sum over j in pairs <i,j>
  print(str(q_point_no)+': From all '+str(self.nbnd_el)+' bands detected in elph calc. only bands ',ell_structure.bands_num,' cross EF and will be written in frmsf')
  self.fermi_nbnd_el=len(ell_structure.bands_num)
#  ELPH=[[[ [] for k in range(self.nkp)] for j in range(self.nbnd_el)] for i in range(self.fermi_nbnd_el)] #stores elph[k][ibnd][jbnd][nmode]
  ELPH=np.zeros(shape=(self.fermi_nbnd_el,self.fermi_nbnd_el,\
                self.nkp,phh_structure.no_of_modes), dtype=complex) #stores elph[k][ibnd][jbnd][nmode]
  #read elph  from all files
  imode=0
  for mode in range(1,len(phh_structure.NONDEG[q_point_no-1])+1):
   #elph
   tree = ET.parse(self.elph_dir+'elph.'+str(q_point_no)+'.'+str(mode)+'.xml')
   root = tree.getroot()
   for country in root.iter('PARTIAL_EL_PHON'):
    if int(country.find('NUMBER_OF_BANDS').text)!=self.nbnd_el: print(str(q_point_no)+'Warning. No of bands!= no of bands from el structure')
    for k in range(1,self.nkp+1):
     for town in country.iter('K_POINT.'+str(k)):
      partial_elph=town.find('PARTIAL_ELPH')
      if k==1: 
       npert0=int(partial_elph.get('size'))/self.nbnd_el/self.nbnd_el
       npert=int(npert0)
     #  print( npert, mode,phh_structure.ORDER_OF_IRR[q_point_no-1],phh_structure.FREQ[q_point_no-1])
       if npert/npert0!=1: 
        print(str(q_point_no)+"WARNING: npert is not int, but is equal to"+str(npert0))
      elph_k=[ complex(float(m.replace(',',' ').split()[0]), 
                       float(m.replace(',',' ').split()[1])) 
                for m in partial_elph.text.split('\n') if len(m.split())>0  ]
     # print('len of elpk',len(elph_k))
      for iipert in range(npert):
       nmode=imode+iipert #phh_structure.ORDER_OF_IRR[q_point_no-1][imode+iipert]
#       elph_k_ii=symmetrize((elph_k_ii))
       for numiband,iband in enumerate(ell_structure.bands_num):
        for numjband,jband in enumerate(ell_structure.bands_num):
         ELPH[numiband][numjband][k-1][nmode]=elph_k[iipert*self.nbnd_el*self.nbnd_el+jband*self.nbnd_el+iband]
  #       ELPH[numiband][numjband][k-1][nmode]+=elph_k[jband*self.nbnd_el*npert+iband*npert+iipert] #--wrong
  #       ELPH[numiband][numjband][k-1][nmode]=elph_k[jband*self.nbnd_el*npert+iipert*self.nbnd_el+iband] #--wrong
#          elph_k[iipert*self.nbnd_el*npert+iband*self.nbnd_el+jband] #--rather wrong

#                 CALL xmlw_writetag( "PARTIAL_ELPH", el_ph_mat_rec_col(:,:,ik,np) )

#              CALL iotk_write_dat(iunpun, "PARTIAL_ELPH", &
#                                         el_ph_mat_rec_col(:,:,ik,:))
#		el_ph_mat_rec_col(nbnd,nbnd,nksqtot,npe)s
# ! el_ph_mat(i,j,k,I)= <\psi(k+q) n_i|dV_{SCF}/du^q_{i a}|\psi(k) n_j>
#elphon.f90:              el_ph_mat (ibnd, jbnd, ik, ipert + imode0) = elphmat (ibnd, jbnd, ipert)
#elphon.f90:              el_ph_mat_rec (ibnd, jbnd, ik, ipert ) = elphmat (ibnd, jbnd, ipert)
#elphon.f90:        el_ph_mat_rec_col => el_ph_mat_rec

   imode=imode+npert   
  return ELPH
 

 def read_elph_file_after_qe7(self,q_point_no,phh_structure,ell_structure,structure):
  tree = ET.parse(self.elph_dir+'elph.'+str(q_point_no)+'.1.xml')
  root = tree.getroot()
  self.nbnd_el=int(root.find('PARTIAL_EL_PHON/NUMBER_OF_BANDS').text)
  #choose only bands which cross EF and sum over j in pairs <i,j>
  print(str(q_point_no)+': From all '+str(self.nbnd_el)+' bands detected in elph calc. only bands ',ell_structure.bands_num,' cross EF and will be written in frmsf')
  self.fermi_nbnd_el=len(ell_structure.bands_num)
#  ELPH=[[[ [] for k in range(self.nkp)] for j in range(self.nbnd_el)] for i in range(self.fermi_nbnd_el)] #stores elph[k][ibnd][jbnd][nmode]
  ELPH=np.zeros(shape=(self.fermi_nbnd_el,self.fermi_nbnd_el,\
                self.nkp,phh_structure.no_of_modes), dtype=complex) #stores elph[k][ibnd][jbnd][nmode]
  #read elph  from all files
  imode=0
  for mode in range(1,len(phh_structure.NONDEG[q_point_no-1])+1):
   #elph
   tree = ET.parse(self.elph_dir+'elph.'+str(q_point_no)+'.'+str(mode)+'.xml')
   root = tree.getroot()
   for country in root.iter('PARTIAL_EL_PHON'):
    if int(country.find('NUMBER_OF_BANDS').text)!=self.nbnd_el: print(str(q_point_no)+'Warning. No of bands!= no of bands from el structure')
    for k in range(1,self.nkp+1):
     for town in country.iter('K_POINT.'+str(k)):
      partial_elph=town.findall('PARTIAL_ELPH')
      if k==1: 
       npert=len(partial_elph)
     #  print( npert, mode,phh_structure.ORDER_OF_IRR[q_point_no-1],phh_structure.FREQ[q_point_no-1])
      for iipert in range(npert):
       nmode=imode+iipert #phh_structure.ORDER_OF_IRR[q_point_no-1][imode+iipert]
#       elph_k_ii=symmetrize((elph_k_ii))
       elph_k=[ complex(float(m.replace(',',' ').split()[0]), 
                       float(m.replace(',',' ').split()[1])) 
                for m in partial_elph[iipert].text.split('\n') if len(m.split())>0  ]
     # print('len of elpk',len(elph_k))
       for numiband,iband in enumerate(ell_structure.bands_num):
        for numjband,jband in enumerate(ell_structure.bands_num):
         ELPH[numiband][numjband][k-1][nmode]=elph_k[iband*self.nbnd_el+jband]
   imode=imode+npert   
  return ELPH
     
 def read_elph_single_q(self,phh_structure,ell_structure,structure): 
  try: self.ELPH=self.read_elph_file_before_qe7(self.q_no,phh_structure,ell_structure,structure)
  except: self.ELPH=self.read_elph_file_after_qe7(self.q_no,phh_structure,ell_structure,structure)
 
  
 def sum_elph_over_jband(self,phh_structure,ell_structure,structure): 
  self.ELPH_sum=np.zeros(shape=(self.fermi_nbnd_el,\
           self.nkp,phh_structure.no_of_modes),dtype=complex)
  ELPH2=np.zeros(shape=(phh_structure.no_of_modes,phh_structure.no_of_modes), dtype=complex) #stores elph[k][ibnd][jbnd][nmode]

  for iband in range(self.fermi_nbnd_el):
    for numk, k in enumerate(self.KPOINTS): 
     for jband in range(self.fermi_nbnd_el):
      weight=el_structure.w0gauss(-ell_structure.ENE_fs[jband][self.KPQ[numk][2]]) #*weight0 #*self.WEIGHTS_OF_K[self.KPQ[numk][1]]
      for iipert in range(phh_structure.no_of_modes):
       for jjpert in range(phh_structure.no_of_modes):
        ELPH2[iipert][jjpert]=np.conjugate(self.ELPH[iband][jband][numk][jjpert])*self.ELPH[iband][jband][numk][iipert] 
     
      ELPH2 =ph_structure.symmetrize(phh_structure.nat,np.array(phh_structure.PATT[self.q_no-1]),
ELPH2 , structure.at,structure.e, structure.SYMM_crystal, self.SYMM_q,self.irt,self.rtau,phh_structure.Q_crystal[self.q_no-1] )
      ELPH2=ELPH2*weight 
      for nu in range(phh_structure.no_of_modes):
       for mu in range(phh_structure.no_of_modes):
        for vu in range(phh_structure.no_of_modes): 
          dyn= (phh_structure.DYN[self.q_no-1][0])
          self.ELPH_sum[iband][numk][nu] += (np.conjugate(dyn[mu][nu])*ELPH2[mu][vu]*dyn[vu][nu])    #--indeksy jak w elphon
 


 def calc_lambda_summed_over_modes(self,structure,phh_structure,ell_structure,l_or_gep):
  self.lambda_or_elph=l_or_gep
  self.COLORS=np.zeros(shape=(self.fermi_nbnd_el, self.nkp),dtype=complex)
  self.COLORS2=np.zeros(shape=(phh_structure.no_of_modes),dtype=complex)
  wagi0=[[el_structure.w0gauss(-ell_structure.ENE_fs[jband][k[4]]) for numk,k in enumerate(self.KPOINTS)] for jband in range(self.fermi_nbnd_el) ]
  wagi=[[wagi0[jband][numk]*self.WEIGHTS_OF_K[numk] for numk in range(len(self.KPOINTS))] for jband in range(self.fermi_nbnd_el) ]
  sumwag=np.sum(wagi)
  if ell_structure.soc!=0: sumwag=np.sum(wagi)/2 #devide by 2 because we are suming not averaging over spin
  if self.lambda_or_elph=='elph':
   for numk in range(self.nkp):
    for jband in range(self.fermi_nbnd_el): 
     self.COLORS[jband][numk]= np.sum(self.ELPH_sum[jband][numk]) #sum over modes)
  else:
   for num in range(phh_structure.no_of_modes): #modes
    if phh_structure.FREQ[self.q_no-1][num]<=20*RY_TO_CM_1: continue
    for numk in range(self.nkp):
     for jband in range(self.fermi_nbnd_el): 
       lam= .5*(self.ELPH_sum[jband][numk][num]) /((phh_structure.FREQ[self.q_no-1][num])**2)
       self.COLORS[jband][numk]+=lam  
       self.COLORS2[num]+=lam *wagi[jband][numk]     
  self.COLORS2=self.COLORS2.real/sumwag
  self.COLORS2=np.round([self.COLORS2[m] for m in phh_structure.ORDER_OF_IRR[self.q_no-1]],4)
  print(self.q_no,self.COLORS2,sum(self.COLORS2),'from main program')  
  lambda_sum0=sum([sum([ self.COLORS[jband][numk]*wagi[jband][numk] for numk in range(len(self.COLORS[jband]))]) for jband in range(len(self.COLORS)) ])
  print (self.q_no,': lambda_sum=',lambda_sum0/sumwag,'sumwag=',sumwag) #,'dos=',dos)
  
 def elph_single_q_in_whole_kgrid(self,structure,phh_structure,ell_structure,l_or_gep):
  wagi0=[[el_structure.w0gauss(-ell_structure.ENE_fs[jband][k[4]]) for numk,k in enumerate(self.KPOINTS)] for jband in range(self.fermi_nbnd_el) ]
  wagi=np.array([[ bnd[k[3]] for k in self.KPOINTS_all] for bnd in wagi0])
  self.COLORS=np.array([[ bnd[k[3]] for k in self.KPOINTS_all] for bnd in self.COLORS])
  lam=np.sum(self.COLORS*wagi)/np.sum(wagi)
  print('whole grid gives',lam)
  h=open('lambda'+str(self.q_no)+'.frmsf','w')
  h.write(str(structure.no_of_kpoints[0])+' '+str(structure.no_of_kpoints[1])+' ' +str(structure.no_of_kpoints[2])+'\n')
  h.write('1\n'+str(len(ell_structure.bands_num))+'\n')
  for i in structure.e:
   for j in i:
    h.write(str(j)+' ')
   h.write('\n')
  for bnd in ell_structure.ENE_fs:
   for k in structure.allk:
    h.write(str(bnd[k[3]])+'\n')
  for bnd in self.COLORS:
   for k in bnd:
    h.write(str(np.abs(k))+'\n')
  h.close()


 def read_lambdas_from_file(self):
   h=open(self.elph_dir2+'/elph.inp_lambda.'+str(self.q_no))
   tmp=h.readlines()
   h.close()
   m=0
   self.lambdas_from_file=[]
   for nj,j in enumerate(tmp):
    if 'Gaussian' in j: m+=1
    if m==4:
     for k in tmp[nj+2:]:
      if 'Gaussian' in k: break
      self.lambdas_from_file.append(np.round(float(k.split()[2]),4))
     break



 def calc_lambda_q_as_elphon_does(self,phh_structure,ell_structure,structure): 
  ELPH2=np.zeros(shape=(phh_structure.no_of_modes,phh_structure.no_of_modes), dtype=complex) #stores elph[k][ibnd][jbnd][nmode]
  self.ELPH_sum2=np.zeros(shape=( phh_structure.no_of_modes),dtype=complex)
  sumweight=0

  for iband in range(self.fermi_nbnd_el):
  #  ELPH2=np.zeros(shape=(len(structure.NONEQ),phh_structure.no_of_modes,phh_structure.no_of_modes), dtype=complex) #stores elph[k][ibnd][jbnd][nmode]
    for numk, k in enumerate(self.KPOINTS): 
     weight0=el_structure.w0gauss(-ell_structure.ENE_fs[iband][k[4]]) *self.WEIGHTS_OF_K[numk]
     sumweight+=weight0
     for jband in range(self.fermi_nbnd_el):
      weight=weight0*el_structure.w0gauss(-ell_structure.ENE_fs[jband][self.KPQ[numk][2]]) 
      for iipert in range(phh_structure.no_of_modes):
       for jjpert in range(phh_structure.no_of_modes):
        ELPH2[iipert][jjpert]+=np.conjugate(self.ELPH[iband][jband][numk][jjpert])*self.ELPH[iband][jband][numk][iipert] *weight 
  ELPH2 =ph_structure.symmetrize(phh_structure.nat,np.array(phh_structure.PATT[self.q_no-1]),
ELPH2 , structure.at,structure.e, structure.SYMM_crystal, self.SYMM_q,self.irt,self.rtau,phh_structure.Q_crystal[self.q_no-1] )
  lam=[]
  for nu in range(phh_structure.no_of_modes):
       for mu in range(phh_structure.no_of_modes):
        for vu in range(phh_structure.no_of_modes): 
          dyn= (phh_structure.DYN[self.q_no-1][0])
          self.ELPH_sum2[nu] += (np.conjugate(dyn[mu][nu])*ELPH2[mu][vu]*dyn[vu][nu])    #--indeksy jak w elphon
       if phh_structure.FREQ[self.q_no-1][nu]>20*RY_TO_CM_1:
        lam.append((.5*(self.ELPH_sum2[nu] ) /((phh_structure.FREQ[self.q_no-1][nu])**2)/sumweight).real)
       else: lam.append(0)
    #   ELPH2[numk] =ph_structure.symmetrize(phh_structure.nat,np.array(phh_structure.PATT[q_point_no-1]),
#ELPH2[numk] , structure.at,structure.e, structure.SYMM_crystal, structure.SYMM_crystal,structure.irt,self.rtau,phh_structure.Q_crystal[q_point_no-1] )
  print('\n',self.q_no,(np.round([lam[m] for m in phh_structure.ORDER_OF_IRR[self.q_no-1]] ,4)),sum(lam),'calc as elphon')
  self.read_lambdas_from_file()
  print(self.q_no,self.lambdas_from_file,sum(self.lambdas_from_file),'from elph_dir')

def cond(ind,nk):
 if ind[0]!=-1 and ind[0]!=nk[0] and ind[1]!=-1 and ind[1]!=nk[1] and ind[2]!=-1 and ind[2]!=nk[2]: return True
 else: return False



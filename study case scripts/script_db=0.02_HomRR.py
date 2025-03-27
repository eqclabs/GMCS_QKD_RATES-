import sys
import os
absolute_path = os.path.abspath("..")
absolute_path=absolute_path.replace(os.sep, '/')
sys.path.append(absolute_path)
#folder_name=os.path.basename(os.path.abspath("."))
script_name = os.path.splitext(os.path.basename(__file__))[0]
db=0.02
data_folder=str(script_name)+'_data'
if not os.path.exists(data_folder):
    os.makedirs(data_folder)
#===========================================
import numpy as np
import random
import Rate_Calculation_module as RCm
###############################################
protocol="hom_RR"
double_mod=False
alpha=8
eta=0.8
vel=0.01
ecor=2**(-32)
eh=2**(-32) 
es=2**(-32)
ePE=2**(-32)
####################################################
t=10**(-db/10)
e=0.01
#------------------------------------
if random.random() <0.8:
             N_sup=4.5+1*np.array([random.random()])
else:
             N_sup=5.5+0.6*np.array([random.random()])
#----------------------------------------
N=int(np.floor(10**N_sup))
p_vec=[4]
V_vec=np.linspace(2,120,150)
Vpe_vec=[3]
RatioPower_vec=np.linspace(np.log10(2),2,10)
r_vec=10**(-RatioPower_vec)
pEC_vec=[0.9]
beta_in=0
#=================================================================================  
RATE=RCm.opt_rate(protocol,double_mod,beta_in,alpha,e,vel,eta,t,N,ecor,eh,es,ePE,
                  p_vec,V_vec,Vpe_vec,r_vec,pEC_vec) 
beta_out=RATE[0]
Rate=RATE[1]
Rsynd=RATE[2]
SNR=RATE[3]
Leak=RATE[4]
V_pe=RATE[5]
V=RATE[6]
pEC=RATE[7]
eps=RATE[8]
Memory=RATE[9]
Mem_sp=RATE[10] 
p=RATE[11]
nval=RATE[12]
rval=RATE[13]
###--------------------------------------------
log_file = open(str(data_folder)+"/"+str(script_name)+"_bksz="+str(N)+"_dB="+str(db)+".txt", "w")#start_date.strftime("%d-%b-%Y (%H.%M.%S.%f)") +
log_file.write('{0:25}  {1:25}\n'.format('Input------',str()))
log_file.write('{0:25}  {1:25}\n'.format('N', str(N)))
log_file.write('{0:25}  {1:25}\n'.format('protocol', str(protocol)))
log_file.write('{0:25}  {1:25}\n'.format('double_Mod', str(double_mod)))
log_file.write('{0:25}  {1:25}\n'.format('ePE', str(ePE)))
log_file.write('{0:25}  {1:25}\n'.format('es', str(es)))
log_file.write('{0:25}  {1:25}\n'.format('eh', str(eh)))
log_file.write('{0:25}  {1:25}\n'.format('ecor', str(ecor)))
log_file.write('{0:25}  {1:25}\n'.format('vel', str(vel)))
log_file.write('{0:25}  {1:25}\n'.format('eta', str(eta)))
log_file.write('{0:25}  {1:25}\n'.format('xi', str(e)))
log_file.write('{0:25}  {1:25}\n'.format('db', str(db)))
log_file.write('{0:25}  {1:25}\n'.format('alpha', str(alpha)))
log_file.write('{0:25}  {1:25}\n'.format('p_vec', str(p_vec)))
log_file.write('{0:25}  {1:25}\n'.format('V_vec', str(V_vec)))
log_file.write('{0:25}  {1:25}\n'.format('Vpe_vec', str(Vpe_vec)))
log_file.write('{0:25}  {1:25}\n'.format('r_vec', str(r_vec)))
log_file.write('{0:25}  {1:25}\n'.format('pEC_vec', str(pEC_vec)))
log_file.write('{0:25}  {1:25}\n'.format('beta_in', str(beta_in)))
log_file.write('{0:25}  {1:25}\n'.format('Output------------------------',str()))
log_file.write('{0:25}  {1:25}\n'.format('Rate', str(Rate)))
log_file.write('{0:25}  {1:25}\n'.format('Rsynd', str(Rsynd)))
log_file.write('{0:25}  {1:25}\n'.format('SNR', str(SNR)))
log_file.write('{0:25}  {1:25}\n'.format('leakage', str(Leak)))
log_file.write('{0:25}  {1:25}\n'.format('V_pe', str(V_pe)))
log_file.write('{0:25}  {1:25}\n'.format('V', str(V)))
log_file.write('{0:25}  {1:25}\n'.format('epsilon', str(eps)))
log_file.write('{0:25}  {1:25}\n'.format('Mem', str(Memory)))
log_file.write('{0:25}  {1:25}\n'.format('Mem_sp', str(Mem_sp)))
log_file.write('{0:25}  {1:25}\n'.format('p', str(p)))
log_file.write('{0:25}  {1:25}\n'.format('beta_out', str(beta_out)))
log_file.write('{0:25}  {1:25}\n'.format('pEC', str(pEC)))
log_file.write('{0:25}  {1:25}\n'.format('nval', str(nval)))
log_file.write('{0:25}  {1:25}\n'.format('rval', str(rval)))

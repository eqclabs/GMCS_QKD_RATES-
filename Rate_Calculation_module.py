import numpy as np
import Leakage_Calculation_module as LCm
import random
###############
def h(x):
    if x>1:
       x0=(x+1)/2
       x1=(x-1)/2
       output= x0*np.log2(x0)-x1*np.log2(x1)
    else:
        output=0
    return output
# --------------------------------------
def rate_leak_bound(protocol,double_mod,beta,alpha,V,V_pe,e,vel,eta,t,p,pe_ratio,N,pEC,ecor,eh,es,ePE):
    if protocol=='het_RR':
                             D=2
    elif protocol=='hom_DR':
                             D=1
    elif protocol=='hom_RR':
                             D=1
    elif protocol=='het_DR':
                             D=2
    #---------------------------------------------
    N=int(np.floor(N))
    M=int(np.floor(pe_ratio*N))
    #----------------------
    w=np.sqrt(2*np.log(ePE**(-1)))
    #--------------------
    XI=eta*t*e 
    scaling_t=(4*t**2/(M*D))
    scaling_XI=(2/M*D)
    if double_mod:
       sigma_t_one_sq=scaling_t*(XI+vel+D)/(eta*t*(V+V_pe))
       sigma_t_two_sq=scaling_t*(eta*t*V+XI+vel+D)/(eta*t*(V_pe))
       sigma_t_sq=1/((1/sigma_t_one_sq)+(1/sigma_t_two_sq))  
       # -------------------------------------------------
       sigma_XI_one_sq=scaling_XI*(XI+vel+D)
       sigma_XI_two_sq=scaling_XI*(eta*t*V+XI+vel+D)
       sigma_XI_sq=1/((1/sigma_XI_one_sq)+(1/sigma_XI_two_sq))
    else:
       sigma_t_sq=scaling_t*(XI+vel+D)/(eta*t*V)
       sigma_XI_sq=scaling_XI*(XI+vel+D)                         
    twc=t-w*np.sqrt(sigma_t_sq)
    XIwc=XI+w*np.sqrt(sigma_XI_sq)
    ewc=XIwc/(eta*twc)
    #------------------------------
    if twc<0 or twc>t or e>ewc or ewc<0:
         R=0
         SNR=0
         rho=0,
         I=0
         leakEC=0
    else:
         mu=V+1
         if twc==1:
                 omega=1
         else:
                 omega = (twc * ewc - twc + 1) / (1 - twc) 
         #------------------------------
         psi = np.sqrt(twc * (omega ** 2 - 1))
         phi = twc * omega + (1 - twc) * mu   
         #----------------------------------------------------
         v_eE = np.array([[omega, 0,     psi,       0],
                     [0,     omega, 0,      -psi],
                     [psi,   0,     phi,       0],
                     [0,    -psi,   0,       phi]])
         Delta_eE=np.linalg.det(v_eE[0:2,0:2])+np.linalg.det(v_eE[2:4,2:4])+2*np.linalg.det(v_eE[0:2,2:4])
         v1=np.sqrt((1/2)*(Delta_eE+np.sqrt(Delta_eE**2-4*np.linalg.det(v_eE))))
         v2=np.sqrt((1/2)*(Delta_eE-np.sqrt(Delta_eE**2-4*np.linalg.det(v_eE))))
         
         if protocol=='het_RR':
             b = twc * eta * (mu + ewc) + 1 - (twc * eta) + vel
             gamma = np.sqrt(eta * (1 - twc) * (omega ** 2 - 1))
             theta = np.sqrt(eta * twc * (1 - twc)) * (omega - mu)
             v_CeE = v_eE - ((b + 1) ** -1) * np.array([[gamma ** 2,   0,             gamma * theta, 0],
                                                  [0,             gamma ** 2,    0,             -gamma * theta],
                                                  [gamma * theta, 0,             theta ** 2,    0],
                                                  [0,             -gamma * theta, 0,             theta ** 2]])
         elif protocol=='hom_DR':
             phi_0=twc * omega + (1 - twc) 
             v_CeE = np.array([[omega, 0,     psi,       0],
                              [0,     omega, 0,      -psi],
                              [psi,   0,     phi_0,       0],
                              [0,    -psi,   0,       phi]])
             
         elif protocol=='hom_RR':
              b = twc * eta * (mu + ewc) + 1 - (twc * eta) + vel
              gamma = np.sqrt(eta * (1 - twc) * (omega ** 2 - 1))
              theta = np.sqrt(eta * twc * (1 - twc)) * (omega - mu)
              v_CeE = v_eE - ((b) ** -1) * np.array([[gamma ** 2,    0,             gamma * theta, 0],
                                                     [0,             0,              0,            0],
                                                     [gamma * theta, 0,             theta ** 2,    0],
                                                     [0,             0,              0,            0]])
         elif protocol=='het_DR':
                phi_0=twc * omega + (1 - twc) 
                v_CeE = np.array([[omega, 0,     psi,       0],
                                 [0,     omega, 0,      -psi],
                                 [psi,   0,     phi_0,       0],
                                 [0,    -psi,   0,       phi_0]])
             
         Delta_CeE=np.linalg.det(v_CeE[0:2,0:2])+np.linalg.det(v_CeE[2:4,2:4])+2*np.linalg.det(v_CeE[0:2,2:4])
         v3=np.sqrt((1/2)*(Delta_CeE+np.sqrt(Delta_CeE**2-4*np.linalg.det(v_CeE))))
         v4=np.sqrt((1/2)*(Delta_CeE-np.sqrt(Delta_CeE**2-4*np.linalg.det(v_CeE))))
         chi=h(v1)+h(v2)-h(v3)-h(v4)
         #---------------------------------------------
         dim=2**(p*D)
         Delta=4*np.log2(np.sqrt(dim)+2)*np.sqrt(np.log2(2/es**2))
         THETA=np.log2(2*ecor*eh**2)
         SNR=V/(e+(D+vel)/(eta*t))
         rho=np.sqrt(SNR/(1+SNR))
         I=(D/2)*np.log2(1+SNR)
         # if  protocol=='het_RR':
                          #Delta_ent=
         #elif  protocol=='hom_RR':
                          #Delta_ent=
         # elif  protocol=='het_DR':
                          #Delta_ent=0
         #elif  protocol=='hom_DR':
                         # Delta_ent=0
            
         if beta==0:
             leakEC=LCm.leak(D*(N-M),pEC*(1-ecor),rho,alpha,p)
             R=D*LCm.Hk(alpha,p)-chi-leakEC/(N-M)-Delta/np.sqrt(N-M)+(THETA)/(N-M)
             beta_out=(D*LCm.Hk(alpha,p)-leakEC/(N-M))/I
         else:
              R=beta*I-chi-Delta/np.sqrt(N-M)+(THETA)/(N-M)
              leakEC=(N-M)*D*LCm.Hk(alpha,p)-beta*I
              beta_out=beta
    return [pEC*((N-M)/N)*R,SNR,leakEC,int(N-M),beta_out]


def opt_rate(protocol,double_mod,beta_in,alpha,e,vel,eta,t,N,ecor,eh,es,ePE,p_vec,V_vec,Vpe_vec,r,pEC_vec): 
    if protocol=='het_RR':
                             D=2
    elif protocol=='hom_DR':
                             D=1
    elif protocol=='hom_RR':
                             D=1
    elif protocol=='het_DR':
                             D=2
    #---------------------------------------------
    # ===============================
    Rate_val=np.zeros((len(V_vec),len(r),len(p_vec),len(pEC_vec),len(Vpe_vec)))
    SNR_val=np.zeros((len(V_vec),len(r),len(p_vec),len(pEC_vec),len(Vpe_vec)))
    leak_val=np.zeros((len(V_vec),len(r),len(p_vec),len(pEC_vec),len(Vpe_vec)))
    n_val=np.zeros((len(V_vec),len(r),len(p_vec),len(pEC_vec),len(Vpe_vec)))
    beta_out_val=np.zeros((len(V_vec),len(r),len(p_vec),len(pEC_vec),len(Vpe_vec)))
    for i in range(len(V_vec)):
           for j in range(len(r)):
              for k in range(len(p_vec)):
                  for l in range(len(pEC_vec)):
                       for m in range(len(Vpe_vec)):
                              Rate=rate_leak_bound(protocol,double_mod,beta_in,alpha,V_vec[i],Vpe_vec[m],e,vel,eta,t,p_vec[k],r[j],N,pEC_vec[l],ecor,eh,es,ePE)
                              Rate_val[i,j,k,l,m]=Rate[0]
                              SNR_val[i,j,k,l,m]=Rate[1]
                              leak_val[i,j,k,l,m]=Rate[2]
                              n_val[i,j,k,l,m]=Rate[3]
                              beta_out_val[i,j,k,l,m]=Rate[4]
                              #print(Rate[0])
    mi,mj,mk,ml,mm=np.unravel_index(np.argmax(Rate_val,axis=None),Rate_val.shape)
    R_synd=leak_val[mi,mj,mk,ml,mm]/(D*p_vec[mk]*n_val[mi,mj,mk,ml,mm])
    eps_tot=es+eh+2*ePE+ecor
    L=int(np.ceil(leak_val[mi,mj,mk,ml,mm]))
     # ------------
    Mem=D*n_val[mi,mj,mk,ml,mm]*L/(8*10**6)
    # -------------
    dc=2
    Mem_sparse=(dc*D*n_val[mi,mj,mk,ml,mm]*(p_vec[mk]+np.ceil(np.log2(n_val[mi,mj,mk,ml,mm])))+((L/p_vec[mk])+1)*np.ceil(np.log2(dc*D*n_val[mi,mj,mk,ml,mm])))/(8*10**6)
    return beta_out_val[mi,mj,mk,ml,mm], Rate_val[mi,mj,mk,ml,mm], R_synd, SNR_val[mi,mj,mk,ml,mm], L/(8*10**6), Vpe_vec[mm],V_vec[mi], pEC_vec[ml], eps_tot, Mem, Mem_sparse, p_vec[mk],n_val[mi,mj,mk,ml,mm],r[mj]

# x=2*10**5*4*2
# y=int(np.ceil(2*10**5*np.ceil(np.log2(10**5))))*2
# z=int(np.ceil((2*0.6667*10**5+1))*np.ceil(np.log2(2*0.6667*10**5+1)))
# print(x, y, z,x+y+z,(x+y+z)/(8*10**6))
# res=(x+y+z)/(8*10**6)
# print(res)


# Notes:-----------------------
# if double_mod:
#     sigma_t_one=np.sqrt(2)*t*M**(-1/2)*np.sqrt(((XI+2+vel)/(eta*t*(V_pe+V))))
#     sigma_t_two=np.sqrt(2)*t*(N-M)**(-1/2)*np.sqrt(((t*V_pe+XI+2+vel)/(eta*t*(V))))
#     sigma_t=1/((1/sigma_t_one)+(sigma_t_two))
#     twc=t-w*sigma_t
#     #------------------------------
#     sigma_XI_one=np.sqrt(((XI+vel+2)**2/M))
#     sigma_XI_two=np.sqrt(((t*V_pe+XI+vel+2)**2/(N-M)))
#     sigma_XI=1/((1/sigma_XI_one)+(sigma_XI_two))
# else:
#     sigma_t=np.sqrt(2)*t*M**(-1/2)*np.sqrt(((XI+2+vel)/(eta*t*(V))))
#     sigma_XI=np.sqrt(((XI+vel+2)**2/M))


# def opt_rate_simple(pEC,protocol,double_mod,beta,alpha,e,vel,eta,t,p,N,ecor,eh,es,ePE): 
#     N=int(np.floor(N))
#     V=np.linspace(2,100,150)#150
#     if double_mod:
#         V_pe=[3]#np.linspace(0,5,1)#10
#     else:
#         V_pe=[0]
#     ratio_power_pe=np.linspace(np.log10(2),2,10)#10
#     r=10**(-ratio_power_pe)
#     # if beta==0:
#     #    pEC=0.9#np.linspace(0.8,0.99,10)#10
#     # elif beta==0.95:
#     #    pEC=[0.9]
#     # elif beta==1:
#     #     pEC==1
#     #-------------------------------
#     Rate_val=np.zeros((len(V),len(r),len(p),len(pEC),len(V_pe)))
#     Vval=np.zeros((len(V),len(r),len(p),len(pEC),len(V_pe)))
#     SNR_val=np.zeros((len(V),len(r),len(p),len(pEC),len(V_pe)))
#     leak_val=np.zeros((len(V),len(r),len(p),len(pEC),len(V_pe)))
#     n_val=np.zeros((len(V),len(r),len(p),len(pEC),len(V_pe)))
#     V_pe_val=np.zeros((len(V),len(r),len(p),len(pEC),len(V_pe)))
#     pvalue=np.zeros((len(V),len(r),len(p),len(pEC),len(V_pe)))
#     for i in range(len(V)):
#            for j in range(len(r)):
#               for k in range(len(p)):
#                   for l in range(len(pEC)):#(beta,alpha,V,V_ratio,e,eta,t,vel,p,pe_ratio,N,pEC,ecor,eh,es,ePE)
#                        for m in range(len(V_pe)):
#                               Rate=rate_het_leak_bound(double_mod,beta,alpha,V[i],V_pe[m],e,vel,eta,t,p[k],r[j],N,pEC[l],ecor,eh,es,ePE)
#                               Rate_val[i,j,k,l,m]=Rate[0]
#                               SNR_val[i,j,k,l,m]=Rate[1]
#                               leak_val[i,j,k,l,m]=Rate[3]
#                               Vval[i,j,k,l,m]=Rate[2]
#                               n_val[i,j,k,l,m]=Rate[4]
#                               V_pe_val[i,j,k,l,m]=Rate[5]
#                               pvalue[i,j,k,l,m]=Rate[6]
#                               print(Rate[0])
#     mi,mj,mk,ml,mm=np.unravel_index(np.argmax(Rate_val,axis=None),Rate_val.shape)
#     R_synd=leak_val[mi,mj,mk,ml,mm]/(2*p[mk])
#     eps_tot=es+eh+2*ePE+ecor
#     Mem=2*n_val[mi,mj,mk,ml,mm]*int(np.ceil(leak_val[mi,mj,mk,ml,mm]))/(8*10**6)
#     wc=2
#     Mem_sparse=2*n_val[mi,mj,mk,ml,mm]*wc*(p[mk]+32)/(8*10**6)
#     return Rate_val[mi,mj,mk,ml,mm], R_synd, SNR_val[mi,mj,mk,ml,mm], leak_val[mi,mj,mk,ml,mm], V_pe_val[mi,mj,mk,ml,mm],Vval[mi,mj,mk,ml,mm], pEC[ml], eps_tot, Mem, Mem_sparse,pvalue[mi,mj,mk,ml,mm]




















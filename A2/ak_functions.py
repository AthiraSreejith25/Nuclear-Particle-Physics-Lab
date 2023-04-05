import numpy as np
import matplotlib.pyplot as pl
from scipy.stats import moyal



class Material:
    def __init__(self,X0,X1,a,C,rho,aA,aZ,m,I):
        
        self.X0 =  X0
        self.X1 = X1
        self.a = a
        self.C = C
        self.rho = rho
        self.aA = aA
        self.aZ = aZ
        self.m = m
        self.I = I



def density_correction(beta, gamma, material):
    X = np.log10(beta*gamma)
    
    if(X < material.X0):
        delta = 0
    elif(X < material.X1):
        delta = 4.6052*X + material.a*(material.X1 - X)**(material.m) + material.C
    elif(material.X1 < X):
        delta = 4.6052*X + material.C
        
    return delta
        
def ionization_loss(E, material):
    m_e = 0.511 #MeV/c^2
    m_u = 105.7 #MeV/c^2

    gamma = E/m_u + 1 
    beta = np.sqrt(1 - 1/gamma**2)
    
    K = 0.307075 #MeV g^-1 cm^2
    
    delta = density_correction(beta, gamma, material)
    
    Em = 2*m_e*(beta*gamma)**2/(1 + 2*gamma*m_e/m_u + (m_e/m_u)**2)
    
    coeff = K * material.aZ/material.aA * 1/(beta**2)
    
    dEdx = coeff*(0.5*np.log((2*m_e*((beta*gamma)**2)*Em)/(material.I**2)) + Em/(8*E) - beta**2 - 0.5*delta)
    
    return material.rho * dEdx

def energy_minimum(E, material, X):
    a, b = ionization_loss(E, material), material.b # MeV/cm, 1/cm
    
    return (a/b)*(np.exp(b*X) - 1)



def loss(E_list,f):
    
        medium = copper
        n = len(E_list)
        #print(E_list)
        #Calculating E_mu
        a_list = []
        E_min_list = []
        for j in E_list:
            a = ionization_loss(j,medium)
            a_list.append(a)

            E_min_list.append(energy_minimum)
        if len(E_list) > 1:
                #E_distri = np.histogram(E_list,n)[0] #Check the bin size
                #pl.plot(E_distri, E_list)
            pl.hist(E_list,n)
            pl.savefig("E_distribution_{}{}muons.png".format(n,f))
            pl.show()

            #print(a_list,E_list)
            a_list.sort()
            pl.plot(a_list,E_list,'o')
            pl.savefig("E_vs_a_for_{}{}muons.png".format(n,f))
            pl.show()

            pl.plot(E_list, [b*E_j for E_j in E_list])
            pl.savefig("E_vs_bE_for_{}{}muons.png".format(n,f))
            pl.show()

            #Energy loss
            dl = 100#Detection length in cm



            #theoretical delta E
            #deltaE_g_list = [0 for l in range(0,10000,10)]
            #deltaE_l_list = [0 for l in range(0,10000,10)]
            #deltaE_t_list = []
            deltaE_t_list = [0 for l in range(0,10000,10)]
            deltaE_g_list = [0 for l in range(0,10000,10)]
            deltaE_l_list = [0 for l in range(0,10000,10)]
            #for i in range(0,10000,10):
                #deltaE_t_list.append(0)
            print(len(deltaE_t_list))

            deltaE_g_slice1 = [[],[]]
            deltaE_l_slice1 = [[],[]]

            if n == 1:
                Ei_t_list = []
                Ei_g_list = []
                Ei_l_list = []
            
            for E_k in range(len(E_list)): #iterating through muons

                E_i_g = E_i_l = E_i_t = E_list[E_k]
                                            
                for i in range(0,10000,10):#iterating through slices
                    i = int(i/10)
                    print("slice",i)
                        

                    if n == 1:
                        Ei_t_list.append(E_i_t)
                        Ei_g_list.append(E_i_g)
                        Ei_l_list.append(E_i_l)

                    #Theoretical
                    if E_i_t > m_mu:
                        #print("yay")
                        print(i)
                        deltaE_t = ionization_loss(E_i_t,medium) + b*E_i_t
                        deltaE_t_list[i] = (deltaE_t_list[i] + deltaE_t)/2 #Adding mean of delta E theoretical
                        E_i_t -= deltaE_t*10
                        

                    #GAussian
                    if E_i_g > m_mu: #muon rest mass energy----------------------------------------
                    
                        deltaE_g = ionization_loss(E_i_g,medium) + b*E_i_g

                        #Experimental delta E
                        deltaE_exp_g = np.random.normal(deltaE_g, 0.1, size=None)
                        deltaE_g_list[i] = (deltaE_g_list[i] + deltaE_exp_g)/2 #Adding mean of delta E gaussain expt for each slice

                        #Update E_i_g
                        E_i_g -= deltaE_exp_g*10


                    #Landau
                    if E_i_l > m_mu: #muon rest mass energy----------------------------------------

                        deltaE = ionization_loss(E_i_l,medium) + b*E_i_l
                        #experimental
                        deltaE_exp_l = moyal.rvs(E_i_l,0.1)
                        deltaE_l_list[i] = (deltaE_l_list[i] + deltaE_exp_l)/2 #Adding mean of delta E landau expt for each slice

                        #upadte E_i_l
                        E_i_l -= deltaE_exp_l*10

                    
                    if n > 1:
                        if i in [1,2]:
                            deltaE_g_slice1[i-1].append(deltaE_exp_g)
                            deltaE_l_slice1[i-1].append(deltaE_exp_l)
                            
                        
                        
            
            #plots
            #print(deltaE_t_list[:10])
            #print(deltaE_g_list[:10])
            #print(deltaE_l_list[:10])
            
            pl.plot(range(0,10000,10), deltaE_t_list,'r')
            pl.plot(range(0,10000,10), deltaE_g_list,'g')
            pl.plot(range(0,10000,10), deltaE_l_list,'b')
            pl.savefig("deltaE_vs_slice_{}{}.png".format(n,f))
            pl.show()

            if n == 1:
                pl.plot(range(0,10000,10), Ei_t_list,'r')
                pl.plot(range(0,10000,10), Ei_t_list,'g')
                pl.plot(range(0,10000,10), Ei_t_list,'b')
                pl.savefig("Ei_vs_slice.png")
                pl.show()

            else:
                pl.hist(deltaE_g_slice1[0],100)
                pl.hist(deltaE_g_slice1[1],100)
                pl.hist(deltaE_l_slice1[0],100)
                pl.hist(deltaE_l_slice1[1],100)
                pl.savefig("E_distri_{}.png".format(f))
                pl.show()




b = 4e-6
E_max = 90e9
E_min = 100e6

c = 3e10

m_e = .511 #MeV
m_mu = 105.7 #MeV

#rho = 1.11 for now
copper = Material(0.20,3,0.119,-4.74,8.96,63.55,29,3.38,27.7)
#Copper
z = 29 #atomic number of medium
i_m = 27.7#i of material
muons = 100#number of muons
"""
#Single muon
E_l = [np.random.uniform(E_max,E_min)]
loss(E_l,1)

#Million muons, fixed E

E_2 = []
a = np.random.uniform(E_max,E_min)
for i in range(muons):
    E_2.append(a)
loss(E_2,2)
"""

#Million muons, different E
E_3 = np.random.uniform(E_max,E_min,muons)
loss(E_3,3)


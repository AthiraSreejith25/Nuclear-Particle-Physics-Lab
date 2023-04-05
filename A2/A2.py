import numpy as np
import matplotlib.pyplot as pl
from scipy.stats import moyal


b = 4e-6
E_max = 90e9
E_min = 100e6
#E = 1000000
   
c = 3e10

m_e = .511 #MeV
m_mu = 105.7 #MeV

X_0 = 0.20
X_1 = 3
a = 0.119
m = 3.38
C = -4.74

def Deltaa(beta,gamma):
        
        X = np.log10(beta*gamma)

        if X < X_0:

                delta = 0

        elif X < X_1:

                delta = 4.6052*X + (a)(X_1+X)*(m) + C

        else:

                delta = 4.6052*X + C

        return(delta)

def bethe(E):
    p = ((E*2-(m_mu**2))**.5)/c

    beta = p/E
    gamma = E/m_mu

    delta = Deltaa(beta,gamma)

    t_max = (2*m_e*(c*beta*gamma)**2)/(1 + ((2*gamma*m_e)/(m_mu + (m_e/m_mu)**2)))
                             
    
    return(0.307075*z*(1/beta**2)*(0.5*np.log10((2*m_e*(c*beta*gamma)**(2)*t_max)/i_m**2) - beta**2 - delta/2))

def loss(E_list,f):
    
        
        n = len(E_list)
        #Calculating E_mu
        a_list = []
        E_min_list = []
        for j in E_list:
            a = bethe(j)
            a_list.append(a)
            E_min_list.append((np.e**(b*1e4) - 1)*a/b)
            if len(E_list) > 1:
                #E_distri = np.histogram(E_list,n)[0] #Check the bin size
                #pl.plot(E_distri, E_list)
                pl.hist(E_list,n)
                pl.savefig("E_distribution_{}muons.png".format(n))
                pl.show()

                pl.plot(E_list,a_list)
                pl.savefig("E_vs_a_for_{}muons.png".format(n))
                pl.show()

                pl.plot(E_list, b*E_list)
                pl.savefig("E_vs_bE_for_{}muons.png".format(n))
                pl.show()

            #Energy loss
            dl = 100#Detection length in cm

            #theoretical delta E
            #deltaE_g_list = [0 for l in range(0,10000,10)]
            #deltaE_l_list = [0 for l in range(0,10000,10)]
            deltaE_t_list = [0 for l in range(0,10000,10)]
            deltaE_g_list = [0 for l in range(0,10000,10)]
            deltaE_l_list = [0 for l in range(0,10000,10)]

            deltaE_g_slice1 = [[],[]]
            deltaE_l_slice1 = [[],[]]

            if n == 1:
                Ei_t_list = []
                Ei_g_list = []
                Ei_l_list = []

            for E_k in range(len(E_list)): #iterating through muons

                E_i_g = E_i_l = E_i_t = E_list[E_k]
                
                
                            
                for i in range(0,10000,10):#iterating through slices
                    #print("slice",i)
                        

                    if n == 1:
                        Ei_t_list.append(E_i_t)
                        Ei_g_list.append(E_i_g)
                        Ei_l_list.append(E_i_l)

                    #Theoretical
                    if E_i_t > m_mu:
                        print("yay")

                        deltaE_t = bethe(E_i_t) + b*E_i_t
                        deltaE_t_list[i] = (deltaE_t_list[i] + deltaE_t)/2 #Adding mean of delta E theoretical
                        E_i_t -= deltaE_t*10
                        

                    #GAussian
                    if E_i_g > m_mu: #muon rest mass energy----------------------------------------
                    
                        deltaE_g = bethe(E_i_g) + b*E_i_g

                        #Experimental delta E
                        deltaE_exp_g = np.random.normal(deltaE_g, 0.1, size=None)
                        deltaE_g_list[i] = (deltaE_g_list[i] + deltaE_exp_g)/2 #Adding mean of delta E gaussain expt for each slice

                        #Update E_i_g
                        E_i_g -= deltaE_exp_g*10


                    #Landau
                    if E_i_l > m_mu: #muon rest mass energy----------------------------------------

                        deltaE = bethe(E_i_l) + b*E_i_l
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



#Copper
z = 29 #atomic number of medium
i_m = 27.7#i of material
#Single muon
E_l = [np.random.uniform(E_max,E_min)]
loss(E_l,1)

#Million muons, fixed E
E_2 = []
a = [np.random.uniform(E_max,E_min)]
for i in range(1000000):
    E_2.append(a)
loss(E_2,2)

#Million muons, different E
E_3 = [np.random.uniform(E_max,E_min) for i in range(10e6)]
loss(E_3,3)
















"""
                    

z - atomic number of medium

p = ((E*2-m2*c4)*.5)/c
beta = p/E
m_e = .511 #MeV
m_mu = 105.7 #MeV
gamma = gamma = E/m_mu
t_max = 2*m_e*(c*beta*gamma)**2/(1 + (2*gamma*m_e)/(m_mu + (m_e/m_mu)**2)
i_m = #i of material
delta =                                 



bethe = 0.307075*z*(1/beta**2)*(0.5*np.log10((2*m_e*(c*beta*gamma)**(2)*t_max)/i_m**2) - beta**2 - delta/2)
    
                    


#finding E_List for each case
for i in range(n):
    E_list = [random.uniform(E_max,E_min) for j in range(E_n)]                    
                    
"""                
                
        
        
    
    

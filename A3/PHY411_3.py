import numpy as np
import matplotlib.pyplot as pl

class Material:
    def __init__(self, I, A, B, C, a, m, X1, X0, aZ, aA, rho):
        
        self.I = I * 13.6056e-6 #parameter for ionization loss

        #parameters for density correction
        self.A, self.B, self.C = A, B, C
        self.a, self.m = a, m 
        self.X0, self.X1 = X0, X1

        
        self.aZ, self.aA = aZ, aA #atomic number, atomic mass

        
        self.rho = rho #material density in gm/cc
        self.b = 4e-6 * rho #MeV/cm


def density_correction(beta, gamma, material):
    X = np.log10(beta*gamma)
    delta = 0
    if(X < material.X0):
        delta = 0
    elif(X < material.X1):
        delta = 4.6052*X + material.a*(material.X1 - X)**(material.m) + material.C
    elif(material.X1 < X):
        delta = 4.6052*X + material.C
    return delta


         
def ionisation_loss(E,mass, material):
    m_e = 0.511 #MeV/c^2
    m_u = 105.7 #MeV/c^2

    gamma = E/mass + 1

    beta = np.sqrt(1 - (1/gamma**2))

    K = 0.307075 #MeV g^-1 cm^2
    
    delta = density_correction(beta, gamma, material)
    Em = 2*m_e*(beta*gamma)**2/(1 + 2*gamma*m_e/m_u + (m_e/m_u)**2)
    coeff = K * material.aZ/material.aA * 1/(beta**2)
    dEdx = coeff*(0.5*np.log((2*m_e*((beta*gamma)**2)*Em)/(material.I**2)) + Em/(8*E) - beta**2 - 0.5*delta)
    
    return material.rho * dEdx


Copper = Material(27.7, 0.0701, 15.09, -4.74, -0.119, 3.38, 3, 0.20, 29, 63.5, 8.92)  

np.random.seed(1)

particles = ['Electron', 'Muon', 'Tauon', 'Pion', 'Kaon', 'Proton', 'Neutron', 'Hydrogen', 'Deuterium', 'Helium']
mass = [0.511,105.658,1776.86,139.57,493.677,938.272,939.565,938.781,1875.6127,3727.379]
n_particles = 10000 #number of particles
E_max = 1
E_min = 10000

thickness = 100 #in mm

E = []
p_all = []



#Plotting the energy spectrum

for j in range(len(particles)):

    E_i = np.random.uniform(E_max, E_min, n_particles)

    p_i = [((1 + (E_i[k]/mass[j]))**2 - 1)*mass[j] for k in range(len(E_i))]    
    E.append(E_i)
    p_all.append(p_i)
'''
    pl.hist(E_i,int(n_particles/100))
    pl.xlabel("Energy in Mev")
    pl.ylabel("Number of particles")
    pl.title("Energy spectrum of {}".format(j))
    pl.savefig("Energy_Spectrum_{}.png".format(j))
    pl.show()
'''


E_loss = []
E_net_loss_all = []


for j in range(len(particles)):
    E_loss_j = []


    for i in range(thickness):
        E_loss_i = []
        
        for k in range(n_particles):
            
            if E[j][k] != 0:
                E_loss_i.append(ionisation_loss(E[j][k],mass[j],Copper)) #Energy loss for all particles
            else:
                E_loss_i.append(0)
                
        E[j] -= E_loss_i[j]

        E_loss_j.append(E_loss_i) 


    #Calculating net energy loss of each particle
    E_net_loss = []
    for each in range(len(E_loss_j[0])):
        net_loss = 0
        for ay in range(len(E_loss_j)):
            net_loss += E_loss_j[ay][each]
        E_net_loss.append(net_loss) 
            
        
    E_net_loss_all.append(E_net_loss)
    E_loss.append(E_loss_j)

    #momentum vs dE/dx plot for each particle
    pl.plot(p_all[j],E_net_loss,".")
    pl.xlabel("Incident momentum")
    pl.ylabel("dE/dx")
    pl.title("{}".format(particles[j]))
    pl.savefig("plot2_{}.png".format(particles[j]))
    pl.show()









    


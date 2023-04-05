import random
import matplotlib.pyplot as pl

def nuclei_track(n,alpha,dt,start,stop):
    time = [i for i in range(start,stop+1,dt)] #time
    print(len(time))
    n_list = [n] #list of number of nuclei at each time step

    for j in range(start,stop,dt):
        for i in range(n):
            if random.random() < alpha*dt: 
                n -= 1 #Nucleus decays, update n
                
        n_list.append(n) #Update the list of nuclei at each time step

    return(n_list,time) 


n0 = [100,9000,6*10**6] #List of initial number of nuceli
a = [0.01,0.05,0.2] #List of alpha values
delt = [1,1,10*60] #List of values of time-step

for k in range(len(n0)):
    nl,tl = nuclei_track(n0[k],a[k],delt[k],0,100*3600)
    pl.plot(tl,nl,marker = 'o')
    pl.xlabel("Time")
    pl.ylabel("Number of nuclei")
    pl.savefig("Nuclei_vs_time_{}n.png".format(k))
    pl.show()



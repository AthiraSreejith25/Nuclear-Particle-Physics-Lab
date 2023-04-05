import random
import matplotlib.pyplot as pl


#Subproblem 1
n = 20000 #The number of coins to be tossed

mean_heads = [] 
mean_tails = []
toss_number = []

def toss(num,toss_num):
    toss_number.append(num)
    heads = 0 #count of the number of heads
    tails = 0 #count of the number of heads

    for i in range(num): #iterating through each coin
        for j in range(toss_num): #iterating through each toss of a single coin
            if random.random() > 0.5: 
                heads += 1 #Head is obtained
            else:
                tails += 1 #Head is obtained

    mean_heads.append(heads/toss_num) #mean number of heads obtained
    mean_tails.append(tails/toss_num) #mean number of tails obtained

    
    num = int(num/2) #Halve the number of coins to be tossed

    if num%2 != 0:
        num += 1 #add 1 if the number of coins to be tossed is odd
        
    if num > 9:  
        toss(num,toss_num) #repeat the experiment if number of coins to be tossed is not a single digit number

toss(n,1000) 
       
pl.plot(toss_number,mean_heads,marker = 'o')
pl.ylabel("Mean number of heads")
pl.xlabel("Number of coins tossed")
pl.grid(True)
pl.savefig("Heads_vs_toss_number{}.png".format(n))
pl.show()

pl.plot(toss_number,mean_tails,marker = 'o')
pl.ylabel("Mean number of tails")
pl.xlabel("Number of coins tossed")
pl.grid(True)
pl.savefig("Tails_vs_toss_number{}.png".format(n))
pl.show()

pl.plot(mean_heads,mean_tails,marker = 'o')
pl.ylabel("Mean number of tails")
pl.xlabel("Mean number of heads")
pl.grid(True)
pl.savefig("Heads vs tails_{}tosses.png".format(n))
pl.show()




        

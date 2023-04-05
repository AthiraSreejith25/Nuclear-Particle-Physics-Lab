import random
import matplotlib.pyplot as pl

num_heads = []
toss_number = []

def toss_2(num):
    toss_number.append(num)
    heads = 0 #count of number of heads
    tails = 0 #count of number of tails

    for i in range(num): #iterating through each coin
        
        if random.random() > 0.5:
                heads += 1 #Head is obtained, update number of heads

    num_heads.append(heads)
    
    if heads > 9: 
        toss_2(heads) #Repeat the experiment if the number of coins to be tossed is not a single digit number

    else:
        pl.plot(toss_number,num_heads,marker = 'o')
        pl.savefig("Q2_Heads_vs_Toss_number.png")
        pl.show()
                
toss_2(100000)

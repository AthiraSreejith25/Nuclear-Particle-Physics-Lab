from scipy.optimize import curve_fit
from scipy import optimize
from scipy.interpolate import interp1d
import csv
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt



def half_min_indices(hm,l):

        max_index = l.index(max(l))

        diff_1 = [abs(hm-l[i]) for i in range(max_index)]
        diff_2 = [abs(hm-l[i]) for i in range(max_index,len(l))]

        return (diff_1.index(min(diff_1)), diff_2.index(min(diff_2)) + (max_index))

def gaussian(x,a,b,c):
    return a*np.exp(-(x-b)**2/(2*c**2))


maxi = []
fwhm = []

with open('vlab4_15-11-21_data.csv', 'r', encoding = "utf-8") as csv_1:
    countie = list(csv.reader(csv_1))
    for i in range(0,len(countie),2):
        x = [i/100 for i in range(150,255,5)]
        print(x)
        
        y = [int(j)/5 for j in countie[i+1]]
        print(y)

        



            
        params, covar = curve_fit(gaussian, x, y)
        mu, sigma, c = params[0], params[1],params[2]

        x_new = np.linspace(min(x),max(x), num=1000)
        y_new = gaussian(x_new, mu, sigma,c).tolist()

        maximum = max(y_new)
        half_max = maximum/2

        y_p = [half_max, maximum, half_max]
        indices = half_min_indices(half_max,y_new)
        x_p = [x_new[indices[0]],x_new[y_new.index(maximum)],x_new[indices[1]]]

        fwhm.append(x_p[2]-x_p[0]/x_new[y_new.index(maximum)])





        plt.title('Operating Voltage = {}V'.format(575 + (25*((i/2)))))
        plt.xlabel('LLD')
        plt.ylabel('Count')
        plt.plot(x, y, 'go', x_p, y_p, 'bo', x_new, y_new, [x_p[0], x_p[2]], [y_p[0], y_p[2]], 'b--')
        plt.legend(['data','peak and half-peak','best-fit','FWHM'],loc='best')
        plt.savefig('{}V.png'.format(575 + (25*(i/2)+1)))


        plt.show()


#--------resolution_vs_voltage--------------
voltage = [k for k in range(575,750,25)]
cubic_ip = interp1d(voltage,fwhm, kind='cubic')
x_new = np.linspace(min(voltage), max(voltage), num=1000)
y_new = cubic_ip(x_new).tolist()

min_index = y_new.index(min(y_new))

#absolute and relative (= (100*FWHM/peak_energy)%)

plt.xlabel('Operating Voltage')
plt.ylabel('Resolution')
plt.plot(voltage, fwhm, 'mo', x_new, y_new, [x_new[min_index]], [y_new[min_index]], 'bo')
plt.legend(['data','interpolated','minimum'],loc='best')
plt.savefig('resolution.png')
plt.show()
print('The optimal operating voltage: {}V'.format(round(x_new[min_index],2)))


fwhm = [j*100 for j in fwhm]
cubic_ip = interp1d(voltage,fwhm, kind='cubic')
x_new = np.linspace(min(voltage), max(voltage), num=1000)
y_new = cubic_ip(x_new).tolist()

min_index = y_new.index(min(y_new))

plt.xlabel('Operating Voltage')
plt.ylabel('Resolution in %')
plt.plot(voltage, fwhm, 'mo', x_new, y_new, [x_new[min_index]], [y_new[min_index]], 'bo')
plt.legend(['data','interpolated','minimum'],loc='best')
plt.savefig('resolution_percent.png')
plt.show()

print('The optimal operating voltage: {}V'.format(round(x_new[min_index],2)))



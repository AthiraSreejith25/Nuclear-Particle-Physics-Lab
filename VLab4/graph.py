from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
import csv
import os

#try gaussian fit as well (instead of interpolation) for spectra

#---------test------------
def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

voltage = list(range(500,1000,50))

#-----half_min_locator-------
def half_min_indices(hm,l):

	max_index = l.index(max(l))

	diff_1 = [abs(hm-l[i]) for i in range(max_index)]
	diff_2 = [abs(hm-l[i]) for i in range(max_index,len(l))]

	return (diff_1.index(min(diff_1)), diff_2.index(min(diff_2)) + (max_index))


#-----------spectra---------------
fwhm = []


'''
os.system('mkdir plots')

files = os.listdir('csv')
files.sort()

voltage = []

for i in files:

	voltage.append(int(i.replace('.csv','')))

	with open(i,'r') as file:

		reader = csv.reader(file)

		x, y = [round(int(j[0])/100,2) for j in reader], [int(j[1]).strip() for j in reader]

'''
for i in voltage:

	x = np.linspace(0, 10, num=15)
	y = gaussian(x, 3, .1*(abs(i-720))**.5)

	'''
	cubic_ip = interp1d(x, y, kind='cubic')
	x_new = np.linspace(min(x), max(x), num=1000)	#linspace range and num might be subject to change
	y_new = cubic_ip(x_new).tolist()

	maximum = max(y_new)
	half_max = maximum/2

	y_p = [half_max, maximum, half_max]
	indices = half_min_indices(half_max,y_new)
	x_p = [x_new[indices[0]],x_new[y_new.index(maximum)],x_new[indices[1]]]

	fwhm.append(x_p[2]-x_p[0])

	plt.title('Spectrum At Operating Voltage = {}V'.format(str(i)))
	plt.xlabel('Energy')
	plt.ylabel('Counts')
	plt.plot(x, y, 'o', x_p, y_p, 'ro', x_new, y_new, [x_p[0], x_p[2]], [y_p[0], y_p[2]], 'r:')
	plt.legend(['data','peak and half-peak','interpolated','FWHM'],loc='best')
	plt.savefig('plots/{}V.png'.format(str(i)))
	plt.clf()
	#plt.show()
	'''
	params, covar = curve_fit(gaussian, x, y)
	mu, sigma = params[0], params[1]

	x_new = np.linspace(min(x),max(x), num=1000)
	y_new = gaussian(x_new, mu, sigma).tolist()

	maximum = max(y_new)
	half_max = maximum/2

	y_p = [half_max, maximum, half_max]
	indices = half_min_indices(half_max,y_new)
	x_p = [x_new[indices[0]],x_new[y_new.index(maximum)],x_new[indices[1]]]

	fwhm.append(x_p[2]-x_p[0])
	


	plt.title('Spectrum At Operating Voltage = {}V'.format(str(i)))
	plt.xlabel('LLD')
	plt.ylabel('Total Counts')
	plt.plot(x, y, 'o', x_p, y_p, 'ro', x_new, y_new, [x_p[0], x_p[2]], [y_p[0], y_p[2]], 'r:')
	plt.legend(['data','peak and half-peak','best-fit','FWHM'],loc='best')
	plt.savefig('plots/{}V.png'.format(str(i)))
#	plt.clf()
	plt.show()



#----------resolution_vs_voltage--------------
cubic_ip = interp1d(voltage,fwhm, kind='cubic')
x_new = np.linspace(min(voltage), max(voltage), num=1000)
y_new = cubic_ip(x_new).tolist()

min_index = y_new.index(min(y_new))
#absolute and relative (= (100*FWHM/peak_energy)%)
plt.title('Finding the Optimal Operating Voltage')
plt.xlabel('Operating Voltage')
plt.ylabel('FWHM')
plt.plot(voltage, fwhm, 'o', x_new, y_new, [x_new[min_index]], [y_new[min_index]], 'ro')
plt.legend(['data','interpolated','minimum'],loc='best')
plt.savefig('plots/resolution.png')
#plt.show()

print('The optimal operating voltage is {}V'.format(round(y_new[min_index],2)))


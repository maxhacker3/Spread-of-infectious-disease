import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap
import os


data_subfolder = 'data/'
os.chdir(data_subfolder)



#check if file exists
if os.path.exists("autocorr.csv"):
	aut = pd.read_csv("autocorr.csv", header = None).to_numpy()
	

	if aut.size != 0: # check if file is not empty
		fig = plt.figure(figsize=(18,7))
		fig.suptitle("Goodness of MT19937 random-number-generator", fontsize=16)
        
		plt.subplot(121)
		x = np.arange(aut.shape[0])
		plt.scatter(x, aut[:,0], s = 7, color = "Darkblue")
		plt.grid()
		plt.xlabel("N")
		plt.ylabel("Random number")
        
		plt.subplot(122)
		auto = aut[1:500,:]
		x = np.arange(auto.shape[0])
		plt.plot(x,auto[:,1], lw = 0.9, color = "Black")
		plt.xlabel("N")
		plt.ylabel(r"Autocorrelation $\rho$")
		plt.grid()
		plt.savefig("../plots/autocorrelation.png", dpi = 400)

		
	else:
	    print("The CSV file for autocorrelation is empty. Skipping plotting")
	    


from mpl_toolkits.mplot3d import Axes3D

if os.path.exists("goodness_random_triple.csv"):
	
	good = pd.read_csv("goodness_random_triple.csv", header = None).to_numpy()
	if not good.size == 0:
	    fig = plt.figure(figsize = (8,5))
	    fig.suptitle("3D-Projection of randomly generated triplets")
	    ax = fig.add_axes([0,0,1,1], projection = '3d')
	    ax.scatter(good[:,0], good[:,1], good[:,2], '.', color = "Red", s = 16, label = "random-generator")
	    plt.savefig("../plots/goodness_3d.png", dpi = 400)
	    
	else:
	    print("The CSV file for goodness of random numbers is empty. Skipping plotting")




#########################################################




# the infection rate as a function of p1 has been calculated for 30 times and the mean and standard deviation determined. In the following the results are plotted
filename = "average_noise_L16_p2_0.60.csv"
fn2 = "infection_rate_L16_p2_0.60.csv"
if os.path.exists(filename) and os.path.exists(fn2):
    dataa = pd.read_csv(filename, header = None).to_numpy()
    data_rate16_p2 = pd.read_csv("infection_rate_L16_p2_0.60.csv", header = None).to_numpy()
    if dataa.size != 0 and (data_rate16_p2.size != 0):
        plt.figure(figsize=(9, 6))
        plt.scatter(dataa[:,0], dataa[:,1], s = 4, label = r"mean infection rate $\overline{\langle I \rangle}$" +" "+ f"(16x16)")
        plt.plot(dataa[:,0], dataa[:,1], alpha = 0.35)
        plt.fill_between(dataa[:,0], dataa[:,1], dataa[:,1] + dataa[:,2], alpha = 0.3, color = "Skyblue", label = "standard-deviation") 
        plt.fill_between(dataa[:,0], dataa[:,1], dataa[:,1] - dataa[:,2], alpha = 0.3, color = "Skyblue") 
        plt.scatter(data_rate16_p2[:,0], data_rate16_p2[:,1], s=4, label=r"infection rate $\overline{\langle I \rangle}$" + f"(16x16)", color = "black")
        plt.plot(data_rate16_p2[:,0], data_rate16_p2[:,1], alpha = 0.35, color = "black")
        plt.legend()
        plt.grid(True)
        plt.xlabel(r"probability of infection $p_1$")
        plt.ylabel(r"av. infection rate $\overline{\langle I \rangle}$ / [1 / T]")
        plt.title(r"average infection rate (mean over 50 iteration) ($p_2 = 0.6$, $p_3 = 0.3$, $p_4 = 0.0$)")
        
        plt.savefig("../plots/average_noise.png", dpi = 400)
    else:
        print("File is empty")
else: 
    print(f"File {filename} does not exist")

plt.show()

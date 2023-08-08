import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap
import os


data_subfolder = 'data'
os.chdir(data_subfolder)

T = 300 #simulation time -> hyperparameter of the simulation


#plotting average infection rate dependent on p_1 with p2 = 0.3, p3 = 0.3, p4 = 0.0

filename = ["infection_rate_L16_p2_0.30.csv","infection_rate_L32_p2_0.30.csv", "infection_rate_L64_p2_0.30.csv", "infection_rate_L150_p2_0.30.csv"]
check = True
for fn in filename:
	if os.path.exists(fn):
		continue
	else:
		print(f"The following file could not be found: {fn}")
		print(f"Skipping plotting\n")
		check = False
#if files exist
if check:
    data_rate16 = pd.read_csv("infection_rate_L16_p2_0.30.csv", header = None).to_numpy()
    data_rate32 = pd.read_csv("infection_rate_L32_p2_0.30.csv", header = None).to_numpy()
    data_rate64 = pd.read_csv("infection_rate_L64_p2_0.30.csv", header = None).to_numpy()
    data_rate150 = pd.read_csv("infection_rate_L150_p2_0.30.csv", header = None).to_numpy()

    if data_rate16.size != 0 and data_rate32.size != 0 and data_rate64.size != 0 and data_rate150.size != 0: #check that files are not empty
        
        files = [data_rate16, data_rate32,data_rate64,data_rate150]
        L = [16,32,64,150]
        plt.figure(figsize=(9,6))
        for i, file in enumerate(files):
            plt.scatter(file[:,0], file[:,1], s=4, label=r"infection rate $\overline{\langle I \rangle}$" + f"({L[i]}x{L[i]})")
            plt.plot(file[:,0], file[:,1], alpha = 0.35)
            plt.legend(loc = 4)
            plt.grid(True)
            plt.xlabel(r"probability of infection $p_1$")
            plt.ylabel(r"av. infection rate $\overline{\langle I \rangle}$ / [1 / T]")
            plt.title(r"average infection rate ($p_2$ = 0.3, $p_3 = 0.3$) for" + " "+ f"T = {T}")
        plt.savefig("../plots/rate_p2_003.png", dpi = 400)

    else:
        print("One of the CSV files for task 2 is empty. Skipping plotting")

else:
    pass



#plotting average infection rate dependent on p_1 with p2 = 0.6, p3 = 0.3, p4 = 0.0

filename = ["infection_rate_L16_p2_0.60.csv","infection_rate_L32_p2_0.60.csv", "infection_rate_L64_p2_0.60.csv", "infection_rate_L150_p2_0.60.csv"]
check = True
for fn in filename:
	if os.path.exists(fn):
		continue
	else:
		print(f"The following file could not be found: {fn}")
		print(f"Skipping plotting\n")
		check = False
		
if check:
    data_rate16_p2 = pd.read_csv("infection_rate_L16_p2_0.60.csv", header = None).to_numpy()
    data_rate32_p2 = pd.read_csv("infection_rate_L32_p2_0.60.csv", header = None).to_numpy()
    data_rate64_p2 = pd.read_csv("infection_rate_L64_p2_0.60.csv", header = None).to_numpy()
    data_rate150_p2 = pd.read_csv("infection_rate_L150_p2_0.60.csv", header = None).to_numpy()
    
    if data_rate16_p2.size != 0 and data_rate32_p2.size != 0 and data_rate64_p2.size != 0 and data_rate150_p2.size != 0:

        files = [data_rate16_p2, data_rate32_p2, data_rate64_p2, data_rate150_p2]
        L = [16, 32, 64, 150]
        plt.figure(figsize=(9, 6))
        for i, file in enumerate(files):
            plt.scatter(file[:,0], file[:,1], s=4, label=r"infection rate $\overline{\langle I \rangle}$" + f" ({L[i]}x{L[i]})")
            plt.plot(file[:,0], file[:,1], alpha=0.35)
            plt.legend()
            plt.grid(True)
            plt.xlabel(r"probability of infection $p_1$")
            plt.ylabel(r"av. infection rate $\overline{\langle I \rangle}$ / [1 / T]")
            plt.title(r"average infection rate ($p_2$ = 0.6, $p_3 = 0.3$) for" + " "+ f"T = {T}")
        plt.savefig("../plots/rate_p2_006.png", dpi=400)

    else:
        print("One of the CSV files of task 2.2 is empty. Skipping plotting")
else:
    pass



#(vaccination) plotting average infection rate dependent on p_4 with p1 = p2 = p3 = 0.5
filename = ["infection_rate_L16_p2_0.50.csv","infection_rate_L32_p2_0.50.csv", "infection_rate_L64_p2_0.50.csv", "infection_rate_L150_p2_0.50.csv"]
check = True
for fn in filename:
	if os.path.exists(fn):
		continue
	else:
		print(f"The following file could not be found: {filename}")
		print(f"Skipping plotting\n")
		check = False
		
if check:
    data_rate16_vac = pd.read_csv("infection_rate_L16_p2_0.50.csv", header = None).to_numpy()
    data_rate32_vac = pd.read_csv("infection_rate_L32_p2_0.50.csv", header = None).to_numpy()
    data_rate64_vac = pd.read_csv("infection_rate_L64_p2_0.50.csv", header = None).to_numpy()
    data_rate150_vac = pd.read_csv("infection_rate_L150_p2_0.50.csv", header = None).to_numpy()

    if data_rate16_vac.size != 0 and data_rate32_vac.size != 0 and data_rate64_vac.size != 0 and data_rate150_vac.size != 0:

            files = [data_rate16_vac, data_rate32_vac,data_rate64_vac,data_rate150_vac]
            
            L = [16,32,64,150]
            plt.figure(figsize=(9,6))
            for i, file in enumerate(files):
            
                plt.scatter(file[:,0], file[:,1], s = 4, label = r"infection rate $\overline{\langle I \rangle}$" +" "+ f"({L[i]}x{L[i]})")
                plt.plot(file[:,0], file[:,1], alpha = 0.35)
                plt.legend()
                plt.grid(True)
                plt.xlabel(r"probability of vaccination $p_4$")
                plt.ylabel(r"av. infection rate $\overline{\langle I \rangle}$ / [1 / T]")
                plt.title(f"average infection rate with vaccination ($p_1 = p_2 = p_3 = 0.5$) for" + " "+ f"T = {T}")
            plt.savefig("../plots/rate_vacc.png", dpi = 400)
           

    else:
        print("One of the CSV files for vaccination is empty. Skipping plotting")
else:
    pass

plt.show()





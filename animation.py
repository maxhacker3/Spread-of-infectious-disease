from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import time
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
import seaborn as sns

print("Loading file...")
t1 = time.time()
data = pd.read_csv("data/disease_dynamics.csv", header = None).to_numpy()
params = pd.read_csv("data/animation_parameter.csv", header = None).to_numpy()
dt = time.time() - t1
print(f"Loading time: {dt} sec.")


palette = ["white", "gray", "black"] #color map
colormap = ListedColormap(palette)

param = params[:,0]
p1 = round(param[0],1)
p2 = round(param[1],1)
p3 = round(param[2],1)
p4 = round(param[3],1)
L = 250

data_e = data.reshape((101, L, L)) #reshape to an 3d-array (depth -> time propagation)

fig, axes = plt.subplots(1,4, figsize = (46,10))
fig.subplots_adjust(wspace=0.01, hspace=0.5)
axes = axes.ravel()

fig.suptitle(r"Spread of disease [250x250, $(p_1,p_2,p_3,p_4)$ = " + f"({p1}, {p2}, {p3}, {p4})]", fontsize = 30)
img_indices = [0, 10, 30, 80]
for i, img_index in enumerate(img_indices):
    ax = axes[i]
    ax.set_title(f"snapshot for t = {img_index} [TU]", fontsize = 26)
    ax.imshow(data_e[img_index,:,:], cmap= 'binary_r') #inferno ist auch cool
    legend_labels = ["Infected", "Recovered", "Susceptible/Vaccinated"]
    patches = [mpatches.Patch(color=palette[k], label=legend_labels[k]) for k in range(len(palette))]
    ax.legend(handles=patches, bbox_to_anchor=(0.45, 1), loc="upper left", edgecolor = "black", frameon = True, fontsize = 19).get_frame().set_edgecolor('black')
    
fig.savefig(f"plots/propagation_p1_{p1}_p2_{p2}_p3_{p3}_p4_{p4}.png", dpi = 300)
  
    
################################################### animation    
import imageio.v2 as imageio

# specify time values T and dt_snap
total_time = 5
dt_snap = 0.05
# get how many snap-datafiles there should be
num_images = (int)(total_time/dt_snap)

print("Animation is created...")
images = []
for n in range(100):
    fig = plt.figure(figsize=(12,10))
    plt.title(r"Spread of disease (250x250, $(p_1,p_2,p_3,p_4)$ = " + f"({p1}, {p2}, {p3}, {p4})")
    plt.imshow(data_e[n,:,:], cmap= 'binary_r')
    legend_labels = ["Infected", "Recovered", "Susceptible/Vaccinated"]
    patches = [mpatches.Patch(color=palette[i], label=legend_labels[i]) for i in range(len(palette))]
    plt.legend(handles=patches, bbox_to_anchor=(0.71, 1), loc="upper left")
    # save result
    fig.savefig("snapshots/snap_{:04d}.jpeg".format(n), dpi=100)
    # and append it to the img array for imageio at the end
    images.append(imageio.imread("snapshots/snap_{:04d}.jpeg".format(n)))
    plt.close()

# create gif from all results
print("starting work on creating animation...")
imageio.mimsave(f"animations/animation_p1_{p1}_p2_{p2}_p3_{p3}_p4_{p4}.gif", images, duration=0.05)
print("Done! Animation has been created! \n")


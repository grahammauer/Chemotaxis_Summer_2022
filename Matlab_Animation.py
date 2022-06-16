##############################################################################################################################
# Modules
##############################################################################################################################

# Necessary packages
import numpy as np

# Plotting data
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import matplotlib.animation as ani

# Runtime analysis
import time
import datetime

# .csv outputs
import pandas as pd

##############################################################################################################################
# Input parameters (copy directly from MATLAB)
##############################################################################################################################

N      = 400 # Number of parts interval [0,1] is broken into
h      = 1/N # Length of eacch interval
ep1    = 0.1 # Coeffificient of diffusion of v
ep2    = ep1 # Coefficient of advection of v
coeff1 = 1   # Logistic term coefficient
T      = 2   # Total system runtime
dt     = h ** 2 * 0.5 # Time step

alpha  = np.arange(0.25,2.1,.25) # Logistic growth exponent
alpha  = np.array([0.12, 0.25, 0.38, 0.50, 0.62, 0.75, 0.88, 1.00, 1.12, 1.25, 1.38, 1.50, 1.62, 1.75, 1.88, 2.00, 10.00, 100.00])
print(f'alpha = {alpha}')

write_to_csv = 1000; # Write the output to a .csv file once every __ outputs

##############################################################################################################################
# Importing .csv files
##############################################################################################################################

# Create the filenames for u and v
filenames_u = []
for i in range(len(alpha)):
    filenames_u.append(f'u_{T:.2f}_{alpha[i]:.2f}.csv')

filenames_v = []
for i in range(len(alpha)):
    filenames_v.append(f'v_{T:.2f}_{alpha[i]:.2f}.csv')


x, y = np.shape(np.genfromtxt(filenames_u[i], delimiter=','))

u_grid = np.zeros([len(filenames_u),x,y])
for i in range(len(filenames_u)):
    u_grid[i] = np.genfromtxt(filenames_u[i], delimiter=',')

x, y = np.shape(np.genfromtxt(filenames_v[i], delimiter=','))

v_grid = np.zeros([len(filenames_v),x,y])
for i in range(len(filenames_v)):
    v_grid[i] = np.genfromtxt(filenames_v[i], delimiter=',')



##############################################################################################################################
# Test 3D Plot
##############################################################################################################################

fig = plt.figure(figsize=(10,5))
ax = plt.axes(projection='3d')

for i in range(len(u_grid)-2):
    zline = u_grid[i,-1]
    xline = np.linspace(0,1,N)
    yline = alpha[i] * np.ones_like(xline)
    ax.plot3D(xline, yline, zline)
    ax.set_xlabel('X')
    ax.set_ylabel('\u03B1')
    ax.set_zlabel('Intensity')
    ax.set_title(f'Test_Plot')
    ax.set_xlim(0,1)
    #ax.set_ylim(np.min(alpha) - 0.1 * (np.max(alpha) - np.min(alpha)), np.max(alpha)+ 0.1 * (np.max(alpha) - np.min(alpha)))
    #ax.set_zlim(0,1)

plt.savefig('Test_plot_u.png', dpi = 170)

fig = plt.figure(figsize=(10,5))
ax = plt.axes(projection='3d')

for i in range(len(v_grid)-2):
    zline = v_grid[i,-1]
    xline = np.linspace(0,1,N+1)
    yline = alpha[i] * np.ones_like(xline)
    ax.plot3D(xline, yline, zline)
    ax.set_xlabel('X')
    ax.set_ylabel('\u03B1')
    ax.set_zlabel('Intensity')
    ax.set_title(f'Test_Plot')
    ax.set_xlim(0,1)
    #ax.set_ylim(np.min(alpha) - 0.1 * (np.max(alpha) - np.min(alpha)), np.max(alpha)+ 0.1 * (np.max(alpha) - np.min(alpha)))
    #ax.set_zlim(0,1)

plt.savefig('Test_plot_v.png', dpi = 170)

print('Starting Animation')

colors = ['tab:blue', 'tab:red', 'tab:green', 'tab:cyan', 'tab:olive', 'tab:orange', 'tab:purple', 'tab:pink', 'blue', 'red', 'green', 'cyan', 'olive', 'orange', 'purple', 'pink']

# Create the animation
writer = ani.FFMpegWriter(fps=45)
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(projection='3d')
with writer.saving(fig, 'Different_Alpha_u_3D.mp4', 300):
    
    # Create a frame for each timestep
    for i in range(len(u_grid[0])):
        print(f'3D Running {i}')
        
        # Clear the plot
        plt.cla()
            
        ax.set_xlabel('X')
        ax.set_ylabel('\u03B1')
        ax.set_zlabel('Intensity')
        ax.set_title(f'Population Density for Different Values of \u03B1 \nt = {i*write_to_csv*dt:0.2f} s')
        ax.set_xlim(0,1)

        # Plot each mass in a different color on the plot
        for j in range(len(u_grid)-2):
            zline = u_grid[j,i]
            xline = np.linspace(0,1,N)
            yline = alpha[j] * np.ones_like(xline)
            ax.plot3D(xline, yline, zline, c=colors[j], label=f'\u03B1 = {alpha[j]}')
            ax.plot3D(xline, yline, np.zeros_like(zline), c=colors[j], alpha=0.5)       
        ax.legend(fontsize='small')
        writer.grab_frame()

plt.figure(figsize=(10,5), dpi=300)
for i in range(len(u_grid)-2):
    plt.plot(np.linspace(0,1,N), u_grid[i,-1], c=colors[i], label=f'\u03B1 = {alpha[i]}')
plt.legend([f'\u03B1 = 0.12',f'\u03B1 = 0.25',f'\u03B1 = 0.38',f'\u03B1 = 0.50'], loc='upper right')
plt.xlabel('X')
plt.ylabel('Intensity')
plt.xlim(0,1)
plt.title(f'Population Density for Different Values of \u03B1 at t = 2.00 s')
plt.savefig('Different_Alpha_u.png', dpi=170)

# Create the animation
writer = ani.FFMpegWriter(fps=45)
fig = plt.figure(figsize=(10,5), dpi=300)
with writer.saving(fig, 'Different_Alpha_u_2D.mp4', 200):
    
    # Create a frame for each timestep
    for i in range(len(u_grid[0])):
        print(f'2D Running {i}')

        # Clear the plot
        plt.cla()
        
        # Plot data
        plt.plot(np.linspace(0,1,N), u_grid[3,i], label='Population Density', c = 'tab:blue')
        plt.plot(np.linspace(0,1,N+1), v_grid[3,i], label='Chemical Concentration', c = 'tab:orange')

        # Label and title axes
        plt.title(f'Population Density and Chemical Concentration for \u03B1 = 1.00 \nt = {i*write_to_csv*dt:0.2f} s')
        plt.xlabel('x')
        plt.ylabel('Magnitude')
        plt.xlim(0, 1)
        plt.ylim(0, 2)
        plt.legend()
                
        # Animator adds frame
        writer.grab_frame()

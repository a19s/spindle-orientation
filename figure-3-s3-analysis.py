'''
Script used for data analysis linked to Figure 3 and S3
'''

# %% General imports and setup

import matplotlib
import os
import csv
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['figure.figsize'] = [6, 4] #control size of figures in the notebook
this_directory = os.getcwd()
main_path  = this_directory
#### define a function which collects data from csv file ####
def load_data_sim(file_name):
    with open(file_name) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
        this_list_c = []
        for row in readCSV:
            this_list_c.append(row)
    return np.array(this_list_c)
#### dictionaries that will contain data####
#### LOAD DATA  ####
c_data_sim = {}
for k in range(1, 13):
    path = main_path + '/data/exp_vs_sim/cell_'+str(k)+'_simulation.csv'
    c_data_sim[k]= load_data_sim(path)
c_data_exp = {}
for k in range(1, 13):
    path = main_path + '/data/exp_vs_sim/cell_'+str(k)+'_experiment.csv'
    c_data_exp[k]= load_data_sim(path)

# A function that will be used to calculate displacements of the centrosome
# L is the size of the domain
# if a displacement is larger than L/2, it means that the point has jumped
# through the periodic boundary condition and L must be added
# or substracted to correct for this

def correct_displacement(this_displacement, L_value):
    """
    correct displacement by assuming that displacements larger
    than L/2 correspond to a jump across the periodic boundary
     """
    if this_displacement > L_value / 2.0:
        return this_displacement - L_value
    elif this_displacement < - L_value / 2.0:
        return this_displacement + L_value
    return this_displacement


SPINDLE_LENGTH = 4.9  # micron
#
# Load additional experimental data: number of data points,
# position of nucleus and centrosome
#
# The position of the nucleus is obtained from the position
# of the centrosome and the spindle angle

x_array_data = {}
time_array_data = {}
DELTA_X = {}
NB_POINTS = {}
L_VALUE = {}
TOTAL_TIME = {}
NB_TIME_POINTS = {}
centrosome_non_periodic_data = {}
nucleus_position_data = {}
for cell_number in range(1,13):
    #load data for x_array
    file_path = os.path.join(main_path,
                             ('data/exp_vs_sim/Data_Monopolar/Cell'
                             + str(cell_number)
                             + '/Cell'
                             +str(cell_number)
                             + '_lgn_init.dat'))
    lgn_init_data_exp = np.loadtxt(file_path,
                                   delimiter=',')
    x_array_data[cell_number] = lgn_init_data_exp[:, 0]
    #extract spatial information from initial array of concentration
    #the spatial discretization step is taken from experimental data
    DELTA_X[cell_number] = lgn_init_data_exp[1, 0] -  lgn_init_data_exp[0, 0]
    NB_POINTS[cell_number] = len(x_array_data[cell_number])
    L_VALUE[cell_number] = DELTA_X[cell_number] * NB_POINTS[cell_number]
    #load data for centrosome position
    file_path = os.path.join(main_path,'data/exp_vs_sim/Data_Monopolar/Cell' + str(cell_number)
                                 + '/Cell' + str(cell_number) + '_Mech.dat')
    mech_data = np.loadtxt(file_path, delimiter=',')
    time_array_data[cell_number], centrosome_position_data, _, spindle_angle_data = np.transpose(mech_data)
    #extract a non periodic version of the centrosome position
    #does this by assuming that displacements larger
    #than L/2 correspond to a jump across the periodic boundary
    centrosome_displacement = [correct_displacement(this_displacement, L_VALUE[cell_number])
                           for this_displacement
                           in np.diff(centrosome_position_data)]
    centrosome_non_periodic_data[cell_number] = np.concatenate(([centrosome_position_data[0]],
                                                    centrosome_position_data[0]
                                                       + np.cumsum(centrosome_displacement)))
    # the nucleus position is obtained from the centrosome position
    # by using the average spindle length and the angle
    # formed by the axis joining the centrosome to the nucleus
    # relative to the trajectory
    nucleus_position_data[cell_number] = (centrosome_non_periodic_data[cell_number]
                                         + SPINDLE_LENGTH
                                         * np.cos(np.pi / 180.0 * spindle_angle_data))
    TOTAL_TIME[cell_number] = time_array_data[cell_number][-1] - time_array_data[cell_number][0]
    NB_TIME_POINTS[cell_number] = len(time_array_data[cell_number])


### Plot kymographs of LGN concentration, with position of nucleus overlaid,
#   for the 12 cells

# PLOT for FIGURE 3B and S3H

#use the l_value of the first graph
#as a reference to have proportional aspect ratios of the plots
l_value_0 = L_VALUE[1]
def plot_kymographs_with_nucleus(array_1, array_2, total_time, l_value, time_array, nucleus_position):
    """ plot kymograph of concentration, array_1 is the experimental spatiotemporal profile,
    array_2 the simulated spatiotemporal profile
    time_array, nucleus_position contains the data on the position of the nucleus
    as a function of time"""
    fig, axes= plt.subplots(nrows=2,
                            ncols=2,
                            gridspec_kw={"height_ratios":[1, 1], "width_ratios":[1, 0.05]})
    plt.subplots_adjust(hspace=0.8, wspace=0.)
    max_concentration_1 = np.max(np.abs(array_1))
    max_concentration_2 = np.max(np.abs(array_2))
    max_concentration = np.max([max_concentration_1, max_concentration_2])
    #plot experimental concentration kymograph
    kymograph_plot=axes[0, 0].imshow(array_1,
                                     clim=[0, max_concentration],
                                     cmap='Reds',
                                     extent=[0, total_time, 0, l_value],
                                     origin='lower',
                                     interpolation='none',
                                     aspect=l_value_0/l_value)
    #plot experimental nucleus position on top of experimental kymograph
    axes[0, 0].scatter(time_array,
                       nucleus_position % l_value,
                       facecolors='none',
                       edgecolors='grey',
                       s=3,
                       marker='o')
    axes[0, 0].set_ylabel('x position (um)', fontsize=13)
    axes[0, 0].set_title('Experiment', fontsize=13)
    axes[0, 0].set_ylim([0, l_value])
    axes[0, 0].set_xlim([0, total_time])
    #plot simulated concentration kymograph
    kymograph_plot=axes[1, 0].imshow(array_2,
                                     clim=[0, max_concentration],
                                     cmap='Reds',
                                     extent=[0, total_time, 0, l_value],
                                     origin='lower',
                                     interpolation='none',
                                     aspect=l_value_0/l_value)
    #plot experimental nucleus position on top of simulation kymograph
    axes[1, 0].scatter(time_array,
                       nucleus_position % l_value,
                       s=3,
                       color='grey',
                       facecolors='none',
                       edgecolors='grey')
    axes[1, 0].set_xlabel('time (min)', fontsize=13)
    axes[1, 0].set_ylabel('x position (um)', fontsize=13)
    axes[1, 0].set_title('Simulation', fontsize=13)
    axes[1, 0].set_ylim([0, l_value])
    axes[1, 0].set_xlim([0, total_time])
    #plot the colorbars
    fig.colorbar(kymograph_plot, cax=axes[0, 1], ax=axes[0,0])
    fig.colorbar(kymograph_plot, cax=axes[1, 1], ax=axes[1,0]) 
#all plots are normalized to spatiotemporal averages
for cell_number in range(1,13):
    data_1=c_data_exp[cell_number]/np.mean(c_data_exp[cell_number])
    data_2=c_data_sim[cell_number]/np.mean(c_data_sim[cell_number])
    plot_kymographs_with_nucleus(np.transpose(data_1),
                                 np.transpose(data_2),
                                 TOTAL_TIME[cell_number],
                                 L_VALUE[cell_number],
                                 time_array_data[cell_number],
                                 nucleus_position_data[cell_number])
    plt.tight_layout()
    plt.show()



### Extract the velocity parameters from fitting to nucleus displacement

### define a function which collects data from csv file ####
def load_data_MT(file_name):
    with open(file_name) as csvfile:
        this_list_distance = {}
        has_header = csv.Sniffer().has_header(csvfile.read(1024))
        csvfile.seek(0)
        readCSV = csv.reader(csvfile, delimiter=',')
        if has_header:
            next(readCSV)  # Skip header row.
        index_cell = 0
        for row in readCSV:
            this_cell = int(row[0]) #cell label
            if not this_cell==index_cell:
                index_cell = this_cell
                this_list_distance[index_cell] = []
            this_list_distance[this_cell].append(float(row[1]))
    return this_list_distance
#### LOAD DATA  ####
path = main_path + '/data/microtubule_end_distribution.csv'
MT_data_distance = load_data_MT(path)

MT_data_all = []
for this_cell in MT_data_distance.keys():
    MT_data_all.extend(MT_data_distance[this_cell])
(MT_values, MT_bin_edges, _) = plt.hist(MT_data_all,bins=30, density=True)
MT_bin_means = (MT_bin_edges[1:] + MT_bin_edges[:-1])/2.0

# Define a function that is used to fit the probability distribution of MT ends

def log_normal(x, mu, sigma):
    """return a variable with log-normal distribution"""
    result= np.zeros_like(x)
    result =   (1.0 /
                     (x * sigma * np.sqrt(2.0 * np.pi))
                     * np.exp(-(np.log(x) - mu) ** 2 / 2.0 / (sigma ** 2)))
    return result

# Perform a fit of the probability distribution of MT end distance
# from the centrosome; the unknown parameters to determine are mu and sigma
# in the log normal function

from scipy.optimize import minimize
#perform the fit to experimental data
res = minimize(lambda x: np.sum(np.square(log_normal(np.array(MT_bin_means), x[0], x[1])
                                          - np.array(MT_values))), [1, 1])
x_min = res.x
print('value of fitted parameters '+ str(x_min))
print('exponential of values of fitted parameters '+ str(np.exp(x_min)))
distance_array_fine_grained = np.linspace(0.00001, 20, 100)
plt.figure()
#plot original data
plt.scatter(MT_bin_means, MT_values)
#plot corresponding fit
plt.plot(distance_array_fine_grained,
         log_normal(distance_array_fine_grained, x_min[0], x_min[1]),
         lw=3,
         color='red')
plt.xlabel('distance from centrosome (um)')
plt.ylabel('MT end probability distribution')
plt.title('MT end probability distribution\n according to distance to centrosome,\n with fitting function')
plt.tight_layout()
plt.show()

# Define the corresponding distribution of MT: here half the MT ends are on
# the positive side, half on the negative side
# The distribution is also signed to take into account the reverse orientation
# of microtubules on each side.

def signed_log_normal(x, mu, sigma):
    """return a variable with log-normal distribution"""
    #some manipulation is required to avoid evaluating the log in 0
    result= np.zeros_like(x)
    mask = (x!=0)
    nonzero_x = x[mask]
    result[mask] =  0.5* (1.0 /
                     (nonzero_x * sigma * np.sqrt(2.0 * np.pi))
                     * np.exp(-(np.log(abs(nonzero_x)) - mu) ** 2 / 2.0 / (sigma ** 2)))
    return result
#values returned by the fit
MT_MU = np.log(10.28)
MT_SIGMA = np.log(1.26)

# For the calculation of force, define a periodic_x_array which contains the center,
# left and right periodic boxes
#
# This is to make sure that when calculating the integral of the force,
# all the relevant concentration data is included

periodic_x_array = {}
nMT_array = {}
normalization_verification = []
for cell_number in range(1,13):
    #define the x_array which contains 3 boxes
    periodic_x_array[cell_number] = np.concatenate((x_array_data[cell_number]
                                                    - L_VALUE[cell_number],
                                                    x_array_data[cell_number],
                                                    x_array_data[cell_number]
                                                    + L_VALUE[cell_number]))
    #corresponding distribution of microtubules, at different times
    nMT_array[cell_number] = []
    for this_time_index in range(len(centrosome_non_periodic_data[cell_number])):
        nMT_array[cell_number].append(signed_log_normal(periodic_x_array[cell_number]
                                                        - centrosome_non_periodic_data[cell_number][this_time_index]%L_VALUE[cell_number],
                                                        MT_MU,
                                                        MT_SIGMA))
    dx= periodic_x_array[cell_number][1]-periodic_x_array[cell_number][0]
    normalization_verification.append(np.sum(abs(nMT_array[cell_number][0]))*dx)
print(normalization_verification)


def plot_c_array_sim_n_MT(cell_number, this_time_index):
    """ plot LGN concentration profile and signed MT density
    for cell number cell_number, at time point this_time_index"""
    c_array =  c_data_sim[cell_number]
    this_c_array = c_array[this_time_index]
    c_3boxes = np.concatenate((this_c_array, this_c_array, this_c_array))
    fig, ax1=plt.subplots()
    #plot concentration profile
    ax1.plot(periodic_x_array[cell_number],
             c_3boxes,
             lw=3,
             label='simulated LGN',
             color='tab:red')
    ax1.set_ylabel('concentration',
                   fontsize=30,
                   color='tab:red')
    ax1.set_xlim([nucleus_position_data[cell_number][this_time_index]%L_VALUE[cell_number]-L_VALUE[cell_number]/2.0,
                  nucleus_position_data[cell_number][this_time_index]%L_VALUE[cell_number]+L_VALUE[cell_number]/2.0])
    ax1.set_ylim([0,
                  1.2*np.max(c_array)])
    ax1.tick_params(axis='x',
                    labelsize=30)
    ax1.tick_params(axis='y',
                    labelcolor='tab:red',
                    labelsize=30)
    ax2 = ax1.twinx()
    #plot MT signed density
    ax2.plot(periodic_x_array[cell_number],
             nMT_array[cell_number][this_time_index],
             lw=3,
             label='signed MT distribution',
             color='green')
    ax2.set_ylabel('pMT',
                   color='tab:green',
                   fontsize=30)
    ax2.tick_params(axis='y',
                    labelcolor='tab:green',
                    labelsize=30)
    #plot a vertical line indicating the centrosome position
    nucleus_position = nucleus_position_data[cell_number][this_time_index]%L_VALUE[cell_number]
    centrosome_position = centrosome_non_periodic_data[cell_number][this_time_index]%L_VALUE[cell_number]
    plt.axvline(x=centrosome_position, color='black')
    plt.axvline(x=nucleus_position, color='gray')
    plt.xlabel('position (um)', fontsize=30)
    plt.title('concentration field and signed MT distribution,\n cell '
              +str(cell_number)
             +', time '
             +str(time_array_data[cell_number][this_time_index])
             +' mins\n',
             fontsize=20)

# PLOT for FIGURE 3C

cell_number = 1
this_time_index = 61
plot_c_array_sim_n_MT(cell_number, this_time_index)
plt.tight_layout()
plt.show()


# Calculated predicted force, from the simulated LGN concentration profile

force_array_sim = {}
for cell_number in range(1,13):
    c_array =  c_data_sim[cell_number]
    force_array_sim[cell_number] = []
    dx = x_array_data[cell_number][1] - x_array_data[cell_number][0]
    for this_time_index in range(len(centrosome_non_periodic_data[cell_number])):
        this_c_array = c_array[this_time_index]
        c_3boxes = np.concatenate((this_c_array, this_c_array, this_c_array))
        force_array_sim[cell_number].append(dx * np.sum(c_3boxes * nMT_array[cell_number][this_time_index]))

# PLOT for FIGURE 3D and S3I

V0_array = []
from scipy.optimize import minimize
n_rows = 4
n_cols = 3
fig, axes = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=(5,8))
fig.subplots_adjust(hspace = 0.4, wspace=0.3)
delta_nucleus_position_exp = {}
delta_nucleus_force_sim = {}
for cell_number in range(1, 13):
    nx = (cell_number-1)//n_cols
    ny = (cell_number-1)%n_cols
    dt = time_array_data[cell_number][1] - time_array_data[cell_number][0]
    #the cumulative value of simulated force starts at 0
    cum_force_array_sim = np.cumsum(np.concatenate(([0],
                                                    force_array_sim[cell_number][:-1])))
    #substract the initial position to obtain the nucleus displacement
    delta_nucleus_position_exp[cell_number] = (nucleus_position_data[cell_number]
                                               -nucleus_position_data[cell_number][0])
    delta_nucleus_force_sim[cell_number] = dt * cum_force_array_sim
    #perform fit
    res = minimize(lambda x: np.sum(np.square(delta_nucleus_position_exp[cell_number]
                                              - x * delta_nucleus_force_sim[cell_number])),
                   [5])
    #store the result of the fit in the array of V0s
    V0_array.append(res.x[0])
    axes[nx][ny].scatter(time_array_data[cell_number],
                         delta_nucleus_position_exp[cell_number],
                         s=1,
                         lw=2,
                         color='blue')
    axes[nx][ny].scatter(time_array_data[cell_number],
                         res.x[0]*delta_nucleus_force_sim[cell_number],
                         s=1,
                         lw=2,
                         color='red')
    axes[nx][ny].set_xlabel('time (mins)',fontsize=8)
    axes[nx][ny].set_ylabel('position (um)',fontsize=8)
    axes[nx][ny].grid(True)
    axes[nx][ny].axhline(y=0, color='black')
    axes[nx][ny].axvline(x=0, color='black')
    axes[nx][ny].set_xticks([0, 100, 200])
    axes[nx][ny].set_xticklabels([0, 100, 200], fontsize=8)
    axes[nx][ny].set_xlim([0, 220])
    axes[nx][ny].set_yticks([-400, 0, 400])
    axes[nx][ny].set_yticklabels([-400, 0, 400], fontsize=8)
    axes[nx][ny].set_ylim([-500, 400])
fig.suptitle('Experimental nucleus position (blue)\n and predicted nucleus position (red) for all cells',
             fontsize=8)
plt.show()
print('array of v0 factors:')
print(V0_array)
print('mean v0 factors: ' + str(np.mean(np.array(V0_array))))
print('std v0 factors: ' + str(np.std(np.array(V0_array))))
print('sem v0 factors: ' + str(np.std(np.array(V0_array)) / np.sqrt(len(V0_array)-1)))


### Analysis of the nucleus motion, with the MT experimental full density

#### define a function which collects data from csv file ####
def load_data_MT(file_name):
    with open(file_name) as csvfile:
        this_list_distance = {}
        this_list_c = {}
        has_header = csv.Sniffer().has_header(csvfile.read(1024))
        csvfile.seek(0)
        readCSV = csv.reader(csvfile, delimiter=',')
        if has_header:
            next(readCSV)  # Skip header row.
        index_cell = 0
        for row in readCSV:
            if row:
                this_cell = int(row[0]) #cell label
                if not this_cell==index_cell:
                    index_cell = this_cell
                    this_list_c[index_cell] = []
                    this_list_distance[index_cell] = []
                this_list_distance[this_cell].append(float(row[1]))
                this_list_c[this_cell].append(float(row[2]))
    return [this_list_distance, this_list_c]
#### dictionaries that will contain data####
#### LOAD DATA  ####
path = main_path + '/data/data_microtubules.csv'
[MT_data_distance, MT_data_c] = load_data_MT(path)

#find the minimum common length to all data
max_length = min([len(MT_data_distance[this_cell]) for this_cell in MT_data_c.keys()])
MT_data_distance_mean = []
MT_data_c_mean = []
MT_data_c_std = []
for k in range(max_length):
    MT_data_distance_mean.append(MT_data_distance[1][k])
    MT_data_c_mean.append(np.mean([MT_data_c[this_cell][k]
                                  for this_cell
                                  in MT_data_c.keys()]))
    MT_data_c_std.append(np.std([MT_data_c[this_cell][k]
                                  for this_cell
                                  in MT_data_c.keys()]))
plt.xlabel('distance from centrosome (um)')
plt.ylabel('normalized MT intensity')
plt.errorbar(MT_data_distance_mean, MT_data_c_mean, MT_data_c_std, lw=3)
plt.title('mean MT intensity according to distance to centrosome')
plt.show()


# Perform a fit to the profile of MT intensit away from the centrosome
#
# Choose a simple function with a constant c0+cinf up to distance d0,
# and then a decaying exponential profile, going to a constant cinf

from scipy.optimize import minimize
def fitted_function(distance, cinf, c0, d0, d):
    func = []
    for this_distance in distance:
        if this_distance < d0:
            func.append(c0 + cinf)
        else:
            func.append(cinf + c0 * np.exp(- (this_distance - d0) / d) )
    return func
res = minimize(lambda x: np.sum(np.square(fitted_function(np.array(MT_data_distance_mean),
                                                          x[0],
                                                          x[1],
                                                          x[2],
                                                          x[3])
                                          - np.array(MT_data_c_mean))),
               [0, 1, 1, 1])
x_min = res.x
print('value of fitted parameters '+ str(x_min))
distance_array_fine_grained = np.linspace(0, 20, 100)
plt.figure()
plt.errorbar(MT_data_distance_mean,
             MT_data_c_mean,
             MT_data_c_std,
             lw=3)
plt.plot(distance_array_fine_grained,
         fitted_function(distance_array_fine_grained,
                         x_min[0],
                         x_min[1],
                         x_min[2],
                         x_min[3]),
        lw=5)
plt.xlabel('distance from centrosome')
plt.ylabel('normalized MT intensity')
plt.title('mean MT intensity according to distance to centrosome, with fitting function')
plt.tight_layout()
plt.show()


# Now use this function to reanalyse the force: use this function
# to have a signed normalized version

def signed_MT_density(distance, d0, d):
    output = []
    for this_distance in distance:
        if(abs(this_distance) < d0):
            output.append(1.0 / 2.0 * np.sign(this_distance) / (d0 + d))
        else:
            output.append(1.0 / 2.0
                          * np.sign(this_distance)
                          * np.exp( - (np.abs(this_distance) - d0) / d)
                          /(d0 + d))
    return output
distance_array_fine_grained = np.linspace(-50, 50, 10000)
signed_MT_density_to_plot = signed_MT_density(distance_array_fine_grained, x_min[2], x_min[3])
#plotting the MT distribution
plt.figure()
plt.plot(distance_array_fine_grained,
         signed_MT_density_to_plot,
         lw = 3)
plt.xlim([-20, 20])
plt.xlabel('distance from centrosome (um)', fontsize=20)
plt.ylabel('normalized MT \nprobability density', fontsize=20)
plt.tight_layout()
plt.show()

#verifying that the integral of the absolute value of the signed probability distribution is 1
integral_value = (np.sum(np.abs(np.array(signed_MT_density_to_plot)))
                *(distance_array_fine_grained[1]- distance_array_fine_grained[0]))
print('the integral is '+str(integral_value))

periodic_x_array = {}
nMT_array_testing = {}
normalization_verification = []
for cell_number in range(1,13):
    periodic_x_array[cell_number] = np.concatenate((x_array_data[cell_number]
                                                    - L_VALUE[cell_number],
                                                    x_array_data[cell_number],
                                                    x_array_data[cell_number]
                                                    + L_VALUE[cell_number]))
    nMT_array_testing[cell_number] = []
    for this_time_index in range(len(centrosome_non_periodic_data[cell_number])):
        nMT_array_testing[cell_number].append(signed_MT_density(periodic_x_array[cell_number]
                                               - centrosome_non_periodic_data[cell_number][this_time_index]%L_VALUE[cell_number] ,
                                                x_min[2],
                                                x_min[3]))
    dx= periodic_x_array[cell_number][1]-periodic_x_array[cell_number][0]
    normalization_verification.append(np.sum(np.abs(np.array(nMT_array_testing[cell_number][0])))*dx)
print(normalization_verification)


force_array_sim_testing = {}
for cell_number in range(1,13):
    c_array =  c_data_sim[cell_number]
    force_array_sim_testing[cell_number] = []
    dx = x_array_data[cell_number][1] - x_array_data[cell_number][0]
    for this_time_index in range(len(centrosome_non_periodic_data[cell_number])):
        this_c_array = c_array[this_time_index]
        c_3boxes = np.concatenate((this_c_array, this_c_array, this_c_array))
        force_array_sim_testing[cell_number].append(dx * np.sum(c_3boxes * nMT_array_testing[cell_number][this_time_index]))

# PLOT for FIGURE S3J

V0_array = []
from scipy.optimize import minimize
n_rows = 4
n_cols = 3
fig, axes = plt.subplots(n_rows, n_cols, sharex=True, sharey=True,figsize=(5,8))
fig.subplots_adjust(hspace = 0.4, wspace=0.3)
for cell_number in range(1, 13):
    #to determine where to plot the cell
    nx=(cell_number-1)//n_cols
    ny=(cell_number-1)%n_cols
    dt=time_array_data[cell_number][1]- time_array_data[cell_number][0]
    #cumulative sum starts from 0
    cum_force_array_sim_testing = np.cumsum(np.concatenate(([0],
                                                            force_array_sim_testing[cell_number][:-1])))
    #take the displacement as the position relative to the position at time 0
    delta_nucleus_position_exp = nucleus_position_data[cell_number]-nucleus_position_data[cell_number][0]
    delta_nucleus_force_sim_testing = dt*cum_force_array_sim_testing
    #actually perform the fit
    res = minimize(lambda x: np.sum(np.square(delta_nucleus_position_exp - x * delta_nucleus_force_sim_testing)), [5])
    V0_array.append(res.x[0])
    axes[nx][ny].scatter(time_array_data[cell_number],
                         delta_nucleus_position_exp,
                         s=1,
                         lw=2,
                         color='blue')
    axes[nx][ny].scatter(time_array_data[cell_number],
                         res.x*delta_nucleus_force_sim_testing,
                         s=1,
                         lw=2,
                         color='red')
    axes[nx][ny].set_xlabel('time (mins)',fontsize=8)
    axes[nx][ny].set_ylabel('position (um)',fontsize=8)
    axes[nx][ny].grid(True)
    axes[nx][ny].axhline(y=0, color='black')
    axes[nx][ny].axvline(x=0, color='black')
    axes[nx][ny].set_xticks([0, 100, 200])
    axes[nx][ny].set_xticklabels([0, 100, 200], fontsize=8)
    axes[nx][ny].set_xlim([0, 220])
    axes[nx][ny].set_yticks([-400, 0, 400])
    axes[nx][ny].set_yticklabels([-400, 0, 400], fontsize=8)
    axes[nx][ny].set_ylim([-500, 400])
fig.suptitle('Experimental nucleus position (blue)\n and predicted nucleus position (red) for all cells',fontsize=8)
plt.tight_layout()
plt.show()
print('array of v0 factors:')
print(V0_array)
print('mean v0 factors: '+str(np.mean(np.array(V0_array))))
print('std v0 factors: '+str(np.std(np.array(V0_array))))
print('sem v0 factors: '+str(np.std(np.array(V0_array))/np.sqrt(len(V0_array)-1.0)))
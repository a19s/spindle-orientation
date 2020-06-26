'''
Script used for data analysis linked to Figure 4 and S4
'''

# %% General imports and setup

import os
import csv
import matplotlib.pyplot as plt
import numpy as np
plt.close('all')
this_directory = os.getcwd()
main_path  = this_directory
#### define a function which collects data from csv file ####
def load_data_shapes(file_name):
    cell_shape_data = []
    with open(file_name) as csvfile:
        has_header = csv.Sniffer().has_header(csvfile.read(1024))
        csvfile.seek(0)
        readCSV = csv.reader(csvfile, delimiter=',')
        if has_header:
            next(readCSV)  # Skip header row.
        cell_id = 1
        this_list_time = []
        this_list_long_axis = []
        this_list_short_axis = []
        for row in readCSV:
            if len(row) == 0: #Skip empty line
                continue
            if not int(row[1]) == cell_id:
                #save current list
                if(len(this_list_time)>0):
                    cell_shape_data.append([np.array(this_list_time),
                                            np.array(this_list_long_axis),
                                            np.array(this_list_short_axis)])
                #create new list
                cell_id = int(row[1])
                this_list_time = []
                this_list_short_axis = []
                this_list_long_axis = []
            this_list_time.append(float(row[0]))
            this_list_long_axis.append(float(row[2])/2.0) ##factor 2 tof obtain semi-long axis
            this_list_short_axis.append(float(row[3])/2.0) ##factor 2 to obtain semi-short axis
         #save final list at the end of the file
        if(len(this_list_time)>0):
            cell_shape_data.append([np.array(this_list_time),
                                    np.array(this_list_long_axis),
                                    np.array(this_list_short_axis)])
    return cell_shape_data
#### dictionaries that will contain data####
#### LOAD DATA  ####
path = main_path + '/data/Micropattern_shape.csv'
cell_shape_data_micropattern = load_data_shapes(path)

# %% Data preparation for Figure S4C

micropattern_pairs = []
for k in range(len(cell_shape_data_micropattern)):
    micropattern_pairs.extend([[this_time, this_long_length, this_short_length]
                                    for this_time, this_long_length, this_short_length
                                    in zip(cell_shape_data_micropattern[k][0],
                                           cell_shape_data_micropattern[k][1],
                                           cell_shape_data_micropattern[k][2])])
from collections import defaultdict
long_axis_micropattern_by_time = defaultdict(list)
short_axis_micropattern_by_time = defaultdict(list)
for this_pair in micropattern_pairs:
    long_axis_micropattern_by_time[int(this_pair[0])].append(this_pair[1])
    short_axis_micropattern_by_time[int(this_pair[0])].append(this_pair[2])


time_array_micropattern= [k for k in range(-21, 16, 3)]
long_axis_micropattern_mean = []
long_axis_micropattern_std = []
short_axis_micropattern_mean = []
short_axis_micropattern_std = []
for this_time in time_array_micropattern:
    long_axis_micropattern_mean.append(np.nanmean(np.array(long_axis_micropattern_by_time[this_time])))
    long_axis_micropattern_std.append(np.nanstd(np.array(long_axis_micropattern_by_time[this_time])))
    short_axis_micropattern_mean.append(np.nanmean(np.array(short_axis_micropattern_by_time[this_time])))
    short_axis_micropattern_std.append(np.nanstd(np.array(short_axis_micropattern_by_time[this_time])))


#perform fit to experimental data of cell shape
from scipy.optimize import minimize
def fitted_function(time, L, Lm0, tm0, taum):
    return L + (Lm0 - L) * (1.0 - np.tanh((time - tm0) / taum))/2.0
res = minimize(lambda x: np.sum(np.square(fitted_function(time_array_micropattern, x[0], x[1], x[2], x[3])
                                          - long_axis_micropattern_mean)
                                +np.square(fitted_function(time_array_micropattern, x[0], x[4], x[5], x[6])
                                          - short_axis_micropattern_mean)
                               )
               , [11, 20, 5, 5,10,5,3])
x_min_long_short = res.x
# Figure S4C
plt.figure()
plt.errorbar(time_array_micropattern,
             long_axis_micropattern_mean,
             long_axis_micropattern_std,
             marker='o',
             capsize=3)
plt.errorbar(time_array_micropattern,
             short_axis_micropattern_mean,
             short_axis_micropattern_std,
             marker='o',
             capsize=3)
time_array_fine_grained = np.linspace(-21, 16, 100)
plt.plot(time_array_fine_grained, fitted_function(time_array_fine_grained,
                                                  x_min_long_short[0],
                                                  x_min_long_short[1],
                                                  x_min_long_short[2],
                                                  x_min_long_short[3]),
         lw=3,
         color='blue')
plt.plot(time_array_fine_grained, fitted_function(time_array_fine_grained,
                                                  x_min_long_short[0],
                                                  x_min_long_short[4],
                                                  x_min_long_short[5],
                                                  x_min_long_short[6]),
         lw=3,
         color='orange')
plt.title('Cell semi-short and semi-long axis\n as a function of time, micropatterns', fontsize=15)
plt.xlabel('Time relative to NEB (minutes)', fontsize=15)
plt.ylabel('Length (um)', fontsize=15)
plt.box(False)
plt.grid(True)
plt.axhline(y=0, color='black')
plt.axvline(x=0, color='black')
plt.xticks(fontsize=15, rotation=0)
plt.yticks(fontsize=15, rotation=0)
plt.tight_layout()
plt.show()

# Data load and preparation for Figure S4D-E

def load_data_shapes(file_name):
    cell_shape_data = []
    with open(file_name) as csvfile:
        has_header = csv.Sniffer().has_header(csvfile.read(1024))
        csvfile.seek(0)
        readCSV = csv.reader(csvfile, delimiter=',')
        if has_header:
            next(readCSV)  # Skip header row.
        cell_id = 1
        this_list_time = []
        this_list_long_axis = []
        this_list_short_axis = []
        this_list_aspect_ratio = []
        for row in readCSV:
            if len(row) == 0: #Skip empty line
                continue
            if not int(row[0]) == cell_id:
                #save current list
                if(len(this_list_time)>0):
                    cell_shape_data.append([np.array(this_list_time),
                                            np.array(this_list_long_axis),
                                            np.array(this_list_short_axis),
                                            np.array(this_list_aspect_ratio)])
                #create new list
                cell_id = int(row[0])
                this_list_time = []
                this_list_short_axis = []
                this_list_long_axis = []
                this_list_aspect_ratio = []
            this_list_time.append(float(row[1]))
            this_list_long_axis.append(float(row[2]) / 2.0) #divide by 2 to obtain the semi-long axis
            this_list_short_axis.append(float(row[3]) / 2.0) #divide by 2 to obtain the semi-short axis
            this_list_aspect_ratio.append(float(row[2])/float(row[3])) #ratio: long over short
         #save final list at the end of the file
        if(len(this_list_time)>0):
            cell_shape_data.append([np.array(this_list_time),
                                    np.array(this_list_long_axis),
                                    np.array(this_list_short_axis),
                                    np.array(this_list_aspect_ratio)])
    return cell_shape_data

path = main_path + '/data/data_dna_lengths.csv'
DNA_shape_data = load_data_shapes(path)

DNA_pairs = []
for k in range(len(DNA_shape_data)):
    DNA_pairs.extend([[this_time, this_long_length, this_short_length, this_aspect_ratio]
                                    for this_time, this_long_length, this_short_length, this_aspect_ratio
                                    in zip(DNA_shape_data[k][0],
                                           DNA_shape_data[k][1],
                                           DNA_shape_data[k][2],
                                           DNA_shape_data[k][3])])
from collections import defaultdict

long_axis_DNA_by_time = defaultdict(list)
short_axis_DNA_by_time = defaultdict(list)
aspect_ratio_DNA_by_time = defaultdict(list)
for this_pair in DNA_pairs:
    long_axis_DNA_by_time[int(this_pair[0])].append(this_pair[1])
    short_axis_DNA_by_time[int(this_pair[0])].append(this_pair[2])
    aspect_ratio_DNA_by_time[int(this_pair[0])].append(this_pair[3])

time_array_DNA= [k for k in range(-15, 16, 3)]
long_axis_DNA_mean = []
long_axis_DNA_std = []
short_axis_DNA_mean = []
short_axis_DNA_std = []
aspect_ratio_DNA_mean = []
aspect_ratio_DNA_std = []
for this_time in time_array_DNA:
    long_axis_DNA_mean.append(np.nanmean(np.array(long_axis_DNA_by_time[this_time])))
    long_axis_DNA_std.append(np.nanstd(np.array(long_axis_DNA_by_time[this_time])))
    short_axis_DNA_mean.append(np.nanmean(np.array(short_axis_DNA_by_time[this_time])))
    short_axis_DNA_std.append(np.nanstd(np.array(short_axis_DNA_by_time[this_time])))
    aspect_ratio_DNA_mean.append(np.nanmean(np.array(aspect_ratio_DNA_by_time[this_time])))
    aspect_ratio_DNA_std.append(np.nanstd(np.array(aspect_ratio_DNA_by_time[this_time])))

# Plot of Figure S4D
plt.figure()
plt.errorbar(time_array_DNA,
             long_axis_DNA_mean,
             long_axis_DNA_std,
             marker='o',
             capsize=3)
plt.errorbar(time_array_DNA,
             short_axis_DNA_mean,
             short_axis_DNA_std,
             marker='o',
             capsize=3)
plt.title('DNA semi-short and semi-long axis as a function of time', fontsize=15)
plt.xlabel('Time relative to NEB (minutes)', fontsize=15)
plt.ylabel('Length (um)', fontsize=15)
plt.box(False)
plt.grid(True)
plt.axhline(y=0, color='black')
plt.axvline(x=0, color='black')
plt.xticks(fontsize=15, rotation=0)
plt.yticks(fontsize=15, rotation=0)
plt.tight_layout()
plt.show()


# Plot of Figure S4E
plt.figure()
plt.errorbar(time_array_DNA,
             aspect_ratio_DNA_mean,
             aspect_ratio_DNA_std,
             marker='o',
             capsize=3)
plt.title('DNA aspect ratio as a function of time', fontsize=15)
plt.xlabel('Time relative to NEB (minutes)', fontsize=15)
plt.ylabel('Aspect ratio', fontsize=15)
plt.box(False)
plt.grid(True)
plt.axhline(y=0, color='black')
plt.axvline(x=0, color='black')
plt.xticks(fontsize=15, rotation=0)
plt.yticks(fontsize=15, rotation=0)
plt.tight_layout()
plt.show()


## %% Data Load and preparation for Figure 4D-E and S4I
#### define a function which collects data from csv file ####
def load_data_LGN(file_name):
    LGN_profile_data = []
    with open(file_name) as csvfile:
        has_header = csv.Sniffer().has_header(csvfile.read(1024))
        csvfile.seek(0)
        readCSV = csv.reader(csvfile, delimiter=',')
        if has_header:
            next(readCSV)  # Skip header row.
        cell_id = 1
        this_list_time = []
        this_list_long_axis = []
        this_list_short_axis = []
        this_list_angles_for_LGN = []
        this_list_LGN = []
        this_list_spindle_angle = []
        for row in readCSV:
            if len(row) == 0: #Skip empty line
                continue
            if not int(row[1]) == cell_id: #when a new cell ID is encountered
                #save current list
                if(len(this_list_time)>0):
                    LGN_profile_data.append([np.array(this_list_time),
                                             np.array(this_list_angles_for_LGN),
                                             np.array(this_list_LGN),
                                             np.array(this_list_long_axis),
                                             np.array(this_list_short_axis),
                                             np.array(this_list_spindle_angle)])
                #create new list for the next cell
                cell_id = int(row[1])
                this_list_time = []
                this_list_long_axis = []
                this_list_short_axis = []
                this_list_angles_for_LGN = []
                this_list_LGN = []
                this_list_spindle_angle = []
            this_list_time.append(float(row[7]))
            this_list_long_axis.append(float(row[4]))
            this_list_short_axis.append(float(row[5]))
            #here angles are in degrees and converted to radians
            this_list_angles_for_LGN.append(np.pi * float(row[8]) / 180.0)
            this_list_LGN.append(float(row[9]))
            #spindle angle; also in degrres
            if(row[11]==''):
                this_list_spindle_angle.append(np.nan)
            else:
                this_list_spindle_angle.append(np.pi * float(row[11]) / 180.0)
        #add the last list at the end of the file
        if(len(this_list_time) > 0):
             LGN_profile_data.append([np.array(this_list_time),
                                     np.array(this_list_angles_for_LGN),
                                     np.array(this_list_LGN),
                                     np.array(this_list_long_axis),
                                     np.array(this_list_short_axis),
                                     np.array(this_list_spindle_angle)])
    return LGN_profile_data
#### dictionaries that will contain data####
#### LOAD DATA  ####
path = main_path + '/data/data_lgn_rounding_bipolar.csv'
LGN_profile_data = load_data_LGN(path)

this_list_angles = []
this_list_LGN = []
all_angles = {}
all_LGN = {}
all_times = {}
all_long_axis = {}
all_short_axis = {}
all_spindle_angle = {}
for k in range(len(LGN_profile_data)):
    all_angles[k] = []
    all_LGN[k] = []
    all_times[k] = []
    all_long_axis[k] = []
    all_short_axis[k] = []
    all_spindle_angle[k] = []
    current_time = LGN_profile_data[k][0][1]
    for this_index, this_time in enumerate(LGN_profile_data[k][0]):
        if not this_time == current_time:
            all_angles[k].append(this_list_angles)
            all_LGN[k].append(this_list_LGN)
            all_times[k].append(current_time)
            all_long_axis[k].append(LGN_profile_data[k][3][this_index-1])#index-1 because one saves the previous encountered value
            all_short_axis[k].append(LGN_profile_data[k][4][this_index-1])
            all_spindle_angle[k].append(LGN_profile_data[k][5][this_index-1])
            this_list_angles = []
            this_list_LGN = []
            current_time = this_time
        this_list_angles.append(LGN_profile_data[k][1][this_index])
        this_list_LGN.append(LGN_profile_data[k][2][this_index])
    all_angles[k].append(this_list_angles)
    all_LGN[k].append(this_list_LGN)
    all_times[k].append(current_time)
    all_long_axis[k].append(LGN_profile_data[k][3][this_index-1])
    all_short_axis[k].append(LGN_profile_data[k][4][this_index-1])
    all_spindle_angle[k].append(LGN_profile_data[k][5][this_index-1])
    this_list_angles = []
    this_list_LGN = []
    current_time = this_time


import copy
all_angles_ordered = copy.deepcopy(all_angles)
all_LGN_ordered = copy.deepcopy(all_LGN)
for this_cell in range(len(all_angles)):
    for this_time in range(len(all_angles[this_cell])):
        this_angles = np.array(all_angles[this_cell][this_time])
        this_LGN = np.array(all_LGN[this_cell][this_time])
        all_angles_ordered[this_cell][this_time] = this_angles[this_angles.argsort()]
        all_LGN_ordered[this_cell][this_time] = this_LGN[this_angles.argsort()]

# Plot of Figure 4C
mean_semi_long_axis = []
mean_semi_short_axis = []
std_semi_long_axis = []
std_semi_short_axis = []
for this_time,_ in enumerate(all_times[0]):
    semi_long_axis_vect=[all_long_axis[k][this_time]/2.0
                         for k in range(len(all_times))]#len(all_times)= number of cells
    semi_short_axis_vect=[all_short_axis[k][this_time]/2.0
                         for k in range(len(all_times))]
    mean_semi_long_axis.append(np.mean(semi_long_axis_vect))
    mean_semi_short_axis.append(np.mean(semi_short_axis_vect))
    std_semi_long_axis.append(np.std(semi_long_axis_vect))
    std_semi_short_axis.append(np.std(semi_short_axis_vect))
plt.errorbar(all_times[0],
             mean_semi_long_axis,
             std_semi_long_axis,
             marker='o',
             capsize=3,
             lw=3)
plt.errorbar(all_times[0],
             mean_semi_short_axis,
             std_semi_short_axis,
             marker='o',
             capsize=3,
             lw=3)
plt.title('Cell semi-short and semi-long axis as a function of time', fontsize=15)
plt.xlabel('Time relative to NEB (minutes)', fontsize=15)
plt.ylabel('Length (um)', fontsize=15)
plt.box(False)
plt.grid(True)
plt.axhline(y=0, color='black')
plt.axvline(x=0, color='black')
plt.xticks(fontsize=15, rotation=0)
plt.yticks(fontsize=15, rotation=0)
# plt.savefig('Graphs/Semi_long_short_axis_LGN_cells.eps', bbox_inches="tight")
plt.tight_layout()
plt.show()

# Plot of Figure S4I
plt.figure()
this_cell_for_plot = 3
for k in range(0, len(all_angles[this_cell_for_plot]), 3):#take a subset of time points
    plt.plot(all_angles_ordered[this_cell_for_plot][k],
             all_LGN_ordered[this_cell_for_plot][k],
             lw=2,
             label=str(all_times[this_cell_for_plot][k])+' mins')
plt.ylim([0, 4.0])
plt.title('Profile of LGN intensity', fontsize=15)
plt.xlabel('Angle relative to long axis (rad)', fontsize=15)
plt.ylabel('Normalized LGN intensity', fontsize=15)
plt.box(False)
plt.grid(True)
plt.axhline(y=0, color='black')
plt.axvline(x=0, color='black')
plt.xticks(fontsize=15, rotation=0)
plt.yticks(fontsize=15, rotation=0)
plt.legend()
plt.show()


# Plot of Figure 4D-E and S4K

def calculate_nematic_Q(list_angles, list_f):
    list_angles_rad = np.array(list_angles)
    list_angles_rad_extended = np.concatenate((list_angles_rad, [(list_angles_rad[0] + 2*np.pi)]))
    list_f_extended = np.concatenate((np.array(list_f), [list_f[0]]))
    Qxx = np.trapz(np.cos(2*list_angles_rad_extended) * list_f_extended,
                   list_angles_rad_extended) / np.trapz(list_f_extended,
                                                        list_angles_rad_extended)
    Qxy = np.trapz(np.sin(2*list_angles_rad_extended ) * list_f_extended,
                   list_angles_rad_extended ) / np.trapz(list_f_extended,
                                                         list_angles_rad_extended)
    return [Qxx, Qxy]
def calculate_nematic_S_phi(list_angles, list_f):
    [Qxx, Qxy] = calculate_nematic_Q(list_angles, list_f)
    S = np.sqrt(Qxx**2 + Qxy**2)
    phi = 1.0/2.0 * np.arctan2(Qxy, Qxx)
    return [S, phi]

#calculate Qxx, S and phi for each cell as a function of time
Qxx_vect = [[] for this_cell in range(len(all_angles_ordered))]
S_vect = [[] for this_cell in range(len(all_angles_ordered))]
phi_vect = [[] for this_cell in range(len(all_angles_ordered))]
a_vect =  [[] for this_cell in range(len(all_angles_ordered))] #alignment vector
for this_cell in range(len(all_angles_ordered)):
    for this_time in range(len(all_angles_ordered[this_cell])):
        [Qxx, Qxy] = calculate_nematic_Q(all_angles_ordered[this_cell][this_time],
                                         all_LGN_ordered[this_cell][this_time])
        S_phi = calculate_nematic_S_phi(all_angles_ordered[this_cell][this_time],
                                        all_LGN_ordered[this_cell][this_time])
        Qxx_vect[this_cell].append(Qxx)
        S_vect[this_cell].append(S_phi[0])
        phi_vect[this_cell].append(S_phi[1])
        a_vect[this_cell].append(np.cos(2.0 * S_phi[1]))
#calculate the mean between different cells
Qxx_vect_mean = []
Qxx_vect_std = []
S_vect_mean = []
S_vect_std = []
phi_vect_mean = []
phi_vect_std = []
a_vect_mean = []
a_vect_std = []
for this_time in range(len(all_angles_ordered[0])):
    Qxx_vect_mean.append(np.mean(np.array([Qxx_vect[this_cell][this_time]
                                         for this_cell
                                         in range(len(all_angles_ordered))])))
    Qxx_vect_std.append(np.std(np.array([Qxx_vect[this_cell][this_time]
                                         for this_cell
                                         in range(len(all_angles_ordered))])))
    S_vect_mean.append(np.mean(np.array([S_vect[this_cell][this_time]
                                         for this_cell
                                         in range(len(all_angles_ordered))])))
    S_vect_std.append(np.std(np.array([S_vect[this_cell][this_time]
                                         for this_cell
                                         in range(len(all_angles_ordered))])))
    phi_vect_mean.append(np.mean(np.array([phi_vect[this_cell][this_time]
                                         for this_cell
                                         in range(len(all_angles_ordered))])))
    phi_vect_std.append(np.std(np.array([phi_vect[this_cell][this_time]
                                         for this_cell
                                         in range(len(all_angles_ordered))])))
    a_vect_mean.append(np.mean(np.array([a_vect[this_cell][this_time]
                                         for this_cell
                                         in range(len(all_angles_ordered))])))
    a_vect_std.append(np.std(np.array([a_vect[this_cell][this_time]
                                         for this_cell
                                         in range(len(all_angles_ordered))])))
#plot the evolution of nematic order S for each cell, and its average
plt.figure()
for this_cell in range(len(all_angles_ordered)):
    plt.plot(all_times[this_cell], S_vect[this_cell])
plt.errorbar(all_times[0], S_vect_mean, S_vect_std,marker='o', color='blue',
             lw=3,
             capsize=3,
             capthick=3)
plt.title('Nematic order parameter S of LGN distribution\n as a function of time\n for different cells', fontsize=15)
plt.xlabel('Time after NEB (minutes)', fontsize=15)
plt.ylabel('S', fontsize=15)
plt.box(False)
plt.grid(True)
plt.axhline(y=0, color='black')
plt.axvline(x=0, color='black')
plt.xticks(fontsize=15, rotation=0)
plt.yticks(fontsize=15, rotation=0)
plt.ylim([0, 0.36])
plt.tight_layout()
plt.show()

plt.figure()
#plot the evolution of nematic angles for each cell
for this_cell in range(len(all_angles_ordered)):
    plt.plot(all_times[this_cell], phi_vect[this_cell])
plt.title('Nematic angle as a function of time', fontsize=15)
plt.xlabel('Time after NEB (minutes)', fontsize=15)
plt.ylabel('phi (rad)', fontsize=15)
plt.box(False)
plt.grid(True)
plt.axhline(y=0, color='black')
plt.axvline(x=0, color='black')
plt.xticks(fontsize=15, rotation=0)
plt.yticks(fontsize=15, rotation=0)
plt.ylim([-np.pi/2.0, np.pi/2.0])
plt.tight_layout()
plt.show()

plt.figure()
#plot the evolution of the LGN nematic angle alignment for each cell, as well as the mean
for this_cell in range(len(all_angles_ordered)):
    plt.plot(all_times[this_cell], a_vect[this_cell])
plt.errorbar(all_times[0], a_vect_mean, a_vect_std,marker='o', color='blue',
             lw=3,
             capsize=3,
             capthick=3)
plt.title('Alignment of LGN nematic angle with cell long axis \n as a function of time', fontsize=15)
plt.xlabel('Time after NEB (minutes)', fontsize=15)
plt.ylabel('alignment', fontsize=15)
plt.box(False)
plt.grid(True)
plt.axhline(y=0, color='black')
plt.axvline(x=0, color='black')
plt.xticks(fontsize=15, rotation=0)
plt.yticks(fontsize=15, rotation=0)
plt.ylim([-1.2, 1.2])
plt.tight_layout()
plt.show()


# %% Data preparation and plots for Figure 4K-M-O-Q and S4F

def load_data_angles(file_name):
    angle_data_control = []
    angle_data_LGN = []
    with open(file_name) as csvfile:
        has_header = csv.Sniffer().has_header(csvfile.read(1024))
        csvfile.seek(0)
        readCSV = csv.reader(csvfile, delimiter=',')
        if has_header:
            next(readCSV)  # Skip header row.
        cell_id = 1
        this_list_time = []
        this_list_angle = []
        this_label = ''
        counter = 0 #to detect the first line that is read
        for row in readCSV:
            if len(row) == 0: #Skip empty line
                continue
            if counter == 0:
                cell_id = int(row[0])
                this_label = str(row[1])
            if not int(row[0]) == cell_id: #when a new cell id is detected
                #save current list
                if this_label == 'micropattern_lgn_rnai':
                    angle_data_LGN.append([np.array(this_list_time),
                                           np.array(this_list_angle)])
                else:
                    angle_data_control.append([np.array(this_list_time),
                                               np.array(this_list_angle)])
                #create new list
                this_label = str(row[1])
                cell_id = int(row[0])
                this_list_time = []
                this_list_angle = []
            this_list_time.append(float(row[2]))
            #modify the angle to have it in radians
            this_list_angle.append(np.pi / 180.0 * float(row[3]))
            counter += 1
        #save final list at the end of the file
        if this_label == 'micropattern_lgn_rnai':
            angle_data_LGN.append([np.array(this_list_time),
                                    np.array(this_list_angle)])
        else:
            angle_data_control.append([np.array(this_list_time),
                                        np.array(this_list_angle)])
    return [angle_data_control, angle_data_LGN]
#### dictionaries that will contain data####
#### LOAD DATA  ####
path = main_path + '/data/data_rnai.csv'
[angle_data_control, angle_data_LGN] = load_data_angles(path)
#offset data by pi or -pi according to initial angle
#This ensures that all angles start in the interval between 0 and pi
for k in range(len(angle_data_control)):
    if angle_data_control[k][1][0] < 0:
        angle_data_control[k][1] = angle_data_control[k][1] + np.pi
    if angle_data_control[k][1][0] > np.pi:
        angle_data_control[k][1] = angle_data_control[k][1] - np.pi
for k in range(len(angle_data_LGN)):
    if angle_data_LGN[k][1][0] < 0:
        angle_data_LGN[k][1] = angle_data_LGN[k][1] + np.pi
    if angle_data_LGN[k][1][0] > np.pi:
        angle_data_LGN[k][1] = angle_data_LGN[k][1] - np.pi

# Plot of Figure 4K
plt.figure()
for k in range(len(angle_data_control)):
    plt.plot(angle_data_control[k][0], angle_data_control[k][1])
plt.title('Angles as a function of time (control)', fontsize=15)
plt.xlabel('Time after spindle formation (minutes)', fontsize=15)
plt.ylabel('Angle(rad)', fontsize=15)
plt.xlim([0, 21])
plt.ylim([-1.4, 4.5])
plt.box(False)
plt.grid(True)
plt.axhline(y=0, color='black')
plt.axvline(x=0, color='black')
plt.xticks(np.arange(0, 21, 5), fontsize=15, rotation=0)
plt.yticks(fontsize=15, rotation=0)
plt.tight_layout()
plt.show()

# Plot of Figure 4O
plt.figure()
for k in range(len(angle_data_LGN)):
    plt.plot(angle_data_LGN[k][0], angle_data_LGN[k][1])
plt.title('Angles as a function of time (LGN RNAi)', fontsize=15)
plt.xlabel('Time after spindle formation (minutes)', fontsize=15)
plt.ylabel('Angle(rad)', fontsize=15)
plt.box(False)
plt.grid(True)
plt.xlim([0, 21])
plt.ylim([-1.4, 4.5])
plt.axhline(y=0, color='black')
plt.axvline(x=0, color='black')
plt.xticks(np.arange(0, 21, 5), fontsize=15, rotation=0)
plt.yticks(fontsize=15, rotation=0)
plt.tight_layout()
plt.show()

# More data preparation

angle_data_control_pairs = []
for k in range(len(angle_data_control)):
    angle_data_control_pairs.extend([[this_time, this_angle]
                                    for this_time, this_angle
                                    in zip(angle_data_control[k][0], angle_data_control[k][1])])

angle_data_LGN_pairs = []
for k in range(len(angle_data_LGN)):
    angle_data_LGN_pairs.extend([[this_time, this_angle]
                                for this_time, this_angle
                                in zip(angle_data_LGN[k][0], angle_data_LGN[k][1])])


from collections import defaultdict
#first for control
angle_data_control_by_time = defaultdict(list)
for this_pair in angle_data_control_pairs:
    angle_data_control_by_time[int(this_pair[0])].append(this_pair[1])
#then for LGN
angle_data_LGN_by_time = defaultdict(list)
for this_pair in angle_data_LGN_pairs:
    angle_data_LGN_by_time[int(this_pair[0])].append(this_pair[1])

max_time_control = max([this_pair[0] for this_pair in angle_data_control_pairs])
time_array_control = [k for k in range(0, int(max_time_control), 3)]#since there is one data point every 3 minutes
alignment_mean_control = []
alignment_std_control = []
for this_time in time_array_control:
    alignment_mean_control.append(np.mean(np.cos(2 * np.array(angle_data_control_by_time[this_time]))))
    alignment_std_control.append(np.std(np.cos(2 * np.array(angle_data_control_by_time[this_time]))))

max_time_LGN = max([this_pair[0] for this_pair in angle_data_LGN_pairs])
time_array_LGN = [k for k in range(0, int(max_time_LGN), 3)]
alignment_mean_LGN = []
alignment_std_LGN = []
for this_time in time_array_LGN:
    alignment_mean_LGN.append(np.mean(np.cos(2 * np.array(angle_data_LGN_by_time[this_time]))))
    alignment_std_LGN.append(np.std(np.cos(2 * np.array(angle_data_LGN_by_time[this_time]))))


# Plot of Figure 4M
plt.figure()
plt.errorbar(time_array_control[0:8],
             alignment_mean_control[0:8],
             alignment_std_control[0:8],
             label='control',
             linewidth=4,
             capsize=6,
             capthick=4)
plt.xlim([-1,22])
plt.ylim([-1.2,1.2])
plt.title('Average alignment as a function of time', fontsize=15)
plt.xlabel('Time after spindle formation (minutes)', fontsize=15)
plt.ylabel('Alignment', fontsize=15)
plt.legend()
plt.box(False)
plt.grid(True)
plt.axhline(y=0, color='black')
plt.axvline(x=0, color='black')
plt.xticks(fontsize=15, rotation=0)
plt.yticks(np.arange(-1, 1.1, 0.5),fontsize=15, rotation=0)
plt.tight_layout()
plt.show()

# Plot of Figure 4Q
plt.figure()
plt.errorbar(time_array_LGN[0:8],
             alignment_mean_LGN[0:8],
             alignment_std_LGN[0:8],
             label='LGN RNAi',
             linewidth=4,
             capsize=6,
             capthick=4)
plt.xlim([-1,22])
plt.ylim([-1.2,1.2])
plt.title('Average alignment as a function of time', fontsize=15)
plt.xlabel('Time after spindle formation (minutes)', fontsize=15)
plt.ylabel('Alignment', fontsize=15)
plt.legend()
plt.box(False)
plt.grid(True)
plt.axhline(y=0, color='black')
plt.axvline(x=0, color='black')
plt.xticks(fontsize=15, rotation=0)
plt.yticks(np.arange(-1, 1.1, 0.5),fontsize=15, rotation=0)
plt.tight_layout()
plt.show()


# Plots of Figure S4F
plt.figure()
font_size = 20
plt.hist(np.cos(2*np.array(angle_data_control_by_time[0])), color='blue', bins=15,alpha=0.5, density=True, range=[-1,1])
plt.xlabel('alignment parameter', fontsize=font_size)
plt.ylabel('probability density', fontsize=font_size)
plt.xticks([-1, -0.5, 0, 0.5, 1],fontsize=font_size, rotation=0)
plt.yticks(fontsize=font_size, rotation=0)
plt.ylim([0,3.2])
plt.title('histogram of initial alignment parameter\n (t=0min, control) ', fontsize=font_size)
plt.tight_layout()
plt.show()

plt.figure()
plt.hist(np.cos(2*np.array(angle_data_control_by_time[21])), color='blue', bins=15,alpha=0.5, density=True, range=[-1, 1])
plt.xlabel('alignment parameter', fontsize=font_size)
plt.ylabel('probability density', fontsize=font_size)
plt.xticks([-1, -0.5, 0, 0.5, 1],fontsize=font_size, rotation=0)
plt.yticks(fontsize=font_size, rotation=0)
plt.title('histogram of final alignment parameter\n (t=21min, control) ', fontsize=font_size)
plt.ylim([0,3.2])
plt.tight_layout()
plt.show()

plt.figure()
plt.hist(np.cos(2*np.array(angle_data_LGN_by_time[0])), color='red', bins=15,alpha=0.5, density=True, range=[-1,1])
plt.xlabel('alignment parameter', fontsize=font_size)
plt.ylabel('probability density', fontsize=font_size)
plt.xticks([-1, -0.5, 0, 0.5, 1],fontsize=font_size, rotation=0)
plt.yticks(fontsize=font_size, rotation=0)
plt.title('histogram of initial alignment parameter\n (t=0min, LGN RNAi) ', fontsize=font_size)
plt.ylim([0,3.2])
plt.tight_layout()
plt.show()

plt.figure()
plt.hist(np.cos(2*np.array(angle_data_LGN_by_time[21])), color='red', bins=15,alpha=0.5, density=True, range=[-1, 1])
plt.xlabel('alignment parameter', fontsize=font_size)
plt.ylabel('probability density', fontsize=font_size)
plt.xticks([-1, -0.5, 0, 0.5, 1],fontsize=font_size, rotation=0)
plt.yticks(fontsize=font_size, rotation=0)
plt.title('histogram of  final alignment parameter\n (t=21min, LGN RNAi) ', fontsize=font_size)
plt.ylim([0,3.2])
plt.tight_layout()
plt.show()

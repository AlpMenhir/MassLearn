# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 10:57:53 2023

@author: Ronan
"""

import Modules.statplot as sp
import Modules.cleaning as cl
import matplotlib.pyplot as plt
from ipywidgets import interact_manual
from scipy.stats import zscore
import pandas as pd
import numpy as np
from fastdtw import fastdtw
from scipy.spatial.distance import euclidean
from scipy.optimize import curve_fit
from scipy.stats import norm, pearsonr



feat = sp.Filter('../feature/fl_ms1_h.csv')
fl = feat.featurelist
feat2 = sp.Filter('../feature/fl_ms2_h.csv')
fl2 = feat2.featurelist
# Create a new column 'ms_level' in df1 and set all values to 'ms1'
fl['ms_level'] = 'ms1'

# Create a new column 'ms_level' in df2 and set all values to 'ms2'
fl2['ms_level'] = 'ms2'

feature = pd.concat([fl, fl2], ignore_index=True)



class TrainingTool:
    def __init__(self, data, rt_tolerance=0.02):
        self.data = data
        self.rt_tolerance = rt_tolerance
        self.labels = [] # The output taking labels information
        self.spectra = {}
        self.metrics = [] # The output taking matching information (m/z diff, rt diff, similarity metric, euclidean distance)
        
        # Group the features based on ms_level and sample
        self.ms1_data = self.data[self.data['ms_level'] == 'ms1']
        self.ms2_data = self.data[self.data['ms_level'] == 'ms2']

        # Create a nested list to store MS1 features and their potential MS2 matches
        self.paired_features_id = [] # we indicate here the paired feature id, ta see how many times a ms2 is matching in the ms1's feature list
        self.paired_features = []
        for _, ms1_feature in self.ms1_data.iterrows():
            matching_ms2 = self.ms2_data[(self.ms2_data['sample'] == ms1_feature['sample']) & 
                                          (np.abs(self.ms2_data['rt'] - ms1_feature['rt']) <= self.rt_tolerance) &
                                         ((self.ms2_data['m/z'] + 1.007) < ms1_feature['m/z'])  &
                                         (self.ms2_data['height'] <= ms1_feature['height'])  ]
            self.paired_features.append((ms1_feature, matching_ms2))
            self.paired_features_id.append((ms1_feature['feature ID'].to_list()[0], matching_ms2['feature ID'].to_list()))
        self.current_index = [0, 0]  # first index for MS1, second index for MS2
        
        for s in data['sample'].unique(): #get samples names to take spectra from mzml files
            spec = cl.Spectra('../mzML/'+s+'.mzML')
            spec.extract_peaks()
            self.spectra[s] = spec
            
    # Function to adjust the length of the ms1 and ms2 feture intensities points for similarity metric and euclidean distance calculation
    def adjust_lengths(self, a, b):
        if len(a) > len(b):
            larger = a
            smaller = b
        else:
            larger = b
            smaller = a

        diff = len(larger) - len(smaller)
        for _ in range(diff):
            # Remove elements alternatively from the beginning (index 0) and end (index -1)
            if len(larger) % 2 == 0:
                larger = np.delete(larger, 0)
            else:
                larger = np.delete(larger, -1)
        return larger, smaller
    
    # Obtain here the different information for the decision tree
    def get_feature_vector(self, ms1_feature, ms2_feature):
        ms1_shape = self.spectra_get_ms1_intensities(ms1_feature), ms1_feature['m/z']
        ms2_shape = self.spectra_get_ms2_intensities(ms2_feature), ms2_feature['m/z']
        
        ms1_i = [i[0] for i in ms1_shape[0]]
        ms2_i = [i[0] for i in ms2_shape[0]]
        
        ms1_i_ad, ms2_i_ad = self.adjust_lengths(ms1_i, ms2_i)
        
        similarity_metric = pearsonr(zscore(ms1_i_ad), zscore(ms2_i_ad))[0]
        euclidean_distance = np.linalg.norm(zscore(ms1_i_ad) - zscore(ms2_i_ad))
        
        mz_diff = ms1_feature['m/z'] - ms2_feature['m/z']
        rt_diff = ms1_feature['rt'] - ms2_feature['rt']
        
        return [mz_diff, rt_diff, similarity_metric, euclidean_distance]

    def prepare_training_data(self):
        self.display_next_pair()
        X = self.metrics
        y = self.labels
        print(f"Size of X: {len(X)}, Size of y: {len(y)}")  # Print the sizes of X and y

        return np.array(X), np.array(y)
    
    def plot_features(self, ms1_shape, ms2_shape):
        ms1_i = [i[0] for i in ms1_shape[0]]
        ms1_rt = [i[1] for i in ms1_shape[0]]
        ms2_i = [i[0] for i in ms2_shape[0]]
        ms2_rt = [i[1] for i in ms2_shape[0]]
        
        ms1_i_ad, ms2_i_ad = self.adjust_lengths(ms1_i, ms2_i)      
        similarity_metric = pearsonr(zscore(ms1_i_ad), zscore(ms2_i_ad))[0]
        euclidean_distance = np.linalg.norm(zscore(ms1_i_ad) - zscore(ms2_i_ad))
        
        fig, axs = plt.subplots(1, 2, figsize=(20, 5))

        # Plot absolute intensities
        axs[0].plot(ms1_rt, ms1_i, label=f'MS1 - {ms1_shape[1]} Da', color='blue')
        axs[0].plot(ms2_rt, ms2_i, label=f'MS2 - {ms2_shape[1]} Da', color='red')
        axs[0].set_xlabel('rt (min)')
        axs[0].set_ylabel('Absolute Intensity')
        axs[0].legend()

        # Plot relative intensities (Z-score normalized)
        axs[1].plot(ms1_rt, zscore(ms1_i), label=f'MS1 - {ms1_shape[1]} Da', color='blue')
        axs[1].plot(ms2_rt, zscore(ms2_i), label=f'MS2 - {ms2_shape[1]} Da', color='red')
        axs[1].set_xlabel('rt (min)')
        axs[1].set_ylabel('Relative Intensity (Z-score)')
        axs[1].set_title(f'{round(similarity_metric, 3)} - {euclidean_distance}')
        axs[1].legend()

        plt.show()
 

    # Function to handle the button click
    def on_button_click(self, label):
    # Before appending the label, also append the current feature values
        self.labels.append(label)  # Map 'Yes' to 1, 'No' to 0, and 'Skip' to None
        self.display_next_pair()   
    

    # Function to obtain intensities for plotting features
    def spectra_get_ms1_intensities(self, ms_feature):
        rt_shift = 0.01 #this is to inrease the observation window of the features
        sample = ms_feature['sample']
        rt_min = float(ms_feature['peak_rt_start'])
        rt_max = float(ms_feature['peak_rt_end'])
        mz_tolerance = 0.002
        mz_min = float(ms_feature['peak_mz_min'])
        mz_max = float(ms_feature['peak_mz_max'])
        spec = self.spectra[sample]       
        scans = [] # scans are minimal and maximal scans index based of the rt range
        for rt_value  in [rt_min - rt_shift, rt_max + rt_shift]:
            match = 100 # match is 100min value to begin with a very high impossible rt value
            for idx, rt in enumerate(spec.rt1): # the idea is to have the lowest rt difference find the corresponding rt value of a scan
                difference = float(rt_value) - float(rt)  # take the difference of rt
                if abs(difference) < abs(match):
                    match = difference
                else: # if the difference increases rather than decreasing we have found the right rt, because the loop is in increase rt duration order
                    scans.append(idx-1)
                    break              
        intensities = []
        for scan in range(scans[0], scans[1]):
            ms_arr = spec.peakarray1[scan]
            df = pd.DataFrame(ms_arr, columns=['mz', 'int'])
            matches = df.loc[(df['mz'] >= (mz_min - mz_tolerance)) & (df['mz'] <= (mz_max + mz_tolerance))] # Perform a match within the m/z range
            if matches.empty == False:
                intensities.append((max(matches['int']), spec.rt1[scan]))
        return intensities
                 
    
    # Function to obtain intensities for plotting features, ms2 (only spec.rt2 changes)
    def spectra_get_ms2_intensities(self, ms_feature):
        rt_shift = 0.01 #this is to inrease the observation window of the features
        sample = ms_feature['sample']
        rt_min = float(ms_feature['peak_rt_start'])
        rt_max = float(ms_feature['peak_rt_end'])
        mz_tolerance = 0.002
        mz_min = float(ms_feature['peak_mz_min'])
        mz_max = float(ms_feature['peak_mz_max'])
        spec = self.spectra[sample]       
        scans = [] # scans are minimal and maximal scans index based of the rt range
        for rt_value  in [rt_min - rt_shift, rt_max + rt_shift]:
            match = 100 # match is 100min value to begin with a very high impossible rt value
            for idx, rt in enumerate(spec.rt2): # the idea is to have the lowest rt difference find the corresponding rt value of a scan
                difference = float(rt_value) - float(rt)  # take the difference of rt
                if abs(difference) < abs(match):
                    match = difference
                else: # if the difference increases rather than decreasing we have found the right rt, because the loop is in increase rt duration order
                    scans.append(idx-1)
                    break              
        intensities = []
        for scan in range(scans[0], scans[1]):
            ms_arr = spec.peakarray1[scan]
            df = pd.DataFrame(ms_arr, columns=['mz', 'int'])
            matches = df.loc[(df['mz'] >= (mz_min - mz_tolerance)) & (df['mz'] <= (mz_max + mz_tolerance))] # Perform a match within the m/z range
            if matches.empty == False:
                intensities.append((max(matches['int']), spec.rt2[scan]))
        return intensities

    # Function to display the next pair
    def display_next_pair(self):
        if self.current_index[0] >= len(self.paired_features):
            print("All pairs have been labeled.")
            return

        ms1_feature, matching_ms2 = self.paired_features[self.current_index[0]] # matching ms2 are all the ms2 feature matching the ms1 feature

        if self.current_index[1] >= len(matching_ms2): # if no matching ms2, len will be 0 and it will be skept
            self.current_index =  [self.current_index[0] + 1, 0]
            self.display_next_pair()
            return

        current_pair_ms2 = matching_ms2.iloc[self.current_index[1]] # take the current ms2 feature
        self.current_index[1] += 1
        
        ms1_shape = self.spectra_get_ms1_intensities(ms1_feature), ms1_feature['m/z']   # Fetch shape data for MS1 feature
        ms2_shape = self.spectra_get_ms2_intensities(current_pair_ms2), current_pair_ms2['m/z']   # Fetch shape data for MS2 feature
        self.metrics.append(self.get_feature_vector)
        self.plot_features(ms1_shape, ms2_shape)
        interact_manual(self.on_button_click, label=['Yes', 'No', 'Skip'])
             
   
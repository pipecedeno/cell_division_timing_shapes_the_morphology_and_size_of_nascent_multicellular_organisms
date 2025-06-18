#!/usr/bin/env python3

import pandas as pd
import numpy as np
from tqdm import tqdm
import os
import warnings
import argparse
warnings.filterwarnings('ignore')

def calc_avg_pairwise_dist(x, y, z, n_pairs):
    """Calculate average pairwise distance between points using random sampling."""
    points = np.column_stack((x, y, z))
    n_points = len(points)
    
    if n_points < 2:
        return np.nan
    
    # Calculate how many pairs we'll actually sample (minimum of possible pairs and requested pairs)
    max_possible_pairs = (n_points * (n_points - 1)) // 2
    n_pairs_to_sample = min(n_pairs, max_possible_pairs)
    
    # Generate random pairs of indices
    idx1 = np.random.randint(0, n_points, n_pairs_to_sample)
    idx2 = np.random.randint(0, n_points, n_pairs_to_sample)
    
    # Ensure we don't calculate distance of a point with itself
    same_point = idx1 == idx2
    idx2[same_point] = (idx2[same_point] + 1) % n_points
    
    # Calculate distances for sampled pairs
    diff = points[idx1] - points[idx2]
    distances = np.sqrt(np.sum(diff**2, axis=1))
    
    return([np.mean(distances), np.std(distances)])

def calc_hopkins_statistic(x, y, z, n_samples=None):
    """Calculate Hopkins statistic for clustering tendency."""
    points = np.column_stack((x, y, z))
    n_points = len(points)
    
    if n_points < 5:
        return np.nan
    
    if n_samples is None:
        n_samples = min(n_points, 50)
    
    # Create bounding box
    ranges = np.vstack([np.min(points, axis=0), np.max(points, axis=0)])
    
    # Generate random points
    random_points = np.random.uniform(
        ranges[0], 
        ranges[1], 
        size=(n_samples, 3)
    )
    
    # Sample points from dataset
    sampled_indices = np.random.choice(n_points, n_samples, replace=False)
    sampled_points = points[sampled_indices]
    
    # Calculate minimum distances for random points
    u_distances = np.zeros(n_samples)
    for i in range(n_samples):
        distances = np.sqrt(np.sum((points - random_points[i])**2, axis=1))
        u_distances[i] = np.min(distances)
    
    # Calculate minimum distances for sampled points
    w_distances = np.zeros(n_samples)
    for i in range(n_samples):
        points_without_sample = np.delete(points, sampled_indices[i], axis=0)
        distances = np.sqrt(np.sum((points_without_sample - sampled_points[i])**2, axis=1))
        w_distances[i] = np.min(distances)
    
    hopkins = np.sum(u_distances) / (np.sum(u_distances) + np.sum(w_distances))
    
    return hopkins

def process_group(group_info, group, n_pairs):
    """Process a single group of data."""
    strain, sim_number, file_num = group_info
    
    try:
        avg_dist, sd_dist = calc_avg_pairwise_dist(
            group['overlap_x'],
            group['overlap_y'],
            group['overlap_z'],
            n_pairs
        )
        
        hopkins = calc_hopkins_statistic(
            group['overlap_x'],
            group['overlap_y'],
            group['overlap_z']
        )
        
        return {
            'strain': strain,
            'sim_number': sim_number,
            'file_num': file_num,
            'avg_pairwise_dist': avg_dist,
            'sd_pairwise_dist': sd_dist,
            'hopkins_stat': hopkins,
            'n_points': len(group)
        }
    except Exception as e:
        print(f"Error processing group {strain}-{sim_number}-{file_num}: {str(e)}")
        return None

# Set up argument parser
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_folder", required=True, help="Folder containing input files")
parser.add_argument("-f", "--file_suffix", required=True, help="Suffix for files of the overlap positions")
parser.add_argument("-n", "--n_pairs", type=int, default=1000, help="Number of pairs to sample for distance calculation")
parser.add_argument("-o", "--output_file", required=True, help="Output file name")
args = parser.parse_args()

# Create empty dataframe
overlap_pos_df = pd.DataFrame()

# Read and combine files
print("\nReading files...")
for strain in ['petite', 'grande']:
    filename = os.path.join(args.input_folder, f"{strain}{args.file_suffix}")
    if os.path.exists(filename):
        temp_df = pd.read_csv(filename)
        temp_df['strain'] = 'Petite' if strain == 'petite' else 'Grande'
        overlap_pos_df = pd.concat([overlap_pos_df, temp_df])
    else:
        print(f"Warning: File not found - {filename}")

# Convert strain to categorical
overlap_pos_df['strain'] = pd.Categorical(
    overlap_pos_df['strain'], 
    categories=['Petite', 'Grande']
)

# Get groups for processing
groups = list(overlap_pos_df.groupby(['strain', 'sim_number', 'file_num']))
total_groups = len(groups)

print(f"\nProcessing {total_groups} groups sequentially...")
print(f"Sampling {args.n_pairs} pairs per group for distance calculation")

# Process groups sequentially with progress bar
results = []
for group_info, group_data in tqdm(groups, total=total_groups, desc="Calculating statistics"):
    result = process_group(group_info, group_data, args.n_pairs)
    if result is not None:
        results.append(result)

# Convert results to dataframe
summ_overlap_pos = pd.DataFrame(results)

# Save results
print("\nSaving results...")
output_path = os.path.join(args.input_folder, args.output_file)
summ_overlap_pos.to_csv(output_path, index=False)
print(f"Results saved to {output_path}")



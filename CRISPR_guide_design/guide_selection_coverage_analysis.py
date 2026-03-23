"""
Guide selection and contig coverage analysis utilities.
"""

# ============================================================
# Imports
# ============================================================
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# ============================================================
# 1) Guide selection across sliding windows
# ============================================================
# function to choose the best guides for each bin (array)
def select_sliding_window_guides(df, bin_number, window_size=100):

    # Filter DataFrame for the given bin
    df_bin = df[df['Bin'] == bin_number]

    # Calculate the global scoring threshold for the given bin (the guide needs to atleast pass this threshold)
    threshold = round(df_bin[df_bin['score']>0]['score'].mean(),3)
    out_rows = []
    # Iterate over each unique contig
    for contig in df_bin['contig'].unique():
        contig_df = df_bin[df_bin['contig'] == contig]

        # Determine the range to scan (based on start and stop columns)
        min_start = contig_df['start'].min()
        max_stop = contig_df['stop'].max()

        # Scan every 100bp window
        for start in range(min_start, max_stop, window_size):
            window_end = start + window_size
            window_df = contig_df[(contig_df['start'] < window_end) & (contig_df['stop'] > start)]

            if not window_df.empty:
                # Apply the decision tree
                if (window_df['score'] > threshold).any():
                    # Select the guide with the highest otCount among those with score > threshold
                    best_guide = window_df[window_df['score'] > threshold].sort_values(by='otCount', ascending=False).iloc[0]
                else:
                    # Select the guide with the highest score
                    best_guide = window_df.sort_values(by='score', ascending=False).iloc[0]

                # Add the selected guide RNA info to the results DataFra
                out_rows.append({
                'bin': bin_number,
                'contig': best_guide['contig'],
                'start': int(best_guide['start']),
                'stop': int(best_guide['stop']),
                'target': best_guide.get('target', None),
                'score': float(best_guide['score']),
                'otCount': int(best_guide['otCount']),
                'window': f"{start}-{window_end}",
            })

    return pd.DataFrame(out_rows)


# ============================================================
# 2) Contig length and guide density calculations
# ============================================================
# calculate some metrics on the coverage over the contigs 
# this enables you to select for contigs that 'dropped out' of the guide design

# Function to calculate contig length
def calculate_contig_length(contig):
    parts = contig.split(':')
    start_stop = parts[1].split('-')
    return int(start_stop[1]) - int(start_stop[0])

# Calculate the density of guides per contig
def calculate_guide_density_over_contig(df):
    # Calculate contig lengths
    df['contig_length'] = df['contig'].apply(calculate_contig_length)

    # Group by contig and calculate density
    density_df = df.groupby('contig').apply(lambda x: len(x) / (x['contig_length'].iloc[0] / 100)).reset_index()
    density_df.columns = ['contig', 'guide_density']

    return density_df


# ============================================================
# 3) Guide coverage plotting
# ============================================================
# function to plot the guide coverage distribution, input is the df and title of the plot
def plot_guide_coverage_distribution(df, title, save=False, filename=None, color='blue'):
    # Plot the distribution of guide coverage
    plt.figure(figsize=(6, 6), dpi=100)
    sns.histplot(data=df, x="guide_density", bins=20, color=color)
    plt.title(title, fontsize=18)
    plt.ylabel('Density', fontsize=14)
    plt.xlabel('Guide Density over Contig', fontsize=14)

    # Calculate mean and its annotation position
    mean_guide_density = df['guide_density'].mean()
    plt.axvline(x=mean_guide_density, color='black', linestyle='--')

    # Get current axis limits
    ax = plt.gca()
    y_lim = ax.get_ylim()

    # Calculate position for the text (e.g., top left corner)
    text_x = mean_guide_density + (ax.get_xlim()[1] - mean_guide_density) * 0.01  # Slightly right from the mean line
    text_y = y_lim[1] * 0.95  # 95% of the way up the y-axis

    plt.text(text_x, text_y, f"Mean = {mean_guide_density:.3f}", verticalalignment='top')

    plt.tight_layout()
    if save:
        plt.savefig(filename, dpi=300)


# ============================================================
# 4) Supplement selected guides with repeat-hit guides
# ============================================================
# function that takes in the original df scored and the array_1 df and the selected bin and takes the array_1 and appends to the dataframe the guides from the scored df in the same bin that have the highest otCount, to get the amount of total guides in array1 to 60,000
def supplement_guides_with_repeat_hits(scored_df, selected_df, bin_number, total_guides=60000):
    scored_df_bin = scored_df[scored_df['Bin'] == bin_number]

    # Calculate the number of guides to select from scored_df, 60k is the max per array 
    num_guides_to_select = total_guides - len(selected_df)

    # Select the guides from scored_df (no scoring threshold)
    selected_guides = scored_df_bin.sort_values(by='otCount', ascending=False).iloc[:num_guides_to_select]
    selected_guides['window'] = 'repeat_guide'
    selected_guides['bin'] = bin_number
    # calculate contig length from contig information by subtracting the two numbers i.e. 'chr1:161540454-161552547' would be 161552547-161540454
    selected_guides['contig_length'] = selected_guides['contig'].apply(calculate_contig_length)

    # Add the selected guides to selected_df
    selected_df = pd.concat([selected_df, selected_guides[['bin', 'contig', 'start', 'stop', 'target', 'score', 'otCount','window']]], ignore_index=True)

    return selected_df

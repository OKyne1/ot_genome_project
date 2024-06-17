import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np

def process_file(file_path):
    # Read the file into a pandas DataFrame
    df = pd.read_csv(file_path, sep='\t', header=None)
    
    # Preprocess column 8 (index 7)
    df[7] = df[7].apply(lambda x: np.nan if pd.isna(x) else (np.mean(list(map(float, x.split(',')))) if isinstance(x, str) else np.nan))
    
    # Scatter plot for columns 9 (index 8) vs 10 (index 9)
    plt.figure(figsize=(10, 6))
    plt.scatter(df.iloc[:, 8], df.iloc[:, 9], alpha=0.5)
    plt.xlabel('Alternative Read Count')
    plt.ylabel('Read Depth')
    plt.title(f'{os.path.basename(file_path)} - AO vs DP')
    scatter_ao_vs_dp_filename = f"{os.path.basename(file_path).split('.')[0]}_AOvsDP.png"
    plt.savefig(scatter_ao_vs_dp_filename)
    plt.close()

    # Scatter plot for columns 8 (index 7) vs 10 (index 9)
    plt.figure(figsize=(10, 6))
    plt.scatter(df.iloc[:, 7], df.iloc[:, 9], alpha=0.5)
    plt.xlabel('snp ratio')
    plt.ylabel('Read Depth')
    plt.title(f'{os.path.basename(file_path)} - AB vs DP')
    scatter_ab_vs_dp_filename = f"{os.path.basename(file_path).split('.')[0]}_ABvsDP.png"
    plt.savefig(scatter_ab_vs_dp_filename)
    plt.close()
    
    return df.iloc[:, 7]


def main(directory):
    all_data = []
    file_names = []

    for file in os.listdir(directory):
        if os.path.isfile(os.path.join(directory, file)):
            file_path = os.path.join(directory, file)
            data = process_file(file_path)
            all_data.append(data)
            file_names.append(os.path.basename(file).split('.')[0])

    # Violin plot for column 7 (index 6) for each file separately
    plt.figure(figsize=(12, 8))
    sns.violinplot(data=all_data, inner=None)
    plt.xticks(range(len(file_names)), file_names, rotation=45)
    plt.xlabel('Sample')
    plt.ylabel('snp ratio')
    plt.title('snp ratio for each sample (+_ RAGE regions)')
    violin_plot_filename = "violin_plot.png"
    plt.savefig(violin_plot_filename)
    plt.close()



if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python plot_vcf.py <directory_path>")
        sys.exit(1)
    
    directory = sys.argv[1]
    main(directory)

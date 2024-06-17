import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def generate_swarm_plot(files):
    data = []
    labels = []

    # Iterate over the provided files
    for file_path in files:
        # Read the file into a DataFrame
        df = pd.read_csv(file_path, sep='\t', header=None)
        
        # Extract values from column 8 (index 7) and append to data list
        if not df.empty and len(df.columns) > 7:
            column_data = df.iloc[:, 7].astype(float).tolist()
            data.extend(column_data)
            labels.extend([file_path.split('/')[-1].split('.')[0]] * len(column_data))

    # Create a swarm plot
    plt.figure(figsize=(10, 6))
    sns.swarmplot(x=labels, y=data)
    plt.xlabel('File')
    plt.ylabel('Value in Column 8')
    plt.title('Swarm Plot')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # Check if file paths are provided as command-line arguments
    if len(sys.argv) < 2:
        print("Usage: python script.py file1.txt file2.txt ...")
        sys.exit(1)

    # Extract file paths from command-line arguments
    files = sys.argv[1:]

    # Generate swarm plot
    generate_swarm_plot(files)

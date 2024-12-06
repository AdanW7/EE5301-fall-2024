#   Adan Wodzinski

#   the following program was generated mostly via chatgpt 4-o 
#   with the following prompt:
#   give me a python program that reads in a csv that has three fields"type, x coord, ycoord" and make a graph representation of it


#pip install matplotlib networkx pandas

import matplotlib.pyplot
import pandas
import sys

# Get circuit name from command line argument
if len(sys.argv) > 1:
    circuitName = sys.argv[1]
else:
    circuitName = ""

# Define file names
pre_spread = f"pre_spread_{circuitName}.csv"
post_spread = f"post_spread_{circuitName}.csv"

# Function to read and process CSV file
def process_csv(file_name):
    try:
        df = pandas.read_csv(file_name)
    except FileNotFoundError:
        print(f"The file {file_name} was not found.")
        exit()

    # Trim spaces from column names
    df.columns = df.columns.str.strip()

    # Check for required columns
    if not all(col in df.columns for col in ['type', 'x coord', 'y coord']):
        print("The CSV file does not contain the required columns: 'type', 'x coord', 'y coord'.")
        print(f"Columns found: {df.columns.tolist()}")
        exit()

    # Map old type names to new ones
    type_mapping = {'cell': 'cell', 'star': 'star', 'pin': 'pin'}
    df['type'] = df['type'].map(type_mapping)

    return df

# Function to plot data from CSV
def plot_data(df, file_name):
    # Define colors and markers for each type
    colors = {'cell': 'red', 'star': 'blue', 'pin': 'green'}
    markers = {'cell': 'o', 'star': 's', 'pin': '^'}

    # Create the plot
    matplotlib.pyplot.figure(figsize=(10, 10))

    # Plot each type with different color and transparency
    for cell_type in df['type'].unique():
        subset = df[df['type'] == cell_type]
        matplotlib.pyplot.scatter(subset['x coord'], subset['y coord'], 
                                 c=colors[cell_type], label=cell_type, alpha=0.5, marker=markers[cell_type])

    # Add legend and labels
    matplotlib.pyplot.legend()
    matplotlib.pyplot.xlabel('X Coord')
    matplotlib.pyplot.ylabel('Y Coord')
    matplotlib.pyplot.title('Locations of Cells and Pins')

    # Save the plot as a PNG file
    matplotlib.pyplot.savefig(file_name)

    # Show the plot (optional)
    # matplotlib.pyplot.show()

# Process pre_spread file
df_pre_spread = process_csv(pre_spread)
plot_data(df_pre_spread, f"pre_spread_cell_pin_{circuitName}.png")

# Process post_spread file
df_post_spread = process_csv(post_spread)
plot_data(df_post_spread, f"post_spread_cell_pin_{circuitName}.png")

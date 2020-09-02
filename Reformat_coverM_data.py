# To use this script, run python Reformat_coverM_data.py Sample-Participant-List MAG-key CoverM_output_directory

# imports --------------------------------------------------------------------------------------------------------------

import pandas as pd
import sys

# set directory to save converted data-sheet ---------------------------------------------------------------------------

output = "/srv/projects3/human_plasmids/georgina/16_Host_Linkage/Test2"

# Import CoverM data (combine multiple runs if necessary), participant to sample and MAG key ---------------------------

print("Importing data...")

p_list = pd.read_csv(open(sys.argv[1], 'r'))
c_data = pd.read_csv(open(sys.argv[3], 'r'),index_col=0)
MAG_key = pd.read_csv(open(sys.argv[2], 'r'))

participants = p_list.Participant.unique()

# Reformat the dataset -------------------------------------------------------------------------------------------------

c_data = c_data.drop(["measure"], axis = 1)

c_data = c_data[["Sample", "Genome", "value"]]

df = c_data.merge(p_list, how='left', on='Sample')

df = df.merge(MAG_key, how='left', on='Genome')

df['Species'].fillna(df['Genome'], inplace = True)

df = df.drop(["Genome"], axis = 1)


# Sum the equivalent MAG values ----------------------------------------------------------------------------------------

df_grouped = df.groupby(['Sample','Species'])["value"].sum().reset_index()

df_grouped.to_csv(output+".csv", index = False, header = True)
# This script returns Clusters that have certain KOs from EnrichM output
# To use this script, run python KO-filter.py EnrichM_KO_Out KO_List


# Load the EnrichM output and the KO list ------------------------------------------------------------------------------

import sys
import pandas as pd

# set directory to save individual participant and plasmid CSVs --------------------------------------------------------

directory = "/srv/projects3/human_plasmids/georgina/9_EnrichM/AMG_Fishing/"

print("Opening data files")

enrichm_ko = pd.read_csv(sys.argv[1], sep='\t', index_col=0)

f = open(sys.argv[2], 'r')
ko_list = [line.rstrip('\n') for line in f]
print("Retrieving {} KOs".format(len(ko_list)))

# Reorganise the dataframe ---------------------------------------------------------------------------------------------

t_df = enrichm_ko.transpose()

t_df2 = t_df[t_df.columns.intersection(ko_list)]

t_df2["Sum"] = t_df2.sum(axis=1)

t_df3 = t_df2[t_df2["Sum"] != 0]

t_df3.to_csv(directory+"AMG_fish_enrichM_onlyARG.csv", index = True, header = True)

print("Complete!")


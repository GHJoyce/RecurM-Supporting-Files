
# To use this script, run python Find-Hosts.py CLR-dataset Sample-Participant-List Plasmid-List

import sys
import pandas as pd
import csv
import numpy
import os

# set directory to save individual participant and plasmid CSVs
directory = "/srv/projects3/human_plasmids/georgina/16_Host_Linkage/4_Propr/output_27-08-20/"

# Proportionality (Rho) equation --------------------------------
# a and b are lists of values (columns from CLR dataframe)

def rho(a, b):
    c = []
    for (i, item) in enumerate(a):
        c.append(item-b[i])
    c_var = numpy.var(c)
    prop = 1 - (c_var/(numpy.var(a) + numpy.var(b)))
    return prop

def func(a):
    x = len(a)
    y = []
    for value in a:
        if value < 0:
            y.append(value)
    z = len(y)/x
    return z

# import CLR data (out from PROPR) and participant to sample key --
print("Opening data files")

clr_df = pd.read_csv(open(sys.argv[1], 'r'))
p_list = pd.read_csv(open(sys.argv[2], 'r'))

participants = p_list.Participant.unique()

print("{} participants to analyse".format(len(participants)))

# List of the plasmids of interest: ------------------------

f3 = open(sys.argv[3], 'r')
plasmids = [line.rstrip('\n') for line in f3]
print("Plasmid list has {} plasmids".format(len(plasmids)))

# Set up the output spreadsheet -----------------------------------

out = open(directory+'Find-Hosts-Out.csv', 'w',  newline='')
w = csv.writer(out)
w.writerow(['Participant','Plasmid','MAG', 'Rho'])

# Clean up CLR data
melt_df = pd.melt(clr_df, id_vars=['Sample'], var_name='Genome', value_name='Value')
print("Melted dataframe.")

# Add participants to dataframe
merge_df = melt_df.merge(p_list, how='left', on='Sample')
print("Merged CLR data with Participant List.")

# Split by participant --------------------------------------------

for p in participants:

    # Subset a df that has one participant and the CLR values for the mags
    df = merge_df[~merge_df.Genome.isin(plasmids) & merge_df.Participant.eq(p)]
    print("Dataframe contains {} samples for participant {}.".format(len(df.Sample.unique()), p))
    mag_list = df.Genome.unique()

    if len(df.Sample.unique()) > 3:

        if not os.path.exists(directory+p):
            os.mkdir(directory + p)

        for plasmid in plasmids:
            # Subset another df that's the CLR values for one plasmid
            p_df = merge_df[merge_df.Genome.eq(plasmid) & merge_df.Participant.eq(p)]

            p_df = p_df.rename(columns = {'Value': plasmid})

            p_df = p_df.drop('Genome', 1)

            for mag in mag_list:
                mag_df = df[df.Genome.eq(mag)]
                mag_df = mag_df.rename(columns={'Value': mag})
                mag_df = mag_df.drop(['Genome', 'Participant'], 1)

                p_df1 = p_df.merge(mag_df, how = 'left', on = 'Sample')

                mag_clr = p_df1[mag].values.tolist()
                p_clr = p_df1[plasmid].values.tolist()

                if p_clr:
                    if func(mag_clr) < 0.4 and func(p_clr) < 0.4:

                        calc = rho(mag_clr, p_clr)
                        if calc > 0.95:
                            w.writerow([p, plasmid, mag, calc])
                            # export dataframe
                            p_df1.to_csv(directory+p+"/CLR_{}_{}_{}.csv".format(p, plasmid, mag), index = False, header = True)

print("Complete!")
from Bio import SeqIO

list_file = "/srv/projects3/human_plasmids/georgina/16_Host_Linkage/1_Blast/Outv2/MAGS_plasmidA_list.txt"
# open list text file
f = open(list_file, 'r')
# set mag directory
mag_dir = "/srv/projects3/human_plasmids/georgina/7_coverm/reformatted_mags/"
# set cluster directory
n = 0
for line in f.readlines():
    subject = line.split("\t")[1].strip("\n\t")
    print("Subject:", subject)
    assemblyname = "f."+subject.split("_")[0]
    print("Assembly:", assemblyname)
    assembly = list(SeqIO.parse("{}/{}".format(mag_dir, assemblyname), 'fasta'))
    print("Opening Assembly {} of length {}".format(assemblyname, len(assembly)))
    new = []
    for seq_record in assembly:
        if seq_record.id != subject:
            new.append(seq_record)
        elif seq_record.id == subject:
            n = n+1
    SeqIO.write(new, "{}/{}".format(mag_dir, assemblyname), 'fasta')
print("Done! Deleted ", n, " contigs.")




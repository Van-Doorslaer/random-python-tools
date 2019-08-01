from Bio import SeqIO
import glob, re, csv, os, sys
from ete3 import Tree

def search_by_size(node, size):
    "Finds nodes with a given number of leaves"
    matches = []
    for n in node.traverse():
       if len(n) == size:
          matches.append(n)
    return matches

def reverse_translate(nucleotide_fasta):
    with open(nucleotide_fasta+".aa.fas","w") as out:
        for record in SeqIO.parse(nucleotide_fasta,"fasta"):
            print >>out, ">"+record.id+"\n"+record.seq.translate()
    os.system( str( "mafft --quiet "+nucleotide_fasta+".aa.fas > "+nucleotide_fasta+".aa.mafft.fas" )) #this can be replaced with any aligner
    with open(nucleotide_fasta+".mafft.fas","w") as out:
        for aa_rec in SeqIO.parse(nucleotide_fasta+".aa.mafft.fas","fasta"):
            for nuc_rec in SeqIO.parse(nucleotide_fasta,"fasta"):
                sequence=[]
                if aa_rec.id == nuc_rec.id:
                    codons = [str(nuc_rec.seq)[i:i+3] for i in range(0, len(str(nuc_rec.seq)), 3)]
                    x=0
                    for aa in list(aa_rec.seq):
                        if aa != "-":
                            sequence.append(codons[x])
                            x=x+1
                        else:
                            sequence.append("---")
                    #print nuc_rec.id, x, len(str(nuc_rec.seq))/3 #">"+nuc_rec.id+"\n"+"".join(sequence)
                    print >>out, ">"+nuc_rec.id+"\n"+"".join(sequence)

def align_2_promals(existing_alignment, csv_file):
        gene = existing_alignment.split(".")[0][-2:]
        print gene
        for r in csv.reader(open(csv_file,"rU")):
            if r[1] == gene:
                if gene == "L1":
                    outfile=r[0]+"_"+gene+".aa.fas"
                    with open(r[0]+"_"+gene+".aa.fas","w") as out:
                        print >>out, ">%s\n%s" %(r[0], r[4])
                    with open(r[0]+"_"+gene+".fas","w") as out:
                        print >>out, ">%s\n%s" %(r[0], r[3])
                    string = "mafft --add %s --reorder %s  > %s" %(outfile,existing_alignment, ".".join(outfile.split(".")[:-1]+["mafft","fas"]))
                    os.system (string)  
                else:
                    outfile=r[0]+"_"+gene+".aa.fas"
                    with open(outfile,"w") as out:
                        print >>out, ">%s\n%s" %(r[0], r[4])
                    string = "mafft --add %s --reorder %s  > %s" %(outfile,existing_alignment, ".".join(outfile.split(".")[:-1]+["mafft","fas"]))
                    os.system (string)

def trimal_gappy(mafft_aligned):
    for aligned in mafft_aligned:
        string = "trimal -in %s -out %s -gappyout" %(aligned, ".".join(aligned.split(".")[:-1]+["gappy","fas"]))
        os.system (string)
        
def concatenate(sequences):  
    viruses=[]
    for r in SeqIO.parse(sequences[0],"fasta"):
        viruses.append( r.id )
    
    with open("E1-E2-L1.fas","w") as out:    
        for v in viruses:
            temp=[]
            for sequence in sequences:
                for r in SeqIO.parse(sequence,"fasta"):
                    if r.id == v:
                        temp.append(str(r.seq))
            print >>out, ">%s\n%s" %(v, "".join(temp))

def find_closest_neighbor(tree, new_virus, L1db):
    with open(tree,"rU") as f:
        for l in f:
            l=l.strip()
            t=Tree(l)
            for clade in search_by_size(t, size=2):
                if new_virus in clade:
                    for leaf in clade:
                        if new_virus not in leaf:
                            neighbor = str(leaf)[3:]
    for r in SeqIO.parse(L1db,"fasta"):
        if r.id == neighbor:
            with open("pairwise.fas","w") as out:
                print >>out, ">%s\n%s" %(r.id, r.seq)
                for rec in SeqIO.parse(new_virus+"_L1.fas","fasta"):
                    print >>out, ">%s\n%s" %(rec.id, rec.seq)

def pairwise_identity(alignment):
    pair=[]
    viruses = []
    identity = 0
    for r in SeqIO.parse(alignment,"fasta"):
        pair.append (list(r.seq.upper()))
        viruses .append(r.id)
    for r in range(len( pair[0] )):
        if pair[0][r] == pair[1][r]:
            identity = identity +1
    percent = 100*(float(identity)/float(len(pair[0])))
    return viruses+[percent]

candidates=[]
PuMA_results = sys.argv[1]
for r in csv.reader(open(PuMA_results,"rU")):
    candidates.append(r[0])
for new_virus in set(candidates):
    #align new sequences to a promals generated guide alignment
    print "#####################################################################################"
    print "aligning %s E1, E2, L1 to existing database" %(new_virus)
    print "#####################################################################################"
    for file in os.listdir("DB"):
        if file.endswith(".profile.fas"):
            profile = os.path.join("DB", file)
            align_2_promals(profile, PuMA_results)
    
    #remove regions that are not well supported
    print "#####################################################################################"
    print "Using TrimAl to remove regions that are not supported"
    print "#####################################################################################"
    trimal_gappy(glob.glob("*.mafft.fas"))
    
    #concatenate the alignments
    print "#####################################################################################"
    print "Concatenating the E1, E2, and L1 alignments"
    print "#####################################################################################"
    concatenate(glob.glob("*.gappy.fas"))
    
    #use FastTree to build a tree
    print "#####################################################################################"
    print "Constructing a ML phylogenetic tree using FastTree"
    print "#####################################################################################"
    os.system( "FastTree E1-E2-L1.fas > E1-E2-L1.tre" )
    
    #find the closest neighbor of the unknown virus
    print "#####################################################################################"
    print "locating the closest neighbor of %s in the phylogeneytic tree" %(new_virus)
    print "#####################################################################################"
    find_closest_neighbor("E1-E2-L1.tre",new_virus, "./DB/L1.nt.fas")
    
    #align the new virus and its neighbor (maintain codons)
    reverse_translate("pairwise.fas")
    
    #calculate pairwise identity
    output = pairwise_identity("pairwise.fas.mafft.fas")
    
    
    print "\n\n\n\n\n"
    
    print "#####################################################################################"
    print "#####################################################################################"
    print "%s and %s are %i %% identitical across the L1 ORF" %(output[0],output[1],output[2])
    print "#####################################################################################"
    print "#####################################################################################"
 
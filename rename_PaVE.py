import re, sys
from Bio import SeqIO
from ete3 import Tree

#pip install ete3
#run as python rename_PaVE.py <file><type>
#kind = fasta, newick

def rename_PaVE(file, kind):
    if kind == "newick":
        with open(file,"rU") as f:
            for l in f:
                l=l.strip()
                t = Tree(l)
                for leaf in t.iter_leaves():
                    leaf.name = (leaf.name).split("|")[3].split("-")[0]
        t.write(format=2, outfile=file+".renamed.tre")
    else:
        with open(file+".renamed.fas","w") as out:
            for r in SeqIO.parse(file,"fasta"):
                print >>out, ">%s\n%s" %(r.id.split("|")[3].split(".")[0],r.seq)
        
    

file, kind = sys.argv[1], sys.argv[2]        
                


rename_PaVE(file, kind)

from Bio import SeqIO
from Bio.Seq import Seq
import os




def reverse_translate(nucleotide_fasta):
    with open(nucleotide_fasta+".aa.fas","w") as out:
        for record in SeqIO.parse(nucleotide_fasta,"fasta"):
            print >>out, ">"+record.id+"\n"+record.seq.translate()
    os.system( str( "/usr/local/bin/mafft --quiet "+nucleotide_fasta+".aa.fas > "+nucleotide_fasta+".aa.mafft.fas" )) #this can be replaced with any aligner
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
                    print nuc_rec.id, x, len(str(nuc_rec.seq))/3 #">"+nuc_rec.id+"\n"+"".join(sequence)
                    print >>out, ">"+nuc_rec.id+"\n"+"".join(sequence)
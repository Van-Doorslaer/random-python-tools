from Bio import SeqIO
import re

def genbank2sequin(genbank):
   for r in SeqIO.parse(genbank,"genbank"):
      with open(r.id+'.fsa', 'w') as fasta, open(r.id+".tbl", 'w') as table:
       print >> fasta, ">%s [organism= %s]\n%s" %(r.id, r.id, r.seq)
       print >>table, ">Feature "+r.id 
       for (index, feature) in enumerate(r.features):
        if feature.type ==  'CDS':
            m=re.search('\[(\d*)\:(\d*)\]',str(feature.location))
            if m:
               start = int(m.groups()[0])
               stop= int(m.groups()[1])
            print >>table, "%s\t%s\tgene\n\t\t\tgene\t%s" %(start+1, stop,feature.qualifiers['gene'][0])
            print >>table, "%s\t%s\tCDS\n\t\t\tproduct\t%s\n\t\t\tgene\t%s\n\t\t\tcodon_start\t%s" %(int(m.groups()[0])+1, m.groups()[1], feature.qualifiers['gene'][0], feature.qualifiers['gene'][0], "1")
        if feature.type ==  'mRNA':
               m=re.search('\[(\d*)\:(\d*)\].*\[(\d*)\:(\d*)\]',str(feature.location))
               if m:
                  start, SD, SA, stop = m.groups()[0],m.groups()[1],m.groups()[2],m.groups()[3]
                  print >>table, "%s\t%s\tmRNA\n%s\t%s\n\t\t\tproduct\t%s" %(int(start)+1, SD, int(SA)+1, stop,feature.qualifiers['gene'][0])
        if "misc" in feature.type:
            m=re.search('\[(\d*)\:(\d*)\]',str(feature.location))
            if m:
               start = int(m.groups()[0])
               stop= int(m.groups()[1])
            print >>table, "%s\t%s\tmisc_feature\n\t\t\tnote\t%s" %(start+1, stop, feature.qualifiers['note'][0])

genbank2sequin("example.gb")

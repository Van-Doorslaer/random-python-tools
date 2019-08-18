This script will help to convert a GenBank formated file for submission to GenBank

The python script turns the GenBank  files into the fasta and *.tbl files needed.
The template.sbt.txt file can be generated through the NCBI website and should be 
downloaded to the same folder.

Place the python script in the same folder as the GenBank files you want to convert
In the terminal type: "python genbank2sequin.py”
Make sure you also save the “template file” in the same folder
In terminal type:"tbl2asn -t template.sbt.txt -p . -M n -Z discrep -r ./sequin"


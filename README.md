# single_copy_gene_analysis_pangenome
a single copy gene analysis pangenome which allows for the orthology already computed and allows for the alignment and phylogeny building using three alignment approaches and also allows for the GTRCAT and GTRGAMA model with 1000 bootstrap. To check the file conversion it iterates the first files at the start and see if the number matches to the number of the genes present in each files. If the number doesnt matches then it will exit such as
```
#! /usr/bin/bash 
# Universitat Potsdam
# Author Gaurav Sablok
# date: 2024-2-16
# an end to end workflow for the complete analysis of the single copy orthologs from the pangenome
echo "starting analyzing the pangenome for the single copy orthologs"
read -r -p "enter the directory path to the single copy orthologs:" dirpath
read -r -p "enter the number of the species involved in the orthology search:" species
read -r -p "enter the number of the threads:" threads
directorypath="${dirpath}"
listedspecies="${species}"
cd "${directorypath}"
echo "checking whethere all the files are correcting"
echo "implementing a regaular expression for the file search and the number enumeration"
if [[ $(for i in *.fa; do grep ">" -c $i; done | head -n 1) == "${species}" ]]
 then 
        echo "files are all correct and the number of the species in all the files are present"
fi
```
Gaurav Sablok,\
Academic Staff Member,\
Bioinformatics,\
Institute for Biochemistry and Biology,\
University of Potsdam,\
Potsdam,Germany

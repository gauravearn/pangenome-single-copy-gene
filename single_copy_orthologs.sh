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
echo "starting the analysis"
declare -a alignmenttools=("macse", "prank", "muscle")
for i in alignmentools \
do \
  if [[ "${i}" == "macse" ]]
  then
        echo "alignining the pangenome using the pangenome using the macse" 
        for j in "${directorypath}"/*.fa
                java -jar -Xmx100g macse -prog alignSequences \
                                    -gc_def 12 -seq "$j" -out_AA "${i}%".AA -out_NT \
                                                        "${i}%".NT > "${i}".macse.run.log.txt
   fi \ 
        for i in *.NT; do mv "${i}" "${i%.*}".ntaligned.fasta; done
        echo "renamed the aligned files as ntaligned"
        for i in *.ntaligned.fasta
        do 
                trimal -in "${i}" -out "${i%.*}".trimmed.fasta -nogaps
        done
        wget https://github.com/marekborowiec/AMAS/AMAS.py 
        chmod 755 AMAS.py
        python3 AMAS.py -in *.fasta -f fasta -d dna
        mv concatenated.out macsealignmentconcatenated.fasta
        mv partitions.txt macsealignmentpartitions.txt
        raxmlHPC-PTHREADS -s macsealignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
                                                 -p 12345 -n macsephylogeny_GAMMA -T "${threads}" -N 50
        raxmlHPC-PTHREADS -s macsealignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
                                                            -p 12345 -n macsephylogeny_GTRCAT -T "${threads}" -N 50 -b 1000 
        echo "finishing up the analysis"   
   fi
   if [[ "${i}" == "prank" ]]
   then
        echo "aligning the pangenome using the prank probabilistic alignment"
        for j in "${directorypath}"/*.fa
        do       
                sudo apt-get install prank
                prank -d "${i}" -o "${i%.*}".prankaligned.fasta
        done
        for i in *.prankaligned.fasta
        do 
                trimal -in "${i}" -out "${i%.*}".pranktrimmed.fasta -nogaps
        done
        wget https://github.com/marekborowiec/AMAS/AMAS.py 
        chmod 755 AMAS.py
        python3 AMAS.py -in *.fasta -f fasta -d dna
        mv concatenated.out prankalignmentconcatenated.fasta
        mv partitions.txt prankalignmentpartitions.txt
        raxmlHPC-PTHREADS -s prankalignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
                                                               -p 12345 -n prankphylogeny_GAMMA -T "${threads}" -N 50
        raxmlHPC-PTHREADS -s prankalignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
                                                               -p 12345 -n prankphylogeny_GTRCAT -T "${threads}" -N 50 -b 1000
   fi
   if [[ "${i}" == "muscle" ]]
   then
        echo "aligning the pangenome using the muscle alignment"
        for j in "${directorypath}"/*.fa
        do   
                sudo apt-get install muscle
                muscle -in "${i}" -out "${i%.*}".musclealigned.fasta
        done
        for i in *.musclealigned.fasta
        do 
                trimal -in "${i}" -out "${i%.*}".muscletrimmed.fasta -nogaps
        done
        wget https://github.com/marekborowiec/AMAS/AMAS.py 
        chmod 755 AMAS.py
        python3 AMAS.py -in *.fasta -f fasta -d dna
        mv concatenated.out musclealignmentconcatenated.fasta
        mv partitions.txt musclealignmentpartitions.txt
        raxmlHPC-PTHREADS -s musclealignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
                                                                     -p 12345 -n musclephylogeny_GAMMA -T "${threads}" -N 50
        raxmlHPC-PTHREADS -s musclelignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
                                                                    -p 12345 -n musclephylogeny_GTRCAT -T "${threads}" -N 50 -b 1000
   fi
done

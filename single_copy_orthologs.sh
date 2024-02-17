#! /usr/bin/bash
# Universitat Potsdam
# Author Gaurav Sablok
# date: 2024-2-16
# updated 2024-2-17
# an end to end workflow for the complete analysis of the single copy orthologs from the pangenome
# it uses three alignment tools and then trims the alignments and concatenates the alignments and builds the phylogeny
echo "starting analyzing the pangenome for the single copy orthologs nucleotide sequences"
echo "alignment tools options are macse, prank, muscle or all"
read -r -p "enter the directory path to the single copy orthologs:" dirpath
read -r -p "enter the number of the species involved in the orthology search:" species
read -r -p "enter the number of the threads:" threads
read -r -p "enter the alignment tools:" tool
directorypath="${dirpath}"
listedspecies="${species}"
cd "${directorypath}"
echo "checking whethere all the files are correcting"
echo "implementing a regaular expression for the file search and the number enumeration"
if [[ $(for i in *.fa; do grep ">" -c $i; done | head -n 1) == "${species}" ]]; then
        echo "files are all correct and the number of the species in all the files are present"
fi
echo "starting the analysis"

if [[ "${tool}" == "macse" ]]; then
        echo "alignining the pangenome using the pangenome using the macse"
        for j in "${directorypath}"/*.fa; do
                java -jar -Xmx100g macse -prog alignSequences \
                        -gc_def 12 -seq "${j}" -out_AA "${j}%".AA -out_NT \
                        "${j}%".NT >"${j}".macse.run.log.txt
        done
        for i in *.NT; do
                mv "${i}" "${i%.*}".ntaligned.fasta
        done
        echo "renamed the aligned files as ntaligned"
        for i in *.ntaligned.fasta; do
                trimal -in "${i}" -out "${i%.*}".trimmed.fasta -nogaps
        done
        wget https://github.com/marekborowiec/AMAS/AMAS.py
        chmod 755 AMAS.py
        python3 AMAS.py -in *.fasta -f fasta -d dna
        mv concatenated.out macsealignmentconcatenated.fasta
        mv partitions.txt macsealignmentpartitions.txt
        iqtree --seqtype DNA -s macsealignmentconcatenated.fasta --alrt 1000 -b 1000 -T "${threads}"
        raxmlHPC-PTHREADS -s macsealignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
                -p 12345 -n macsephylogeny_GAMMA -T "${threads}" -N 50
        raxmlHPC-PTHREADS -s macsealignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
                -p 12345 -n macsephylogeny_GTRCAT -T "${threads}" -N 50 -b 1000
        echo "finishing up the analysis"
fi

if [[ "${tool}" == "prank" ]]; then
        echo "aligning the pangenome using the prank probabilistic alignment"
        for i in "${directorypath}"/*.fa; do
                sudo apt-get install prank
                prank -d "${i}" -o "${i%.*}".prankaligned.fasta
        done
        for i in *.prankaligned.fasta; do
                trimal -in "${i}" -out "${i%.*}".pranktrimmed.fasta -nogaps
        done
        wget https://github.com/marekborowiec/AMAS/AMAS.py
        chmod 755 AMAS.py
        python3 AMAS.py -in *.fasta -f fasta -d dna
        mv concatenated.out prankalignmentconcatenated.fasta
        mv partitions.txt prankalignmentpartitions.txt
        iqtree --seqtype DNA -s prankalignmentconcatenated.fasta --alrt 1000 -b 1000 -T "${threads}"
        raxmlHPC-PTHREADS -s prankalignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
                -p 12345 -n prankphylogeny_GAMMA -T "${threads}" -N 50
        raxmlHPC-PTHREADS -s prankalignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
                -p 12345 -n prankphylogeny_GTRCAT -T "${threads}" -N 50 -b 1000
fi

if [[ "${tool}" == "muscle" ]]; then
        echo "aligning the pangenome using the muscle alignment"
        for i in "${directorypath}"/*.fa; do
                sudo apt-get install muscle
                muscle -in "${i}" -out "${i%.*}".musclealigned.fasta
        done
        for i in *.musclealigned.fasta; do
                trimal -in "${i}" -out "${i%.*}".muscletrimmed.fasta -nogaps
        done
        wget https://github.com/marekborowiec/AMAS/AMAS.py
        chmod 755 AMAS.py
        python3 AMAS.py -in *.fasta -f fasta -d dna
        mv concatenated.out musclealignmentconcatenated.fasta
        mv partitions.txt musclealignmentpartitions.txt
        iqtree --seqtype DNA -s musclealignmentconcatenated.fasta --alrt 1000 -b 1000 -T "${threads}"
        raxmlHPC-PTHREADS -s musclealignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
                -p 12345 -n musclephylogeny_GAMMA -T "${threads}" -N 50
        raxmlHPC-PTHREADS -s musclelignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
                -p 12345 -n musclephylogeny_GTRCAT -T "${threads}" -N 50 -b 1000
fi

if [[ "${tool}" == "all" ]]; then

        echo "alignining the pangenome using the pangenome using the macse"
        for j in "${directorypath}"/*.fa; do
                java -jar -Xmx100g macse -prog alignSequences \
                        -gc_def 12 -seq "${j}" -out_AA "${j}%".AA -out_NT \
                        "${j}%".NT >"${j}".macse.run.log.txt
        done
        for i in *.NT; do
                mv "${i}" "${i%.*}".ntaligned.fasta
        done
        echo "renamed the aligned files as ntaligned"
        for i in *.ntaligned.fasta; do
                trimal -in "${i}" -out "${i%.*}".trimmed.fasta -nogaps
        done
        wget https://github.com/marekborowiec/AMAS/AMAS.py
        chmod 755 AMAS.py
        python3 AMAS.py -in *.fasta -f fasta -d dna
        mv concatenated.out macsealignmentconcatenated.fasta
        mv partitions.txt macsealignmentpartitions.txt
        iqtree --seqtype DNA -s macsealignmentconcatenated.fasta --alrt 1000 -b 1000 -T "${threads}"
        raxmlHPC-PTHREADS -s macsealignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
                -p 12345 -n macsephylogeny_GAMMA -T "${threads}" -N 50
        raxmlHPC-PTHREADS -s macsealignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
                -p 12345 -n macsephylogeny_GTRCAT -T "${threads}" -N 50 -b 1000
        echo "finishing up the analysis"
        echo "aligning the pangenome using the prank probabilistic alignment"
        for i in "${directorypath}"/*.fa; do
                sudo apt-get install prank
                prank -d "${i}" -o "${i%.*}".prankaligned.fasta
        done
        for i in *.prankaligned.fasta; do
                trimal -in "${i}" -out "${i%.*}".pranktrimmed.fasta -nogaps
        done
        wget https://github.com/marekborowiec/AMAS/AMAS.py
        chmod 755 AMAS.py
        python3 AMAS.py -in *.fasta -f fasta -d dna
        mv concatenated.out prankalignmentconcatenated.fasta
        mv partitions.txt prankalignmentpartitions.txt
        iqtree --seqtype DNA -s prankalignmentconcatenated.fasta --alrt 1000 -b 1000 -T "${threads}"
        raxmlHPC-PTHREADS -s prankalignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
                -p 12345 -n prankphylogeny_GAMMA -T "${threads}" -N 50
        raxmlHPC-PTHREADS -s prankalignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
                -p 12345 -n prankphylogeny_GTRCAT -T "${threads}" -N 50 -b 1000
        echo "aligning the pangenome using the muscle alignment"
        for i in "${directorypath}"/*.fa; do
                sudo apt-get install muscle
                muscle -in "${i}" -out "${i%.*}".musclealigned.fasta
        done
        for i in *.musclealigned.fasta; do
                trimal -in "${i}" -out "${i%.*}".muscletrimmed.fasta -nogaps
        done
        wget https://github.com/marekborowiec/AMAS/AMAS.py
        chmod 755 AMAS.py
        python3 AMAS.py -in *.fasta -f fasta -d dna
        mv concatenated.out musclealignmentconcatenated.fasta
        mv partitions.txt musclealignmentpartitions.txt
        iqtree --seqtype DNA -s musclealignmentconcatenated.fasta --alrt 1000 -b 1000 -T "${threads}"
        raxmlHPC-PTHREADS -s musclealignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
                -p 12345 -n musclephylogeny_GAMMA -T "${threads}" -N 50
        raxmlHPC-PTHREADS -s musclelignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
                -p 12345 -n musclephylogeny_GTRCAT -T "${threads}" -N 50 -b 1000
fi

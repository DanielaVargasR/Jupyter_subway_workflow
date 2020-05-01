#Assembly of pair-end sequences
ls *.fastq.gz | perl -pe 's/_R.*.fastq.gz//g' | sort | uniq >lista #generate list
ln -s scritps/assemblyPANDA.sh assemblyPANDA.sh #symbolic link
bash assemblyPANDA.sh NOMBRE_TRABAJO 
for N in `ls *.scr`; do qsub $N; done # se mandan al cluster los trabajos de ensamblado

#Rename and enumerate sequences per sample
perl scripts/header.fasta.numbers.pl PREFIX nombre_del_archivo.fasta 

#Concatenate all samples 
cat *.numbered.fasta >estudio_completo.fas

#Clustering
qsub -N NOMBRE_TRABAJO -b y -j y -cwd -V "cd-hit-est -c 0.97 -T 25 -M 0 -i estudio_completo.fas -o output.clstr"

#Edit
perl -pne 's/\t//g;s/^.*,//g;s/\.\.\..*$//g;s/\n/\t/g; s/\>Cluster\ /\n/g;s/\>//g; eof && do{chomp; print "$_ \n"; exit}' output.clstr.clstr >estudio.otu

#Representative sequences
pick_rep_set.py -i estudio.otu -f estudio_completo.fas -o rep_set.fna

##Remove no 16S sequences. Blast against gg_otus-13_8-release/rep_set/70_otus.fasta and remove those that did not aligned
parallel_assign_taxonomy_blast.py -i rep_set1.fna -o no16S_screen -r /qiime/gg_otus-13_8-release/rep_set/70_otus.fasta -t /qiime/gg_otus-13_8-release/taxonomy/70_otu_taxonomy.txt

cat no16S_screen/rep_set_tax_assignments.txt | grep -c "No blast hit"

#File screened
cat no16S_screen/rep_set_tax_assignments.txt | grep -v "No blast hit" | cut -f1 >ids_screened.txt

#Sequences removed
cat no16S_screen/rep_set_tax_assignments.txt | grep "No blast hit" | cut -f1 >ids_REMOVE_biom.txt

#Extracts sequences matching with 16S and generates a new file with representative sequences 
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' ids_screened.txt rep_set.fna >rep_set.screened.fna 

#Taxonomic assignment at 97%
parallel_assign_taxonomy_blast.py -i rep_set.screened.fna -o taxonomy -r /qiime/gg_otus-13_8-release/rep_set/97_otus.fasta -t /qiime/gg_otus-13_8-release/taxonomy/97_otu_taxonomy.txt

#OTU table
make_otu_table.py -i estudio.otu -t taxonomy/rep_set.screened_tax_assignments.txt -o estudio.biom 

#Remove OTUs with no hit or singletons
filter_otus_from_otu_table.py -i estudio.biom -e ids_REMOVE_biom.txt -o estudio_screened.biom -n2 ; mv estudio_screened.biom estudio.biom

#Chimera identification
parallel_align_seqs_pynast.py -i rep_set.screened.fna -o chimalign -X estudio

#Chimera removal
parallel_identify_chimeric_seqs.py -m blast_fragments -i rep_set.screened.fna -a chimalign/rep_set.screened_aligned.fasta -o estudio.chimera.txt -X chimerablast --id_to_taxonomy_fp /qiime/gg_otus-13_8-release/taxonomy/97_otu_taxonomy.txt -r /qiime/gg_otus-13_8-release/rep_set/97_otus.fasta

#Chimera removal from OTU table
filter_otus_from_otu_table.py -i estudio.biom -e estudio.chimera.txt -o estudio_chimera.biom; mv estudio_chimera.biom estudio.biom

#Generate tabular OTU table
biom convert --to-tsv -i estudio.biom -o estudio.biom.tsv --table-type "Taxon table" --header-key=taxonomy

#Phylogenetic tree
qsub -pe completenode 40 -N arbolote -b y -j y -cwd -V "export OMP_NUM_THREADS=40; FastTreeMP -nt -gtr chimalign/rep_set.screened_aligned.fasta >tree_daniela.nwk"

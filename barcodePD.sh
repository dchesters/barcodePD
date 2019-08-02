
###################################################################################################################################
#
#
# 
#	barcodePD: Phylogenetic Diversity from DNA barcodes 	  
#
#    	Copyright (C) 2019  Douglas Chesters
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#	contact address dc0357548934@live.co.uk
#
#
#
#
##########################################################################################################################################
#
#	
#	
#	
#	Pipeline and Perl scripts available at https://github.com/dchesters/barcodePD
#
#
#
#	CITATIONS
#
# 	For the pipeline as a whole: 
# 	  Wang MQ. et al. 
#	  Multiple components of plant diversity loss determine herbivore phylogenetic diversity in a subtropical forest experiment. 
#	  Journal of Ecology.
#
#	For making backbone constraint tree:
#	  Chesters D. 2017. 
#	  Construction of a Species-Level Tree-of-Life for the Insects and Utility in Taxonomic Profiling. 
#	  Systematic Biology 66: 426–439. 
#
#	For high throughput alignment:
#	  Chesters, D. 2019 
#	  aligner.pl. Computer software available at https://github.com/dchesters/aligner.
#
#
#
#	Change Log.
#
#	2019-AUG-02:
#		ALPHA version, under active development. please check back soon for updated versions.
# 		repository initiated at https://github.com/dchesters/barcodePD
#
#
#
#
#
#
#
#
##########################################################################################################################################


# Construct reference phylogeny, which is made from taxonomically identified DNA barcodes and constrained according to a published backbone phylogeny.
# Processed DNA barcodes in file named processed_refs.fas, and already aligned.
# Second input is backbone phylogeny. Note backbone phylogeny will probably need some processing, 
# For example we used Heikkilä" et al., which required conersion to Newick format, 
# Processing terminals such as removing some single quotes around terminal IDs and replacing whitespaces with underscores, 

# backbone constraining is in two parts, relational constraints and taxon/monophyly constraints (see Chesters 2017)

 # here instead should use newick method ... 
backbone=ditrysia_fig2.process_newick0.subtree_processed
sp_matrix=processed_refs.fas
# [script renamed, backbone_constraints_newick.pl to relational_constraints.pl]
perl relational_constraints.pl -node 7088 -seqfile $sp_matrix -treefile $backbone -outfile_prefix $backbone.constraint1 -backbone_terminal_format 2 -constrain_ranks suborder infraorder superfamily family subfamily genus subgenus
# output is relational_constraints.ft_format

########################

perl taxon_constraints.pl -node 7088 -seqfile $sp_matrix -treefile $intree -ignore_IDs_with_no_tax_assignment -print_taxon_constraints_only
mv taxon_constraints_only Coleoptera.taxon_constraints
# output will omit only order labeled data, although will include all insects.
# output is taxon_constraints_only, looks like: 
# Aaaba_fossicollis        0000000000001000000000000000000000000000000000000000000000000000000000000011000000000000000000000000000000000000000000000000001000000000000000000000000000000001000000000000000001000000000000000000000000000000000000000000000000000000010001000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000001000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
# Aaroniella_badonneli        0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001

# concatenate constraints, filter anything not of order under study, and entries labeled only to order level.
rm constraints_concat Coleoptera.constraints_concat
perl ~/usr_scripts/integrate_constraints.pl -constraint_files Coleoptera.taxon_constraints relational_constraints.ft_format -fasta processed_refs.fas -node 7088

# Use FastTree for initial topology, as it is one of few software with flexible constraining method:
FastTree_buffer80000 -log fasttree_logfile -gtr -gamma -nt -constraints constraints_concat < processed_refs.fas > processed_refs.ft_tax_and_rel_constraints.ftOUT

# Resolve any polytomies:
raxml-ng-0.6.0 --seed 123 --threads 2 --redo --msa processed_refs.fas --model DNA --tree-constraint processed_refs.ft_tax_and_rel_constraints.ftOUT --prefix processed_refs.ft_tax_and_rel_constraints.ftOUT
mv processed_refs.ft_tax_and_rel_constraints.ftOUT.raxml.bestTree refs.TR.nwk

# Root your reference tree here.

###################################################


# Inferring MOTU from plot level DNA barcodes.
# First is PTP, this is most involved, requires its own alignment and phylogeny inference step.
# Plot barcodes are processed and saved as file named update_data_caterpillars.fas
# Reference alignment (doublecheck it is aligned) saved in fasta format as processed_refs.fas

coi_file=update_data_caterpillars.fas

rm BlastOUT
blastn-2.2.28+-64bit -task blastn -query $coi_file -subject processed_refs.fas -out BlastOUT -word_size 10 -perc_identity 60 -max_target_seqs 1 -evalue 1e-8 -dust no -strand plus -outfmt '6 qseqid sseqid evalue pident length'
perl ~/usr_scripts/aligner.pl BlastOUT processed_refs.fas $coi_file unfiltered.alignerOUT

# Here add an outgroup to unfiltered.alignerOUT, perhaps labelled outgroup_Druceiella,
# save as unfiltered.alignerOUT2
# Note, here conducting an unconstrained tree search of All data (not just OTU)

# This step is too slow on literally all data, remove identical seqs
rm update_data_caterpillars.bcOUT
blastclust -i update_Leps.fas -o update_data_caterpillars.bcOUT -S 99.8 -p F -e F -a 4 -L .98
rm update_data_caterpillars.dereplicated
perl ~/usr_scripts/remove_fasta_entries.pl update_Leps.fas update_data_caterpillars.bcOUT update_data_caterpillars.dereplicated
rm update_data_caterpillars.dereplicated.phy
perl ~/usr_scripts/format_conversion.pl update_data_caterpillars.dereplicated update_data_caterpillars.dereplicated.phy fasta phylip
rm R*update_data_caterpillars.dereplicated
raxmlHPC-8.2.4 -s update_data_caterpillars.dereplicated.phy -n update_data_caterpillars.dereplicated -m GTRCAT -c 4 -p 123 -o outgroup_Druceiella

# Finally, now have rooted phylogeny of plot data, run PTP
rm update_data_caterpillars.dereplicated.ptp
python ~/software/ptp/SpeciesCounting/PTP.py -t RAxML_bestTree.update_data_caterpillars.dereplicated -p -pvalue 0.001 > update_data_caterpillars.dereplicated.ptp
perl ~/usr_scripts/parse_ptp_output.pl update_data_caterpillars.dereplicated.ptp ptp_results


#############################


# Probably the simplest method for making MOTU, and is robust software.
# ALthough criteria it uses is rudimentary.
# Default parameters: 97.8 percent identity over 98 percent of seq. 
blastclust -i $coi_file -o $coi_file.bcOUT -S 97.8 -p F -e F -a 4 -L .98
# Running time for our data ~10 minutes

# scan output, how many MOTU:
wc -l $coi_file.bcOUT

# Get sequence of single member of each OTU
# set this variable in the script:$rename_OTUs = 1;
rm OTU.fas
coi_file=moth_COI.fas
otu_results=moth_COI.fas.bcOUT
perl ~/usr_scripts/filter_otus.pl $coi_file $otu_results OTU.fas
grep ">" OTU.fas | wc -l

#########################################


# Now align MOTU representatives
reference_alignment=processed_refs.fas
rm BlastOUT alignerOUT
blastn-2.2.28+-64bit -task blastn -query OTU.fas -subject processed_refs.fas -out BlastOUT -word_size 10 -perc_identity 60 -max_target_seqs 1 -evalue 1e-8 -dust no -strand plus -outfmt '6 qseqid sseqid evalue pident length'
perl ~/usr_scripts/aligner.pl BlastOUT processed_refs.fas OTU.fas alignerOUT

# combine references and queries
cat alignerOUT processed_refs.fas > OTU_and_REFs.fas


# here add outgroup, perhaps a non-Ditrysia Neolepidoptera, Druceiella sp. BOLD:AAI6074
>outgroup_Druceiella
ATTGGCGATGATCAAATCTATAATGTTATCGTAACAGCTCATGCTTTCATTATAATTTTTTTTATAGTTATACCAATTATAATTGGAGGATTTGGAAATTGATTAGTTCCTTTGATATTAGGAGCCCCTGATATAGCTTTCCCACGAATAAATAATATGAGATTTTGATTATTACCACCATCATTAATACTACTAATTTCAAGAAGAATTGTAGAAAATGGAGCAGGTACAGGATGAACTGTTTATCCCCCACTATCATCAAACATTGCACATACAGGAAGATCTGTAGATTTAGCTATTTTTTCACTACATTTAGCTGGTATCTCATCTATTTTAGGTGCAGTAAATTTTATTACTACTGTAATTAATATACGAACAGAAGGAATATCTTTTGATCGCATACCACTATTTGTATGAAGAGTTGCTATTACTGCTTTATTATTATTATTATCCTTACCGGTATTAGCCGGAGCTATTACTATATTATTAACAGATCGAAATTTAAATACCTCATTTTTTGACCCTGCTGGGGGTGGTGATCCAATTTTATATCAACATTTATTT---

# Here consider removing problematic OTU,
# Those with very long branches, which can indicate sequencing error, alignment error,
# We removed OTU401, OTU384, OTU462
# saved as OTU_and_REFs2.fas

#  This is the main tree search of OTU and references
raxmlHPC-8.2.4 -s OTU_and_REFs.phy -n OTU_and_REFs -m GTRCAT -c 4 -p 123 -g refs.TR.nwk -o outgroup_Druceiella


# Here, to strengthen results, try incorporating uncertainty by calculating PDs on bootstrapped trees
# Note, regular raxml bootstrap is of no use, as need branch lengths for calculating PD

# create bootstrapped alignments:
raxmlHPC-8.2.4 -f j -b 123 -# 100 -s OTU_and_REFs.phy -n OTU_and_REFs.boot_alignments -m GTRCAT -c 4 -p 123 -g refs.TR.nwk -o outgroup_Druceiella

# parallelized:
no_cpus=8
count=0
for i in {0..99};do echo "$i run RAxML";
  ##################
  raxmlHPC-8.2.4 -D -e 1.0 -s OTU_and_REFs.phy.BS$i -n OTU_and_REFs.phy.BS$i -m GTRCAT -c 4 -p 123 -g refs.TR.nwk -o outgroup_Druceiella &
  ##################
  let count+=1
  [[ $((count%no_cpus)) -eq 0 ]] && wait
done

# note, Ape has a bug reading branch lengths if multiple trees in one file. better read them seperatly

# following files not needed. bestTrees retained.
rm RAxML_log.OTU_and_REFs.phy.BS* RAxML_result.OTU_and_REFs.phy.BS* RAxML_info.OTU_and_REFs.phy.BS*






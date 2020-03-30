

# 
# barcodePD_v2.sh 
#	Douglas Chesters, Institute of Zoology, Chinese Academy of Sciences.
# 	Available at github.com/dchesters/barcodePD
# 
# citation:
#	Wang, Li, Chesters, Bruelheide, Ma, Zhou, Staab, Zhu, Schuldt.
# 	Host functional and phylogenetic composition rather than host diversity structure plant-herbivore networks. 
#	in revision.
# 
# requirements:
# backbone trees, processed; we used phylogenies from Kawahara et al. 2019, Heikkila et al. 2015, Regier et al. 2016.
# tree|stake scripts taxon_table.pl, relational_constraints.pl, taxon_constraints.pl, see github.com/dchesters/treestake for latest versions
# ITIS and NCBI taxonomic databases (see treestake documentation)
# Perl scripts parse_hits.pl, unalign.pl, species_filter.pl, format_conversion.pl, concatenate_v2.pl, process_BOLD_API_seqIDs.pl, 
#  aligner.pl (github.com/dchesters/aligner)
# installed software raxmlHPC, Emboss, blastn.
# DNA barcodes of plot data, here named OTUs_MQ20200322.fas 
# and optional, plot reference data, which is DNA barcodes with taxonomic labels, here named processed_refs_MQ20200325.fas
# These last 2 two files can be obtained from the authors for replicating this analysis.



# First, download all DNA barcodes for families expected in the plot data. 
# WARNING! check each of these has actually downloaded, this may not be apparent by file size. 
# Sometimes html of several kb will be returned instead of sequence, if for example bold server is being maintained.
# From some locations, connection can be unreliable, if this the case beter do these one at a time and check screen output for each.
# Any file suspected not or only partially downloaded, just try again.

wget --output-document=Crambidae.fas http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Crambidae
wget --output-document=Drepanidae.fas http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Drepanidae
wget --output-document=Erebidae.fas http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Erebidae
wget --output-document=Euteliidae.fas http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Euteliidae
wget --output-document=Gelechiidae.fas http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Gelechiidae
wget --output-document=Geometridae.fas http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Geometridae
wget --output-document=Gracillariidae.fas http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Gracillariidae
wget --output-document=Lasiocampidae.fas http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Lasiocampidae
wget --output-document=Limacodidae.fas http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Limacodidae
wget --output-document=Lycaenidae.fas http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Lycaenidae
wget --output-document=Noctuidae.fas http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Noctuidae
wget --output-document=Nolidae.fas http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Nolidae
wget --output-document=Notodontidae.fas http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Notodontidae
wget --output-document=Nymphalidae.fas http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Nymphalidae
wget --output-document=Papilionidae.fas http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Papilionidae
wget --output-document=Psychidae.fas http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Psychidae
wget --output-document=Pyralidae.fas http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Pyralidae
wget --output-document=Saturniidae.fas http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Saturniidae
wget --output-document=Sphingidae.fas http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Sphingidae
wget --output-document=Tortricidae.fas http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Tortricidae
wget --output-document=Zygaenidae.fas http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Zygaenidae

# Combine all to a single file
cat Crambidae.fas Drepanidae.fas Erebidae.fas Euteliidae.fas Gelechiidae.fas Geometridae.fas Gracillariidae.fas Lasiocampidae.fas Limacodidae.fas Lycaenidae.fas Noctuidae.fas Nolidae.fas Notodontidae.fas Nymphalidae.fas Papilionidae.fas Psychidae.fas Pyralidae.fas Saturniidae.fas Sphingidae.fas Tortricidae.fas Zygaenidae.fas > BOLD_Lep.fas
# Following is bash command used regularly for counting entries in fasta file,
# it is a good idea to keep notes of these, to track degree of changes upon adjusting various parameters
grep ">" BOLD_Lep.fas | wc -l
# count entries, as of 2020-MARCH, 716660


# Get 5 prime COI barcodes only (sometimes there are other loci).
# Note there are format inconsistencies, sometimes has accession 
#  (for which there would be an additional delimiter and this would be used: |COI-5P|)
# But roughly half have no accession after marker_code.
grep -A 1 "|COI-5P" BOLD_Lep.fas > BOLD_Lep.COI

grep ">" BOLD_Lep.COI | wc -l
# count 5' COIs only: 700,673

# Extract only entries which are labelled fully to species level. 
# Partially labelled data can be difficult/impossible to parse (eg >BLPDV11202-19|zygJanzen01 Janzen01|COI-5P), 
# difficult to assign taxonomy (due to rank variation), and often cannot be constrained in the phylogenetic analysis.
# Strict binomial parse:
grep -A 1 "|[A-Z][a-z]\+ [a-z]\+|" BOLD_Lep.COI > BOLD_Lep.COI.fas

grep ">" BOLD_Lep.COI.fas | wc -l
# count: 488083

# Processes the BOLD id format into preffered format which is taxonomic name folowed by accession
# nb Good practice to delete previous output files as a step immediatly prior to anything that produces them.
# Also, from here, change filepaths according to your setup
rm BOLD_Lep.COI.parsed
perl ~/usr_scripts/process_BOLD_API_seqIDs.pl BOLD_Lep.COI.fas BOLD_Lep.COI.parsed

# Dereplicate, leaving a single exemplar for each species.
# [Script arguments probably not required: -trim_accession  -identified_species_only]
rm BOLD_Lep.COI.parsed.ID_filtered
perl ~/usr_scripts/species_filter.pl -in BOLD_Lep.COI.parsed -format 1 -filter_method 2
# count:36071

# Dataset will be a mix of sequences aligned and unaligned, thus they will be re-aligned. 
# In this case, always safer to remove any previously inserted gaps, due to default behaviour of some alignment software.
rm BOLD_Lep.COI.parsed.ID_filtered.unaligned
perl ~/usr_scripts/alignment_processing/unalign.pl BOLD_Lep.COI.parsed.ID_filtered

# processing BOLD barcodes complete.

###############################################################################################################





# Using the plot data as queries, search the BOLD database for anything broadly similar.
# The option -outfmt 6 allows to specify exactly the fields which are output.
# Columns output are aligned with an old blast output parser used next, but you can change if wanting further details
# Either here or in the following step, you might want to experiment with parameters, particularly -perc_identity.
# I tried several, -perc_identity of 95 left too many of the plot entities with nothing of broad similarity in the reference data,
# thus often placing in random locations. perc_identity of 90 gave quite a huge and unwieldy reference dataset, the tree search became very slow.
# Thus, default = 92.
# Other parameters maybe to look at -evalue, but best to keep stringent, since often even a small number of spurious alignments can cause large headaches later.
# and max_target_seqs, similarly to perc_identity, can be increased for a more comprehensive reference phylogeny, 
# but the extra members wont help much in query placement.
rm blast_against_BOLD
blastn-2.2.28+-64bit -task blastn -subject BOLD_Lep.COI.parsed.ID_filtered.unaligned -out blast_against_BOLD -dust no -perc_identity 92 -evalue 1e-10 -query OTUs_MQ20200322.fas -num_threads 1 -max_target_seqs 40 -outfmt '6 qseqid sseqid evalue pident length sstart send qframe sframe'

# Parse the blast output.
rm blast_against_BOLD.retreived
perl ~/usr_scripts/parse_hits.pl BOLD_Lep.COI.parsed.ID_filtered.unaligned blast_against_BOLD 90 1 450 2 blastdbcmd-2.2.28+-64bit
grep ">" blast_against_BOLD.retreived | wc -l
# Key parameter here is length of alignment between query and reference, less than which is not incorporated.
# Do not worry too much about allowing alignments quite a bit less than average full length. 
# A reference aligned at something quite low such as 450 bp, is still perfectly capable of identifying to very usefull taxonomic level. 
# and will be valuable if barcodes are sparse for that taxa.
# There was a specific case of this in the current dataset, ie if alignment threshold is increased to 500,
# one key reference species will be ommited, with the result that one OTU cannot be placed and jumps into some random location.




# Sequence alignment. Would be reasonable to assume to be trivial to align a length invarient marker. 
# this is not the case for a number of reasons (wont get into here).
# Alignment is conducted against a 'profile alignment' which must be prealigned, visually verified, 
# and current setup requires all sequences in the profile to be without gaps, including terminal gaps.
# NW parameters used here have been tested, were optimized on a Apoidae dataset previously
# aligner.pl requires Emboss to be installed on your system.
reference_alignment=processed_refs_MQ.aligner_profile
unaligned=blast_against_BOLD.retreived
rm BlastOUT newBOLD_aligned
blastn-2.2.28+-64bit -task blastn -query $unaligned -subject $reference_alignment -out BlastOUT -word_size 10 -perc_identity 60 -max_target_seqs 1 -evalue 1e-8 -dust no -strand plus -outfmt '6 qseqid sseqid evalue pident length'
perl ~/usr_scripts/aligner.pl -blast_results BlastOUT -references $reference_alignment -subject $unaligned -outfile newBOLD_aligned -gapopen 20.0 -gapextend 0.5 -endweight Y

# finally , append new bold data 
# Sequence files here are the plot OTU, plot references, BOLD references
# Making 2 new files, one with and one without the OTU, as these are integrated at different steps.
cat OTUs_MQ20200322.fas processed_refs_MQ20200325.fas newBOLD_aligned > OTU_BEFCNmoths_and_BOLD_references
cat processed_refs_MQ20200325.fas newBOLD_aligned > BEFCNmoths_and_BOLD_references

grep ">" OTU_BEFCNmoths_and_BOLD_references | wc -l
# 4119




# this section is based on github.com/dchesters/treestake. 
# modified as per Wang et al., with 3 backbone trees incorporated.
# Backbone trees probably require some processing (not shown), but as a quick guide for the current study,
# the 1KITE tree, available in Newick format in their supplemental data, had higher taxonomic ranks prepended to terminal 
# labels and needed to be removed. the Heikkilä tree was recieved from the authors in Nexus format,
# thus convered to Newick in Figtree, and some modification to IDs 
# (eg quoted labeled with taxonomic acronyms 'Jana_nr._eurymas' simpified to Jana_sp)


# make a taxon table for Leps using both ITIS and NCBI databases (see treestake documentation)
ITIS_database=/home/douglas/scripted_analyses/sDiv/data/taxa_3842024.txt
NCBI_names=~/scripted_analyses/sDiv/data/names.dmp
NCBI_nodes=~/scripted_analyses/sDiv/data/nodes.dmp
taxon_of_interest=Lepidoptera
rm Taxon_Table
perl ~/usr_scripts/taxon_table.pl -node 7088 -ITIS $ITIS_database $taxon_of_interest -NCBI $NCBI_names $NCBI_nodes -outfile Taxon_Table


# infer relational constraints, these are applied as a newick constraint tree input to raxml
# We chose to use the 1KITE tree, for reasons of taxonomic overlap, you might opt for a different one.
# Format for backbone tree is plain Newick, terminal labels are Linnean binomials in format Genus_species. 
# A given species should have no more than one terminal. Limited species ambiguity is allowed, e.g. Genus_sp1 
# Branch lengths and node support ignored if present.
backbone=1KITE2019
sp_matrix=BEFCNmoths_and_BOLD_references

taxon_table=Taxon_Table
rm constraint_phylogeny_pruned2.nwk constr2.1KITE
perl ~/usr_scripts/relational_constraints.pl -seqfile $sp_matrix -treefile $backbone -outfile_prefix $backbone.constraint1 -backbone_terminal_format 0 -taxon_table $taxon_table
mv constraint_phylogeny_pruned2.nwk constr2.1KITE
# file constr2.1KITE will be used as contraint tree input into Raxml.



#####################################################################################################################################################



# Infer taxon constraints for both backbone trees, these are applied as mock characters
# Note, when applied this way, these are not *strictly* constrained.
# This will be noticable in the presence of conflicting constraints,
# for instance, Regier et al find Noctudidae monophyletic, and this is incorporated as characters into the matrix here,
# However, in the 1KITE topology Noctudidae are paraphyletic with regard to Lymantridae,
# and thus as the 1KITE tree is selected for use as relational constraints, 
# Noctudidae will not be monophyletic in the resulting output (even if you try increasing weights in taxon_constraints)

# As before, always have a step for deleting previous files.
# since it is not always obvious a process has crashed and you inadvertently use some old files.
rm taxon_constraints.fas taxon_constraints.1KITE.fas taxon_constraints.Heikkila.fas taxon_constraints.Regier2016.fas

sp_matrix=BEFCNmoths_and_BOLD_references
taxon_table=Taxon_Table

# If subsequent phylogenetic analysis step are running fast, you can use higher character weighting to ensure (non-conflicting) taxon are constrained.
# But these can make unwieldy matrices
perl ~/usr_scripts/taxon_constraints.pl -character_weight 6 -seqfile $sp_matrix -treefile ditrysia_fig2.process_newick0.subtree_processed -print_taxon_constraints_only -backbone_terminal_format 0 -taxon_table $taxon_table
mv taxon_constraints.fas taxon_constraints.Heikkila.fas
mv taxon_constraints_INFO taxon_constraints_INFO.Heikkila; mv taxon_constraints_genus_level_table taxon_constraints_genus_level_table.Heikkila; mv taxon_constraints_output_LOG taxon_constraints_output_LOG.Heikkila

perl ~/usr_scripts/taxon_constraints.pl -character_weight 6 -seqfile $sp_matrix -treefile 1KITE2019 -print_taxon_constraints_only -backbone_terminal_format 0 -taxon_table $taxon_table
mv taxon_constraints.fas taxon_constraints.1KITE.fas
mv taxon_constraints_INFO taxon_constraints_INFO.1KITE; mv taxon_constraints_genus_level_table taxon_constraints_genus_level_table.1KITE; mv taxon_constraints_output_LOG taxon_constraints_output_LOG.1KITE

perl ~/usr_scripts/taxon_constraints.pl -character_weight 6 -seqfile $sp_matrix -treefile Regier2016 -print_taxon_constraints_only -backbone_terminal_format 0 -taxon_table $taxon_table
mv taxon_constraints.fas taxon_constraints.Regier2016.fas
mv taxon_constraints_INFO taxon_constraints_INFO.Regier2016; mv taxon_constraints_genus_level_table taxon_constraints_genus_level_table.Regier2016; mv taxon_constraints_output_LOG taxon_constraints_output_LOG.Regier2016
# note There are several useful output files produced in addition to the constraint matrix.



################################################


# concatenate the dna sequences to the character constraints
# if you want to incorporate another backbone tree, then add another ? to required data (and the other way round).
rm current_supermatrix2 dna_and_char_constraints12
perl ~/usr_scripts/concatenate_v2.pl -missing_data_char ? -remove_accession 0 -required_data ???? -matrices $sp_matrix taxon_constraints.1KITE.fas taxon_constraints.Heikkila.fas taxon_constraints.Regier2016.fas
mv current_supermatrix2 dna_and_char_constraints12

# here make a partition file, worth doing manually just to check and think about the matrix that will be used. 
# Save as partitionfile.1
# Partition file should probably look like:
DNA, barcode1 = 1-567\3
DNA, barcode2 = 2-567\3
DNA, barcode3 = 3-567\3
BIN, constraint_binary = 568-1671


# Here is a key step, the constrained tree search of reference data.
# Dont skip branch length optimization; although best tree (output with branchlengths) not used in analysis, 
#  it is useful for visualizing where constraints applied
rm dna_and_char_constraints12.phy
perl ~/usr_scripts/format_conversion.pl dna_and_char_constraints12 dna_and_char_constraints12.phy fasta phylip
constraint_tree=constr2.1KITE
rm R*constrained12
raxmlHPC-8.2.4 -q partitionfile.1 -g $constraint_tree -p 123 -s dna_and_char_constraints12.phy -n constrained12 -m GTRCAT -c 24 -o outgroup_Druceiella


# This is not an essential step, but if you want to plot the reference tree only, with 'actual' branchlengths.
#  ... optimize branch lengths on sequence data only
perl ~/usr_scripts/format_conversion.pl $sp_matrix $sp_matrix.phy fasta phylip
rm R*constrained1b
raxmlHPC-8.2.4 -g RAxML_result.constrained1 -p 12345 -s $sp_matrix.phy -n constrained1b -m GTRCAT -c 24 -o outgroup_Druceiella


# Final step is place plot level OTU onto the tree .... 
# Optional to regular tree search as done below, is the phylogenetic placement function;
#	which you might consider if your query data is a lot more substantial (such as NGS generated)
# Also consider doing a bootstrap to reveal otu unreliably placed.
rm OTU_BEFCNmoths_and_BOLD_references.phy
perl ~/usr_scripts/format_conversion.pl OTU_BEFCNmoths_and_BOLD_references OTU_BEFCNmoths_and_BOLD_references.phy fasta phylip
rm R*OTU_and_REFs.costrained_phylo12
raxmlHPC-8.2.4 -g RAxML_result.constrained12 -p 123 -s OTU_BEFCNmoths_and_BOLD_references.phy -n OTU_and_REFs.costrained_phylo12 -m GTRCAT -c 24 -o outgroup_Druceiella




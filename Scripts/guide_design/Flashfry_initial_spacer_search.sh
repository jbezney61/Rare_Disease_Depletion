#!/usr/bin/bash
#SBATCH --job-name=Flashfry_hits
#SBATCH --output=Flashfry_hits.%j.out
#SBATCH --error=Flashfry_hits.%j.err

eval "$(conda shell.bash hook)"
conda activate FLASHFRY

########################
# Build indexed database 
########################

# targeted contigs 
java -Xmx200g -jar ,/FLASHFRY/FlashFry-assembly-1.15.jar index \
	--tmpLocation tmp \
	--database Targeted_contigs_cas9ngg_database \
	--reference FASTA_targeted_contigs_GREGOR_10k_depletion.fa \
	--enzyme spcas9ngg

#non-targeted, preserved contigs 
java -Xmx200g -jar ./FLASHFRY/FlashFry-assembly-1.15.jar index \
	--tmpLocation tmp \
	--database NON_targeted_preserved_contigs_cas9ngg_database \
	--reference FASTA_NON_targeted_contigs_GREGOR_preserved.fa \
	--enzyme spcas9ngg

################
# On-target hits
################

#determine hits, had to splice the file to account for memory 
for i in {a..m}; do
    for j in {a..z}; do
        file="x$i$j"
        if [[ -f $file ]]; then
            java -Xmx300g -jar ./FLASHFRY/FlashFry-assembly-1.15.jar discover \
                --database Targeted_contigs_cas9ngg_database \
                --fasta $file \
                --output ${file}_Targeted_contigs_guide_design.output \
                --forceLinear --maxMismatch 3 --maximumOffTargets 2000000 \
                --flankingSequence 10 --minGC 0.35 --maxGC 0.80
        fi
    done
done

################
# Off-target hits
################

for i in {a..m}; do
    for j in {a..z}; do
        file="x$i$j"
        if [[ -f $file ]]; then
            java -Xmx300g -jar ./FLASHFRY/FlashFry-assembly-1.15.jar discover \
                --database NON_targeted_preserved_contigs_cas9ngg_database \
                --fasta $file \
                --output ${file}_NON_targeted_contigs_guide_design.output \
                --forceLinear --maxMismatch 3 \
                --flankingSequence 10 --minGC 0.35 --maxGC 0.80
        fi
    done
done

# guides were then removed if they had poor IVT production
# and if they were identified in targeting the preserved contigs 
# the remaining were then scored

#score the guide hits
java -Xmx300g -jar ./FLASHFRY/FlashFry-assembly-1.15.jar score \
	-input 3million_filtered_guides_for_scoring.out \
	--output GUIDES_3m_selection_design_0_offtargets.output.scored \
	--scoringMetrics doench2014ontarget,moreno2015,dangerous \
	--maxMismatch 3 --database NON_targeted_preserved_contigs_cas9ngg_database \
	--inputAnnotationBed Targeted_contigs_bed_intervals.bed

 
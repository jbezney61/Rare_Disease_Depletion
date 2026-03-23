# Rare_Disease_Depletion
Improving the diagnostic capability of RNA-seq for rare disease 


Differential change in splice junctions from Leafcutter
To visualize for every gene the change in intron clusters

1) leverage the prepared leafviz_results.RData in the leafviz repository
2) install leafviz locally 
3) Open the shiny app locally 

# running shiny appy locally 
    ## in R:
    install.packages("remotes")
    remotes::install_github("jackhump/leafviz")

    ## on the command line:
    git clone https://github.com/jackhump/leafviz.git

    # now run shiny app from command line 
    R -q -e "library(leafviz);options(shiny.host='127.0.0.1', shiny.port=8787, shiny.launch.browser=TRUE);leafviz::leafviz('leafviz_results.RData')"

    # now go to browse splice differences
    http://localhost:8787

LeafViz results were prepared following the LeafViz tutorial https://github.com/jackhump/leafviz

1) Gencode V48 annotation was prepared using leafviz gtf2leafcutter.pl
2) Results were prepared using prepare_results.R with:
	- metadata indicating the 96 standard versus the 96 depleted names 
	- _perind_numers.counts.gz 
	- _cluster_significance.txt 
	- _effect_sizes.txt 
	- gencode_v48 annotations


Fraser analysis 

1) The Fraser analysis was run using the pipelines documented here:
https://github.com/maurermaggie/Transcriptome_Wide_Splicing_Analysis/tree/main


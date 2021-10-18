## Genome-wide analysis of *Schistosoma mansoni* population structure and praziquantel drug selection pressure within Ugandan hot-spot communities.

This is code to accompany the above paper, and regenerate the included figures.


### egg reduction rate and GWAS of ERR vs. genome-wide variation
#

The script [ERR_estimates.R](ERR_estimates.R) takes the [egg count data](Eggs_88 - Eggs_88.tsv) and generates Egg reduction rates per individual and marginal estimates per village and per treatment arm, and draws the figures. It uses the generalised linear mixed-effect model developed by Martin Walker (Martin Walker (mwalker@rvc.ac.uk), and is very closely based on code by Tom Crellen (thomas.crellen@bdi.ox.ac.uk), and Martin, plus a tiny bit by James Cotton (jc17@sanger.ac.uk. For an explanation of the mode, see:

Martin Walker, Thomas S. Churcher, María-Gloria Basáñez,
Models for measuring anthelmintic drug efficacy for parasitologists,
*Trends in Parasitology* 30:528-537 (2014).

 There is a bit more information in [a separate github repository](https://github.com/jacotton/ERR_mixed_model)

The ERR estimates per individual (written as file JV.ERR phenotypes) are then used as input to run plink v1.90b6, with command line:

> plink --vcf Bi_174_new.recode.vcf --pheno JV.ERR.phenotypes --assoc --allow-extra-chr --allow-no-sex --qt-means

Script [ERR_assoc_plots.R](ERR_assoc_plots.R) generates the final Manhattan plot and QQplot figures.
This also reads Schistosoma_mansoni_v7.fa.fai, which is generated with 

> samtools faidx Schistosoma_mansoni_v7.fa

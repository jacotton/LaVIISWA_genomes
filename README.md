## Genome-wide analysis of *Schistosoma mansoni* population structure and praziquantel drug selection pressure within Ugandan hot-spot communities.

This is code to accompany the above paper, and regenerate the included figures.


### egg reduction rate and GWAS of ERR vs. genome-wide variation
#

The script [ERR_estimates.R](ERR_estimates.R) takes the [egg count data](Eggs_88 - Eggs_88.tsv) and generates Egg reduction rates per individual and marginal estimates per village and per treatment arm, and draws the figures. It uses the generalised linear mixed-effect model developed by Martin Walker (Martin Walker (mwalker@rvc.ac.uk), and is very closely based on code by Tom Crellen (thomas.crellen@bdi.ox.ac.uk), and Martin, plus a tiny bit by James Cotton (jc17@sanger.ac.uk). For an explanation of the mode, see:

Martin Walker, Thomas S. Churcher, María-Gloria Basáñez,
Models for measuring anthelmintic drug efficacy for parasitologists,
*Trends in Parasitology* 30:528-537 (2014).

 There is a bit more information in [a separate github repository](https://github.com/jacotton/ERR_mixed_model)

The ERR estimates per individual (written as file JV.ERR phenotypes) are then used as input to run plink v1.90b6, with command line:

> plink --vcf Bi_174_new.recode.vcf --pheno JV.ERR.phenotypes --assoc --allow-extra-chr --allow-no-sex --qt-means

this genrates a file called plink.qassoc. Script [ERR_assoc_plots.R](ERR_assoc_plots.R) generates the final Manhattan plot and QQplot figures.
This also reads [Schistosoma_mansoni_v7.fa.fai](Schistosoma_mansoni_v7.fa.fai), which is generated with 

> samtools faidx Schistosoma_mansoni_v7.fa

### analysis of between-site Fst and gravity model

original data is [fst_villages_3.csv](fst_villages_3.csv)

Code to fit these models and some brief explanation, together with code to reproduce figure 4 and S2 figure in the manuscript are in [GravityModel.R](GravityModel.R)

#### mapping

code to redraw map in figure 1 is self-contained in Mapping.R, with links to layer data included in that file.

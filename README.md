# PPI-ID: Protein-Protein Interaction Predictor
## About
This is a structural bioinformatics tool that aids domain-domain or domain-motif interaction prediction between proteins. This tool can help researchers refine inputted amino acid sequences when predicting protein-protein interactions using AlphaFold Multimer, and it can also help researchers identify interacting regions as found between 2 proteins in a PDB file. To use, provide protein domain or small linear motif (SLiM) information generated by either InterPro or the ELM database. Interaction predictions, as sorted in increasing order of combined probability of each domain/SLiM sequence occurring by chance, can be downloaded as a tsv by the user. 

After running predictions on AlphaFold, the user is able to upload the resulting PDB file and interact with the data frame of predicted interactions. Please note that in order for data frame and the 3D molecular model to properly interact, the sequences submitted to InterPro/ELM must be the exact same as the sequences that were folded. Also note that the order of protein sequences folded in AlphaFold must match that of Protein 1 and Protein 2 tsv's uploaded to the tool. 

This tool takes advantage of a compiled dataset of domain-domain interactions from the 3did (2022 release) and DOMINE databases. As a result, domains are identified by their Pfam ID, and domain-SLiM interactions are provided by the Eukaryotic Linear Motif (ELM) Database. Potential interactions are determined according to the appropriate algorithm implemented in the R script.

The web version of this tool can be accessed at the following link, http://ppi-id.biosci.utexas.edu:7215/.

## To Use
One is able to run this Shiny app from RStudios. To use, simply download the ppid.R script as well as the compiled_interactions.csv and interaction_id.tsv files. Copy and paste the file path name at the appropriate spot in the script, which will be annotated in-script for ease of identification. Afterwards, one is able to run the entire script to use the tool. I am currently working on getting the tool published on a server so that anyone can access the app directly from a link, without having to execute any code on RStudios.

## References
Raghavachari B, Tasneem A, Przytycka T, and Jothi R. (2008). "DOMINE: A database of protein domain interactions." _Nucl. Acids Res., 36 (Database Issue)_, D656-661, doi: 10.1093/nar/gkm761.

Rego N, Koes D. (2015). "3Dmol.js: molecular visualization with WebGL." _Bioinformatics, 31_(8):1322–1324, doi: 10.1093/bioinformatics/btu829.

Roberto Mosca, Arnaud Ceol, Amelie Stein, Roger Olivella & Patrick Aloy. (2014) "3did: a catalogue of domain-based interactions of known three-dimensional structure." _Nucl. Acids Res., 42_(D1):D374-D379, doi: 10.1093/nar/gkt887.

Yellaboina S, Tasneem A, Zaykin DV, Raghavachari B, and Jothi R. (2011). "DOMINE: A comprehensive collection of known and predicted domain-domain interactions." _Nucl. Acids Res., 39 (Database Issue)_, D730-735, doi: 10.1093/nar/gkq1229.

## Collaborators:
Haley Vy Goodwin

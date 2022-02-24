# Location of Coral aDNA Pipeline Files

Comp Methods folder is the location of all processing scripts used.
 - **aDNA_remote_workflow.sh**: software run remotely on supercomputing cluster to
 create required files.
 - **chr_namer.R, cmdJack_Scaffold.R**: Rscripts to be run remotely used in the
  workflow. These create jackknifed files at the scaffold level for the
  projection-based pipeline.
- **putatively_ancient_microbes.R**: Returns putatively ancient microbes to visually
  check. Requires mapdamage output for each microbe mapping.
 Runs locally, though could be edited to run remotely.
- **Projection-Based-PCA.R**: local script using projection pipeline output to
create PCA of samples.
- **pcangsd-pipeline.R**: local script using pcangsd output to create PCA of samples.
- **plot_admixture.R**: local script, plots NGSadmix output as barplot
- **final_metagenome_plot.R**: local script, creates dot plot from paper. requires
significant editing after the fact for labelling, etc.
- **Create_Sample_Map.R**: map study area and ocean bathymetry, plots reefs cores
were sourced from.


Data_for_scripts contains some of the metadata and other dependencies of the comp
methods folder. You probably shouldn't need these - all of these outputs are created
by the aDNA_remote_workflow pipeline. Note, the metadata from Kitchen et al (2019)
had some typos. If recreating that plot, use the Kitchen_meta_edit.csv for a fixed
version.

Lab Methods contains the published protocols used to extract DNA and UDG-half
treat aDNA libraries. 

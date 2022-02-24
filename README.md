# Location of Coral aDNA Pipeline Files

Comp Methods folder is the location of all processing scripts used.
 - aDNA_remote_workflow.sh: software run remotely on supercomputing cluster to
 create required files.
 - chr_namer.R, cmdJack_Scaffold.R: Rscripts to be run remotely used in the
  workflow. These create jackknifed files at the scaffold level for the
  projection-based pipeline.
- putatively_ancient_microbes: Returns putatively ancient microbes to visually
  check. Requires mapdamage output for each microbe mapping.
 Runs locally, though could be edited to run remotely.
- 

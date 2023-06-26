# Spatial Discriminability Analysis
Code for Comms Bio paper, Example processed spikes from single subject is provided.

Dependencies from other toolboxes:
* cSPIKE (https://thomas-kreuz.complexworld.net/source-codes/spiky)
* npy-matlab (https://github.com/kwikteam/npy-matlab)
* spikes (https://github.com/cortex-lab/spikes)
* sortingQuality toolbox (https://github.com/cortex-lab/sortingQuality)
* sigstar (https://github.com/raacampbell/sigstar/tree/master)

-------------
# Main scripts to run

To get discriminability from spiking data:
- spgrid_loop_analysis (runs spgrid_Phy_to_results for each subject)

For population analysis:
- finalPopAnalysis (experimental group)
- aggregateControlData (control group)
- doColdspotStats (experimental group)

For supplementary results (comparing PV-Arch and control data):
- compareLaserOnsetDiff
- compareArchvsControlPerf


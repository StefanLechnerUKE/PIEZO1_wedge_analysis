# PIEZO1 MINFLUX analysis
 Matlab scripts and datasets for "State-dependent binding of the wedge domain controls inactivation of the mechanosensitive ion channel PIEZO1"

by Stefan Lechner, Clement Verkest, Lucas Roettger and Nadja Zeitzschel 
Contact: s.lechner@uke.de / c.verkest@uke.de

This set of Matlab scripts are used to visualize and explore PIEZO1 DNA-PAINT Minflux data presented in an upcoming article (Roettger et al., 2026)
The scripts and data are provided for academic and visualization purposes only. Commercial usage of the provided Minflux data (reproduction outside of the publication, etc. ) is forbidden.
Please contact the authors for further inquiries.


# Requirements

Scripts generated and tested with Matlab R2023b and R2024a. Require additional Matlab features such as the Signal_Processing_Toolbox or Image_Processing_Toolbox

Required Main scripts:

-   PIEZO1_ALFA_wedge_analysis_2.m
-   PIEZO1_ALFA_TrimerInPlaneProjection_wedge.m
-   FitDistributionInterbladeDist_wedge.m

and associated subroutines:
-   DBSCAN
-   dbscan2
-   PlotRawData
-   CalcTrimerAngle1
-   CalculateTraceMean
-   PIEZO1Superparticle
-   PlotTrimerAnalysisResult
-   SetALFAanalysisParameters_V2


Raw data provided:

-   confocal images as tif files, with region selected for Minflux scan as ROI
-   Matlab structure files with raw, unfiltered minflux localizations and traces

	->PIEZO1mGL_CTL_RawData, PIEZO1_ALFAmGL_WedgeDel_RawData. Contain data array (X rows, 8 cols) with raw loc and traces, organized as follows:
	col1: X coord /col2: Y coord / col3: Z coord / col4: TID (trace identification number) / col5: efo (Hz) / col6: cfr / col7: background (Hz) / col8: time (s)
	

# How to use 

-To install/run the code deposited here, download all folder and subfolders and add them to your Matlab path. Then open the files 'PIEZO1_ALFA_wedge_analysis_2.m' and 'PIEZO1_ALFA_TrimerInPlaneProjection_wedge.m'


 
-Run PIEZO1_ALFA_wedge_analysis_2.m	

	-This script is used to plot most panels in figure 4, and some supplementary Data Fig.
	-Select options in 'CHOOSE OPTIONS' section
	-Select the data source in the 'load data' section (0 = PIEZO1 Control, 1 = PIEZO1 WedgeDel)
	-It will generate a graph with some of the plots (interblade distance etc.)
	-The main results are stored in a structure called "results"

-Run PIEZO1_ALFA_TrimerInPlaneProjection_wedge.m, only after running PIEZO1_ALFA_wedge_analysis_2.m to plot the selected trimers displayed in the figure ('IndivTrimersRAW' is required for this script to work)

-Run FitDistributionInterbladeDist.m after running PIEZO1_ALFA_wedge_analysis_2.m to perform multi gauss fit on the interblade distances distribution ('results' structure is required)



# Citation
=======
Ongoing publication in Nature Communications
otherwise:
Roettger et al., 2026, biorxiv
"State-dependent binding of the wedge domain controls inactivation of the mechanosensitive ion channel PIEZO1"
doi: https://doi.org/10.64898/2026.01.12.699060

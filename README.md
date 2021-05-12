# Congyu Luo's 2021 degree project.
Here recorded the code and files related to the 2021 master's thesis <Cellular deconvolution of spatial transcriptomics in ALS mice>.


-Introduction


Amyotrophic lateral sclerosis(ALS) is a progressive neurodegenerative disease that affects motor neurons causing loss of muscle control, eventually causes muscular atrophy, paralysis and death.

Currently, researches main focused on nerve cells, and more evidence has showed that non-neuronal cells also play a role in the early stages of disease progressionthis. This project was aimed to find out which non-neural cell types were early responders before nerve cells degenerate during onset and progression of ALS. And also used single-cell transcriptomics to spatio-temporal interpret a large number of gene expression values.

The source spatial transcription data were from Maniatis's lab[1], which the data were derived from the spinal cords of two mouse models, and were classified according to 11 anatomical areas and 4 time points. This project mainly used EWCE package[2] and spatial transcriptomics methods.


-Overview


Input files folder: 
- 'aav9776_data' folder: Included tables were separated by 4 time points from differential expression results between 2 mouse models (SOD1-WT and SOD1-G93A) of Maniatis's lab. Each table lists the Bayesian factor, posterior means and standard deviation of the β distribution for comparison.

Output files folder: 
- 'BootstrapPlots_aav' folder: Included figures were plotted by '250genes_bootstrap_plot.R' script. These figures included bootstrap diagrams and corresponding Q-Q diagrams of 250 up-regulated/down-regulated genes in 10 cell types at 11 anatomical sites on 4 time points. 

- ' EWCE input files' folder: Included Excels were the output files of 'ewce_prepare.R' script. The file consisted of data statistics for 11 anatomical areas on 4 time points and used as input files for bootstrap diagrams plotting and EWCE analysis.

-' EWCE results' folder: Included the results of the EWCE analysis (figures and tables) by 'EWCE_process.R' script. Results were classified by anatomical areas and time points and also used as inputs of 'timeplot_11.R' script.

Python scripts folder:
-'ST_plot.ipynb' script: For plotting the distribution of genes expression in the spinal cord, the Jupyter notebook is required and enter the gene name in the target gene list.

R scripts folder:
- Contained the relevant R script: '250genes_bootstrap_plot.R', 'EWCE_process.R', 'ewce_prepare.R', 'timeplot_11.R', 'Bootstrap_plots_aav.R', and 'ewce_expression_data_aav.r'. Where 'Bootstrap_plots_aav.R', and 'ewce_expression_data_aav.r' are downloaded from NathanSkene's github and modified. And their functions called by '250genes_bootstrap_plot.R'.


Shell files folder:
- Contained the relevant shell files which used to get file names which names based on 11 anatomic areas and 4 time points, in order to let 'timeplot_11.R' and '250genes_bootstrap_plot.R' read file names automatically instead of manually modified the input file names.


-Citation
1. Maniatis, S. and Äijö, T.Spatiotemporal dynamics of molecular pathology in amyotrophic lateral sclerosis. Science (2019). doi: 10.1126/science.aav9776
2. Skene, N. & Grant, S. Identification of vulnerable cell types in major brain disorders using single cell transcriptomes and expression weighted cell type enrichment. Frontiers in Neuroscience (2016). doi:10.3389/fnins.2016.00016


Steps for running wgcna + documentation of beta coefficients picked for each tissue
Remember to also adjust the paths in each script

1. Soft thresholding
	- run soft_thresholding_daniel.R to generate values and plots
	- run soft_thresholding.sh to confirm values from the plots (if deemed necessary)
	- decide on coefficient to use (see below section on choosing coefficients)
2. Generate adjacency tom
	- Edit adjacency_tom.sbatch to use the coefficient decided on in the previous soft thresholding step
	- run adjacency_tom.sbatch
3. Run wgcna_post_processing.sh

Beta Coefficients - look at 'WGCNA soft thresholding (no design)' powerpoint for how they were chosen

Some examples of beta coefficients (see the powerpoint described above for all tissues):
breast_mammary_tissue: 6
esophagus_mucosa: 6
liver: 7
lung: 6
nerve_tibial: 5
whole_blood: 8
 

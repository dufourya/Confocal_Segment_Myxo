# segment_myxo_confocal

Code for segmentation and analysis of confocal images of Myxococcus xanthus mound formation.

Citation:  
**Cell density, alignment, and orientation correlate with C-signal–dependent gene expression during Myxococcus xanthus development**  
**Y Hoang, Joshua L. Franklin, Yann S. Dufour, Lee Kroos**  
**Proceedings of the National Academy of Sciences Nov 2021, 118 (45) e2111706118; DOI: 10.1073/pnas.2111706118**  

For the analysis code to run properly the folder structure including empty folder needs to be recreated:  
   .\csv  
   .\csv\arrangement  
   .\csv\correlations  
   .\csv\fluorescence_hypotheses  
   .\csv\fluorescence_hypotheses\max_distances  
   .\csv\fluorescence_hypotheses\max_distances_total_fluor  
   .\csv\hypothesis  
   .\csv\hypothesis\figure_3  
   .\csv\hypothesis\figure_3\maximum_radial_proportions  
   .\csv\hypothesis\figure_3\val  
   .\csv\raw_data  
   .\figure_2  
   .\figure_3  
   .\figure_4  
   .\figure_5  
   .\figure_supplemental  
   .\rds  
   .\rds\model_rds  
   .\rds\model_rds\radius_angle_align  
   .\rds\model_rds\radius_density  
   .\rds\model_rds\radius_intensity  
   .\rds\model_rds\radius_label  
   .\val  

The data tables in '.\raw_data' contain object parameters generated from the segmentation and analyses of confocal fluorescence imaging of cell mounds, such as size, shape, position, orientation, or fluorescence intensity. Each file represents one confocal image stack.  

| Fruiting body | Strain description |
| ------ | ------ |
| 1,2,3 | YH8; Pvan-tdTomato; collected with setting optimized for each timepoints (Fig. 2,3, S4, S5) |
| 13,14,15 | YH7; Pvan-mNeongreen; collected with the same setting over time (Fig. S9) |
| 16,17,18 | YH8; Pvan-tdTomato; collected with the same setting over time (Fig. S9) |
| 19,20,21 | YH14; Pdev-tdTomato + Pvan-mNeongreen (Fig. 2, 3, 4, 5, S4, S5, S8) |
| 22,23,24 | YH15; PfmgE-tdTomato + Pvan-mNeonGreen (Fig. 2, 3, 4, 5, S4, S5, S8) |

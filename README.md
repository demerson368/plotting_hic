# README #

Plotting example using lib5c package

### How do I get set up? ###

sh setup_venv.sh

### Guidelines ###

### Step 1 ###
plotting example for domain calls on heatmap:

run plot_domain.py    
(example input provided in plot_domain.sh. Run as "sh plot_domain.sh")    

Prior to running this code raw Hi-C matrices should be balanced. Chromosome balanced sparse matrix is
used as input.   

#### input ####

heatmap = chromosome sparse matrix .npz in input_heatmaps directory     
domains = domain calls bed file in input_bed directory     
chromosome = chromosome of plot      
start = start coordinate of plot     
end = end coordinate of plot      
resolution = resolution basepair of heatmap   
max_color = maximum color of plot     

#### Output ####
output placed in output/ subdirectory    
*{domains}*'_' + args.chr  + '_' + str(start)  + '_' + str(end) + '_' +  str(args.max_color) + '.png'


![Scheme](output/plots_domain_calls_example_chr21_chr21_30000000_33000000_100.0.png)


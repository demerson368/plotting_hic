import matplotlib as mpl
mpl.use('Agg')

from lib5c.plotters.extendable import ExtendableHeatmap
#from lib5c.parsers.genes import load_gene_table
#from lib5c.contrib.pybigwig.bigwig import BigWig
import lib5c.plotters
import scipy
import numpy as np
import argparse
import glob

def load_PCA(PCA_file):

    """Assume that PCA file is 4 column .bed file
       
    """

    PCA_x = {}
    PCA_y = {}

    input = open(PCA_file, 'r')

    for line in input:
        if line.startswith("chrom"):
            continue
        else:
            chrom = str(line.strip().split('\t')[0])
            start = int(line.strip().split('\t')[1])
            end = int(line.strip().split('\t')[2])
            value = float(line.strip().split('\t')[3])
            if chrom not in PCA_x:
                PCA_x[chrom] = []
                PCA_y[chrom] = []
            PCA_x[chrom].append(start)
            PCA_y[chrom].append(value)
    input.close()

    return PCA_x,PCA_y


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('heatmap_dir',type=str,help="heatmap dir")
    parser.add_argument('heatmap_name',type=str,help="name of heatmap")
    parser.add_argument('PCA_track',type=str,help="PCA track")
    parser.add_argument('min_PCA',type=float,help="minimum PCA value")
    parser.add_argument('max_PCA',type=float,help="maximum PCA value")
    parser.add_argument('chr', type=str, help = "chromosome")
    parser.add_argument('start', type=int, help = "start")
    parser.add_argument('end', type=int, help = "end")
    parser.add_argument('resolution',type=int, help = "resolution")
    parser.add_argument('max_color',default = 100,type=float,help="maximum color scale")

    args = parser.parse_args()

    total_heatmap_file = args.heatmap_dir + args.heatmap_name

    matrix_pre = scipy.sparse.load_npz(total_heatmap_file)
    matrix_pre_csr = matrix_pre.tocsr()
    size = matrix_pre_csr.shape[0]

    PCA_x,PCA_y = load_PCA(args.PCA_track)

    PCA_pos_signal = np.array(PCA_y[args.chr]).copy()
    PCA_neg_signal = np.array(PCA_y[args.chr]).copy()

    PCA_pos_signal[PCA_pos_signal <= 0] = 0 #np.nan
    PCA_neg_signal[PCA_neg_signal > 0] = 0 #np.nan

    max_color = float(args.max_color)  #-10*np.log2(args.max_color)

    start =  args.start
    end = args.end

    
    matrix = matrix_pre_csr[int(start/args.resolution):int(end/args.resolution),int(start/args.resolution):int(end/args.resolution)].todense()
    some_square_matrix = np.triu(matrix) + np.triu(matrix).T - np.diag(np.diag(matrix))

    h = ExtendableHeatmap(array=some_square_matrix,grange_x={'chrom': args.chr, 'start': start,'end': end},colorscale=(0, max_color),colormap='Reds')

    width = 1
    h.add_ax('PCA')
    h['PCA'].bar(np.array(PCA_x[args.chr])/args.resolution,PCA_pos_signal,width,color='g')
    h['PCA'].bar(np.array(PCA_x[args.chr])/args.resolution,PCA_neg_signal,width,color='r')
    h['PCA'].set_xlim((start/args.resolution,end/args.resolution))
    h['PCA'].set_ylim((args.min_PCA,args.max_PCA))
    h['PCA'].set_xticklabels([])
    h['PCA'].axis('off')

    h.save(args.heatmap_name[:-4] + '_' + args.chr  + '_' + str(start)  + '_' + str(end) + '_' +  str(args.max_color) + '.png')


if __name__ == "__main__":
    main()
     

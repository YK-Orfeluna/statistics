# coding: utf-8

import argparse
import logging
from os import makedirs
from os.path import splitext
from pdb import set_trace

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from factor_analyzer import FactorAnalyzer
from factor_analyzer.factor_analyzer import calculate_bartlett_sphericity, calculate_kmo

ROTATIONS = ("varimax", "promax", "oblimin", "oblimax", "quartimin", "quartimax", "equamax")
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("inf", type=str, \
        help="Input file that is csv or tsv only.")
    parser.add_argument("outd", type=str, \
        help="Output files will saves in $outd.")
    parser.add_argument("N_adj", type=int, \
        help="Number of adjective pairs.")
    parser.add_argument("--inf_header", action="store_true", default=False, \
        help="If $inf has header, to activate this argument.")
    parser.add_argument("--inf_index", action="store_true", default=False, 
        help="If $inf has index, to activate this argument.")
    parser.add_argument("--inf_drop_NA", action="store_true", default=False, \
        help="If activate, 'NA's in $inf are dropped from $inf.")
    parser.add_argument("--p_alpha", type=float, default=0.05, \
        help="Significance level for chi square test; default is 0.05.")
    parser.add_argument("--FA_rotation", type=str, choices=ROTATIONS, default=None, \
        help="To choose rotation type of FactorAnalysis; default is None.")
    parser.add_argument("--FA_method", type=str, choices=("minres", "ml", "principal"), default="minres", \
        help="To choose method of FactorAnalysis; default is 'minres'.")
    parser.add_argument("--FA_n_factor", type=int, default=3, \
        help="Number of factors in FactorAnalysis; default is 3. However, if activate $kaiser_guttman, this argument is ignored.")
    parser.add_argument("--kaiser_guttman", action="store_true", default=False, \
        help="If activate, number of factors is decided based on Kaiser-Guttman Method ($FA_n_factor is ignored.)")
    parser.add_argument("--fig_ext", type=str, choices=("png", "pdf", "eps"), default="png", \
        help="To choose output figures' extention; default is 'png'.")
    parser.add_argument("--fig_dpi", type=int, default=300, \
        help="Output figures' dpi; default is 300.")
    parser.add_argument("--heatmap_figsize_x", type=int, default=5, \
        help="Heatmap's x-size; default is 5.")
    parser.add_argument("--heatmap_figsize_y", type=int, default=6, \
        help="Heatmap's y-size; default is 6.")
    parser.add_argument("--heatmap_sort_factor", type=int, default=1, \
        help="To choose sort criteria of heatmap; default is 1 (i.e., Factor 1).")
    parser.add_argument("--heatmap_cmap", type=str, default=None, \
        help="To choose heatmap's color map. To see matplotlib's color map.")
    return parser.parse_args()

class Factor_Analysis:
    def __init__(self, inf, outd, N_adj, \
                    inf_header=False, \
                    inf_index=False, \
                    inf_dropna=False, \
                    N_factor=3, \
                    kaiser_guttman=False, \
                    method="minres", \
                    rotation=None, \
                    p_alpha=0.05, \
                    dpi=300, \
                    fig_ext="png", \
                    heatmap_cmap=None, \
                    heatmap_size=(5, 6), \
                    heatmap_sort_factor=1):
        self.outd = outd
        self.logger = self.logger_init()
        self.dataset = self.read_dataset(inf, inf_header, inf_index, dropna=inf_dropna)
        self.N_adj = N_adj
        self.fig_dpi = dpi
        self.fig_ext = fig_ext
        self.p_alpha = p_alpha
        self.rotation = rotation
        self.method = method
        self.N_factor = N_factor
        self.kaiser_guttman = kaiser_guttman
        self.loadigs = None
        self.heatmap_cmap = heatmap_cmap
        self.heatmap_size = heatmap_size
        self.heatmap_sort_factor = heatmap_sort_factor
        
    def logger_init(self):
        logger = logging.getLogger("Factor_Analysis")
        logger.setLevel(logging.DEBUG)
        fh = logging.FileHandler("%s/rslt.txt" %self.outd)
        logger.addHandler(fh)
        sh = logging.StreamHandler()
        logger.addHandler(sh)
        return logger

    def read_dataset(self, inf, header=False, index=False, dropna=False):
        ext = splitext(inf)[1].lower()
        if ext==".csv":
            sep = ","
        elif ext==".tsv":
            sep = "\t"
        else:
            assert False, "This script supports only csv or tsv file as $inf."
        header = 0 if header else None
        index = 0 if index else None
        df = pd.read_csv(inf, sep=sep, header=header, index_col=index)
        if dropna:
            df.dropna(inplace=True)
        self.logger.debug("Successfull: %s is read." %inf)
        return df

    def bartlett_sphericity(self):
        chi_square_value,p_value=calculate_bartlett_sphericity(self.dataset)
        out = "chi^2 = %.3f" %chi_square_value
        if p_value<0.001:
            out += ", p<0.0001"
        else:
            out += ", p=%.3f" %p_value
        self.logger.info(out)
        if p_value<self.p_alpha:
            self.logger.info("It is not an identity matrix.")
        else:
            self.logger.info("It is an identity matrix.")
        return 0

    def KMO_test(self):
        kmo_all, kmo_model=calculate_kmo(self.dataset)
        if kmo_model>=0.5:
            kmo_rslt = "marvelous"
        elif kmo_model>=0.6:
            kmo_rslt = "mediocre"
        elif kmo_model>=0.7:
            kmo_rslt = "middling"
        elif kmo_model>=0.8:
            kmo_rslt = "meritorious"
        elif kmo_model>=0.9:
            kmo_rslt = "marvelous"
        else:
            kmo_rslt = "unacceptable"
        self.logger.info("KMO measure = %.3f (%s)" %(kmo_model, kmo_rslt))

    def KaiserGuttman(self):
        fa = FactorAnalyzer(self.N_adj, rotation=None)
        fa.fit(self.dataset)
        ev, v = fa.get_eigenvalues()
        N_factor = np.where(ev>=1)[0].shape[0]
        # To draw eigen values
        plt.figure()
        plt.plot()
        plt.plot(range(1, ev.shape[0]+1), ev, label="eigen value", marker="o")
        plt.axhline(y=1, color='red', linestyle="dotted", label=r"$y=1$")
        plt.legend()
        plt.xlim(0.5, ev.shape[0]+0.6)
        plt.xlabel("factor")
        plt.ylabel("eigen value")
        outf = "%s/eigen_value.%s" %(self.outd, self.fig_ext)
        plt.savefig(outf, dpi=args.fig_dpi)
        self.logger.info("%s is saved." %outf)
        return N_factor

    def factor_analysis(self):
        fa = FactorAnalyzer(n_factors=self.N_factor, rotation=self.rotation, method=self.method) 
        score = fa.fit_transform(self.dataset)
        header = ["Factor_%s" %i for i in range(1, self.N_factor+1)]
        ### 因子負荷量
        self.loadings = fa.loadings_
        outf = "%s/factor_loadings.tsv" %self.outd
        df = pd.DataFrame(fa.loadings_, columns=header)
        df.to_csv(outf, sep="\t", index=False)
        self.logger.info("Facotr loadings are saved as %s." %outf)
        ### 因子得点
        outf = "%s/factor_score.tsv" %self.outd
        df = pd.DataFrame(score, columns=header)
        df.to_csv(outf, sep="\t", index=False)
        self.logger.info("Facotr scores are saved as %s." %outf)
        return 0

    def draw_loadings(self):
        if self.heatmap_sort_factor>self.N_factor:
            axis = self.N_factor - 1
            self.logger.debug("$heatmap_sort_factor is adjusted to %d." %self.N_factor)
        else:
            axis = max(0, self.heatmap_sort_factor - 1)
        sort = np.argsort(self.loadings[:, axis])[::-1]
        df = pd.DataFrame(self.loadings[sort], columns=["Factor_%s" %i for i in range(1, self.N_factor+1)])
        fig, ax = plt.subplots(figsize=(self.heatmap_size))
        sns.heatmap(df, cmap=self.heatmap_cmap, annot=True, vmax=1.0, vmin=-1.0, center=0.0, yticklabels=sort, fmt=".2f")
        ax.set_ylim(self.N_adj, 0)
        plt.ylabel("Adj. pairs")
        plt.yticks(rotation=0)
        outf = "%s/factor_loadings_heatmap.%s" %(self.outd, self.fig_ext)
        plt.savefig(outf, dpi=self.fig_dpi)
        self.logger.info("Factor loadings are saved as %s." %outf)
        return 0

    def main(self):
        self.bartlett_sphericity()
        self.KMO_test()
        if self.kaiser_guttman:
            self.N_factor = self.KaiserGuttman()
        self.factor_analysis()
        self.draw_loadings()
        return 0

    def __call__(self):
        self.main()

if __name__=="__main__":
    args = get_args()
    makedirs(args.outd, exist_ok=True)
    fa = Factor_Analysis(\
        args.inf, args.outd, args.N_adj, \
        inf_header=args.inf_header, \
        inf_index=args.inf_index, \
        inf_dropna=args.inf_drop_NA, \
        N_factor=args.FA_n_factor, \
        method=args.FA_method, \
        kaiser_guttman=args.kaiser_guttman, \
        rotation=args.FA_rotation, \
        p_alpha=args.p_alpha, \
        dpi=args.fig_dpi, \
        fig_ext=args.fig_ext, \
        heatmap_size=(args.heatmap_figsize_x, args.heatmap_figsize_y), \
        heatmap_sort_factor=args.heatmap_sort_factor, \
        heatmap_cmap=args.heatmap_cmap, \
        )
    fa()

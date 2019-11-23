# -*- coding: utf-8 -*-

import argparse
from os import makedirs
from os.path import splitext, basename

import numpy as np
from scipy.stats import ttest_rel, ttest_ind, ttest_1samp

from ftest import ftest
from pandas import read_csv

def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("inf", help="only CSV or TSV file.")
	parser.add_argument("outd", help="the result is saved to 'outd' (directory).")
	parser.add_argument("--alpha", type=float, default=0.05, \
		help="significance level (default = '0.05').")
	parser.add_argument("--type", choices=["independent", "pair", "one"], default="independent", \
		help="t-test' type: 'independent' = independent t-test, 'pair': paired t-test, 'one': one-sample t-test (default is 'independent').")
	parser.add_argument("--population-mean", type=float, default=0.0, \
		help="to use one-sample t-test, the population mean is set to this argment.")
	parser.add_argument("--effect-size", choices=["d", "g", "r"], default="d", \
		help="type of effect size (default = 'd').")
	parser.add_argument("--g-adjust", choices=["yes", "no"], default="yes", \
		help="if 'yes' and two sample sizes are same, effect size 'g' is adjusted (default = 'yes')")
	return parser.parse_args()

class Ttest :
	def __init__(self, inf, outd, alpha=0.05, type_ttest="independent", ef_flag="d", g_adj="yes", population=0) :
		self.population = population			# 1標本t検定を行うときの母集団の値
		
		self.inf = inf
		self.outd = outd
		self.alpha = alpha
		self.x = np.array([])
		self.y = np.array([])

		self.type_ttest = type_ttest			# 対応のあるt検定を行う場合はTrue

		self.ef_flag = ef_flag					# 効果量を指定する
		self.ef = None
		self.ef_level = ""
		self.cid = []
		self.adj = True if g_adj=="yes" else False

		self.f_value = 0 						# F値
		self.f_p = 1.0							# Ftestのp値（帰無仮説に従い，初期値を1.0としておく）
		
		self.t_value = 0						# T値
		self.t_p = 1.0							# Ttestのp値（帰無仮説に従い，初期値を1.0としておく）
		
	def read_data(self) :						# CSVファイルの読み込み
		if splitext(self.inf)[1].lower()==".csv" :					# csvならそのまま関数実行
			delimiter = ","
		elif splitext(self.inf)[1].lower()==".tsv":
			delimiter = "\t"
		else:
			exit("Error: this script supports CSV or TSV file only.\nYour chose file is not CSV or TSV.")
		df = read_csv(self.inf, header=None, index_col=None, delimiter=delimiter)
		data = df.values.T
		self.x = data[0]
		if self.type_ttest!="one":
			self.y = data[1]
		else:
			self.y = np.zeros_like(self.x)
			self.y.fill(self.population)
		self.mean_x, self.sd_x, self.sem_x = self.compute_stats(self.x)
		self.mean_y, self.sd_y, self.sem_y = self.compute_stats(self.y)

	def compute_stats(self, arr):
		mean = np.mean(arr)
		std = np.std(arr, ddof=1)
		sem = std / arr.shape[0]
		return mean, std, sem

	def ttest(self) :
		if self.type_ttest=="one":    # one-sample t-test
			self.t_value, self.t_p = ttest_1samp(self.x, self.population)
		elif self.type_ttest=="independent":    # independent t-test
			self.f_value, self.f_p = ftest(self.x, self.y)
			equal = False if self.f_p<self.alpha else True
			self.t_value, self.t_p = ttest_ind(self.x, self.y, equal_var=equal)
		elif self.type_ttest=="pair":    # paired t-test
			self.t_value, self.t_p = ttest_rel(self.x, self.y)			
		self.ef, self.ef_level, self.cid = self.effect_size()
		return 0

	def effect2level(self, effect_size) :							# 効果量の解釈を行う
		if self.ef_flag == "r" :
			l, m, s = 0.5, 0.3, 0.1
		elif self.ef_flag in ("d", "g") :
			l, m, s = 0.8, 0.5, 0.2

		if effect_size> l :
			return "Large"				# 大
		elif effect_size> m :
			return "Medium"			# 中
		elif effect_size> s :
			return "Small"				# 小
		else :
			return"Rare"				# ほどんどなし

	def effect_size(self) :		# 効果量を求める（クラスの引数により求める効果量の種類を決定する）
		if self.ef_flag == "r" :
			effect_size = np.sqrt(self.t_value**2 / (self.t_value**2 + (self.x.shape[0] + self.y.shape[0] - 2)))
		elif self.ef_flag in ("d", "g"):
			ddof = 0 if self.ef_flag=="d" else 1
			n1, n2 = self.x.shape[0] - ddof, self.y.shape[0] - ddof
			s1, s2 = np.var(self.x, ddof=ddof), np.var(self.y, ddof=ddof)
			sp = np.sqrt((n1*s1 + n2*s2) / (n1 + n2))
			effect_size = np.abs(self.mean_x - self.mean_y) / (sp + 1e-7)
			if n1==n2 and self.adj and self.ef_flag=="g":
				coef = 1 - (3 / (4 * (n1 + n2) - 9))
				effect_size *= coef
			else:
				self.adj = False
			semd1 = (n1 + n2) / (n1 * n2)
			semd2 = effect_size**2 / ((n1 + n2 - 2) * 2)
			semd = np.sqrt(semd1 + semd2)						# d/gの標準誤差
			cid = [effect_size - 1.96 * semd, effect_size+ 1.96 * semd]	# d/gの95%信頼区間
		ef_level = self.effect2level(effect_size)
		return effect_size, ef_level, cid

	def write(self) :											# txtファイルへの書き出し
		n = "\n"
		outf = "%s/%s.txt" %(self.outd, splitext(basename(self.inf))[0])
		out  =""
		with open(outf, "w") as fd :
			if self.type_ttest=="one":
				out += "One-Sample T-test\n"
				out += "population: mean = %.3f\n" %self.population
				out += "sample: mean = %.3f, stdev = %.3f, N = %d\n\n" %(self.mean_x, self.sd_x, self.x.shape[0])
			else :
				if self.type_ttest=="paried":
					out += "Two-sample paired t-test\n"
				else:
					out += "Two-Sample independent t-test\n"
				out += "x: mean = %.3f, stdev = %.3f, N = %d\n" %(self.mean_x, self.sd_x, self.x.shape[0])
				out += "y: mean = %.3f, stdev = %.3f, N = %d\n\n" %(self.mean_y, self.sd_y, self.y.shape[0])
				if self.type_ttest=="independent":
					out += "F-Test: F = %.3f, p = %.3f %s\n" %(self.f_value, self.f_p, "[sig]" if self.f_p<self.alpha else "")
			out += "T-Test: t = %.3f, p = %.3f%s\n" %(self.t_value, self.t_p, "[sig]" if self.t_p<self.alpha else "")
			out += "Effect Sise: %s%s = %.3f [%s]\n" %(self.ef_flag, "[adj]" if self.adj else "", self.ef, self.ef_level)
			if self.ef_flag in ("g", "d"):
				out += "95per. CI of %s: [%.3f, %.3f]\n" %(self.ef_flag, self.cid[0], self.cid[1])
			fd.write(out)
	def main(self) :
		self.read_data()			# ファイルを読み込む
		self.ttest()
		self.write()
	def __call__(self):
		self.main()

if __name__ == "__main__" :
	args = get_args()
	makedirs(args.outd, exist_ok=True)
	ttest = Ttest(args.inf, args.outd, args.alpha, args.type, args.effect_size, args.g_adjust, args.population_mean)
	ttest()
	exit("done: process")


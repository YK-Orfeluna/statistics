# -*- coding: utf-8 -*-

from pandas import read_csv
import numpy as np
from scipy.stats import chi2_contingency, chisquare
from p2ast import *
from os.path import splitext, basename
import codecs

DEBUG = True
CHI = u"\u03c7"
V = "Cram%sr's V" %u"\u00e9"
W = "Cohen's w"

class Chi2 :
	def __init__(self, filename, prob=None) :
		self.filename = filename		# 読み込み対象のファイル名

		self.data = np.array([])		# 実測値
		self.n = 0						# N（sample数）
		self.c = 1						# 表の縦
		self.r = 1						# 表の横

		self.prob = prob				# 理論確立
		if type(prob) != np.ndarray :
			prob = np.array(prob, dtype=np.float64)

		self.chi2 = 0.0					# カイ二乗値
		self.p = 1.0					# p値（帰無仮説に従い，初期値は1.0とする）
		self.df = 0						# 自由度
		self.expected = np.array([])	# 期待値

		self.ef = None

	def read_data(self) :								# データの読み込み
		if splitext(self.filename)[1] == ".csv" :
			df = read_csv(self.filename, index_col=None, header=None)
		elif splitext(self.filename)[1] == ".tsv" :
			df = read_csv(self.filename, index_col=None, header=None, delimiter="\t")
		else :
			exit("Error: This script only suppot CSV or TSV file.")

		self.data = df.values
		self.n = np.sum(self.data)						# サンプル数

		if self.data.shape[0] > 1 :
			self.r = self.data.shape[0]

		self.c = self.data.shape[1]

		if DEBUG :
			print("data: \n", self.data)


	def test(self) :									# 検定実施
		if self.c == 1 :								# 1 x j　のカイ二乗検定
			if self.prob == None :						# 理論確立の入力がなかった場合，自動で等分になるように理論確立を設定する
				self.prob = np.ones(self.r)

			self.expected = self.prob * self.n			# 理論確立とNから期待値を求める

			self.chi2, self.p = chisquare(self.data, self.expected)
			self.df = self.r - 1
	
		else :											# i x j のカイ二乗検定
			self.chi2, self.p, self.df, self.expected = chi2_contingency(self.data)

		self.effect_size()								# 効果量を求める

		if DEBUG :
			print("%s^2(%d): %f" %(CHI, self.df, self.chi2))
			print("p-value: %f" %self.p + p2ast(self.p))
			print("expected frequencies: \n", self.expected)
			if self.c == 1 :
				print("effect size %s: %f" %(W, self.ef))
			else :
				print("effect size %s: %f" %(V, self.ef))

	def effect_size(self) :
		if self.c == 1 :
			p0 = self.prob
			p1 = self.data / self.data.sum()

			self.ef = np.sqrt(np.sum((p0 - p1)**2 / p0))							# Cohen's w

		else :
			self.ef = np.sqrt(self.chi2 / (min(self.r-1, self.c-1) * self.n))		# Cramer's V

	def write(self) :
		n = "\n"
		outname = splitext(self.filename)[0]
		outname += "_chi2.txt"

		with codecs.open(outname, "w", "utf-8") as fd :
			fd.write("Chi-squared test" + n + n)

			fd.write("Measured value: " + n)
			fd.write("%s" %self.data + n + n)

			fd.write("expected frequencies: " + n)
			fd.write("%s" %self.expected + n + n)

			fd.write("%s^2(%d): %f" %(CHI, self.df, self.chi2) + n)
			fd.write("p-value: %f" %self.p + p2ast(self.p) + n + n)
			
			if self.c == 1 :
				fd.write("effect size %s: %f" %(W, self.ef))
			else  :
				fd.write("effect size %s: %f" %(V, self.ef))

	def run(self) :
		self.read_data()

		self.test()

		self.write()

		
if __name__ == "__main__" :
	from sys import argv
	if len(argv) < 1 :
		exit("Error: arg1 is filename(CSV or TSV).")
	elif len(argv) == 2 :
		prob = None
	elif len(argv) == 3 :
		prob = argv[2]
	elif len(argv) > 3 :
		exit("Error: many args.")
	
	app = Chi2(argv[1], prob)
	app.run()
	

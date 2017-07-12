# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency, chisquare
from p2ast import *

DEBUG = True
CHI = u"\u03c7"
V = "Cram%sr's V" %u"\u00e9"

class Chi2 :
	def __init__(self, filename, prob=None) :
		self.filename = filename		# 読み込み対象のファイル名

		self.data = np.array([])		# 実測値
		self.n = 0						# N（被験者数）
		self.c = 1						# 表の縦
		self.r = 1						# 表の横

		if prob != None :				# 1*c
			self.prob = prob
			if type(prob) != np.ndarray :
				prob = np.array(prob, dtype=np.float64)

		self.chi2 = 0.0					# カイ二乗値
		self.p = 1.0					# p値
		self.df = 0						# 自由度
		self.expected = np.array([])	# 期待値

		self.ef = None

	def read_data(self) :
		df = pd.read_csv(self.filename, index_col=0)

		self.data = df.values
		self.n = np.sum(self.data)

		self.c = self.data.shape[0]
		if len(self.data.shape) != 1 :
			self.r = self.data.shape[1]

		print("data: ")
		print(self.data)

	def test(self) :
		if self.c == 1 :
			if self.prob == None :
				self.prob = np.ones(self.r)
			self.prob /= np.sum(self.prob)			# 理論確立
			self.expected = prob * self.n

			self.chi2, self.p = chisquare(self.data, self.expected)
			self.df = self.r - 1
	
		else :
			self.chi2, self.p, self.df, self.expected = chi2_contingency(self.data)
			self.effect_size()

		if DEBUG :
			print("%s2(%d): %f" %(CHI, self.df, self.chi2))
			print("p-value: %f" %self.p + p2ast(self.p))
			print("expected frequencies: ")
			print(self.expected)
			if self.ef != None :
				print("effect size %s: %f" %(V, self.ef))

	def effect_size(self) :
		self.ef = np.sqrt(self.chi2 / (min(self.r-1, self.c-1) * self.n))

	def write(self) :
		n = "\n"
		outname = self.filename.rstrip(".csv")
		outname += "_chi2.txt"

		with open(outname, "w") as txt :
			txt.write("Chi-squared test" + n + n)

			txt.write("Measured value: " + n)
			txt.write("%s" %self.data + n + n)

			txt.write("expected frequencies: " + n)
			txt.write("%s" %self.expected + n + n)

			txt.write("chi^2(%d): %f" %(self.df, self.chi2) + n)
			txt.write("p-value: %f" %self.p + p2ast(self.p) + n + n)
			
			if self.ef != None :
				txt.write("effect size Crame'r's V: %f" %self.ef)

	def run(self) :
		self.read_data()

		self.test()

		self.write()

		
if __name__ == "__main__" :
	app = Chi2("data_chi2.csv")
	app.run()
	

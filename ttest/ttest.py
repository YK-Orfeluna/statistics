# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
from scipy.stats import ttest_rel, ttest_ind, ttest_1samp, sem
from ftest import ftest
import pandas as pd

DEBUG = True
#DEBUG = False

def p2ast(p) :
		if p < 0.001 :
			ast = "****"
		elif p < 0.005 :
			ast = "***"
		elif p < 0.01 :
			ast = "**"
		elif p < 0.05 :
			ast = "*"
		else :
			ast = "nan sig."

		return ast

class Ttest :
	def __init__(self, filename, pair=True, ef_flag="r", sample=2, population=0) :
		self.sample = sample					# t検定の標本数
		self.population = population			# 1標本t検定を行うときの母集団の値
		
		self.filename = filename
		self.x = np.array([])
		self.y = np.array([])

		self.pair = pair						# 対応のあるt検定を行い場合はTrue

		self.ef_flag = ef_flag
		self.ef = None
		self.ef_level = ""
		self.cid = []
		self.adj = False

		self.equal = True						# 等分散かどうか（帰無仮説に従い，初期値をTrueとしておく）

		self.f_value = 0 						# F値
		self.f_p = 1.0							# Ftestのp値（帰無仮説に従い，初期値を1.0としておく）
		
		self.t_value = 0						# T値
		self.t_p = 1.0							# Ttestのp値（帰無仮説に従い，初期値を1.0としておく）
		
	def read_data(self) :						# CSVファイルの読み込み
		df = pd.read_csv(self.filename, index_col=None)
		if DEBUG :
			print(df)

		data = df.values.T
		self.x = data[0]

		if self.sample == 2 :
			self.y = data[1]
		elif self.sample == 1 :
			self.y = self.x.copy()
			self.y.fill(self.population)

		self.data_features()

	def data_features(self) :
		self.meanx = np.mean(self.x)			# 平均
		self.meany = np.mean(self.y)

		self.varx = np.var(self.x)				# 分散
		self.vary = np.var(self.y)

		self.sdx = np.std(self.x, ddof=1)		# SD（標準偏差（不偏））
		self.sdy = np.std(self.y, ddof=1)

		self.semx = sem(self.x)					# SEM（平均の標準誤差）
		self.semy = sem(self.y)

		self.dfx = self.x.shape[0] - 1			
		self.dfy = self.y.shape[0] - 1
		self.df = self.dfx + self.dfy			# 自由度

	def ttest(self) :
		if self.sample == 1 :					# 1標本t検定（サンプルデータがある値より有意に大きいか小さいか）
			self.t_value, self.t_p = ttest_1samp(self.x, self.population)

		elif self.pair :						# 対応のあるt検定
			self.t_value, self.t_p = ttest_rel(self.x, self.y)

		else :									# 対応のないt検定（F検定により等分散か否か調査してから，その結果に対応したt検定を実行する）
			self.f_value, self.f_p = ftest(self.x, self.y)
			if self.f_p < 0.05 :
				self.equal = False				# F検定の結果，等分散でない場合，t検定の引数を変更する

			self.t_value, self.t_p = ttest_ind(self.x, self.y, equal_var=self.equal)

		self.effect_size()

		if DEBUG :
			if self.sample == 1 :
				print("One-Sample T-test")
				print("population: %d" %self.population)
				print("data-mean: %f" %self.meanx)

			else :
				print("Two-Sample T-test")
				print("x_mean, y-mean: %f, %f" %(self.meanx, self.meany))

				if self.pair == False :
					print("F-value: ", self.f_value)
					print("p-value of F-test: ", self.f_p, p2ast(self.f_p))

			print("t-value: %f" %self.t_value)
			print("p_value: %f %s" %(self.t_p, p2ast(self.t_p)))

			if self.adj :
				print("effect_size %s(adj): %f" %(self.ef_flag, self.ef))
			else :
				print("effect_size %s: %f" %(self.ef_flag, self.ef))
			if self.ef_flag == "d" or self.ef_flag == "g" :
				print("95per. CI of ", self.ef_flag, ": ", self.cid)
			print("effect size level: ", self.ef_level)

	def level(self) :							# 効果量の解釈を行う
		if self.ef_flag == "r" :
			l, m, s = 0.5, 0.3, 0.1
		elif self.ef_flag == "g" or self.ef_flag == "d" :
			l, m, s = 0.8, 0.5, 0.2

		if self.ef > l :
			self.ef_level = "Large"
		elif self.ef > m :
			self.ef_level = "Medium"
		elif self.ef > s :
			self.ef_level = "Small"
		else :
			self.ef_level ="Rare"

	def effect_size(self) :		# 効果量を求める（クラスの引数により求める効果量の種類を決定する）
		if self.ef_flag == "r" :
			self.ef = np.abs(np.sqrt(self.t_value**2 / (self.t_value**2 + self.df)))

		# elif self.pair :
		# 	self.ef_flag = "d"
		# 	xd = self.x - self.y
		# 	self.ef = np.mean(xd) / np.std(xd, ddof=1)

		elif self.ef_flag == "d" :		# Cohen's d
			sp = np.sqrt(((self.dfx+1) * np.std(self.x)**2 + (self.dfy+1) * np.std(self.y)**2) / (self.df + 2))
			self.ef = np.abs((self.meanx - self.meany) / sp)

		elif self.ef_flag == "g" :		# Hedges' g
			sp = np.sqrt((self.dfx * self.sdx**2 + self.dfy * self.sdy**2) / self.df)
			self.ef = np.abs((self.meanx - self.meany) / sp)
				
			if self.dfx == self.dfy :							# gに対して補正を行う
				correction = 1 - (3 / (4 * (self.df + 2) - 9))
				self.ef *= correction
				self.adj = True

			semd1 = (self.df + 2) / ((self.dfx+1) * (self.dfy+1))
			semd2 = (self.ef**2) / (2 * self.df)
			semd = np.sqrt(semd1 + semd2)						# d/gの標準誤差
			self.cid =[self.ef-1.96*semd, self.ef+1.96*semd]	# d/gの95%信頼区間

		if self.ef != None :
			self.level()

	def write(self) :											# txtファイルへの書き出し
		n = "\n"
		with open("ttest.txt", "w") as txt :
			if self.sample == 1 :
				txt.write("One-Sample T-test" + n)
				txt.write("population: %d" %self.population + n)
				txt.write("sample_mean: %f" %self.meanx + n)
				txt.write("samle_var: %f" %self.varx + n)
				txt.write("sample_SD: %f" %self.sdx + n)
				txt.write("sample SEM: %f" %self.semx + n + n)

			else :
				if self.pair :
					txt.write("Two-sample Paired T-test"+ n)
				else :
					txt.write("Two-Sample Independent T-test" + n)

				txt.write("x_mean: %f" %self.meanx + n)
				txt.write("x_var: %f" %self.varx + n)
				txt.write("x_SD: %f" %self.sdx + n)
				txt.write("x_SEM: %f" %self.semx + n + n)

				txt.write("y_mean: %f" %self.meany + n)
				txt.write("y_var: %f" %self.vary + n)
				txt.write("y_SD: %f" %self.sdy + n)
				txt.write("y_SEM: %f" %self.semy + n + n)

				if self.pair :
					txt.write("F-Test Result" + n)
					txt.write("F-value: %f" %self.f_value + n)
					txt.write("p-valueof F-test: %f" %self.f_p + p2ast(self.f_p) + n + n)

			txt.write("t_value: %f" %self.t_value + n)
			txt.write("p_value of T-test: %f" %self.t_p + p2ast(self.t_p) + n + n)

			if self.adj :
				txt.write("effect size %s(adj): %f" %(self.ef_flag, self.ef))
			else : 
				txt.write("effect size %s: %f" %(self.ef_flag, self.ef))
			if self.ef_flag == "d" or self.ef_flag == "g" :
				txt.write("95per. CI of %s: [%f, %f]" %(self.ef_flag, self.cid[0], self.cid[1]))
			txt.write("effect size level: %s" %self.ef_level)


	def run(self) :
		self.read_data()			# CSVファイルを読み込む

		self.ttest()

		self.write()

if __name__ == "__main__" :
	### Pythonのバージョンに合わせたtkinterのimport
	import sys
	v = sys.version_info[0]
	if v == 2 :
		import Tkinter as tkinter
		import tkMessageBox as messagebox
		import tkFileDialog as filedialog
	elif v == 3 :
		import tkinter
		from tkinter import filedialog
	else :
		exit("*This script only supports Python2.x or 3.x.\nSorry, we can not support your Python.")


	### GUI用のおまじない
	root = tkinter.Tk()
	root.option_add('*font', ('FixedSys', 14))
	fTyp=[('csvファイル','*.csv')]
	iDir='.'


	### ファイル選択
	lb=tkinter.Label(root, text="Chose answer-file",width=20)
	lb.pack()
	filename = filedialog.askopenfilename(filetypes=fTyp,initialdir=iDir)


	t = Ttest(filename, pair=True, ef_flag="g")
	t.run()
	exit("System Exit")


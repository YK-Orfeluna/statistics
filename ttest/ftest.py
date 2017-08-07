# -*- coding: utf-8 -*-

import numpy as np
from scipy.stats import f
from p2ast import *

def ftest(x, y) :							# F検定を行って，等分散かどうかを判定する（p>.05なら等分散）
	dfx = len(x) - 1
	dfy = len(y) - 1

	varx = np.var(x)						# 2つのデータの分散
	vary = np.var(y)

	if varx > vary :						# F値=分散比を計算（大きい方を小さい方で必ず割る）
		f_value = varx / vary
	else :
		f_value = vary / varx

	p_value = f.cdf(f_value,dfx,dfy)		# F検定を実施して，p値を取得

	return f_value, p_value


if __name__ == "__main__" :
	from sys import argv
	from pandas import read_csv
	from os.path import splitext, basename
	from scipy.stats import sem
	import codecs

	if len(argv) != 2 :
		exit("Error: arg is filename(CSV or TSV)")

	if splitext(argv[1])[1] == ".csv" :
		df = read_csv(argv[1], header=None, index_col=None)
	elif splitext(argv[1])[1] == ".tsv" :
		df = read_csv(argv[1], header=None, index_col=None, delimiter="\t")
	else :
		exit("Error: your chosen file is not CSV or TSV.")

	data = df.values.T

	x = data[0]
	y = data[1]

	f_value, p_value = ftest(x, y)

	print("f value: ", f_value)
	print("p value:", p_value, p2ast(p_value))
	
	if p_value < 0.05 :
		print("These data are not Equal Variance")
	else :
		print("These data are Equal Variance")

	with codecs.open(splitext(argv[1])[0]+"_F.txt", "w", "utf-8") as fd :
		n = "\n"

		fd.write("x_mean: %f" %np.mean(x) + n)
		fd.write("x_var: %f" %np.var(x) + n)
		fd.write("x_SD: %f" %np.std(x, ddof=1) + n)
		fd.write("x_SEM: %f" %sem(x) + n + n)

		fd.write("y_mean: %f" %np.mean(y) + n)
		fd.write("y_var: %f" %np.var(y) + n)
		fd.write("y_SD: %f" %np.std(y, ddof=1) + n)
		fd.write("y_SEM: %f" %sem(y) + n + n)

		fd.write("f-value: %f" %f_value + n)
		fd.write("p-value: %f" %p_value + p2ast(p_value) + n)

		if p_value < 0.05 :
			fd.write("These data are not Equal Variance")
		else :
			fd.write("These data are Equal Variance")

	exit("System Exit")

# -*- coding: utf-8 -*-
import numpy as np
from scipy.stats import f

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

def ftest(x, y) :							# F検定を行って，等分散かどうかを判定する（p>.05なら等分散）
	dfx = len(x) - 1
	dfy = len(y) - 1

	varx = np.var(x)						# 2つのデータの分散
	vary = np.var(y)

	if varx > vary :						# 分散比を計算（大きい方を小さい方で必ず割る）
		f_value = varx / vary
	else :
		f_value = vary / varx

	p_value = f.cdf(f_value,dfx,dfy)		# F検定を実施して，p値を取得

	return f_value, p_value


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

	import pandas as pd
	from scipy.stats import sem

	df = pd.read_csv(filename, index_col=None)
	print(df)
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

	with open("ftest.txt", "w") as txt :
		n = "\n"

		txt.write("x_mean: %f" %np.mean(x) + n)
		txt.write("x_var: %f" %np.var(x) + n)
		txt.write("x_SD: %f" %np.std(x, ddof=1) + n)
		txt.write("x_SEM: %f" %sem(x) + n + n)

		txt.write("y_mean: %f" %np.mean(y) + n)
		txt.write("y_var: %f" %np.var(y) + n)
		txt.write("y_SD: %f" %np.std(y, ddof=1) + n)
		txt.write("y_SEM: %f" %sem(y) + n + n)

		txt.write("f-value: %f" %f_value + n)
		txt.write("p-value: %f" %p_value + p2ast(p_value) + n)

		if p_value < 0.05 :
			txt.write("These data are not Equal Variance")
		else :
			txt.write("These data are Equal Variance")

	exit("System Exit")

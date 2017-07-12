import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency

DEBUG = True
CHI = u"\u03c7"
V = "Cram%sr's V" %u"\u00e9"

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

class Chi2 :
	def __init__(self, filename) :
		self.filename = filename

		self.data = np.array([])
		self.n = 0
		self.r = 1
		self.c = 1

		self.chi2 = 0.0
		self.p = 1.0
		self.df = 0
		self.expected = np.array([])

		self.ef = 0.0

	def read_data(self) :
		df = pd.read_csv(self.filename, index_col=0)

		self.data = df.values
		self.n = np.sum(self.data)

		self.r = self.data.shape[0]
		self.c = self.data.shape[1]

		print("data: ")
		print(self.data)

	def test(self) :
		self.chi2, self.p, self.df, self.expected = chi2_contingency(self.data)
		self.effect_size()

		if DEBUG :
			print("%s2(%d): %f" %(CHI, self.df, self.chi2))
			print("p-value: %f" %self.p + p2ast(self.p))
			print("expected frequencies: ")
			print(self.expected)
			print("effect size %s: %f" %(V, self.ef))

	def effect_size(self) :
		self.ef = np.sqrt(self.chi2 / (min(self.r-1, self.c-1) * self.n))


	def run(self) :
		self.read_data()

		self.test()

		
if __name__ == "__main__" :
	app = Chi2("data_chi2.csv")
	app.run()
	

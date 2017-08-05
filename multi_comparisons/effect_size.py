# -*- coding: utf-8 -*-
import numpy as np

class EffectSize :
	def __init__(self, dis=None, meanx=None, meany=None, sem=None, sd=None, N=None, index="r") :
		self.dis = dis
		if dis == None :
			self.dis = np.abs(meanx - meany)

		self.sd = sd
		if sd == None :
			self.sd = sem * np.sqrt(N)


		self.t = np.sqrt(self.dis / (self.sd / N))
		self.df = N * 2 - 2

		self.index = index

		self.ef = 0

	def run(self) :
		if self.index == "r" :
			r = np.abs(np.sqrt(self.t**2 / (self.t**2 + self.df)))
		print("effect size r: ", r)


if __name__ == "__main__" :
	ef = EffectSize(dis=6, sem=2.249, N=33)
	ef.run()
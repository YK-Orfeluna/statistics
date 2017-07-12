# -*- coding: utf-8 -*-

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
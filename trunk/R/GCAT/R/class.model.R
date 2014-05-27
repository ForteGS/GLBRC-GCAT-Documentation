#Copyright 2012 The Board of Regents of the University of Wisconsin System.
#Contributors: Jason Shao, James McCurdy, Enhai Xie, Adam G.W. Halstead, 
#Michael H. Whitney, Nathan DiPiazza, Trey K. Sato and Yury V. Bukhman
#
#This file is part of GCAT.
#
#GCAT is free software: you can redistribute it and/or modify
#it under the terms of the GNU Lesser General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#GCAT is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Lesser General Public License for more details.
#
#You should have received a copy of the GNU Lesser General Public License  
#along with GCAT.  If not, see <http://www.gnu.org/licenses/>.

########################################################################
#                                                                      #
#  <model> class definition and functions. Objects contain equations   #
#  and other information for parameterized growth curve models.        #
#                                                                      #
########################################################################

setClass("model", representation(name = "character",
					  expression = "expression",
					  formula = "formula", 
					  guess = "function"))
# Slots:
#   name - a simple description of the model.
#   expression - an object of class "expression" that evaluates the response (transformed OD) with respect to the variable Time. 
#   formula - same as expression, but with y as the response.
#   guess - a function that computes initial guesses for the parameters given a well object with a valid "screen.data" slot 
#           containing useable OD values and slope estimates


# --------------------------------------------------------------------
# Function to create a new model
model = function(name, expression, formula, guess){
	new("model", name = name, expression = expression, formula = formula, guess = guess)
	}

########################################################################
#     Create the logistic 4-parameter model (the extra fourth para-    #
#      meter allows estimating the baseline as well)                   #
########################################################################

logistic.g =  function(well){
	slopes = slopes(well)
	data = data.from(well)
	growth = data[,2]
	time = data[,1]


	min.index = min(which(growth == min(growth)))
	max.index = max(which(growth == max(growth)))
	g50 = (min(growth) + max(growth)) / 2
	g50.dist = abs(growth - g50) 
	g50.index = min(which(g50.dist == min(g50.dist)))


	#  Get the slopes at min, max, and half-way points
	#  The highest positive slope is used as the halfway slope
	min.slope = slopes[min.index]
	max.slope = slopes[max.index]
	g50.slope = slopes[g50.index]
	
	#  Compute an initial guess 
	a = min(growth)
	b = max(growth)
	t50 = (time[g50.index])
	th = (max(growth) - min(growth)) / (4 * (g50.slope))

	c(a = a, b = b, t50 = t50 - 1, th = th)

	}

logistic.e = expression(a + (b-a)/(1 + exp((t50 - Time)/th)))
logistic.f = formula(y~ a + (b-a)/(1 + exp((t50 - Time)/th)))
logistic = model("logistic 4-par.", logistic.e, logistic.f, logistic.g)

remove(logistic.e, logistic.f, logistic.g)

########################################################################
#     Create the Richards 5-parameter model                            #
########################################################################

richards.g =  function(well){
	slopes = slopes(well)
	data = data.from(well)
	growth = data[,2]
	time = data[,1]


	min.index = min(which(growth == min(growth)))
	max.index = max(which(growth == max(growth)))
	g50 = (min(growth) + max(growth)) / 2
	g50.dist = abs(growth - g50) 
	g50.index = min(which(g50.dist == min(g50.dist)))


	#  Get the slopes at min, max, and half-way points
	#  The highest positive slope is used as the halfway slope
	min.slope = slopes[min.index]
	max.slope = slopes[max.index]
	g50.slope = slopes[g50.index]
	
	#  Compute an initial guess 
	a = min(growth)
	b = max(growth)
	t50 = (time[g50.index])
	th = (max(growth) - min(growth)) / (4 * (g50.slope))

	c(a = a, b = b, t50 = t50 - 1, th = th, c = 1)

	}

richards.e = expression(a + (b-a)*(1 + c * exp((t50 - Time)/th)) ^(-1/c))
richards.f = formula(y~ a + (b-a)*(1 + c * exp((t50 - Time)/th)) ^(-1/c))
richards = model("richards 5-par.", richards.e, richards.f, richards.g)

remove(richards.e, richards.f, richards.g)

########################################################################
#     Create the Gompertz model (might be useful as a                  #
#     limiting case of Richards model when c = 0)                      #
########################################################################

gompertz.g =  function(well){
	slopes = slopes(well)
	data = data.from(well)
	growth = data[,2]
	time = data[,1]


	min.index = min(which(growth == min(growth)))
	max.index = max(which(growth == max(growth)))
	g50 = (min(growth) + max(growth)) / 2
	g50.dist = abs(growth - g50) 
	g50.index = min(which(g50.dist == min(g50.dist)))


	#  Get the slopes at min, max, and half-way points
	#  The highest positive slope is used as the halfway slope
	min.slope = slopes[min.index]
	max.slope = slopes[max.index]
	g50.slope = slopes[g50.index]
	
	#  Compute an initial guess 
	a = min(growth)
	b = max(growth)
	t50 = (time[g50.index])
	th = (max(growth) - min(growth)) / (4 * (g50.slope))

	c(a = a, b = b, t50 = t50 - 1, th = th)

	}

gompertz.e = expression(a + (b-a)/exp(exp((t50 - Time)/th)))
gompertz.f = formula(y~ a + (b-a)/exp(exp((t50 - Time)/th)))
gompertz = model("gompertz 4-par.", gompertz.e, gompertz.f, gompertz.g)

remove(gompertz.e, gompertz.f, gompertz.g)

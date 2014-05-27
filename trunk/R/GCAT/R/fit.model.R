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
#  Fit a parameterized model to the growth data in a well object.      #
#                                                                      #
########################################################################
#
# growth.model and backup.growth.model: objects of class "model" (see class.model.R) to use in non-linear regression
#
# This function will use the function stored in the "guess" slot of <growth.model> to calculate initial guesses for growth.model parameters
#   then it will use the "formula" slot with <nls> to fit a non-linear least squares growth.model to the data. 
# If fitting with <growth.model> fails the function will attempt to use <backup.growth.model> (which should be a simpler growth.model with fewer parameters for best results)
#
# fit.if.no.growth: should the function attempt to fit a well even if there was no growth detected? default is F 
# silent: output back to R console? 

fit.model = function(input.well, growth.model, backup.growth.model = NULL, fit.if.no.growth = F, silent = T){

	# Change all relevant slots to <NA> or blank values
	input.well@model.name = "<NA>"
  input.well@fit.par = list()
  input.well@equation = expression() 
  
  # Get OD vs. time data from well
 	input.data = data.from(input.well, na.rm = T) 
  
  # Skip well if <no.growth> in slot "curve.par" is set to true, and <fit.if.no.growth> is false. 
	if(!fit.if.no.growth & lacks.growth(input.well)){
		input.well@fit.info = "skipped - no growth in well."
		if (!silent)
			cat(plate.name(input.well), well.name(input.well), ":", input.well@fit.info, "\n")
		return(input.well)
		}
	# Skip well if there are fewer than 5 data points left in the analysis. 
	if (length(input.data$Time) < 5){
		input.well@fit.info = "skipped - not enough points."
		if (!silent)
			cat(plate.name(input.well), well.name(input.well), ":", input.well@fit.info, "\n")
		return(input.well)
		}
	
	# Change column headers of input.data to the more general "Time" vs. "y"
  names(input.data) = c("Time", "y")

  # Extract the model formula from <growth.model> (slot "formula")
  # Use the function from slot "guess" to calculate initial guesses for model parameters based on slope estimates in <input.well>
  # Attempt to fit a nonlinear least squares odel using <nls> 
   
	fit = try(nls(formula = growth.model@formula, data = input.data, start = growth.model@guess(input.well)), silent = T)
	
	# If using a Richards 5 parameter model, report an error if the fitted value for c is negative,			
  #   if the baseline is below the plateau, or if the fitted theta is negative. 
	
	if(growth.model@name == "richards 5-par."){
		if(class(fit) == "nls"){
			fit.par = as.list(coef(fit))       
			if (fit.par$c < 0 | fit.par$th < 0 | fit.par$b < fit.par$a)
				class(fit) = "try-error"
			}
		}

  # If no error was reported by the model fitting, report successful fit in slot "fit.info", add coefficients to slot "fit.par", 
  #   add the model equation and name to slots "equation" and "growth.model.name"
	if (class(fit) == "nls"){
		input.well@fit.info = paste("Model fit successfully.")
		input.well@fit.par = as.list(coef(fit))
		input.well@equation = growth.model@expression
		input.well@model.name = growth.model@name
		}
	# If the model fitting returned an error, do the same thing as above except using the backup model. 
	else{
		fit = try(nls(formula = backup.growth.model@formula, data = input.data, start = backup.growth.model@guess(input.well)), silent = T)
    if (class(fit) == "nls"){
      # If using a logistic 4-parameter model, report an error if the fitted baseline is below the plateau, or if the fitted theta is 
      fit.par = as.list(coef(fit))  
      if(backup.growth.model@name == "logistic 4-par." & (fit.par$th < 0 | fit.par$b < fit.par$a))
  			input.well@fit.info = "Model fitting failed."		
 			else{
		    input.well@fit.info = paste("Model fit successfully.")
  			input.well@fit.par = as.list(coef(fit))
  			input.well@equation = backup.growth.model@expression
  			input.well@model.name = backup.growth.model@name
        }	
  		}
 		# If both models failed, report a failure to slot "fit.info." leave the other slots blank or NA. 
		else
			input.well@fit.info = "Model fitting failed."		
		}
  # Output to console
	if (!silent)
		cat(plate.name(input.well), well.name(input.well), ":", input.well@fit.info, growth.model@name, "\n")
	return(input.well)
	}

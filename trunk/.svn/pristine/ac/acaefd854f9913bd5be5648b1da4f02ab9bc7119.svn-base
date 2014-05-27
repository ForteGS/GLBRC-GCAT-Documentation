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
# Functions to calculate various things about wells based on fit model #
#                                                                      #
########################################################################
#
#  Common arguments:
#   fitted.well - should be a well containing the results of <fit.model>, most functions will return NA if well has not been fit yet.
#   unlog - should the value be returned on the linear scale as opposed to the log-transformed scale?
#   constant.added - for returning values on the linear scale, what was the constant added before the log transform?
#   digits - passed to the <round> function, default is no rounding (infinity digits)

well.eval = function(fitted.well, Time = NULL, unlog = T, constant.added = 0){
  ########################################################################
  #   Evaluate estimated OD at any timepoints using the fitted model     #
  ########################################################################

  # If no timepoints are provided, use the ones collected in the experiment itself.
	if(!is.numeric(Time))
		Time = data.from(fitted.well)$Time

  # Attempt to use <eval> with the fitted equation and parameters to get estimates for OD at the given timepoints.
	output = try(eval(fitted.well@equation, fitted.well@fit.par), silent = T)

  # Return values on log scale or linear scale. If OD evaluation failed for any reason, return NULL.
  if (is.numeric(output)){
    if(unlog)
  		return(exp(output) - constant.added)
  	else
  		return(output)
    }
   else
    return(NULL)
	}

model.residuals = function(fitted.well, unlog = F){
  ########################################################################
  #   Evaluate model residuals using the measured vs. fitted OD values   #
  ########################################################################
  # Get measured OD values from well using normal or log scale
  if (unlog)
		fitted.well@use.log = F
	else
		fitted.well@use.log = T
	measured.OD = data.from(fitted.well)[,2]

	# Use <well.eval> with no Time argument to get fitted OD values at measured timepoints.
	predicted.OD = well.eval(fitted.well, unlog = unlog)

	# If all values are valid, return the differences
	if (!is.numeric(predicted.OD))
		return(NA)
	else
    return(measured.OD - predicted.OD)
	}


model.good.fit = function(fitted.well, digits = Inf, unlog = F){
  ########################################################################
  #   Calculate a metric for fit accuracy using squared residuals        #
  ########################################################################

  # Sum of squared residuals
	rss = sum(model.residuals(fitted.well)^2, unlog = unlog)

	if (unlog)
		fitted.well@use.log = F
	else
		fitted.well@use.log = T

  # Variance in x and y of the measured timepoints
	varx = var(data.from(fitted.well)[,1])
	vary = var(data.from(fitted.well)[,2])

  # Negative log(10) of the residual sum of squares divided by variance on either axis.
	return(round(-log10(rss/varx/vary),digits = digits))
	}

parameter.text = function(fitted.well){
  ########################################################################
  #           Output a string with values of fitted parameters           #
  ########################################################################

  # Get a list of fitted parameters
  fit.par = fitted.well@fit.par

  # Return nothing if the list is empty. Otherwise, concatenate the terms in the list with the parameter names.
	if(!is.list(fit.par))
		return()
  else{
  	output = ""
  	i = 1
  	while(i <= length(fit.par)){
  		output = paste(output, names(fit.par)[i], "=", round(as.numeric(fit.par[i]),3), "; ", sep = "")
  		i = i + 1
  		}
  	output
  	}
	}

specific.growth = function(fitted.well, digits = Inf){
  ########################################################################
  #        Calculate specific growth from fitted parameters              #
  ########################################################################

  # If there are no fitted OD values, return NA
  if (is.null(well.eval(fitted.well)))
		return(NA)

	Time = fitted.well@fit.par$t50

  # If the timepoint at which maximum specific growth was reached according to fitted model was before inoculation, also return NA
	if (Time < 0)
		return(NA)

  # Evalute the derivative of the model equation with respect to time, at t50 (maximum slope)
	slope.max = eval(deriv(fitted.well@equation, "Time"), fitted.well@fit.par)
	slope.max = attr(slope.max, "gradient")

  # Return the maximum slope, rounded as specified by <digits> (default: do not round)
	if(is.nan(slope.max[1]) | is.na(slope.max[1]))
		return(NA)
	else
		return(round(slope.max[1], digits))
	}


plateau = function(fitted.well, unlog = T, constant.added = 1.0, digits = Inf){
  ########################################################################
  #        Calculate plateau OD from fitted parameters                   #
  ########################################################################

	plat = fitted.well@fit.par$b
	if (!is.numeric(plat))
		plat = NA
	else{
    if(!unlog)
		  plat = round(plat, digits)
    else
		  plat = round(exp(plat) - constant.added, digits)
		}
	return(plat)
	}

baseline = function(fitted.well, unlog = T, constant.added = 1.0, digits = Inf){
  ########################################################################
  #        Calculate baseline OD from fitted parameters                  #
  ########################################################################

  base = fitted.well@fit.par$a

  # If b (plateau OD) is invalid, return NA.
	if (!is.numeric(fitted.well@fit.par$b))
		base = NA
  # If a (baseline OD) is invalid but plateau OD was valid, return zero.
  else if (!is.numeric(base))
		base = 0
	else{
    if(!unlog)
		  base = round(base, digits)
    else
		  base = round(exp(base) - constant.added, digits)
		}
	return(base)
	}



inoc.OD = function(fitted.well, unlog = T, constant.added = 1.0, digits = Inf){
  ########################################################################
  #        Calculate OD at inoculation from fitted parameters            #
  ########################################################################

  # Evaluated the fitted model at the inoculation timepoint (should be zero from using <start.times> from table2wells.R)
	if (is.null(well.eval(fitted.well)))
		return(NA)
  else{
    return(round(well.eval(fitted.well, 0, unlog = unlog), digits) - constant.added)
    }
	}


end.OD = function(fitted.well, unlog = T, constant.added = 1.0, digits = Inf){
  ########################################################################
  #        Calculate OD at end of experiment from fitted parameters      #
  ########################################################################

  # Evaluated the fitted model at the final timepoint (just the last valid timepoint in the experiment)
	if (is.null(well.eval(fitted.well)))
		return(NA)
  else{
  	time.fin = max(data.from(fitted.well)$Time)
  	return(round(well.eval(fitted.well, time.fin, unlog = unlog), digits) - constant.added)
    }
  }


  ########################################################################
  #   Calculate total growth: plateau minus the inoculated OD            #
  ########################################################################

total.growth = function(fitted.well,...)
	plateau(fitted.well,...) - inoc.OD(fitted.well,...)



reach.plateau = function(fitted.well, cutoff = 0.75, unlog = T){
  ########################################################################
  # Did the curve come close to the plateau OD during the experiment?    #
  ########################################################################

  plat = plateau(fitted.well, unlog=unlog)
  inoc = inoc.OD(fitted.well, unlog=unlog)
  final = end.OD(fitted.well, unlog=unlog)

	if (!is.na(final)){
    # If the plateau is the same as the OD at inoculation, return TRUE
    if ((plat - inoc) == 0)
      return(T)
     # If the difference between the final OD and inoculation OD is at least a certain proportion
     #  <cutoff> of the difference between the plateau and inoculated ODs, return TRUE.
    else
      return((final - inoc) / (plat - inoc) > cutoff)
		}
	else
		return(T)
		# If no final OD was calculated (if curve was not fit properly) just return T.
	}


lag.time = function(fitted.well, digits = Inf){
  ########################################################################
  #              Calculate the lag time from the fitted OD               #
  ########################################################################

  # Return NA if the well has no fitted OD values or a specific growth could not be calculated.
  if (is.null(well.eval(fitted.well)))
		return(NA)
	max.slope = specific.growth(fitted.well)
	if (is.na(max.slope))
		return(0)

  # Get the inoculated, transformed OD from the well (use 0 if it could not be calculated).
	inoc = inoc.OD(fitted.well, unlog = F, constant.added=0)
	if (!is.numeric(inoc))
		inoc = 0

  # Evaluate OD at the point of maximum slope
  x = fitted.well@fit.par$t50
	y = well.eval(fitted.well, x, unlog = F)

  # Calculate lag time by extending the tangent line from the point of maximum slope (in the log-transformed data)
  #  and calculating the time at the intercept with OD at inoculation, then subtracting the time at inoculation.

	output = x - (y - inoc) / max.slope

  # If the calculation returned an erroneous or negative value, return zero.
	if(is.nan(output) | is.na(output))
		return(0)
	if (output < 0)
		return(0)
	else
    return(round(output, digits))
	}




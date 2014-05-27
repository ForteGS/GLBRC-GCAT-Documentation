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
#  <well> class definition and functions. Objects contain raw          #
#  data from screening runs on single wells from 96-well plates, and   #
#  other slots for processing and model-fitting details.               #
#                                                                      #
########################################################################

# Windows OS compatibility
Sys.setlocale(locale="C")
#require(RExcelXML)

setClass("well", representation(position = "character",
					  well.info = "list",
					  screen.data = "data.frame", 
					  use.log = "logical",
					  norm = "numeric", 
					  curve.par = "list", 
					  fit.par = "list",
					  equation = "expression",
					  model.name = "character",
					  fit.info = "character",
					  add.info = "character"))

# Slots:
#   position - 3-member vector containing identifying information for the well: row (letters), column (numbers) and plate ID. 
#   well.info - a list containing strain and media names if provided
#   screen.data - a data frame with Time and raw OD values. This is the only slot that is filled upon creation of a well object. 
#                 as different functions are run on the well the data frame gets filled with additional columns. 
#   use.log - a single logical value denoting whether to return log-transformed values when data is requested from the well
#   norm - a value to subtract from all OD values before returning data. filled by <normalize.ODs> (see normalize.and.transform.R)
#   curve.par - a list of parameters that denote whether the well is empty, whether it contains ODs indicating a viable culture, whether it tanks at a certain timepoint.  

#   if model fitting using <fit.model> is successful:
#     fit.par - will be a list containing the fitted model parameters
#     equation - will contain an expression for evaluating the successfully fitted model 
#     model.name - will contain the name of the successfully fit model

#   fit.info - a message with info about whether the fit was successful, failed, or skipped. 
#   add.info - a message with info about whether jumps in OD were detected or removed, or if ODs were detected below the blank OD.  


# --------------------------------------------------------------------
# Function to create a new well (requires only Time and OD vectors, which will fill slot "screen.data")
well = function(Time = NULL, OD = NULL){
	new("well", screen.data = data.frame(Time, OD, stringsAsFactors=F))
	}

########################################################################
# Some miscellaneous functions to extract info from well objects       #
# Most of these return a single value from the well.                   #
########################################################################
#
#   Since many of these need to be applied to all wells over an array, while conserving the dimensions of 
#   that array, this file includes a wrapper function <aapply> (see bottom of file).

plate.name = function(well)
	well@position[1]

# Return the full alphanumeric well name (with leading zeros if applicable)
well.name = function(well){
	row = well@position[2]
	col = as.numeric(well@position[3])
	if (col>9)
		col = as.character(col)
	else
		col = paste("0", col, sep = "")

	paste(row,col,sep = "")
	}

is.empty = function(well)
	well@curve.par$empty.well

lacks.growth = function(well)
	well@curve.par$no.growth

tanking.start = function(well)
	well@curve.par$tanking.start

removed.points = function(well)
	(1:length(well))[well@screen.data$Remove]

remaining.points = function(well,...){
	as.numeric(rownames(data.from(well,...)))
	}

strain.name = function(well){
  if(is.null(well@well.info$Strain))
    return("<NA>")
  else
    return(well@well.info$Strain)
  }
media.name = function(well){
  if(is.null(well@well.info$Media))
    return("<NA>")
  else
    return(well@well.info$Media)
  }

raw.data = function(well)
	data.from(well, remove.tanking = F, remove = F, na.rm = F, raw.data = T)

contains.fit = function(well)
	length(well@fit.par) > 0

setMethod("length", signature(x = "well"), function(x) length(x@screen.data[,1]))

#   The <data.from> function has some options: by default it returns a two-column data frame with time and OD 
#   (or log OD if the <use.log> slot is true in the object), after normalization to the value specified in <norm> slot.    
#   - With <remove> set to true the rows specified in the <remove> column of the <screen.data> slot are not returned. 
#   - With <remove.tanking> set to true all the rows after the <tanking.start> index are removed. 
#   - Setting <raw.data> to true overrides all these settings and just returns 2 columns with Time and Raw OD.

data.from = function(well, remove = T, remove.tanking = T, raw.data = F, na.rm = F){
	
	if (length(well@use.log) == 0)
		OD.column = "OD"
	else if (well@use.log)
		OD.column = "log.OD"
	else
		OD.column = "OD"
	
	if (raw.data){
		OD.column = "OD"
		norm = 0
		}
	else if (!well@use.log)
		norm = well@norm
	else
		norm = 0

	if(remove.tanking & is.numeric(tanking.start(well)))
		well = remove.points(well, (tanking.start(well)):length(well))
	if (!remove | is.null(well@screen.data$Remove))
		output = well@screen.data[c("Time", OD.column)]
	else
		output = well@screen.data[!well@screen.data$Remove ,c("Time", OD.column)]

	output[,2] = output[,2] - norm

	if (!raw.data){
		if (!length(well@use.log))
			names(output)[2] = "Corrected.OD"
		if (!well@use.log)
			names(output)[2] = "Corrected.OD"
		}

	if (na.rm)
		output[!is.na(output[,2]),]	 
	else
		output
	}


# Functions much like <data.from> but gives a single vector containing the 
# slope at each point. Has a parameter allowing removal of NA values. 

slopes = function(well, remove = T, remove.tanking = T, na.rm = F){

	if(remove.tanking & is.numeric(tanking.start(well)))
		well = remove.points(well, (tanking.start(well)):length(well))
	if (!remove | is.null(well@screen.data$Remove))
		output = well@screen.data$Slope
	else
		output = well@screen.data$Slope[!well@screen.data$Remove]

	if (na.rm)
		output[!is.na(output)]	 
	else
		output
	}

# -----------------------------------------------------------------------
# Well array functions: these must be used on entire arrays of well objects
# instead of single ones. 

plate.names = function(well.array)
	dimnames(well.array)[[3]]

tanking.start.values = function(well.array, array = F){
	if (array)
		aapply(well.array, function(well) tanking.start(well))
	else
		sapply(well.array, function(well) tanking.start(well))
	}

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

# Notes by Jason
# 3/15/2011

# Note: package RExcelXML is generally not supported on Windows.   
#require(RExcelXML)

########################################################################
#                                                                      #
#    Function for loading experiment data from .xlsx format            #
#                                                                      #
########################################################################
#
#  This function handles loading data from excel spreadsheet format (.xlsx format) 
#  and returns an array of well objects, each filled with raw Time vs. OD data.  
#  It can handle multiple worksheets, and calls <gcat.load.data> separately on each one. 
#  It also extracts plate layout information if a correctly labeled worksheet exists 


# Arguments : see <HTS.fit.R>

load.xlsx = function(file.name = NULL, input.data = NULL, blank.value = NULL, single.plate.ID = "", single.plate = F, silent = T){
   
  ########################################################################
  #   Using RExcelXML, read from .xlsx into a list of data frame objects #
  ########################################################################
  # 
	
  if(is.null(input.data)){
    if (!silent)
  		cat("reading the excel spreadsheet...\n")
  		
    # Use <fix.xlsx> (below) to remove formatting issues with the <read.xlsx> function
    xlsx.in = fix.xlsx(read.xlsx(file.name, as.list = T))
    } else
    xlsx.in = input.data
    
  # Read plate layout and media definitions information, if it is included	
  layout.sheet = which(tolower(names(xlsx.in)) == tolower(xlsx.layout.sheet))

  if(length(layout.sheet) > 0) plate.layout = xlsx.in[[layout.sheet[1]]]
  else plate.layout = NULL
  
  if (!silent){
    if(is.null(plate.layout))
      cat("\tPlate layout info not found, assuming all wells are full\n")
    else
      cat("\tPlate layout info read successfully\n")	
  } 
  
	# Determine the first part of the plate ID for single plate format data, if not specified. 
  if (is.null(single.plate.ID)){
    # Split the file name by "." and discard the last member (file extension).
    single.plate.ID = strsplit(basename(file.name),"\\.")[[1]]
    single.plate.ID = paste(single.plate.ID[-length(single.plate.ID)],collapse=".")
    }
	 
  ########################################################################
  #          Load data from each sheet into an array                     #
  ########################################################################
  #	
  # Read all worksheets from the excel file containing actual data i.e. with headers matching those specified in <xlsx.data.headers>	
  xlsx.data.headers = xlsx.headers #constant assigned in rails
  well.array = NULL
  
	for(i in 1:length(xlsx.in)){ 
	  if (names(xlsx.in[[i]])[1] %in% xlsx.data.headers){
     	if (!silent)
	    	cat("Reformatting sheet", i, "\n")
	    
      # Call <gcat.load.data> on data from each sheet.	
      
      # Use <plate.layout> from above, if it exists.
      
      data.sheet = gcat.load.data(input.data = xlsx.in[[i]], blank.value = blank.value, single.plate = single.plate, 
        single.plate.ID = paste(single.plate.ID, names(xlsx.in)[i], sep=""), plate.layout = plate.layout)
   
      # Return an error if there is a problem with file loading.  	
      if (class(data.sheet) == "try-error")
      	stop("Error in <gcat.load.data>")
      	
      # Start a <well.array> object or add the new sheet to existing data 
  		if(is.null(well.array))
  			well.array = data.sheet
  		else
  			well.array = gcat.append.arrays(well.array, data.sheet)	
		  }
    }
  if (is.null(well.array))
    stop("No data found in workbook(check cell A1 for the correct heading?)")
  return(well.array)
  }


########################################################################
#                                                                      #
#     Function for fixing format issues in Excel spreadsheets          #
#                                                                      #
########################################################################
#
#  This function is no longer being worked on - instead we should create a template file for user input 
#    to avoid unforseen formatting problems with the excel spreadsheet format.
#
#  xlsx.in should be a list of data frames (the output from <read.xlsx>)

fix.xlsx = function(xlsx.in){
    
# Cycle through input from each sheet 
  for(i in 1:length(xlsx.in)){
    data = xlsx.in[[i]]
    # Check whether sheet contains any data
    if (length(data) > 0){
      # Cycle through each column and reformat any factors as character vectors
      for(j in 1:ncol(data))
        data[,j] = as.character(unlist(data[,j]))
      # Use the first row for column headers  
      names(data) = data[1,]
      data = data[-1,]
      # Cycle through each column again, replacing NAs with blanks and reformatting numeric columns as such  
      for(j in 1:ncol(data)){
        data[is.na(data[,j]),j] = ""
        if (all(!is.na(as.numeric(data[,j]))))
          data[,j] = as.numeric(data[,j])
        }
      # Delete any empty rows      
      delete.list = c()     
      for(j in 1:nrow(data)){
        if (all (data[j,] == ""))
          delete.list = c(delete.list, j)
        }
      if (length(delete.list) > 0)   
       data = data[-delete.list,]
      # Replace spaces in column names with dots
      names(data) = sub(pattern=" ",replacement=".",x=names(data))
      # Replace the spreadsheet with reformatted one, repeat     
      xlsx.in[[i]] = data
      }
    }
    return(xlsx.in)
  }
                                 

\name{gcat.fit.main}
\alias{gcat.fit.main}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Main analysis function for GCAT}
\description{
This is the main function that handles all the analyses for data files in both single and multiple plate formats. 
It is called by the top level function \code{gcat.analysis.main} along with \code{gcat.output.main}.
}
\usage{
gcat.fit.main(file.name, input.data = NULL, load.type = \"xlsx\", layout.file = NULL, out.dir = getwd(), 
graphic.dir = paste(out.dir, \"/pics\", sep = \"\"), single.plate = F, blank.value = NULL, start.index = 2, 
time.format = \"\%Y-\%m-\%d \%H:\%M:\%S\", normalize.method = \"default\", 
add.constant = 1, use.log = T, points.to.remove = 0, fall.cutoff = -0.0025, growth.cutoff = 0.05, remove.jumps = F, 
growth.model = richards, backup.growth.model = logistic, silent = F, verbose = F)
}
\arguments{

\item{file.name}{
Complete path and file name of an excel spreadsheet (.xlsx) or comma-separated values 
(.csv) file containing growth curve data in the multiple-plate (long) format.  
See GCAT documentation for specific description of formats. 
}
  \item{input.data}{
A list of tables representing input files read with \code{read.table}.
Used to save time in cases of running multiple analyses on the same dataset.
If used, the function will ignore \code{file.name} entirely. 
}
  \item{load.type}{
A character specifying the type of input. Can be:
  %\item{\"csv\"}{
%Comma-separated values (.csv) file. GCAT uses \code{\link{read.csv}}.}
  %\item{\"xlsx\"}{
%Spreadsheet from Microsoft Excel 2007 or later. Must be in the .xlsx format.}
  %\item{\"...\"}{
  %Any other argument will default to tab-delimited files. GCAT uses \code{\link{read.table}}.}
}
  \item{layout.file}{
If \code{load.type} is \"csv\", this argument specifies the location of a layout file containing
identifying information. See GCAT documentation for specific description of formats.  
}
 \item{out.dir}{
Name a directory to which GCAT will export text results.  
Defaults to working directory. 
}
  \item{graphic.dir}{
Name a directory to which GCAT will export graphics results.  
Defaults to a \"pics\" subdirectory of \code{out.dir}. 
Note: directories will be created if not existant. 
}
  \item{single.plate}{
Is the input data in the single-plate (wide) format? If yes, GCAT calls 
\code{reorganize.single.plate.data} to convert it to the multiple-plate (long) format. 
}
  \item{blank.value}{
Blank OD measurement for uninoculated wells. By default (NULL), 
the value of the first OD measurement in each well is used.   
}
\item{start.index}{
Which timepoint should be used as the first one after inoculation?  
}
  \item{time.format}{
Either a character describing the format used to convert timestamps 
in the input to numbers representing number of seconds (see \code{\link{strptime}}), or a 
factor to divide entries in the \code{Time} column by to get the number of hours. 
}  
  \item{normalize.method}{
Describes the method used by \code{normalize.ODs} to normalize cell
density values using blank reads. Possible values:
  \code{\"first\"} Subtracts \code{blank.value} from all ODs
  \code{\"...\"} any other argument does nothing, and assumes no blank.
  \code{\"average.first\"} forces all filled wells on each plate to match the 
  average value at \code{start.index} (after subtracting \code{blank.value})
}
  \item{add.constant}{
A value for n in the log(OD + n) transformation.  
}
  \item{points.to.remove}{
A vector of integers specifying which timepoints should be removed 
across all wells. By default (0), none are marked for removal.  
}
  \item{fall.cutoff}{
A cutoff used by \code{check.slopes} to decide on thresholds for jumps and tanking. 
}
  \item{growth.cutoff}{
A threshold used by \code{check.growth} to decide whether a well displays growth. 
}
  \item{remove.jumps}{
Should jumps in OD detected by the subfunction \code{check.slopes}? 
}
  \item{growth.model}{
The primary \code{model} used to fit growth trends.
}
  \item{backup.growth.model}{
The backup \code{model} used in case of failure of the algorithm to fit the primary model. 
}
  \item{silent}{
Should all messages returned to the R console 
(by this function and its sub-functions) be suppressed?  
}
  \item{verbose}{
Should sub-functions return lengthy messages to the console? 
}
}

\value{
An array of multiple \code{well} objects containing 
the original data sorted into wells and model fits.
}
\author{
J. Shao
}
\seealso{
\code{single.plate.analysis},  \code{multiple.plate.analysis}, 
}
\examples{

}

#----------------------------------------------
# Routine to get data from a binary file
#----------------------------------------------
dffbinary<- function(FILENAME, NROW, NAMES)
{
  #Estimate the size of the file
  newdata = file(FILENAME, "rb")
  size = 1;
  while (length(a <- readBin(newdata, double(), n = 1000)) > 0) {
    size = size +1;
  }
  size = 1000*(size+1)/NROW;
  close(newdata)
  
  #Read the data in binary form
  newdata = file(FILENAME, "rb")
  datavals = readBin(newdata, double(), n = NROW * size)
  close(newdata)
  #Columns names
  dimnames <- list(
    row = seq(1,length(datavals)/NROW),
    name = NAMES
  )
  #Create matrix from array
  matavals = matrix(datavals, ncol = NROW, byrow = TRUE, dimnames = dimnames)
  #Create data.frame from matrix
  imap = data.frame(matavals)
  #Return data frame
  return(imap)
}
args <- commandArgs(TRUE)
path <- as.character(args[1])
outname <- as.character(args[2])

##Read files names
files <- list.files(path=path, pattern="*.txt")
print(files)
# using perl to manpulate file names by trimming file extension
labs <- paste("", gsub("\\.txt", "", files, perl=TRUE), sep="")
cov <- list()
for (i in labs) {
  filepath <- file.path(path,paste(i,".txt",sep=""))
  cov[[i]] <- read.table(filepath,sep = "\t", header=F, stringsAsFactors=FALSE)
  colnames(cov[[i]]) <- c("Geneid", i)
}
df <-Reduce(function(x,y) merge(x = x, y = y, by ="Geneid"), cov)
write.table(df,outname, sep="\t", quote= F, row.names = F)


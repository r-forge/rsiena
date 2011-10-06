if(!nzchar(Sys.getenv("RSIENA_TESTING"))) q("no")

runone <- function(f)
{
    message("  Running ", sQuote(basename(f)))
    outfile <- paste(basename(f), "out", sep = "")
    failfile <- paste(outfile, "fail", sep=".")
    unlink(c(outfile, failfile))
    cmd <- paste(shQuote(file.path(R.home(), "bin", "R")),
                 "CMD BATCH --no-save",
                 shQuote(f), shQuote(outfile))
    res <- system(cmd)
    if (res) {
        cat(tail(readLines(outfile), 20), sep="\n")
        file.rename(outfile, failfile)
        return(1L)
    }
    0L
}


##library(RSienaTest)
## get the data files
datafiles <- system.file("examples", package="RSiena")
files1 <- list.files(datafiles, pattern="\\.dat$", full.names=TRUE)
file.copy(files1, ".", overwrite=TRUE)

## write some initialisation data
unlink("scriptfile.R")
writeLines(c("options(error=NULL)", "set.seed(1)"), "scriptfile.R")

## now concatenate the scripts
dd <- system.file("scripts", package="RSiena")
files <- list.files(dd, pattern="\\.R$", full.names=TRUE)
for (f in files)
{
	file.append("scriptfile.R", f)
}

## now run it
res <- 0L
runone("scriptfile.R")

## now look at the differences

system("diff scriptfile.Rout scriptFile.Rout.save")

proc.time()

if (res > 0)
{
	stop("scripts failed")
}

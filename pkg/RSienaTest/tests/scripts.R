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
writeLines(c("options(error=NULL)", "set.seed(1)",
			 "options(help_type='text')"), "scriptfile.R")
if (.Platform$OS.type == "windows")
{
	cat("options(pager='console')\n", file="scriptfile.R", append=TRUE)
}
## now concatenate the scripts
dd <- system.file("scripts", package="RSiena")
files <- list.files(dd, pattern="\\.R$", full.names=TRUE)
for (f in files)
{
	if (!grepl("SNA", f))
	{
		file.append("scriptfile.R", f)
	}
}

## now run it
res <- 0L
runone("scriptfile.R")

## now look at the differences in extra braces because of executing in batch
{if (.Platform$OS.type == "windows")
{
	previousFile <- "scriptFile.Rout.win"
}
else
{
	previousFile <- "scriptFile.Rout.save"
}}
#system("diff scriptfile.Rout scriptFile.Rout.save")
library(tools)
Rdiff("scriptfile.Rout", previousFile)

proc.time()

if (res > 0)
{
	stop("scripts failed")
}

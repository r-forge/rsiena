#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snidjers/siena
# *
# * File: effectsDocumentation.r
# *
# * Description: This module contains a function for documenting the shortNames
# * and other fields of an effects object.
# *****************************************************************************/

##@effectsDocumentation Documentation
effectsDocumentation <- function(type="html", display=type=="html",
                                 filename="effects")
{
    x <- allEffects[, c("effectGroup", "effectName", "shortName",
                        "endowment", "interaction1", "interaction2",
                        "parm", "interactionType")]
    storage.mode(x$parm) <- "integer"
    names(x)[4] <- "endow?"
    names(x)[5] <- "inter1"
    names(x)[6] <- "inter2"
    names(x)[8] <- "ego?"
    x$row <- as.integer(row.names(x))
    x <- x[, c(9, 1:8)]

    myorder <- c("nonSymmetricRate",
                 "covarNonSymmetricRate",

                 "symmetricRate",
                 "covarSymmetricRate",

                 "bipartiteRate",
                 "covarBipartiteRate",

                 "behaviorRate",
                 "behaviorOneModeRate",
                 "behaviorBipartiteRate",
                 "covarBehaviorRate",

                 "nonSymmetricObjective",
                 "dyadObjective",
                 "covarNonSymmetricObjective",
                 "unspecifiedNetInteraction",
                 "nonSymmetricNonSymmetricObjective",
                 "nonSymmetricSymmetricObjective",
                 "nonSymmetricBipartiteObjective",
                 "covarNetNetObjective",

                 "symmetricObjective",
                 "dyadObjective",
                 "covarSymmetricObjective",
                 "unspecifiedNetInteraction",

                 "bipartiteObjective",
                 "dyadObjective",
                 "covarBipartiteObjective",
                 "unspecifiedNetInteraction",
                 "bipartiteNonSymmetricObjective",
                 "bipartiteSymmetricObjective",
                 "bipartiteBipartiteObjective",
                 "covarNetNetObjective",

                 "behaviorObjective",
                 "behaviorOneModeObjective",
                 "behaviorBipartiteObjective",
                 "covarBehaviorObjective",
                 "behaviorOneModeObjective2",
                 "behaviorBipartiteObjective2",
                 "unspecifiedBehaviorInteraction")

    mytab <- table(allEffects[,1])

    addtorowPos <- cumsum(c(0, mytab[myorder]))[1:37]
    addtorowText <- names(mytab[myorder])
    if (type=="latex")
    {
        addtorowText <- paste(" \\hline \\multicolumn{4}{l}{",
                              addtorowText, "}\\\\ \\hline")
        addtorowText[1] <- paste(" \\endhead ", addtorowText[1], collapse="")
    }
    else
    {
        x[is.na(x)] <- "FALSE" ## endow? field
       x[x==""] <- "<br>"
         addtorowText <- paste(' <TR> <TD colspan="8" >',
                              addtorowText, "</TD> </TR>")
    }
    add.to.row  <-  NULL
    add.to.row$pos <- lapply(addtorowPos, function(x)x)
    add.to.row$command <- as.vector(sapply(addtorowText, function(x)x))

    order2 <- match(myorder, x[, 2])
    order3 <- as.vector(mytab[myorder])

    order4 <- unlist(apply(cbind(order2, order3), 1,
                           function(x)x[1]:(x[1] + x[2] -1)))
    y <- x[order4, -2]
    row.names(y) <- 1:nrow(y)

    if (type =="latex")
    {
        filename2 <- paste(filename,".tex", sep="", collapse="")
        includefile <- paste(filename,".include.tex", sep="", collapse="")
        includepart <- paste(filename,".include", sep="", collapse="")
        print(xtable::xtable(y), add.to.row=add.to.row, file=includefile,
              tabular.environment="longtable", hline.after=c(-1),
              floating=FALSE, include.rownames=FALSE)

        cat(file=filename2, "\\documentclass[12pt,a4paper]{article}\n",
            "\\usepackage[pdftex,dvipsnames]{color}\n",
            "\\usepackage[landscape]{geometry}\n",
            "\\usepackage{longtable}\n",
            "\\pagestyle{empty}\n",
            "\\textheight=7.25in\n",
            "\\topmargin=-1in\n",
            "\\evensidemargin=-0in\n",
            "\\oddsidemargin=-0in\n",
            "\\marginparwidth=-0in\n",
            "\\begin{document}\n",
            "\\include{", includepart, "}\n",
            "\\end{document}\n", sep=""
            )
    }
    else
    {
        filename <- paste(filename,".html", sep="", collapse="")
        print(xtable::xtable(y), add.to.row=add.to.row,
              file=filename,
              type="html", hline.after=c(-1),
              sanitize.text.function=function(x){x},
              floating=FALSE, include.rownames=FALSE,
              html.table.attributes="border=1 cellpadding=3 cellspacing=0")
        if (display)
        {
            browseURL(paste("file://", getwd(), "/", filename,
                            collapse="", sep=""))
        }
    }
}

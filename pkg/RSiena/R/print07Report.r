#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snidjers/siena
# *
# * File: print07report.r
# *
# * Description: This module contains the code to produce the report on a
# * siena07 model fit.
# *****************************************************************************/
##@PrintReport siena07 Print report
PrintReport <- function(z, x)
{
    Report('\n\n', outf)
    Heading(2, outf, "Estimation Results.")
    if (!z$OK)
   {
       Report("Error end of estimation algorithm", outf)
   }
   else
       {
           Report("Regular end of estimation algorithm.\n", outf)
           Report(c("Total of", z$n, "iteration steps.\n\n"), outf)
           Report(c("Total of", z$n, "iteration steps.\n\n"), bof)
           Heading(3, outf, "Estimates and standard errors")
           Heading(3, bof, "Estimates and standard errors")
           if (z$cconditional) ## deal with rate parameter
           {
               Report('Rate parameters:\n', outf)
               Report('Rate parameters:\n', bof)
               if (length(z$rate) == 1)
               {
                   if (length(attr(z$f,'netnames')) == 1)
                   {
                       Report(format(' 0. Rate parameter', width = 43), outf)
                       Report(format(' 0. Rate parameter', width = 43), bof)
                   }
                   else
                   {
                       Report(format(' 0. Rate parameter of conditioning variable',
                                     width = 43), outf)
                       Report(format(' 0. Rate parameter of conditioning variable',
                                     width = 43), bof)
                   }
                   Report(format(round(z$rate[1], digits = 4), width = 9), outf)
                   Report(format(round(z$rate[1], digits = 4), width = 9), bof)
                   Report(c('  (', format(round(z$vrate[1], digits = 4),
                                        width = 9), ')\n'), sep = '', outf)
                   Report(c('  (', format(round(z$vrate[1], digits = 4),
                                        width = 9), ')\n'), sep = '', bof)
               }
               else ## observations > 2
               {
                   nn <- length(z$rate)
                   nnstr <- format(as.character(1:nn),width=2,justify='left')
                   if (z$f$nDepvars == 1 || is.null(z$f$nDepvars))
                   {
                       tmp <- paste(' 0.', nnstr, ' Rate parameter period ',
                                    1:nn, '              ',
                                    format(round(z$rate,4),width=9),
                                    '  (',format(round(z$vrate,4),width=9),
                                    ')\n', sep = '')
                   }                   else{
                       tmp <- paste(' 0.', nnstr,
                                    'Rate parameter cond. variable period ',
                                    1:nn, '              ',
                                    format(round(z$rate,4),width=9),
                                    '  (',format(round(z$vrate,4),width=9),
                                    ')\n',   sep='')
                   }
                   Report(tmp, outf, sep='')
                   Report(tmp, bof, sep='')
                   Report('\nOther parameters:\n', outf)
                   Report('\nOther parameters:\n', bof)
               }
           }
           nBehavs <- length(z$f[[1]]$behavs)
           nOneModes <- length(z$f[[1]]$nets)
           if (nBehavs > 0 && nOneModes > 0)
           {
               Report("Network Dynamics\n", outf)
           }
           ses <- ifelse(diag(z$covtheta) >= 0.0 | z$fixed,
                         paste('  (', sprintf("%9.4f",sqrt(diag(z$covtheta))),
                                             ')', sep=''), '        ---')
           if (!all(z$fixed))
           {
               ses[z$fixed] <- '  (  fixed  )'
           }
           theta <- ifelse(diag(z$covtheta) >= 0.0 | z$fixed,
                           sprintf("%9.4f",z$theta),
                           '       ---')
           if (nBehavs > 0)
           {
               behEffects <- z$effects[z$effects$netType == 'behavior',]
               behNames <- unique(behEffects$name)
               if (nBehavs > 1)
               {
                   behEffects$effectName <- paste('<',
                                             (1:nBehavs)[match(behEffects$name,
                                                                    behNames)],
                                                  '> ', behEffects$effectName,
                                                  sep='')
                   z$effects$effectName[z$effects$netType=='behavior'] <-
                       behEffects$effectName
               }
           }
           typesp <- ifelse (z$effects$type== "endow", ": ", ":  ")
           tmp <- paste(sprintf("%2d", 1:length(z$effects$effectName)),
                        '. ',format(paste(z$effects$type,
                        typesp, z$effects$effectName, sep = ''), width=50),
                         theta, ses, '\n', sep='', collapse = '')
           if (nBehavs > 0 && nOneModes > 0)
           {
               nOneModeEff <- nrow(z$effects) - nrow(behEffects)
               tmpstr <- paste(nOneModeEff + 1, '. ', sep='')
               tmpsub <- regexpr(tmpstr, tmp, fixed=TRUE)
               tmp1 <- substring(tmp, 1, tmpsub - 2)
               tmp2 <- substring(tmp, tmpsub - 1, nchar(tmp))
               Report(tmp1, outf)
               Report(tmp1, bof)
               Report('\nBehavior Dynamics\n', outf)
               Report(tmp2, outf)
               Report(tmp2, bof)
           }
           else
           {
               Report(tmp, outf)
               Report(tmp, bof)
           }
           Report('\n', outf)
           Report('\n', bof)
           if (z$cconditional && length(attr(z$f, 'netnames')) > 1)
           {
               Report(c('For conditional estimation, ',
                        'the standard errors of rate parameters\n',
                        'not used for conditioning are unreliable.'), outf)
           }
           Heading(3, outf, "Covariance matrices")
           if (any(z$fixed))
           {
               Report(c('(Values of the covariance matrix of estimates\n',
                        ' are meaningless for the fixed parameters.)\n\n'),
                      outf)
           }

           Report(c("Covariance matrix of estimates",
                    "(correlations below diagonal):\n"), outf)
           covcor <- z$covtheta
           correl <- z$covtheta/sqrt(diag(z$covtheta))[row(z$covtheta)]/
               sqrt(diag(z$covtheta))[col(z$covtheta)]
           covcor[lower.tri(covcor)] <- correl[lower.tri(correl)]
           PrtOutMat(format(round(covcor, digits = 3), width = 10), outf)
           Report(c('Derivative matrix of expected statistics X by',
                    'parameters and\n'), outf)
           Report(c("covariance/correlation matrix of X can be found using\n",
                    "summary(ans) within R,",
           " or by using the 'verbose' option in Siena07.\n "), sep = "", outf)
           Report(c('Derivative matrix of expected statistics X by',
                    'parameters:\n\n '), lf)
           PrtOutMat(z$dfrac, lf)
           Report('Covariance matrix of X (correlations below the diagonal):\n',
                  lf)
           covcor <- z$msf
           correl <- z$msf/sqrt(diag(z$msf))[row(z$msf)]/
               sqrt(diag(z$msf))[col(z$msf)]
           covcor[lower.tri(covcor)] <- correl[lower.tri(correl)]
           PrtOutMat(format(round(covcor, digits = 3), width = 10), lf)
           Report('\n', outf)
           Report('\n', lf)
      }

}

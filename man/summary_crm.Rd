% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/applied_crm.R
\name{summary_crm}
\alias{summary_crm}
\title{Provide a summary of applied_crm output}
\usage{
summary_crm(x)
}
\arguments{
\item{x}{An object assigned to be the output from applied_crm.}
}
\value{
Dataframe of the summary of the output from applied_crm.
}
\description{
summary_crm is used to return a dataframe of the summary of the
    output from applied_crm.
}
\details{
This function takes an object of class "mtd" and produces a
    dataframe containing a summary of information within the
    object. Specifically it shows the dose levels, prior probabilities, number
    of evaluable patients, number of DLTs and the posterior probability
    estimates along with confidence/probability intervals if estimated in the
    underlying object.
}
\examples{
prior  <- c(0.1, 0.3, 0.5)
target <- 0.2
tox    <- c(0, 0, 1, 0, 1, 1)
level  <- c(1, 1, 1, 2, 2, 2)

crm_obj <- applied_crm(prior, target, tox, level, no_skip_esc = TRUE, no_skip_deesc = TRUE,
                       global_coherent_esc = TRUE, stop_func = NULL)

summary_crm(crm_obj)

}
\references{
O'Quigley, J. O., Pepe, M., and Fisher, L. (1990). Continual
    reassessment method: A practical design for phase I clinical trials in
    cancer. Biometrics 46:33-48.

Cheung, Y. K. (2011). Dose Finding by the Continual Reassessment Method. New
York: Chapman & Hall/CRC Press.
}
\keyword{CRM}
\keyword{dtpcrm}

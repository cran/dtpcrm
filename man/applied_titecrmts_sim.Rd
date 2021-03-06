% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_titecrmts.R
\name{applied_titecrmts_sim}
\alias{applied_titecrmts_sim}
\title{Simulate TITE-CRM trials using specified design options}
\usage{
applied_titecrmts_sim(true_tox, prior, target, max_sample_size,
    num_sims, cohort_size = 1, obswin, minfu, recrate, initdes,
    dose_func = applied_titecrm, ...)
}
\arguments{
\item{true_tox}{A vector of 'true' underlying rates of toxicity for each of
the dose levels.}

\item{prior}{A vector of prior estimates of toxicity probabilties for the
dose levels.}

\item{target}{The target DLT rate.}

\item{max_sample_size}{The maximum number of subjects to be recruited in any
simulation.}

\item{num_sims}{The total number of simulations to be run.}

\item{cohort_size}{The size of the cohorts. Default is 1.}

\item{obswin}{The observation period for total subject follow up.}

\item{minfu}{The minimum amount of follow-up required for each subject.}

\item{recrate}{The number of subjects recruited per obswin.}

\item{initdes}{A vector specifying the doses to be assisned to subjects as
per the initial design.}

\item{dose_func}{The function to be employed in executing the CRM. Default is
applied_titecrm.}

\item{...}{Any other arguements detailed in dtp::applied_titecrm.}
}
\value{
A list containg two further lists. The first of these lists contains
    the operating charateristics of the design, the second contains the
    underlying data for each of the simulation iterations.
}
\description{
applied_titecrmts_sim is used to simulate trials using the
    two-stage time-to-event continual reassessment method with specified
    design options to determine the operating characteristics.
}
\examples{
# It may take quite long for large num_sims
prior  <- c(0.1, 0.3, 0.5)
target <- 0.2
true_tox <- c(0.05, 0.2, 0.35)
first_dose <- 1
num_sims <- 5  # recommend doing 5000 simulations for the final design
obswin = 80

applied_titecrmts_sim(true_tox = true_tox, prior = prior, target = target,
                      max_sample_size = 21, num_sims = num_sims,
                      cohort_size = 3, obswin = obswin, minfu = 20,
                      recrate = 3, initdes = c(rep(1, 3), rep(2, 3), rep(3, 15)),
                      dose_func = applied_titecrm)
                      
}
\references{
O'Quigley, J. O., Pepe, M., and Fisher, L. (1990). Continual
    reassessment method: A practical design for phase I clinical trials in
    cancer. Biometrics 46:33-48.

Cheung, Y. K. (2011). Dose Finding by the Continual Reassessment Method. New
York: Chapman & Hall/CRC Press.
}
\keyword{CRM}
\keyword{Simulations}
\keyword{TITE-CRM}
\keyword{dtpcrm}

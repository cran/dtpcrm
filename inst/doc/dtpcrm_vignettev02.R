## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(dtpcrm)

## ------------------------------------------------------------------------
number.doses= 7

start.dose.level = 3

max.sample.size  = 21

target.DLT = 0.2

cohort.size = 3


## ------------------------------------------------------------------------
prior.DLT = c(0.03, 0.07, 0.12, 0.20, 0.30, 0.40, 0.52)
prior.var = 0.75

## ------------------------------------------------------------------------
no_skip_esc = TRUE
no_skip_deesc = FALSE
global_coherent_esc = TRUE

## ------------------------------------------------------------------------
y = 0.72
x = 0.1
tox_lim = target.DLT + x

## ------------------------------------------------------------------------
n.mtd= 12     

## ------------------------------------------------------------------------
stop_func <- function(x) {
  y = stop_for_excess_toxicity_empiric(x, tox_lim = tox_lim, prob_cert = y, dose = 1, nsamps = 100000)
  if(y$stop){
    x <- y
  } else {x = stop_for_consensus_reached(x, req_at_mtd = n.mtd)}
}

## ------------------------------------------------------------------------
S1=c(0.03, 0.05, 0.07, 0.20, 0.36, 0.45, 0.55)
S2=c(0.005, 0.01, 0.03, 0.05, 0.07, 0.20, 0.36)
S3=c(0.45, 0.55, 0.65, 0.7, 0.75, 0.80, 0.85)

## ---- echo=FALSE---------------------------------------------------------
scenarios<-rbind.data.frame(S1,
                            S2,
                            S3)
dose.labels= seq(-2,4,by=1)
dimnames(scenarios)[[2]]<-paste("d(", dose.labels,")", sep="")
scenarios<-cbind.data.frame("Scenario" = c("S1", "S2", "S3"), scenarios)
knitr::kable(scenarios)

## ------------------------------------------------------------------------
set.seed(1000)
No.of.simulations = 50

## ----message=FALSE, warning=FALSE, results= "hide"-----------------------
sim1_1 <- applied_crm_sim(true_tox = S1, prior = prior.DLT, target = target.DLT, max_sample_size=max.sample.size, first_dose=start.dose.level, num_sims = No.of.simulations, cohort_size = cohort.size, dose_func = applied_crm, scale = sqrt(prior.var), stop_func = stop_func)

## ------------------------------------------------------------------------
stop.fun<-function(x, matchtxt) {
   prob<-round(ifelse( !is.null(grep(matchtxt, names(x))>0) , x[grep(matchtxt, names(x))], 0),2)
   prob[is.na(prob)]<-0
   return(prob)
}

out<- cbind.data.frame(
           rbind.data.frame("S1 True Tox" = sim1_1$summary$true_tox,
                            "S1 Prob Select"= sim1_1$summary$mtd,
                            "S1 No.Subjects"= sim1_1$summary$doses_given),
      Stop.Tox  =   c("", stop.fun(sim1_1$summary$prob_stop, "Tox"), ""),
      Stop.nMTD =   c("", stop.fun(sim1_1$summary$prob_stop, "Con"), "")
)
dimnames(out)[[2]][1:length(dose.labels)]<-paste("d(", dose.labels,")", sep="")
dimnames(out)[[1]] <- paste(c(paste("S1", " True DLT rate", sep = ""), paste("S1",  " Selection Probability", sep = ""), paste("S1",  " Mean Number of Subjects", sep = "")))


## ---- echo= FALSE, results="hide"----------------------------------------
sim1_2 <- applied_crm_sim(true_tox = S2, prior = prior.DLT, target = target.DLT, max_sample_size=max.sample.size, first_dose=start.dose.level, num_sims = No.of.simulations, cohort_size = cohort.size, dose_func = applied_crm, scale = sqrt(prior.var), stop_func = stop_func)

## ---- echo= FALSE, results="hide"----------------------------------------
sim1_3 <- applied_crm_sim(true_tox = S3, prior = prior.DLT, target = target.DLT, max_sample_size=max.sample.size, first_dose=start.dose.level, num_sims = No.of.simulations, cohort_size = cohort.size, dose_func = applied_crm, scale = sqrt(prior.var), stop_func = stop_func)

## ---- echo=FALSE---------------------------------------------------------
out2<- cbind.data.frame(
           rbind.data.frame("S2 True Tox" = sim1_2$summary$true_tox,
                            "S2 Prob Select"= sim1_2$summary$mtd,
                            "S2 No.Subjects"= sim1_2$summary$doses_given),
      Stop.Tox  =   c("", stop.fun(sim1_2$summary$prob_stop, "Tox"), ""),
      Stop.nMTD =   c("", stop.fun(sim1_2$summary$prob_stop, "Con"), "")
)
dimnames(out2)[[2]][1:length(dose.labels)]<-paste("d(", dose.labels,")", sep="")
dimnames(out2)[[1]] <- paste(c(paste("S2", " True DLT rate", sep = ""), paste("S2", " Selection Probability", sep = ""), paste("S2", " Mean Number of Subjects", sep = "")))

out3<- cbind.data.frame(
           rbind.data.frame("S3 True Tox" = sim1_3$summary$true_tox,
                            "S3 Prob Select"= sim1_3$summary$mtd,
                            "S3 No.Subjects"= sim1_3$summary$doses_given),
      Stop.Tox  =   c("", stop.fun(sim1_3$summary$prob_stop, "Tox"), ""),
      Stop.nMTD =   c("", stop.fun(sim1_3$summary$prob_stop, "Con"), "")
)
dimnames(out3)[[2]][1:length(dose.labels)]<-paste0("d(", dose.labels,")")
dimnames(out3)[[1]] <- paste(c(paste("S3",  " True DLT rate", sep = ""), paste("S3",  " Selection Probability", sep = ""), paste("S3",  " Mean Number of Subjects", sep = "")))


out<- rbind.data.frame(
      "Prior DLT" = c(prior.DLT, "", ""),
      out, 
      out2,
      out3
)      

## ---- echo=FALSE---------------------------------------------------------
knitr::kable(out)

## ------------------------------------------------------------------------
start.dose.level<-3   #(eg. 1,2,3 etc)

viola_dtp <- calculate_dtps(next_dose = start.dose.level, cohort_sizes = c(cohort.size, cohort.size, cohort.size), dose_func = applied_crm, prior = prior.DLT, target = target.DLT,
stop_func = stop_func, scale = sqrt(prior.var),
no_skip_esc = no_skip_esc, no_skip_deesc = no_skip_deesc,
global_coherent_esc = global_coherent_esc)

## ------------------------------------------------------------------------
# Using dose labels                         
viola_dtp[seq(1, ncol(viola_dtp), by=2)] <- viola_dtp[seq(1, ncol(viola_dtp), by=2)] - start.dose.level

# Indicate when stopping early occurs
indSTOP<-is.na(viola_dtp[seq(1,ncol(viola_dtp), by=2)])
viola_dtp.pretty<-viola_dtp
viola_dtp.pretty[seq(1,ncol(viola_dtp.pretty), by=2)][indSTOP]<-"STOP"
viola_dtp.pretty<-cbind.data.frame(1:nrow(viola_dtp.pretty), viola_dtp.pretty)
dimnames(viola_dtp.pretty)[[2]]<-c('Pathway', 'C1 Dose', 'C1 No.DLT',  'C2 Dose', 'C2 No.DLT','C3 Dose', 'C3 No.DLT', 'C4 Dose')


## ---- echo=FALSE---------------------------------------------------------
knitr::kable(viola_dtp.pretty)

## ---- fig.width=8, fig.height=6------------------------------------------
dtpflow(viola_dtp)

## ------------------------------------------------------------------------
DLT.outcomes<-c(0,0,0,0,0,0)
dose.level = c(3,3,3,4,4,4)

## ------------------------------------------------------------------------
C2 <- applied_crm(prior=prior.DLT, target = target.DLT, tox = DLT.outcomes, level = dose.level, 
                              stop_func = stop_func, 
                              no_skip_esc = no_skip_esc, 
                              no_skip_deesc = no_skip_deesc,
                              global_coherent_esc = global_coherent_esc,
                              scale = sqrt(prior.var))
results<-cbind.data.frame(Dose=paste0("d(", dose.labels, ")"), summary_crm(C2)[2:ncol(summary_crm(C2))])
dimnames(results)[[2]][2:ncol(results)]<-c("Prior DLT", "No.subjects", "No.DLT", "Posterior DLT (90% PI)")

## ------------------------------------------------------------------------
C2$mtd 

## ------------------------------------------------------------------------
paste0("d(", C2$mtd-start.dose.level, ")"  )

## ---- echo=FALSE---------------------------------------------------------
knitr::kable(results)

## ---- fig.width=6, fig.height=6------------------------------------------
plot_crm(C2, dose_labels = paste0("d(",dose.labels,")"), cohort_sizes = c(3,3), dose_func =applied_crm, 
        no_skip_esc = no_skip_esc, 
        no_skip_deesc = no_skip_deesc,
        global_coherent_esc = global_coherent_esc, 
        scale = sqrt(prior.var),
        ylim = c(0, 0.7), lwd = 2, cex.axis = 1.5, cex.lab = 1.4, cex = 1, cohort.last = T
        )

## ----message=FALSE, warning=FALSE, results="hide"------------------------
dtp_c3 <- calculate_dtps(next_dose = 5, cohort_sizes = c(3, 3, 3),
                         prev_tox = DLT.outcomes, prev_dose = dose.level,
                         prior = prior.DLT, target = target.DLT, 
                         no_skip_esc = no_skip_esc, 
                         no_skip_deesc = no_skip_deesc,
                         global_coherent_esc = global_coherent_esc, 
                         stop_func = stop_func, 
                         scale = sqrt(prior.var))

# Using dose labels                         
dtp_c3[seq(1, ncol(dtp_c3), by=2)] <- dtp_c3[seq(1, ncol(dtp_c3), by=2)] - start.dose.level
dtp_c3_pretty<-cbind.data.frame(1:nrow(dtp_c3), dtp_c3)
dimnames(dtp_c3_pretty)[[2]]<-c('Pathway','C3 Dose', 'C3 No.DLT',  'C4 Dose', 'C4 No.DLT','C5 Dose', 'C5 No.DLT', 'C6 Dose')

## ----echo=FALSE, warning=FALSE-------------------------------------------
knitr::kable(dtp_c3_pretty)

## ---- fig.width=8, fig.height=6------------------------------------------
dtpflow(dtp_c3, cohort.labels = c('C3', 'C4', 'C5') )



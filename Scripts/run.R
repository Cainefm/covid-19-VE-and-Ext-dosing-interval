## Script name:  run.R
## Purpose of script: Effectiveness of extending the interdoses interval for BNT162b2
## Author: FAN Min
## Date Created: 2022-11
## Email: minfan@connect.hku.hk
##
## Notes:
## Adapt the matching process into density sampling.
# cases define: 2022-1-1 ormicron? draw plot
# To do list
# sensi: add rat
# 1.  ppl without extending may over
# 2. cases after 2022-1-1 over
# 3. square testing for non-linear 
# 4. interaction of sex

options(scipen = 99999)
library(readxl)
library(survival)
library(tidyverse)
library(tableone)
library(data.table)

# functions ---------------------------------------------------------------
covar_dz <- function(ct, name, icd,dx_casectrl) {
  hx <- left_join(dx_casectrl, cc[, .(patient_pssn, index_date)],by="patient_pssn")[date_past_event < index_date  &
                                                                                      grepl(icd, codes, ignore.case = T), 
                                                                                    .(patient_pssn,index_date)][,unique(.SD)][,c(name) := 1]
  ct <- merge(ct,hx,by=c("patient_pssn","index_date"),all.x=T,allow.cartesian=TRUE)
  ct[is.na(get(name)),c(name):=0]
  cat(hx[get(name)==1,.N], "\n")
  return(ct)
}
bl_dz <- function(ct, dz_def,dx_casectrl) {
  for (i in 1:nrow(dz_def)) {
    cat(dz_def[i, Name], "...")
    ct <- covar_dz(ct, dz_def[i, Name], dz_def[i, Regex],dx_casectrl)
  }
  return(ct)
}
table1 <- function(cohort, strata="case") {;
  cohort_ <- copy(cohort)
  cohort_[, (grep("^(dx|px|ip|rx|hx)",names(cohort_),value=T)) := lapply(.SD, as.logical), .SDcols=grep("^(dx|px|ip|rx|hx)",names(cohort_),value=T)]
  tb1_def <- setDT(read_excel("./documents/codes.xlsx", sheet="paed"))
  t1 <- CreateTableOne(c("age_new","sex_new",tb1_def[!is.na(Name), Name]), c(strata), data=cohort_)
  t1_print <- print(t1, smd=T, test=F)
  return(t1_print)
}
as.data.table.TableOne <- function(t) {
  tb1_def <- setDT(read_excel("./documents/codes.xlsx", sheet="paed"))[!is.na(Name), .(Name, Description)]
  tb1_def <- rbind(list(c("age_new","sex_new"),c("Age","Sex")),tb1_def)
  tb1_def <- rbind(list(Name="n", Description="n"), tb1_def)
  t <- as.data.frame(print(t, test=F, dropEqual=T, noSpaces=T))
  varlabels <- rownames(t)
  t$Name = sub("^([a-zA-Z0-9._]+).*$", "\\1", varlabels)
  t <- merge(as.data.table(t), tb1_def, by="Name", all.x=T, sort=F)
  t$Description = paste0(t$Description, sub("^([a-zA-Z0-9._]+)", "", varlabels))
  t[!is.na(Description), Name:=Description][, Description:=NULL]
  return(t)
}
as.data.table.clogit <- function(r, n=3) {
  t <- cbind(OR=exp(r$coefficients), exp(confint(r)))
  t <- as.data.table(t, keep.rownames=T)
  #t$rn <- gsub("vaccinated.brand", "Vaccinated: ", t$rn)
  t[, `OR (95% CI)` := paste0(round(OR,n), " (", round(`2.5 %`,n), " - ", round(`97.5 %`,n), ")")]
  return(t[, .(`gp`=rn, `OR (95% CI)`)])
}
obt_past_hx <- function() {
  # past disease hx ---------------------------------------------------------
  dx_casectrl <- past_history_gopsopipaepx[patient_pssn %in% cc$patient_pssn]
  codes <- setDT(readxl::read_excel("Documents/codes.xlsx",sheet = "paed"))[!is.na(Regex)]
  codes$Name <- gsub("paed","dx",codes$Name)
  cc <- bl_dz(cc, codes,dx_casectrl)
  
  # immusupp drugs ----------------------------------------------------------
  rx_immunsupp <- merge(cc[,.(patient_pssn,index_date)],rx_2022_latest_cld,by = "patient_pssn")[presc_start_date_ymd < index_date &
                                                                                                  presc_end_date_ymd  >= index_date-90 &
                                                                                                  grepl("^8.2",bnfno_p),.(patient_pssn,index_date,rx.immunsupp=1)][,unique(.SD)]
  cc <- merge(cc,rx_immunsupp,by=c("patient_pssn","index_date"),all.x = T)
  cc[,rx.immunsupp:=ifelse(is.na(rx.immunsupp),0,rx.immunsupp)]
  # cc[, score.cci := (
  #   dx.mi + dx.chf + dx.pvd + dx.cbd + dx.copd + dx.dementia + dx.paralysis +
  #     (dx.dm_com0 &
  #        !dx.dm_com1) + dx.dm_com1 * 2 + dx.crf * 2 + (dx.liver_mild &
  #                                                        !dx.liver_modsev) + dx.liver_modsev * 3 + dx.ulcers + dx.ra + dx.aids *
  #     6 + dx.cancer * 2 + dx.cancer_mets * 6
  # )]
  # cc[,.N,score.cci]
  return(cc)
}
obt_expo <- function(cutoff_expo=28) {
  cc[, interdose_interval := as.numeric(Date.of.vaccination_2nd.dose - Date.of.vaccination_1st.dose)]
  cc[, expo:=factor(ifelse(interdose_interval<cutoff_expo,"normal","ext"),c("normal","ext"))]
  cc[, onset_time:=as.numeric(index_date-Date.of.vaccination_2nd.dose)]
  
  cc[,.N,case][order(case),message("case:control --- ",paste0(N,collapse = ":"))]
  cc[,.N,expo][order(expo),message("normal:ext --- ",paste0(N,collapse = ":"))]
  cc[,.N,.(case,expo)]
}
table1 <- function(cc,exp="expo",out="case") {
  cols <- as.data.table(readxl::read_excel("Documents/codes.xlsx",sheet = "paed"))[!is.na(Name),Name]
  colsd <- as.data.table(readxl::read_excel("Documents/codes.xlsx",sheet = "paed"))[!is.na(Name),Description]
  cc[,(cols):=lapply(.SD, function(x) factor(x,levels=c("0","1"))),.SDcols=cols]
  labelled::var_label(cc)[cols] <- colsd
  labelled::var_label(cc)[exp] <- "exposure"
  labelled::var_label(cc)[c("age_new","sex_new")] <- c("age","sex")
  cc[,case:=factor(case,labels =  c("control","exposure"))]
  
  a <- CreateTableOne(vars = c(exp,"onset_time","age_new","sex_new",cols),
                      factorVars = c(exp,cols,"sex_new"),
                      strata = out,smd = T,
                      cc)
  return(print(a,varLabels = T,test=F,dropEqual = T,nonnormal=c("onset_time","age"),printToggle = F,smd=T))
}

matching <- function(cases,control_pool,sex=NA,sen_cut=NA,test=NA) {
  cases <- merge(cases,cohort,by="patient_pssn") # 6843
  if(!is.na(sex)){
    cases <- cases[sex_new==sex]
    control_pool <- control_pool[sex==sex]
    message(sex)
  }
  if(!is.na(sen_cut)){ # sensitivity analysis -
    cases <- cases[index_date-Date.of.vaccination_2nd.dose<=sen_cut]
  }
  message("cases within the cohort: ",nrow(cases)," (",cases[,uniqueN(patient_pssn)],")")
  cases <- cases[index_date>Date.of.vaccination_2nd.dose & 
                   (index_date<Date.of.vaccination_3rd.dose | is.na(Date.of.vaccination_3rd.dose))] # remove ppl have 1st covid before 2nd dose : 1408
  message("cases within the 2nd and 3rd dose: ",nrow(cases)," (",cases[,uniqueN(patient_pssn)],")")
  # matching process --------------------------------------------------------
  setorder(control_pool,patient_pssn,sex,age,Date.of.vaccination_1st.dose,Date.of.vaccination_2nd.dose,obed_date)
  setorder(cases,patient_pssn,index_date)
  ctrl <- data.table()
  j <- 0
  set.seed(475)
  for(i in 1:nrow(cases)){ #nrow(cases)
    if(i %% floor(nrow(cases)/10)==0) cat(i,"/",nrow(cases),"\n")
    temp_case <- cases[i,]
    temp_control <- control_pool[sex==temp_case$sex_new & 
                                   temp_case$index_year-dob_y==temp_case$age_new & 
                                   (is.na(obed_date)|obed_date > temp_case$index_date) &
                                   Date.of.vaccination_2nd.dose < temp_case$index_date
    ]
    if(!is.na(sen_cut)){
      temp_control <- control_pool[temp_case$index_date-Date.of.vaccination_2nd.dose<=sen_cut]
    }
    
    if(nrow(temp_control)==0){
      message(temp_case$patient_pssn)
      j <- j+1
      # }else if(nrow(temp_control)<4){
      
      #   selected_control <- temp_control$patient_pssn[sample(1:nrow(temp_control),size=nrow(temp_control))]
      #   control <- data.table(patient_pssn=selected_control,index_date=temp_case$index_date,case=0,
      #                         index_year=year(temp_case$index_date),matched_id=temp_case$patient_pssn)
    }else{
      if(nrow(temp_control)<4) message(temp_case$patient_pssn,"---",nrow(temp_control))
      selected_control <- temp_control[sample(1:nrow(temp_control),size = min(4,nrow(temp_control)),replace = F),patient_pssn]
      control <- data.table(patient_pssn=selected_control,index_date=temp_case$index_date,case=0,
                            index_year=year(temp_case$index_date),matched_id=temp_case$patient_pssn)
      if(!is.na(test) & temp_case$patient_pssn==test){
        print(temp_case$patient_pssn)
        print(nrow(temp_control))
        print(head(.Random.seed))
        # print(temp_control$patient_pssn)
        return(temp_control)
      }
      ctrl <- rbind(ctrl,control) 
    }
  }
  cc <- rbind(cases,
              merge(ctrl, cohort, by = "patient_pssn", all.x = T)) # covid_infection:5396:21584
  message("case: control:",cc[,.N,case][order(case,decreasing = T),paste0(N,collapse = ":")])
  return(cc)
}

as.logistic <- function(model) {
  temp <- cbind(exp(coef(model)),exp(confint.default(model)))[1,]
  return(temp)
}

as.logistic.two <- function(crude_est, adjust_est, analysis_type) {
  temp <- as.logistic(crude_est)
  temp1 <- as.logistic(adjust_est)
  data.table(
    `Analysis` = analysis_type,
    `crude Odds Ratio (95% CI)`=paste0(round(temp[1],3)," (",round(temp[2],3),", ",round(temp[3],3),")"),
    `adjusted Odds Ratio (95% CI)`=paste0(round(temp1[1],3)," (",round(temp1[2],3),", ",round(temp1[3],3),")"),
    `VE (95% CI)`=paste0(round(1-temp1[1],3)," (",round(1-temp1[3],3),", ",round(1-temp1[2],3),")")
  )
}



# load data ---------------------------------------------------------------

load("Data/RDS dataset/sccs_env_202210xx.Rdata")
past_history_gopsopipaepx <- readRDS("Data/RDS dataset/cleaned/past_history_gopsopipaepx.rds")

# underlying cohort -------------------------------------------------------
cohort <- Demographic_cleaned_removed  # 7492083
cohort <- cohort[!is.na(patient_pssn)]  # 4755069
cohort <- cohort[Vaccine.Brand_1st.dose=="BioNTech/Fosun" & Vaccine.Brand_2nd.dose=="BioNTech/Fosun"] #2055313
cohort <- cohort[!grepl("01",status)] #2054969
cohort <- cohort[age_new<18 ] #136367
cohort[death_date_ymd=="", death_date_ymd := NA]

# Lab processing ----------------------------------------------------------
# covid_all <- covid_lab[result=="detected",setnames(.SD,"date_ymd","covid_date")] # change names
covid_cpr_1st <- covid_lab[result == "detected", setnames(.SD, "date_ymd", "covid_date")][order(patient_pssn, covid_date)][, .SD[1], patient_pssn] # keep each one's first positive covid records
covid_omi <- covid_cpr_1st[covid_date >= as.Date("2022-01-01")] # keep only covid records after 2022-1-1, few are in 2018

# covid cases -------------------------------------------------------------
cases_covid <- covid_omi[,.(patient_pssn,index_date=covid_date,case=1)]
cases_covid$matched_id <- cases_covid$patient_pssn
cases_covid$index_year <- year(cases_covid$index_date)

dt_model <- list()

# change outcome ----------------------------------------------------------
cases <- copy(cases_covid)
control_pool <- merge(cohort[,.(patient_pssn,sex=sex_new,age=age_new,dob_y,death_date_ymd,
                                Date.of.vaccination_1st.dose,Date.of.vaccination_2nd.dose,Date.of.vaccination_3rd.dose)],
                      covid_cpr_1st[,.(patient_pssn,covid_date)],by="patient_pssn",all.x=T
)[,.(patient_pssn,sex,age,dob_y,
     Date.of.vaccination_1st.dose,
     Date.of.vaccination_2nd.dose,
     obed_date=pmin(death_date_ymd,covid_date,Date.of.vaccination_3rd.dose,na.rm=T))]

# main --------------------------------------------------------------------
cc <- matching(cases,control_pool)
cc <- obt_past_hx()
obt_expo()
dt_model$cov$primary <- cc
crude_est <- clogit(case ~ expo + strata(matched_id)+ onset_time, data = cc) 
adjust_est <- clogit(case ~ expo + strata(matched_id)+onset_time+
                       dx.asthma+
                       dx.dm+
                       dx.epilepsy+
                       rx.immunsupp, data = cc) 
out_primary <- as.logistic.two(crude_est,adjust_est,"covid-primary")


# sen1 --------------------------------------------------------------------
obt_expo(cutoff_expo = 56)
dt_model$cov$sen1_56 <- cc
crude_est <- clogit(case ~ expo + strata(matched_id)+ onset_time, data = cc) 
adjust_est <- clogit(case ~ expo + strata(matched_id)+onset_time+ dx.asthma+ dx.dm+ dx.epilepsy+ rx.immunsupp, data = cc) 
out_sen1_56 <- as.logistic.two(crude_est,adjust_est,"covid-sen1")


# sen2 --------------------------------------------------------------------
id_case_124 <- merge(cases,cohort,by="patient_pssn")[Date.of.vaccination_2nd.dose-Date.of.vaccination_1st.dose<=124,patient_pssn]
id_ctrl_124 <- control_pool[Date.of.vaccination_2nd.dose-Date.of.vaccination_1st.dose<=124,patient_pssn]
cc <- matching(cases[patient_pssn %in% id_case_124],control_pool[patient_pssn %in% id_ctrl_124])
cc <- obt_past_hx()
obt_expo()
dt_model$cov$sen2_124 <- cc
crude_est <- clogit(case ~ expo + strata(matched_id)+ onset_time, data = cc) 
adjust_est <- clogit(case ~ expo + strata(matched_id)+onset_time+ dx.asthma+ dx.dm+ dx.epilepsy+ rx.immunsupp, data = cc) 
out_sen2_124 <- as.logistic.two(crude_est,adjust_est,"covid-sen2")


# subgroup male -----------------------------------------------------------
cc <- matching(cases,control_pool,"M")
cc <- obt_past_hx()
obt_expo()
dt_model$cov$male <- cc
crude_est <- clogit(case ~ expo + strata(matched_id)+ onset_time, data = cc) 
adjust_est <- clogit(case ~ expo + strata(matched_id)+onset_time+ dx.asthma+ dx.dm+ dx.epilepsy+ rx.immunsupp, data = cc) 
out_male <- as.logistic.two(crude_est,adjust_est,"covid-male")

cc <- matching(cases,control_pool,"F")
cc <- obt_past_hx()
obt_expo()
dt_model$cov$female <- cc
crude_est <- clogit(case ~ expo + strata(matched_id)+ onset_time, data = cc) 
adjust_est <- clogit(case ~ expo + strata(matched_id)+onset_time+ dx.asthma+ dx.dm+ dx.epilepsy+ rx.immunsupp, data = cc) 
out_female <- as.logistic.two(crude_est,adjust_est,"covid-female")


# subgroup age ------------------------------------------------------------
# 
# id_case_age11 <- merge(cases,cohort,by="patient_pssn")[age_new<=11,patient_pssn]
# id_ctrl_age11 <- control_pool[age<=11,patient_pssn]
# cc <- matching(cases[patient_pssn %in% id_case_age11],control_pool[patient_pssn %in% id_ctrl_age11])
# cc <- obt_past_hx()
# obt_expo()
# dt_model$cov$age11 <- cc
# crude_est <- clogit(case ~ expo + strata(matched_id)+ onset_time, data = cc) 
# adjust_est <- clogit(case ~ expo + strata(matched_id)+onset_time+dx.seizure+dx.asthma+dx.dm, data = cc) 
# out_age11 <- as.logistic.two(crude_est,adjust_est,"covid-age11")
# 
# id_case_age12 <- merge(cases,cohort,by="patient_pssn")[age_new>11,patient_pssn]
# id_ctrl_age12 <- control_pool[age>11,patient_pssn]
# cc <- matching(cases[patient_pssn %in% id_case_age12],control_pool[patient_pssn %in% id_ctrl_age12])
# cc <- obt_past_hx()
# obt_expo()
# dt_model$cov$age12 <- cc
# crude_est <- clogit(case ~ expo + strata(matched_id)+ onset_time, data = cc) 
# adjust_est <- clogit(case ~ expo + strata(matched_id)+onset_time+dx.seizure+dx.asthma+dx.dm, data = cc) 
# out_age12 <- as.logistic.two(crude_est,adjust_est,"covid-age12")


# sen3 --------------------------------------------------------------------
cases_rat <- merge(readRDS("D:/R_batch7/LAB_CHP_COVID_0815.RDS"),
                   as.data.table(readRDS("D:/R_batch7/doc_no.RDS")),
                   by.x="PseudoID",by.y="pseudo_hkid")[,date_ymd:=pmin(ymd(date),ymd(report.date),na.rm = T)][,.(patient_pssn,date_ymd,order=99,result="detected")][]
covid_rat_1st <- cases_rat[result == "detected", setnames(.SD, "date_ymd", "covid_date")][order(patient_pssn, covid_date)][, .SD[1], patient_pssn][covid_date >= as.Date("2022-01-01")]

covid_rat_1st <- covid_rat_1st[,.(patient_pssn,index_date=covid_date,case=1)]
covid_rat_1st$matched_id <- covid_rat_1st$patient_pssn
covid_rat_1st$index_year <- year(covid_rat_1st$index_date)

cases <- copy(covid_rat_1st)
control_pool <- merge(cohort[,.(patient_pssn,sex=sex_new,age=age_new,dob_y,death_date_ymd,
                                Date.of.vaccination_1st.dose,Date.of.vaccination_2nd.dose,Date.of.vaccination_3rd.dose)],
                      covid_rat_1st[,.(patient_pssn,covid_date=index_date)],by="patient_pssn",all.x=T
)[,.(patient_pssn,sex,age,dob_y,
     Date.of.vaccination_1st.dose,Date.of.vaccination_2nd.dose,
     obed_date=pmin(death_date_ymd,covid_date,Date.of.vaccination_3rd.dose,na.rm=T))]
cc <- matching(cases,control_pool)
cc <- obt_past_hx()
obt_expo()
dt_model$cov$sen3_rat <- cc
crude_est <- clogit(case ~ expo + strata(matched_id)+ onset_time, data = cc) 
adjust_est <- clogit(case ~ expo + strata(matched_id)+onset_time+ dx.asthma+ dx.dm+ dx.epilepsy+ rx.immunsupp, data = cc) 
out_sen3_rat <- as.logistic.two(crude_est,adjust_est,"covid-sen3")




# hos ---------------------------------------------------------------------
# covid hospitalization cases ---------------------------------------------
cases_hosp <- merge(covid_omi[, .(patient_pssn, covid_date)],
                    ip[, .(patient_pssn, adm_date_ymd = date_ymd, dischg_date_ymd)],
                    by = "patient_pssn")[adm_date_ymd - covid_date >= 0 &
                                           adm_date_ymd - covid_date <= 28, unique(.SD)
                    ][order(patient_pssn,adm_date_ymd),.SD[1],patient_pssn]
cases_hosp <- merge(cases_hosp,cohort[,.(patient_pssn,sex=sex_new)],by="patient_pssn")
cases_hosp <- cases_hosp[,.(patient_pssn,index_date=adm_date_ymd,case=1,matched_id=patient_pssn,index_year=year(adm_date_ymd))]

# reviewers comments
cases_hosp <- cases_hosp[,.(patient_pssn,index_date=covid_date,case=1,matched_id=patient_pssn,index_year=year(adm_date_ymd))]


cases <- copy(cases_hosp[,sex:=NULL])
control_pool <- merge(cohort[,.(patient_pssn,sex=sex_new,age=age_new,dob_y,death_date_ymd,
                                Date.of.vaccination_1st.dose,Date.of.vaccination_2nd.dose,Date.of.vaccination_3rd.dose)],
                      cases_hosp[,.(patient_pssn,covid_date=index_date)],by="patient_pssn",all.x=T
)[,.(patient_pssn,sex,age,dob_y,
     Date.of.vaccination_1st.dose,Date.of.vaccination_2nd.dose,
     obed_date=pmin(death_date_ymd,covid_date,Date.of.vaccination_3rd.dose,na.rm=T))]

# main --------------------------------------------------------------------
cc <- matching(cases,control_pool)
cc <- obt_past_hx()
obt_expo()
dt_model$covhos$primary <- cc
crude_est <- clogit(case ~ expo + strata(matched_id)+ onset_time, data = cc) 
adjust_est <- clogit(case ~ expo + strata(matched_id)+onset_time+ dx.asthma+ dx.dm+ dx.epilepsy+ rx.immunsupp, data = cc) 
out_covhos_primary <- as.logistic.two(crude_est,adjust_est,"covidhos-primary")


# sen1 --------------------------------------------------------------------
obt_expo(cutoff_expo = 56)
dt_model$covhos$sen1_56 <- cc
crude_est <- clogit(case ~ expo + strata(matched_id)+ onset_time, data = cc) 
adjust_est <- clogit(case ~ expo + strata(matched_id)+onset_time+ dx.asthma+ dx.dm+ dx.epilepsy+ rx.immunsupp, data = cc) 
out_covhos_sen1_56 <- as.logistic.two(crude_est,adjust_est,"covid-sen1")


# sen2 --------------------------------------------------------------------
id_case_124 <- merge(cases,cohort,by="patient_pssn")[Date.of.vaccination_2nd.dose-Date.of.vaccination_1st.dose<=124,patient_pssn]
id_ctrl_124 <- control_pool[Date.of.vaccination_2nd.dose-Date.of.vaccination_1st.dose<=124,patient_pssn]
cc <- matching(cases[patient_pssn %in% id_case_124],control_pool[patient_pssn %in% id_ctrl_124])
cc <- obt_past_hx()
obt_expo()
dt_model$covhos$sen2_124 <- cc
crude_est <- clogit(case ~ expo + strata(matched_id)+ onset_time, data = cc) 
adjust_est <- clogit(case ~ expo + strata(matched_id)+onset_time+ dx.asthma+ dx.dm+ dx.epilepsy+ rx.immunsupp, data = cc) 
out_covhos_sen2_124 <- as.logistic.two(crude_est,adjust_est,"covid-sen2")

# subgroup male -----------------------------------------------------------
cc <- matching(cases,control_pool,"M")
cc <- obt_past_hx()
obt_expo()
dt_model$covhos$male <- cc
crude_est <- clogit(case ~ expo + strata(matched_id)+ onset_time, data = cc) 
adjust_est <- clogit(case ~ expo + strata(matched_id)+onset_time+ dx.asthma+ dx.dm+ dx.epilepsy+ rx.immunsupp, data = cc) 
out_covhos_male <- as.logistic.two(crude_est,adjust_est,"covid-hos-male")

cc <- matching(cases,control_pool,"F")
cc <- obt_past_hx()
obt_expo()
dt_model$covhos$female <- cc
crude_est <- clogit(case ~ expo + strata(matched_id)+ onset_time, data = cc) 
adjust_est <- clogit(case ~ expo + strata(matched_id)+onset_time+ dx.asthma+ dx.dm+ dx.epilepsy+ rx.immunsupp, data = cc) 
out_covhos_female <- as.logistic.two(crude_est,adjust_est,"covid-hos-female")



# subgroup age ------------------------------------------------------------
# 
# id_case_age11 <- merge(cases,cohort,by="patient_pssn")[age_new<=11,patient_pssn]
# id_ctrl_age11 <- control_pool[age<=11,patient_pssn]
# cc <- matching(cases[patient_pssn %in% id_case_age11],control_pool[patient_pssn %in% id_ctrl_age11])
# cc <- obt_past_hx()
# obt_expo()
# dt_model$covhos$age11 <- cc
# crude_est <- clogit(case ~ expo + strata(matched_id)+ onset_time, data = cc) 
# adjust_est <- clogit(case ~ expo + strata(matched_id)+onset_time+dx.seizure+dx.asthma+dx.dm, data = cc) 
# out_covhos_age11 <- as.logistic.two(crude_est,adjust_est,"covid-hos-age11")
# 
# id_case_age12 <- merge(cases,cohort,by="patient_pssn")[age_new>11,patient_pssn]
# id_ctrl_age12 <- control_pool[age>11,patient_pssn]
# cc <- matching(cases[patient_pssn %in% id_case_age12],control_pool[patient_pssn %in% id_ctrl_age12])
# cc <- obt_past_hx()
# obt_expo()
# dt_model$covhos$age12 <- cc
# crude_est <- clogit(case ~ expo + strata(matched_id)+ onset_time, data = cc) 
# adjust_est <- clogit(case ~ expo + strata(matched_id)+onset_time+dx.seizure+dx.asthma+dx.dm, data = cc) 
# out_covhos_age12 <- as.logistic.two(crude_est,adjust_est,"covid-hos-age12")





# sen3 --------------------------------------------------------------------
cases_rat <- merge(readRDS("D:/R_batch7/LAB_CHP_COVID_0815.RDS"),
                   as.data.table(readRDS("D:/R_batch7/doc_no.RDS")),
                   by.x="PseudoID",by.y="pseudo_hkid")[,date_ymd:=pmin(ymd(date),ymd(report.date),na.rm = T)][,.(patient_pssn,date_ymd,order=99,result="detected")][]
covid_rat_1st <- cases_rat[result == "detected", setnames(.SD, "date_ymd", "covid_date")][order(patient_pssn, covid_date)][, .SD[1], patient_pssn][covid_date >= as.Date("2022-01-01")]

cases_hosp <- merge(covid_rat_1st[, .(patient_pssn, covid_date)],
                    ip[, .(patient_pssn, adm_date_ymd = date_ymd, dischg_date_ymd)],
                    by = "patient_pssn")[adm_date_ymd - covid_date >= 0 &
                                           adm_date_ymd - covid_date <= 28, unique(.SD)
                    ][order(patient_pssn,adm_date_ymd),.SD[1],patient_pssn]
cases_hosp <- merge(cases_hosp,cohort[,.(patient_pssn,sex=sex_new)],by="patient_pssn")
cases_hosp <- cases_hosp[,.(patient_pssn,index_date=adm_date_ymd,case=1,matched_id=patient_pssn,index_year=year(adm_date_ymd))]


cases <- copy(cases_hosp)
control_pool <- merge(cohort[,.(patient_pssn,sex=sex_new,age=age_new,dob_y,death_date_ymd,
                                Date.of.vaccination_1st.dose,Date.of.vaccination_2nd.dose,Date.of.vaccination_3rd.dose)],
                      covid_rat_1st[,.(patient_pssn,covid_date)],by="patient_pssn",all.x=T
)[,.(patient_pssn,sex,age,dob_y,
     Date.of.vaccination_1st.dose,Date.of.vaccination_2nd.dose,
     obed_date=pmin(death_date_ymd,covid_date,Date.of.vaccination_3rd.dose,na.rm=T))]
cc <- matching(cases,control_pool)
cc <- obt_past_hx()
obt_expo()
dt_model$covhos$sen3_rat <- cc
crude_est <- clogit(case ~ expo + strata(matched_id)+ onset_time, data = cc) 
adjust_est <- clogit(case ~ expo + strata(matched_id)+onset_time+ dx.asthma+ dx.dm+ dx.epilepsy+ rx.immunsupp, data = cc) 
out_covhos_sen3_rat <- as.logistic.two(crude_est,adjust_est,"covid-sen3")






xlsx::write.xlsx(table1(dt_model$cov$primary),file="model estimation.xlsx",sheet="table 1",append = T)
xlsx::write.xlsx(table1(dt_model$covhos$primary),file="model estimation.xlsx",sheet="table 1_hos",append = T)
xlsx::write.xlsx(rbind(out_primary,out_sen1_56,out_sen2_124,out_sen3_rat),file = "model estimation.xlsx",sheet="infection",append = T,row.names = F)
xlsx::write.xlsx(rbind(out_covhos_primary,out_covhos_sen1_56,out_covhos_sen2_124,out_covhos_sen3_rat),file = "model estimation.xlsx",sheet="hosp",append = T,row.names = F)
xlsx::write.xlsx(rbind(out_male,out_female,out_covhos_male,out_covhos_female),file = "model estimation.xlsx",sheet="subgroup",append = T,row.names = F)
# xlsx::write.xlsx(rbind(out_age11,out_age12,out_covhos_age11,out_covhos_age12),file = "model estimation.xlsx",sheet="subgroup_age",append = T,row.names = F)

saveRDS(dt_model,"Data/RDS dataset/out/VE Ext_20221213.rds")


dt_model <- readRDS("Data/RDS dataset/out/VE Ext_20221213.rds")
cc <- dt_model$covhos$primary
# cc[,.N,expo]
# cc[, expo:=factor(ifelse(interdose_interval<28,"normal",
#                          ifelse(interdose_interval < 56,"28-55",ifelse(interdose_interval<84,"56-83",ifelse(interdose_interval<112,"84-111","ext-112 and above")))),
#                   c("normal","28-55","56-83","84-111","ext-112 and above"))]
# crude_est <- clogit(case ~ expo + strata(matched_id)+ onset_time, data = cc)
# adjust_est <- clogit(case ~ expo + strata(matched_id)+onset_time+
#                        dx.asthma+
#                        dx.dm+
#                        dx.epilepsy+
#                        rx.immunsupp, data = cc)

as.logistic_twogrp <- function(model) {
  temp <- cbind(exp(coef(model)),exp(confint.default(model)))[c(1,2,3),]
  return(temp)
}
as.logistic.two_twogrp <- function(crude_est, adjust_est, analysis_type) {
  temp <- apply(as.logistic_twogrp(crude_est),1,function(x) paste0(round(x[1],3)," (",round(x[2],3),"-",round(x[3],3),")"))
  temp1 <- apply(as.logistic_twogrp(adjust_est),1,function(x) paste0(round(x[1],3)," (",round(x[2],3),"-",round(x[3],3),")"))
  out <- as.data.table(cbind(temp,temp1),keep.rownames = T)
  setnames(out,c("Analysis","crude Odds Ratio (95% CI)","adjusted Odds Ratio (95% CI)"))
  print(out)
  # data.table(
  #   `Analysis` = analysis_type,
  #   `crude Odds Ratio (95% CI)`=paste0(round(temp[1],3)," (",round(temp[2],3),", ",round(temp[3],3),")"),
  #   `adjusted Odds Ratio (95% CI)`=paste0(round(temp1[1],3)," (",round(temp1[2],3),", ",round(temp1[3],3),")")
  # )
}
# as.logistic.two_twogrp(crude_est,adjust_est,"covid-sen3")

cc[, expo:=factor(
  ifelse(interdose_interval < 28,"normal",
         ifelse(interdose_interval < 56,"28-55",
                ifelse(interdose_interval < 84,"56-83","ext-84 and above"))),
  c("normal","28-55","56-83","ext-84 and above"))]
crude_est <- clogit(case ~ expo + strata(matched_id)+ onset_time, data = cc) 
adjust_est <- clogit(case ~ expo + strata(matched_id)+onset_time+
                       dx.asthma+
                       dx.dm+
                       dx.epilepsy+
                       rx.immunsupp, data = cc) 
as.logistic.two_twogrp(crude_est,adjust_est,"covid-sen3")



writexl::write_xlsx(list(table1=as.data.table.TableOne(table1(cc,strata = "case")),
                         crude_est=as.data.table.clogit(crude_est),
                         adjust_est=as.data.table.clogit(adjust_est)),
                    "table 1_VE and ext interval.xlsx")

writexl::write_xlsx(list(crude_est_hos=as.data.table.clogit(crude_est),
                         adjust_est_hos=as.data.table.clogit(adjust_est)),
                    "table 1_VE and ext interval_hosp.xlsx")


cc[,expo_gp:=floor(interdose_interval/7)]
as.data.table.clogit(clogit(case ~ expo_gp + onset_time + strata(matched_id), data = cc))


library(ggplot2)
ggplot(data=cohort[,.N,as.numeric(Date.of.vaccination_2nd.dose-Date.of.vaccination_1st.dose)],aes(x = as.numeric, y = N))+
  geom_bar(stat = "identity", color = "black", fill = "grey") +
  labs(title = "Frequency by interdose interval in the adolescents who finished priming doses of BNT162b2\n", x = "\ninterdose interval", y = "Frequency\n") +
  ggthemes::theme_stata()+
  scale_x_continuous(breaks=0:20*50)



cohort[,quantile(as.numeric(Date.of.vaccination_2nd.dose-Date.of.vaccination_1st.dose),0.80)] #


cc[, interdose_interval := as.numeric(Date.of.vaccination_2nd.dose - Date.of.vaccination_1st.dose)]
cc[, expo:=factor(ifelse(interdose_interval<56,"normal","ext"),c("normal","ext"))]
cc[, onset_time:=as.numeric(index_date-Date.of.vaccination_2nd.dose)]
cc[,.N,case]
cc[,.N,.(case,expo)]
crude_est <- clogit(case ~ expo + strata(matched_id)+ onset_time, data = cc) 
table1(cc,strata = "case")
adjust_est <- clogit(case ~ expo + strata(matched_id)+onset_time+dx.seizure+dx.asthma+dx.dm, data = cc) 
as.data.table.clogit(crude_est)
as.data.table.clogit(adjust_est)


cc[, interdose_interval := as.numeric(Date.of.vaccination_2nd.dose - Date.of.vaccination_1st.dose)]
cc[, expo:=factor(ifelse(interdose_interval<28,"normal","ext"),c("normal","ext"))]
cc[, onset_time:=as.numeric(index_date-Date.of.vaccination_2nd.dose)]
cc <- cc[interdose_interval<124]
cc[,.N,case]
cc[,.N,.(case,expo)]
crude_est <- clogit(case ~ expo + strata(matched_id)+ onset_time, data = cc) 
adjust_est <- clogit(case ~ expo + strata(matched_id)+onset_time+dx.seizure+dx.asthma+dx.dm, data = cc) 
as.data.table.clogit(crude_est)
as.data.table.clogit(adjust_est)


testing <- copy(cohort)
testing[, interdose_interval := as.numeric(Date.of.vaccination_2nd.dose - Date.of.vaccination_1st.dose)]
testing[, expo:=factor(ifelse(interdose_interval<56,"normal","ext"),c("normal","ext"))]
testing[, before_announcement:= Date.of.vaccination_2nd.dose < as.Date("2022-3-1")]

testing[,.N,.(expo,before_announcement)]




a<- dcast(cc[,.N,.(age_new,case)],age_new~case)
b<- dcast(dt_model$covhos$primary[,.N,.(age_new,case)],age_new~case)
setnames(a,c("0","1"),c("covid-19-control","covid_19-cases"))
setnames(b,c("0","1"),c("covid-19-hos-control","covid_19-hos-cases"))
ab<- merge(a,b)

write.csv(ab,"reply_review_age_distribution.csv")
sum(ab$`covid-19-control`)
sum(ab$`covid_19-cases`)
sum(ab$`covid-19-hos-control`)
sum(ab$`covid_19-hos-cases`)

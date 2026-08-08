
# =============================================================================
# Discordance between 10-year and 30-year cardiovascular risk estimates
#   among US adults with chronic kidney disease (NHANES 1999-2018)
#
# Full replication script: cohort construction -> PREVENT scoring ->
#   survey-weighted discordance, reclassification & 95% CIs (Taylor
#   linearization) under BOTH threshold definitions -> Figure 1 (vector PDF).
#
# Run:  Rscript prevent_ckd_replication.R
# Input:  data/ckm_3.27_v3.csv   (NHANES-derived analytic file, provided)
# Output: results/*.csv, figures/figure1.pdf (+ .png)
#
# Packages: preventr (CRAN, PREVENT equations), survey, tidyverse, arrow,
#           ggplot2, patchwork.  preventr version used: 0.12.0
# =============================================================================

suppressPackageStartupMessages({
  library(preventr); library(tidyverse); library(purrr)
  library(survey);   library(arrow);     library(patchwork)
})
options(survey.lonely.psu = "adjust")

DATA_FILE <- "data/ckm_3.27_v3.csv"
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
message("preventr version: ", as.character(packageVersion("preventr")))

# ---- PREVENT input clamping ranges (per equation documentation) -------------
PREVENT_RANGES <- list(
  sbp=list(lo=90,hi=180), tc=list(lo=130,hi=320), hdl=list(lo=20,hi=100),
  bmi=list(lo=18.5,hi=39.9), egfr=list(lo=15,hi=140), uacr=list(lo=1,hi=9800))
clamp <- function(x,lo,hi) pmax(lo,pmin(hi,x))

# 10-year clinical action thresholds
THRESH_10 <- 10; THRESH_ASCVD_10 <- 7.5; THRESH_HF_10 <- 5
# 30-year high-risk fixed threshold
THRESH_30 <- 20

# =============================================================================
# 1. COHORT CONSTRUCTION
# =============================================================================
dat_raw <- read_csv(DATA_FILE, show_col_types = FALSE)

prep <- dat_raw |> mutate(
  sex_chr=RIAGENDR, sex_lc=tolower(RIAGENDR),
  diabetes=as.integer(any_dmrx_==1 | (!is.na(LBXGH) & LBXGH>=6.5)),
  bp_tx=as.integer(any_htnrx_==1), statin=as.integer(statin_==1),
  smoking=as.integer(currsmoke==1), wt=WTMEC8YR,
  sbp_sc =clamp(sbp,    PREVENT_RANGES$sbp$lo, PREVENT_RANGES$sbp$hi),
  tc_sc  =clamp(LBXTC,  PREVENT_RANGES$tc$lo,  PREVENT_RANGES$tc$hi),
  hdl_sc =clamp(hdl,    PREVENT_RANGES$hdl$lo, PREVENT_RANGES$hdl$hi),
  bmi_sc =clamp(BMXBMI, PREVENT_RANGES$bmi$lo, PREVENT_RANGES$bmi$hi),
  egfr_sc=clamp(egfr_cr,PREVENT_RANGES$egfr$lo,PREVENT_RANGES$egfr$hi),
  uacr_sc=clamp(uacr,   PREVENT_RANGES$uacr$lo,PREVENT_RANGES$uacr$hi),
  ckd_mech=case_when(
    uacr>=30 & egfr_cr>=60 ~ "UACR-driven",
    uacr< 30 & egfr_cr< 60 ~ "eGFR-driven",
    uacr>=30 & egfr_cr< 60 ~ "Both", TRUE ~ "Other"))

# CKD definition: UACR>=30 OR eGFR<60; no prior CVD; complete PREVENT inputs
apply_ckd_filter <- function(df) df |> filter(
  !is.na(uacr) & (uacr>=30 | egfr_cr<60), cvdhist==0,
  !is.na(LBXTC), !is.na(hdl), !is.na(sbp), !is.na(BMXBMI), !is.na(egfr_cr),
  !is.na(sex_chr), !is.na(wt), wt>0)

age_breaks <- c(29,34,39,44,49,54,59)
age_labels <- c("30-34","35-39","40-44","45-49","50-54","55-59")

dat30 <- prep |> filter(RIDAGEYR>=30, RIDAGEYR<=59) |> apply_ckd_filter() |>
  mutate(age_group=cut(RIDAGEYR,age_breaks,age_labels,right=TRUE),
         ageband=cut(RIDAGEYR,c(29,39,49,59),c("30-39","40-49","50-59"),right=TRUE))

# Reference population for age/sex percentile cut-offs:
#   all NHANES adults 30-59 without prior CVD & complete inputs (regardless of CKD)
ref_pop <- prep |> filter(RIDAGEYR>=30, RIDAGEYR<=59, cvdhist==0,
  !is.na(LBXTC),!is.na(hdl),!is.na(sbp),!is.na(BMXBMI),!is.na(egfr_cr),
  !is.na(sex_chr),!is.na(wt),wt>0) |>
  mutate(age_group=cut(RIDAGEYR,age_breaks,age_labels,right=TRUE))

message(sprintf("CKD cohort N=%d (weighted %.0f); reference N=%d",
                nrow(dat30), sum(dat30$wt), nrow(ref_pop)))

# =============================================================================
# 2. ATTRITION FLOW
# =============================================================================
elig <- prep |> filter(RIDAGEYR>=30, RIDAGEYR<=59, cvdhist==0)
miss_inputs <- elig |> filter(is.na(LBXTC)|is.na(hdl)|is.na(sbp)|is.na(BMXBMI)|
                                is.na(egfr_cr)|is.na(uacr))
attrition <- tibble(
  step=c("Aged 30-59, no prior CVD (eligible)",
         "Excluded: missing >=1 PREVENT input",
         "Complete PREVENT inputs",
         "Met CKD definition (UACR>=30 or eGFR<60)"),
  n=c(nrow(elig), nrow(miss_inputs), nrow(elig)-nrow(miss_inputs), nrow(dat30)))
write_csv(attrition, "results/attrition_flow.csv")

# =============================================================================
# 3. PREVENT SCORING  (per-row call to preventr::estimate_risk)
#    Base = without UACR; +UACR = albuminuria-enhanced equation. x100 -> %.
# =============================================================================
score_cohort_df <- function(dat, time_val, include_uacr=FALSE) {
  pmap_dfr(dat, function(RIDAGEYR,sex_lc,sbp_sc,bp_tx,tc_sc,hdl_sc,
                         statin,diabetes,smoking,egfr_sc,bmi_sc,uacr_sc,...) {
    args <- list(age=RIDAGEYR,sex=sex_lc,sbp=sbp_sc,bp_tx=bp_tx,
                 total_c=tc_sc,hdl_c=hdl_sc,statin=statin,dm=diabetes,
                 smoking=smoking,egfr=egfr_sc,bmi=bmi_sc,time=time_val,quiet=TRUE)
    if (include_uacr) args$uacr <- uacr_sc
    tryCatch({
      r <- do.call(estimate_risk,args)
      tibble(total_cvd=r$total_cvd*100, ascvd=r$ascvd*100, heart_failure=r$heart_failure*100)
    }, error=function(e) tibble(total_cvd=NA_real_,ascvd=NA_real_,heart_failure=NA_real_))
  })
}

message("Scoring CKD cohort (30-yr & 10-yr, base & +UACR)...")
s30b<-score_cohort_df(dat30,30,FALSE); s30u<-score_cohort_df(dat30,30,TRUE)
s10b<-score_cohort_df(dat30,10,FALSE); s10u<-score_cohort_df(dat30,10,TRUE)
dat30 <- dat30 |> mutate(
  risk30_base=s30b$total_cvd, risk30_ascvd_base=s30b$ascvd, risk30_hf_base=s30b$heart_failure,
  risk30_uacr=s30u$total_cvd, risk30_ascvd_uacr=s30u$ascvd, risk30_hf_uacr=s30u$heart_failure,
  risk10_base=s10b$total_cvd, risk10_ascvd_base=s10b$ascvd, risk10_hf_base=s10b$heart_failure,
  risk10_uacr=s10u$total_cvd, risk10_ascvd_uacr=s10u$ascvd, risk10_hf_uacr=s10u$heart_failure)

message("Scoring reference population (30-yr base) for percentile cut-offs...")
r30b <- score_cohort_df(ref_pop,30,FALSE)
ref_pop <- ref_pop |> mutate(ref_cvd=r30b$total_cvd, ref_ascvd=r30b$ascvd, ref_hf=r30b$heart_failure)

write_parquet(dat30, "results/scored_cohort.parquet")

# =============================================================================
# 4. SURVEY DESIGN + Taylor-linearization CI helpers
# =============================================================================
svy   <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE,
                   data=dat30 |> filter(!is.na(wt), wt>0))
svyref<- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE,
                   data=ref_pop |> filter(!is.na(wt), wt>0))

# proportion + logit 95% CI
ci_prop <- function(design, var){
  f <- as.formula(paste0("~",var))
  r <- svyciprop(f, design, method="logit", level=0.95)
  ci <- as.numeric(confint(r))
  c(est=as.numeric(r)*100, lo=ci[1]*100, hi=ci[2]*100)
}
# weighted count + Wald 95% CI
ci_count <- function(design, var){
  t <- svytotal(as.formula(paste0("~",var)), design, na.rm=TRUE)
  est<-as.numeric(t); se<-as.numeric(SE(t))
  c(est=est, lo=est-1.96*se, hi=est+1.96*se)
}

# =============================================================================
# 5. AGE x SEX 75th-PERCENTILE CUT-OFFS (from reference population)
# =============================================================================
p75 <- expand_grid(.sex=c("Female","Male"), .age=age_labels) |>
  pmap_dfr(function(.sex, .age){
    idx <- svyref$variables$sex_chr==.sex & as.character(svyref$variables$age_group)==.age
    sub <- svyref[idx, ]
    q <- function(v) as.numeric(coef(svyquantile(as.formula(paste0("~",v)),
                                sub, 0.75, na.rm=TRUE, ci=FALSE)))
    tibble(sex_chr=.sex, age_group=.age,
           p75_cvd=q("ref_cvd"), p75_ascvd=q("ref_ascvd"), p75_hf=q("ref_hf"))
  })
write_csv(p75, "results/percentile_cutoffs.csv")

# attach percentile cut-offs to the CKD cohort & rebuild design
dat30 <- dat30 |> left_join(p75, by=c("sex_chr","age_group"))
svy   <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE,
                   data=dat30 |> filter(!is.na(wt), wt>0))

# endpoint metadata: (key, base10, base30, uacr30, t10, p75col)
ENDPOINTS <- list(
  CVD  =list(r10="risk10_base",       r30="risk30_base",       u30="risk30_uacr",       t10=THRESH_10,       p="p75_cvd"),
  ASCVD=list(r10="risk10_ascvd_base", r30="risk30_ascvd_base", u30="risk30_ascvd_uacr", t10=THRESH_ASCVD_10, p="p75_ascvd"),
  HF   =list(r10="risk10_hf_base",    r30="risk30_hf_base",    u30="risk30_hf_uacr",    t10=THRESH_HF_10,    p="p75_hf"))

STRATA <- list(
  Overall=quote(rep(TRUE,nrow(V))), Females=quote(V$sex_chr=="Female"), Males=quote(V$sex_chr=="Male"),
  `Aged 30-39`=quote(V$ageband=="30-39"), `Aged 40-49`=quote(V$ageband=="40-49"),
  `Aged 50-59`=quote(V$ageband=="50-59"))
# subset a survey design by a logical expression evaluated against V=design$variables
svy_sub <- function(design, expr){ V<-design$variables; design[eval(expr), ] }

# ---- 6. DISCORDANCE (low 10-y AND high 30-y) under both thresholds ----------
disc_rows <- list()
for(oc in names(ENDPOINTS)){ s<-ENDPOINTS[[oc]]
  for(thr in c("20","p75")){
    V<-svy$variables
    lo10 <- V[[s$r10]] <  s$t10
    hi30 <- if(thr=="20") V[[s$r30]]>=THRESH_30 else V[[s$r30]]>=V[[s$p]]
    svy$variables$disc <- as.integer(lo10 & hi30)
    for(st in names(STRATA)){
      sub<-svy_sub(svy, STRATA[[st]]); ci<-ci_prop(sub,"disc")
      disc_rows[[length(disc_rows)+1]]<-tibble(outcome=oc,threshold=thr,stratum=st,
        pct=ci["est"], pct_lo=ci["lo"], pct_hi=ci["hi"])
    }}}
discordance <- bind_rows(disc_rows); write_csv(discordance,"results/discordance_CI.csv")

# ---- 7. RECLASSIFICATION when UACR added (gross upward) under both thresholds
recl_rows <- list()
for(oc in names(ENDPOINTS)){ s<-ENDPOINTS[[oc]]
  for(thr in c("20","p75")){
    V<-svy$variables
    cutv <- if(thr=="20") THRESH_30 else V[[s$p]]
    svy$variables$newup <- as.integer(V[[s$r30]]<cutv & V[[s$u30]]>=cutv)
    for(st in names(STRATA)){
      sub<-svy_sub(svy, STRATA[[st]])
      cN<-ci_count(sub,"newup"); cP<-ci_prop(sub,"newup")
      recl_rows[[length(recl_rows)+1]]<-tibble(outcome=oc,threshold=thr,stratum=st,
        N=cN["est"], N_lo=cN["lo"], N_hi=cN["hi"],
        pct=cP["est"], pct_lo=cP["lo"], pct_hi=cP["hi"])
    }}}
reclass <- bind_rows(recl_rows); write_csv(reclass,"results/reclassification_CI.csv")

# ---- 8. RISK-CATEGORY COMPOSITION (both lower / discordant / higher-10y) -----
comp_rows <- list()
for(oc in names(ENDPOINTS)){ s<-ENDPOINTS[[oc]]
  for(thr in c("20","p75")){
    V<-svy$variables
    hi10<-V[[s$r10]]>=s$t10
    hi30<-if(thr=="20") V[[s$r30]]>=THRESH_30 else V[[s$r30]]>=V[[s$p]]
    svy$variables$cat3<-factor(ifelse(hi10,"Higher 10-y",
                        ifelse(hi30,"Discordant","Both lower")),
                        levels=c("Both lower","Discordant","Higher 10-y"))
    for(st in names(STRATA)){
      sub<-svy_sub(svy, STRATA[[st]]); pp<-svymean(~cat3,sub,na.rm=TRUE)
      comp_rows[[length(comp_rows)+1]]<-tibble(outcome=oc,threshold=thr,stratum=st,
        both_lower=as.numeric(pp["cat3Both lower"])*100,
        discordant=as.numeric(pp["cat3Discordant"])*100,
        higher10 =as.numeric(pp["cat3Higher 10-y"])*100)
    }}}
composition <- bind_rows(comp_rows); write_csv(composition,"results/composition_both.csv")

# ---- 9. CIRCULARITY: reclassification by CKD mechanism -----------------------
svy$variables$subgroup <- ifelse(svy$variables$ckd_mech=="eGFR-driven","eGFR-only","UACR-enriched")
circ_rows <- list()
for(oc in names(ENDPOINTS)){ s<-ENDPOINTS[[oc]]
  for(thr in c("20","p75")){
    V<-svy$variables
    cutv <- if(thr=="20") THRESH_30 else V[[s$p]]
    svy$variables$newup <- as.integer(V[[s$r30]]<cutv & V[[s$u30]]>=cutv)
    for(g in c("UACR-enriched","eGFR-only")){
      sub<-svy[svy$variables$subgroup==g, ]; ci<-ci_prop(sub,"newup")
      circ_rows[[length(circ_rows)+1]]<-tibble(subgroup=g,outcome=oc,threshold=thr,
        pct=ci["est"], pct_lo=ci["lo"], pct_hi=ci["hi"])
    }}}
circularity <- bind_rows(circ_rows); write_csv(circularity,"results/circularity_CI.csv")

# ---- 10. MEDIAN RISK by stratum (with Taylor-linearized 95% CI) --------------
# svyquantile(ci=TRUE) inverts a Taylor-linearized test on the CDF for the CI.
med_ci <- function(design, v){
  q <- tryCatch(svyquantile(as.formula(paste0("~",v)), design, 0.5,
                            na.rm=TRUE, ci=TRUE), error=function(e) NULL)
  if(is.null(q)) return(c(est=NA,lo=NA,hi=NA))
  est <- as.numeric(coef(q)); ci <- as.numeric(confint(q))
  c(est=est, lo=ci[1], hi=ci[2])
}
med_rows <- list()
for(oc in names(ENDPOINTS)){ s<-ENDPOINTS[[oc]]
  for(st in names(STRATA)){
    sub<-svy_sub(svy, STRATA[[st]])
    a<-med_ci(sub,s$r10); b<-med_ci(sub,s$r30)
    med_rows[[length(med_rows)+1]]<-tibble(outcome=oc,stratum=st,
      med10=a["est"], med10_lo=a["lo"], med10_hi=a["hi"],
      med30=b["est"], med30_lo=b["lo"], med30_hi=b["hi"])
  }}
medians <- bind_rows(med_rows); write_csv(medians,"results/median_risk.csv")

message("Analyses complete. Headline check (Overall, fixed 20%):")
print(discordance |> filter(stratum=="Overall", threshold=="20"))

# =============================================================================
# 11. FIGURE 1  ->  vector PDF (ggplot2 + patchwork)
#     Palette matches the original submission figure.
# =============================================================================
LBLUE<-"#88b8d8"; NAVY<-"#103858"; REDO<-"#e05028"; BLUE<-"#3868a0"
PURPLE<-"#785090"; GREY<-"#d9d9d9"; THRLINE<-"#981830"
strata_lev <- c("Overall","Females","Males","Aged 30-39","Aged 40-49","Aged 50-59")
ep_lev     <- c("CVD","ASCVD","HF")
ep_lab     <- c(CVD="Total CVD", ASCVD="ASCVD", HF="Heart failure")
base_thm <- theme_classic(base_size=18) +
  theme(strip.background=element_blank(), strip.text=element_text(face="bold",size=19),
        plot.title=element_text(face="bold",size=19),
        axis.title=element_text(size=16),
        axis.text=element_text(colour="black",size=15), legend.position="bottom",
        legend.text=element_text(size=15),
        legend.title=element_blank(), panel.spacing=unit(1.2,"lines"))
LBL <- 5.2   # data-label text size (geom_text)

# ---- Panel A: median 10-y vs 30-y risk --------------------------------------
medA <- bind_rows(
  medians |> transmute(outcome,stratum,horizon="10-year",
                       risk=med10, lo=med10_lo, hi=med10_hi),
  medians |> transmute(outcome,stratum,horizon="30-year",
                       risk=med30, lo=med30_lo, hi=med30_hi)) |>
  mutate(stratum=factor(stratum,rev(strata_lev)),
         outcome=factor(outcome,ep_lev,ep_lab[ep_lev]))
pA <- ggplot(medA, aes(risk, stratum, fill=horizon)) +
  geom_col(position=position_dodge(width=.7), width=.6) +
  geom_errorbarh(aes(xmin=lo, xmax=hi), position=position_dodge(width=.7),
                 height=.18, colour="grey30", linewidth=.4) +
  geom_text(aes(x=hi, label=ifelse(horizon=="30-year", sprintf("%.0f", risk), "")),
            position=position_dodge(width=.7), hjust=-0.35, size=LBL,
            fontface="bold", colour=NAVY) +
  geom_vline(xintercept=THRESH_30, linetype="dashed", colour=THRLINE) +
  expand_limits(x=max(medA$hi,na.rm=TRUE)*1.18) +
  facet_wrap(~outcome, nrow=1) +
  scale_fill_manual(values=c("10-year"=LBLUE,"30-year"=NAVY)) +
  labs(title="A  Median 10-year vs 30-year predicted risk", x="Median predicted risk, %", y=NULL) +
  base_thm

# ---- Panel B: composition, both thresholds ----------------------------------
compB <- composition |>
  pivot_longer(c(both_lower,discordant,higher10),names_to="cat",values_to="pct") |>
  mutate(stratum=factor(stratum,rev(strata_lev)),
         outcome=factor(outcome,ep_lev,ep_lab[ep_lev]),
         cat=factor(cat,c("both_lower","discordant","higher10"),
                    c("Both lower","Discordant (low 10-y, high 30-y)","Higher 10-y & higher 30-y")),
         thr=factor(threshold,c("20","p75"),c("Fixed 20%","75th pctile")))
# Centre a label in EVERY wide-enough segment. ggplot stacks in reverse factor
# order, so compute the running position in that same (descending) order.
compB_lab <- compB |> group_by(outcome,threshold,stratum) |>
  arrange(desc(cat)) |>
  mutate(cum=cumsum(pct), center=cum-pct/2) |> ungroup() |>
  filter(pct>=6) |>
  mutate(txtcol=ifelse(cat=="Both lower","grey10","white"))
pB <- ggplot(compB, aes(pct, interaction(thr,stratum,sep="  "), fill=cat)) +
  geom_col(width=.72, colour="white", linewidth=.3) +
  geom_text(data=compB_lab, inherit.aes=FALSE, fontface="bold", size=LBL-1.2,
            aes(x=center, y=interaction(thr,stratum,sep="  "),
                label=sprintf("%.0f",pct), colour=txtcol)) +
  scale_colour_identity() +
  facet_wrap(~outcome, nrow=1) +
  scale_fill_manual(values=c("Both lower"=GREY,
        "Discordant (low 10-y, high 30-y)"=REDO,"Higher 10-y & higher 30-y"=BLUE)) +
  labs(title="B  Risk-category composition (top=fixed 20%, bottom=75th pctile)",
       x="US adults, %", y=NULL) + base_thm +
  theme(axis.text.y=element_text(size=12))

# ---- Panel C: reclassification counts (millions), both thresholds -----------
reclC <- reclass |> filter(stratum %in% c("Overall","Females","Males")) |>
  mutate(stratum=recode(stratum,Females="Women",Males="Men"),
         stratum=factor(stratum,rev(c("Overall","Women","Men"))),
         outcome=factor(outcome,ep_lev,ep_lab[ep_lev]),
         thr=factor(threshold,c("20","p75"),c("Fixed 20% threshold","75th-percentile threshold")))
pC <- ggplot(reclC, aes(N/1e6, stratum, fill=thr)) +
  geom_col(position=position_dodge(width=.7), width=.6) +
  geom_errorbarh(aes(xmin=pmax(N_lo,0)/1e6, xmax=N_hi/1e6),
                 position=position_dodge(width=.7), height=.2, colour="grey30") +
  geom_text(aes(x=N_hi/1e6, label=sprintf("%.2fM", N/1e6)),
            position=position_dodge(width=.7), hjust=-0.2, size=LBL-0.6, fontface="bold") +
  expand_limits(x=0) +
  scale_x_continuous(expand=expansion(mult=c(0,0.22))) +
  facet_wrap(~outcome, nrow=1, scales="free_x") +
  scale_fill_manual(values=c("Fixed 20% threshold"=BLUE,"75th-percentile threshold"=PURPLE)) +
  labs(title="C  Population reclassified upward when UACR added",
       x="No. of US adults reclassified upward (millions)", y=NULL) + base_thm

# ---- Panel D: circularity by CKD mechanism, both thresholds -----------------
circD <- circularity |>
  mutate(subgroup=factor(subgroup,c("UACR-enriched","eGFR-only"),
                         c("UACR-enriched\n(87%)","eGFR-only\n(13%)")),
         outcome=factor(outcome,ep_lev,ep_lab[ep_lev]),
         thr=factor(threshold,c("20","p75"),c("Fixed 20%","75th pctile")))
pD <- ggplot(circD, aes(subgroup, pct, fill=thr)) +
  geom_col(position=position_dodge(width=.8), width=.7) +
  geom_errorbar(aes(ymin=pmax(pct_lo,0), ymax=pct_hi),
                position=position_dodge(width=.8), width=.2, colour="grey30") +
  geom_text(aes(y=pct_hi, label=sprintf("%.1f%%", pct)),
            position=position_dodge(width=.8), vjust=-0.6, size=LBL-0.6, fontface="bold") +
  expand_limits(y=max(circD$pct_hi)*1.15) +
  facet_wrap(~outcome, nrow=1) +
  scale_fill_manual(values=c("Fixed 20%"=BLUE,"75th pctile"=PURPLE)) +
  labs(title="D  UACR reclassification by CKD mechanism",
       x=NULL, y="Reclassified upward, %") + base_thm

fig1 <- pA / pB / pC / pD + plot_layout(heights=c(1,1.5,1,0.9)) +
  plot_annotation(
    title="Discordance between 10- and 30-year PREVENT cardiovascular risk in US adults with CKD, ages 30-59",
    subtitle="NHANES 1999-2018 | N=2326 (9.6 million weighted) | Both thresholds: fixed 30-year >=20% and age/sex 75th percentile",
    theme=theme(plot.title=element_text(face="bold",size=22),
                plot.subtitle=element_text(size=16,colour="grey30")))

ggsave("figures/figure1.pdf", fig1, width=16, height=21, device=cairo_pdf)
ggsave("figures/figure1.png", fig1, width=16, height=21, dpi=300)
message("Figure written: figures/figure1.pdf (vector) + figures/figure1.png")

# =============================================================================
# 12. TABLE 1  - survey-weighted baseline characteristics of the CKD cohort
# =============================================================================
suppressPackageStartupMessages({ library(gt); library(htmltools) })

# weighted continuous: mean (SD) using survey mean & population variance
w_meansd <- function(design, var){
  m  <- as.numeric(coef(svymean(as.formula(paste0("~",var)), design, na.rm=TRUE)))
  v  <- as.numeric(coef(svyvar (as.formula(paste0("~",var)), design, na.rm=TRUE)))
  sprintf("%.1f (%.1f)", m, sqrt(v))
}
w_medIQR <- function(design, var){
  q <- svyquantile(as.formula(paste0("~",var)), design, c(.25,.5,.75), na.rm=TRUE, ci=FALSE)
  q <- as.numeric(coef(q)); sprintf("%.0f (%.0f-%.0f)", q[2], q[1], q[3])
}
# weighted % for a level of a factor/indicator
w_pct <- function(design, expr){
  design$variables$.ind <- as.integer(eval(expr, design$variables))
  p <- as.numeric(coef(svymean(~.ind, design, na.rm=TRUE)))*100
  sprintf("%.1f", p)
}
V <- svy$variables
tab1 <- tribble(
  ~Characteristic, ~Value,
  "Participants, unweighted N",                as.character(nrow(V)),
  "Represented US adults, weighted",           format(round(sum(V$wt)), big.mark=","),
  "Age, y - mean (SD)",                        w_meansd(svy,"RIDAGEYR"),
  "Female, %",                                 w_pct(svy, quote(sex_chr=="Female")),
  "Systolic BP, mm Hg - mean (SD)",            w_meansd(svy,"sbp"),
  "Total cholesterol, mg/dL - mean (SD)",      w_meansd(svy,"LBXTC"),
  "HDL cholesterol, mg/dL - mean (SD)",        w_meansd(svy,"hdl"),
  "Body-mass index, kg/m2 - mean (SD)",        w_meansd(svy,"BMXBMI"),
  "eGFR, mL/min/1.73m2 - mean (SD)",           w_meansd(svy,"egfr_cr"),
  "UACR, mg/g - median (IQR)",                 w_medIQR(svy,"uacr"),
  "Diabetes, %",                               w_pct(svy, quote(diabetes==1)),
  "Antihypertensive treatment, %",             w_pct(svy, quote(bp_tx==1)),
  "Statin use, %",                             w_pct(svy, quote(statin==1)),
  "Current smoker, %",                         w_pct(svy, quote(smoking==1)),
  "CKD mechanism: UACR-driven (UACR>=30, eGFR>=60), %", w_pct(svy, quote(ckd_mech=="UACR-driven")),
  "CKD mechanism: eGFR-driven (UACR<30, eGFR<60), %",   w_pct(svy, quote(ckd_mech=="eGFR-driven")),
  "CKD mechanism: Both (UACR>=30, eGFR<60), %",         w_pct(svy, quote(ckd_mech=="Both")))
write_csv(tab1, "results/table1_baseline.csv")

# =============================================================================
# 13. EXPORT ALL TABLES TO HTML (one file per table + a combined index)
# =============================================================================
lbl_thr <- function(x) recode(x, "20"="Fixed 20%", "p75"="75th percentile")
fmt <- function(df) df |> mutate(across(where(is.numeric), ~round(.,1)))

gt_disc <- discordance |> mutate(threshold=lbl_thr(threshold),
    `Discordance % (95% CI)`=sprintf("%.1f (%.1f-%.1f)",pct,pct_lo,pct_hi)) |>
  select(Endpoint=outcome, Threshold=threshold, Stratum=stratum, `Discordance % (95% CI)`)
gt_recl <- reclass |> mutate(threshold=lbl_thr(threshold),
    `N reclassified (95% CI)`=sprintf("%s (%s-%s)",
        format(round(N),big.mark=","),format(round(N_lo),big.mark=","),format(round(N_hi),big.mark=",")),
    `% (95% CI)`=sprintf("%.1f (%.1f-%.1f)",pct,pct_lo,pct_hi)) |>
  select(Endpoint=outcome, Threshold=threshold, Stratum=stratum, `N reclassified (95% CI)`, `% (95% CI)`)
gt_circ <- circularity |> mutate(threshold=lbl_thr(threshold),
    `Reclassified % (95% CI)`=sprintf("%.2f (%.2f-%.2f)",pct,pct_lo,pct_hi)) |>
  select(Subgroup=subgroup, Endpoint=outcome, Threshold=threshold, `Reclassified % (95% CI)`)
gt_comp <- composition |> mutate(threshold=lbl_thr(threshold)) |>
  transmute(Endpoint=outcome, Threshold=threshold, Stratum=stratum,
    `Both lower %`=round(both_lower,1), `Discordant %`=round(discordant,1),
    `Higher 10-y & 30-y %`=round(higher10,1))
gt_med <- medians |> transmute(Endpoint=outcome, Stratum=stratum,
    `Median 10-y risk %`=round(med10,1), `Median 30-y risk %`=round(med30,1))
gt_attr <- attrition |> rename(Step=step, N=n)

tables <- list(
  "Table 1. Baseline characteristics (survey-weighted)" = list(tab1, "table1_baseline"),
  "Table 2. Discordance (low 10-year, high 30-year risk)" = list(gt_disc, "table2_discordance"),
  "Table 3. Reclassification when UACR added" = list(gt_recl, "table3_reclassification"),
  "Table 4. Risk-category composition" = list(gt_comp, "table4_composition"),
  "Table 5. Reclassification by CKD mechanism (circularity)" = list(gt_circ, "table5_circularity"),
  "Table 6. Median predicted risk by stratum" = list(gt_med, "table6_median_risk"),
  "Table 7. Cohort attrition flow" = list(gt_attr, "table7_attrition"))

subt <- "US adults with CKD, ages 30-59, NHANES 1999-2018 | survey-weighted, 95% CIs by Taylor linearization"
for(ttl in names(tables)){
  df <- tables[[ttl]][[1]]; fn <- tables[[ttl]][[2]]
  g <- df |> gt() |> tab_header(title=ttl, subtitle=subt) |>
    opt_table_font(font="Helvetica") |>
    tab_options(table.font.size=px(13), heading.title.font.size=px(15),
                heading.subtitle.font.size=px(11), data_row.padding=px(4))
  gtsave(g, paste0("results/", fn, ".html"))
}

# combined single-page HTML index of all tables
parts <- lapply(names(tables), function(ttl){
  df <- tables[[ttl]][[1]]
  as.character(df |> gt() |> tab_header(title=ttl, subtitle=subt) |>
    opt_table_font(font="Helvetica") |>
    tab_options(table.font.size=px(13)) |> as_raw_html())
})
index_html <- paste0(
  "<!DOCTYPE html><html><head><meta charset='utf-8'><title>PREVENT-CKD Tables</title>",
  "<style>body{font-family:Helvetica,Arial,sans-serif;max-width:1000px;margin:24px auto;padding:0 16px}",
  "h1{font-size:20px}hr{border:none;border-top:1px solid #ddd;margin:28px 0}</style></head><body>",
  "<h1>Discordance between 10- and 30-year PREVENT cardiovascular risk in US adults with CKD</h1>",
  "<p>All tables, survey-weighted (NHANES 1999-2018). preventr ",
  as.character(packageVersion("preventr")), ".</p><hr>",
  paste(parts, collapse="<hr>"), "</body></html>")
writeLines(index_html, "results/all_tables.html")

message("Tables written: results/table1_baseline.csv + 7 HTML tables + all_tables.html")
message("DONE.")

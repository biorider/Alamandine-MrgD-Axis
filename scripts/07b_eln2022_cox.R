# 07b / Step 3: ELN 2022 sensitivity for gate 0.4. NOT the primary analysis (frozen primary = ELN2017 / RISK_MOLECULAR).
# Same GSVA scores (gate04_scores.csv), same Cox specification; only the risk covariate changes; unclassifiable excluded.
suppressPackageStartupMessages({library(survival); library(data.table)})
sc <- fread("results/tables/gate04_scores.csv")
fitc <- function(d, label, cohort) {
  d$score_sd <- as.numeric(scale(d$score)); f <- coxph(Surv(time, event) ~ score_sd + age + risk, data=d)
  s <- summary(f); z <- cox.zph(f)
  data.frame(cohort=cohort, model=label, risk_var="ELN2022_derived", n=nrow(d), events=sum(d$event),
             HR_per_SD=s$conf.int["score_sd",1], lo=s$conf.int["score_sd",3], hi=s$conf.int["score_sd",4],
             p=s$coefficients["score_sd","Pr(>|z|)"], zph_score_p=z$table["score_sd","p"])
}
b <- fread("data/interim/g04_beat_clin.csv"); e <- fread("results/tables/step3_beat_eln2022.csv")
b <- merge(b, e[, .(dbgap_rnaseq_sample=id, ELN2022_derived)], by="dbgap_rnaseq_sample")
b <- b[ELN2022_derived %in% c("Favorable","Intermediate","Adverse") & vitalStatus %in% c("Alive","Dead") & overallSurvival > 0]
bd <- data.frame(id=b$dbgap_rnaseq_sample, time=b$overallSurvival/30.4375, event=as.integer(b$vitalStatus=="Dead"),
                 age=as.numeric(b$ageAtDiagnosis), risk=factor(b$ELN2022_derived, levels=c("Favorable","Intermediate","Adverse")))
t <- fread("data/interim/g04_tcga_clin.csv"); te <- fread("results/tables/step3_tcga_eln2022.csv")
t <- merge(t, te[, .(id, ELN2022_derived)], by="id")
t <- t[ELN2022_derived %in% c("Favorable","Intermediate","Adverse") & OS_MONTHS > 0]
td <- data.frame(id=t$id, time=t$OS_MONTHS, event=as.integer(grepl("^1", t$OS_STATUS)), age=as.numeric(t$AGE),
                 risk=factor(t$ELN2022_derived, levels=c("Favorable","Intermediate","Adverse")))
out <- list()
for (nm in c("RAS8","RAS7")) {
  s1 <- sc[cohort=="BeatAML", .(id, score=get(nm))]; s2 <- sc[cohort=="TCGA-LAML", .(id, score=get(nm))]
  out[[paste0("B",nm)]] <- fitc(merge(bd, s1, by="id"), nm, "BeatAML")
  out[[paste0("T",nm)]] <- fitc(merge(td, s2, by="id"), nm, "TCGA-LAML")
}
res <- rbindlist(out)
for (nm in c("RAS8","RAS7")) res[model==nm, q_BH_2cohorts := p.adjust(p, "BH")]
fwrite(res, "results/tables/step3_eln2022_cox.csv"); print(res, digits=3)
# risk covariate itself should be prognostic (sanity check that the derivation is meaningful)
cat("ELN2022 risk LRT p (Beat, TCGA):", anova(coxph(Surv(time,event)~age, data=bd), coxph(Surv(time,event)~age+risk, data=bd))$`Pr(>|Chi|)`[2],
    anova(coxph(Surv(time,event)~age, data=td), coxph(Surv(time,event)~age+risk, data=td))$`Pr(>|Chi|)`[2], "\n")

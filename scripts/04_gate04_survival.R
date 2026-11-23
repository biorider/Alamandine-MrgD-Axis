# 04: Gate 0.4 — RAS score vs overall survival. Frozen protocol v1 (sha256 466cf2ff...6844), Sections 3.3, 3.4, C (0.4), F.6-F.7.
# DECLARED DEVIATION from kickoff §3.4 (frozen addendum F.7): transplant status is NOT a covariate. Neither cohort has
# transplant dates; a baseline ever/never transplant term would introduce immortal-time bias. Beat AML has no transplant
# field; TCGA TRANSPLANT_TYPE is tabulated descriptively below only. Risk = ELN2017 (Beat) / RISK_MOLECULAR (TCGA), not ELN 2022 (F.6).
# NOTE (Windows/conda R): if BiocSingular crashes on load (Mingw-w64 pseudo-relocation), install env/BiocSingular_stub (see README).
suppressPackageStartupMessages({library(GSVA); library(survival); library(data.table)})
cat("TCGA TRANSPLANT_TYPE (descriptive only):\n"); print(table(fread("data/interim/g04_tcga_clin.csv")$TRANSPLANT_TYPE, useNA="ifany"))
RAS8 <- c("AGT","REN","ACE","ACE2","AGTR1","AGTR2","MAS1","MRGPRD")
RAS7 <- setdiff(RAS8, "MRGPRD")

score_cohort <- function(mat) {
  out <- list(); log <- list()
  for (nm in c("RAS8","RAS7")) {
    gs <- get(nm)
    present <- intersect(gs, rownames(mat))
    nonconst <- present[apply(mat[present,,drop=FALSE], 1, sd) > 0]
    log[[nm]] <- data.frame(set=nm, gene=gs, present=gs %in% present, nonconstant=gs %in% nonconst)
    s <- tryCatch({
      p <- gsvaParam(mat, setNames(list(gs), nm), kcdf="Gaussian", minSize=1, maxDiff=TRUE)
      as.numeric(gsva(p, verbose=FALSE)[nm, ])
    }, error=function(e) { message("GSVA FAILED for ", nm, ": ", conditionMessage(e), " -> pre-specified z-score fallback"); NULL })
    method <- "GSVA"
    if (is.null(s)) { z <- t(scale(t(mat[nonconst,,drop=FALSE]))); s <- colMeans(z, na.rm=TRUE); method <- "zscore_fallback" }
    out[[nm]] <- s; attr(out[[nm]], "method") <- method
  }
  list(scores=out, log=do.call(rbind, log))
}

fit_cox <- function(d, label) {
  d$score_sd <- as.numeric(scale(d$score))
  f <- coxph(Surv(time, event) ~ score_sd + age + risk, data=d)
  sm <- summary(f)$coefficients["score_sd", ]
  ci <- summary(f)$conf.int["score_sd", c("lower .95","upper .95")]
  z <- cox.zph(f)
  data.frame(model=label, n=nrow(d), events=sum(d$event), HR_per_SD=exp(sm["coef"]), lo=ci[1], hi=ci[2],
             p=sm["Pr(>|z|)"], zph_score_p=z$table["score_sd","p"], zph_global_p=z$table["GLOBAL","p"],
             unadj_HR=exp(coef(coxph(Surv(time, event) ~ score_sd, data=d))), row.names=NULL)
}

run <- function(cohort, mfile, cfile, prep) {
  m <- as.matrix(fread(mfile), rownames=1)
  cl <- prep(fread(cfile))
  m <- m[, cl$id]
  sc <- score_cohort(m)
  res <- list(); sco <- data.frame(id=cl$id)
  for (nm in names(sc$scores)) {
    d <- cbind(cl, score=sc$scores[[nm]]); sco[[nm]] <- sc$scores[[nm]]
    r <- fit_cox(d, nm); r$cohort <- cohort; r$score_method <- attr(sc$scores[[nm]], "method"); res[[nm]] <- r
  }
  list(res=do.call(rbind, res), log=cbind(cohort=cohort, sc$log), scores=cbind(cohort=cohort, sco), cor=cor(sco$RAS8, sco$RAS7))
}

prep_tcga <- function(x) {
  x <- x[RISK_MOLECULAR %in% c("Good","Intermediate","Poor") & OS_MONTHS > 0]
  data.frame(id=x$id, time=x$OS_MONTHS, event=as.integer(grepl("^1", x$OS_STATUS)), age=as.numeric(x$AGE),
             risk=factor(x$RISK_MOLECULAR, levels=c("Good","Intermediate","Poor")))
}
prep_beat <- function(x) {
  x <- x[ELN2017 %in% c("Favorable","Intermediate","Adverse") & vitalStatus %in% c("Alive","Dead") & overallSurvival > 0]
  data.frame(id=x$dbgap_rnaseq_sample, time=x$overallSurvival/30.4375, event=as.integer(x$vitalStatus=="Dead"),
             age=as.numeric(x$ageAtDiagnosis), risk=factor(x$ELN2017, levels=c("Favorable","Intermediate","Adverse")))
}

B <- run("BeatAML", "data/interim/g04_beat_log2tpm.csv.gz", "data/interim/g04_beat_clin.csv", prep_beat)
T <- run("TCGA-LAML", "data/interim/g04_tcga_log2tpm.csv.gz", "data/interim/g04_tcga_clin.csv", prep_tcga)
res <- rbind(B$res, T$res)
for (nm in c("RAS8","RAS7")) { i <- res$model == nm; res$q_BH_2cohorts[i] <- p.adjust(res$p[i], "BH") }
fwrite(res, "results/tables/gate04_cox.csv"); fwrite(rbind(B$log, T$log), "results/tables/gate04_gsva_gene_log.csv"); fwrite(rbind(B$scores, T$scores), "results/tables/gate04_scores.csv")
print(res, digits=3); cat("cor(RAS8,RAS7): Beat", round(B$cor,4), " TCGA", round(T$cor,4), "\n")
print(rbind(B$log, T$log))

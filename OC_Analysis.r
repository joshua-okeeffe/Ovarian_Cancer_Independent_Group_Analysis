##Biostatistics: Assignment 1
##Immune checkpoint proteins, treatment arms, and ER status

## ------------------------------------------------------------------
## 0. Setup and data import
## ------------------------------------------------------------------

#Library prep
library(dunn.test)
library(FSA)

#Check current working directory
getwd()

#Read in the dataset
Ass1 <- read.csv("OC_Dataset.csv")

#Quick data exploration
colnames(Ass1)   # variable names
summary(Ass1)    # basic summaries
class(Ass1)      # object class

#Attach for convenience
attach(Ass1)

## ------------------------------------------------------------------
## 1. Difference in pre-treatment expression between treatment arms
##    for each immune checkpoint protein separately
## ------------------------------------------------------------------

#Markers (pre-treatment): PD.1_pre, PD.L1_pre, CTLA.4_pre, SIRP_pre
#Treatment_Arm has levels: "Rucaparib" and "niraparib"

### 1a. Check normality with QQ-plots

qqnorm(PD.1_pre)
qqline(PD.1_pre, col = "red")

qqnorm(PD.L1_pre)
qqline(PD.L1_pre, col = "red")

qqnorm(CTLA.4_pre)
qqline(CTLA.4_pre, col = "red")

qqnorm(SIRP_pre)
qqline(SIRP_pre, col = "red")

### 1b. Split each marker by treatment arm and test normality

#PD1
PD1_ruc <- PD.1_pre[Treatment_Arm == "Rucaparib"]
PD1_nir <- PD.1_pre[Treatment_Arm == "niraparib"]
shapiro.test(PD1_ruc)
shapiro.test(PD1_nir)

#PD-L1
PDL1_ruc <- PD.L1_pre[Treatment_Arm == "Rucaparib"]
PDL1_nir <- PD.L1_pre[Treatment_Arm == "niraparib"]
shapiro.test(PDL1_ruc)
shapiro.test(PDL1_nir)

#CTLA-4
CTLA4_ruc <- CTLA.4_pre[Treatment_Arm == "Rucaparib"]
CTLA4_nir <- CTLA.4_pre[Treatment_Arm == "niraparib"]
shapiro.test(CTLA4_ruc)
shapiro.test(CTLA4_nir)

#SIRP
SIRP_ruc <- SIRP_pre[Treatment_Arm == "Rucaparib"]
SIRP_nir <- SIRP_pre[Treatment_Arm == "niraparib"]
shapiro.test(SIRP_ruc)
shapiro.test(SIRP_nir)

### 1c. Parametric comparison (t test)

#Check variances for PD-1
var(PD1_ruc)
var(PD1_nir)
var.test(PD1_ruc, PD1_nir)   #test equality of variances

#Check variances for PD-L1
var(PDL1_ruc)
var(PDL1_nir)
var.test(PDL1_ruc, PDL1_nir)

#Boxplot for PD-1 and PD-L1 by arm
boxplot(
  PD1_ruc, PD1_nir, PDL1_ruc, PDL1_nir,
  names = c("PD1_ruc", "PD1_nir", "PDL1_ruc", "PDL1_nir"),
  main  = "Protein expression level by treatment arm",
  ylab  = "Expression level",
  col   = c("red", "blue", "red", "blue")
)

#Two-sample t tests
t.test(PD1_ruc, PD1_nir, var.equal = TRUE)
t.test(PDL1_ruc, PDL1_nir, var.equal = TRUE)

### 1d. Non-parametric comparison (Wilcoxon) for non-normal markers

wilcox.test(CTLA4_ruc, CTLA4_nir)
wilcox.test(SIRP_ruc,  SIRP_nir)

#Boxplot for CTLA-4 and SIRP by arm
boxplot(
  CTLA4_ruc, CTLA4_nir, SIRP_ruc, SIRP_nir,
  names = c("CTLA4_ruc", "CTLA4_nir", "SIRP_ruc", "SIRP_nir"),
  main  = "Protein expression level by treatment arm",
  ylab  = "Expression level",
  col   = c("red", "blue", "red", "blue")
)

## ------------------------------------------------------------------
## 2. Association between ER status and treatment response
##    (responders vs non-responders)
## ------------------------------------------------------------------

#Question: Are ER-positive samples more likely to be non-responders?

#Compute counts by ER status and treatment response
n_pos_nonresp <- length(ER_status[ER_status == "positive" & Treatment_Response == "nonresponders"])
n_neg_nonresp <- length(ER_status[ER_status == "negative" & Treatment_Response == "nonresponders"])
n_pos_resp    <- length(ER_status[ER_status == "positive" & Treatment_Response == "responders"])
n_neg_resp    <- length(ER_status[ER_status == "negative" & Treatment_Response == "responders"])

n_pos_nonresp
n_neg_nonresp
n_pos_resp
n_neg_resp

#2x2 contingency table:
#rows: ER status (pos/neg)
#cols: Treatment response (responders/nonresponders)
contable <- matrix(
  c(n_pos_resp, n_pos_nonresp,
    n_neg_resp, n_neg_nonresp),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    "ER status" = c("pos", "neg"),
    "Response"  = c("responders", "nonresponders")
  )
)

print(contable)

#Fisher's exact test (appropriate for small cell counts)
fisher.test(contable)

#Mosaic plot for visualising association
mosaicplot(
  contable,
  main  = "ER status vs treatment response",
  color = TRUE
)

## ------------------------------------------------------------------
## 3. Differences between the four checkpoint proteins
##    within the Rucaparib arm (pre-treatment)
## ------------------------------------------------------------------

# Question: Is there a difference in expression between PD-1, PD-L1,
#           CTLA-4, and SIRP within the Rucaparib arm pre-treatment,
#           and where does that difference lie?

#Combine values for all four markers in the Rucaparib arm
AllVal <- c(PD1_ruc, PDL1_ruc, CTLA4_ruc, SIRP_ruc)

#Check group sizes (should all be equal)
length(PD1_ruc)
length(PDL1_ruc)
length(CTLA4_ruc)
length(SIRP_ruc)

#Create a factor indicating the protein for each value
PD    <- rep("PD1",   times = length(PD1_ruc))
PDL1  <- rep("PDL1",  times = length(PDL1_ruc))
CTLA4 <- rep("CTLA4", times = length(CTLA4_ruc))
SIRP  <- rep("SIRP",  times = length(SIRP_ruc))

AllPro <- c(PD, PDL1, CTLA4, SIRP)

#Boxplot of expression by protein (Rucaparib arm only)
boxplot(
  PD1_ruc, PDL1_ruc, CTLA4_ruc, SIRP_ruc,
  names = c("PD.1", "PD.L1", "CTLA.4", "SIRP"),
  main  = "Expression of proteins in Rucaparib arm (pre-treatment)",
  ylab  = "Expression level",
  col   = c("red", "blue", "green", "yellow")
)

#Create combined long-format data frame
combined <- data.frame(
  Protein = factor(AllPro, levels = c("PD1", "PDL1", "CTLA4", "SIRP")),
  Value   = AllVal
)

#Because normality is doubtful, use Kruskal–Wallis test
kruskal.test(Value ~ Protein, data = combined)

#Post-hoc multiple comparisons using Dunn's test with BH correction
result <- dunnTest(Value ~ Protein, data = combined, method = "bh")
result

## ------------------------------------------------------------------
## 4. Save workspace
## ------------------------------------------------------------------

save(Ass1, file = "OC_Analysis.Rdata")

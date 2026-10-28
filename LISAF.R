############################################################
# LISA v7
# Gastric epithelial identity × biological state
#
# CORE PRINCIPLE
#
# LINEAGE = WHAT THE CELL IS
#
# STATE = WHAT THE CELL IS DOING
#
# DEDIFFERENTIATION = HOW STRONGLY LINEAGE IDENTITY IS
# BEING RELAXED / LOST WHILE MALIGNANT-PLASTICITY
# PROGRAMS BECOME ENRICHED
#
# IMPORTANT:
#
# LISA does NOT assume that dedifferentiation means
# complete disappearance of lineage identity.
#
# Therefore:
#
#   Pit + EMT
#   Pit + Tumor_Associated
#   Chief + Glycolysis
#   Goblet + Hypoxia
#   Parietal + Migration
#
# are all biologically valid.
#
# Dedifferentiation is treated primarily as a CONTINUUM.
#
# Only cells showing:
#
#   1. reduced lineage specificity
#   AND
#   2. strong malignancy / plasticity
#
# are classified as:
#
#   Dedifferentiated_Epithelial
#
# All other cells retain their best epithelial lineage.
############################################################


############################################################
# 0. PACKAGES
############################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(patchwork)


############################################################
# 1. LOAD EPITHELIAL OBJECT
############################################################

Epithelial <- readRDS("C:/Users/somay/Desktop/PhD thesis/Data/single cell/Epithelial.After.Umap.rds")

DefaultAssay(Epithelial) <- "RNA"

Idents(Epithelial) <- "RNA_snn_res.0.6"


cat(
  "\n============================================\n",
  "LISA v7\n",
  "============================================\n"
)

cat(
  "Cells:",
  ncol(Epithelial),
  "\n"
)

cat(
  "Genes:",
  nrow(Epithelial),
  "\n"
)


############################################################
# 2. HOUSEKEEPING / NON-INFORMATIVE GENES
############################################################

housekeeping_genes <- unique(
  c(
    grep("^RPS", rownames(Epithelial), value = TRUE),
    grep("^RPL", rownames(Epithelial), value = TRUE),
    grep("^MT-", rownames(Epithelial), value = TRUE),
    
    "GAPDH",
    "ACTB",
    "ACTG1",
    "TUBB",
    "TUBA1A",
    "TUBA1B",
    "B2M",
    "MALAT1",
    "UBC",
    "UBA52",
    "EEF1A1",
    "HSP90AA1",
    "HSPA8",
    "HSP90AB1"
  )
)

housekeeping_genes <-
  intersect(
    housekeeping_genes,
    rownames(Epithelial)
  )

cat(
  "\nHousekeeping genes excluded:",
  length(housekeeping_genes),
  "\n"
)


############################################################
# 3. EPITHELIAL LINEAGE MARKERS
############################################################

lineage_markers <- list(
  
  ##########################################################
  # PIT
  ##########################################################
  
  Pit = c(
    "MUC5AC",
    "TFF1",
    "GKN1",
    "GKN2"
  ),
  
  ##########################################################
  # CHIEF
  ##########################################################
  
  Chief = c(
    "PGA3",
    "PGA4",
    "PGC",
    "LIPF"
  ),
  
  ##########################################################
  # PARIETAL
  ##########################################################
  
  Parietal = c(
    "ATP4A",
    "ATP4B",
    "GIF",
    "KCNQ1"
  ),
  
  ##########################################################
  # ENTEROENDOCRINE
  ##########################################################
  
  Enteroendocrine = c(
    "CHGA",
    "CHGB",
    "PCSK1",
    "PCSK2",
    "TPH1",
    "TAC1",
    "NEUROG3"
  ),
  
  ##########################################################
  # GOBLET
  ##########################################################
  
  Goblet = c(
    "MUC2",
    "TFF3",
    "SPINK4",
    "CLCA1"
  ),
  
  ##########################################################
  # ENTEROCYTE
  ##########################################################
  
  Enterocyte = c(
    "FABP1",
    "CA1",
    "VIL1",
    "ANPEP",
    "ALPI",
    "APOA1"
  ),
  
  ##########################################################
  # ISTHMUS / STEM
  #
  # IMPORTANT:
  #
  # MKI67 is NOT lineage.
  # Proliferation is a state.
  ##########################################################
  
  Isthmus_Stem = c(
    "LGR5",
    "OLFM4",
    "SMOC2",
    "ASCL2",
    "SOX9"
  ),
  
  ##########################################################
  # GASTRIC MUCOSAL / NECK
  ##########################################################
  
  Gastric_Mucosal = c(
    "MUC6",
    "TFF2",
    "GKN3",
    "MUC1"
  )
)


############################################################
# 4. MALIGNANCY / PLASTICITY PROGRAM
############################################################

dedifferentiated_markers <- c(
  
  # Tumor-associated epithelial phenotype
  "CEACAM5",
  "CEACAM6",
  "FUT3",
  "ST3GAL3",
  
  # Plasticity
  "CD44",
  "LAMC2",
  "SERPINE1",
  
  # Invasion / mesenchymal transition
  "VIM",
  "FN1",
  "ITGA5",
  "ITGAV",
  "MMP14",
  "SNAI1",
  "SNAI2",
  "ZEB1",
  "ZEB2",
  
  # Proliferative malignant phenotype
  "MKI67",
  "TOP2A"
)


############################################################
# 5. BIOLOGICAL STATE PROGRAMS
############################################################

state_markers <- list(
  
  ##########################################################
  # TUMOR-ASSOCIATED
  ##########################################################
  
  Tumor_Associated = c(
    "CEACAM5",
    "CEACAM6",
    "FUT3",
    "FUT2",
    "ST3GAL3",
    "CD44",
    "LAMC2",
    "SERPINE1"
  ),
  
  ##########################################################
  # HYPOXIA
  ##########################################################
  
  Hypoxia = c(
    "CA9",
    "VEGFA",
    "BNIP3",
    "PGK1",
    "ENO1",
    "ALDOA",
    "PDK1",
    "ANGPTL4",
    "ADM",
    "LOX",
    "NDRG1"
  ),
  
  ##########################################################
  # GLYCOLYSIS / LACTATE
  ##########################################################
  
  Glycolysis_Lactate = c(
    "SLC2A1",
    "HK2",
    "PFKM",
    "ALDOA",
    "ENO1",
    "PKM",
    "LDHA",
    "LDHB",
    "SLC16A1",
    "SLC16A3",
    "HCAR1",
    "PDK1"
  ),
  
  ##########################################################
  # OXIDATIVE PHOSPHORYLATION
  ##########################################################
  
  Oxidative_Phosphorylation = c(
    "NDUFA1",
    "NDUFA2",
    "NDUFB8",
    "SDHA",
    "SDHB",
    "UQCRC1",
    "UQCRC2",
    "COX4I1",
    "COX5A",
    "ATP5F1A",
    "ATP5F1B",
    "ATP5MC1"
  ),
  
  ##########################################################
  # DNA REPLICATION
  ##########################################################
  
  DNA_Replication = c(
    "MKI67",
    "PCNA",
    "MCM2",
    "MCM3",
    "MCM4",
    "MCM5",
    "MCM6",
    "MCM7",
    "POLA1",
    "POLD1",
    "POLE",
    "RPA1",
    "FEN1",
    "LIG1",
    "PRIM1",
    "TOP2A",
    "RRM1",
    "RRM2",
    "TK1",
    "TYMS"
  ),
  
  ##########################################################
  # MIGRATION / METASTASIS
  ##########################################################
  
  Migration_Metastasis = c(
    "VIM",
    "ITGA5",
    "ITGB1",
    "ITGAV",
    "FN1",
    "MMP2",
    "MMP9",
    "MMP14",
    "CXCR4",
    "SNAI1",
    "SNAI2",
    "TWIST1",
    "ZEB1",
    "ZEB2"
  ),
  
  ##########################################################
  # EMT / PLASTICITY
  ##########################################################
  
  EMT_Plasticity = c(
    "VIM",
    "FN1",
    "SNAI1",
    "SNAI2",
    "TWIST1",
    "TWIST2",
    "ZEB1",
    "ZEB2",
    "CDH2",
    "ITGA5",
    "ITGB1",
    "ITGAV",
    "LAMC2",
    "SERPINE1",
    "CD44"
  ),
  
  ##########################################################
  # ER STRESS / UPR
  ##########################################################
  
  Stress_UPR = c(
    "HSPA5",
    "XBP1",
    "ATF4",
    "ATF6",
    "DDIT3",
    "ERN1",
    "EIF2AK3",
    "DNAJB9",
    "PDIA4",
    "PDIA6",
    "CALR",
    "CANX",
    "HERPUD1",
    "DERL1",
    "SEL1L"
  ),
  
  ##########################################################
  # INFLAMMATION / INTERFERON
  ##########################################################
  
  Inflammation_Interferon = c(
    "STAT1",
    "STAT2",
    "IRF1",
    "IRF7",
    "ISG15",
    "IFIT1",
    "IFIT2",
    "IFIT3",
    "IFITM1",
    "IFITM2",
    "IFITM3",
    "MX1",
    "OAS1",
    "OAS2",
    "OASL",
    "BST2",
    "CXCL10",
    "CXCL11"
  )
)


############################################################
# 6. REMOVE HOUSEKEEPING / MISSING GENES
############################################################

lineage_markers <-
  lapply(
    lineage_markers,
    function(x) {
      intersect(
        setdiff(
          x,
          housekeeping_genes
        ),
        rownames(Epithelial)
      )
    }
  )


state_markers <-
  lapply(
    state_markers,
    function(x) {
      intersect(
        setdiff(
          x,
          housekeeping_genes
        ),
        rownames(Epithelial)
      )
    }
  )


dedifferentiated_markers <-
  intersect(
    setdiff(
      dedifferentiated_markers,
      housekeeping_genes
    ),
    rownames(Epithelial)
  )


############################################################
# 7. MARKER QC
############################################################

lineage_summary <-
  data.frame(
    lineage = names(lineage_markers),
    detected = lengths(lineage_markers)
  )

state_summary <-
  data.frame(
    state = names(state_markers),
    detected = lengths(state_markers)
  )

print(lineage_summary)
print(state_summary)

cat(
  "\nDedifferentiation markers detected:",
  length(dedifferentiated_markers),
  "\n"
)

if (
  any(
    lengths(lineage_markers) == 0
  )
) {
  stop(
    "At least one lineage marker set has zero detected genes."
  )
}

if (
  any(
    lengths(state_markers) == 0
  )
) {
  stop(
    "At least one state marker set has zero detected genes."
  )
}


############################################################
# 8. REMOVE OLD LISA METADATA
############################################################

old_lisa <-
  grep(
    "^LISA_",
    colnames(Epithelial@meta.data),
    value = TRUE
  )

if (
  length(old_lisa) > 0
) {
  
  Epithelial@meta.data[
    ,
    old_lisa
  ] <- NULL
  
  cat(
    "\nRemoved",
    length(old_lisa),
    "previous LISA columns.\n"
  )
}


############################################################
# 9. LINEAGE MODULE SCORES
############################################################

for (
  nm in names(lineage_markers)
) {
  
  genes <-
    lineage_markers[[nm]]
  
  tmp_name <-
    paste0(
      "LISA_Lineage_",
      nm
    )
  
  Epithelial <-
    AddModuleScore(
      object = Epithelial,
      features = list(genes),
      name = tmp_name,
      assay = "RNA",
      search = FALSE
    )
  
  colnames(
    Epithelial@meta.data
  )[
    colnames(Epithelial@meta.data) ==
      paste0(tmp_name, "1")
  ] <- tmp_name
}


############################################################
# 10. STATE MODULE SCORES
############################################################

for (
  nm in names(state_markers)
) {
  
  genes <-
    state_markers[[nm]]
  
  tmp_name <-
    paste0(
      "LISA_State_",
      nm
    )
  
  Epithelial <-
    AddModuleScore(
      object = Epithelial,
      features = list(genes),
      name = tmp_name,
      assay = "RNA",
      search = FALSE
    )
  
  colnames(
    Epithelial@meta.data
  )[
    colnames(Epithelial@meta.data) ==
      paste0(tmp_name, "1")
  ] <- tmp_name
}


############################################################
# 11. SCORE MATRICES
############################################################

lisa_lineage_scores <-
  grep(
    "^LISA_Lineage_",
    colnames(Epithelial@meta.data),
    value = TRUE
  )

lisa_state_scores <-
  grep(
    "^LISA_State_",
    colnames(Epithelial@meta.data),
    value = TRUE
  )


lineage_score_matrix <-
  as.matrix(
    Epithelial@meta.data[
      ,
      lisa_lineage_scores,
      drop = FALSE
    ]
  )

colnames(lineage_score_matrix) <-
  sub(
    "^LISA_Lineage_",
    "",
    colnames(lineage_score_matrix)
  )


state_score_matrix <-
  as.matrix(
    Epithelial@meta.data[
      ,
      lisa_state_scores,
      drop = FALSE
    ]
  )

colnames(state_score_matrix) <-
  sub(
    "^LISA_State_",
    "",
    colnames(state_score_matrix)
  )


############################################################
# 12. RAW LINEAGE ASSIGNMENT
############################################################

top_idx <-
  max.col(
    lineage_score_matrix,
    ties.method = "first"
  )

Epithelial$LISA_lineage_raw <-
  colnames(lineage_score_matrix)[top_idx]

Epithelial$LISA_lineage_score <-
  lineage_score_matrix[
    cbind(
      seq_len(
        nrow(lineage_score_matrix)
      ),
      top_idx
    )
  ]


############################################################
# 13. SECOND-BEST LINEAGE
############################################################

sorted_lineage_scores <-
  t(
    apply(
      lineage_score_matrix,
      1,
      sort,
      decreasing = TRUE
    )
  )

Epithelial$LISA_second_lineage_score <-
  sorted_lineage_scores[, 2]

Epithelial$LISA_lineage_margin <-
  Epithelial$LISA_lineage_score -
  Epithelial$LISA_second_lineage_score


############################################################
# 14. LINEAGE SPECIFICITY
#
# IMPORTANT:
#
# We now explicitly distinguish:
#
#   lineage identity
#   lineage specificity
#
# A cell can have a lineage assignment while having
# progressively weaker lineage specificity.
############################################################

Epithelial$LISA_lineage_specificity <-
  apply(
    lineage_score_matrix,
    1,
    function(x) {
      
      max(x) -
        median(x)
    }
  )


############################################################
# 15. LINEAGE MARKER SUPPORT
############################################################

all_lineage_genes <-
  unique(
    unlist(
      lineage_markers
    )
  )

lineage_expr <-
  FetchData(
    Epithelial,
    vars = all_lineage_genes
  )


lineage_support <-
  sapply(
    names(lineage_markers),
    function(nm) {
      
      genes <-
        lineage_markers[[nm]]
      
      rowMeans(
        lineage_expr[
          ,
          genes,
          drop = FALSE
        ] > 0
      )
    }
  )

colnames(lineage_support) <-
  names(lineage_markers)


############################################################
# 16. SUPPORT OF ASSIGNED LINEAGE
############################################################

Epithelial$LISA_lineage_marker_fraction <-
  lineage_support[
    cbind(
      seq_len(
        nrow(lineage_support)
      ),
      match(
        Epithelial$LISA_lineage_raw,
        colnames(lineage_support)
      )
    )
  ]


############################################################
# 17. NUMBER OF LINEAGE MARKERS DETECTED
############################################################

lineage_detected <-
  sapply(
    names(lineage_markers),
    function(nm) {
      
      genes <-
        lineage_markers[[nm]]
      
      rowSums(
        lineage_expr[
          ,
          genes,
          drop = FALSE
        ] > 0
      )
    }
  )

colnames(lineage_detected) <-
  names(lineage_markers)


Epithelial$LISA_lineage_markers_detected <-
  lineage_detected[
    cbind(
      seq_len(
        nrow(lineage_detected)
      ),
      match(
        Epithelial$LISA_lineage_raw,
        colnames(lineage_detected)
      )
    )
  ]


############################################################
# 18. LINEAGE CONFIDENCE
############################################################

Epithelial$LISA_lineage_confidence <-
  case_when(
    
    Epithelial$LISA_lineage_margin >=
      0.25 ~
      "High",
    
    Epithelial$LISA_lineage_margin >=
      0.10 ~
      "Intermediate",
    
    TRUE ~
      "Low"
  )


############################################################
# 19. MALIGNANCY / PLASTICITY COMPONENTS
############################################################

tumor_score <-
  Epithelial$LISA_State_Tumor_Associated

emt_score <-
  Epithelial$LISA_State_EMT_Plasticity

migration_score <-
  Epithelial$LISA_State_Migration_Metastasis


############################################################
# 20. STANDARDIZE MALIGNANCY COMPONENTS
############################################################

zscore <- function(x) {
  
  x_sd <-
    sd(
      x,
      na.rm = TRUE
    )
  
  if (
    is.na(x_sd) ||
    x_sd == 0
  ) {
    return(
      rep(
        0,
        length(x)
      )
    )
  }
  
  as.numeric(
    scale(x)
  )
}


tumor_z <-
  zscore(tumor_score)

emt_z <-
  zscore(emt_score)

migration_z <-
  zscore(migration_score)


############################################################
# 21. DEDIFFERENTIATION SCORE
#
# Dedifferentiation is now a CONTINUOUS AXIS.
#
# It combines:
#
#   tumor-associated phenotype
#   EMT/plasticity
#   migration/invasion
#
# without requiring all three to be simultaneously extreme.
############################################################

Epithelial$LISA_dedifferentiation_score <-
  rowMeans(
    cbind(
      tumor_z,
      emt_z,
      migration_z
    ),
    na.rm = TRUE
  )


############################################################
# 22. DEDIFFERENTIATION MARKER SUPPORT
############################################################

dediff_expr <-
  FetchData(
    Epithelial,
    vars = dedifferentiated_markers
  )


Epithelial$LISA_dedifferentiation_marker_fraction <-
  rowMeans(
    dediff_expr > 0
  )


############################################################
# 23. DEDIFFERENTIATION COMPONENTS
############################################################

Epithelial$LISA_dediff_tumor_z <-
  tumor_z

Epithelial$LISA_dediff_emt_z <-
  emt_z

Epithelial$LISA_dediff_migration_z <-
  migration_z


############################################################
# 24. QUANTILE-BASED DEDIFFERENTIATION CATEGORIES
#
# IMPORTANT:
#
# These categories describe the distribution of the
# continuous dedifferentiation axis.
#
# They are NOT equivalent to biological lineage labels.
############################################################

dediff_q33 <-
  quantile(
    Epithelial$LISA_dedifferentiation_score,
    0.33,
    na.rm = TRUE
  )

dediff_q67 <-
  quantile(
    Epithelial$LISA_dedifferentiation_score,
    0.67,
    na.rm = TRUE
  )

Epithelial$LISA_dedifferentiation_state <-
  case_when(
    
    Epithelial$LISA_dedifferentiation_score >=
      dediff_q67 ~
      "High",
    
    Epithelial$LISA_dedifferentiation_score >=
      dediff_q33 ~
      "Intermediate",
    
    TRUE ~
      "Low"
  )


############################################################
# 25. DEFINE LINEAGE-POOR CELLS
#
# We deliberately use a less restrictive threshold than v6.
#
# Bottom 30% of lineage specificity.
#
# This is NOT itself dedifferentiation.
############################################################

lineage_specificity_threshold <-
  quantile(
    Epithelial$LISA_lineage_specificity,
    0.30,
    na.rm = TRUE
  )


Epithelial$LISA_lineage_poor_flag <-
  ifelse(
    Epithelial$LISA_lineage_specificity <=
      lineage_specificity_threshold,
    "Yes",
    "No"
  )


############################################################
# 26. STRONG DEDIFFERENTIATION THRESHOLD
#
# Top 30% of the continuous dedifferentiation score.
############################################################

dediff_high_threshold <-
  quantile(
    Epithelial$LISA_dedifferentiation_score,
    0.70,
    na.rm = TRUE
  )


############################################################
# 27. FINAL DEDIFFERENTIATED CRITERION
#
# A cell is classified as:
#
#   Dedifferentiated_Epithelial
#
# only when:
#
#   A. lineage specificity is low
#
# AND
#
#   B. dedifferentiation score is high
#
# This is deliberately simpler and more biologically
# interpretable than the v6 three-way intersection.
############################################################

dediff_candidate <-
  
  Epithelial$LISA_lineage_specificity <=
  lineage_specificity_threshold &
  
  Epithelial$LISA_dedifferentiation_score >=
  dediff_high_threshold


############################################################
# 28. FINAL LINEAGE
############################################################

Epithelial$LISA_lineage_final <-
  Epithelial$LISA_lineage_raw


Epithelial$LISA_lineage_final[
  dediff_candidate
] <-
  "Dedifferentiated_Epithelial"


############################################################
# 29. DEDIFFERENTIATION FLAG
############################################################

Epithelial$LISA_dedifferentiated_flag <-
  ifelse(
    dediff_candidate,
    "Yes",
    "No"
  )


############################################################
# 30. TOP BIOLOGICAL STATE
#
# IMPORTANT:
#
# This is only a descriptive "dominant state".
#
# A cell can simultaneously have multiple active states.
############################################################

state_idx <-
  max.col(
    state_score_matrix,
    ties.method = "first"
  )


Epithelial$LISA_top_state <-
  colnames(state_score_matrix)[state_idx]


Epithelial$LISA_top_state_score <-
  state_score_matrix[
    cbind(
      seq_len(
        nrow(state_score_matrix)
      ),
      state_idx
    )
  ]


############################################################
# 31. SECOND-BEST STATE
############################################################

sorted_state_scores <-
  t(
    apply(
      state_score_matrix,
      1,
      sort,
      decreasing = TRUE
    )
  )


Epithelial$LISA_second_state_score <-
  sorted_state_scores[, 2]


Epithelial$LISA_state_margin <-
  Epithelial$LISA_top_state_score -
  Epithelial$LISA_second_state_score


############################################################
# 32. STATE CONFIDENCE
############################################################

Epithelial$LISA_state_confidence <-
  case_when(
    
    Epithelial$LISA_state_margin >=
      0.20 ~
      "High",
    
    Epithelial$LISA_state_margin >=
      0.10 ~
      "Intermediate",
    
    TRUE ~
      "Low"
  )


############################################################
# 33. LINEAGE × STATE SUMMARY
############################################################

lineage_state_long <-
  
  Epithelial@meta.data %>%
  
  dplyr::select(
    LISA_lineage_final,
    dplyr::all_of(
      lisa_state_scores
    )
  ) %>%
  
  tidyr::pivot_longer(
    cols =
      dplyr::all_of(
        lisa_state_scores
      ),
    names_to = "state",
    values_to = "score"
  ) %>%
  
  dplyr::mutate(
    state =
      sub(
        "^LISA_State_",
        "",
        state
      )
  ) %>%
  
  dplyr::group_by(
    LISA_lineage_final,
    state
  ) %>%
  
  dplyr::summarise(
    
    mean_score =
      mean(
        score,
        na.rm = TRUE
      ),
    
    median_score =
      median(
        score,
        na.rm = TRUE
      ),
    
    n =
      dplyr::n(),
    
    .groups = "drop"
  )


############################################################
# 34. LINEAGE × CLUSTER
############################################################

lineage_cluster <-
  
  Epithelial@meta.data %>%
  
  dplyr::count(
    RNA_snn_res.0.6,
    LISA_lineage_final,
    name = "n"
  ) %>%
  
  dplyr::group_by(
    RNA_snn_res.0.6
  ) %>%
  
  dplyr::mutate(
    percent =
      n /
      sum(n) *
      100
  ) %>%
  
  dplyr::ungroup()


############################################################
# 35. STATE × CLUSTER
############################################################

state_cluster <-
  
  Epithelial@meta.data %>%
  
  dplyr::count(
    RNA_snn_res.0.6,
    LISA_top_state,
    name = "n"
  ) %>%
  
  dplyr::group_by(
    RNA_snn_res.0.6
  ) %>%
  
  dplyr::mutate(
    percent =
      n /
      sum(n) *
      100
  ) %>%
  
  dplyr::ungroup()


############################################################
# 36. LINEAGE × STATE × CLUSTER
############################################################

lineage_state_cluster <-
  
  Epithelial@meta.data %>%
  
  dplyr::count(
    RNA_snn_res.0.6,
    LISA_lineage_final,
    LISA_top_state,
    name = "n"
  ) %>%
  
  dplyr::group_by(
    RNA_snn_res.0.6
  ) %>%
  
  dplyr::mutate(
    percent =
      n /
      sum(n) *
      100
  ) %>%
  
  dplyr::ungroup()


############################################################
# 37. LINEAGE × CONDITION × CLUSTER
############################################################

lineage_condition_cluster <-
  
  Epithelial@meta.data %>%
  
  dplyr::count(
    RNA_snn_res.0.6,
    LISA_lineage_final,
    condition,
    name = "n"
  ) %>%
  
  dplyr::group_by(
    RNA_snn_res.0.6,
    LISA_lineage_final
  ) %>%
  
  dplyr::mutate(
    percent =
      n /
      sum(n) *
      100
  ) %>%
  
  dplyr::ungroup()


############################################################
# 38. STATE × CONDITION × CLUSTER
############################################################

state_condition_cluster <-
  
  Epithelial@meta.data %>%
  
  dplyr::count(
    RNA_snn_res.0.6,
    LISA_top_state,
    condition,
    name = "n"
  ) %>%
  
  dplyr::group_by(
    RNA_snn_res.0.6,
    LISA_top_state
  ) %>%
  
  dplyr::mutate(
    percent =
      n /
      sum(n) *
      100
  ) %>%
  
  dplyr::ungroup()


############################################################
# 39. CONDITION × CLUSTER
############################################################

condition_cluster <-
  
  Epithelial@meta.data %>%
  
  dplyr::count(
    RNA_snn_res.0.6,
    condition,
    name = "n"
  ) %>%
  
  dplyr::group_by(
    RNA_snn_res.0.6
  ) %>%
  
  dplyr::mutate(
    percent =
      n /
      sum(n) *
      100
  ) %>%
  
  dplyr::ungroup()


############################################################
# 40. CONDITION WIDE
############################################################

condition_cluster_wide <-
  
  condition_cluster %>%
  
  dplyr::select(
    RNA_snn_res.0.6,
    condition,
    percent
  ) %>%
  
  tidyr::pivot_wider(
    names_from = condition,
    values_from = percent,
    values_fill = 0
  )


############################################################
# 41. LINEAGE SCORE BY CLUSTER
############################################################

lineage_score_by_cluster <-
  
  FetchData(
    Epithelial,
    vars = c(
      "RNA_snn_res.0.6",
      lisa_lineage_scores
    )
  ) %>%
  
  tidyr::pivot_longer(
    cols =
      dplyr::all_of(
        lisa_lineage_scores
      ),
    names_to = "lineage",
    values_to = "score"
  ) %>%
  
  dplyr::mutate(
    lineage =
      sub(
        "^LISA_Lineage_",
        "",
        lineage
      )
  ) %>%
  
  dplyr::group_by(
    RNA_snn_res.0.6,
    lineage
  ) %>%
  
  dplyr::summarise(
    
    mean_score =
      mean(
        score,
        na.rm = TRUE
      ),
    
    median_score =
      median(
        score,
        na.rm = TRUE
      ),
    
    max_score =
      max(
        score,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


############################################################
# 42. RANK LINEAGES WITHIN CLUSTERS
############################################################

lineage_cluster_rank <-
  
  lineage_score_by_cluster %>%
  
  dplyr::group_by(
    RNA_snn_res.0.6
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(mean_score),
    .by_group = TRUE
  ) %>%
  
  dplyr::mutate(
    
    rank =
      dplyr::row_number(),
    
    relative_to_top =
      mean_score /
      max(mean_score)
  ) %>%
  
  dplyr::ungroup()


############################################################
# 43. DEDIFFERENTIATION BY CLUSTER
############################################################

dediff_cluster <-
  
  Epithelial@meta.data %>%
  
  dplyr::group_by(
    RNA_snn_res.0.6
  ) %>%
  
  dplyr::summarise(
    
    n_cells =
      dplyr::n(),
    
    mean_dediff_score =
      mean(
        LISA_dedifferentiation_score,
        na.rm = TRUE
      ),
    
    median_dediff_score =
      median(
        LISA_dedifferentiation_score,
        na.rm = TRUE
      ),
    
    high_dediff_cells =
      sum(
        LISA_dedifferentiation_state ==
          "High"
      ),
    
    dedifferentiated_cells =
      sum(
        LISA_dedifferentiated_flag ==
          "Yes"
      ),
    
    dedifferentiated_percent =
      100 *
      dedifferentiated_cells /
      n_cells,
    
    .groups = "drop"
  )


############################################################
# 44. LINEAGE × CONDITION
############################################################

lineage_condition <-
  
  Epithelial@meta.data %>%
  
  dplyr::count(
    condition,
    LISA_lineage_final,
    name = "n"
  ) %>%
  
  dplyr::group_by(
    condition
  ) %>%
  
  dplyr::mutate(
    percent =
      n /
      sum(n) *
      100
  ) %>%
  
  dplyr::ungroup()


############################################################
# 45. STATE × LINEAGE × CONDITION
############################################################

state_lineage_condition <-
  
  Epithelial@meta.data %>%
  
  dplyr::select(
    condition,
    LISA_lineage_final,
    dplyr::all_of(
      lisa_state_scores
    )
  ) %>%
  
  tidyr::pivot_longer(
    cols =
      dplyr::all_of(
        lisa_state_scores
      ),
    names_to = "state",
    values_to = "score"
  ) %>%
  
  dplyr::mutate(
    state =
      sub(
        "^LISA_State_",
        "",
        state
      )
  ) %>%
  
  dplyr::group_by(
    condition,
    LISA_lineage_final,
    state
  ) %>%
  
  dplyr::summarise(
    
    mean_score =
      mean(
        score,
        na.rm = TRUE
      ),
    
    median_score =
      median(
        score,
        na.rm = TRUE
      ),
    
    n =
      dplyr::n(),
    
    .groups = "drop"
  )


############################################################
# 46. DEDIFFERENTIATION TRAJECTORY TABLE
#
# This is particularly important for your hypothesis.
#
# It allows you to examine whether lineage specificity
# decreases as dedifferentiation score increases.
############################################################

dedifferentiation_trajectory <-
  
  Epithelial@meta.data %>%
  
  dplyr::select(
    LISA_lineage_raw,
    LISA_lineage_score,
    LISA_lineage_specificity,
    LISA_lineage_margin,
    LISA_lineage_marker_fraction,
    LISA_dedifferentiation_score,
    LISA_dedifferentiation_state,
    LISA_dedifferentiation_marker_fraction,
    LISA_State_Tumor_Associated,
    LISA_State_EMT_Plasticity,
    LISA_State_Migration_Metastasis
  )


############################################################
# 47. UMAP — FINAL LINEAGE
############################################################

p_lineage_umap <-
  
  DimPlot(
    Epithelial,
    reduction = "umap",
    group.by = "LISA_lineage_final",
    label = TRUE,
    repel = TRUE,
    raster = FALSE,
    pt.size = 0.5
  ) +
  
  labs(
    title =
      "Gastric Epithelial Lineage"
  ) +
  
  theme_classic() +
  
  theme(
    plot.title =
      element_text(
        face = "bold",
        hjust = 0.5,
        size = 15
      ),
    
    legend.title =
      element_blank()
  )

print(p_lineage_umap)


############################################################
# 48. UMAP — LINEAGE CONFIDENCE
############################################################

p1 <-
  
  DimPlot(
    Epithelial,
    reduction = "umap",
    group.by = "LISA_lineage_confidence",
    raster = FALSE,
    pt.size = 0.5
  ) +
  
  scale_color_manual(
    values = c(
      "High" = "darkblue",
      "Intermediate" = "grey",
      "Low" = "pink"
    )
  ) +
  
  labs(
    title =
      "LISA Lineage Confidence"
  ) +
  
  theme_classic() +
  
  theme(
    plot.title =
      element_text(
        face = "bold",
        hjust = 0.5
      ),
    
    legend.title =
      element_blank()
  )

print(p1)


############################################################
# 49. UMAP — DEDIFFERENTIATION SCORE
############################################################

p2 <-
  
  FeaturePlot(
    Epithelial,
    features =
      "LISA_dedifferentiation_score",
    reduction = "umap",
    raster = FALSE,
    pt.size = 0.5,
    cols =
      c(
        "white",
        "darkblue"
      )
  ) +
  
  labs(
    title =
      "LISA Dedifferentiation Continuum"
  ) +
  
  theme_classic() +
  
  theme(
    plot.title =
      element_text(
        face = "bold",
        hjust = 0.5
      )
  )

print(p2)


############################################################
# 50. UMAP — DEDIFFERENTIATION STATE
############################################################

p3 <-
  
  DimPlot(
    Epithelial,
    reduction = "umap",
    group.by =
      "LISA_dedifferentiation_state",
    raster = FALSE,
    pt.size = 0.5
  ) +
  
  labs(
    title =
      "LISA Dedifferentiation State"
  ) +
  
  theme_classic() +
  
  theme(
    plot.title =
      element_text(
        face = "bold",
        hjust = 0.5
      ),
    
    legend.title =
      element_blank()
  )

print(p3)


############################################################
# 51. UMAP — FINAL DEDIFFERENTIATED CELLS
############################################################
p4 <-
  DimPlot(
    Epithelial,
    reduction = "umap",
    group.by = "LISA_dedifferentiated_flag",
    cols = c(
      "Yes" = "#D73027",
      "No"  = "grey"
    ),
    raster = FALSE,
    pt.size = 0.5
  ) +
  
  labs(
    title =
      "Strong Lineage-poor Malignant/Plasticity Phenotype"
  ) +
  
  theme_classic() +
  
  theme(
    plot.title =
      element_text(
        face = "bold",
        hjust = 0.5
      ),
    
    legend.title =
      element_blank()
  )

print(p4)
############################################################

p1+p2+p3+p4
############################################################
# 52. UMAP — BIOLOGICAL STATES
############################################################

p_states <-
  FeaturePlot(
    Epithelial,
    features = lisa_state_scores,
    reduction = "umap",
    ncol = 3,
    raster = FALSE,
    cols = c(
      "lightgrey",
      "pink",
      "darkred"
    ),
    pt.size = 0.5
  )

print(p_states)


############################################################
# 53. LINEAGE × STATE HEATMAP
############################################################

p_lineage_state <-
  
  ggplot(
    lineage_state_long,
    aes(
      x = state,
      y = LISA_lineage_final,
      fill = mean_score
    )
  ) +
  
  geom_tile(
    color = "white",
    linewidth = 0.4
  ) +
  
  geom_text(
    aes(
      label =
        paste0(
          "n=",
          format(
            n,
            big.mark = ","
          )
        )
    ),
    size = 3
  ) +
  
  scale_fill_gradient2(
    low = "navy",
    mid = "white",
    high = "firebrick",
    midpoint = 0
  ) +
  
  labs(
    x = "Biological state",
    y = "Epithelial lineage",
    fill = "Mean LISA\nscore",
    title =
      "LISA: epithelial identity × biological state"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1,
        face = "bold"
      ),
    
    axis.text.y =
      element_text(
        face = "bold"
      ),
    
    plot.title =
      element_text(
        face = "bold",
        hjust = 0.5
      )
  )

print(p_lineage_state)


############################################################
# 54. KEY PLOT:
# LINEAGE SPECIFICITY × DEDIFFERENTIATION
############################################################

p_dediff_relationship <-
  
  ggplot(
    Epithelial@meta.data,
    aes(
      x =
        LISA_dedifferentiation_score,
      y =
        LISA_lineage_specificity
    )
  ) +
  
  geom_point(
    alpha = 0.25,
    size = 0.4
  ) +
  
  geom_smooth(
    method = "loess",
    se = TRUE
  ) +
  
  labs(
    x =
      "Dedifferentiation score",
    
    y =
      "Lineage specificity",
    
    title =
      "Relationship between lineage specificity and dedifferentiation"
  ) +
  
  theme_classic() +
  
  theme(
    plot.title =
      element_text(
        face = "bold",
        hjust = 0.5
      )
  )

print(p_dediff_relationship)
############################################################
cor.test(
  Epithelial$LISA_dedifferentiation_score,
  Epithelial$LISA_lineage_specificity,
  method = "spearman"
)
Epithelial$dediff_group <- cut(
  Epithelial$LISA_dedifferentiation_score,
  breaks = c(
    -Inf,
    0,
    1,
    2,
    Inf
  ),
  labels = c(
    "Low",
    "Intermediate",
    "High",
    "Very_high"
  )
)

library(ggplot2)

# Prepare plotting data
plot_df <- data.frame(
  Dedifferentiation = Epithelial$LISA_dedifferentiation_score,
  Lineage_specificity = Epithelial$LISA_lineage_specificity,
  Dediff_group = Epithelial$dediff_group
)

# Remove missing values
plot_df <- plot_df[complete.cases(plot_df), ]

# Spearman correlation
cor_res <- cor.test(
  plot_df$Dedifferentiation,
  plot_df$Lineage_specificity,
  method = "spearman"
)

ggplot(
  plot_df,
  aes(
    x = Dedifferentiation,
    y = Lineage_specificity,
    color = Dediff_group
  )
) +
  geom_point(
    alpha = 0.35,
    size = 0.8
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0(
      "Spearman \u03c1 = ",
      round(cor_res$estimate, 3),
      "\nP = ",
      format.pval(cor_res$p.value, digits = 3, eps = 1e-10)
    ),
    hjust = 1.1,
    vjust = 1.5,
    size = 4
  ) +
  labs(
    x = "LISA dedifferentiation score",
    y = "LISA lineage specificity",
    color = "Dedifferentiation group"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )


############################################################
# 55. LINEAGE × DEDIFFERENTIATION
############################################################

lineage_dediff <-
  
  Epithelial@meta.data %>%
  
  dplyr::group_by(
    LISA_lineage_raw
  ) %>%
  
  dplyr::summarise(
    
    n =
      dplyr::n(),
    
    mean_dediff =
      mean(
        LISA_dedifferentiation_score,
        na.rm = TRUE
      ),
    
    median_dediff =
      median(
        LISA_dedifferentiation_score,
        na.rm = TRUE
      ),
    
    high_dediff_percent =
      mean(
        LISA_dedifferentiation_state ==
          "High"
      ) *
      100,
    
    .groups = "drop"
  )

print(lineage_dediff)
############################################################


library(ggplot2)
library(patchwork)

p1 <- ggplot(
  lineage_dediff,
  aes(
    x = reorder(LISA_lineage_raw, mean_dediff),
    y = mean_dediff
  )
) +
  geom_col(width = 0.7) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Mean dedifferentiation score"
  ) +
  theme_classic(base_size = 13)

p2 <- ggplot(
  lineage_dediff,
  aes(
    x = reorder(LISA_lineage_raw, high_dediff_percent),
    y = high_dediff_percent
  )
) +
  geom_col(width = 0.7) +
  coord_flip() +
  labs(
    x = "Epithelial lineage",
    y = "High dedifferentiation (%)"
  ) +
  theme_classic(base_size = 13)

p1 + p2
############################################################
# 56. SAVE OBJECT
############################################################

saveRDS(Epithelial,"C:/Users/somay/Desktop/PhD thesis/Data/single cell/Epithelial.LISA.v7.rds")


############################################################
# 57. SAVE TABLES
############################################################

write.csv(
  lineage_summary,
  "E:/LISA_v7_lineage_marker_summary.csv",
  row.names = FALSE
)

write.csv(
  state_summary,
  "E:/LISA_v7_state_marker_summary.csv",
  row.names = FALSE
)

write.csv(
  lineage_cluster,
  "E:/LISA_v7_lineage_cluster.csv",
  row.names = FALSE
)

write.csv(
  state_cluster,
  "E:/LISA_v7_state_cluster.csv",
  row.names = FALSE
)

write.csv(
  lineage_state_cluster,
  "E:/LISA_v7_lineage_state_cluster.csv",
  row.names = FALSE
)

write.csv(
  lineage_condition_cluster,
  "E:/LISA_v7_lineage_condition_cluster.csv",
  row.names = FALSE
)

write.csv(
  state_condition_cluster,
  "E:/LISA_v7_state_condition_cluster.csv",
  row.names = FALSE
)

write.csv(
  condition_cluster_wide,
  "E:/LISA_v7_condition_cluster.csv",
  row.names = FALSE
)

write.csv(
  lineage_score_by_cluster,
  "E:/LISA_v7_lineage_scores_by_cluster.csv",
  row.names = FALSE
)

write.csv(
  lineage_cluster_rank,
  "E:/LISA_v7_lineage_cluster_rank.csv",
  row.names = FALSE
)

write.csv(
  lineage_state_long,
  "E:/LISA_v7_lineage_state_summary.csv",
  row.names = FALSE
)

write.csv(
  lineage_condition,
  "E:/LISA_v7_lineage_condition.csv",
  row.names = FALSE
)

write.csv(
  state_lineage_condition,
  "E:/LISA_v7_state_lineage_condition.csv",
  row.names = FALSE
)

write.csv(
  dediff_cluster,
  "E:/LISA_v7_dedifferentiation_by_cluster.csv",
  row.names = FALSE
)

write.csv(
  dedifferentiation_trajectory,
  "E:/LISA_v7_dedifferentiation_trajectory.csv",
  row.names = FALSE
)

write.csv(
  lineage_dediff,
  "E:/LISA_v7_lineage_dedifferentiation.csv",
  row.names = FALSE
)


############################################################
# 58. FINAL QC REPORT
############################################################

cat(
  "\n\n============================================\n",
  "LISA v7 FINAL QC REPORT\n",
  "============================================\n"
)

cat(
  "\nTotal cells:",
  ncol(Epithelial),
  "\n"
)

cat(
  "\nFinal lineages:\n"
)

print(
  sort(
    table(
      Epithelial$LISA_lineage_final
    ),
    decreasing = TRUE
  )
)

cat(
  "\nFinal lineage percentages:\n"
)

print(
  round(
    prop.table(
      table(
        Epithelial$LISA_lineage_final
      )
    ) *
      100,
    2
  )
)

cat(
  "\nLineage confidence:\n"
)

print(
  table(
    Epithelial$LISA_lineage_confidence
  )
)

cat(
  "\nDedifferentiation state:\n"
)

print(
  table(
    Epithelial$LISA_dedifferentiation_state
  )
)

cat(
  "\nStrong dedifferentiated epithelial cells:",
  sum(
    Epithelial$LISA_dedifferentiated_flag ==
      "Yes"
  ),
  "\n"
)

cat(
  "\nTop biological states:\n"
)

print(
  sort(
    table(
      Epithelial$LISA_top_state
    ),
    decreasing = TRUE
  )
)

cat(
  "\n============================================\n",
  "LISA v7 COMPLETE\n",
  "============================================\n"
)
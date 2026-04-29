############################################################
# Broad stage I-IV cell-type composition analysis
# Dataset: GSE131907 lung cancer metadata
############################################################

proj_dir <- "/Users/rongrong/Desktop/JHU/Academics/Course/大三/大三下/Monte Carlo Methods/Project/"
setwd(proj_dir)

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

############################################################
# 1. Load corrected metadata
############################################################

meta_file <- file.path(
  proj_dir,
  "GSE131907_cell_annotation_with_clinical_metadata.csv"
)

meta <- read.csv(
  meta_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Metadata dimensions:\n")
print(dim(meta))

cat("\nColumn names:\n")
print(colnames(meta))

############################################################
# 2. Create broad stage label: I, II, III, IV
############################################################

detailed_stage_levels <- c("IA", "IA2", "IA3", "IB", "IIA", "IIB", "IIIA", "IV")
broad_stage_levels <- c("I", "II", "III", "IV")

meta$Stage <- factor(meta$Stage, levels = detailed_stage_levels)

meta <- meta %>%
  mutate(
    Stage_broad = case_when(
      Stage %in% c("IA", "IA2", "IA3", "IB") ~ "I",
      Stage %in% c("IIA", "IIB") ~ "II",
      Stage %in% c("IIIA") ~ "III",
      Stage %in% c("IV") ~ "IV",
      TRUE ~ NA_character_
    )
  )

meta$Stage_broad <- factor(meta$Stage_broad, levels = broad_stage_levels)

cat("\nCell counts by detailed Stage:\n")
print(table(meta$Stage, useNA = "ifany"))

cat("\nCell counts by broad Stage:\n")
print(table(meta$Stage_broad, useNA = "ifany"))

cat("\nPatient counts by broad Stage:\n")
patient_stage_broad <- meta %>%
  distinct(Patient, Stage_broad) %>%
  count(Stage_broad, name = "n_patients")

print(patient_stage_broad)

############################################################
# 3. Save metadata with broad stage
############################################################

write.csv(
  meta,
  file = file.path(
    proj_dir,
    "GSE131907_cell_annotation_with_clinical_metadata_stage_broad.csv"
  ),
  row.names = FALSE
)

cat("\nSaved metadata with broad stage:\n")
cat("GSE131907_cell_annotation_with_clinical_metadata_stage_broad.csv\n")

############################################################
# 4. Choose cell-type column
############################################################

# In your metadata, Cell_type was detected successfully.
# You can change this to "Cell_type.refined" or "Cell_subtype"
# if you want a more detailed annotation later.

celltype_col <- "Cell_type"

if (!(celltype_col %in% colnames(meta))) {
  stop("Cell_type column not found. Please check colnames(meta).")
}

cat("\nUsing cell-type column:\n")
print(celltype_col)

############################################################
# 5. Cell counts and proportions by broad stage
############################################################

broad_stage_celltype_counts <- meta %>%
  count(Stage_broad, .data[[celltype_col]], name = "n_cells") %>%
  group_by(Stage_broad) %>%
  mutate(
    stage_total_cells = sum(n_cells),
    proportion_within_stage = n_cells / stage_total_cells
  ) %>%
  ungroup()

write.csv(
  broad_stage_celltype_counts,
  file = file.path(
    proj_dir,
    "broad_stage_celltype_counts_and_proportions.csv"
  ),
  row.names = FALSE
)

cat("\nSaved:\n")
cat("broad_stage_celltype_counts_and_proportions.csv\n")

############################################################
# 6. Patient-level cell-type proportions using broad stage
############################################################

patient_broad_celltype_counts <- meta %>%
  count(Patient, Stage_broad, Tissue, .data[[celltype_col]], name = "n_cells") %>%
  group_by(Patient) %>%
  mutate(
    patient_total_cells = sum(n_cells),
    proportion_within_patient = n_cells / patient_total_cells
  ) %>%
  ungroup()

write.csv(
  patient_broad_celltype_counts,
  file = file.path(
    proj_dir,
    "patient_broad_stage_celltype_counts_and_proportions.csv"
  ),
  row.names = FALSE
)

cat("patient_broad_stage_celltype_counts_and_proportions.csv\n")

############################################################
# 7. Tissue-aware patient-level proportions
############################################################

patient_tissue_broad_celltype_counts <- meta %>%
  count(Patient, Stage_broad, Tissue, .data[[celltype_col]], name = "n_cells") %>%
  group_by(Patient, Tissue) %>%
  mutate(
    patient_tissue_total_cells = sum(n_cells),
    proportion_within_patient_tissue = n_cells / patient_tissue_total_cells
  ) %>%
  ungroup()

write.csv(
  patient_tissue_broad_celltype_counts,
  file = file.path(
    proj_dir,
    "patient_tissue_broad_stage_celltype_counts_and_proportions.csv"
  ),
  row.names = FALSE
)

cat("patient_tissue_broad_stage_celltype_counts_and_proportions.csv\n")

############################################################
# 8. Patient-level matrix for modeling
############################################################
# Rows = patients
# Columns = cell-type proportions
# Label = broad stage

patient_celltype_matrix <- patient_broad_celltype_counts %>%
  select(
    Patient,
    Stage_broad,
    Tissue,
    cell_type = .data[[celltype_col]],
    proportion_within_patient
  ) %>%
  group_by(Patient, Stage_broad, cell_type) %>%
  summarise(
    proportion_within_patient = sum(proportion_within_patient),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = cell_type,
    values_from = proportion_within_patient,
    values_fill = 0
  )

write.csv(
  patient_celltype_matrix,
  file = file.path(
    proj_dir,
    "patient_celltype_proportion_matrix_broad_stage.csv"
  ),
  row.names = FALSE
)

cat("patient_celltype_proportion_matrix_broad_stage.csv\n")

############################################################
# 9. Patient-by-tissue matrix for modeling
############################################################
# This keeps tissue origin separate.
# Useful because tissue origin and stage are strongly related.

patient_tissue_celltype_matrix <- patient_tissue_broad_celltype_counts %>%
  select(
    Patient,
    Stage_broad,
    Tissue,
    cell_type = .data[[celltype_col]],
    proportion_within_patient_tissue
  ) %>%
  pivot_wider(
    names_from = cell_type,
    values_from = proportion_within_patient_tissue,
    values_fill = 0
  )

write.csv(
  patient_tissue_celltype_matrix,
  file = file.path(
    proj_dir,
    "patient_tissue_celltype_proportion_matrix_broad_stage.csv"
  ),
  row.names = FALSE
)

cat("patient_tissue_celltype_proportion_matrix_broad_stage.csv\n")

############################################################
# 10. Plot: cell-type composition by broad stage
############################################################

p1 <- ggplot(
  broad_stage_celltype_counts,
  aes(
    x = Stage_broad,
    y = proportion_within_stage,
    fill = .data[[celltype_col]]
  )
) +
  geom_col(position = "stack") +
  labs(
    title = "Cell-type composition by broad cancer stage",
    x = "Broad cancer stage",
    y = "Proportion of cells within stage",
    fill = "Cell type"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(
    proj_dir,
    "plot_celltype_composition_by_broad_stage.png"
  ),
  plot = p1,
  width = 8,
  height = 6,
  dpi = 300
)

cat("\nSaved plot:\n")
cat("plot_celltype_composition_by_broad_stage.png\n")

############################################################
# 11. Plot: patient-level composition by broad stage
############################################################

p2 <- ggplot(
  patient_broad_celltype_counts,
  aes(
    x = Patient,
    y = proportion_within_patient,
    fill = .data[[celltype_col]]
  )
) +
  geom_col(position = "stack") +
  facet_grid(. ~ Stage_broad, scales = "free_x", space = "free_x") +
  labs(
    title = "Patient-level cell-type composition by broad stage",
    x = "Patient",
    y = "Proportion of cells within patient",
    fill = "Cell type"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(
    proj_dir,
    "plot_patient_celltype_composition_by_broad_stage.png"
  ),
  plot = p2,
  width = 12,
  height = 6,
  dpi = 300
)

cat("plot_patient_celltype_composition_by_broad_stage.png\n")

############################################################
# 12. Plot: tissue origins by broad stage
############################################################

broad_stage_tissue_counts <- meta %>%
  distinct(Patient, Stage_broad, Tissue) %>%
  count(Stage_broad, Tissue, name = "n_patient_tissue_samples")

write.csv(
  broad_stage_tissue_counts,
  file = file.path(
    proj_dir,
    "broad_stage_tissue_patient_sample_counts.csv"
  ),
  row.names = FALSE
)

p3 <- ggplot(
  broad_stage_tissue_counts,
  aes(
    x = Stage_broad,
    y = n_patient_tissue_samples,
    fill = Tissue
  )
) +
  geom_col(position = "stack") +
  labs(
    title = "Tissue origins represented within each broad stage",
    x = "Broad cancer stage",
    y = "Number of patient-tissue samples",
    fill = "Tissue"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(
    proj_dir,
    "plot_tissue_origin_by_broad_stage.png"
  ),
  plot = p3,
  width = 8,
  height = 5,
  dpi = 300
)

cat("broad_stage_tissue_patient_sample_counts.csv\n")
cat("plot_tissue_origin_by_broad_stage.png\n")

############################################################
# 13. Quick sanity checks
############################################################

cat("\nSanity check: patient-level matrix dimensions:\n")
print(dim(patient_celltype_matrix))

cat("\nSanity check: patient-tissue matrix dimensions:\n")
print(dim(patient_tissue_celltype_matrix))

cat("\nPatient counts by broad stage from matrix:\n")
print(table(patient_celltype_matrix$Stage_broad, useNA = "ifany"))

cat("\nDone.\n")
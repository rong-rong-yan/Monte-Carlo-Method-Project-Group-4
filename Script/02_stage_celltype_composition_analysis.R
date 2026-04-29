############################################################
# Stage-level and patient-level cell-type composition analysis
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

meta_file <- file.path(proj_dir, "GSE131907_cell_annotation_with_clinical_metadata.csv")

meta <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)

cat("Metadata dimensions:\n")
print(dim(meta))

cat("\nColumn names:\n")
print(colnames(meta))

############################################################
# 2. Set detailed stage order
############################################################

stage_levels <- c("IA", "IA2", "IA3", "IB", "IIA", "IIB", "IIIA", "IV")

meta$Stage <- factor(meta$Stage, levels = stage_levels)

cat("\nCell counts by Stage:\n")
print(table(meta$Stage, useNA = "ifany"))

cat("\nPatient counts by Stage:\n")
patient_stage <- meta %>%
  distinct(Patient, Stage) %>%
  count(Stage, name = "n_patients")

print(patient_stage)

############################################################
# 3. Identify cell-type column
############################################################
# Check the printed column names first.
# Common possibilities in this dataset may include:
# "Cell_type", "CellType", "cell_type", "Major_cell_type", etc.

possible_celltype_cols <- c(
  "Cell_type",
  "CellType",
  "cell_type",
  "Major_cell_type",
  "MajorCellType",
  "Annotation",
  "annotation",
  "Cell_subtype",
  "Cell_subtype_1"
)

celltype_col <- possible_celltype_cols[possible_celltype_cols %in% colnames(meta)][1]

if (is.na(celltype_col)) {
  stop("Could not automatically find cell-type column. Please check colnames(meta) and set celltype_col manually.")
}

cat("\nUsing cell-type column:\n")
print(celltype_col)

############################################################
# 4. Cell counts by stage and cell type
############################################################

stage_celltype_counts <- meta %>%
  count(Stage, .data[[celltype_col]], name = "n_cells") %>%
  group_by(Stage) %>%
  mutate(
    stage_total_cells = sum(n_cells),
    proportion_within_stage = n_cells / stage_total_cells
  ) %>%
  ungroup()

write.csv(
  stage_celltype_counts,
  file = file.path(proj_dir, "stage_celltype_counts_and_proportions.csv"),
  row.names = FALSE
)

cat("\nSaved: stage_celltype_counts_and_proportions.csv\n")

############################################################
# 5. Patient-level cell-type proportions
############################################################
# This is more important than raw cell-level counts because patients
# contribute different numbers of cells.

patient_celltype_counts <- meta %>%
  count(Patient, Stage, Tissue, .data[[celltype_col]], name = "n_cells") %>%
  group_by(Patient) %>%
  mutate(
    patient_total_cells = sum(n_cells),
    proportion_within_patient = n_cells / patient_total_cells
  ) %>%
  ungroup()

write.csv(
  patient_celltype_counts,
  file = file.path(proj_dir, "patient_celltype_counts_and_proportions.csv"),
  row.names = FALSE
)

cat("Saved: patient_celltype_counts_and_proportions.csv\n")

############################################################
# 6. Tissue-aware patient-level cell-type proportions
############################################################
# This is useful because Stage and Tissue are strongly related in this dataset.

patient_tissue_celltype_counts <- meta %>%
  count(Patient, Stage, Tissue, .data[[celltype_col]], name = "n_cells") %>%
  group_by(Patient, Tissue) %>%
  mutate(
    patient_tissue_total_cells = sum(n_cells),
    proportion_within_patient_tissue = n_cells / patient_tissue_total_cells
  ) %>%
  ungroup()

write.csv(
  patient_tissue_celltype_counts,
  file = file.path(proj_dir, "patient_tissue_celltype_counts_and_proportions.csv"),
  row.names = FALSE
)

cat("Saved: patient_tissue_celltype_counts_and_proportions.csv\n")

############################################################
# 7. Plot: cell-type composition by detailed stage
############################################################

p1 <- ggplot(
  stage_celltype_counts,
  aes(x = Stage, y = proportion_within_stage, fill = .data[[celltype_col]])
) +
  geom_col(position = "stack") +
  labs(
    title = "Cell-type composition by detailed cancer stage",
    x = "Cancer stage",
    y = "Proportion of cells within stage",
    fill = "Cell type"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(proj_dir, "plot_celltype_composition_by_stage.png"),
  plot = p1,
  width = 9,
  height = 6,
  dpi = 300
)

cat("Saved: plot_celltype_composition_by_stage.png\n")

############################################################
# 8. Plot: patient-level cell-type proportions
############################################################

p2 <- ggplot(
  patient_celltype_counts,
  aes(x = Patient, y = proportion_within_patient, fill = .data[[celltype_col]])
) +
  geom_col(position = "stack") +
  facet_grid(. ~ Stage, scales = "free_x", space = "free_x") +
  labs(
    title = "Patient-level cell-type composition, grouped by detailed stage",
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
  filename = file.path(proj_dir, "plot_patient_celltype_composition_by_stage.png"),
  plot = p2,
  width = 14,
  height = 6,
  dpi = 300
)

cat("Saved: plot_patient_celltype_composition_by_stage.png\n")

############################################################
# 9. Plot: tissue origin by stage
############################################################

stage_tissue_counts <- meta %>%
  distinct(Patient, Stage, Tissue) %>%
  count(Stage, Tissue, name = "n_patient_tissue_samples")

write.csv(
  stage_tissue_counts,
  file = file.path(proj_dir, "stage_tissue_patient_sample_counts.csv"),
  row.names = FALSE
)

p3 <- ggplot(
  stage_tissue_counts,
  aes(x = Stage, y = n_patient_tissue_samples, fill = Tissue)
) +
  geom_col(position = "stack") +
  labs(
    title = "Tissue origins represented within each detailed stage",
    x = "Cancer stage",
    y = "Number of patient-tissue samples",
    fill = "Tissue"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(proj_dir, "plot_tissue_origin_by_stage.png"),
  plot = p3,
  width = 8,
  height = 5,
  dpi = 300
)

cat("Saved: plot_tissue_origin_by_stage.png\n")

############################################################
# 10. Optional: create broad stage label too, but keep detailed stage
############################################################

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

meta$Stage_broad <- factor(meta$Stage_broad, levels = c("I", "II", "III", "IV"))

write.csv(
  meta,
  file = file.path(proj_dir, "GSE131907_cell_annotation_with_clinical_metadata_stage_broad.csv"),
  row.names = FALSE
)

cat("Saved: GSE131907_cell_annotation_with_clinical_metadata_stage_broad.csv\n")

############################################################
# 11. Done
############################################################

cat("\nDone.\n")
# data-raw/make-narps-dcm.R
#
# Build the narps_dcm dataset: subject-level DCM posterior summaries for the
# 48-subject NARPS group analysis, plus participant covariates. This is the
# real-data input for the robust + sparse group DCM (rsdcm).
#
# Source files (from the thesis Results/Presented folder, not shipped):
#   res2.rds         - list with $eta (48 x 22), $Cp_list (48 x [22 x 22]),
#                      $parameter_names, $subject_names
#   participants.tsv - participant_id, group, gender, age
#
# The full res2.rds also carries PEB_list / Ce_list, which are not needed here
# and are dropped to keep the shipped object small (~0.3 MB xz-compressed).
#
# Usage (from the package root, adjusting SRC if needed):
#   Rscript data-raw/make-narps-dcm.R

SRC <- Sys.getenv("NARPS_SRC",
                  "../../Results/Presented")   # relative to package root

res2  <- readRDS(file.path(SRC, "res2.rds"))
parts <- read.delim(file.path(SRC, "participants.tsv"),
                    stringsAsFactors = FALSE)

eta <- res2$eta
Cp  <- res2$Cp_list
pnames <- res2$parameter_names
N <- nrow(eta); p <- ncol(eta)
stopifnot(length(Cp) == N, ncol(Cp[[1]]) == p)

# Align covariates: subject id like "PEB_sub-001_runwise" -> "sub-001"
sid  <- sub(".*(sub-[0-9]+).*", "\\1", res2$subject_names)
idx  <- match(sid, parts$participant_id)
if (any(is.na(idx)))
  stop("Unmatched subjects: ", paste(sid[is.na(idx)], collapse = ", "))
meta <- parts[idx, ]

narps_dcm <- list(
  eta             = eta,                       # 48 x 22 posterior means
  Cp              = Cp,                         # list of 48, each 22 x 22
  parameter_names = pnames,                     # 16 A(i,j) + 6 B(i,j,input)
  subject         = sid,                        # de-identified NARPS ids
  covariates      = data.frame(group  = meta$group,
                               gender = meta$gender,
                               age    = meta$age,
                               stringsAsFactors = FALSE),
  regions         = c("vmPFC", "vStr", "amyg", "aIns")
)

if (!dir.exists("data")) dir.create("data")
save(narps_dcm, file = file.path("data", "narps_dcm.rda"),
     compress = "xz", compression_level = 9)

fsize <- file.info(file.path("data", "narps_dcm.rda"))$size
message(sprintf("Wrote data/narps_dcm.rda (%.1f KB): %d subjects, %d parameters.",
                fsize / 1024, N, p))

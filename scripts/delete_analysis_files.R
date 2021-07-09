## DELETE ALL FILES IN THESE FOLDERS ## 

# stage <- "incubation"
stage <- "chick_rearing"

## FILTER STEPS ##
#1#
# Speed filtered
do.call(file.remove, list(list.files(paste0("data/analysis/speed_filtered/", stage), full.names=T)))
# trip split
do.call(file.remove, list(list.files(paste0("data/analysis/trip_split/", stage), full.names=T)))
# trip summaries
do.call(file.remove, list(list.files(paste0("data/analysis/trip_summary/", stage), full.names=T)))
# interpolated
do.call(file.remove, list(list.files(paste0("data/analysis/interpolated/", stage), full.names=T)))
# individual UDs
do.call(file.remove, list(list.files(paste0("data/analysis/ind_UDs/", stage), full.names=T)))
# inter-annual UDs
do.call(file.remove, list(list.files(paste0("data/analysis/interannual_UDs/", stage), full.names=T)))
# yearly UDs
do.call(file.remove, list(list.files(paste0("data/analysis/yearly_UDs/", stage), full.names=T)))
# inter-annual % home ranges
do.call(file.remove, list(list.files(paste0("data/analysis/interannaul_HRs/", stage, "/polygon/50"), full.names=T)))
do.call(file.remove, list(list.files(paste0("data/analysis/interannaul_HRs/", stage, "/polygon/95"), full.names=T)))
# yearly % home ranges
do.call(file.remove, list(list.files(paste0("data/analysis/yearly_HRs/", stage, "/polygon/50"), full.names=T)))
do.call(file.remove, list(list.files(paste0("data/analysis/yearly_HRs/", stage, "/polygon/95"), full.names=T)))

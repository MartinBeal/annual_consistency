## Delete analysis files ##

## Choose breeding stage # ----------------------------------------------------
# stage  <- "incubation"
stage  <- "chick_rearing"

## Choose h parameter method # ------------------------------------------------
htype <- "mag"
htype <- "href1"
htype <- "href2"


#------------------------------------------------------------------------------
## individual UDs ## ----------------------------------------------------------
folder <- "data/analysis/ind_UDs/"
files <- list.files(paste0(folder, stage), pattern = fixed(htype), full.names = T)
do.call(file.remove, list(files, full.names=T))
do.call(file.remove, list(files[!str_detect(files, pattern = "href2")], full.names=T))

## interannual UDs ## ----------------------------------------------------------
folder <- "data/analysis/interannual_UDs/"
files <- list.files(paste0(folder, stage), pattern = fixed(htype), full.names = T)
do.call(file.remove, list(files, full.names=T))
# do.call(file.remove, list(files[!str_detect(files, pattern = "href2")], full.names=T))

## interannual HRs ## ----------------------------------------------------------
folder <- "data/analysis/interannaul_HRs/"
folders1 <- list.files(paste0(folder, stage), full.names = T)

# files <- list.files(paste0(folder, stage), pattern = fixed(htype), full.names = T)
# do.call(file.remove, list(files, full.names=T)))

## yearly UDs ## ----------------------------------------------------------
folder  <- "data/analysis/yearly_UDs/"
folders <- list.files(paste0(folder, stage), full.names = T)
do.call(file.remove, 
        lapply(folders, function(x) list.files(x, full.names = T, pattern=fixed(htype)))
)

# do.call(file.remove,
#         lapply(folders, function(x) {
#   xx <- list.files(x, full.names = T)
#   xxx <- xx[!str_detect(xx, pattern = "href2")]
#   })
# )


## yearly HRs ## ----------------------------------------------------------
folder <- "data/analysis/yearly_HRs/"
folders1 <- list.files(paste0(folder, stage), full.names = T)
files <- list.files(folders1, full.names = T)

# do.call(file.remove, 
#         lapply(folders, function(x) list.files(x, full.names = T, pattern=fixed(htype)))
# )


### tracking data (e.g. interpolated data) ##

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

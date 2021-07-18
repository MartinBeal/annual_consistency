### Generate arbitrary-start DateTime sequence that has exact same timestep 
## durations and total extent as original data, but ignores individuals and 
# separates them by 15 min interval 

## get time step interval btwn points per ind ##
ts <- lapply(split(one, one$ID), function(x){
  timestep <- c(as.vector(difftime(x$DateTime[2:(length(x$DateTime))], 
                       x$DateTime[1:(nrow(x)-1)], units="secs")), NA)
})

one$timestep <- unlist(ts)

## make time steps btwn individuals an even 15 min ##
one$timestep <- ifelse(is.na(one$timestep), 15*60, one$timestep)
## last point in dataset should be NA ##
one$timestep[nrow(one)] <- NA

## generate artificial time series across IDs based on (mostly) real time steps ##
# make the desired time sequence
dati_start <- as.POSIXct("1970-01-01 00:00:01", "UTC")
dati_end   <- dati_start + sum(na.omit(timestep))*60

# datiseq <- seq(from=dati_start, to=dati_end, by=1)

dts <- vector(mode = "list", length=nrow(one))
dts[[1]] <- dati_start

for(i in 2:nrow(one)){
  
  ts <- one$timestep[i-1]
  
  dts[[i]] <- dts[[i-1]] + ts
  
}

dts_all <- do.call(rbind, dts)
attributes(dts_all) <- attributes(dts[[1]]) ## add POSIXct info back

one$genDT <- dts_all

## check that first and last date times make sense ## 
one %>% group_by(ID) %>% summarise(first(DateTime), last(DateTime), first(genDT), last(genDT))

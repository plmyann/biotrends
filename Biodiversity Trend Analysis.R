setwd('D:\\Rspace\\')

library(sp)
library(dplyr)
library(rgdal)

#########load the BBS routes of polyline and create a connected column for further combination with other information

rtlocation <- readOGR('/BBSroutes.shp')##load BBS routes of United States

rtlocation@data$ID <- 1:length(rtlocation)

###notes: the rtlocation included a column of RTENO as the ID, which is the combination of state ID and Route ID
###e.g., 14001 represents the data coming from Route 001 of State 14. But the info in the BBS observations was stored in separate
###columns (i.e., StateNum, Route)

############Create columns of CountryNum, StateNum and Route for further combination with BBS observations

rtlocation@data$RTENO <- as.numeric(levels(rtlocation@data$RTENO))[rtlocation@data$RTENO]

rtlocation@data$CountryNum <- 840##840 is the Code for United States

rtlocation@data$StateNum <- gsub('\\d{3}$','',rtlocation@data$RTENO)##the last three charaters as the StateNum

rtlocation@data$Route <- as.numeric(stringr::str_sub(rtlocation@data$RTENO,start = -3))##exclue the last three letters

rtlocation@data$RTENO <- paste(840, rtlocation@data$RTENO, sep = '')##Add 840 at the beginning as the country ID

############## extract the land cover data to routes

lcfolder <- '/MCD12Q1'

lcfile <- list.files(lcfolder, pattern = '.tif$', full.names = T)

lclist2 <- lapply(lcfile,function(lc){ ## loop for each year 2001-2018
  
  lcrst <- raster::raster(lc)
  
  lcdf <- raster::extract(lcrst, rtlocation, df=T, buffer = 5000) ##extract the land cover around the routes with 5km buffer 
  
  lcdf <- data.frame(lcdf)
  
  colnames(lcdf)[2] <- 'Value'
  
  lcdf1 <- lcdf %>% group_by(ID, Value) %>% summarise(Freq = n()) %>%
    
    mutate(proc = (Freq/sum(Freq) * 100)) ##calculate the proportion of each land cover type
  
  lc <- sub("/MCD12Q1/",
                   replacement = "", lc)
  lcdf1$Year <- sub(".tif$", replacement = "", lc)
  
  cat(lc, '\n')
  
  return(data.frame(lcdf1))
  
})

lcdf.bbs_modis <- bind_rows(lclist2)

lcdf.bbs_modis$X <- NULL

stpyp <- reshape2::dcast(lcdf.bbs_modis,  ID + Year ~ Value, fill = 0, value.var = 'proc' ) ##redistribute the data by proportion

### weather2019.csv stored all information of each survey of along each route at all years.
weather <- read.csv('/weather2019.csv')

###create a RTENO column for weather to join with rtlocaiton

weather$RTENO <- paste(weather$CountryNum, weather$StateNum, sprintf("%03d", weather$Route), sep = '')

weather <- weather[weather$Year > 2000 & weather$Year < 2019,]##remain the data from 2001-2018

###stpyp only include route ID and the habitat composition, which need a combination with rtlocation to get the RTENO

route.lc.proc <- merge(stpyp, rtlocation, by = 'ID')

###merge the route habitat composition with survey info and only include the columns of RTENO, Year, proportion and RPID
###RPID represent the survey protocol with 101 as a standard protocol, 501 as a modified protocol and NA as no survey.

route.lc1.proc <- merge(route.lc.proc, weather, by=c('RTENO','Year'), all.y = F, all.x = T)[,c(1:2,4:19,41)]

route.lc1.proc <- data.frame(route.lc1.proc)

route.lc1.proc$Year <- as.integer(route.lc1.proc$Year)

#####################summarise the habitat change

####by threshold of the habitat composition

route.lcrt.proc <- route.lc1.proc

route.lcrt.proc <- route.lcrt.proc[with(route.lcrt.proc, order(RTENO,Year)),]

lcnames <- colnames(route.lcrt.proc)[3:18]##get the colnames to calculate the difference of habitat composition 

route.lcrt.proc1 <- route.lcrt.proc %>% group_by(RTENO) %>% ###calculate the difference of habitat in 2 years.
  
  mutate_at(lcnames, list(~c(0, diff(.))))

route.lcrt.proc1 <- route.lcrt.proc1[!duplicated(route.lcrt.proc1[c(1:2)]),]###remove the duplicates with same year and RTENO

####creat a new column to indicate any of the habitat type changed more than 5%, 10%...

##set 1 to the new column of habitat change index when the absolute value of columns beginning 
##with 'X' (column of the habitat proportion change) is over 3

### 5% change

d.tmp <- route.lcrt.proc1 %>% filter_at(vars(starts_with('X')), any_vars(abs(.) >= 5)) %>% mutate(HCIndex5 = 1) 
## merge the filtered data.frame with the previous one
route.lcrt.proc1 <- merge(route.lcrt.proc1, d.tmp, all.x =T)
##set 0 to NAs in the habitat change index column
route.lcrt.proc1$HCIndex5[is.na(route.lcrt.proc1$HCIndex5)] <- 0

### 10% change
d.tmp <- route.lcrt.proc1 %>% filter_at(vars(starts_with('X')), any_vars(abs(.) >= 10)) %>% mutate(HCIndex10 = 1) 
route.lcrt.proc1 <- merge(route.lcrt.proc1, d.tmp, all.x =T)
route.lcrt.proc1$HCIndex10[is.na(route.lcrt.proc1$HCIndex10)] <- 0

### 15% change
d.tmp <- route.lcrt.proc1 %>% filter_at(vars(starts_with('X')), any_vars(abs(.) >= 15)) %>% mutate(HCIndex15 = 1) 
route.lcrt.proc1 <- merge(route.lcrt.proc1, d.tmp, all.x =T)
route.lcrt.proc1$HCIndex15[is.na(route.lcrt.proc1$HCIndex15)] <- 0

### 20% change
d.tmp <- route.lcrt.proc1 %>% filter_at(vars(starts_with('X')), any_vars(abs(.) >= 20)) %>% mutate(HCIndex20 = 1) 
route.lcrt.proc1 <- merge(route.lcrt.proc1, d.tmp, all.x =T)
route.lcrt.proc1$HCIndex20[is.na(route.lcrt.proc1$HCIndex20)] <- 0

### 25% change
d.tmp <- route.lcrt.proc1 %>% filter_at(vars(starts_with('X')), any_vars(abs(.) >= 25)) %>% mutate(HCIndex25 = 1) 
route.lcrt.proc1 <- merge(route.lcrt.proc1, d.tmp, all.x =T)
route.lcrt.proc1$HCIndex25[is.na(route.lcrt.proc1$HCIndex25)] <- 0

### 30% change
d.tmp <- route.lcrt.proc1 %>% filter_at(vars(starts_with('X')), any_vars(abs(.) >= 30)) %>% mutate(HCIndex30 = 1) 
route.lcrt.proc1 <- merge(route.lcrt.proc1, d.tmp, all.x =T)
route.lcrt.proc1$HCIndex30[is.na(route.lcrt.proc1$HCIndex30)] <- 0

### remove the survey routes with whole NAs in 18 years
route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RTENO2 = case_when(length(RPID[is.na(RPID)])==18 ~ 0, TRUE ~ 1)) 

route.lcrt.proc1 <- route.lcrt.proc1[route.lcrt.proc1$RTENO2==1,]

route.lcrt.proc1$RTENO2 <- NULL

###build a funciton to set the index for survey stop

foo <- function(x, test) {
  
  xlist <- lapply(1:length(which(test)), function(i){
    ###the loop accouts for multiple rows of an individual route that match the test
  ind_v <- which(test)[i]###the ordinal number of row that match the test firstly.
  x[seq_along(x) <= ind_v] <- NA ##Set NAs to all the row before the first matched row in the column of stop index 
  ind_non_na = which(!is.na(x))[1] ##the ordinal number of the row that match the test secondly
  x[seq_along(x) >= ind_non_na] <- NA ## set NAs to all the row after the second matched row
  x[ind_v] = max(-1, ind_non_na, na.rm = TRUE) ##set the larger value between -1 and ordinal number of the second matched row
                                               ####to the first matched row
  return(x)
  
  })
  
  xlist <- bind_cols(xlist) ##column bind for all the rows match the test and sum all the columns
  
  xlist <- rowSums(xlist,na.rm = T)
  
  ###so there will be three different types indicating the stop index: NA for no change, 0 for stop until the end of
  ###the survey, and a ordinal value showing when matching the test again.
  
  xlist[xlist==0] <- NA
  xlist[xlist==-1] <- 0
  
  return(xlist)
  
}

tests <- c("((is.na(RPID) & !is.na(lag(RPID,default = RPID[1])))|(!is.na(RPID) & is.na(lead(RPID,default = RPID[18]))))",
           "(is.na(RPID) & !is.na(lag(RPID,default = RPID[1])))", "(!is.na(RPID) & is.na(lead(RPID,default = RPID[18])))")
 
route.lcrt.list <- list() 
rndlist <- list()
bbslist <- list()

for (t in 1:3) {
  

tt <- tests[t]  

if(t==3){
  
##To test a habitat change happens associated with a survey stop
route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndex5 = foo(RPID, test = HCIndex5!=0 & lead(HCIndex5==0,default = HCIndex5[18]) & eval(parse(text = tt)))) %>%

mutate(RSIndex5 = case_when(RSIndex5==0 ~ 1, TRUE ~0)) ##t+1 

##hcindex10
route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndex10 = foo(RPID, test = HCIndex10!=0 & lead(HCIndex10==0,default = HCIndex10[18]) & eval(parse(text = tt)))) %>%
  
  mutate(RSIndex10 = case_when(RSIndex10==0 ~ 1, TRUE ~0))

##hcindex15
route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndex15 = foo(RPID, test = HCIndex15!=0 & lead(HCIndex15==0,default = HCIndex15[18]) & eval(parse(text = tt)))) %>%
  
  mutate(RSIndex15 = case_when(RSIndex15==0 ~ 1, TRUE ~0))
##hcindex20
route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndex20 = foo(RPID, test = HCIndex20!=0 & lead(HCIndex20==0,default = HCIndex5[18]) & eval(parse(text = tt)))) %>%
  
  mutate(RSIndex20 = case_when(RSIndex20==0 ~ 1, TRUE ~0))
##hcindex25
route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndex25 = foo(RPID, test = HCIndex25!=0 & lead(HCIndex25==0,default = HCIndex25[18]) & eval(parse(text = tt)))) %>%
  
  mutate(RSIndex25 = case_when(RSIndex25==0 ~ 1, TRUE ~0))
##hcindex30
route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndex30 = foo(RPID, test = HCIndex30!=0 & lead(HCIndex30==0,default = HCIndex30[18]) & eval(parse(text = tt)))) %>%
  
  mutate(RSIndex30 = case_when(RSIndex30==0 ~ 1, TRUE ~0))
} else {
  
  ##hcindex 5
  route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndex5 = foo(RPID, test = HCIndex5!=0 & eval(parse(text = tt)))) %>%
    
    mutate(RSIndex5 = case_when(RSIndex5==0 ~ 1, TRUE ~0))##t and both
  
  ##hcindex10
  route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndex10 = foo(RPID, test = HCIndex10!=0 & eval(parse(text = tt)))) %>%
    
    mutate(RSIndex10 = case_when(RSIndex10==0 ~ 1, TRUE ~0))
  
  ##hcindex15
  route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndex15 = foo(RPID, test = HCIndex15!=0 & eval(parse(text = tt)))) %>%
    
    mutate(RSIndex15 = case_when(RSIndex15==0 ~ 1, TRUE ~0))
  ##hcindex20
  route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndex20 = foo(RPID, test = HCIndex20!=0 & eval(parse(text = tt)))) %>%
    
    mutate(RSIndex20 = case_when(RSIndex20==0 ~ 1, TRUE ~0))
  ##hcindex25
  route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndex25 = foo(RPID, test = HCIndex25!=0 & eval(parse(text = tt)))) %>%
    
    mutate(RSIndex25 = case_when(RSIndex25==0 ~ 1, TRUE ~0))
  ##hcindex30
  route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndex30 = foo(RPID, test = HCIndex30!=0 & eval(parse(text = tt)))) %>%
    
    mutate(RSIndex30 = case_when(RSIndex30==0 ~ 1, TRUE ~0))
  
}

route.lcrt.list[[t]] <- route.lcrt.proc1

hcroute.r5 <- NA
hcroute.r10 <- NA
hcroute.r15 <- NA
hcroute.r20 <- NA
hcroute.r25 <- NA
hcroute.r30 <- NA
stroute.r5 <- NA
stroute.r10 <- NA
stroute.r15 <- NA
stroute.r20 <- NA
stroute.r25 <- NA
stroute.r30 <- NA


for(i in 1:1000){###calculate the random for n=1000 times
  
  if(t==3){
    
    set.seed(i)
    
    ### random distribute the habitat change index of 5% 
    ###random 5
    route.lcrt.proc1$HCIndexr5dm <- sample(route.lcrt.proc1$HCIndex5, length(route.lcrt.proc1$HCIndex5))
    
    route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndexr5dm = foo(RPID, test = HCIndexr5dm!=0 & lead(HCIndexr5dm==0,default = HCIndexr5dm[18]) & eval(parse(text = tt)))) %>%
      
      mutate(RSIndexr5dm = case_when(RSIndexr5dm==0 ~ 1, TRUE ~0))
    
    ###random 10
    route.lcrt.proc1$HCIndexr10dm <- sample(route.lcrt.proc1$HCIndex10, length(route.lcrt.proc1$HCIndex10))
    
    route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndexr10dm = foo(RPID, test = HCIndexr10dm!=0 & lead(HCIndexr10dm==0,default = HCIndexr10dm[18]) & eval(parse(text = tt)))) %>%
      
      mutate(RSIndexr10dm = case_when(RSIndexr10dm==0 ~ 1, TRUE ~0))
    
    ###random 15
    route.lcrt.proc1$HCIndexr15dm <- sample(route.lcrt.proc1$HCIndex15, length(route.lcrt.proc1$HCIndex15))
    
    route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndexr15dm = foo(RPID, test = HCIndexr15dm!=0 & lead(HCIndexr15dm==0,default = HCIndexr15dm[18]) & eval(parse(text = tt)))) %>%
      
      mutate(RSIndexr15dm = case_when(RSIndexr15dm==0 ~ 1, TRUE ~0))
    
    ###random 20
    route.lcrt.proc1$HCIndexr20dm <- sample(route.lcrt.proc1$HCIndex20, length(route.lcrt.proc1$HCIndex20))
    
    route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndexr20dm = foo(RPID, test = HCIndexr20dm!=0 & lead(HCIndexr20dm==0,default = HCIndexr20dm[18]) & eval(parse(text = tt)))) %>%
      
      mutate(RSIndexr20dm = case_when(RSIndexr20dm==0 ~ 1, TRUE ~0))
    
    ###random 25
    route.lcrt.proc1$HCIndexr25dm <- sample(route.lcrt.proc1$HCIndex25, length(route.lcrt.proc1$HCIndex25))
    
    route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndexr25dm = foo(RPID, test = HCIndexr25dm!=0 & lead(HCIndexr25dm==0,default = HCIndexr25dm[18]) & eval(parse(text = tt)))) %>%
      
      mutate(RSIndexr25dm = case_when(RSIndexr25dm==0 ~ 1, TRUE ~0))
    
    ###random 30
    route.lcrt.proc1$HCIndexr30dm <- sample(route.lcrt.proc1$HCIndex30, length(route.lcrt.proc1$HCIndex30))
    
    route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndexr30dm = foo(RPID, test = HCIndexr30dm!=0 & lead(HCIndexr30dm==0,default = HCIndexr30dm[18]) & eval(parse(text = tt)))) %>%
      
      mutate(RSIndexr30dm = case_when(RSIndexr30dm==0 ~ 1, TRUE ~0))
  
  } else {
 
    
    set.seed(i)
    
    ### random distribute the habitat change index of 5% 
    ###random 5
    route.lcrt.proc1$HCIndexr5dm <- sample(route.lcrt.proc1$HCIndex5, length(route.lcrt.proc1$HCIndex5))
    
    route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndexr5dm = foo(RPID, test = HCIndexr5dm!=0 & eval(parse(text = tt)))) %>%
      
      mutate(RSIndexr5dm = case_when(RSIndexr5dm==0 ~ 1, TRUE ~0))
    
    ###random 10
    route.lcrt.proc1$HCIndexr10dm <- sample(route.lcrt.proc1$HCIndex10, length(route.lcrt.proc1$HCIndex10))
    
    route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndexr10dm = foo(RPID, test = HCIndexr10dm!=0 & eval(parse(text = tt)))) %>%
      
      mutate(RSIndexr10dm = case_when(RSIndexr10dm==0 ~ 1, TRUE ~0))
    
    ###random 15
    route.lcrt.proc1$HCIndexr15dm <- sample(route.lcrt.proc1$HCIndex15, length(route.lcrt.proc1$HCIndex15))
    
    route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndexr15dm = foo(RPID, test = HCIndexr15dm!=0 & eval(parse(text = tt)))) %>%
      
      mutate(RSIndexr15dm = case_when(RSIndexr15dm==0 ~ 1, TRUE ~0))
    
    ###random 20
    route.lcrt.proc1$HCIndexr20dm <- sample(route.lcrt.proc1$HCIndex20, length(route.lcrt.proc1$HCIndex20))
    
    route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndexr20dm = foo(RPID, test = HCIndexr20dm!=0 & eval(parse(text = tt)))) %>%
      
      mutate(RSIndexr20dm = case_when(RSIndexr20dm==0 ~ 1, TRUE ~0))
    
    ###random 25
    route.lcrt.proc1$HCIndexr25dm <- sample(route.lcrt.proc1$HCIndex25, length(route.lcrt.proc1$HCIndex25))
    
    route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndexr25dm = foo(RPID, test = HCIndexr25dm!=0 & eval(parse(text = tt)))) %>%
      
      mutate(RSIndexr25dm = case_when(RSIndexr25dm==0 ~ 1, TRUE ~0))
    
    ###random 30
    route.lcrt.proc1$HCIndexr30dm <- sample(route.lcrt.proc1$HCIndex30, length(route.lcrt.proc1$HCIndex30))
    
    route.lcrt.proc1 <- route.lcrt.proc1 %>% group_by(RTENO) %>% mutate(RSIndexr30dm = foo(RPID, test = HCIndexr30dm!=0 & eval(parse(text = tt)))) %>%
      
      mutate(RSIndexr30dm = case_when(RSIndexr30dm==0 ~ 1, TRUE ~0)) 
    
}

HC.sum <- route.lcrt.proc1 %>% group_by(RTENO) %>% 
  
  summarise(rtr5stops = length(RSIndexr5dm[RSIndexr5dm!=0]),
            rtr10stops = length(RSIndexr10dm[RSIndexr10dm!=0]),
            rtr15stops = length(RSIndexr15dm[RSIndexr15dm!=0]),
            rtr20stops = length(RSIndexr20dm[RSIndexr20dm!=0]),
            rtr25stops = length(RSIndexr25dm[RSIndexr25dm!=0]),
            rtr30stops = length(RSIndexr30dm[RSIndexr30dm!=0]),
            hcr5times = length(HCIndexr5dm[HCIndexr5dm!=0]),
            hcr10times = length(HCIndexr10dm[HCIndexr10dm!=0]),
            hcr15times = length(HCIndexr15dm[HCIndexr15dm!=0]),
            hcr20times = length(HCIndexr20dm[HCIndexr20dm!=0]),
            hcr25times = length(HCIndexr25dm[HCIndexr25dm!=0]),
            hcr30times = length(HCIndexr30dm[HCIndexr30dm!=0]))


hcroute.r5[i] <- length(HC.sum$hcr5times[HC.sum$hcr5times!=0])
stroute.r5[i] <- length(HC.sum$rtr5stops[HC.sum$rtr5stops!=0])

hcroute.r10[i] <- length(HC.sum$hcr10times[HC.sum$hcr10times!=0])
stroute.r10[i] <- length(HC.sum$rtr10stops[HC.sum$rtr10stops!=0])

hcroute.r15[i] <- length(HC.sum$hcr15times[HC.sum$hcr15times!=0])
stroute.r15[i] <- length(HC.sum$rtr15stops[HC.sum$rtr15stops!=0])

hcroute.r20[i] <- length(HC.sum$hcr20times[HC.sum$hcr20times!=0])
stroute.r20[i] <- length(HC.sum$rtr20stops[HC.sum$rtr20stops!=0])

hcroute.r25[i] <- length(HC.sum$hcr25times[HC.sum$hcr25times!=0])
stroute.r25[i] <- length(HC.sum$rtr25stops[HC.sum$rtr25stops!=0])

hcroute.r30[i] <- length(HC.sum$hcr30times[HC.sum$hcr30times!=0])
stroute.r30[i] <- length(HC.sum$rtr30stops[HC.sum$rtr30stops!=0])

cat(i, '\n')

}

rndrouts <- data.frame(hcroute.r5, hcroute.r10, hcroute.r15, hcroute.r20, hcroute.r25, hcroute.r30, 
                       stroute.r5, stroute.r10, stroute.r15, stroute.r20, stroute.r25, stroute.r30)

rndrouts$changepro.r5 <- rndrouts$stroute.r5/rndrouts$hcroute.r5*100
rndrouts$changepro.r10 <- rndrouts$stroute.r10/rndrouts$hcroute.r10*100
rndrouts$changepro.r15 <- rndrouts$stroute.r15/rndrouts$hcroute.r15*100
rndrouts$changepro.r20 <- rndrouts$stroute.r20/rndrouts$hcroute.r20*100
rndrouts$changepro.r25 <- rndrouts$stroute.r25/rndrouts$hcroute.r25*100
rndrouts$changepro.r30 <- rndrouts$stroute.r30/rndrouts$hcroute.r30*100

rndlist[[t]] <- rndrouts

###routes that stopped at least once

consecutive <- function(t){
  
  t <- sort(as.numeric(unique(t)))
  
  b <- cumsum(c(1, diff(t)>1))
  
  c <- rle(b)
  
  return(c)
}

route.lcrt.year.na <-  route.lcrt.proc1
route.lcrt.year.na$Year[is.na(route.lcrt.year.na$RPID)] <- NA

rstop.sum <- route.lcrt.year.na %>% group_by(RTENO) %>%
  
  summarise(rtstops = length(consecutive(Year)$lengths)-1)

rtstops <- length(rstop.sum$rtstops[rstop.sum$rtstops!=0])

####generate the whole output.

HC.sum <- route.lcrt.proc1 %>% group_by(RTENO) %>% 
  
  summarise(hc5times = length(HCIndex5[HCIndex5!=0]),
            hc10times = length(HCIndex10[HCIndex10!=0]),
            hc15times = length(HCIndex15[HCIndex15!=0]),
            hc20times = length(HCIndex20[HCIndex20!=0]),
            hc25times = length(HCIndex25[HCIndex25!=0]),
            hc30times = length(HCIndex30[HCIndex30!=0]),
            
            rt5stops = length(RSIndex5[RSIndex5!=0]),
            rt10stops = length(RSIndex10[RSIndex10!=0]),
            rt15stops = length(RSIndex15[RSIndex15!=0]),
            rt20stops = length(RSIndex20[RSIndex20!=0]),
            rt25stops = length(RSIndex25[RSIndex25!=0]),
            rt30stops = length(RSIndex30[RSIndex30!=0]))

output <- data.frame(rstops = rtstops,
                     hcroute = c(length(HC.sum$hc5times[HC.sum$hc5times!=0]),length(HC.sum$hc10times[HC.sum$hc10times!=0]),
                                 length(HC.sum$hc15times[HC.sum$hc15times!=0]), length(HC.sum$hc20times[HC.sum$hc20times!=0]),
                                 length(HC.sum$hc25times[HC.sum$hc25times!=0]),length(HC.sum$hc30times[HC.sum$hc30times!=0])),
                     
                     stroute = c(length(HC.sum$rt5stops[HC.sum$rt5stops!=0]),length(HC.sum$rt10stops[HC.sum$rt10stops!=0]),
                                 length(HC.sum$rt15stops[HC.sum$rt15stops!=0]),length(HC.sum$rt20stops[HC.sum$rt20stops!=0]),
                                 length(HC.sum$rt25stops[HC.sum$rt25stops!=0]),length(HC.sum$rt30stops[HC.sum$rt30stops!=0])),
                     changeindex = c(5,10,15,20,25,30)
                     
                  )


output$changepro <- output$stroute/output$hcroute*100

bbslist[[t]] <- output

}###end of the t loop


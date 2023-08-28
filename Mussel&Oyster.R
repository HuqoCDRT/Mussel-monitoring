#### Start ####


setwd("P:/Hugo/Data/data_chris")


#### Packages ####
library(patchwork)
library(exactRankTests)
library(ggplot2)
library(stats)
library(FSA)
library(readxl)
library(rstatix)
library(viridis)
library(TropFishR)
library(cohorts)
library(tidyverse)
library(rcompanion)
library(gridExtra)
library(tidyr)
library(dplyr)
library(ggsignif)
library(ggpubr)
library(pwr)
library(lsr)
library(exactRankTests)
library(reshape2)
library(cluster)    
library(factoextra) 
library(dendextend)
library(purrr)
library(lubridate)
library(ggsignif)

length_2022 <- read_excel("length_2022.xlsx")
View(length_2022)

lemu=length_2022[length_2022$species=="mussel",]
oys=length_2022[length_2022$species=="oyster",]






#Mussel
names(lemu)[5]="dates"  
lemu$dates=as.Date(lemu$dates,format="%Y%m%d")




#Oyster
names(oys)[5]="dates"  
oys$dates=as.Date(oys$dates,format="%Y%m%d")






############################## Mussel & Oyster measurements #########################################



######### Stats for 12 and 8 samples #########

######### 1. Mussel #########
##### 1.a Length #####

### 12 sample
sdm12<-aggregate(lemu$length,by=list(lemu$id),FUN=sd) #standard deviation

mean12=aggregate(lemu$length,by=list(lemu$id),FUN=mean) # mean of the length
res_lemu<-merge(mean12,sdm12,by="Group.1",all = T) #add mean and sd column to mussel data frame


### 8a sample
set.seed(1)
lemu8a=lemu %>% filter(sample %in% sample(unique(sample),8)) #8 random sample 

sdm8a<-aggregate(lemu8a$length,by=list(lemu8a$id),FUN=sd)

mean8a=aggregate(lemu8a$length,by=list(lemu8a$id),FUN=mean)
res_lemu8a<-merge(mean8a,sdm8a,by="Group.1",all = T)


### 8b sample
set.seed(12)
lemu8b=lemu %>% filter(sample %in% sample(unique(sample),8))

sdm8b<-aggregate(lemu8b$length,by=list(lemu8b$id),FUN=sd)

mean8b=aggregate(lemu8b$length,by=list(lemu8b$id),FUN=mean)
res_lemu8b<-merge(mean8b,sdm8b,by="Group.1",all = T)



### 8c sample
set.seed(123)
lemu8c=lemu %>% filter(sample %in% sample(unique(sample),8))

sdm8c<-aggregate(lemu8c$length,by=list(lemu8c$id),FUN=sd)

mean8c=aggregate(lemu8c$length,by=list(lemu8c$id),FUN=mean)
res_lemu8c<-merge(mean8c,sdm8c,by="Group.1",all = T)


list_df8<-list(res_lemu,res_lemu8a,res_lemu8b,res_lemu8c)
res_lemu12<-list_df8%>%
  reduce(inner_join,by="Group.1") #join all data frames by column Group.1


colnames(res_lemu12)<-c("id","mean_L12","sd_L12","mean_L8a","sd_L8a","mean_L8b","sd_L8b","mean_L8c","sd_L8c")



list_df8<-list(mean12,mean8a,mean8b,mean8c)
mean_length8<-list_df8%>%
  reduce(inner_join,by="Group.1")

colnames(mean_length8)<-c("id","meanL12","meanL8a","meanL8b","meanL8c")


melt8<-melt(mean_length8,id="id")
colnames(melt8)<-c("id","replicates","mean_length")

melt8%>% group_by(replicates) %>% shapiro_test(mean_length) # normally distributed

melt8 %>% t_test(mean_length~replicates,paired = T) #non significant














##### 1.b Density #####


banktype_list <- read_excel("banktype_list.xlsx")
names(banktype_list)[names(banktype_list) == 'STATION_NAME'] <- 'bank'




res_lemu12$bank <- as.vector(str_trim(str_extract(res_lemu12$id,"([:upper:][:upper:][:number:][:number:])"))) #create new column with something of id column
res_lemu12 <- merge(x = res_lemu12, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE) #join data frames



#add sub samp area column
require(dplyr)
res_lemu12 <- res_lemu12%>%
  mutate(sub_samp_area = case_when(grepl("Miesmuschelbank",res_lemu12$BANK_TYPE) ~ 176.7145867644259,
                                   grepl("Austernbank",res_lemu12$BANK_TYPE) ~ 625))










g1=group_by(length_2022,id,sample)
densityall=summarise(g1,length(sample)) #calculate number of ind for each sample






names(densityall)[names(densityall) == 'length(sample)'] <- 'n_mussels'


densityall$bank <- as.vector(str_trim(str_extract(densityall$id,"([:upper:][:upper:][:number:][:number:])")))
densityall <- merge(x = densityall, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)
require(dplyr)
densityall <- densityall%>%
  mutate(sub_samp_area = case_when(grepl("Miesmuschelbank",densityall$BANK_TYPE) ~ 176.7145867644259,
                                   grepl("Austernbank",densityall$BANK_TYPE) ~ 625))










# #add density column
densityall=densityall%>%
  mutate(density=n_mussels/sub_samp_area)




#add mean density to res_lemu
mean_D=aggregate(densityall$density,by=list(densityall$id),FUN=mean)
names(mean_D)[names(mean_D)=="Group.1"]<-"id"
densityall<-merge(mean_D,densityall,by="id",all = T)
names(densityall)[names(densityall)=="x"]<-"mean_D"






plotml8<- melt(res_lemu12, id.vars = c("id","bank"))






###density for 12 samples###
g=group_by(lemu,id,sample)
density12=summarise(g,length(sample)) #calculate number of ind for each sample



names(density12)[names(density12) == 'length(sample)'] <- 'n_mussels12'


density12$bank <- as.vector(str_trim(str_extract(density12$id,"([:upper:][:upper:][:number:][:number:])")))
density12 <- merge(x = density12, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)
require(dplyr)
density12 <- density12%>%
  mutate(sub_samp_area = case_when(grepl("Miesmuschelbank",density12$BANK_TYPE) ~ 176.7145867644259,
                                   grepl("Austernbank",density12$BANK_TYPE) ~ 625))






# #add density column
density12=density12%>%
  mutate(density=n_mussels12/sub_samp_area)




#add mean density to res_lemu
mean_D12=aggregate(density12$density,by=list(density12$id),FUN=mean)
names(mean_D12)[names(mean_D12)=="Group.1"]<-"id"
res_lemu12<-merge(mean_D12,res_lemu12,by="id",all = T)
names(res_lemu12)[names(res_lemu12)=="x"]<-"mean_D12"





#add sd density to res_lemu
sd_D12<-aggregate(density12$density,by=list(density12$id),FUN=sd)
names(sd_D12)[names(sd_D12)=="Group.1"]<-"id"
res_lemu12<-merge(sd_D12,res_lemu12,by="id",all = T)
names(res_lemu12)[names(res_lemu12)=="x"]<-"sd_D12"









###density for 8a samples###
set.seed(1)
g=group_by(lemu8a,id,sample)
density8a=summarise(g,length(sample))



names(density8a)[names(density8a) == 'length(sample)'] <- 'n_mussels8a'


density8a$bank <- as.vector(str_trim(str_extract(density8a$id,"([:upper:][:upper:][:number:][:number:])")))
density8a <- merge(x = density8a, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)
require(dplyr)
density8a <- density8a%>%
  mutate(sub_samp_area = case_when(grepl("Miesmuschelbank",density8a$BANK_TYPE) ~ 176.7145867644259,
                                   grepl("Austernbank",density8a$BANK_TYPE) ~ 625))





#add density column
density8a=density8a%>%
  mutate(density=n_mussels8a/sub_samp_area)



#add mean density to res_lemu
mean_D8a=aggregate(density8a$density,by=list(density8a$id),FUN=mean)
names(mean_D8a)[names(mean_D8a)=="Group.1"]<-"id"
res_lemu12<-merge(mean_D8a,res_lemu12,by="id",all = T)
names(res_lemu12)[names(res_lemu12)=="x"]<-"mean_D8a"





#add sd density to res_lemu
sd_D8a<-aggregate(density8a$density,by=list(density8a$id),FUN=sd)
names(sd_D8a)[names(sd_D8a)=="Group.1"]<-"id"
res_lemu12<-merge(sd_D8a,res_lemu12,by="id",all = T)
names(res_lemu12)[names(res_lemu12)=="x"]<-"sd_D8a"









###density for 8b samples###
set.seed(12)
g=group_by(lemu8b,id,sample)
density8b=summarise(g,length(sample))



names(density8b)[names(density8b) == 'length(sample)'] <- 'n_mussels8b'


density8b$bank <- as.vector(str_trim(str_extract(density8b$id,"([:upper:][:upper:][:number:][:number:])")))
density8b <- merge(x = density8b, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)
require(dplyr)
density8b <- density8b%>%
  mutate(sub_samp_area = case_when(grepl("Miesmuschelbank",density8b$BANK_TYPE) ~ 176.7145867644259,
                                   grepl("Austernbank",density8b$BANK_TYPE) ~ 625))





#add density column
density8b=density8b%>%
  mutate(density=n_mussels8b/sub_samp_area)




#add mean density to res_lemu
mean_D8b=aggregate(density8b$density,by=list(density8b$id),FUN=mean)
names(mean_D8b)[names(mean_D8b)=="Group.1"]<-"id"
res_lemu12<-merge(mean_D8b,res_lemu12,by="id",all = T)
names(res_lemu12)[names(res_lemu12)=="x"]<-"mean_D8b"





#add sd density to res_lemu
sd_D8b<-aggregate(density8b$density,by=list(density8b$id),FUN=sd)
names(sd_D8b)[names(sd_D8b)=="Group.1"]<-"id"
res_lemu12<-merge(sd_D8b,res_lemu12,by="id",all = T)
names(res_lemu12)[names(res_lemu12)=="x"]<-"sd_D8b"









###density for 8c samples###
set.seed(13)
g=group_by(lemu8c,id,sample)
density8c=summarise(g,length(sample))



names(density8c)[names(density8c) == 'length(sample)'] <- 'n_mussels8c'


density8c$bank <- as.vector(str_trim(str_extract(density8c$id,"([:upper:][:upper:][:number:][:number:])")))
density8c <- merge(x = density8c, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)
require(dplyr)
density8c <- density8c%>%
  mutate(sub_samp_area = case_when(grepl("Miesmuschelbank",density8c$BANK_TYPE) ~ 176.7145867644259,
                                   grepl("Austernbank",density8c$BANK_TYPE) ~ 625))





#add density column
density8c=density8c%>%
  mutate(density=n_mussels8c/sub_samp_area)




#add mean density to res_lemu
mean_D8c=aggregate(density8c$density,by=list(density8c$id),FUN=mean)
names(mean_D8c)[names(mean_D8c)=="Group.1"]<-"id"
res_lemu12<-merge(mean_D8c,res_lemu12,by="id",all = T)
names(res_lemu12)[names(res_lemu12)=="x"]<-"mean_D8c"





#add sd density to res_lemu
sd_D8c<-aggregate(density8c$density,by=list(density8c$id),FUN=sd)
names(sd_D8c)[names(sd_D8c)=="Group.1"]<-"id"
res_lemu12<-merge(sd_D8c,res_lemu12,by="id",all = T)
names(res_lemu12)[names(res_lemu12)=="x"]<-"sd_D8c"




list_dfD8<-list(mean_D12,mean_D8a,mean_D8b,mean_D8c)
mean_density8<-list_dfD8%>%
  reduce(inner_join,by="id")

colnames(mean_density8)<-c("id","meanD12","meanD8a","meanD8b","meanD8c")

melt8d<-melt(mean_density8,id="id")
colnames(melt8d)<-c("id","replicates","mean_density")

melt8d%>% group_by(replicates) %>% shapiro_test(mean_density) # non normally distributed

melt8d %>% wilcox_test(mean_density~replicates,paired = T) #non significant




















######### 2. Oyster ##########
##### 2.a Length #####

### 12 sample
sdo12<-aggregate(oys$length,by=list(oys$id),FUN=sd)

mean12o=aggregate(oys$length,by=list(oys$id),FUN=mean)
res_oys<-merge(mean12o,sdo12,by="Group.1",all = T)


### 8a sample
set.seed(1)
oys8a=oys %>% filter(sample %in% sample(unique(sample),8))

sdo8a<-aggregate(oys8a$length,by=list(oys8a$id),FUN=sd)

mean8ao=aggregate(oys8a$length,by=list(oys8a$id),FUN=mean)
res_oys8a<-merge(mean8ao,sdo8a,by="Group.1",all = T)


### 8b sample
set.seed(12)
oys8b=oys %>% filter(sample %in% sample(unique(sample),8))

sdo8b<-aggregate(oys8b$length,by=list(oys8b$id),FUN=sd)

mean8bo=aggregate(oys8b$length,by=list(oys8b$id),FUN=mean)
res_oys8b<-merge(mean8bo,sdo8b,by="Group.1",all = T)



### 8c sample
set.seed(13)
oys8c=oys %>% filter(sample %in% sample(unique(sample),8))

sdo8c<-aggregate(oys8c$length,by=list(oys8c$id),FUN=sd)

mean8co=aggregate(oys8c$length,by=list(oys8c$id),FUN=mean)
res_oys8c<-merge(mean8co,sdo8c,by="Group.1",all = T)


list_dfO<-list(res_oys,res_oys8a,res_oys8b,res_oys8c)
res_oys12<-list_dfO%>%
  reduce(inner_join,by="Group.1")


colnames(res_oys12)<-c("id","mean_Lo12","sd_Lo12","mean_Lo8a","sd_Lo8a","mean_Lo8b","sd_Lo8b","mean_Lo8c","sd_Lo8c")




list_dfO<-list(mean12o,mean8ao,mean8bo,mean8co)
mean_lengthO<-list_dfO%>%
  reduce(inner_join,by="Group.1")

colnames(mean_lengthO)<-c("id","mean12","mean8a","mean8b","mean8c")



melt8o<-melt(mean_lengthO,id="id")
colnames(melt8o)<-c("id","replicates","mean_length")

melt8o%>% group_by(replicates) %>% shapiro_test(mean_length) # not normally distributed

melt8o %>% wilcox_test(mean_length~replicates,paired = T) #non significant











##### 2.b Density #####



banktype_list <- read_excel("banktype_list.xlsx")
names(banktype_list)[names(banktype_list) == 'STATION_NAME'] <- 'bank'



res_oys12$bank <- as.vector(str_trim(str_extract(res_oys12$id,"([:upper:][:upper:][:number:][:number:])")))
res_oys12 <- merge(x = res_oys12, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)



#add sub samp area column
require(dplyr)
res_oys12 <- res_oys12%>%
  mutate(sub_samp_area = case_when(grepl("Miesmuschelbank",res_oys12$BANK_TYPE) ~ 176.7145867644259,
                                   grepl("Austernbank",res_oys12$BANK_TYPE) ~ 625))








###density for 12 samples###
g=group_by(oys,id,sample)
density12o=summarise(g,length(sample))



names(density12o)[names(density12o) == 'length(sample)'] <- 'n_Oysters12'


density12o$bank <- as.vector(str_trim(str_extract(density12o$id,"([:upper:][:upper:][:number:][:number:])")))
density12o <- merge(x = density12o, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)
require(dplyr)
density12o <- density12o%>%
  mutate(sub_samp_area = case_when(grepl("Miesmuschelbank",density12o$BANK_TYPE) ~ 176.7145867644259,
                                   grepl("Austernbank",density12o$BANK_TYPE) ~ 625))






# #add density column
density12o=density12o%>%
  mutate(density=n_Oysters12/sub_samp_area)




#add mean density to res_oys12
mean_Do12=aggregate(density12o$density,by=list(density12o$id),FUN=mean)
names(mean_Do12)[names(mean_Do12)=="Group.1"]<-"id"
res_oys12<-merge(mean_Do12,res_oys12,by="id",all = T)
names(res_oys12)[names(res_oys12)=="x"]<-"mean_Do12"





#add sd density to res_oys12
sd_Do12<-aggregate(density12o$density,by=list(density12o$id),FUN=sd)
names(sd_Do12)[names(sd_Do12)=="Group.1"]<-"id"
res_oys12<-merge(sd_Do12,res_oys12,by="id",all = T)
names(res_oys12)[names(res_oys12)=="x"]<-"sd_Do12"








###density for 8a samples###
set.seed(1)
g=group_by(oys8a,id,sample)
density8ao=summarise(g,length(sample))



names(density8ao)[names(density8ao) == 'length(sample)'] <- 'n_Oysters8a'


density8ao$bank <- as.vector(str_trim(str_extract(density8ao$id,"([:upper:][:upper:][:number:][:number:])")))
density8ao <- merge(x = density8ao, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)
require(dplyr)
density8ao <- density8ao%>%
  mutate(sub_samp_area = case_when(grepl("Miesmuschelbank",density8ao$BANK_TYPE) ~ 176.7145867644259,
                                   grepl("Austernbank",density8ao$BANK_TYPE) ~ 625))





#add density column
density8ao=density8ao%>%
  mutate(density=n_Oysters8a/sub_samp_area)



#add mean density to res_oys12
mean_Do8a=aggregate(density8ao$density,by=list(density8ao$id),FUN=mean)
names(mean_Do8a)[names(mean_Do8a)=="Group.1"]<-"id"
res_oys12<-merge(mean_Do8a,res_oys12,by="id",all = T)
names(res_oys12)[names(res_oys12)=="x"]<-"mean_Do8a"





#add sd density to res_oys12
sd_Do8a<-aggregate(density8ao$density,by=list(density8ao$id),FUN=sd)
names(sd_Do8a)[names(sd_Do8a)=="Group.1"]<-"id"
res_oys12<-merge(sd_Do8a,res_oys12,by="id",all = T)
names(res_oys12)[names(res_oys12)=="x"]<-"sd_Do8a"







###density for 8b samples###
set.seed(12)
g=group_by(oys8b,id,sample)
density8bo=summarise(g,length(sample))



names(density8bo)[names(density8bo) == 'length(sample)'] <- 'n_Oysters8b'


density8bo$bank <- as.vector(str_trim(str_extract(density8bo$id,"([:upper:][:upper:][:number:][:number:])")))
density8bo <- merge(x = density8bo, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)
require(dplyr)
density8bo <- density8bo%>%
  mutate(sub_samp_area = case_when(grepl("Miesmuschelbank",density8bo$BANK_TYPE) ~ 176.7145867644259,
                                   grepl("Austernbank",density8bo$BANK_TYPE) ~ 625))





#add density column
density8bo=density8bo%>%
  mutate(density=n_Oysters8b/sub_samp_area)




#add mean density to res_oys12
mean_Do8b=aggregate(density8bo$density,by=list(density8bo$id),FUN=mean)
names(mean_Do8b)[names(mean_Do8b)=="Group.1"]<-"id"
res_oys12<-merge(mean_Do8b,res_oys12,by="id",all = T)
names(res_oys12)[names(res_oys12)=="x"]<-"mean_Do8b"





#add sd density to res_oys12
sd_Do8b<-aggregate(density8bo$density,by=list(density8bo$id),FUN=sd)
names(sd_Do8b)[names(sd_Do8b)=="Group.1"]<-"id"
res_oys12<-merge(sd_Do8b,res_oys12,by="id",all = T)
names(res_oys12)[names(res_oys12)=="x"]<-"sd_Do8b"







###density for 8c samples###
set.seed(13)
g=group_by(oys8c,id,sample)
density8co=summarise(g,length(sample))



names(density8co)[names(density8co) == 'length(sample)'] <- 'n_Oysters8c'


density8co$bank <- as.vector(str_trim(str_extract(density8co$id,"([:upper:][:upper:][:number:][:number:])")))
density8co <- merge(x = density8co, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)
require(dplyr)
density8co <- density8co%>%
  mutate(sub_samp_area = case_when(grepl("Miesmuschelbank",density8co$BANK_TYPE) ~ 176.7145867644259,
                                   grepl("Austernbank",density8co$BANK_TYPE) ~ 625))





#add density column
density8co=density8co%>%
  mutate(density=n_Oysters8c/sub_samp_area)




#add mean density to res_oys12
mean_Do8c=aggregate(density8co$density,by=list(density8co$id),FUN=mean)
names(mean_Do8c)[names(mean_Do8c)=="Group.1"]<-"id"
res_oys12<-merge(mean_Do8c,res_oys12,by="id",all = T)
names(res_oys12)[names(res_oys12)=="x"]<-"mean_Do8c"





#add sd density to res_oys12
sd_Do8c<-aggregate(density8co$density,by=list(density8co$id),FUN=sd)
names(sd_Do8c)[names(sd_Do8c)=="Group.1"]<-"id"
res_oys12<-merge(sd_Do8c,res_oys12,by="id",all = T)
names(res_oys12)[names(res_oys12)=="x"]<-"sd_Do8c"







list_dfDO<-list(mean_Do12,mean_Do8a,mean_Do8b,mean_Do8c)
mean_densityO<-list_dfDO%>%
  reduce(inner_join,by="id")

colnames(mean_densityO)<-c("id","meanD12","meanD8a","meanD8b","meanD8c")



melt8do<-melt(mean_densityO,id="id")
colnames(melt8do)<-c("id","replicates","mean_density")

melt8do%>% group_by(replicates) %>% shapiro_test(mean_density) # non normally distributed

melt8do %>% wilcox_test(mean_density~replicates,paired = T) #non significant








######### Stats for 12 and 6 samples #########

######### 1. Mussel #########
##### 1.a Length #####
set.seed(1)
lemu6a=lemu %>% filter(sample %in% sample(unique(sample),6))
set.seed(12)
lemu6b=lemu %>% filter(sample %in% sample(unique(sample),6))
set.seed(13)
lemu6c=lemu %>% filter(sample %in% sample(unique(sample),6))


## 6a sample
sdm6a<-aggregate(lemu6a$length,by=list(lemu6a$id),FUN=sd)

mean6a=aggregate(lemu6a$length,by=list(lemu6a$id),FUN=mean)
res_lemu6a<-merge(mean6a,sdm6a,by="Group.1",all = T)



## 6b sample
sdm6b<-aggregate(lemu6b$length,by=list(lemu6b$id),FUN=sd)

mean6b=aggregate(lemu6b$length,by=list(lemu6b$id),FUN=mean)
res_lemu6b<-merge(mean6b,sdm6b,by="Group.1",all = T)



## 6c sample
sdm6c<-aggregate(lemu6c$length,by=list(lemu6c$id),FUN=sd)

mean6c=aggregate(lemu6c$length,by=list(lemu6c$id),FUN=mean)
res_lemu6c<-merge(mean6c,sdm6c,by="Group.1",all = T)


list_df6<-list(res_lemu,res_lemu6a,res_lemu6b,res_lemu6c)
res_lemu6<-list_df6%>%
  reduce(inner_join,by="Group.1") #join all data frames by column Group.1


colnames(res_lemu6)<-c("id","mean_L12","sd_L12","mean_L6a","sd_L6a","mean_L6b","sd_L6b","mean_L6c","sd_L6c")






list_df6bis<-list(mean12,mean6a,mean6b,mean6c)
mean_length6<-list_df6bis%>%
  reduce(inner_join,by="Group.1")

colnames(mean_length6)<-c("id","meanL12","meanL6a","meanL6b","meanL6c")

melt6<-melt(mean_length6,id="id")
colnames(melt6)<-c("id","replicates","mean_length")

melt6%>% group_by(replicates) %>% shapiro_test(mean_length) # not normally distributed

melt6 %>% t_test(mean_length~replicates,paired = T) #non significant






##### 1.b Density #####


###density for 6a samples###
set.seed(1)
g=group_by(lemu6a,id,sample)
density6a=summarise(g,length(sample))



names(density6a)[names(density6a) == 'length(sample)'] <- 'n_mussels6a'


density6a$bank <- as.vector(str_trim(str_extract(density6a$id,"([:upper:][:upper:][:number:][:number:])")))
density6a <- merge(x = density6a, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)
require(dplyr)
density6a <- density6a%>%
  mutate(sub_samp_area = case_when(grepl("Miesmuschelbank",density6a$BANK_TYPE) ~ 176.7145667644259,
                                   grepl("Austernbank",density6a$BANK_TYPE) ~ 625))





#add density column
density6a=density6a%>%
  mutate(density=n_mussels6a/sub_samp_area)



#add mean density to res_lemu6
mean_D6a=aggregate(density6a$density,by=list(density6a$id),FUN=mean)
names(mean_D6a)[names(mean_D6a)=="Group.1"]<-"id"
res_lemu6<-merge(mean_D6a,res_lemu6,by="id",all = T)
names(res_lemu6)[names(res_lemu6)=="x"]<-"mean_D6a"





#add sd density to res_lemu6
sd_D6a<-aggregate(density6a$density,by=list(density6a$id),FUN=sd)
names(sd_D6a)[names(sd_D6a)=="Group.1"]<-"id"
res_lemu6<-merge(sd_D6a,res_lemu6,by="id",all = T)
names(res_lemu6)[names(res_lemu6)=="x"]<-"sd_D6a"





###density for 6b samples###
set.seed(12)
g=group_by(lemu6b,id,sample)
density6b=summarise(g,length(sample))



names(density6b)[names(density6b) == 'length(sample)'] <- 'n_mussels6b'


density6b$bank <- as.vector(str_trim(str_extract(density6b$id,"([:upper:][:upper:][:number:][:number:])")))
density6b <- merge(x = density6b, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)
require(dplyr)
density6b <- density6b%>%
  mutate(sub_samp_area = case_when(grepl("Miesmuschelbank",density6b$BANK_TYPE) ~ 176.7145667644259,
                                   grepl("Austernbank",density6b$BANK_TYPE) ~ 625))




#add density column
density6b=density6b%>%
  mutate(density=n_mussels6b/sub_samp_area)




#add mean density to res_lemu6
mean_D6b=aggregate(density6b$density,by=list(density6b$id),FUN=mean)
names(mean_D6b)[names(mean_D6b)=="Group.1"]<-"id"
res_lemu6<-merge(mean_D6b,res_lemu6,by="id",all = T)
names(res_lemu6)[names(res_lemu6)=="x"]<-"mean_D6b"





#add sd density to res_lemu6
sd_D6b<-aggregate(density6b$density,by=list(density6b$id),FUN=sd)
names(sd_D6b)[names(sd_D6b)=="Group.1"]<-"id"
res_lemu6<-merge(sd_D6b,res_lemu6,by="id",all = T)
names(res_lemu6)[names(res_lemu6)=="x"]<-"sd_D6b"






###density for 6c samples###
set.seed(13)
g=group_by(lemu6c,id,sample)
density6c=summarise(g,length(sample))



names(density6c)[names(density6c) == 'length(sample)'] <- 'n_mussels6c'


density6c$bank <- as.vector(str_trim(str_extract(density6c$id,"([:upper:][:upper:][:number:][:number:])")))
density6c <- merge(x = density6c, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)
require(dplyr)
density6c <- density6c%>%
  mutate(sub_samp_area = case_when(grepl("Miesmuschelbank",density6c$BANK_TYPE) ~ 176.7145667644259,
                                   grepl("Austernbank",density6c$BANK_TYPE) ~ 625))





#add density column
density6c=density6c%>%
  mutate(density=n_mussels6c/sub_samp_area)




#add mean density to res_lemu6
mean_D6c=aggregate(density6c$density,by=list(density6c$id),FUN=mean)
names(mean_D6c)[names(mean_D6c)=="Group.1"]<-"id"
res_lemu6<-merge(mean_D6c,res_lemu6,by="id",all = T)
names(res_lemu6)[names(res_lemu6)=="x"]<-"mean_D6c"
res_lemu6<-merge(mean_D12,res_lemu6,by="id",all = T)
names(res_lemu6)[names(res_lemu6)=="x"]<-"mean_D12"




#add sd density to res_lemu6
sd_D6c<-aggregate(density6c$density,by=list(density6c$id),FUN=sd)
names(sd_D6c)[names(sd_D6c)=="Group.1"]<-"id"
res_lemu6<-merge(sd_D6c,res_lemu6,by="id",all = T)
names(res_lemu6)[names(res_lemu6)=="x"]<-"sd_D6c"
res_lemu6<-merge(sd_D12,res_lemu6,by="id",all = T)
names(res_lemu6)[names(res_lemu6)=="x"]<-"sd_D12"








list_dfD6bis<-list(mean_D12,mean_D6a,mean_D6b,mean_D6c)
mean_density6<-list_dfD6bis%>%
  reduce(inner_join,by="id")

colnames(mean_density6)<-c("id","meanD12","meanD6a","meanD6b","meanD6c")



melt6d<-melt(mean_density6,id="id")
colnames(melt6d)<-c("id","replicates","mean_density")

melt6d%>% group_by(replicates) %>% shapiro_test(mean_density) # non normally distributed

melt6d %>% wilcox_test(mean_density~replicates,paired = T) #non significant



















######### 2. Oyster ##########
##### 2.a Length #####

set.seed(1)
oys6a=oys %>% filter(sample %in% sample(unique(sample),6))
set.seed(12)
oys6b=oys %>% filter(sample %in% sample(unique(sample),6))
set.seed(13)
oys6c=oys %>% filter(sample %in% sample(unique(sample),6))


## 6a sample
sdo6a<-aggregate(oys6a$length,by=list(oys6a$id),FUN=sd)

mean6ao=aggregate(oys6a$length,by=list(oys6a$id),FUN=mean)
res_oys6a<-merge(mean6ao,sdo6a,by="Group.1",all = T)



## 6b sample
sdo6b<-aggregate(oys6b$length,by=list(oys6b$id),FUN=sd)

mean6bo=aggregate(oys6b$length,by=list(oys6b$id),FUN=mean)
res_oys6b<-merge(mean6bo,sdo6b,by="Group.1",all = T)



## 6c sample
sdo6c<-aggregate(oys6c$length,by=list(oys6c$id),FUN=sd)

mean6co=aggregate(oys6c$length,by=list(oys6c$id),FUN=mean)
res_oys6c<-merge(mean6co,sdo6c,by="Group.1",all = T)


list_dfO6<-list(res_oys,res_oys6a,res_oys6b,res_oys6c)
res_oys6<-list_dfO6%>%
  reduce(inner_join,by="Group.1")


colnames(res_oys6)<-c("id","mean_Lo12","sd_Lo12","mean_Lo6a","sd_Lo6a","mean_Lo6b","sd_Lo6b","mean_Lo6c","sd_Lo6c")






list_dfO6<-list(mean12o,mean6ao,mean6bo,mean6co)
mean_lengthO6<-list_dfO6%>%
  reduce(inner_join,by="Group.1")

colnames(mean_lengthO6)<-c("id","mean12","mean6a","mean6b","mean6c")




melt6o<-melt(mean_lengthO6,id="id")
colnames(melt6o)<-c("id","replicates","mean_length")

melt6o%>% group_by(replicates) %>% shapiro_test(mean_length) # not normally distributed

melt6o %>% wilcox_test(mean_length~replicates,paired = T) #non significant









##### 2.b Density #####



###density for 6a samples###
set.seed(1)
g=group_by(oys6a,id,sample)
density6ao=summarise(g,length(sample))



names(density6ao)[names(density6ao) == 'length(sample)'] <- 'n_Oysters6a'


density6ao$bank <- as.vector(str_trim(str_extract(density6ao$id,"([:upper:][:upper:][:number:][:number:])")))
density6ao <- merge(x = density6ao, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)
require(dplyr)
density6ao <- density6ao%>%
  mutate(sub_samp_area = case_when(grepl("Miesmuschelbank",density6ao$BANK_TYPE) ~ 176.7145667644259,
                                   grepl("Austernbank",density6ao$BANK_TYPE) ~ 625))





#add density column
density6ao=density6ao%>%
  mutate(density=n_Oysters6a/sub_samp_area)



#add mean density to res_oys6
mean_Do6a=aggregate(density6ao$density,by=list(density6ao$id),FUN=mean)
names(mean_Do6a)[names(mean_Do6a)=="Group.1"]<-"id"
res_oys6<-merge(mean_Do6a,res_oys6,by="id",all = T)
names(res_oys6)[names(res_oys6)=="x"]<-"mean_Do6a"





#add sd density to res_oys6
sd_Do6a<-aggregate(density6ao$density,by=list(density6ao$id),FUN=sd)
names(sd_Do6a)[names(sd_Do6a)=="Group.1"]<-"id"
res_oys6<-merge(sd_Do6a,res_oys6,by="id",all = T)
names(res_oys6)[names(res_oys6)=="x"]<-"sd_Do6a"









###density for 6b samples###
set.seed(12)
g=group_by(oys6b,id,sample)
density6bo=summarise(g,length(sample))



names(density6bo)[names(density6bo) == 'length(sample)'] <- 'n_Oysters6b'


density6bo$bank <- as.vector(str_trim(str_extract(density6bo$id,"([:upper:][:upper:][:number:][:number:])")))
density6bo <- merge(x = density6bo, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)
require(dplyr)
density6bo <- density6bo%>%
  mutate(sub_samp_area = case_when(grepl("Miesmuschelbank",density6bo$BANK_TYPE) ~ 176.7145667644259,
                                   grepl("Austernbank",density6bo$BANK_TYPE) ~ 625))





#add density column
density6bo=density6bo%>%
  mutate(density=n_Oysters6b/sub_samp_area)




#add mean density to res_oys6
mean_Do6b=aggregate(density6bo$density,by=list(density6bo$id),FUN=mean)
names(mean_Do6b)[names(mean_Do6b)=="Group.1"]<-"id"
res_oys6<-merge(mean_Do6b,res_oys6,by="id",all = T)
names(res_oys6)[names(res_oys6)=="x"]<-"mean_Do6b"





#add sd density to res_oys6
sd_Do6b<-aggregate(density6bo$density,by=list(density6bo$id),FUN=sd)
names(sd_Do6b)[names(sd_Do6b)=="Group.1"]<-"id"
res_oys6<-merge(sd_Do6b,res_oys6,by="id",all = T)
names(res_oys6)[names(res_oys6)=="x"]<-"sd_Do6b"









###density for 6c samples###
set.seed(13)
g=group_by(oys6c,id,sample)
density6co=summarise(g,length(sample))



names(density6co)[names(density6co) == 'length(sample)'] <- 'n_Oysters6c'


density6co$bank <- as.vector(str_trim(str_extract(density6co$id,"([:upper:][:upper:][:number:][:number:])")))
density6co <- merge(x = density6co, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)
require(dplyr)
density6co <- density6co%>%
  mutate(sub_samp_area = case_when(grepl("Miesmuschelbank",density6co$BANK_TYPE) ~ 176.7145667644259,
                                   grepl("Austernbank",density6co$BANK_TYPE) ~ 625))





#add density column
density6co=density6co%>%
  mutate(density=n_Oysters6c/sub_samp_area)




#add mean density to res_oys6
mean_Do6c=aggregate(density6co$density,by=list(density6co$id),FUN=mean)
names(mean_Do6c)[names(mean_Do6c)=="Group.1"]<-"id"
res_oys6<-merge(mean_Do6c,res_oys6,by="id",all = T)
names(res_oys6)[names(res_oys6)=="x"]<-"mean_Do6c"

res_oys6<-merge(mean_Do12,res_oys6,by="id",all = T)
names(res_oys6)[names(res_oys6)=="x"]<-"mean_Do12"



#add sd density to res_oys6
sd_Do6c<-aggregate(density6co$density,by=list(density6co$id),FUN=sd)
names(sd_Do6c)[names(sd_Do6c)=="Group.1"]<-"id"
res_oys6<-merge(sd_Do6c,res_oys6,by="id",all = T)
names(res_oys6)[names(res_oys6)=="x"]<-"sd_Do6c"

res_oys6<-merge(sd_Do12,res_oys6,by="id",all = T)
names(res_oys6)[names(res_oys6)=="x"]<-"sd_Do12"



list_dfDO6<-list(mean_Do12,mean_Do6a,mean_Do6b,mean_Do6c)
mean_densityO6<-list_dfDO6%>%
  reduce(inner_join,by="id")

colnames(mean_densityO6)<-c("id","meanD12","meanD6a","meanD6b","meanD6c")



melt6do<-melt(mean_densityO6,id="id")
colnames(melt6do)<-c("id","replicates","mean_density")

melt6do%>% group_by(replicates) %>% shapiro_test(mean_density) # non normally distributed

melt6do %>% wilcox_test(mean_density~replicates,paired = T) #non significant









######### Stats for 12 and 5 samples ##########

######### 1. Mussel #########
##### 1.a Length #####

set.seed(1)
lemu5a=lemu %>% filter(sample %in% sample(unique(sample),5))
set.seed(12)
lemu5b=lemu %>% filter(sample %in% sample(unique(sample),5))
set.seed(13)
lemu5c=lemu %>% filter(sample %in% sample(unique(sample),5))


## 5a sample
sdm5a<-aggregate(lemu5a$length,by=list(lemu5a$id),FUN=sd)

mean5a=aggregate(lemu5a$length,by=list(lemu5a$id),FUN=mean)
res_lemu5a<-merge(mean5a,sdm5a,by="Group.1",all = T)



## 5b sample
sdm5b<-aggregate(lemu5b$length,by=list(lemu5b$id),FUN=sd)

mean5b=aggregate(lemu5b$length,by=list(lemu5b$id),FUN=mean)
res_lemu5b<-merge(mean5b,sdm5b,by="Group.1",all = T)



## 5c sample
sdm5c<-aggregate(lemu5c$length,by=list(lemu5c$id),FUN=sd)

mean5c=aggregate(lemu5c$length,by=list(lemu5c$id),FUN=mean)
res_lemu5c<-merge(mean5c,sdm5c,by="Group.1",all = T)


list_df5<-list(res_lemu,res_lemu5a,res_lemu5b,res_lemu5c)
res_lemu5<-list_df5%>%
  reduce(inner_join,by="Group.1") #join all data frames by column Group.1


colnames(res_lemu5)<-c("id","mean_L12","sd_L12","mean_L5a","sd_L5a","mean_L5b","sd_L5b","mean_L5c","sd_L5c")






list_df5bis<-list(mean12,mean5a,mean5b,mean5c)
mean_length5<-list_df5bis%>%
  reduce(inner_join,by="Group.1")

colnames(mean_length5)<-c("id","meanL12","meanL5a","meanL5b","meanL5c")



melt5<-melt(mean_length5,id="id")
colnames(melt5)<-c("id","replicates","mean_length")

melt5%>% group_by(replicates) %>% shapiro_test(mean_length) # not normally distributed

melt5 %>% t_test(mean_length~replicates,paired = T) #non significant








##### 1.b Density #####


###density for 5a samples###
set.seed(1)
g=group_by(lemu5a,id,sample)
density5a=summarise(g,length(sample))



names(density5a)[names(density5a) == 'length(sample)'] <- 'n_mussels5a'


density5a$bank <- as.vector(str_trim(str_extract(density5a$id,"([:upper:][:upper:][:number:][:number:])")))
density5a <- merge(x = density5a, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)
require(dplyr)
density5a <- density5a%>%
  mutate(sub_samp_area = case_when(grepl("Miesmuschelbank",density5a$BANK_TYPE) ~ 175.7155557555259,
                                   grepl("Austernbank",density5a$BANK_TYPE) ~ 525))





#add density column
density5a=density5a%>%
  mutate(density=n_mussels5a/sub_samp_area)



#add mean density to res_lemu5
mean_D5a=aggregate(density5a$density,by=list(density5a$id),FUN=mean)
names(mean_D5a)[names(mean_D5a)=="Group.1"]<-"id"
res_lemu5<-merge(mean_D5a,res_lemu5,by="id",all = T)
names(res_lemu5)[names(res_lemu5)=="x"]<-"mean_D5a"





#add sd density to res_lemu5
sd_D5a<-aggregate(density5a$density,by=list(density5a$id),FUN=sd)
names(sd_D5a)[names(sd_D5a)=="Group.1"]<-"id"
res_lemu5<-merge(sd_D5a,res_lemu5,by="id",all = T)
names(res_lemu5)[names(res_lemu5)=="x"]<-"sd_D5a"









###density for 5b samples###
set.seed(12)
g=group_by(lemu5b,id,sample)
density5b=summarise(g,length(sample))



names(density5b)[names(density5b) == 'length(sample)'] <- 'n_mussels5b'


density5b$bank <- as.vector(str_trim(str_extract(density5b$id,"([:upper:][:upper:][:number:][:number:])")))
density5b <- merge(x = density5b, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)
require(dplyr)
density5b <- density5b%>%
  mutate(sub_samp_area = case_when(grepl("Miesmuschelbank",density5b$BANK_TYPE) ~ 175.7155557555259,
                                   grepl("Austernbank",density5b$BANK_TYPE) ~ 525))





#add density column
density5b=density5b%>%
  mutate(density=n_mussels5b/sub_samp_area)




#add mean density to res_lemu5
mean_D5b=aggregate(density5b$density,by=list(density5b$id),FUN=mean)
names(mean_D5b)[names(mean_D5b)=="Group.1"]<-"id"
res_lemu5<-merge(mean_D5b,res_lemu5,by="id",all = T)
names(res_lemu5)[names(res_lemu5)=="x"]<-"mean_D5b"





#add sd density to res_lemu5
sd_D5b<-aggregate(density5b$density,by=list(density5b$id),FUN=sd)
names(sd_D5b)[names(sd_D5b)=="Group.1"]<-"id"
res_lemu5<-merge(sd_D5b,res_lemu5,by="id",all = T)
names(res_lemu5)[names(res_lemu5)=="x"]<-"sd_D5b"









###density for 5c samples###
set.seed(13)
g=group_by(lemu5c,id,sample)
density5c=summarise(g,length(sample))



names(density5c)[names(density5c) == 'length(sample)'] <- 'n_mussels5c'


density5c$bank <- as.vector(str_trim(str_extract(density5c$id,"([:upper:][:upper:][:number:][:number:])")))
density5c <- merge(x = density5c, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)
require(dplyr)
density5c <- density5c%>%
  mutate(sub_samp_area = case_when(grepl("Miesmuschelbank",density5c$BANK_TYPE) ~ 175.7155557555259,
                                   grepl("Austernbank",density5c$BANK_TYPE) ~ 525))





#add density column
density5c=density5c%>%
  mutate(density=n_mussels5c/sub_samp_area)




#add mean density to res_lemu5
mean_D5c=aggregate(density5c$density,by=list(density5c$id),FUN=mean)
names(mean_D5c)[names(mean_D5c)=="Group.1"]<-"id"
res_lemu5<-merge(mean_D5c,res_lemu5,by="id",all = T)
names(res_lemu5)[names(res_lemu5)=="x"]<-"mean_D5c"
res_lemu5<-merge(mean_D12,res_lemu5,by="id",all = T)
names(res_lemu5)[names(res_lemu5)=="x"]<-"mean_D12"



#add sd density to res_lemu5
sd_D5c<-aggregate(density5c$density,by=list(density5c$id),FUN=sd)
names(sd_D5c)[names(sd_D5c)=="Group.1"]<-"id"
res_lemu5<-merge(sd_D5c,res_lemu5,by="id",all = T)
names(res_lemu5)[names(res_lemu5)=="x"]<-"sd_D5c"
res_lemu5<-merge(sd_D12,res_lemu5,by="id",all = T)
names(res_lemu5)[names(res_lemu5)=="x"]<-"sd_D12"




list_dfD5bis<-list(mean_D12,mean_D5a,mean_D5b,mean_D5c)
mean_density5<-list_dfD5bis%>%
  reduce(inner_join,by="id")

colnames(mean_density5)<-c("id","meanD12","meanD5a","meanD5b","meanD5c")



melt5d<-melt(mean_density5,id="id")
colnames(melt5d)<-c("id","replicates","mean_density")

melt5d%>% group_by(replicates) %>% shapiro_test(mean_density) # non normally distributed

melt5d %>% wilcox_test(mean_density~replicates,paired = T) #significant




















####### 2. Oyster #######
#### 2.a Length ####

set.seed(1)
oys5a=oys %>% filter(sample %in% sample(unique(sample),5))
set.seed(12)
oys5b=oys %>% filter(sample %in% sample(unique(sample),5))
set.seed(13)
oys5c=oys %>% filter(sample %in% sample(unique(sample),5))


## 5a sample
sdo5a<-aggregate(oys5a$length,by=list(oys5a$id),FUN=sd)

mean5ao=aggregate(oys5a$length,by=list(oys5a$id),FUN=mean)
res_oys5a<-merge(mean5ao,sdo5a,by="Group.1",all = T)



## 5b sample
sdo5b<-aggregate(oys5b$length,by=list(oys5b$id),FUN=sd)

mean5bo=aggregate(oys5b$length,by=list(oys5b$id),FUN=mean)
res_oys5b<-merge(mean5bo,sdo5b,by="Group.1",all = T)



## 5c sample
sdo5c<-aggregate(oys5c$length,by=list(oys5c$id),FUN=sd)

mean5co=aggregate(oys5c$length,by=list(oys5c$id),FUN=mean)
res_oys5c<-merge(mean5co,sdo5c,by="Group.1",all = T)


list_dfO5<-list(res_oys,res_oys5a,res_oys5b,res_oys5c)
res_oys5<-list_dfO5%>%
  reduce(inner_join,by="Group.1")


colnames(res_oys5)<-c("id","mean_Lo12","sd_Lo12","mean_Lo5a","sd_Lo5a","mean_Lo5b","sd_Lo5b","mean_Lo5c","sd_Lo5c")








list_dfO5<-list(mean12o,mean5ao,mean5bo,mean5co)
mean_lengthO5<-list_dfO5%>%
  reduce(inner_join,by="Group.1")

colnames(mean_lengthO5)<-c("id","mean12","mean5a","mean5b","mean5c")





melt5o<-melt(mean_lengthO5,id="id")
colnames(melt5o)<-c("id","replicates","mean_length")

melt5o%>% group_by(replicates) %>% shapiro_test(mean_length) # not normally distributed

melt5o %>% wilcox_test(mean_length~replicates,paired = T) #non significant











#### 2.b Density ####


###density for 5a samples###
set.seed(1)
g=group_by(oys5a,id,sample)
density5ao=summarise(g,length(sample))



names(density5ao)[names(density5ao) == 'length(sample)'] <- 'n_Oysters5a'


density5ao$bank <- as.vector(str_trim(str_extract(density5ao$id,"([:upper:][:upper:][:number:][:number:])")))
density5ao <- merge(x = density5ao, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)
require(dplyr)
density5ao <- density5ao%>%
  mutate(sub_samp_area = case_when(grepl("Miesmuschelbank",density5ao$BANK_TYPE) ~ 175.7155557555259,
                                   grepl("Austernbank",density5ao$BANK_TYPE) ~ 525))





#add density column
density5ao=density5ao%>%
  mutate(density=n_Oysters5a/sub_samp_area)



#add mean density to res_oys5
mean_Do5a=aggregate(density5ao$density,by=list(density5ao$id),FUN=mean)
names(mean_Do5a)[names(mean_Do5a)=="Group.1"]<-"id"
res_oys5<-merge(mean_Do5a,res_oys5,by="id",all = T)
names(res_oys5)[names(res_oys5)=="x"]<-"mean_Do5a"





#add sd density to res_oys5
sd_Do5a<-aggregate(density5ao$density,by=list(density5ao$id),FUN=sd)
names(sd_Do5a)[names(sd_Do5a)=="Group.1"]<-"id"
res_oys5<-merge(sd_Do5a,res_oys5,by="id",all = T)
names(res_oys5)[names(res_oys5)=="x"]<-"sd_Do5a"









###density for 5b samples###
set.seed(12)
g=group_by(oys5b,id,sample)
density5bo=summarise(g,length(sample))



names(density5bo)[names(density5bo) == 'length(sample)'] <- 'n_Oysters5b'


density5bo$bank <- as.vector(str_trim(str_extract(density5bo$id,"([:upper:][:upper:][:number:][:number:])")))
density5bo <- merge(x = density5bo, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)
require(dplyr)
density5bo <- density5bo%>%
  mutate(sub_samp_area = case_when(grepl("Miesmuschelbank",density5bo$BANK_TYPE) ~ 175.7155557555259,
                                   grepl("Austernbank",density5bo$BANK_TYPE) ~ 525))





#add density column
density5bo=density5bo%>%
  mutate(density=n_Oysters5b/sub_samp_area)




#add mean density to res_oys5
mean_Do5b=aggregate(density5bo$density,by=list(density5bo$id),FUN=mean)
names(mean_Do5b)[names(mean_Do5b)=="Group.1"]<-"id"
res_oys5<-merge(mean_Do5b,res_oys5,by="id",all = T)
names(res_oys5)[names(res_oys5)=="x"]<-"mean_Do5b"





#add sd density to res_oys5
sd_Do5b<-aggregate(density5bo$density,by=list(density5bo$id),FUN=sd)
names(sd_Do5b)[names(sd_Do5b)=="Group.1"]<-"id"
res_oys5<-merge(sd_Do5b,res_oys5,by="id",all = T)
names(res_oys5)[names(res_oys5)=="x"]<-"sd_Do5b"









###density for 5c samples###
set.seed(13)
g=group_by(oys5c,id,sample)
density5co=summarise(g,length(sample))



names(density5co)[names(density5co) == 'length(sample)'] <- 'n_Oysters5c'


density5co$bank <- as.vector(str_trim(str_extract(density5co$id,"([:upper:][:upper:][:number:][:number:])")))
density5co <- merge(x = density5co, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)
require(dplyr)
density5co <- density5co%>%
  mutate(sub_samp_area = case_when(grepl("Miesmuschelbank",density5co$BANK_TYPE) ~ 175.7155557555259,
                                   grepl("Austernbank",density5co$BANK_TYPE) ~ 525))





#add density column
density5co=density5co%>%
  mutate(density=n_Oysters5c/sub_samp_area)




#add mean density to res_oys5
mean_Do5c=aggregate(density5co$density,by=list(density5co$id),FUN=mean)
names(mean_Do5c)[names(mean_Do5c)=="Group.1"]<-"id"
res_oys5<-merge(mean_Do5c,res_oys5,by="id",all = T)
names(res_oys5)[names(res_oys5)=="x"]<-"mean_Do5c"

res_oys5<-merge(mean_Do12,res_oys5,by="id",all = T)
names(res_oys5)[names(res_oys5)=="x"]<-"mean_Do12"



#add sd density to res_oys5
sd_Do5c<-aggregate(density5co$density,by=list(density5co$id),FUN=sd)
names(sd_Do5c)[names(sd_Do5c)=="Group.1"]<-"id"
res_oys5<-merge(sd_Do5c,res_oys5,by="id",all = T)
names(res_oys5)[names(res_oys5)=="x"]<-"sd_Do5c"

res_oys5<-merge(sd_Do12,res_oys5,by="id",all = T)
names(res_oys5)[names(res_oys5)=="x"]<-"sd_Do12"



list_dfDO5<-list(mean_Do12,mean_Do5a,mean_Do5b,mean_Do5c)
mean_densityO5<-list_dfDO5%>%
  reduce(inner_join,by="id")

colnames(mean_densityO5)<-c("id","meanD12","meanD5a","meanD5b","meanD5c")



melt5do<-melt(mean_densityO5,id="id")
colnames(melt5do)<-c("id","replicates","mean_density")

melt5do%>% group_by(replicates) %>% shapiro_test(mean_density) # non normally distributed

melt5do %>% wilcox_test(mean_density~replicates,paired = T) #significant













#### Plot Mussel Length ####


### 8 samples ###
meltLm<-melt(res_lemu12,id.vars =("id"),measure.vars = c("mean_L12","mean_L8a","mean_L8b","mean_L8c"))
meltLm_sd<-melt(res_lemu12,id.vars =("id"),measure.vars = c("sd_L12","sd_L8a","sd_L8b","sd_L8c"))
meltLm$variable=gsub("mean_","",as.character(meltLm$variable))
meltLm_sd$variable=gsub("sd_","",as.character(meltLm_sd$variable))

meltLenMuss=merge(meltLm,meltLm_sd,by=c("id","variable"))

colnames(meltLenMuss)<-c("id","sample","mean","sd")


Lm8<-ggplot(meltLenMuss[grepl("^LT",meltLenMuss$id),],aes(x=id,y=mean))+
  geom_point(position=position_dodge(0.5),aes(color=sample,shape=sample))+
  labs(y="Mean length",x="Bank and date",title= "Mussel mean length")+
  coord_flip()+
  geom_vline(xintercept=seq(1.5, length(unique(meltLenMuss$id))-0.5, 1), 
             lwd=1, colour="grey",linetype="longdash")+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd,color=sample),position=position_dodge(0.5),width=.2)



### 6 samples ###
meltLm6<-melt(res_lemu6,id.vars =("id"),measure.vars = c("mean_L12","mean_L6a","mean_L6b","mean_L6c"))
meltLm6_sd<-melt(res_lemu6,id.vars =("id"),measure.vars = c("sd_L12","sd_L6a","sd_L6b","sd_L6c"))
meltLm6$variable=gsub("mean_","",as.character(meltLm6$variable))
meltLm6_sd$variable=gsub("sd_","",as.character(meltLm6_sd$variable))

meltLenMuss6=merge(meltLm6,meltLm6_sd,by=c("id","variable"))

colnames(meltLenMuss6)<-c("id","sample","mean","sd")


Lm6<-ggplot(meltLenMuss6[grepl("^LT",meltLenMuss6$id),],aes(x=id,y=mean))+
  geom_point(position=position_dodge(0.5),aes(color=sample,shape=sample))+
  labs(y="Mean length",x="Bank and date")+
  coord_flip()+
  geom_vline(xintercept=seq(1.5, length(unique(meltLenMuss6$id))-0.5, 1), 
             lwd=1, colour="grey",linetype="longdash")+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd,color=sample),position=position_dodge(0.5),width=.2)




### 5 samples ###
meltLm5<-melt(res_lemu5,id.vars =("id"),measure.vars = c("mean_L12","mean_L5a","mean_L5b","mean_L5c"))
meltLm5_sd<-melt(res_lemu5,id.vars =("id"),measure.vars = c("sd_L12","sd_L5a","sd_L5b","sd_L5c"))
meltLm5$variable=gsub("mean_","",as.character(meltLm5$variable))
meltLm5_sd$variable=gsub("sd_","",as.character(meltLm5_sd$variable))

meltLenMuss5=merge(meltLm5,meltLm5_sd,by=c("id","variable"))

colnames(meltLenMuss5)<-c("id","sample","mean","sd")


Lm5<-ggplot(meltLenMuss5[grepl("^LT",meltLenMuss5$id),],aes(x=id,y=mean))+
  geom_point(position=position_dodge(0.5),aes(color=sample,shape=sample))+
  labs(y="Mean length",x="Bank and date")+
  coord_flip()+
  geom_vline(xintercept=seq(1.5, length(unique(meltLenMuss5$id))-0.5, 1), 
             lwd=1, colour="grey",linetype="longdash")+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd,color=sample),position=position_dodge(0.5),width=.2)


lay<-rbind(c(1,1,2,2),c(NA,3,3,NA))
grid.arrange(Lm8,Lm6, Lm5,layout_matrix=lay)







#### Plot Mussel Density ####


### 8 samples ###
meltDm<-melt(res_lemu12,id.vars =("id"),measure.vars = c("mean_D12","mean_D8a","mean_D8b","mean_D8c"))
meltDm_sd<-melt(res_lemu12,id.vars =("id"),measure.vars = c("sd_D12","sd_D8a","sd_D8b","sd_D8c"))
meltDm$variable=gsub("mean_","",as.character(meltDm$variable))
meltDm_sd$variable=gsub("sd_","",as.character(meltDm_sd$variable))

meltDensMuss=merge(meltDm,meltDm_sd,by=c("id","variable"))

colnames(meltDensMuss)<-c("id","sample","mean","sd")


Dm8<-ggplot(meltDensMuss[grepl("^LT",meltDensMuss$id),],aes(x=id,y=mean))+
  geom_point(position=position_dodge(0.5),aes(color=sample,shape=sample))+
  labs(y="Mean density",x="Bank and date",title= "Mussel mean density")+
  coord_flip()+
  geom_vline(xintercept=seq(1.5, length(unique(meltDensMuss$id))-0.5, 1), 
             lwd=1, colour="grey",linetype="longdash")+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd,color=sample),position=position_dodge(0.5),width=.2)



### 6 samples ###
meltDm6<-melt(res_lemu6,id.vars =("id"),measure.vars = c("mean_D12","mean_D6a","mean_D6b","mean_D6c"))
meltDm6_sd<-melt(res_lemu6,id.vars =("id"),measure.vars = c("sd_D12","sd_D6a","sd_D6b","sd_D6c"))
meltDm6$variable=gsub("mean_","",as.character(meltDm6$variable))
meltDm6_sd$variable=gsub("sd_","",as.character(meltDm6_sd$variable))

meltDensMuss6=merge(meltDm6,meltDm6_sd,by=c("id","variable"))

colnames(meltDensMuss6)<-c("id","sample","mean","sd")


Dm6<-ggplot(meltDensMuss6[grepl("^LT",meltDensMuss6$id),],aes(x=id,y=mean))+
  geom_point(position=position_dodge(0.5),aes(color=sample,shape=sample))+
  labs(y="Mean density",x="Bank and date")+
  coord_flip()+
  geom_vline(xintercept=seq(1.5, length(unique(meltDensMuss6$id))-0.5, 1), 
             lwd=1, colour="grey",linetype="longdash")+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd,color=sample),position=position_dodge(0.5),width=.2)

  


### 5 samples ###
meltDm5<-melt(res_lemu5,id.vars =("id"),measure.vars = c("mean_D12","mean_D5a","mean_D5b","mean_D5c"))
meltDm5_sd<-melt(res_lemu5,id.vars =("id"),measure.vars = c("sd_D12","sd_D5a","sd_D5b","sd_D5c"))
meltDm5$variable=gsub("mean_","",as.character(meltDm5$variable))
meltDm5_sd$variable=gsub("sd_","",as.character(meltDm5_sd$variable))

meltDensMuss5=merge(meltDm5,meltDm5_sd,by=c("id","variable"))

colnames(meltDensMuss5)<-c("id","sample","mean","sd")


Dm5<-ggplot(meltDensMuss5[grepl("^LT",meltDensMuss5$id),],aes(x=id,y=mean))+
  geom_point(position=position_dodge(0.5),aes(color=sample,shape=sample))+
  labs(y="Mean density",x="Bank and date")+
  coord_flip()+
  geom_vline(xintercept=seq(1.5, length(unique(meltDensMuss5$id))-0.5, 1), 
             lwd=1, colour="grey",linetype="longdash")+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd,color=sample),position=position_dodge(0.5),width=.2)

lay<-rbind(c(1,1,2,2),c(NA,3,3,NA))
grid.arrange(Dm8,Dm6, Dm5,layout_matrix=lay)










#### Plot Oyster Length ####


### 8 samples ###
meltLmo<-melt(res_oys12,id.vars =("id"),measure.vars = c("mean_Lo12","mean_Lo8a","mean_Lo8b","mean_Lo8c"))
meltLmo_sd<-melt(res_oys12,id.vars =("id"),measure.vars = c("sd_Lo12","sd_Lo8a","sd_Lo8b","sd_Lo8c"))
meltLmo$variable=gsub("mean_","",as.character(meltLmo$variable))
meltLmo_sd$variable=gsub("sd_","",as.character(meltLmo_sd$variable))

meltLenOys=merge(meltLmo,meltLmo_sd,by=c("id","variable"))

colnames(meltLenOys)<-c("id","sample","mean","sd")


Lmo8<-ggplot(meltLenOys[grepl("^LT",meltLenOys$id),],aes(x=id,y=mean))+
  geom_point(position=position_dodge(0.5),aes(color=sample,shape=sample))+
  labs(y="Mean length",x="Bank and date",title= "Oyster mean length")+
  coord_flip()+
  geom_vline(xintercept=seq(1.5, length(unique(meltLenOys$id))-0.5, 1), 
             lwd=1, colour="grey",linetype="longdash")+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd,color=sample),position=position_dodge(0.5),width=.2)



### 6 samples ###
meltLmo6<-melt(res_oys6,id.vars =("id"),measure.vars = c("mean_Lo12","mean_Lo6a","mean_Lo6b","mean_Lo6c"))
meltLmo6_sd<-melt(res_oys6,id.vars =("id"),measure.vars = c("sd_Lo12","sd_Lo6a","sd_Lo6b","sd_Lo6c"))
meltLmo6$variable=gsub("mean_","",as.character(meltLmo6$variable))
meltLmo6_sd$variable=gsub("sd_","",as.character(meltLmo6_sd$variable))

meltLenOys6=merge(meltLmo6,meltLmo6_sd,by=c("id","variable"))

colnames(meltLenOys6)<-c("id","sample","mean","sd")


Lmo6<-ggplot(meltLenOys6[grepl("^LT",meltLenOys6$id),],aes(x=id,y=mean))+
  geom_point(position=position_dodge(0.5),aes(color=sample,shape=sample))+
  labs(y="Mean length",x="Bank and date")+
  coord_flip()+
  geom_vline(xintercept=seq(1.5, length(unique(meltLenOys6$id))-0.5, 1), 
             lwd=1, colour="grey",linetype="longdash")+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd,color=sample),position=position_dodge(0.5),width=.2)




### 5 samples ###
meltLmo5<-melt(res_oys5,id.vars =("id"),measure.vars = c("mean_Lo12","mean_Lo5a","mean_Lo5b","mean_Lo5c"))
meltLmo5_sd<-melt(res_oys5,id.vars =("id"),measure.vars = c("sd_Lo12","sd_Lo5a","sd_Lo5b","sd_Lo5c"))
meltLmo5$variable=gsub("mean_","",as.character(meltLmo5$variable))
meltLmo5_sd$variable=gsub("sd_","",as.character(meltLmo5_sd$variable))

meltLenOys5=merge(meltLmo5,meltLmo5_sd,by=c("id","variable"))

colnames(meltLenOys5)<-c("id","sample","mean","sd")


Lmo5<-ggplot(meltLenOys5[grepl("^LT",meltLenOys5$id),],aes(x=id,y=mean))+
  geom_point(position=position_dodge(0.5),aes(color=sample,shape=sample))+
  labs(y="Mean length",x="Bank and date")+
  coord_flip()+
  geom_vline(xintercept=seq(1.5, length(unique(meltLenOys5$id))-0.5, 1), 
             lwd=1, colour="grey",linetype="longdash")+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd,color=sample),position=position_dodge(0.5),width=.2)


lay<-rbind(c(1,1,2,2),c(NA,3,3,NA))
grid.arrange(Lmo8,Lmo6, Lmo5,layout_matrix=lay)









#### Plot Oyster Density ####


### 8 samples ###
meltDmo<-melt(res_oys12,id.vars =("id"),measure.vars = c("mean_Do12","mean_Do8a","mean_Do8b","mean_Do8c"))
meltDmo_sd<-melt(res_oys12,id.vars =("id"),measure.vars = c("sd_Do12","sd_Do8a","sd_Do8b","sd_Do8c"))
meltDmo$variable=gsub("mean_","",as.character(meltDmo$variable))
meltDmo_sd$variable=gsub("sd_","",as.character(meltDmo_sd$variable))

meltDensOys=merge(meltDmo,meltDmo_sd,by=c("id","variable"))

colnames(meltDensOys)<-c("id","sample","mean","sd")


Dmo8<-ggplot(meltDensOys[grepl("^LT",meltDensOys$id),],aes(x=id,y=mean))+
  geom_point(position=position_dodge(0.5),aes(color=sample,shape=sample))+
  labs(y="Mean density",x="Bank and date",title= "Oyster mean density")+
  coord_flip()+
  geom_vline(xintercept=seq(1.5, length(unique(meltDensOys$id))-0.5, 1), 
             lwd=1, colour="grey",linetype="longdash")+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd,color=sample),position=position_dodge(0.5),width=.2)



### 6 samples ###
meltDmo6<-melt(res_oys6,id.vars =("id"),measure.vars = c("mean_Do12","mean_Do6a","mean_Do6b","mean_Do6c"))
meltDmo6_sd<-melt(res_oys6,id.vars =("id"),measure.vars = c("sd_Do12","sd_Do6a","sd_Do6b","sd_Do6c"))
meltDmo6$variable=gsub("mean_","",as.character(meltDmo6$variable))
meltDmo6_sd$variable=gsub("sd_","",as.character(meltDmo6_sd$variable))

meltDensOys6=merge(meltDmo6,meltDmo6_sd,by=c("id","variable"))

colnames(meltDensOys6)<-c("id","sample","mean","sd")


Dmo6<-ggplot(meltDensOys6[grepl("^LT",meltDensOys6$id),],aes(x=id,y=mean))+
  geom_point(position=position_dodge(0.5),aes(color=sample,shape=sample))+
  labs(y="Mean density",x="Bank and date")+
  coord_flip()+
  geom_vline(xintercept=seq(1.5, length(unique(meltDensOys6$id))-0.5, 1), 
             lwd=1, colour="grey",linetype="longdash")+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd,color=sample),position=position_dodge(0.5),width=.2)




### 5 samples ###
meltDmo5<-melt(res_oys5,id.vars =("id"),measure.vars = c("mean_Do12","mean_Do5a","mean_Do5b","mean_Do5c"))
meltDmo5_sd<-melt(res_oys5,id.vars =("id"),measure.vars = c("sd_Do12","sd_Do5a","sd_Do5b","sd_Do5c"))
meltDmo5$variable=gsub("mean_","",as.character(meltDmo5$variable))
meltDmo5_sd$variable=gsub("sd_","",as.character(meltDmo5_sd$variable))

meltDensOys5=merge(meltDmo5,meltDmo5_sd,by=c("id","variable"))

colnames(meltDensOys5)<-c("id","sample","mean","sd")


Dmo5<-ggplot(meltDensOys5[grepl("^LT",meltDensOys5$id),],aes(x=id,y=mean))+
  geom_point(position=position_dodge(0.5),aes(color=sample,shape=sample))+
  labs(y="Mean density",x="Bank and date")+
  coord_flip()+
  geom_vline(xintercept=seq(1.5, length(unique(meltDensOys5$id))-0.5, 1), 
             lwd=1, colour="grey",linetype="longdash")+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd,color=sample),position=position_dodge(0.5),width=.2)

lay<-rbind(c(1,1,2,2),c(NA,3,3,NA))
grid.arrange(Dmo8,Dmo6, Dmo5,layout_matrix=lay)















#### Cluster plots on density and length ####

#Mussel

colnames(mean12)<-c("id","length")
colnames(mean_D12)<-c("id","density")
meanDL12<-merge(mean12,mean_D12,by="id",all=T)

meanDL12<-meanDL12 %>% remove_rownames %>% column_to_rownames(var="id")

kms<-kmeans(meanDL12$length,9,nstart = 100)


fviz_nbclust(meanDL12,kmeans,method = "silhouette")
fviz_cluster(kms,data = meanDL12)






#Oyster

colnames(mean12o)<-c("id","length")
colnames(mean_Do12)<-c("id","density")
meanDL12o<-merge(mean12o,mean_Do12,by="id",all=T)

meanDL12o<-meanDL12o %>% remove_rownames %>% column_to_rownames(var="id")

kms<-kmeans(meanDL12o$length,2,nstart = 100)


fviz_nbclust(meanDL12o,kmeans,method = "silhouette")
fviz_cluster(kms,data = meanDL12o)










#Miesmuchelbank

mean12$bank <- as.vector(str_trim(str_extract(mean12$id,"([:upper:][:upper:][:number:][:number:])")))
mean12 <- merge(x = mean12, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)
mean_D12$bank <- as.vector(str_trim(str_extract(mean_D12$id,"([:upper:][:upper:][:number:][:number:])")))
mean_D12 <- merge(x = mean_D12, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)




meanmmbD12<-mean_D12[mean_D12$BANK_TYPE=="Miesmuschelbank",]
meanmmbL12<-mean12[mean12$BANK_TYPE=="Miesmuschelbank",]

meanmmb12<-merge(meanmmbL12,meanmmbD12,by="id",all=T)

meanmmb12<-meanmmb12 %>% remove_rownames %>% column_to_rownames(var="id")

meanmmb12<-meanmmb12[,-c(1,3,4,6)]


kms<-kmeans(meanmmb12$length,9,nstart = 100)


fviz_nbclust(meanmmb12,kmeans,method = "silhouette")
fviz_cluster(kms,data = meanmmb12)






#Austernbank

meanausD12<-mean_D12[mean_D12$BANK_TYPE=="Austernbank",]
meanausL12<-mean12[mean12$BANK_TYPE=="Austernbank",]

meanaus12<-merge(meanausL12,meanausD12,by="id",all=T)

meanaus12<-meanaus12 %>% remove_rownames %>% column_to_rownames(var="id")

meanaus12<-meanaus12[,-c(1,3,4,6)]


kms<-kmeans(meanaus12$length,2,nstart = 100)

fviz_nbclust(meanaus12,kmeans,method = "silhouette")
fviz_cluster(kms,data = meanaus12)







#### 12 and 8 samples ####
### 12 sample
length_2022 <- length_2022 %>%
  mutate(size_class = case_when(species == "mussel" &  length >19.99 & length <= 29.99 ~ "m20",
                                species == "mussel" &length >29.99 & length <= 39.99 ~ "m30",
                                species == "mussel" & length >39.99 & length <= 49.99 ~ "m40",
                                species == "mussel" &length >49.99 & length <= 59.99 ~ "m50",
                                species == "mussel" & length >59.99 & length <= 69.99 ~ "m60",
                                species == "mussel" & length >69.99 & length <= 79.99 ~ "m70",
                                species == "oyster" &  length >24.99 & length <= 49.99 ~ "o25",
                                species == "oyster" &length >49.99 & length <= 74.99 ~ "o50",
                                species == "oyster" & length >74.99 & length <= 99.99 ~ "o75",
                                species == "oyster" & length >99.99 & length <= 124.99 ~ "o100",
                                species == "oyster" & length >124.99 & length <= 149.99 ~ "o125",
                                species == "oyster" & length >149.99 & length <= 174.99 ~ "o150"))


distrib12<-length_2022 %>% group_by(id, size_class) %>% summarise(length(length))
distrib12<-na.omit(distrib12)



distrib12$bank <- as.vector(str_trim(str_extract(distrib12$id,"([:upper:][:upper:][:number:][:number:])")))
distrib12 <- merge(x = distrib12, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)


names(distrib12)[names(distrib12) == 'length(length)'] <- 'N12'


#rearrange size_class order 
distrib12$size_class<-factor(distrib12$size_class,
                             levels = c("m20", "m30", "m40",  "m50" , "m60",  "m70","o25","o50" , "o75","o100" ,"o125", "o150"))


### 8a sample

length_2022_8a=length_2022 %>% filter(sample %in% sample(unique(sample),8))

length_2022_8a <- length_2022_8a %>%
  mutate(size_class = case_when(species == "mussel" &  length >19.99 & length <= 29.99 ~ "m20",
                                species == "mussel" &length >29.99 & length <= 39.99 ~ "m30",
                                species == "mussel" & length >39.99 & length <= 49.99 ~ "m40",
                                species == "mussel" &length >49.99 & length <= 59.99 ~ "m50",
                                species == "mussel" & length >59.99 & length <= 69.99 ~ "m60",
                                species == "mussel" & length >69.99 & length <= 79.99 ~ "m70",
                                species == "oyster" &  length >24.99 & length <= 49.99 ~ "o25",
                                species == "oyster" &length >49.99 & length <= 74.99 ~ "o50",
                                species == "oyster" & length >74.99 & length <= 99.99 ~ "o75",
                                species == "oyster" & length >99.99 & length <= 124.99 ~ "o100",
                                species == "oyster" & length >124.99 & length <= 149.99 ~ "o125",
                                species == "oyster" & length >149.99 & length <= 174.99 ~ "o150"))

distrib8a<-length_2022_8a %>% group_by(id, size_class) %>% summarise(length(length))
distrib8a<-na.omit(distrib8a)


distrib8a$bank <- as.vector(str_trim(str_extract(distrib8a$id,"([:upper:][:upper:][:number:][:number:])")))
distrib8a <- merge(x = distrib8a, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)


names(distrib8a)[names(distrib8a) == 'length(length)'] <- 'N8a'


distrib8a$size_class<-factor(distrib8a$size_class,
                             levels = c("m20", "m30", "m40",  "m50" , "m60",  "m70","o25","o50" , "o75","o100" ,"o125", "o150"))







### 8b sample


length_2022_8b=length_2022 %>% filter(sample %in% sample(unique(sample),8))

length_2022_8b <- length_2022_8b %>%
  mutate(size_class = case_when(species == "mussel" &  length >19.99 & length <= 29.99 ~ "m20",
                                species == "mussel" &length >29.99 & length <= 39.99 ~ "m30",
                                species == "mussel" & length >39.99 & length <= 49.99 ~ "m40",
                                species == "mussel" &length >49.99 & length <= 59.99 ~ "m50",
                                species == "mussel" & length >59.99 & length <= 69.99 ~ "m60",
                                species == "mussel" & length >69.99 & length <= 79.99 ~ "m70",
                                species == "oyster" &  length >24.99 & length <= 49.99 ~ "o25",
                                species == "oyster" &length >49.99 & length <= 74.99 ~ "o50",
                                species == "oyster" & length >74.99 & length <= 99.99 ~ "o75",
                                species == "oyster" & length >99.99 & length <= 124.99 ~ "o100",
                                species == "oyster" & length >124.99 & length <= 149.99 ~ "o125",
                                species == "oyster" & length >149.99 & length <= 174.99 ~ "o150"))

distrib8b<-length_2022_8b %>% group_by(id, size_class) %>% summarise(length(length))
distrib8b<-na.omit(distrib8b)


distrib8b$bank <- as.vector(str_trim(str_extract(distrib8b$id,"([:upper:][:upper:][:number:][:number:])")))
distrib8b <- merge(x = distrib8b, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)


names(distrib8b)[names(distrib8b) == 'length(length)'] <- 'N8b'




distrib8b$size_class<-factor(distrib8b$size_class,
                             levels = c("m20", "m30", "m40",  "m50" , "m60",  "m70","o25","o50" , "o75","o100" ,"o125", "o150"))





### 8c sample
length_2022_8c=length_2022 %>% filter(sample %in% sample(unique(sample),8))

length_2022_8c <- length_2022_8c %>%
  mutate(size_class = case_when(species == "mussel" &  length >19.99 & length <= 29.99 ~ "m20",
                                species == "mussel" &length >29.99 & length <= 39.99 ~ "m30",
                                species == "mussel" & length >39.99 & length <= 49.99 ~ "m40",
                                species == "mussel" &length >49.99 & length <= 59.99 ~ "m50",
                                species == "mussel" & length >59.99 & length <= 69.99 ~ "m60",
                                species == "mussel" & length >69.99 & length <= 79.99 ~ "m70",
                                species == "oyster" &  length >24.99 & length <= 49.99 ~ "o25",
                                species == "oyster" &length >49.99 & length <= 74.99 ~ "o50",
                                species == "oyster" & length >74.99 & length <= 99.99 ~ "o75",
                                species == "oyster" & length >99.99 & length <= 124.99 ~ "o100",
                                species == "oyster" & length >124.99 & length <= 149.99 ~ "o125",
                                species == "oyster" & length >149.99 & length <= 174.99 ~ "o150"))

distrib8c<-length_2022_8c %>% group_by(id, size_class) %>% summarise(length(length))
distrib8c<-na.omit(distrib8c)


distrib8c$bank <- as.vector(str_trim(str_extract(distrib8c$id,"([:upper:][:upper:][:number:][:number:])")))
distrib8c <- merge(x = distrib8c, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)


names(distrib8c)[names(distrib8c) == 'length(length)'] <- 'N8c'



distrib8c$size_class<-factor(distrib8c$size_class,
                             levels = c("m20", "m30", "m40",  "m50" , "m60",  "m70","o25","o50" , "o75","o100" ,"o125", "o150"))














require(reshape2)
distrib12_cast <- dcast(distrib12, bank+id+BANK_TYPE ~ size_class, value.var = "N12")





distrib12 <- distrib12 %>%
  complete(id,size_class, fill = list(N12 = 0))
table(distrib12$id,distrib12$size_class)

# order df (relevant for the subsequent "filling" process)
distrib12 <- distrib12[order(distrib12$id,distrib12$size_class),]

# missing information in columns PROBE_NUMBER, METHOD_ID will be filled in the direction down (grouped by column SAMPLE_ID)
#fill is part of tidyr
distrib12 <- distrib12 %>% 
  group_by(id) %>% 
  fill(bank, BANK_TYPE) 






distrib8a <- distrib8a %>%
  complete(id,size_class, fill = list(N8a = 0))
table(distrib8a$id,distrib8a$size_class)

# order df (relevant for the subsequent "filling" process)
distrib8a <- distrib8a[order(distrib8a$id,distrib8a$size_class),]

# missing information in columns PROBE_NUMBER, METHOD_ID will be filled in the direction down (grouped by column SAMPLE_ID)
#fill is part of tidyr
distrib8a <- distrib8a %>% 
  group_by(id) %>% 
  fill(bank, BANK_TYPE) 






distrib8b <- distrib8b %>%
  complete(id,size_class, fill = list(N8b = 0))
table(distrib8b$id,distrib8b$size_class)

# order df (relevant for the subsequent "filling" process)
distrib8b <- distrib8b[order(distrib8b$id,distrib8b$size_class),]

# missing information in columns PROBE_NUMBER, METHOD_ID will be filled in the direction down (grouped by column SAMPLE_ID)
#fill is part of tidyr
distrib8b <- distrib8b %>% 
  group_by(id) %>% 
  fill(bank, BANK_TYPE) 






distrib8c <- distrib8c %>%
  complete(id,size_class, fill = list(N8c = 0))
table(distrib8c$id,distrib8c$size_class)

# order df (relevant for the subsequent "filling" process)
distrib8c <- distrib8c[order(distrib8c$id,distrib8c$size_class),]

# missing information in columns PROBE_NUMBER, METHOD_ID will be filled in the direction down (grouped by column SAMPLE_ID)
#fill is part of tidyr
distrib8c <- distrib8c %>% 
  group_by(id) %>% 
  fill(bank, BANK_TYPE) 













distrib12_Nclass<-read.csv("P:/Hugo/Data/data_chris/distrib_Nclass.csv")
distrib12_Nclass3<-read.csv("P:/Hugo/Data/data_chris/distrib_Nclass3.csv")

distrib8a_Nclass3<-read.csv("P:/Hugo/Data/data_chris/distrib8a_Nclass3.csv")
distrib8b_Nclass3<-read.csv("P:/Hugo/Data/data_chris/distrib8b_Nclass3.csv")
distrib8c_Nclass3<-read.csv("P:/Hugo/Data/data_chris/distrib8c_Nclass3.csv")


colnames(distrib12_Nclass3) <- c("id","state_class","n_class")





library(dplyr)
compilation <- left_join(distrib12_Nclass3, distrib8a_Nclass3, by=c("id","state_class")) %>%
  left_join(., distrib8b_Nclass3,  by=c("id","state_class"))%>%
  left_join(., distrib8c_Nclass3,  by=c("id","state_class")) 

colnames(compilation) <- c("id","state_class","n_class12","n_class8a","n_class8b","n_class8c")






shapiro.test(compilation$n_class12) #non normal
shapiro.test(compilation$n_class8a) #non normal
shapiro.test(compilation$n_class8b) #non normal
shapiro.test(compilation$n_class8c) #non normal




compilgood<-compilation[compilation$state_class=="Good",]


wilcox.exact(compilgood$n_class12,compilgood$n_class8a,paired = T) #***
wilcox.exact(compilgood$n_class12,compilgood$n_class8b,paired = T) #***
wilcox.exact(compilgood$n_class12,compilgood$n_class8c,paired = T) #***




ggplot(compilgood)+
  geom_boxplot(aes(x="12",y=n_class12))+
  geom_boxplot(aes(x="8a",y=n_class8a))+
  geom_boxplot(aes(x="8b",y=n_class8b))+
  geom_boxplot(aes(x="8c",y=n_class8c))+
  labs(x="Number of samples",y="Number of good classes")


compilmedium<-compilation[compilation$state_class=="Medium",]


wilcox.exact(compilmedium$n_class12,compilmedium$n_class8a,paired = T) #ns
wilcox.exact(compilmedium$n_class12,compilmedium$n_class8b,paired = T) #ns
wilcox.exact(compilmedium$n_class12,compilmedium$n_class8c,paired = T) #ns





ggplot(compilmedium)+
  geom_boxplot(aes(x="12",y=n_class12))+
  geom_boxplot(aes(x="8a",y=n_class8a))+
  geom_boxplot(aes(x="8b",y=n_class8b))+
  geom_boxplot(aes(x="8c",y=n_class8c))+
  labs(x="Number of samples",y="Number of medium classes")




compilbad<-compilation[compilation$state_class=="Bad",]


wilcox.exact(compilbad$n_class12,compilbad$n_class8a,paired = T) #***
wilcox.exact(compilbad$n_class12,compilbad$n_class8b,paired = T) #**
wilcox.exact(compilbad$n_class12,compilbad$n_class8c,paired = T) #**




ggplot(compilbad)+
  geom_boxplot(aes(x="12",y=n_class12))+
  geom_boxplot(aes(x="8a",y=n_class8a))+
  geom_boxplot(aes(x="8b",y=n_class8b))+
  geom_boxplot(aes(x="8c",y=n_class8c))+
  labs(x="Number of samples",y="Number of bad classes")






kruskal_test(compilation,compilation$n_class12~compilation$state_class)
dunn_test(compilation,n_class12~state_class)




kruskal_test(compilation,compilation$n_class8a~compilation$state_class)
dunn_test(compilation,n_class8a~state_class)



kruskal_test(compilation,compilation$n_class8b~compilation$state_class)
dunn_test(compilation,n_class8b~state_class)



kruskal_test(compilation,compilation$n_class8c~compilation$state_class)
dunn_test(compilation,n_class8c~state_class)




require(reshape2)
melted_comp <- melt(compilation, id.vars = c("id","state_class"))


ggplot(melted_comp)+
  geom_bar(aes(state_class,value,fill=state_class),stat = "identity")+
  geom_text(aes(x="Bad",y=255,label = "ns"))+
  geom_text(aes(x="Good",y=240,label = "ns"))+
  geom_text(aes(x="Medium",y=45,label = "****"))+
  scale_fill_manual(values = c("red3","chartreuse3","orange3"))+
  facet_grid(.~variable)



kruskal_test(melted_comp,value~variable)











#### 12 and 6 samples ####


### 6a sample

length_2022_6a=length_2022 %>% filter(sample %in% sample(unique(sample),6))

length_2022_6a <- length_2022_6a %>%
  mutate(size_class = case_when(species == "mussel" &  length >19.99 & length <= 29.99 ~ "m20",
                                species == "mussel" &length >29.99 & length <= 39.99 ~ "m30",
                                species == "mussel" & length >39.99 & length <= 49.99 ~ "m40",
                                species == "mussel" &length >49.99 & length <= 59.99 ~ "m50",
                                species == "mussel" & length >59.99 & length <= 69.99 ~ "m60",
                                species == "mussel" & length >69.99 & length <= 79.99 ~ "m70",
                                species == "oyster" &  length >24.99 & length <= 49.99 ~ "o25",
                                species == "oyster" &length >49.99 & length <= 74.99 ~ "o50",
                                species == "oyster" & length >74.99 & length <= 99.99 ~ "o75",
                                species == "oyster" & length >99.99 & length <= 124.99 ~ "o100",
                                species == "oyster" & length >124.99 & length <= 149.99 ~ "o125",
                                species == "oyster" & length >149.99 & length <= 174.99 ~ "o150"))

distrib6a<-length_2022_6a %>% group_by(id, size_class) %>% summarise(length(length))
distrib6a<-na.omit(distrib6a)


distrib6a$bank <- as.vector(str_trim(str_extract(distrib6a$id,"([:upper:][:upper:][:number:][:number:])")))
distrib6a <- merge(x = distrib6a, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)


names(distrib6a)[names(distrib6a) == 'length(length)'] <- 'N6a'


distrib6a$size_class<-factor(distrib6a$size_class,
                             levels = c("m20", "m30", "m40",  "m50" , "m60",  "m70","o25","o50" , "o75","o100" ,"o125", "o150"))







### 6b sample


length_2022_6b=length_2022 %>% filter(sample %in% sample(unique(sample),6))

length_2022_6b <- length_2022_6b %>%
  mutate(size_class = case_when(species == "mussel" &  length >19.99 & length <= 29.99 ~ "m20",
                                species == "mussel" &length >29.99 & length <= 39.99 ~ "m30",
                                species == "mussel" & length >39.99 & length <= 49.99 ~ "m40",
                                species == "mussel" &length >49.99 & length <= 59.99 ~ "m50",
                                species == "mussel" & length >59.99 & length <= 69.99 ~ "m60",
                                species == "mussel" & length >69.99 & length <= 79.99 ~ "m70",
                                species == "oyster" &  length >24.99 & length <= 49.99 ~ "o25",
                                species == "oyster" &length >49.99 & length <= 74.99 ~ "o50",
                                species == "oyster" & length >74.99 & length <= 99.99 ~ "o75",
                                species == "oyster" & length >99.99 & length <= 124.99 ~ "o100",
                                species == "oyster" & length >124.99 & length <= 149.99 ~ "o125",
                                species == "oyster" & length >149.99 & length <= 174.99 ~ "o150"))

distrib6b<-length_2022_6b %>% group_by(id, size_class) %>% summarise(length(length))
distrib6b<-na.omit(distrib6b)


distrib6b$bank <- as.vector(str_trim(str_extract(distrib6b$id,"([:upper:][:upper:][:number:][:number:])")))
distrib6b <- merge(x = distrib6b, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)


names(distrib6b)[names(distrib6b) == 'length(length)'] <- 'N6b'




distrib6b$size_class<-factor(distrib6b$size_class,
                             levels = c("m20", "m30", "m40",  "m50" , "m60",  "m70","o25","o50" , "o75","o100" ,"o125", "o150"))





### 6c sample
length_2022_6c=length_2022 %>% filter(sample %in% sample(unique(sample),6))

length_2022_6c <- length_2022_6c %>%
  mutate(size_class = case_when(species == "mussel" &  length >19.99 & length <= 29.99 ~ "m20",
                                species == "mussel" &length >29.99 & length <= 39.99 ~ "m30",
                                species == "mussel" & length >39.99 & length <= 49.99 ~ "m40",
                                species == "mussel" &length >49.99 & length <= 59.99 ~ "m50",
                                species == "mussel" & length >59.99 & length <= 69.99 ~ "m60",
                                species == "mussel" & length >69.99 & length <= 79.99 ~ "m70",
                                species == "oyster" &  length >24.99 & length <= 49.99 ~ "o25",
                                species == "oyster" &length >49.99 & length <= 74.99 ~ "o50",
                                species == "oyster" & length >74.99 & length <= 99.99 ~ "o75",
                                species == "oyster" & length >99.99 & length <= 124.99 ~ "o100",
                                species == "oyster" & length >124.99 & length <= 149.99 ~ "o125",
                                species == "oyster" & length >149.99 & length <= 174.99 ~ "o150"))

distrib6c<-length_2022_6c %>% group_by(id, size_class) %>% summarise(length(length))
distrib6c<-na.omit(distrib6c)


distrib6c$bank <- as.vector(str_trim(str_extract(distrib6c$id,"([:upper:][:upper:][:number:][:number:])")))
distrib6c <- merge(x = distrib6c, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)


names(distrib6c)[names(distrib6c) == 'length(length)'] <- 'N6c'



distrib6c$size_class<-factor(distrib6c$size_class,
                             levels = c("m20", "m30", "m40",  "m50" , "m60",  "m70","o25","o50" , "o75","o100" ,"o125", "o150"))










require(reshape2)


distrib6a <- distrib6a %>%
  complete(id,size_class, fill = list(N6a = 0))
table(distrib6a$id,distrib6a$size_class)

# order df (relevant for the subsequent "filling" process)
distrib6a <- distrib6a[order(distrib6a$id,distrib6a$size_class),]

# missing information in columns PROBE_NUMBER, METHOD_ID will be filled in the direction down (grouped by column SAMPLE_ID)
#fill is part of tidyr
distrib6a <- distrib6a %>% 
  group_by(id) %>% 
  fill(bank, BANK_TYPE) 






distrib6b <- distrib6b %>%
  complete(id,size_class, fill = list(N6b = 0))
table(distrib6b$id,distrib6b$size_class)

# order df (relevant for the subsequent "filling" process)
distrib6b <- distrib6b[order(distrib6b$id,distrib6b$size_class),]

# missing information in columns PROBE_NUMBER, METHOD_ID will be filled in the direction down (grouped by column SAMPLE_ID)
#fill is part of tidyr
distrib6b <- distrib6b %>% 
  group_by(id) %>% 
  fill(bank, BANK_TYPE) 






distrib6c <- distrib6c %>%
  complete(id,size_class, fill = list(N6c = 0))
table(distrib6c$id,distrib6c$size_class)

# order df (relevant for the subsequent "filling" process)
distrib6c <- distrib6c[order(distrib6c$id,distrib6c$size_class),]

# missing information in columns PROBE_NUMBER, METHOD_ID will be filled in the direction down (grouped by column SAMPLE_ID)
#fill is part of tidyr
distrib6c <- distrib6c %>% 
  group_by(id) %>% 
  fill(bank, BANK_TYPE) 






distrib6a_Nclass <- distrib6a %>%
  mutate(state_class = case_when(N6a <= 4 ~ "Bad",
                                 N6a <= 9 & N6a > 4 ~ "Medium",
                                 N6a >= 10 ~ "Good"))

distrib6b_Nclass <- distrib6b %>%
  mutate(state_class = case_when(N6b <= 4 ~ "Bad",
                                 N6b <= 9 & N6b > 4 ~ "Medium",
                                 N6b >= 10 ~ "Good"))

distrib6c_Nclass <- distrib6c %>%
  mutate(state_class = case_when(N6c <= 4 ~ "Bad",
                                 N6c <= 9 & N6c > 4 ~ "Medium",
                                 N6c >= 10 ~ "Good"))







distrib6a_Nclass3<-distrib6a_Nclass %>%
  group_by(id)%>%
  summarise(Good=sum(state_class=="Good"))

distrib6a_Nclass3.2<-distrib6a_Nclass %>%
  group_by(id)%>%
  summarise(Medium=sum(state_class=="Medium"))


distrib6a_Nclass3.3<-distrib6a_Nclass %>%
  group_by(id)%>%
  summarise(Bad=sum(state_class=="Bad"))

distrib6a_Nclass3 <- merge(x = distrib6a_Nclass3, y = distrib6a_Nclass3.2, by = "id", all.x = TRUE)
distrib6a_Nclass3 <- merge(x = distrib6a_Nclass3, y = distrib6a_Nclass3.3, by = "id", all.x = TRUE)

distrib6a_Nclass3<-melt(distrib6a_Nclass3, id.vars = c("id"))





distrib6b_Nclass3<-distrib6b_Nclass %>%
  group_by(id)%>%
  summarise(Good=sum(state_class=="Good"))

distrib6b_Nclass3.2<-distrib6b_Nclass %>%
  group_by(id)%>%
  summarise(Medium=sum(state_class=="Medium"))


distrib6b_Nclass3.3<-distrib6b_Nclass %>%
  group_by(id)%>%
  summarise(Bad=sum(state_class=="Bad"))

distrib6b_Nclass3 <- merge(x = distrib6b_Nclass3, y = distrib6b_Nclass3.2, by = "id", all.x = TRUE)
distrib6b_Nclass3 <- merge(x = distrib6b_Nclass3, y = distrib6b_Nclass3.3, by = "id", all.x = TRUE)

distrib6b_Nclass3<-melt(distrib6b_Nclass3, id.vars = c("id"))





distrib6c_Nclass3<-distrib6c_Nclass %>%
  group_by(id)%>%
  summarise(Good=sum(state_class=="Good"))

distrib6c_Nclass3.2<-distrib6c_Nclass %>%
  group_by(id)%>%
  summarise(Medium=sum(state_class=="Medium"))


distrib6c_Nclass3.3<-distrib6c_Nclass %>%
  group_by(id)%>%
  summarise(Bad=sum(state_class=="Bad"))

distrib6c_Nclass3 <- merge(x = distrib6c_Nclass3, y = distrib6c_Nclass3.2, by = "id", all.x = TRUE)
distrib6c_Nclass3 <- merge(x = distrib6c_Nclass3, y = distrib6c_Nclass3.3, by = "id", all.x = TRUE)
distrib6c_Nclass3<-melt(distrib6c_Nclass3, id.vars = c("id"))


colnames(distrib6a_Nclass3) <- c("id","state_class","n_class")
colnames(distrib6b_Nclass3) <- c("id","state_class","n_class")
colnames(distrib6c_Nclass3) <- c("id","state_class","n_class")





library(dplyr)
compilation2 <- left_join(distrib12_Nclass3, distrib6a_Nclass3, by=c("id","state_class")) %>%
  left_join(., distrib6b_Nclass3,  by=c("id","state_class"))%>%
  left_join(., distrib6c_Nclass3,  by=c("id","state_class")) 

colnames(compilation2) <- c("id","state_class","n_class12","n_class6a","n_class6b","n_class6c")






shapiro.test(compilation2$n_class12) #non normal
shapiro.test(compilation2$n_class6a) #non normal
shapiro.test(compilation2$n_class6b) #non normal
shapiro.test(compilation2$n_class6c) #non normal




compilgood2<-compilation2[compilation2$state_class=="Good",]


wilcox.exact(compilgood2$n_class12,compilgood2$n_class6a,paired = T) #***
wilcox.exact(compilgood2$n_class12,compilgood2$n_class6b,paired = T) #***
wilcox.exact(compilgood2$n_class12,compilgood2$n_class6c,paired = T) #***




ggplot(compilgood2)+
  geom_boxplot(aes(x="12",y=n_class12))+
  geom_boxplot(aes(x="6a",y=n_class6a))+
  geom_boxplot(aes(x="6b",y=n_class6b))+
  geom_boxplot(aes(x="6c",y=n_class6c))+
  labs(x="Number of samples",y="Number of good classes")


compilmedium2<-compilation2[compilation2$state_class=="Medium",]


wilcox.exact(compilmedium2$n_class12,compilmedium2$n_class6a,paired = T) #ns
wilcox.exact(compilmedium2$n_class12,compilmedium2$n_class6b,paired = T) #ns
wilcox.exact(compilmedium2$n_class12,compilmedium2$n_class6c,paired = T) #ns





ggplot(compilmedium2)+
  geom_boxplot(aes(x="12",y=n_class12))+
  geom_boxplot(aes(x="6a",y=n_class6a))+
  geom_boxplot(aes(x="6b",y=n_class6b))+
  geom_boxplot(aes(x="6c",y=n_class6c))+
  labs(x="Number of samples",y="Number of medium classes")




compilbad2<-compilation2[compilation2$state_class=="Bad",]


wilcox.exact(compilbad2$n_class12,compilbad2$n_class6a,paired = T) #***
wilcox.exact(compilbad2$n_class12,compilbad2$n_class6b,paired = T) #**
wilcox.exact(compilbad2$n_class12,compilbad2$n_class6c,paired = T) #**




ggplot(compilbad2)+
  geom_boxplot(aes(x="12",y=n_class12))+
  geom_boxplot(aes(x="6a",y=n_class6a))+
  geom_boxplot(aes(x="6b",y=n_class6b))+
  geom_boxplot(aes(x="6c",y=n_class6c))+
  labs(x="Number of samples",y="Number of bad classes")






kruskal_test(compilation2,compilation2$n_class12~compilation2$state_class)
dunn_test(compilation2,n_class12~state_class)




kruskal_test(compilation2,compilation2$n_class6a~compilation2$state_class)
dunn_test(compilation2,n_class6a~state_class)



kruskal_test(compilation2,compilation2$n_class6b~compilation2$state_class)
dunn_test(compilation2,n_class6b~state_class)



kruskal_test(compilation2,compilation2$n_class6c~compilation2$state_class)
dunn_test(compilation2,n_class6c~state_class)




require(reshape2)
melted_comp2 <- melt(compilation2, id.vars = c("id","state_class"))


ggplot(melted_comp2)+
  geom_bar(aes(state_class,value,fill=state_class),stat = "identity")+
  geom_text(aes(x="Bad",y=255,label = "ns"))+
  geom_text(aes(x="Good",y=240,label = "ns"))+
  geom_text(aes(x="Medium",y=45,label = "****"))+
  scale_fill_manual(values = c("red3","chartreuse3","orange3"))+
  facet_grid(.~variable)



kruskal_test(melted_comp2,value~variable)












##### Reduce sample size #####


#### MUSSEL & OYSTERBANK ####
lemu$id_sample<- paste(lemu$id,lemu$sample)





### Add BANK_TYPE column
lemu$bank <- as.vector(str_trim(str_extract(lemu$id,"([:upper:][:upper:][:number:][:number:])")))
lemu <- merge(x = lemu, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)


### Subset with only Austernbank
lemu<-lemu[lemu$BANK_TYPE=="Austernbank",]








##### 64% random samples (3 times) #####
set.seed(1)

lemured<-lemu %>% 
  group_by(id_sample) %>%
  sample_frac(.64)



lemured2<-lemu %>% 
  group_by(id_sample) %>%
  sample_frac(.64)

lemured2b<-lemu %>% 
  group_by(id_sample) %>%
  sample_frac(.64)



### Mean length comparison
sd100<-aggregate(lemu$length,by=list(lemu$id),FUN=sd)

mean100=aggregate(lemu$length,by=list(lemu$id),FUN=mean)
res_lemu100<-merge(mean100,sd100,by="Group.1",all = T)





## First try
sd64_1<-aggregate(lemured$length,by=list(lemured$id),FUN=sd)

mean64_1=aggregate(lemured$length,by=list(lemured$id),FUN=mean)
res_lemu64_1<-merge(mean64_1,sd64_1,by="Group.1",all = T)



## Second try
sd64_2<-aggregate(lemured2$length,by=list(lemured2$id),FUN=sd)

mean64_2=aggregate(lemured2$length,by=list(lemured2$id),FUN=mean)
res_lemu64_2<-merge(mean64_2,sd64_2,by="Group.1",all = T)


##Third try
sd64_3<-aggregate(lemured2b$length,by=list(lemured2b$id),FUN=sd)

mean64_3=aggregate(lemured2b$length,by=list(lemured2b$id),FUN=mean)
res_lemu64_3<-merge(mean64_3,sd64_3,by="Group.1",all = T)



list_df64<-list(mean64_1,mean64_2,mean64_3,mean100)
mean_length64<-list_df64%>%
  reduce(inner_join,by="Group.1")

colnames(mean_length64)<-c("id","mean64_1","mean64_2","mean64_3","mean100")




melt64<-melt(mean_length64,id="id")
colnames(melt64)<-c("id","replicates","length")

melt64 %>% group_by(replicates) %>% shapiro_test(length) # normally distributed

melt64 %>% t_test(length~replicates,paired = T) # no differences between replicates







###### 50% random samples (3 times) #####
set.seed(1)
lemured3<-lemu %>% 
  group_by(id_sample) %>%
  sample_frac(.5)



lemured4<-lemu %>% 
  group_by(id_sample) %>%
  sample_frac(.5)


lemured4b<-lemu %>% 
  group_by(id_sample) %>%
  sample_frac(.5)




### Mean length comparison


## First try
sd50_1<-aggregate(lemured3$length,by=list(lemured3$id),FUN=sd)

mean50_1=aggregate(lemured3$length,by=list(lemured3$id),FUN=mean)
res_lemu50_1<-merge(mean50_1,sd50_1,by="Group.1",all = T)



## Second try
sd50_2<-aggregate(lemured4$length,by=list(lemured4$id),FUN=sd)

mean50_2=aggregate(lemured4$length,by=list(lemured4$id),FUN=mean)
res_lemu50_2<-merge(mean50_2,sd50_2,by="Group.1",all = T)


##Third try
sd50_3<-aggregate(lemured4b$length,by=list(lemured4b$id),FUN=sd)

mean50_3=aggregate(lemured4b$length,by=list(lemured4b$id),FUN=mean)
res_lemu50_3<-merge(mean50_3,sd50_3,by="Group.1",all = T)



list_df50<-list(mean50_1,mean50_2,mean50_3,mean100)
mean_length50<-list_df50%>%
  reduce(inner_join,by="Group.1")

colnames(mean_length50)<-c("id","mean50_1","mean50_2","mean50_3","mean100")




melt50<-melt(mean_length50,id="id")
colnames(melt50)<-c("id","replicates","length")

melt50 %>% group_by(replicates) %>% shapiro_test(length) # normally distributed

melt50 %>% t_test(length~replicates,paired = T) # no differences between replicates








##### 25% random samples (3 times) #####
set.seed(1)
lemured5<-lemu %>% 
  group_by(id_sample) %>%
  sample_frac(.25)



lemured6<-lemu %>% 
  group_by(id_sample) %>%
  sample_frac(.25)


lemured6b<-lemu %>% 
  group_by(id_sample) %>%
  sample_frac(.25)




### Mean length comparison


## First try
sd25_1<-aggregate(lemured5$length,by=list(lemured5$id),FUN=sd)

mean25_1=aggregate(lemured5$length,by=list(lemured5$id),FUN=mean)
res_lemu25_1<-merge(mean25_1,sd25_1,by="Group.1",all = T)



## Second try
sd25_2<-aggregate(lemured6$length,by=list(lemured6$id),FUN=sd)

mean25_2=aggregate(lemured6$length,by=list(lemured6$id),FUN=mean)
res_lemu25_2<-merge(mean25_2,sd25_2,by="Group.1",all = T)



##Third try
sd25_3<-aggregate(lemured6b$length,by=list(lemured6b$id),FUN=sd)

mean25_3=aggregate(lemured6b$length,by=list(lemured6b$id),FUN=mean)
res_lemu25_3<-merge(mean25_3,sd25_3,by="Group.1",all = T)



list_df25<-list(mean25_1,mean25_2,mean25_3,mean100)
mean_length25<-list_df25%>%
  reduce(inner_join,by="Group.1")

colnames(mean_length25)<-c("id","mean25_1","mean25_2","mean25_3","mean100")



melt25<-melt(mean_length25,id="id")
colnames(melt25)<-c("id","replicates","length")


melt25 %>% group_by(replicates) %>% shapiro_test(length) # normally distributed

melt25 %>% t_test(length~replicates,paired = T) # no differences between replicates












# 100%
ggplot(lemu,aes(x=length))+
  geom_density()+
  facet_wrap(vars(bank))



# 64%
ggplot(lemured,aes(x=length))+
  geom_density()+
  facet_wrap(vars(bank))


ggplot(lemured2,aes(x=length))+
  geom_density()+
  facet_wrap(vars(bank))


mean64stack<-cbind(mean_length64[1],stack(mean_length64[2:5]))

ggplot(mean64stack,aes(x=ind,y=values,fill=ind))+
  geom_boxplot()+
  labs(y="Length (mm)",x="Mean",title = "Average length of blue mussels in 36% reduced samples and base samples")+
  geom_signif(comparisons = list(c("mean64_1","mean100"),c("mean64_2","mean100"),c("mean64_3","mean100")),
                                  map_signif_level=TRUE, y_position = c(43, 41,39))+
  scale_fill_manual(values = c("#6E8B3D","#6E8B3D","#6E8B3D","#009ACD"))




# 50%
ggplot(lemured3,aes(x=length))+
  geom_density()+
  facet_wrap(vars(bank))


ggplot(lemured4,aes(x=length))+
  geom_density()+
  facet_wrap(vars(bank))


mean50stack<-cbind(mean_length50[1],stack(mean_length50[2:5]))

ggplot(mean50stack,aes(x=ind,y=values,fill=ind))+
  geom_boxplot()+
  labs(y="Length (mm)",x="Mean",title = "Average length of blue mussels in 50% reduced samples and base samples")+
  geom_signif(comparisons = list(c("mean50_1","mean100"),c("mean50_2","mean100"),c("mean50_3","mean100")),
              map_signif_level=TRUE, y_position = c(43, 41,39))+
  scale_fill_manual(values = c("#B452CD","#B452CD","#B452CD","#009ACD"))




# 25%
ggplot(lemured5,aes(x=length))+
  geom_density()+
  facet_wrap(vars(bank))


ggplot(lemured6,aes(x=length))+
  geom_density()+
  facet_wrap(vars(bank))



mean25stack<-cbind(mean_length25[1],stack(mean_length25[2:5]))

ggplot(mean25stack,aes(x=ind,y=values,fill=ind))+
  geom_boxplot()+
  labs(y="Length (mm)",x="Mean",title = "Average length of blue mussels in 75% reduced samples and base samples")+
  geom_signif(comparisons = list(c("mean25_1","mean100"),c("mean25_2","mean100"),c("mean25_3", "mean100")),
              map_signif_level=TRUE, y_position = c(43, 41,39))+
  scale_fill_manual(values = c("#CD5555","#CD5555","#CD5555","#009ACD"))















#### OYSTER & OYSTERBANK ####

oys$id_sample<- paste(oys$id,oys$sample)





### Add BANK_TYPE column
oys$bank <- as.vector(str_trim(str_extract(oys$id,"([:upper:][:upper:][:number:][:number:])")))
oys <- merge(x = oys, y = banktype_list[,c("bank", "BANK_TYPE")], by = "bank", all.x = TRUE)


### Subset with only Austernbank
oys<-oys[oys$BANK_TYPE=="Austernbank",]








##### 64% random samples (3 times) #####
set.seed(1)
oysred<-oys %>% 
  group_by(id_sample) %>%
  sample_frac(.64)



oysred2<-oys %>% 
  group_by(id_sample) %>%
  sample_frac(.64)



oysred2b<-oys %>% 
  group_by(id_sample) %>%
  sample_frac(.64)



### Mean length comparison
sd100o<-aggregate(oys$length,by=list(oys$id),FUN=sd)

mean100o=aggregate(oys$length,by=list(oys$id),FUN=mean)
res_oys100<-merge(mean100o,sd100o,by="Group.1",all = T)





## First try
sd64_1o<-aggregate(oysred$length,by=list(oysred$id),FUN=sd)

mean64_1o=aggregate(oysred$length,by=list(oysred$id),FUN=mean)
res_oys64_1<-merge(mean64_1o,sd64_1o,by="Group.1",all = T)



## Second try
sd64_2o<-aggregate(oysred2$length,by=list(oysred2$id),FUN=sd)

mean64_2o=aggregate(oysred2$length,by=list(oysred2$id),FUN=mean)
res_oys64_2<-merge(mean64_2o,sd64_2o,by="Group.1",all = T)


## Third try
sd64_3o<-aggregate(oysred2b$length,by=list(oysred2b$id),FUN=sd)

mean64_3o=aggregate(oysred2b$length,by=list(oysred2b$id),FUN=mean)
res_oys64_3<-merge(mean64_3o,sd64_3o,by="Group.1",all = T)



list_df64o<-list(mean64_1o,mean64_2o,mean64_3o,mean100o)
mean_length64o<-list_df64o%>%
  reduce(inner_join,by="Group.1")

colnames(mean_length64o)<-c("id","mean64_1","mean64_2","mean64_3","mean100")




melt64o<-melt(mean_length64o,id="id")
colnames(melt64o)<-c("id","replicates","length")

melt64o %>% group_by(replicates) %>% shapiro_test(length) # not normally distributed

melt64o %>% wilcox_test(length~replicates,paired = T) # no differences between replicates







##### 50% random samples (3 times) #####
set.seed(1)
oysred3<-oys %>% 
  group_by(id_sample) %>%
  sample_frac(.5)



oysred4<-oys %>% 
  group_by(id_sample) %>%
  sample_frac(.5)


oysred4b<-oys %>% 
  group_by(id_sample) %>%
  sample_frac(.5)




### Mean length comparison


## First try
sd50_1o<-aggregate(oysred3$length,by=list(oysred3$id),FUN=sd)

mean50_1o=aggregate(oysred3$length,by=list(oysred3$id),FUN=mean)
res_oys50_1<-merge(mean50_1o,sd50_1o,by="Group.1",all = T)



## Second try
sd50_2o<-aggregate(oysred4$length,by=list(oysred4$id),FUN=sd)

mean50_2o=aggregate(oysred4$length,by=list(oysred4$id),FUN=mean)
res_oys50_2<-merge(mean50_2o,sd50_2o,by="Group.1",all = T)




## Third try
sd50_3o<-aggregate(oysred4b$length,by=list(oysred4b$id),FUN=sd)

mean50_3o=aggregate(oysred4b$length,by=list(oysred4b$id),FUN=mean)
res_oys50_3<-merge(mean50_3o,sd50_3o,by="Group.1",all = T)



list_df50o<-list(mean50_1o,mean50_2o,mean50_3o,mean100o)
mean_length50o<-list_df50o%>%
  reduce(inner_join,by="Group.1")

colnames(mean_length50o)<-c("id","mean50_1","mean50_2","mean50_3","mean100")




melt50o<-melt(mean_length50o,id="id")
colnames(melt50o)<-c("id","replicates","length")

melt50o %>% group_by(replicates) %>% shapiro_test(length) # not normally distributed

melt50o %>% wilcox_test(length~replicates,paired = T) # no differences between replicates








##### 25% random samples (3 times) #####
set.seed(1)
oysred5<-oys %>% 
  group_by(id_sample) %>%
  sample_frac(.25)



oysred6<-oys %>% 
  group_by(id_sample) %>%
  sample_frac(.25)


oysred6b<-oys %>% 
  group_by(id_sample) %>%
  sample_frac(.25)




### Mean length comparison


## First try
sd25_1o<-aggregate(oysred5$length,by=list(oysred5$id),FUN=sd)

mean25_1o=aggregate(oysred5$length,by=list(oysred5$id),FUN=mean)
res_oys25_1<-merge(mean25_1o,sd25_1o,by="Group.1",all = T)



## Second try
sd25_2o<-aggregate(oysred6$length,by=list(oysred6$id),FUN=sd)

mean25_2o=aggregate(oysred6$length,by=list(oysred6$id),FUN=mean)
res_oys25_2<-merge(mean25_2o,sd25_2o,by="Group.1",all = T)


## Third try
sd25_3o<-aggregate(oysred6b$length,by=list(oysred6b$id),FUN=sd)

mean25_3o=aggregate(oysred6b$length,by=list(oysred6b$id),FUN=mean)
res_oys25_3<-merge(mean25_3o,sd25_3o,by="Group.1",all = T)



list_df25o<-list(mean25_1o,mean25_2o,mean25_3o,mean100o)
mean_length25o<-list_df25o%>%
  reduce(inner_join,by="Group.1")

colnames(mean_length25o)<-c("id","mean25_1","mean25_2","mean25_3","mean100")



melt25o<-melt(mean_length25o,id="id")
colnames(melt25o)<-c("id","replicates","length")


melt25o %>% group_by(replicates) %>% shapiro_test(length) # not normally distributed

melt25o %>% wilcox_test(length~replicates,paired = T) # no differences between replicates




# 100%
ggplot(oys,aes(x=length))+
  geom_density()+
  facet_wrap(vars(bank))


# 64%
ggplot(oysred,aes(x=length))+
  geom_density()+
  facet_wrap(vars(bank))


ggplot(oysred2,aes(x=length))+
  geom_density()+
  facet_wrap(vars(bank))


mean64stacko<-cbind(mean_length64o[1],stack(mean_length64o[2:5]))

ggplot(mean64stacko,aes(x=ind,y=values,fill=ind))+
  geom_boxplot()+
  labs(y="Length (mm)",x="Mean",title = "Average length of Pacific oysters in 36% reduced samples and base samples")+
  geom_signif(comparisons = list(c("mean64_1","mean100"),c("mean64_2","mean100"),c("mean64_3", "mean100")),
              map_signif_level=TRUE, y_position = c(98, 92,85))+
  scale_fill_manual(values = c("#458B74","#458B74","#458B74","#8EE5EE"))




# 50%
ggplot(oysred3,aes(x=length))+
  geom_density()+
  facet_wrap(vars(bank))


ggplot(oysred4,aes(x=length))+
  geom_density()+
  facet_wrap(vars(bank))




mean50stacko<-cbind(mean_length50o[1],stack(mean_length50o[2:5]))

ggplot(mean50stacko,aes(x=ind,y=values,fill=ind))+
  geom_boxplot()+
  labs(y="Length (mm)",x="Mean",title = "Average length of Pacific oysters in 50% reduced samples and base samples")+
  geom_signif(comparisons = list(c("mean50_1","mean100"),c("mean50_2","mean100"),c("mean50_3", "mean100")),
              map_signif_level=TRUE, y_position = c(98, 92,85))+
  scale_fill_manual(values = c("#9F79EE","#9F79EE","#9F79EE","#8EE5EE"))




# 25%
ggplot(oysred5,aes(x=length))+
  geom_density()+
  facet_wrap(vars(bank))


ggplot(oysred6,aes(x=length))+
  geom_density()+
  facet_wrap(vars(bank))



mean25stacko<-cbind(mean_length25o[1],stack(mean_length25o[2:5]))

ggplot(mean25stacko,aes(x=ind,y=values,fill=ind))+
  geom_boxplot()+
  labs(y="Length (mm)",x="Mean",title = "Average length of Pacific oyster in 75% reduced samples and base samples")+
  geom_signif(comparisons = list(c("mean25_1","mean100"),c("mean25_2","mean100"),c("mean25_3", "mean100")),
              map_signif_level=TRUE, y_position = c(88, 85,82))+
  scale_fill_manual(values = c("#CD2626","#CD2626","#CD2626","#8EE5EE"))















#### Difference of density over the years ####
Oybed<-read.table("HCA_Oysterbed_length.txt",sep = ",")




## Better view of the dataset
colnames(Oybed)<-c("bank","date","BANK_TYPE","method","sample","species","length")
Oybed$species[Oybed$species=="2225"]<-"mussel"
Oybed$species[Oybed$species=="1353"]<-"oyster"
Oybed$date<-format(as.Date(Oybed$date),"%Y")





oy_mod <- Oybed
oy_mod$yea_mon <- paste(format(as.Date(oy_mod$date),"%Y"),format(as.Date(oy_mod$date),"%m"),sep = "_")

oy_mod$year <- format(as.Date(oy_mod$date),"%Y")
oy_mod$month <- format(as.Date(oy_mod$date),"%m")


oy_mod_sel <- oy_mod[oy_mod$month=="04"|oy_mod$month=="05"|oy_mod$month=="09"|oy_mod$month=="10",]

oy_mod_sel <- oy_mod_sel%>%
  mutate(season = case_when(oy_mod_sel$month =="04"|oy_mod_sel$month =="05"~ "spring",
                            oy_mod_sel$month =="09"|oy_mod_sel$month =="10" ~ "autumn"))

library(dplyr)                                        
oy_mod_sel_unique <- oy_mod_sel %>% distinct(bank, yea_mon, .keep_all = TRUE)

overview_oy <- as.data.frame(unclass(table(paste(oy_mod_sel_unique$year,oy_mod_sel_unique$season,sep = "_"),oy_mod_sel_unique$bank)))














#add sub samp area column
require(dplyr)
Oybed <- Oybed%>%
  mutate(sub_samp_area = case_when(grepl("325",Oybed$method) ~ 176.7145867644259,
                                   grepl("332",Oybed$method) ~ 625,
                                   grepl("333",Oybed$method) ~ 625,
                                   grepl("326",Oybed$method) ~ 314.1593,
                                   grepl("331",Oybed$method) ~ 2500,
                                   grepl("336",Oybed$method) ~ 113.0973))





#### 1. Mussel ####
oybedmuss<-subset(Oybed,species=="mussel")

oybedmuss<-oybedmuss%>%mutate_all(~replace(., is.na(.), 625))

densoybedmuss=oybedmuss %>% count(sample, bank,date,sub_samp_area)





##add density column
densoybedmuss=densoybedmuss%>%
  mutate(density=n/sub_samp_area)



#Calculate mean density 

#sd_Dmuss=  densoybedmuss %>% group_by(id) %>% summarize(mean_dens = mean(density),sd_dens = sd(density))



#### 1.a Plot of some banks ####
densoybedmuss$year <- format(as.Date(densoybedmuss$date),"%Y")
densoybedmuss$month <- format(as.Date(densoybedmuss$date),"%m")

densoybedmuss2 <- densoybedmuss[densoybedmuss$month=="04"|densoybedmuss$month=="05"|densoybedmuss$month =="06"|densoybedmuss$month =="08"|densoybedmuss$month=="09"|densoybedmuss$month=="10",]






densoybedmuss2 <- densoybedmuss2%>%
  mutate(season = case_when(densoybedmuss2$month =="04"|densoybedmuss2$month =="05"|densoybedmuss2$month =="06"~ "spring",
                            densoybedmuss2$month =="08"|densoybedmuss2$month =="09"|densoybedmuss2$month =="10" ~ "autumn"))





spr<-ggplot(densoybedmuss2[densoybedmuss2$season=="spring"&densoybedmuss2$bank=="LT01",], aes(x = year, y = density))+
  geom_boxplot()+
  facet_wrap(.~bank, scales = "free")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  stat_summary(fun = mean,
               geom = "line",
               aes(group = 1),
               col = "red",linewidth=1)+
  labs(title = "Spring",y="density (ind/cm²)")+
  ylim(0,0.4)


aut<-ggplot(densoybedmuss2[densoybedmuss2$season=="autumn"&densoybedmuss2$bank=="LT01",], aes(x = year, y = density))+
  geom_boxplot()+
  facet_wrap(.~bank, scales = "free")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.title.y=element_blank())+ 
  stat_summary(fun = mean,
               geom = "line",
               aes(group = 1),
               col = "red",linewidth=1)+
  labs(title = "Autumn")+
  ylim(0,0.4)

spr+aut


#### 1.b Statistics tests for some years and some banks ####
sel_densoymuss <- densoybedmuss[densoybedmuss$bank=="LT01"|densoybedmuss$bank=="LT02"|densoybedmuss$bank=="LT13",]


LT01m<-sel_densoymuss[sel_densoymuss$bank=="LT01",]
LT02m<-sel_densoymuss[sel_densoymuss$bank=="LT02",]
#LT07m<-sel_densoymuss[sel_densoymuss$bank=="LT07",]
LT13m<-sel_densoymuss[sel_densoymuss$bank=="LT13",]
#NA23m<-sel_densoymuss[sel_densoymuss$bank=="NA23",]
#NA38m<-sel_densoymuss[sel_densoymuss$bank=="NA38",]



## Means and sd ##

LT01myears %>%
  group_by(year) %>%
  summarise_at(vars(density), list(mean = mean))

LT01myears %>%
  group_by(year) %>%
  summarise_at(vars(density), list(sd = sd))



LT02myears %>%
  group_by(year) %>%
  summarise_at(vars(density), list(mean = mean))

LT02myears %>%
  group_by(year) %>%
  summarise_at(vars(density), list(sd = sd))



LT13myears %>%
  group_by(year) %>%
  summarise_at(vars(density), list(mean = mean))

LT13myears %>%
  group_by(year) %>%
  summarise_at(vars(density), list(sd = sd))


### Select years ###
LT01myears<-LT01m[LT01m$year=="2005"|LT01m$year=="2010"|LT01m$year=="2015"|LT01m$year=="2020",]
LT02myears<-LT02m[LT02m$year=="2005"|LT02m$year=="2010"|LT02m$year=="2015"|LT02m$year=="2020",]
#LT07myears<-LT07m[LT07m$year=="2005"|LT07m$year=="2010"|LT07m$year=="2015"|LT07m$year=="2020",]
LT13myears<-LT13m[LT13m$year=="2005"|LT13m$year=="2010"|LT13m$year=="2015"|LT13m$year=="2020",]
#NA23myears<-NA23m[NA23m$year=="2005"|NA23m$year=="2010"|NA23m$year=="2015"|NA23m$year=="2020",]
#NA38myears<-NA38m[NA38m$year=="2005"|NA38m$year=="2010"|NA38m$year=="2015"|NA38m$year=="2020",]


shapiro.test(LT01myears$density) #Non normal
shapiro.test(LT02myears$density) #Normal
#shapiro.test(LT07myears$density) #Non normal
shapiro.test(LT13myears$density) #Non normal
#shapiro.test(NA23myears$density) #Non normal
#shapiro.test(NA38myears$density) #Non normal



pairwise.wilcox.test(LT01myears$density,LT01myears$year,p.adjust.method = "bonferroni")
pairwise.wilcox.test(LT02myears$density,LT02myears$year,p.adjust.method = "bonferroni")
#pairwise.wilcox.test(LT07myears$density,LT07myears$year,p.adjust.method = "bonferroni")
pairwise.wilcox.test(LT13myears$density,LT13myears$year,p.adjust.method = "bonferroni")
#pairwise.wilcox.test(NA23myears$density,NA23myears$year,p.adjust.method = "bonferroni")
#pairwise.wilcox.test(NA38myears$density,NA38myears$year,p.adjust.method = "bonferroni")






### Plot only for some years ###
sel_densoymuss_years<- sel_densoymuss[sel_densoymuss$year=="2005"|sel_densoymuss$year=="2010"|sel_densoymuss$year=="2015"|sel_densoymuss$year=="2020",]






### For LT banks
ggplot(sel_densoymuss_years[grep("^LT", sel_densoymuss_years$bank), ], aes(x = year, y = density,fill=year))+
  geom_boxplot()+
  geom_signif(comparisons = list(c("2005", "2020"),c("2010","2020"),c("2015","2020")),
              map_signif_level=TRUE, y_position = c(.7, .65,.6))+
  facet_wrap(.~bank, scales = "free")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title = "Mussel density over time on oyster beds LT01, LT02 and LT13",x="year",y="density (ind/cm²)")+
  scale_fill_brewer(palette="Set3")+
  ylim(0,.75)







### For NA23 bank
# ggplot(sel_densoymuss_years[sel_densoymuss_years$bank=="NA23",], aes(x = date, y = density,fill=date))+
#   geom_boxplot()+
#   geom_signif(comparisons = list(c("2015", "2010"),c("2010","2020"),c("2015","2020")),
#               map_signif_level=TRUE, y_position = c(0.65, 0.55,0.45))+
#   facet_wrap(.~bank, scales = "free")+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   labs(title = "Mussel density over time on oyster bed NA23",x="Date (year)",y="Density (ind/cm2)")+
#   scale_fill_brewer(palette="Set3")





### For NA38 bank
# ggplot(sel_densoymuss_years[sel_densoymuss_years$bank=="NA38",], aes(x = date, y = density,fill=date))+
#   geom_boxplot()+
#   geom_signif(comparisons = list(c("2010","2020")),
#               map_signif_level=TRUE, y_position = 0.35)+
#   facet_wrap(.~bank, scales = "free")+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   labs(title = "Mussel density over time on oyster bed NA38",x="Date (year)",y="Density (ind/cm2)")+
#   scale_fill_brewer(palette="Set3")
# 
# 







#### 2. Oyster ####
oybedoys<-subset(Oybed,species=="oyster")

oybedoys<-oybedoys%>%mutate_all(~replace(., is.na(.), 625))

densoybedoys=oybedoys %>% count(sample, bank,date,sub_samp_area)






##add density column
densoybedoys=densoybedoys%>%
  mutate(density=n/sub_samp_area)



#Calculate mean density 

#sd_Dmuss=  densoybedmuss %>% group_by(id) %>% summarize(mean_dens = mean(density),sd_dens = sd(density))



#### 1.a Plot of some banks ####
densoybedoys$year <- format(as.Date(densoybedoys$date),"%Y")
densoybedoys$month <- format(as.Date(densoybedoys$date),"%m")

densoybedoys2 <- densoybedoys[densoybedoys$month=="04"|densoybedoys$month=="05"|densoybedoys$month =="06"|densoybedoys$month =="08"|densoybedoys$month=="09"|densoybedoys$month=="10",]





densoybedoys2 <- densoybedoys2%>%
  mutate(season = case_when(densoybedoys2$month =="04"|densoybedoys2$month =="05"|densoybedoys2$month =="06"~ "spring",
                            densoybedoys2$month =="08"|densoybedoys2$month =="09"|densoybedoys2$month =="10" ~ "autumn"))







spr2<-ggplot(densoybedoys2[densoybedoys2$season=="spring"&densoybedoys2$bank=="LT01",], aes(x = year, y = density))+
  geom_boxplot()+
  facet_wrap(.~bank, scales = "free")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  stat_summary(fun = mean,
               geom = "line",
               aes(group = 1),
               col = "red",linewidth=1)+
  labs(title = "Spring",y="density (ind/cm²)")+
  ylim(0,0.4)


aut2<-ggplot(densoybedoys2[densoybedoys2$season=="autumn"&densoybedoys2$bank=="LT01",], aes(x = year, y = density))+
  geom_boxplot()+
  facet_wrap(.~bank, scales = "free")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.title.y=element_blank())+ 
  stat_summary(fun = mean,
               geom = "line",
               aes(group = 1),
               col = "red",linewidth=1)+
  labs(title = "Autumn")+
  ylim(0,0.4)

spr2+aut2



#### 2.b Statistics tests for some years and some banks ####
sel_densoyoys <- densoybedoys[densoybedoys$bank=="LT01"|densoybedoys$bank=="LT02"|densoybedoys$bank=="LT13",]


LT01o<-sel_densoyoys[sel_densoyoys$bank=="LT01",]
LT02o<-sel_densoyoys[sel_densoyoys$bank=="LT02",]
#LT07o<-sel_densoyoy[sel_densoyoy$bank=="LT07",]
LT13o<-sel_densoyoys[sel_densoyoys$bank=="LT13",]
#NA23o<-sel_densoyoy[sel_densoyoy$bank=="NA23",]
#NA38o<-sel_densoyoy[sel_densoyoy$bank=="NA38",]



shapiro.test(LT01o$density) #Non normal
shapiro.test(LT02o$density) #Non normal
#shapiro.test(LT07o$density) #Non normal
shapiro.test(LT13o$density) #Non normal
#shapiro.test(NA23o$density) #Non normal
#shapiro.test(NA38o$density) #Non normal




### Select years ###
LT01oyears<-LT01o[LT01o$year=="2005"|LT01o$year=="2010"|LT01o$year=="2015"|LT01o$year=="2020",]
LT02oyears<-LT02o[LT02o$year=="2005"|LT02o$year=="2010"|LT02o$year=="2015"|LT02o$year=="2020",]
#LT07oyears<-LT07o[LT07o$year=="2005"|LT07o$year=="2010"|LT07o$year=="2015"|LT07o$year=="2020",]
LT13oyears<-LT13o[LT13o$year=="2005"|LT13o$year=="2010"|LT13o$year=="2015"|LT13o$year=="2020",]
#NA23oyears<-NA23o[NA23o$year=="2005"|NA23o$year=="2010"|NA23o$year=="2015"|NA23o$year=="2020",]
#NA38oyears<-NA38o[NA38o$year=="2005"|NA38o$year=="2010"|NA38o$year=="2015"|NA38o$year=="2020",]


## Means and sd ##

LT01oyears %>%
  group_by(year) %>%
  summarise_at(vars(density), list(mean = mean))

LT01oyears %>%
  group_by(year) %>%
  summarise_at(vars(density), list(sd = sd))



LT02oyears %>%
  group_by(year) %>%
  summarise_at(vars(density), list(mean = mean))

LT02oyears %>%
  group_by(year) %>%
  summarise_at(vars(density), list(sd = sd))



LT13oyears %>%
  group_by(year) %>%
  summarise_at(vars(density), list(mean = mean))

LT13oyears %>%
  group_by(year) %>%
  summarise_at(vars(density), list(sd = sd))




pairwise.wilcox.test(LT01oyears$density,LT01oyears$year)
pairwise.wilcox.test(LT02oyears$density,LT02oyears$year)
#pairwise.wilcox.test(LT07oyears$density,LT07oyears$year)
pairwise.wilcox.test(LT13oyears$density,LT13oyears$year)
#pairwise.wilcox.test(NA23oyears$density,NA23oyears$year)
#pairwise.wilcox.test(NA38oyears$density,NA38oyears$year)






### Plot only for some years and some banks ###
sel_densoyoy_years<- sel_densoyoys[sel_densoyoys$year=="2005"|sel_densoyoys$year=="2010"|sel_densoyoys$year=="2015"|sel_densoyoys$year=="2020",]





### For LT banks
ggplot(sel_densoyoy_years[grep("^LT", sel_densoyoy_years$bank), ], aes(x = year, y = density,fill=year))+
  geom_boxplot()+
  geom_signif(comparisons = list(c("2005", "2020"),c("2010","2020"),c("2015","2020")),
              map_signif_level=TRUE, y_position = c(.55, .5,.45))+
  facet_wrap(.~bank, scales = "free")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title = "Oyster density over time on oyster beds LT01, LT02 and LT13",x="year",y="density (ind/cm²)")+
  scale_fill_brewer(palette="Set3")+
  ylim(0,.6)






### For NA23 bank
# ggplot(sel_densoyoy_years[sel_densoyoy_years$bank=="NA23",], aes(x = date, y = density,fill=date))+
#   geom_boxplot()+
#   geom_signif(comparisons = list(c("2015", "2010"),c("2010","2020"),c("2015","2020")),
#               map_signif_level=TRUE, y_position = c(0.25, 0.2,0.15))+
#   facet_wrap(.~bank, scales = "free")+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   labs(title = "Oyster density over time on oyster bed NA23",x="Date (year)",y="Density (ind/cm2)")+
#   scale_fill_brewer(palette="Set3")




### For NA38 bank
# ggplot(sel_densoyoy_years[sel_densoyoy_years$bank=="NA38",], aes(x = date, y = density,fill=date))+
#   geom_boxplot()+
#   geom_signif(comparisons = list(c("2010","2020")),
#               map_signif_level=TRUE, y_position = 0.2)+
#   facet_wrap(.~bank, scales = "free")+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   labs(title = "Oyster density over time on oyster bed NA38",x="Date (year)",y="Density (ind/cm2)")+
#   scale_fill_brewer(palette="Set3")
# 
# 
















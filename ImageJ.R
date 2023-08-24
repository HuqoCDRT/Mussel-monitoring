library(tidyverse)
library(rstatix)
library(ggpubr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
library(readxl)

Mussel_ImageJ_HCA_1 <- read_excel("Mussel_ImageJ_HCA_1.xlsx")
Test_ImageJ2 <- read_excel("Test_ImageJ2.xlsx")
Results2 <- read_csv("Results2.csv")
Results3 <- read_csv("Results3.csv")

### First try ###


res<- Results2%>% 
  mutate(r_number = case_when(Y >= 395 &  Y <= 407 ~ 1,
                              Y >= 347 &  Y <= 358 ~ 2,
                              Y >= 298 &  Y <= 308 ~ 3,
                              Y >= 249 &  Y <= 258 ~ 4,
                              Y >= 197 &  Y <= 211 ~ 5,
                              Y >= 146 &  Y <= 158 ~ 6,
                              Y >= 92 &  Y <= 102 ~ 7,
                              Y >= 42 &  Y <= 52 ~ 8))



res$r_number = factor(res$r_number, levels=c(1,2,3,4,5,6,7,8))


res<-res %>% arrange(r_number, desc(-X))

res_norubbish <- res[!(res$...1==63|is.na(res$r_number)),]

res_norubbish$real_order <- 1:104


table_mussel<-merge(Mussel_ImageJ_HCA_1,res_norubbish,by.x = "mussel",by.y = "real_order")




test<-melt(table_mussel,id.vars="mussel",measure.vars = c("Hugo","Laurens","Chris","Feret"))

colnames(test)<-c("mussel","method","length")








### Test statistic ###
test%>%
  group_by(method)%>%
  shapiro_test(length) #normally distributed


test %>% t_test(length~method) #no differences



ggplot(test, aes(x = mussel,color=method)) + 
  geom_point(aes(y = length,shape=method)) + 
  geom_line(aes(y = length,linetype=method)) + 
  labs(title="Comparison of length measurements between ImageJ and 3 hand measurements",
       x="Number of mussel",y="Length (mm)")+ 
  scale_x_continuous(n.breaks = 10)+
  scale_y_continuous(n.breaks = 10)+
  scale_color_brewer(palette = "Set1")



### Second try ###


res2<- Results3%>% 
  mutate(r_number = case_when(Y >= 395 &  Y <= 407 ~ 1,
                              Y >= 342 &  Y <= 356 ~ 2,
                              Y >= 290 &  Y <= 303 ~ 3,
                              Y >= 242 &  Y <= 254 ~ 4,
                              Y >= 193 &  Y <= 201 ~ 5,
                              Y >= 142 &  Y <= 149 ~ 6,
                              Y >= 91 &  Y <= 97 ~ 7,
                              Y >= 34 &  Y <= 45 ~ 8))



res2$r_number = factor(res2$r_number, levels=c(1,2,3,4,5,6,7,8))


res2<-res2 %>% arrange(r_number, desc(-X))



res2$real_order <- 1:99


new.row <- data.frame(A=11, B="K", stringsAsFactors=F)
res2 <- rbind.fill(res2, new.row)


table_mussel2<-merge(Test_ImageJ2,res2,by.x = "mussel",by.y = "real_order")

names(table_mussel2)[2] <- "Hugo"

test2<-melt(table_mussel2,id.vars="mussel",measure.vars = c("Hugo","Feret"))

colnames(test2)<-c("mussel","method","length")







### Test statistic ###

test2%>%
  group_by(method)%>%
  shapiro_test(length) #non normally distributed


test2 %>% wilcox_test(length~method) #no differences


ggplot(test2, aes(x = mussel,color=method)) + 
  geom_point(aes(y = length,shape=method)) + 
  geom_line(aes(y = length,linetype=method)) + 
  labs(title="Comparison of length measurements between ImageJ and 1 hand measurement",
       x="Number of mussel",y="Length (mm)")+ 
  scale_x_continuous(n.breaks = 10)+
  scale_y_continuous(n.breaks = 10)+
  scale_color_brewer(palette = "Set1")




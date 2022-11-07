library(readxl)
library(Matrix)
library(tidyverse)
library(brms)
library(tidybayes)
library(magrittr)
library(gridExtra)

#Data setup
dat<-read_excel("Range_Test.xlsx", sheet = "Sheet1")
V16<-filter(dat, Tag=="V16")
V13<-filter(dat, Tag=="V13")

############################################################
########### General Summaries 
############################################################
#Visualize the relationships
ggplot(V16, aes(x=Distance, y=Detections/20, col=(Wind_avg), shape=Depth, group=Depth)) + 
  geom_point(size=2, position=position_jitter(width=0.5))+
  geom_smooth(method="loess")+
  facet_wrap("Receiver")

ggplot(V13, aes(x=Distance, y=Detections/108, col=(Wind_avg), shape=Depth, group=Depth)) + 
  geom_point(size=2, position=position_jitter(width=0.5))+
  geom_smooth(method="loess")+
  facet_wrap("Receiver")

##Summaries
V16 %>% 
  group_by(Depth, Receiver, Distance) %>% 
  summarize(M=mean(Detections), 
            n=length(Detections))

V13 %>% 
  group_by(Depth, Receiver, Distance) %>% 
  summarize(M=mean(Detections),
            n=length(Detections), 
            Prob=M/108, 
            SE=sd(Detections)/n)

#########################################################################
############### Model preparation
########################################################################
#Determine the necessary priors
get_prior(Detections|trials(Possible)~Distance*Depth+Wind_avg*Depth+Receiver,
          data=V16,
          family = binomial(link="logit"))

#########################################################################################
########  V16
#########################################################################################
V16_mod<-brm(Detections|trials(Possible)~Distance*Depth+Wind_avg*Depth+Receiver,
             data=V16,
             family = binomial(link="logit"), 
             prior = c (prior(normal(3,0.5),class="Intercept"),
                        prior(normal(0,1), coef="DepthShallow"), 
                        prior(normal(-1,.5), class="b"), 
                        prior(normal(1,1), coef="ReceiverPipe")),
             #sample_prior = "only",
             file="V16_model.rds",
             chains=4, iter=2000, cores=4)

#View model results and checks
print(V16_mod, digits=4)
pp16 <- pp_check(V16_mod, ndraws=50) 
plot(V16_mod)
bayes_R2(V16_mod)

#Sample from the posterior
new.dat<-tibble(expand_grid(Receiver=c("Midriver", "Pipe"),
                Distance=c(100, 200, 300),
                Wind_avg=c(10,25,40),
                Depth=c("Shallow","Deep"), 
                Possible=20))

draw<-add_epred_draws(V16_mod, newdata=new.dat)

#Plot results
aa<-draw %>% 
  group_by(Distance, Depth, Receiver, Wind_avg) %>% 
  mean_qi(.epred) %>%   
  mutate(Wind_avg=as.factor(Wind_avg), 
         Mean_det_prob = .epred/20,
         Lower95_det_prob = .lower/20,
         Upper95_det_prob = .upper/20)
aa %>% 
  ggplot(aes(x=Distance, y=Mean_det_prob, col=Depth))+
  #geom_point(data = V16, aes(x=Distance, y=Detections/20), position=position_jitter(width=1) )+
  geom_ribbon(data=filter(aa, Depth=="Shallow", Receiver=="Midriver"),aes(ymin=Lower95_det_prob, ymax=Upper95_det_prob, fill=Wind_avg), 
              alpha=0.4, col=NA)+
  geom_ribbon(data=filter(aa, Depth=="Deep", Receiver=="Midriver"),aes(ymin=Lower95_det_prob, ymax=Upper95_det_prob, fill=Wind_avg), 
              alpha=0.4, col=NA)+
  geom_ribbon(data=filter(aa, Depth=="Shallow", Receiver=="Pipe"),aes(ymin=Lower95_det_prob, ymax=Upper95_det_prob, fill=Wind_avg), 
              alpha=0.4, col=NA)+
  geom_ribbon(data=filter(aa, Depth=="Deep", Receiver=="Pipe"),aes(ymin=Lower95_det_prob, ymax=Upper95_det_prob, fill=Wind_avg), 
              alpha=0.4, col=NA)+
  theme_classic()+
  scale_fill_manual(values=c("gray70", "gray70", "gray70"))+
  scale_color_manual(values=c("black", "gray50" ))+
  geom_line(aes(linetype=Wind_avg), size=1)+
  ylab("Detection probability")+
  xlab("Distance (m)")+
  facet_wrap(vars(Receiver),scales='free')+
  scale_y_continuous(limits=c(0,1))+
  theme(strip.background = element_blank(),
        axis.title = element_text(size=16), 
        axis.text = element_text(size=11),
        strip.text = element_text(size=16))

#sensitivity analysis
V16_sens<-brm(Detections|trials(Possible)~Distance*Depth+Wind_avg*Depth+Receiver,
             data=V16,
             family = binomial(link="logit"), 
             prior = c (prior(normal(3,1),class="Intercept"),
                        prior(normal(0,2), coef="DepthShallow"), 
                        prior(normal(-1,1), class="b"), 
                        prior(normal(1,2), coef="ReceiverPipe")),
             #sample_prior = "only",
             file="V16_sens.rds",
             chains=4, iter=2000, cores=4)

#View model results and checks
print(V16_sens, digits=4)
pp_check(V16_sens, ndraws=50)
plot(V16_sens)

#########################################################################################
########  V13
#########################################################################################
#model
V13_mod<-brm(Detections|trials(Possible)~Distance*Depth+Wind_avg*Depth+Receiver,
             data=V13,
             family = binomial(link="logit"), 
             prior = c (prior(normal(3,0.5),class="Intercept"),
                        prior(normal(0,1), coef="DepthShallow"), 
                        prior(normal(-1,.5), class="b"), 
                        prior(normal(1,1), coef="ReceiverPipe")),
             #sample_prior = "only",
             file="V13_model.rds",
             chains=4, iter=2000, cores=4)

#View model results and checks
print(V13_mod, digits=5)
pp13 <- pp_check(V13_mod, ndraws=50)
plot(V13_mod)
bayes_R2(V13_mod)

#Sample from the posterior
new.dat2<-tibble(expand_grid(Receiver=c("Midriver", "Pipe"),
                            Distance=c(100, 200, 300),
                            Wind_avg=c(10,25,40),
                            Depth=c("Shallow","Deep"), 
                            Possible=108))

draw2<-add_epred_draws(V13_mod, newdata=new.dat2)

#Plot results
bb<-draw2 %>% 
  group_by(Distance, Depth, Receiver, Wind_avg) %>% 
  mean_qi(.epred) %>%   
  mutate(Wind_avg=as.factor(Wind_avg), 
         Mean_det_prob = .epred/108,
         Lower95_det_prob = .lower/108,
         Upper95_det_prob = .upper/108)
bb %>% 
  ggplot(aes(x=Distance, y=Mean_det_prob, col=Depth))+
  #geom_point(data = V13, aes(x=Distance, y=Detections/108), position=position_jitter(width=1) )+
  geom_ribbon(data=filter(bb, Depth=="Shallow", Receiver=="Midriver"),aes(ymin=Lower95_det_prob, ymax=Upper95_det_prob, fill=Wind_avg), 
              alpha=0.4, col=NA)+
  geom_ribbon(data=filter(bb, Depth=="Deep", Receiver=="Midriver"),aes(ymin=Lower95_det_prob, ymax=Upper95_det_prob, fill=Wind_avg), 
              alpha=0.4, col=NA)+
  geom_ribbon(data=filter(bb, Depth=="Shallow", Receiver=="Pipe"),aes(ymin=Lower95_det_prob, ymax=Upper95_det_prob, fill=Wind_avg), 
              alpha=0.4, col=NA)+
  geom_ribbon(data=filter(bb, Depth=="Deep", Receiver=="Pipe"),aes(ymin=Lower95_det_prob, ymax=Upper95_det_prob, fill=Wind_avg), 
              alpha=0.4, col=NA)+
  theme_classic()+
  scale_fill_manual(values=c("gray70", "gray70", "gray70"))+
  scale_color_manual(values=c("black", "gray50" ))+
  geom_line(aes(linetype=Wind_avg), size=1)+
  ylab("Detection probability")+
  xlab("Distance (m)")+
  facet_wrap(vars(Receiver),scales='free')+
  scale_y_continuous(limits=c(0,1))+
  theme(strip.background = element_blank(),
        axis.title = element_text(size=16), 
        axis.text = element_text(size=11),
        strip.text = element_text(size=16))

#sensitivity analysis
V13_sens<-brm(Detections|trials(Possible)~Distance*Depth+Wind_avg*Depth+Receiver,
              data=V13,
              family = binomial(link="logit"), 
              prior = c (prior(normal(3,1),class="Intercept"),
                         prior(normal(0,2), coef="DepthShallow"), 
                         prior(normal(-1,1), class="b"), 
                         prior(normal(1,2), coef="ReceiverPipe")),
              #sample_prior = "only",
              file="V13_sens.rds",
              chains=4, iter=2000, cores=4)

#View model results and checks
print(V13_sens, digits=5)
pp_check(V13_mod)
plot(V13_mod)

#################################################################
########### Model pp_check final plot
#################################################################
pp16 <- pp16 + labs(title = "V16") + theme(plot.title = element_text(hjust=0.5))
pp13 <- pp13 +  labs(title = "V13") + xlab("Number of detections") + theme(plot.title = element_text(hjust=0.5))

(final <- grid.arrange(pp16, pp13))
ggsave("S1. Posterior predicitve check plot.jpeg", final, dpi=600, width=14, height=16, units="cm")

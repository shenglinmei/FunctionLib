



dat=data.frame('ERBB3'=t1,'NRG1'=t2)

p=ggplot(dat, aes(x=ERBB3, y=NRG1)) +
  geom_point(size=2, shape=23)+
  ggtitle(title)
p



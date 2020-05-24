#just playing with shape:



plot(1:150,cumsum( sim(params = set_constant_params(c(3,0.0,-7,0,0),
                                                    parameters = c("beta0","beta1","gama0","gama1","gama2"))
                       , N = 150)),ylab="cumulative",xlab = "year")

#seems like the range of values for beta1 is limited

plostlist <- vector(mode = "list")
i <- 1
for (beta1 in seq(-0.1,0.1,by = 0.01)){
  for (gama0 in seq(-3,3,by = 0.5)){
    plostlist[[i]] <- data.frame(beta1 = rep(beta1,200),gama0 = rep(gama0,200),sim = sim(params = set_constant_params(c(-1,beta1,gama0,0,0))
                                      , N = 200))
    i <- i + 1
  }
}



plotlist <- lapply(plostlist, function(x) x %>% as.tibble() %>% mutate(year = 1:200) %>% {ggplot(.)+aes(x=year,y = cumsum(sim))+geom_point()+
    ggtitle(paste("beta1 = ",.$beta1,"; gama0 = ",.$gama0,sep = ""))})


shape_plot <- gridExtra::arrangeGrob(grobs = plotlist,ncol = 20)
ggsave(shape_plot,filename = "shape.png",
       width = 60, height =  40,limitsize = F)




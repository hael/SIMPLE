################################################################################
#Script written by Frederic Bonnet on the 7th of March 2014                    #
#                                                                              #
#Scrit to evaluate the distribution of the generated random numbers from       #
#the fortran code.                                                             #
#                                                                              #
################################################################################

require(zoo)
require(quadprog)
require(tseries)
require(stat)

data<-read.table("dfxdfy_box_stat.asc")
summary(data)

microg<-data[[1]]
dfx<-data[[2]]
dfy<-data[[3]]
angast<-data[[4]]

add_phase_shift<-data[[5]]
ctf_res<-data[[6]]

df_ave = ( (dfx + dfy ) / 2 ) *(0.0001)


x11()
#postscript("R_output_fig/df_ave.eps")
plot(df_ave,main="Plot of the averaged defocus",xlab="n sample",ylab="(dfx+dfy)/2",col="blue",axes=TRUE)
#dev.off()

x11()
#postscript("R_output_fig/hist_df_ave.eps")
hist(df_ave,nclass=15,col="steelblue",prob=TRUE,main="Averaged defocus distribution (micron meter)");

#dev.off()

x11()
#postscript("R_output_fig/box_df_ave.eps")
boxplot(df_ave,ctf_res)

#dev.off()

x11()
#postscript("R_output_fig/hist_ctf_res.eps")
hist(ctf_res,nclass=15,col="steelblue",prob=TRUE,main="CTF resolution distribution");

#dev.off()

x11()
#postscript("R_output_fig/ctf_res.eps")
plot(ctf_res,main="Plot of the CTF resolution",xlab="n sample",ylab="(dfx+dfy)/2",col="blue",axes=TRUE)

#dev.off()

summary(ctf_res)

x11()
#postscript("R_output_fig/qqnorm_ctf_res.eps")
qqnorm(ctf_res,col="blue")
qqline(ctf_res)
#dev.off()

x11()
#postscript("R_output_fig/ecdf_ctf_res.eps")
plot(ecdf(ctf_res),main="Empirical cummulattive distribution")
#dev.off()

require(fGarch)
spec = garchSpec()
spec

x11()
fit=garchFit(~garch(1,1),data=ctf_res)
print(fit)
#postscript("R_output_fig/ctf_res_garch11.eps")
plot(fit)
1
#2
#3
#4
#5
#6
#7
#8
#9
#10
#11
#12
#13
#0
#dev.off()
sink("R_output_fig/ctf_res_garch11.asc",append=TRUE)
summary(fit)
sink()

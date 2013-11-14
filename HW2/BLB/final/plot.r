beta.se=matrix(unlist(read.table(file='./blb_lin_reg_data_s5_r50_SE.txt',header=TRUE)))
plot(1:1000,beta.se,main='Stand error of beta',
xlab='i-th coefficient',ylab='SE of beta',type='l')
savePlot(filename='se_beta',type='png')
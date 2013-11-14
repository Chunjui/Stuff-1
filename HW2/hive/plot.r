va=matrix(unlist(read.table('./var.txt')))
avg=matrix(unlist(read.table('./ave.txt')))
plot(avg,va,main='Means vs. Variance',xlab='within-group means',
ylab='within-group variance')
savePlot(filename='withinplot',type='png')
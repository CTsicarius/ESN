#1 VARIABLE PLOT
library(ggplot2)
library(reshape)
library(RColorBrewer)
row = 1
x_for_plot = y_train[1, ]
x_index = c(1:ncol(y_train))
df = data.frame('C1' = x_index, 'C2' = x_for_plot)
gg = ggplot(data=df, aes(x= C1, y = C2)) + theme_bw() +
            theme(text = element_text(size=14)) + theme(axis.text=element_text(size=rel(1.2))) + 
            geom_line() + xlab('Index') + ylab('Value')
plot(gg)
ggsave('ex1_output.pdf',plot = gg, device = 'pdf')

from = 1
library(ggplot2)
library(reshape)
library(RColorBrewer)
it_array = Nx_array
final = nrow(MSE_array_esn)
n_t = number_of_train[from:final]
df1 <- as.data.frame(cbind(n_t, cbind(MSE_array_var[from:final, ], MSE_array_esn[from:final, ])))
#df <- as.data.frame(cbind(number_of_train, MSE_array_esn))
df <- melt(df1, id.vars = 'n_t')
gg <- ggplot(data=df,
             aes(x=n_t, y=value, colour=variable)) +
  theme_bw() + theme(text = element_text(size=14)) + theme(axis.text=element_text(size=rel(1.2))) +
  geom_line() + geom_point() + scale_x_log10() + xlab('Size of trainig') + ylab('MSE error') +
  scale_color_manual(name = 'Nx', labels = c('VAR5 model', it_array),
                     values = c('blue', brewer.pal(length(it_array), "Paired")))
plot(gg)


#PARA PLOT SIN VAR
library(ggplot2)
library(reshape)
library(RColorBrewer)
it_array = Nx_array
df <- as.data.frame(cbind(number_of_train, MSE_array_esn))
#df <- as.data.frame(cbind(number_of_train, MSE_array_esn))
df <- melt(df, id.vars = 'number_of_train')
gg <- ggplot(data=df,
             aes(x=number_of_train, y=value, colour=variable)) +
  theme_bw() + theme(text = element_text(size=14)) + theme(axis.text=element_text(size=rel(1.2))) +
  geom_line() + geom_point() + scale_x_log10() + xlab('Size of trainig') + ylab('MSE error') +
  scale_color_manual(name = 'Nx', labels = it_array,
                     values = c(brewer.pal(length(it_array), "Paired")))
plot(gg)

#FOR SAVE
ggsave('sqrtbeta10.pdf',plot = gg, device = 'pdf')

# Import the dataset
signal <- read.delim('RORDEPs_signal_intensity.txt', header = 1, row.names = 1)
cor.test(signal$RORDEP1, signal$RORDEP2, method = 'spearman')


p_RORDEP1_RORDEP2 <- ggplot(signal, aes(x = log10(signal$RORDEP1), y = log10(signal$RORDEP2)))+
  geom_point(size = 3, shape = 1) + geom_smooth(method = lm, linetype="dashed",
                                                color="darkred", fill="blue")+
  theme(axis.text = element_text(size = 14, colour = 'black'), 
        axis.title = element_text(size = 16, color = 'black' ),
        plot.title = element_text(size = 14, color = 'black'))+
  xlab( "RORDEP1 intensity (log10 transformed)")+
  ylab("RORDEP2 intensity (log10 transformed)")+
  # labs(title = "Spearman Rho = 0.42, adjusted p = 3.3e-03")+
  xlim(5.8, 6.3)




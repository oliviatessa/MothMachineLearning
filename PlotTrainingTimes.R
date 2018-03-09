install.packages('ggplot2')
newDF = data.frame(Method = c("Sklearn_CPU", "Keras_GPU", "Keras_CPU"),                    Time = c(155.6,10, 26.22)/100)
newDF$Method = factor(newDF$Method)
newDF$Method = relevel(newDF$Method, ref = "Keras_GPU")


require(ggplot2)
version

theme_set(theme_classic())

ggplot(newDF, aes(x = Method, y = Time)) + 
  geom_bar(stat = "identity", width = 0.2) + 
  ylab( "Training Time per epoch (s)") + 
  coord_flip()

ggsave("E:\\Dropbox\\mothMachineLearning_dataAndFigs\\Figs\\TrainingTimes.pdf", width = 4, height = 1.5)


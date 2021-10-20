library(tidyverse)
library(ggplot2)
library(palmerpenguins)
ggplot(data=penguins,aes(x=bill_depth_mm))+
  geom_histogram()

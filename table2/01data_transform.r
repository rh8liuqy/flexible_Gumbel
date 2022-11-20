library(tidyverse)
library(lubridate)
library(readxl)
df1 <- read_xlsx("data.xlsx")
df1$datetime <- as_date(df1$datetime)

df1 <- df1 %>% 
  group_by(datetime) %>%
  summarise(change = (max(Y) - min(Y))*(which.max(Y) - which.min(Y)))

write.csv(df1,"data_toy.csv",row.names = FALSE)
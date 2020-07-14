# extra: Figure out why there is a systematic difference for BCoV copies between master mixes (706)
# going to do scatter plots for Ct and copues.ul and copies/l across weeks to figure this out

# pre run the weekly comparision.Rmd (dont knit, run chunks)

# scatter plot of Ct
results_abs %>% filter(Target == 'BCoV_M') %>%  ggplot(aes(x = Week, y = Ct)) + geom_jitter(width = .1)

df1 <- results_abs %>% 
  filter(Target == 'BCoV_M') %>%
  mutate(half_week = str_match(Tube_ID, '[:digit:]*(?= )') %>% {as.numeric(.) - as.numeric(as.character(Week)) } %>% as.factor(.) ) 

# Ct coloured by half samples per week (example : 706, 707)
df1 %>% ggplot(aes(x = Week, y = Ct, colour = half_week)) + geom_jitter(width = .1)

# Copies/ul RNA BCoV
plt1 <- df1 %>% ggplot(aes(x = Week, y = `Copies/ul RNA`, colour = half_week)) + geom_violin() + geom_jitter(width = .1) + ggtitle('BCoV data across weeks')
print(plt1)

# interactive copies/ul RNA BCoV
plotly::ggplotly(plt1)
# plotly::plot_ly(df1, x = ~Week, y = ~`Copies/ul RNA`, color = ~half_week) 

plt2 <- df1 %>% ggplot(aes(x = Week, y = `Copies/l WW`, colour = half_week)) + geom_violin() + geom_jitter(width = .1) + ggtitle('BCoV data across weeks') 
plotly::ggplotly(plt2)


# save plots manually
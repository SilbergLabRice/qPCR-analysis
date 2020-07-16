# extra: Figure out why there is a systematic difference for BCoV copies between master mixes (706)
# going to do scatter plots for Ct and copues.ul and copies/l across weeks to figure this out

# pre run the weekly comparision.Rmd (dont knit, run chunks)
library(plotly)

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

# interactive copies/l WW BCoV
plt2 <- df1 %>% ggplot(aes(x = Week, y = `Copies/l WW`, colour = half_week)) + geom_violin() + geom_jitter(width = .1) + ggtitle('BCoV data across weeks') 
plotly::ggplotly(plt2)

# Old vs new master mix plots
old_mm_707 <- str_c('707', ' [A-E]', '[1-3]') 
df2_by_mm <- df1 %>% mutate_when(str_detect(Week, '608|615|622|629'), list(half_week = factor(0)), str_detect(Tube_ID, old_mm_707), list(half_week = factor(0))) %>% 
  mutate('master mix' = str_replace_all(half_week, c('3' = '1', '0' = 'old (PCRbio)', '1' = 'New (Fastviral)')))

# Scatter plots
plt3 <- df2_by_mm %>% ggplot(aes(x = Week, y = `Copies/l WW`, colour = `master mix`)) + geom_violin() + geom_jitter(width = .1) + ggtitle('BCoV data across weeks') 
plotly::ggplotly(plt3)

# WWTP labelled plots
plt4 <- df2_by_mm %>% ggplot(aes(x = WWTP, y = `Copies/l WW`, fill = Week, shape = `master mix`)) + geom_jitter(width = .1) + ggtitle('BCoV data across weeks') + scale_shape_manual(values = c(21, 16))
plotly::ggplotly(plt4)

# WWTP labelled plots : Recovery fraction
plt5 <- df2_by_mm %>% filter(!str_detect(Week, '608|615')) %>% ggplot(aes(x = WWTP, y = `Recovery fraction`, colour = Week, shape = `master mix`)) + geom_jitter(width = .1) + ggtitle('BCoV data across weeks') + scale_shape_manual(values = c(21, 16))
plotly::ggplotly(plt5)

# Scatter : recovery fraction
plt6 <- df2_by_mm %>% ggplot(aes(x = Week, y = `Recovery fraction`, colour = `master mix`)) + geom_violin() + geom_jitter(width = .1) + ggtitle('BCoV data across weeks') 
ggplotly(plt6)

# Scatter : CT
df2_by_mm %>% ggplot(aes(x = Week, y = Ct, colour = `master mix`)) + geom_violin() + geom_jitter(width = .1) + ggtitle('BCoV data across weeks') 

# save plots manually


# Comparing old stuff with mixed master mix vials (RT qPCR and qPCR probe from PCRbio)
sheet_dump <- 'https://docs.google.com/spreadsheets/d/1ouk-kCJHERRhOMNP07lXfiC3aGB4wtWXpnYf5-b2CI4/edit#gid=730074403'

mixed_dat <- read_sheet(sheet_dump, sheet = 'WW12-608 part2a_BCoV')
dat_608 <- read_sheet(sheet_dump, sheet = 'WW14_608_part2')

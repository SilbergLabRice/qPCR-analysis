# Read in the file and do manual analysis and plotting
# Author: Prashant Kalvapalle;  October 10 2018

# The last part of the script contains important codes from the console that were used for specific situations : These will be commented

source('./general_functions.R') # Source the general_functions file before running this

# User inputs ----
# choose file name, title for plots and experiment mode (file name starts in the same directory as Rproject) 

flnm <- 'WW8_608_N1 N2 BCoV'  # set the filename
flpath <- str_c('excel files/',flnm,'.xls') # this completes the file path
plate_template_raw <- read_sheet('https://docs.google.com/spreadsheets/d/19oRiRcRVS23W3HqRKjhMutJKC2lFOpNK8aNUkC-No-s/edit#gid=478762118', sheet = 'Plate import setup', range = 'G1:S9')

title_name <-'6/8_Data1'

errorbar_width = 0.1; # width of errorbars - emperically change

# Spike in and concentration details
HA_concentration_factor <- 500 # concentration factor of wastewater -> RNA
spike_virus_conc <- 1124 # copies/ul viral suspension spiked in
spike_virus_volume <- 50 # ul of viral suspenses spiked in x ml WW; (x ~ 350 - 450 and varies for each sample)


# Assay mode features (choose if you want absolute quantification)
std_par <- tibble(                       # Input the slope and intercept from standard curve of various primer pairs/targets here - Target should match Target field (provided in excel sheet - Sample input reference.csv) 
  target = c('BRSV_N', 'BCoV_M', 'N1_CoV2', 'N2_CoV2', 'N1_multiplex',  'N2_multiplex'),
  slope =  c(-3.62, -3.49, -3, -3.12, -3.09, -3.1),
  intercept = c(39, 39, 39, 40, 39, 40) # values for various targets
)

# Input the data ----

# reading in file and polishing
fl <- readqpcr(flpath) # read excel file exported by Quantstudio

sample_order = columnwise_index(fl) # this gives a vector to order the samples columnwise in the PCR plate or strip (by default : data is shown row-wise) => This command will enable plotting column wise order
results_relevant <- fl$Results %>% select(`Well Position`, `Sample Name`, CT, starts_with('Tm'),`Target Name`) %>% rename(Target = `Target Name`) %>%  .[sample_order,] # select only the results used for plotting, calculations etc. and arrange them according to sample order

plate_template <- read_plate_to_column(plate_template_raw, 'Sample Name') # convert plate template (sample names) into a single vector, columnwise
results_relevant %<>% select(-`Sample Name`) %>% right_join(plate_template, by = 'Well Position') %>%  # Incorporate samples names from the google sheet
  mutate(Target = str_replace(Target, 'BSRV', 'BRSV'))  # correcting mis-spelled name of BRSV target

rm(fl, plate_template_raw)  # remove old data for sparsity

# Data cleaning ----
# (facetted by Sample category; naming: 'Sample Name'_variable primer pair)

# Separate the sample name into columns and make factors in the right order for plotting (same order as the plate setup)

# isolate the primer pair and assay_variable into 3 columns : Sample name, assay variable and primer pair 
results_relevant %<>% separate(`Sample Name`,c(NA, 'Sample Name'),'-') %>% 
  separate(`Sample Name`,c('Sample Name','Tube ID'),'_') %>% 
  mutate(`Tube ID` = if_else(`Sample Name` == 'NTC', '0', `Tube ID`)) %>% 
  separate(`Tube ID`, c('assay_variable', 'biological_replicates'), remove = F) %>%  # seperating biological replicates 1.1, 1.2, 1.3
  arrange(`Well Position`) # re-arrange the results in same order as the above factors (columnwise order of the plate)

# Computing copy number from standard curve linear fit information
results_abs <- results_relevant %>% group_by(Target)  %>% do(., absolute_backcalc(., std_par)) %>%  # iteratively calculates copy #'s from standard curve parameters of each Target
  unite('Biobot_ID', c('Sample Name', 'Tube ID'), sep = "", remove = F)

# WWTP identifiers ----
# Get ID to WWTP list from google sheet
biobot_lookup <- map_df(c('Week 9 (6/8)') , ~ read_sheet('https://docs.google.com/spreadsheets/d/1ghb_GjTS4yMFbzb65NskAlm-2Gb5M4SNYi4FHE4YVyI/edit#gid=233791008', sheet = .x, range = 'G:I')) %>% 
  mutate(Biobot_ID = str_remove(Biobot_ID,'\\.'), WWTP = as.character(SYMBOL), SYMBOL = NULL)

results_abs %<>% left_join(biobot_lookup, by = 'Biobot_ID') %>%  # join the results with the WWTP identifiers and names
  mutate(WWTP = if_else(is.na(WWTP), assay_variable, WWTP)) 

# Recovery calculations ----

# volumes_data <- read_sheet('https://docs.google.com/spreadsheets/d/1mJcCt1wMiOuBic6sRlBZJf8KSNu2y-B5PjzCUu7jPM8/edit#gid=521099478', sheet = 'Concentrated samples', range = 'A1:D38') %>% 
#   rename('WW_vol' = `Total WW vol (ml)`, 'Biobot_ID' = `Biobot/other ID`) %>% 
#   select(WW_vol, Biobot_ID) %>% 
#   fill(WW_vol) %>% distinct() %>% 
#   mutate(Biobot_ID = str_remove(Biobot_ID, " "))
#   
# results_abs %<>% left_join(volumes_data, by = 'Biobot_ID') %>%   # join the results with the WWTP identifiers and names
#   mutate(`Actual spike-in` = spike_virus_conc * spike_virus_volume / (WW_vol * 1e-3), Recovered = `Copy #` * 1e6/HA_concentration_factor, `Recovery fraction` = 100 * Recovered/`Actual spike-in`)

# Plotting ----


summary_results_abs <- results_abs %>%  group_by(`Sample Name`, Target, assay_variable, WWTP) %>% summarise_at(vars(`Copy #`), funs(mean(.,na.rm = T), sd(., na.rm = T))) # find mean and SD of individual copy #s for each replicate
results_abs %<>% replace_na(replace = list('Copy #' = 0, Recovered = 0)) # Make unamplified values zero: for plotting - not counted in the mean calculations

plt <- summary_results_abs %>% ggplot(aes(x = WWTP, y = mean, color = Target)) + ylab('Copies/ul RNA extract')    # Specify the plotting variables 

plt <- plt + geom_errorbar(aes(ymin = mean -sd, ymax = mean + sd, width = errorbar_width)) + # plot errorbars if mean and SD are desired
  geom_jitter(data = results_abs, aes(x = WWTP, y = `Copy #`, colour = Target), size = 1, width = .2) + # plot raw data
  geom_point(size = 2) +
  facet_grid(~`Sample Name`, scales = 'free_x', space = 'free_x') # plot points and facetting

# Formatting plot
plt.formatted <- plt %>% format_classic(., title_name, 'WWTP or Sample name') %>% format_logscale() # formatting plot, axes labels, title and logcale plotting

print(plt.formatted)

# write_sheet(results_abs,'https://docs.google.com/spreadsheets/d/1ouk-kCJHERRhOMNP07lXfiC3aGB4wtWXpnYf5-b2CI4/edit#gid=0', sheet = title_name) # save results to a google sheet
# ggsave('qPCR analysis/WW1_Baylor-bovine_pilot.png', plot = plt.formatted, width = 5, height = 4)
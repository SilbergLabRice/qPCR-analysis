# Read in the linregpcr processed data for data formatting and plotting (RAM paper)
# Author: Prashant Kalvapalle;  Sep 27 2021

source('./0-general_functions_main.R') # Source the general_functions file before running this

# User inputs  ----
# choose file name, title for plots (file name starts in the same directory as Rproject)

flnms <- c('q16c_Uxx-positives_15-10-21', 
           'q016b_Uxx-negatives_15-10-21')  # raw filename, without the "-processed" keyword

title_name <-'q16_Uxx variants'

# options

# x-axis label
axislabel.assay_variable <- 'Template name' # printed on the x axis of the graph


# Labelling translators ----

assay_var_translation <- c('295' = 'cat-RNA design 3', # U64 // U1376 # regex based translation to change x-axis labels
                           '297' = 'cat-RNA design 1',   #  U1 // U17
                           '298' = '(8)', # U8  --- Removed from further analysis
                           '299' = 'cat-RNA design 2',  # U34 // U891
                           '300' = 'cat-RNA design 4', # U73 // 1490
                           
                           '186' = 'split-gfp design', # gfp:gfp
                           
                           '103' = 'Empty vector', # previously : (-)
                           'NTC' = 'ntc') # use full form for easy understanding

target_name_translation <- c('34' = '2', '64' = '3', '73' = '4', # change to design numbers 
                             'U(.*)$' = 'rRNA-barcoded primer pair \\1', # add pre-fix
                             
                             'gfp:gfp' = 'gfp-complemented primer pair',
                             '16s' = '16S rRNA primer pair', # regex to change the target names for publication
                             'gfpbarcode' = 'cat-RNA barcode primer pair')

target_name_to_design <- c('34' = '2', '64' = '3', '73' = '4', # change to design numbers 
                           'U(.*)$' = 'design \\1', # add pre-fix
                           'gfp:gfp' = 'split-gfp design')

# Input the data ----

# reading in file and polishing
.df <- map_dfr(flnms, 
               # read linregpcr + R processed csv
               ~ read_csv(str_c('excel files/linregpcr_w_metadata/', .x , '-linreg-processed.csv'), 
                          col_types = 'ncccccnnnnnnnnnnnccnnn')
)


# Processing ---- 
 
# any custom processing goes here
forplot_reduced.data <- 
  
  # Removing U8 data due to false positives ~ mispriming
  filter(.df, !str_detect(Target_name, 'U8'), !str_detect(assay_variable, '298')) %>% 
  select(-1) %>% # remove first row which has numbering
  
  # translate to meaningful x axis labels
  mutate(across(assay_variable, ~ str_replace_all(.x, assay_var_translation)),
  
         # translate target names
         original_target_name = Target_name, # take backup
         across(Target_name, ~ str_replace_all(.x, target_name_translation)))
         
  
polished_data_for_output <- 
  
  mutate(forplot_reduced.data,
         
         # Make a design field for both template and primer set
         
         # Design = str_replace_all(original_target_name, target_name_to_design),          
         Design = if_else(str_detect(assay_variable, 'Empty vector|ntc') &
                            !str_detect(original_target_name, 'gfpbarcode|16s'), # for negative controls
                          # transform the original target name into design x numbering
                          str_replace_all(original_target_name, target_name_to_design),
                          str_remove(assay_variable, 'cat-RNA ')),

         .before = 1) %>%
  
  # # Refine target to be unique for each spliced design
  # mutate(across(Target_name, ~ if_else(.x == 'Spliced', # for spliced samples
  #                                      str_c(.x, ' (', Design, ')'), # attach the design name
  #                                      .x))) %>% 
  
  # Rearrange as factors in custom order
  
  # Make an ordered vector Sample from assay_variable
  mutate(Sample = fct_relevel(assay_variable, 'Empty vector') %>% # put empty vector first
           fct_relevel('ntc', after = Inf)) %>%  # put no template controls last

  # rearrange Target_name and Design  # Bring Spliced on the top
  mutate(across(Target_name, ~ fct_relevel(.x, '16S rRNA primer pair', 'cat-RNA barcode primer pair', after = Inf))) %>% 
  mutate(across(Design, ~ fct_relevel(.x, "Empty vector", 'ntc', after = Inf))) %>% 
  
  arrange(Design, Sample, Target_name) %>%  # arrange all rows in systematic order for output
  
  # :: Arrange columns
  select(Design, Sample, Target_name, contains('Copies'), everything()) %>% 
  select(-Sample_category, -matches('assay_var')) # remove redundant columns to avoid confusion



# summarize mean and stdev for replicates
polished.summary <- 
  group_by(polished_data_for_output, # group replicates
           Design, Sample, Target_name) %>% 
  
  summarize(mean_Copies.per.ul.template = mean(Copies.per.ul.template, na.rm = TRUE),
            stdev_Copies.per.ul.template = sd(Copies.per.ul.template, na.rm = TRUE))


# Save data ----

# save raw data
write_csv(polished_data_for_output,
          str_c('excel files/paper_data/2C_Uxx_variants', '-raw', '.csv'),
          na = 'nd')

# save summary data
write_csv(polished.summary,
          str_c('excel files/paper_data/2C_Uxx_variants', '-summary', '.csv'),
          na = 'nd')




# Plotting ----

plt.cq <- {plot_facetted_assay(.data = forplot_reduced.data,  # plotting function call with defaults
                               .yvar_plot = 40 - Cq,
                               .xvar_plot = assay_variable,
                               .colourvar_plot = Target_name, 
                               .facetvar_plot = NULL,
                               points_plt.style = 'jitter')  
  
  # geom_boxplot(aes(y = mean_Copies.per.ul.template), show.legend = FALSE) # plot mean
  
} %>% 
  
  # format_logscale_y() %>% # show logscale
  print()


# Plot Calibrated copies-linregpcr
plt.copies_w.mean <- {plot_facetted_assay(.data = forplot_reduced.data, 
                               .yvar_plot = Copies.per.ul.template,
                               .xvar_plot = assay_variable,
                               .xaxis.label.custom = axislabel.assay_variable,
                               .colourvar_plot = Target_name, 
                               .facetvar_plot = NULL,
                               points_plt.style = 'jitter') +  
  
  geom_boxplot(aes(y = mean_Copies.per.ul.template), show.legend = FALSE) # plot mean
} %>% 
  format_logscale_y() %>% # show logscale
  print()


# Troubleshooting: Plot N0-linregpcr
plt.N0 <- {plot_facetted_assay(.data = forplot_reduced.data, 
                                          .yvar_plot = N0,
                                          .xvar_plot = assay_variable,
                                          .xaxis.label.custom = axislabel.assay_variable,
                                          .colourvar_plot = Target_name, 
                                          .facetvar_plot = NULL,
                                          points_plt.style = 'jitter')  
    
    # geom_boxplot(aes(y = mean_Copies.per.ul.template), show.legend = FALSE) # plot mean
} %>% 
  format_logscale_y() %>% # show logscale
  print()


# Ratio analysis -----

# Calculating spliced fraction of the Ribozyme/16s

simple_target_translation <- c('U.*' = 'Spliced',
                             'gfp:gfp' = 'Spliced')


wider_reduced.dat <- 
  
  forplot_reduced.data %>% 
  
  # rename target for easier calculations
  mutate(Target_name = str_replace_all(original_target_name, simple_target_translation)) %>% 
  
  # select only required columns and rows
  select(Sample_category, Target_name, Copies.per.ul.template, assay_variable, biological_replicates) %>%  
  filter(Sample_category != 'negative') %>% # remove negative rows, with non-unique targets (pivoting issues)
  
  pivot_wider(names_from = Target_name, 
              values_from = Copies.per.ul.template, 
              names_prefix = 'Copies_') %>%  # each target to a col
  
  # calculate ratios of each targets per replicate
  mutate(Splicing_fraction_ribozyme = Copies_Spliced / Copies_gfpbarcode,
         Splicing_fraction_16s = Copies_Spliced / Copies_16s, 
         expression_ratio = Copies_gfpbarcode / Copies_16s) %>% 
  
  # find the mean across biological replicates for all numeric cols
  group_by(assay_variable) %>% 
  mutate(across(where(is.numeric), 
                ~ mean(.x, ignore.na = TRUE),
                .names = "mean_{.col}")) 



# Plotting ratios ----
# All code below does not work yet ------------------------

# plotting splicing fractions


plt.splicing_ratio_ribozyme <- 
  {plot_facetted_assay(.data = wider_reduced.dat,  # plotting function call with defaults
                       .yvar_plot = Splicing_fraction_ribozyme,
                       .colourvar_plot = NULL, # remove colour: since only 1 Sample_category
                       .facetvar_plot = organism,
                       points_plt.style = 'jitter') + 
      
      geom_boxplot(aes(y = mean_Splicing_fraction_ribozyme), show.legend = FALSE)  # plot mean
    # ggh4x::facet_nested(cols = vars(Target_name, organism), scales = 'free_x', space = 'free_x')   # redo facets
  } %>% 
  
  format_logscale_y() %>% # show logscale
  print()


plt.splicing_ratio_16s <- 
  {plot_facetted_assay(.data = wider_reduced.dat,  # plotting function call with defaults
                       .yvar_plot = Splicing_fraction_16s,
                       .colourvar_plot = NULL, # remove colour: since only 1 Sample_category
                       .facetvar_plot = organism,
                       points_plt.style = 'jitter') + 
    
    geom_boxplot(aes(y = mean_Splicing_fraction_16s), show.legend = FALSE)  # plot mean
  # ggh4x::facet_nested(cols = vars(Target_name, organism), scales = 'free_x', space = 'free_x') # redo facets
} %>% 
  
  format_logscale_y() %>% # show logscale
  print()

# plot expression ratio of ribozyme

plt.expression_ratio_ribozyme <- 
  {plot_facetted_assay(.data = wider_reduced.dat,  # plotting function call with defaults
                       .yvar_plot = expression_ratio,
                       .colourvar_plot = NULL, # remove colour: since only 1 Sample_category
                       .facetvar_plot = organism,
                       points_plt.style = 'jitter') + 
      
      geom_boxplot(aes(y = mean_expression_ratio), show.legend = FALSE)  # plot mean
    # ggh4x::facet_nested(cols = vars(Target_name, organism), scales = 'free_x', space = 'free_x')   # redo facets
  } %>% 
  
  format_logscale_y() %>% # show logscale
  print()



# Save plots ----

# Cq plot
ggsave(str_c('qPCR analysis/Archive/', title_name, '-Cq-linreg.png'),
       plt.cq,
       width = 5,
       height = 4)

ggsave(str_c('qPCR analysis/Archive/', title_name, '-thin.png'),
       plt.copies_w.mean,
       width = 5,
       height = 4)



# #  @@@@@@@   @@@@@@@              @@@@@@@  @@@@@@@@  @@@@@@@   @@@@@@    @@@@@@
# #  @@@@@@@@  @@@@@@@             @@@@@@@@  @@@@@@@@  @@@@@@@  @@@@@@@   @@@@@@@@
# #  @@!  @@@    @@!               !@@       @@!         @@!    !@@       @@!  @@@
# #  !@!  @!@    !@!               !@!       !@!         !@!    !@!       !@!  @!@
# #  @!@!!@!     @!!    @!@!@!@!@  !@!       @!!!:!      @!!    !!@@!!    @!@!@!@!
# #  !!@!@!      !!!    !!!@!@!!!  !!!       !!!!!:      !!!     !!@!!!   !!!@!!!!
# #  !!: :!!     !!:               :!!       !!:         !!:         !:!  !!:  !!!
# #  :!:  !:!    :!:               :!:       :!:         :!:        !:!   :!:  !:!
# #  ::   :::     ::                ::: :::   :: ::::     ::    :::: ::   ::   :::
# # NonParametric Multiparameter Analysis of CETSA/RT-CETSA Experimental Sets
# #
# # Written by: Michael Ronzetti {NIH/NCATS/UMD}
# # Patents: PCT/US21/45184, HHS E-022-2022-0-US-01
# # Main Analysis

library(tidyverse)
library(writexl)
source('functions.R')

# EXPERIMENTAL PARAMETERS AND SETUP
#
# Input experiment parameters here
startTemp <- 37
endTemp <- 95
plate_format <- 384
control <- 'vehicle'
pc <- 'control'

# Prepare the MatLab file for MoltenProt processing
moltprot_df <-
  prepMatLabforMolt(file_loc = './data/example_plate.xlsx',
                    start_temp = startTemp,
                    end_temp = endTemp)

# Read in the processed MoltProt data and prepare the data frames.
full_param <- retrieveMoltenData(model = 'standard')
curve_df <-
  retrieve_FittedCurves(model = 'baseline-fit',
                        start_temp = startTemp,
                        end_temp = endTemp)
full_df <- full_param %>%
  plate_assignment(., './data/platemap.xlsx')
full_df <- bind_fulldf(full_df, curve_df) %>%
  kelToCel(.)
full_df <- calculate_auc(full_df)

# Perform some preliminary control group analysis of variability
control_df <-
  control_grouping(full_df, control, pc) # Pull out control compound datapoints
control_var <-
  control_variability(control_df) # Read out the control group variability
controlPlot <-
  control_analysis(
    full_df,
    nc = 'vehicle',
    pc = 'control',
    output = 'plot',
    controlDF = control_var
  )
print(controlPlot)

#Calculate melting parameter difference for each well from MoltenProt
# full_df <- calculate_meltingparams(full_df) %>%
#   calculate_zscore() %>%
#   convert_zscore

#Derive RSS values for null and alternate model for each compound from full_df
rss <- compute.rss.models(full_df, rssPlot = TRUE, drPlot = TRUE, plotModel = TRUE)

#Perform dose-response for each thermogram parameter
parameters <- compute_parameter.rssmodel(full_df, plotModel = TRUE)

#Merge these plots for further analysis
signif_df <- merge(rss, parameters)
colnames(signif_df)[9] <- 'mannwhit.pval'
signif_df <- determineSig(signif_df)
signif_df <- rankOrder(signif_df)

# Volcano plots comparing the different parameters of analysis against the NPARC RSS Difference
# Colored by significance test and whether the compound passes any.
#plot_volcanos(signif_df)
# Plot of RSS Differences vs. p-values for NPARC
#rss.pval.plot(signif_df, savePlot = FALSE)
#Heatmap of compounds vs. different measurement styles.
parameter_heatmaps(signif_df, plotHeat = TRUE)

# Pull out raw data for reference #
raw_param <- retrieveMoltenData(model = 'standard')
raw_curves <- retrieve_FittedCurves(model = "raw_data",
                                  start_temp = startTemp,
                                  end_temp = endTemp)
raw_df <- raw_param %>% plate_assignment(., './data/platemap.xlsx')
raw_df <- cbind(raw_param, raw_curves) %>%
  plate_assignment(., './data/platemap.xlsx') %>%
  kelToCel(.)
raw_df <- calculate_auc(raw_df)

# Single temp analysis
tmFit_df <- findClosestTmColumn(raw_df)
tmFit_fits <- fitRawTmData(tmFit_df)

# Write out all data from analysis
write_xlsx(list(
  "Raw data" = raw_df,
  "Raw data fits" = tmFit_fits,
  "NPARC" = full_df,
  "Significance tests" = signif_df
), path = "./data/results.xlsx")
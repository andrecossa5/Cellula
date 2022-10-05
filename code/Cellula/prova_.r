# GBC CBC classification workflow from Roman Hillje

#########################################################################################

# Libraries
library(tidyverse)


CBC_GBC_combinations <- read_tsv('/Users/IEO5505/Desktop/CBC_GBC_combinations.tsv.gz')
CBC_GBC_combinations


dividing_line_slope <- 0.9                  #0.9
dividing_line_intercept <- -0.8                # -1.1

p <-
  ggplot(CBC_GBC_combinations, aes(read_count, UMI_count)) +
  geom_abline(intercept = dividing_line_intercept, slope = dividing_line_slope, color = 'red') +
  geom_point(size = 0.01, alpha = 0.05, color = 'black') +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  annotation_logticks(sides = 'bl') +
  labs(x = 'Read count', y = 'UMI count') +
  theme_bw()
p


CBC_GBC_combinations <- CBC_GBC_combinations %>%
  mutate(
    value = log10(read_count) * dividing_line_slope + dividing_line_intercept,
    below_line = value > log10(UMI_count)
  )

CBC_GBC_combinations$below_line %>% summary()



CBC_with_single_GBC <- CBC_GBC_combinations %>%
  filter(
    read_count >= 50,
    UMI_count >= 3,
    below_line == TRUE) %>%
  group_by(CBC) %>%
  summarize(GBC_count = n()) %>%
  filter(GBC_count == 1) %>%
  pull(CBC) %>%
  unique()

length(CBC_with_single_GBC)



CBC_with_multiple_GBC <- CBC_GBC_combinations %>%
  filter(
    read_count >= 50,
    UMI_count >= 3,
    below_line == TRUE
  ) %>%
  group_by(CBC) %>%
  summarize(GBC_count = n()) %>%
  filter(GBC_count > 1) %>%
  pull(CBC) %>%
  unique()

length(CBC_with_multiple_GBC)



CBC_GBC_combinations <- CBC_GBC_combinations %>%
  mutate(
    value = log10(read_count) * dividing_line_slope + dividing_line_intercept,
    below_line = value > log10(UMI_count),
    class = case_when(
      read_count < 50 | UMI_count < 3 | below_line == FALSE ~ 'fail',
      CBC %in% CBC_with_single_GBC & below_line == TRUE ~ 'single',
      CBC %in% CBC_with_multiple_GBC & below_line == TRUE ~ 'multiple',
    ),
    class = factor(class, levels = c('fail','multiple','single'))
  ) %>%
  arrange(class)

CBC_GBC_combinations

table(CBC_GBC_combinations$class)


p <- 
  ggplot(CBC_GBC_combinations) +
  geom_point(aes(read_count, UMI_count, color = class), size = 0.01, alpha = 0.5) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  scale_color_manual(
    name = 'Class',
    values = c('single' = '#2980b9', 'multiple' = '#27ae60', 'fail' = '#c0392b'),
    guide = guide_legend(override.aes = list(size = 1, alpha = 1), reverse = TRUE)
  ) +
  labs(x = 'Read count', y = 'UMI count') +
  annotation_logticks(sides = 'lb') +
  theme_bw() +
  theme(legend.position = c(0.2,0.8))
p


CBC_GBC_combinations %>%
filter(class == 'single') %>%
select(CBC, GBC, read_count, UMI_count) %>%
write.csv('/Users/IEO5505/Desktop/summary_sheet_cells.csv')



CBC_GBC_combinations %>%
filter(class == 'single') %>%
group_by(GBC) %>%
summarize(
  cell_count = n(),
  cells = paste(CBC, collapse = ','),
  mean_read_count = mean(read_count),
  mean_UMI_count = mean(UMI_count),
  median_read_count = median(read_count),
  median_UMI_count = median(UMI_count)
) %>%
arrange(-cell_count) %>%
write.csv('/Users/IEO5505/Desktop/summary_sheet_clones.csv')



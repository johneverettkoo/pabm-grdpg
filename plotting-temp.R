clustering.df %>%
  dplyr::group_by(n, K) %>%
  dplyr::summarise(
    med.err = median(error),
    first.q = quantile(error, .25),
    third.q = quantile(error, .75),
    med.err.ssc = median(error.ssc),
    first.q.ssc = quantile(error.ssc, .25),
    third.q.ssc = quantile(error.ssc, .75),
    med.err.ep = median(error.ep, na.rm = TRUE),
    first.q.ep = quantile(error.ep, .25, na.rm = TRUE),
    third.q.ep = quantile(error.ep, .75, na.rm = TRUE),
    med.err.louvain = median(error.louvain),
    first.q.louvain = quantile(error.louvain, .25),
    third.q.louvain = quantile(error.louvain, .75)
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::inner_join(
    ssc.df %>% 
      dplyr::group_by(n, K) %>% 
      dplyr::summarise(med.err.ssc.A = median(error.ssc2),
                       first.q.ssc.A = quantile(error.ssc2, .25),
                       third.q.ssc.A = quantile(error.ssc2, .75)) %>% 
      dplyr::ungroup()
  ) %>% 
  ggplot() +
  theme_bw() + 
  theme(text = element_text(size = 20)) + 
  scale_x_log10(breaks = c(128, 256, 512, 1024, 2048, 4096)) +
  # scale_x_continuous(breaks = c(128, 256, 512, 1024, 2048, 4096)) + 
  scale_y_log10() +
  labs(y = 'community detection error count', 
       colour = NULL, shape = NULL, size = NULL) +
  geom_line(aes(x = n, y = med.err * n,
                colour = 'OSC')) +
  geom_point(aes(x = n, y = med.err * n,
                 colour = 'OSC', shape = 'OSC'), size = 3) +
  geom_errorbar(aes(x = n, ymin = first.q * n, ymax = third.q * n,
                    colour = 'OSC'), width = .1) + 
  geom_line(aes(x = n, y = med.err.ssc * n,
                colour = 'SSC-ASE')) +
  geom_point(aes(x = n, y = med.err.ssc * n,
                 colour = 'SSC-ASE', shape = 'SSC-ASE'), size = 3) +
  geom_errorbar(aes(x = n, ymin = first.q.ssc * n, ymax = third.q.ssc * n,
                    colour = 'SSC-ASE'), width = .1) +
  geom_line(aes(x = n, y = med.err.ssc.A * n,
                colour = 'SSC-A')) +
  geom_point(aes(x = n, y = med.err.ssc.A * n,
                 colour = 'SSC-A', shape = 'SSC-A'), size = 3) +
  geom_errorbar(aes(x = n, ymin = first.q.ssc.A * n, ymax = third.q.ssc.A * n,
                    colour = 'SSC-A'), width = .1) +
  geom_line(aes(x = n, y = med.err.louvain * n,
                colour = 'MM-Louvain')) + 
  geom_point(aes(x = n, y = med.err.louvain * n,
                 colour = 'MM-Louvain', shape = 'MM-Louvain'), size = 3) + 
  geom_errorbar(aes(x = n, ymin = first.q.louvain * n, ymax = third.q.louvain * n,
                    colour = 'MM-Louvain'), width = .1) + 
  scale_colour_brewer(palette = 'Set1') + 
  facet_wrap(~ K, labeller = 'label_both')

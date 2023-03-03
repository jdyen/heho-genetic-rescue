# function to extract parameter draws from fitted Stan models
extract_parameters <- function(x, pars, model, transform = identity, ...) {
  
  # grab the parameters we want
  tmp <- as.matrix(x)[, pars, drop = FALSE]
  
  # add the intercept to all other terms if included
  if (any(grepl("Intercept", colnames(tmp))) & ncol(tmp > 1)) {
    if (ncol(tmp) == 2)
      tmp[, 2] <- tmp[, 2] + tmp[, 1]
    else
      tmp[, 2:ncol(tmp)] <- sweep(tmp[, 2:ncol(tmp)], 1, tmp[, 1], "+") 
    
    # and make an intercept name that is informative
    colnames(tmp)[1] <- "SEX_SPECIFIC_PAIR_TYPECC"
    
  }
  
  # now flatten and transform it back to response scale
  tmp <- data.frame(
    iter = seq_len(length(c(tmp))),
    value = transform(c(tmp)),
    group = rep(colnames(tmp), each = nrow(tmp)),
    model = model
  )
  
  # return
  tmp
  
}

# function to compare cross types based on the expected number of
#   independent chicks (overall summary result)
compare_cross_type <- function(pars, var = "expected_independent_chicks") {
  
  # filter to target variable
  tmp <- pars %>% filter(name == var)
  
  # clean up names of cross types
  tmp <- tmp %>%mutate(
    group = gsub("SEX_SPECIFIC_PAIR_TYPE", "", group)
  )
  
  # grab values for each cross type
  tmp <- data.frame(
    cc = tmp %>% filter(group == "CC") %>% pull(value),
    cff1m = tmp %>% filter(group == "CfF1m") %>% pull(value),
    cfgm = tmp %>% filter(group == "CfGm") %>% pull(value),
    f1f1 = tmp %>% filter(group == "F1F1") %>% pull(value),
    f1fcm = tmp %>% filter(group == "F1fCm") %>% pull(value),
    gfcm = tmp %>% filter(group == "GfCm") %>% pull(value) 
  )
  
  # calculate proportion of values greater than the CC cross type
  out <- rep(NA, ncol(tmp) - 1)
  for (i in 2:ncol(tmp))
    out[i - 1] <- mean(tmp[, i] > tmp[, 1])
  
  # return
  out
  
}

# plot parameters from the model on link or observation scales
plot_parameters <- function(
    x,
    pars,
    prob = c(0.8, 0.95),
    title = "",
    transform = identity
) {
  
  # grab the parameters we want
  tmp <- as.matrix(x)[, pars, drop = FALSE]
  
  # add a zero line to the plot
  int_line <- 0
  
  # add the intercept to all other terms if included
  if (any(grepl("Intercept", colnames(tmp))) & ncol(tmp) > 1) {
    if (ncol(tmp) == 2)
      tmp[, 2] <- tmp[, 2] + tmp[, 1]
    else
      tmp[, 2:ncol(tmp)] <- sweep(tmp[, 2:ncol(tmp)], 1, tmp[, 1], "+") 
    
    # shift ther zero line to the CC intercept if plotting on obs scale
    int_line <- median(transform(tmp[, 1]))
    
  }
  
  # now flatten so we can summarise and plot it (transform if required)
  tmp <- data.frame(
    value = transform(c(tmp)),
    group = rep(colnames(tmp), each = nrow(tmp))
  )
  
  # summarise time
  eps <- 1 - prob
  probs <- c(eps[2] / 2, eps[1] / 2, 1 - eps[1] / 2, 1 - eps[2] / 2)
  tmp <- tmp %>%
    dplyr::group_by(group) %>% 
    dplyr::mutate(
      group = gsub("SEX_SPECIFIC_PAIR_TYPE", "", group),
      group = gsub("\\(Intercept\\)", "CC", group)
    ) %>%
    summarise(
      lowest = quantile(value, probs = probs[1]),
      lower = quantile(value, probs = probs[2]),
      mid = median(value),
      upper = quantile(value, probs = probs[3]),
      upperest = quantile(value, probs = probs[4])
    ) %>%
    ungroup()
  
  # set plot limits
  x_lim <- range(c(tmp$lowest, tmp$upperest))
  x_range <- diff(x_lim)
  x_lim[1] <- x_lim[1] - 0.05 * x_range
  x_lim[2] <- x_lim[2] + 0.05 * x_range
  
  # plot parameter estimates
  p <- tmp %>%
    ggplot() +
    geom_vline(aes(xintercept = int_line), size = 1, col = "gray") +
    geom_point(aes(x = mid, y = group), size = 4) +
    geom_segment(aes(x = lower, xend = upper, y = group, yend = group), size = 2) +
    geom_segment(aes(x = lowest, xend = upperest, y = group, yend = group), size = 0.5) +
    labs(title = title) +
    xlim(x_lim) +
    scale_y_discrete(limits = unique(rev(tmp$group))) +
    yaxis_text(face = "bold") +
    yaxis_title(FALSE) +
    yaxis_ticks(size = 1) +
    xaxis_title(FALSE)
  
  # return
  p
  
}

# function to plot the expected number of independent chicks by cross type
plot_expected_independent_chicks <- function(
    pars, 
    prob = c(0.8, 0.95), 
    subset = NULL,
    var = "expected_independent_chicks",
    log_x = FALSE,
    ...
) {
  
  # pull out the overall outcome variable
  pars <- pars %>% filter(name == var)
  
  # summarise these values to return the quantiles used in the plot
  eps <- 1 - prob
  probs <- c(eps[2] / 2, eps[1] / 2, 1 - eps[1] / 2, 1 - eps[2] / 2)
  tmp <- pars %>%
    mutate(group = gsub("SEX_SPECIFIC_PAIR_TYPE", "", group)) %>%
    group_by(group) %>% 
    summarise(
      lowest = quantile(value, probs = probs[1]),
      lower = quantile(value, probs = probs[2]),
      mid = median(value),
      upper = quantile(value, probs = probs[3]),
      upperest = quantile(value, probs = probs[4])
    ) %>%
    ungroup()
  
  # plot parameter estimates
  p <- tmp %>%
    ggplot() +
    geom_vline(aes(xintercept = 0), size = 1, col = "gray") +
    geom_point(aes(x = mid, y = group), size = 3, position = position_dodge(0.8)) +
    geom_linerange(aes(xmin = lower, xmax = upper, y = group), size = 2, position = position_dodge(0.8)) +
    geom_linerange(aes(xmin = lowest, xmax = upperest, y = group), size = 0.5, position = position_dodge(0.8)) +
    scale_y_discrete(limits = rev(unique(tmp$group))) +
    xlab("Expected number of independent chicks") +
    yaxis_text(face = "bold") +
    yaxis_title(FALSE) +
    yaxis_ticks(size = 1) +
    theme(legend.position = "none")
  
  # scale x axis if needed
  if (log_x)
    p <- p + scale_x_log10(labels = scales::number_format(accuracy = 0.001))
  
  # return
  p
  
}

# function to extract random effects from fitted models
extract_random_effects <- function(x, term, model, ...) {
  
  # put all samples into a single matrix 
  tmp <- as.matrix(x)
  
  # grab the term we want
  tmp <- tmp[, grepl(term, colnames(tmp)), drop = FALSE]
  tmp <- tmp[, !grepl("Sigma", colnames(tmp)), drop = FALSE]
  
  # clean up the group names
  group <- rep(colnames(tmp), each = nrow(tmp))
  group <- gsub("b\\[|\\]|\\(|\\)|Intercept| ", "", group)
  
  # split up group and term names
  term <- strsplit(group, ":")
  group <- sapply(term, function(x) x[2])
  term <- sapply(term, function(x) x[1])
  
  # flatten and return
  data.frame(
    iter = seq_len(length(c(tmp))),
    value = c(tmp),
    group = group,
    term = term
  )
  
}

# function to plot random effects
plot_random_effects <- function(pars, prob = c(0.8, 0.95), subset = NULL) {
  
  # summarise values to calculate quantiles for plotting
  eps <- 1 - prob
  probs <- c(eps[2] / 2, eps[1] / 2, 1 - eps[1] / 2, 1 - eps[2] / 2)
  tmp <- pars %>%
    group_by(group, term) %>% 
    summarise(
      lowest = quantile(value, probs = probs[1]),
      lower = quantile(value, probs = probs[2]),
      mid = median(value),
      upper = quantile(value, probs = probs[3]),
      upperest = quantile(value, probs = probs[4])
    ) %>%
    ungroup()
  
  # plot parameter estimates
  p <- tmp %>%
    ggplot() +
    geom_vline(aes(xintercept = 0), size = 1, col = "gray") +
    geom_point(aes(x = mid, y = group), size = 3, position = position_dodge(0.8)) +
    geom_linerange(aes(xmin = lower, xmax = upper, y = group), size = 2, position = position_dodge(0.8)) +
    geom_linerange(aes(xmin = lowest, xmax = upperest, y = group), size = 0.5, position = position_dodge(0.8)) +
    facet_wrap( ~ term, scales = "free_y") +
    yaxis_text(face = "bold") +
    yaxis_title(FALSE) +
    yaxis_ticks(size = 1) +
    xaxis_title(FALSE)
  
  # return
  p
  
}

# function to calculate the probability that reproductive success is higher
#   than the CC cross type
probability_greater_than_others <- function(
    x, 
    pars = c(
      "(Intercept)",
      "SEX_SPECIFIC_PAIR_TYPECfF1m", 
      "SEX_SPECIFIC_PAIR_TYPECfGm",
      "SEX_SPECIFIC_PAIR_TYPEF1F1",
      "SEX_SPECIFIC_PAIR_TYPEF1fCm",
      "SEX_SPECIFIC_PAIR_TYPEGfCm"
    )
) {
  
  # grab the parameters we want
  draws <- as.matrix(x)[, pars, drop = FALSE]
  
  # add the intercept to all other terms if included
  if (any(grepl("Intercept", colnames(draws))) & ncol(draws) > 1) {
    if (ncol(draws) == 2)
      draws[, 2] <- draws[, 2] + draws[, 1]
    else
      draws[, 2:ncol(draws)] <- sweep(draws[, 2:ncol(draws)], 1, draws[, 1], "+") 
  }
  
  # calculate probability that a parameter is greater than each
  #   other cross type
  out <- matrix(NA, nrow = ncol(draws), ncol = ncol(draws))
  for (i in seq_len(ncol(draws))) {
    tmp <- sweep(draws, 1, draws[, i], "-")
    out[i, ] <- apply(tmp, 2, function(x) mean(x > 0))
  }
  
  # remove self-comparisons
  diag(out) <- NA
  
  # and some informative labels
  group_ids <- c("CC", "CfF1m", "CfGm", "F1F1", "F1fCm", "GfCm")
  out <- data.frame(
    group = paste0(rep(group_ids, each = nrow(out)), " > ", rep(group_ids, times = nrow(out))),
    value = c(out)
  )
  
  # and return
  out
  
}

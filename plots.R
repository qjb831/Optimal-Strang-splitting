library(plotly)
library(tidyr)
library(dplyr)
library(MASS)
library(expm)
library(ggnewscale)
library(patchwork)


plot_2d_curve <- function(df,line_size  = 0.2,palette    = "viridis") {
  if (!all(c("X1", "X2") %in% names(df))) {
    stop("df must contain columns 'X1' and 'X2'")
  }
  
  df$idx <- seq_len(nrow(df))
  cols <- viridis::viridis(256, option = palette)

  
  ggplot(df, aes(x = X1, y = X2, color = idx)) +
    geom_path(linewidth = line_size, lineend = "round") +
    scale_color_gradientn(colors = cols, name = "index", guide = "colourbar") +
    labs(x = "X", y = "Y") +
    theme_minimal()
}

plot_FHN_phase_portrait <- function(par=c(0.1,0.5,1.5,1.4),
                             x_range = c(-1.4, 1.4),
                             y_range = c(-0.3, 1.2),
                             n_points = 40,
                             arrow_scale = 0.12,
                             show_nullclines = TRUE,
                             show_equilibria = TRUE) {
  
  epsilon   <- par[1]; alpha <- par[2]; gamma <- par[3]; beta  <- par[4]
  f <- function(x,y) {  
    (-x^3+x+alpha-y)/epsilon
  }
  g <- function(x,y) {    # y' = (2 - 0.4*x)*x + (0.3 + 0.1*x)
    gamma*x+beta-y
  }
  # build grid
  x <- seq(x_range[1], x_range[2], length.out = n_points)
  y <- seq(y_range[1], y_range[2], length.out = n_points)
  grid <- expand.grid(x = x, y = y, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  
  # evaluate f and g on the grid 
  grid$xdot <- f(grid$x, grid$y)
  grid$ydot <- g(grid$x, grid$y)
  
  # magnitude and normalization
  grid$magnitude <- sqrt(grid$xdot^2 + grid$ydot^2)
  grid$xdot_n <- grid$xdot / grid$magnitude
  grid$ydot_n <- grid$ydot / grid$magnitude
  
  # compute arrow endpoints scaled to grid spacing
  dx <- (x_range[2] - x_range[1]) / max(1, (n_points - 1))
  dy <- (y_range[2] - y_range[1]) / max(1, (n_points - 1))
  base_step <- min(dx, dy)
  len <- arrow_scale * base_step
  grid$xend <- grid$x + len * grid$xdot_n
  grid$yend <- grid$y + len * grid$ydot_n
  
  # plot
  p <- ggplot(grid, aes(x = x, y = y)) +
    geom_raster(aes(fill = magnitude), interpolate = TRUE) +
    geom_segment(aes(xend = xend, yend = yend),
                 arrow = arrow(length = unit(0.08, "inches")), color = "black", alpha = 0.6) +
    labs(x = "X", y = "Y", fill = "||F(X,Y)||") +
    theme_minimal(base_size = 17)
  
  p <- p + scale_fill_viridis_c(option = "plasma",trans='log10')

  return(p)
}

plot_estimates <- function(est_df, true_params, param_labels,
                           subtitle = "", center_names = splitting_labels,
                           title_text = NULL) {
  # Create a data frame of true values repeated for each method
  true_df <- tibble(param = param_labels, true_value = true_params)
  true_df_expanded <- expand_grid(param = param_labels, center = center_names) %>%
    left_join(true_df, by = "param") %>%
    mutate(center = factor(center, levels = center_names))
  
  # Ensure est_df$center is factor with correct levels
  est_df <- est_df %>% mutate(center = factor(center, levels = center_names))
  est_df <- est_df %>%
    group_by(param) %>%
    filter(
      estimate >= quantile(estimate, 0.01, na.rm = TRUE),
      estimate <= quantile(estimate, 0.99, na.rm = TRUE)
    )
  # Main boxplot
  p <- ggplot(est_df, aes(x = center, y = estimate, fill = center)) +
    geom_boxplot(outlier.size = 0.3, show.legend = FALSE) +
    geom_point(data = true_df_expanded,
               aes(x = center, y = true_value),
               color = "gold", size = 1, inherit.aes = FALSE) +
    geom_hline(data = true_df,
               aes(yintercept = true_value),
               color = "gold", linewidth = 1,
               inherit.aes = FALSE,linetype = "dashed")+
    facet_wrap(~ param, scales = "free_y", ncol = 3, labeller = label_parsed) +
    labs(title = title_text,
         subtitle = subtitle,
         x = "Method", y = "Estimate") +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  
  
  return(p)
}

plot_2d_curve_lorenz <- function(df,
                                 axisnames = c('X', 'Y'),
                                 line_size = 0.2,
                                 palette = "viridis",
                                 hide_legend=TRUE) {
  if (!all(c("X1", "X2") %in% names(df))) {
    names(df) <- c('X1', 'X2')
  }
  
  df$idx <- seq_len(nrow(df))
  cols <- viridis::viridis(256, option = palette)
   
  
  p <- ggplot(df, aes(x = X1, y = X2, color = idx)) +
    geom_path(linewidth = line_size, lineend = "round") +
    scale_color_gradientn(colors = cols, name = "step") +
    labs(
      x = axisnames[1],
      y = axisnames[2]
    ) +
    theme_bw(base_size = 17)
      
    if (hide_legend) p <- p + guides(color = "none") 
  return(p)
}



# Kræver df med d+1 columns X1,X2,...,Xd
plot_variable_vs_time <- function(df,times, to_time = 50) {
  df$times <- times
  h <- times[2]-times[1]
  steps <- to_time / h
  
  df <- df[1:steps,]
  
  df_long <- df %>%
    pivot_longer(cols = starts_with("X"), 
                 names_to = "variable", 
                 values_to = "value")

  ggplot(df_long, aes(x = times, y = value, color = variable)) +
    geom_line() +
    labs(x = "Time", y = "Value", title = "Each Xi vs Time")
}

plot_steps <- function(df_path, EM_steps, df_steps, b,
                       title = "Linearization around fixpoint",
                       add_nullclines = TRUE,
                       xlim = c(-2, 2), n = 400, plane = c('X','Y')) {
  
  # Helper to ensure x,y columns exist
  ensure_xy <- function(df) {
    df <- as.data.frame(df)
    if (all(c("x", "y") %in% names(df))) return(df)
    if (all(c("X1", "X2") %in% names(df))) {
      names(df)[names(df) == "X1"] <- "x"
      names(df)[names(df) == "X2"] <- "y"
      return(df)
    }
    names(df)[1:2] <- c("x", "y")
    return(df)
  }
  
  df_path  <- ensure_xy(df_path)
  EM_steps <- ensure_xy(EM_steps)
  df_steps <- as.data.frame(df_steps)
  
  if (!"label" %in% names(df_steps)) stop("df_steps must contain a 'label' column (e.g. 'S1','S2','S3').")
  if (!all(c("x0_1","x0_2") %in% names(df_steps))) stop("df_steps must contain 'x0_1' and 'x0_2' columns.")
  
  # Prepare plot data
  path_df  <- transform(df_path, which = "Path realization")
  df_S3    <- transform(subset(df_steps, label == "S3"), which = "Distribution of Strang step")
  EM_df    <- transform(EM_steps, which = "True distribution (EM simulated)")
  init_pts <- unique(df_steps[, c("x0_1", "x0_2")])
  names(init_pts) <- c("x", "y"); init_pts$which <- "Initial points"
  
  # ----- handle b (centers) : allow single point or multiple -----
  # Accept: numeric length-2 c(x,y), data.frame/matrix with two columns, or named list/data.frame with x/y
  make_centers_df <- function(b) {
    if (is.null(b)) return(NULL)
    if (is.numeric(b) && length(b) == 2) {
      return(data.frame(x = b[1], y = b[2], which = "Center of linearization"))
    }
    bdf <- as.data.frame(b)
    # If columns named x and y exist, use them; otherwise use first two columns
    if (all(c("x","y") %in% names(bdf))) {
      centers <- bdf[, c("x","y")]
    } else if (all(c("X1","X2") %in% names(bdf))) {
      centers <- bdf[, c("X1","X2")]
      names(centers) <- c("x","y")
    } else {
      centers <- bdf[, 1:2]
      names(centers) <- c("x","y")
    }
    centers$which <- "Center of linearization"
    # drop NA rows (if any)
    centers <- centers[!is.na(centers$x) & !is.na(centers$y), , drop = FALSE]
    rownames(centers) <- NULL
    return(centers)
  }
  
  center_df <- make_centers_df(b)
  
  # color map
  color_map <- c(
    "Path realization" = "grey80",
    "Distribution of Strang step" = "purple",
    "True distribution (EM simulated)" = "gold",
    "Initial points" = "black",
    "Center of linearization" = "red",
    "Nullclines" = "darkgrey"
  )
  
  # base plot
  p <- ggplot() +
    geom_path(data = path_df, aes(x = x, y = y, color = which),
              linewidth = 0.3, alpha = 0.7) +
    geom_point(data = df_S3, aes(x = x, y = y, color = which), size = 0.1) +
    geom_point(data = EM_df, aes(x = x, y = y, color = which), size = 0.1) +
    geom_point(data = init_pts, aes(x = x, y = y, color = which), size = 1)
  
  # add centers if provided (can be multiple rows)
  if (!is.null(center_df) && nrow(center_df) > 0) {
    p <- p + geom_point(data = center_df, aes(x = x, y = y, color = which),
                        size = 2.5, shape = 21, stroke = 0.6)
  }
  
  # Add nullclines (separate dashed lines, same legend)
  if (add_nullclines) {
    # note: 'par' must exist in calling environment (keeps same behavior as original)
    eps   <- par[1]; alpha <- par[2]; gamma <- par[3]; beta <- par[4]
    xs <- seq(xlim[1], xlim[2], length.out = n)
    df_xnull <- data.frame(x = xs, y = xs - xs^3 + alpha, which = "Nullclines")
    df_ynull <- data.frame(x = xs, y = gamma * xs + beta, which = "Nullclines")
    
    p <- p +
      geom_line(data = df_xnull, aes(x = x, y = y, color = which, linetype = which),
                linewidth = 0.6, alpha = 0.8) +
      geom_line(data = df_ynull, aes(x = x, y = y, color = which, linetype = which),
                linewidth = 0.6, alpha = 0.8) +
      scale_linetype_manual(values = c("Nullclines" = "dashed"))
  }
  
  p <- p +
    scale_color_manual(
      name = "",
      values = color_map,
      breaks = c(
        "Initial points",
        "Center of linearization",
        "Distribution of Strang step",
        "True distribution (EM simulated)"
      )
    ) +
    labs(x = plane[1], y = plane[2], title = title) +
    theme_bw() +
    guides(
      color = guide_legend(override.aes = list(size = 3)),
      linetype = "none"
    )
  
  return(p)
}

plot_multiple_steps <- function(df_path, EM_steps,
                                df_steps_list, b_list,
                                titles = NULL,
                                ncol = NULL,  # number of columns in the layout; if NULL a sensible default is chosen
                                xlim = NULL, ylim = NULL,
                                legend_pos = "bottom",
                                legend_text_size=20) {
  
  # validate inputs
  if (!is.list(df_steps_list)) stop("df_steps_list must be a list of data.frames (one per linearization).")
  n_plots <- length(df_steps_list)
  
  # Helper to normalize one b-element (can be numeric length2, data.frame/matrix with 2 cols, NULL)
  normalize_b_element <- function(be) {
    if (is.null(be)) return(NULL)
    # numeric vector length 2 -> single center
    if (is.numeric(be) && length(be) == 2) return(as.numeric(be)[1:2])
    # data.frame or matrix -> return 2-col data.frame (rows = centers)
    if (is.data.frame(be) || is.matrix(be)) {
      bdf <- as.data.frame(be)
      if (ncol(bdf) < 2) stop("b element data.frame/matrix must have at least two columns (x,y).")
      bdf <- bdf[, 1:2, drop = FALSE]
      names(bdf) <- c("x", "y")
      # keep as data.frame (possibly multiple rows)
      return(bdf)
    }
    # list: try to convert to data.frame or numeric
    if (is.list(be)) {
      # if it looks like c(x=..., y=...) or list of numeric pairs
      # try to coerce to data.frame
      try_df <- try(as.data.frame(be), silent = TRUE)
      if (!inherits(try_df, "try-error") && ncol(try_df) >= 2) {
        df2 <- try_df[, 1:2, drop = FALSE]
        names(df2) <- c("x", "y")
        return(df2)
      }
      # if it's numeric-ish vector inside list
      if (length(unlist(be)) == 2 && all(sapply(unlist(be), is.numeric))) {
        return(as.numeric(unlist(be))[1:2])
      }
    }
    stop("Unsupported b_list element type. Each element must be NULL, numeric length-2, or a 2-column data.frame/matrix (rows = centers).")
  }
  
  # If b_list provided as matrix/data.frame with nrow == n_plots, treat each row as a single center
  if (is.matrix(b_list) || is.data.frame(b_list)) {
    if (nrow(b_list) != n_plots) {
      stop("When b_list is a matrix/data.frame it must have one row per df_steps_list element (nrow(b_list) == length(df_steps_list)).")
    }
    # convert rows to list of numeric pairs
    b_list <- lapply(seq_len(nrow(as.data.frame(b_list))), function(i) {
      row <- as.numeric(as.data.frame(b_list)[i, 1:2])
      row
    })
  } else {
    # if not a data.frame/matrix, coerce to list if needed
    if (!is.list(b_list)) {
      b_list <- as.list(b_list)
    }
    if (length(b_list) != n_plots) stop("Length of b_list must match length of df_steps_list.")
  }
  
  # Normalize each element (can become numeric length2 or a data.frame of centers or NULL)
  b_list_norm <- vector("list", n_plots)
  for (i in seq_len(n_plots)) {
    b_list_norm[[i]] <- normalize_b_element(b_list[[i]])
  }
  
  if (is.null(titles)) {
    titles <- paste("Linearization", seq_len(n_plots))
  } else {
    if (length(titles) != n_plots) stop("Length of titles must match number of items in df_steps_list.")
  }
  
  # sensible default for ncol
  if (is.null(ncol)) {
    ncol <- if (n_plots == 1) 1 else if (n_plots == 2) 2 else ceiling(sqrt(n_plots))
  }
  
  # build individual plots
  plots <- vector("list", n_plots)
  for (i in seq_len(n_plots)) {
    # pass b_list_norm[[i]] directly to plot_steps(); plot_steps handles single or multiple centers
    plots[[i]] <- plot_steps(df_path, EM_steps, df_steps_list[[i]], b_list_norm[[i]], title = titles[[i]])
    # apply axis limits if provided
    if (!is.null(xlim) || !is.null(ylim)) {
      plots[[i]] <- plots[[i]] + coord_cartesian(xlim = xlim, ylim = ylim)
    }
  }
  
  # combine with patchwork, collect legends and set legend position
  combined <- wrap_plots(plots, ncol = ncol) +
    plot_layout(guides = "collect") &
    theme(legend.position = legend_pos,
          legend.text  = element_text(size = legend_text_size))
  
  return(combined)
}


plot_top_k_b_choices <- function(paths, x0, b_results, par, h, k = 5,
                                 xlim = c(-2, 2), ylim = c(-0.2, 2.1)) {
  library(ggplot2)
  
  # Helper: ensure x,y column names
  ensure_xy <- function(df) {
    df <- as.data.frame(df)
    if (all(c("x", "y") %in% names(df))) return(df)
    if (all(c("X1", "X2") %in% names(df))) {
      names(df)[names(df) == "X1"] <- "x"
      names(df)[names(df) == "X2"] <- "y"
      return(df)
    }
    names(df)[1:2] <- c("x", "y")
    return(df)
  }
  
  # Convert x0 into a data frame with proper column names
  x0_df <- as.data.frame(t(x0))
  names(x0_df) <- c("x", "y")
  x0_df$which <- "Initial point"
  
  # Nullcline data
  xs <- seq(xlim[1], xlim[2], length.out = 400)
  bs <- seq(xlim[1], xlim[2], length.out = 400)
  df_xnull <- data.frame(x = xs, y = xs - xs^3 + par[2], which = "Nullclines")
  df_ynull <- data.frame(x = xs, y = par[3] * xs + par[4], which = "Nullclines")
  
  # Polynomial overlay using x = x0[1]
  x_val <- x0[1]
  y_val <- x0[2]
  eps <- par[1]
  alpha <- par[2]
  gamma <- par[3]
  beta <- par[4]
  roots <- polyroot(c(alpha-y_val,1,0,-1))
  if (y_val > -(1/sqrt(3))^3+1/sqrt(3)+alpha){
    x_tilde <- min(Re(roots[abs(Im(roots))<0.00001]))
  }else if (y_val < -(-1/sqrt(3))^3-1/sqrt(3)+alpha) {
    x_tilde <- max(Re(roots[abs(Im(roots))<0.00001]))
  }else{
    
    real_roots <- Re(roots[abs(Im(roots)) < 0.00001])
    x_tilde <- real_roots[which.min(abs(real_roots - x_val))]
  }
  x_tilde_min <- min(Re(roots[abs(Im(roots))<0.00001]))
  x_tilde_med <- median(Re(roots[abs(Im(roots))<0.00001]))
  x_tilde_max <- max(Re(roots[abs(Im(roots))<0.00001]))
  
  
  
  # n <- length(xs)
  # 
  # b2_vals <- numeric(n)
  # fh_vals <- numeric(n)
  # 
  # for (i in seq_along(xs)) {
  #   
  #   b1 <- xs[i]
  #   
  #   find_b2 <- function(b2){
  #     fh <- fh_rcpp(c(x_val,y_val), -0.07/2,par, center=c(b1,b2), method='custom')
  #     abs(-3*b1^2*fh[1] - b1 + 3*b1^3 + b2 + fh[1]^3 - alpha)
  #   }
  #   
  #   opt <- optimize(find_b2, interval = c(-2,2))
  #   
  #   b2_vals[i] <- opt$minimum
  #     }
  
  # OU <- function(x, b1, b2){
  #   (x - 3*b1^2*x - y_val - b1 + 3*b1^3 + b2) 
  # }
  # 
  # N <- function(x, b1, b2){
  #   (-x^3 + 3*b1^2*x - 3*b1^3 + b1 + 0.5 - b2) 
  # }
  # 
  # n <- length(xs)
  # 
  # b2_vals <- numeric(n)
  # 
  # for (i in seq_along(xs)) {
  #   
  #   b1 <- xs[i]
  #   
  #   find_b2 <- function(b2){
  #     
  #     f1 <- N(x_val, b1, b2) / 2
  #     f2 <- OU(f1, b1, b2)
  #     f3 <- N(f2, b1, b2) / 2
  #     
  #     Fh  <- N(x_val, b1, b2) + OU(x_val, b1, b2)
  #     
  #     abs(f1 - Fh)
  #   }
  #   
  #   opt <- optimize(find_b2, interval = c(-0.7,2))
  #   
  #   b2_vals[i] <- opt$minimum
  # }
  # sigma <- par[5]
  # bias<-function(h,b,x0,y_val){
  #   N <- function(x){(1/eps)*(-x^3+3*b^2*x+b-3*b^3+alpha-y_val)}
  #   dN <- function(x){(1/eps)*(-3*x^2+3*b^2)}
  #   d2N <- function(x){-6*x/eps}
  #   d3N <- function(x){-6/eps}
  #   A <- (1-3*b^2)/eps
  #   E <- function(x){ 
  #     (3/8)*A*dN(x)*N(x) + (1/4)*A^2*N(x) + (1/6)*A^3*(x-b) + (1/6)*d2N(x)*N(x)^2+(1/6)*dN(x)^2*N(x) +  (1/4)*dN(x)*A^2*(x-b) + (1/4)*d2N(x)*A^2*(x-b)^2 + (3/8)*d2N(x)*N(x)*A*(x-b) + (1/8)*dN(x)^2*A*(x-b)+ sigma^2 *( (1/4)*A*d2N(x) + (3/16)*d2N(x)*dN(x)+(1/4)*d3N(x)*A*(x-b)+(3/16)*d3N(x)*N(x) )   }
  #   true_E <- function(x){
  #     (1/3)*A*dN(x)*N(x) + (1/6)*A^2*N(x) + (1/6)*A^3*(x-b) + (1/6)*d2N(x)*N(x)^2+(1/6)*dN(x)^2*N(x) +  (1/3)*dN(x)*A^2*(x-b) + (1/6)*d2N(x)*A^2*(x-b)^2 + (1/3)*d2N(x)*N(x)*A*(x-b) + (1/6)*dN(x)^2*A*(x-b)+ sigma^2 *( (1/4)*A*d2N(x) + (1/4)*d2N(x)*dN(x)+(1/12)*d3N(x)*A*(x-b)+(1/12)*d3N(x)*N(x) ) }
  #   err <- function(x){h^3 * (E(x)-true_E(x))}
  #   return(abs(err(x0)))
  # }
  
  new_roots <- polyroot(c(-(gamma-1)*(x_tilde)+alpha-beta,0,-3*(x_tilde),2))
  b_tilde <- min(Re(new_roots[abs(Im(new_roots))<0.00001]))
  
  F1 <- (x_val - x_val^3 + alpha - y_val) / eps
  F1_prime <- (1-3*x_val^2)/eps
  
  N1 <- (x_val - x_val^3 + alpha - y_val) / eps
  N1_prime <- (1-3*x_val^2)/eps
  
  F2 <- gamma * x_val + beta - y_val
  
  # determinant
  detF <- (-1 + 3*x_val^2 + gamma) / eps
  J <- matrix(c(
    (1-3*x_val^2)/eps,           -1/eps,
    gamma, -1
  ), nrow = 2, byrow = TRUE) 
  # inverse Jacobian entries
  Jinv <- matrix(c(
    -1,           1/eps,
    -gamma, (1 - 3*x_val^2)/eps
  ), nrow = 2, byrow = TRUE) / detF
  eig <- eigen(J)$values
  # compute b
  Fvec <- c(F1, F2)
  speed <- sqrt(F1^2+F2^2)
  speed_mat <- diag(abs(Re(1/(1+eig))))
  
  p <- min(abs(Re(eig[1])),1)
  newton_2d_line_search <- function(x0, y0,
                                    alpha, beta, gamma, eps,
                                    tol = 1e-6,
                                    maxit = 50,
                                    rho = 0.5,
                                    c = 1e-6,
                                    verbose = FALSE) {
    
    x <- c(x0, y0)
    
    for (k in 1:maxit) {

      
      Fvec <- c(
        (x0 - x0^3 + alpha - y0) / eps,
        gamma * x0 + beta - y0
      )
      
      norm_F <- sum(Fvec^2)
      
      if (norm_F < tol) {
        if (verbose) cat("Converged in", k, "iterations\n")
        return(x)
      }
      
      J <- matrix(c(
        (1 - 3*x0^2)/eps,  -1/eps,
        gamma,                -1
      ), nrow = 2, byrow = TRUE)
      
      step <- tryCatch(
        solve(J, Fvec),
        error = function(e) {
          warning("Jacobian singular, fallback step")
          return(Fvec * 0.1)
        }
      )

      step_norm <- sqrt(sum(step^2))
      if (step_norm > 1) {
        step <- step / step_norm
      }
      
      alp <- 1
      
      while (TRUE) {
        x_new <- x - alp * step
        
        F_new <- c(
          (x_new[1] - x_new[1]^3 + alpha - x_new[2]) / eps,
          gamma * x_new[1] + beta - x_new[2]
        )
        
        if (sum(F_new^2) <= (1 - c * alp) * norm_F) {
          break
        }
        
        alp <- rho * alp
        
        if (alp < 1e-3) {
          # accept anyway (important!)
          break
        }
      }
      
      x <- x - alp * step
    }
    
    warning("Did not converge")
    return(x)
  }
  trust_region_step <- function(x, alpha, beta, gamma, eps, Delta) {
    
    x_val <- x[1]
    y_val <- x[2]
    
    # F(x)
    Fvec <- c(
      (x_val - x_val^3 + alpha - y_val) / eps,
      gamma * x_val + beta - y_val
    )
    
    # Jacobian
    J <- matrix(c(
      (1 - 3*x_val^2)/eps,  -1/eps,
      gamma,                -1
    ), nrow = 2, byrow = TRUE)
    
    # Newton step
    step <- tryCatch(
      solve(J, Fvec),
      error = function(e) Fvec * 0.1  # fallback
    )
    
    # Trust region clipping
    step_norm <- sqrt(sum(step^2))
    
    if (step_norm > Delta) {
      step <- (Delta / step_norm) * step
    }
    
    x_new <- x - step
    
    return(x_new)
  }
  trust_region_adaptive <- function(x, alpha, beta, gamma, eps, Delta) {
    
    F_fun <- function(x) {
      c(
        (x[1] - x[1]^3 + alpha - x[2]) / eps,
        gamma * x[1] + beta - x[2]
      )
    }
    
    J_fun <- function(x) {
      matrix(c(
        (1 - 3*x[1]^2)/eps,  -1/eps,
        gamma,               -1
      ), 2, 2, byrow = TRUE)
    }
    
    Fvec <- F_fun(x)
    J <- J_fun(x)
    
    step <- tryCatch(solve(J, Fvec), error = function(e) Fvec * 0.1)
    
    # clip step
    norm_s <- sqrt(sum(step^2))
    if (norm_s > Delta) {
      step <- (Delta / norm_s) * step
    }
    
    x_new <- x - step
    
    # compute rho
    pred <- sum(Fvec^2) - sum((Fvec - J %*% step)^2)
    actual <- sum(Fvec^2) - sum(F_fun(x_new)^2)
    
    rho <- actual / (pred + 1e-10)
    
    # update Delta
    if (rho < 0.25) {
      Delta <- 0.25 * Delta
    } else if (rho > 0.75) {
      Delta <- min(2 * Delta, 10)
    }
    
    return(list(x = x_new, Delta = Delta))
  }
  speed <- sqrt(sum(Fvec^2))
  jac_size <- norm(J, type = "2")
  
  Delta <- min(
    speed,
    1 / (jac_size + 1e-8)
  )
  #b <- c(x_val - p*F1/F1_prime,y_val)
  b_newt <- newton_2d_line_search(x_val,y_val,0.5,1.4,1.5,0.05)
  b_trust <- trust_region_adaptive(c(x_val,y_val),0.5,1.4,1.5,0.05,Delta)$x
  Fb1 <- (b_trust[1] - b_trust[1]^3 + alpha - b_trust[2]) / eps
  Fb2 <- gamma * b_trust[1] + beta - b_trust[2]
  
  # determinant
  detFb <- (-1 + 3*b_trust[1]^2 + gamma) / eps
  
  # inverse Jacobian entries
  Jbinv <- matrix(c(
    -1,           1/eps,
    -gamma, (1 - 3*b_trust[1]^2)/eps
  ), nrow = 2, byrow = TRUE) / detF
  
  # compute b
  Fbvec <- c(Fb1, Fb2)
  speedb <- sqrt(Fb1^2+Fb2^2)
  b_tilde <- c(b_trust[1], b_trust[2]) - Jbinv %*% Fbvec
  
  
  
  
  #N_roots <- Re(polyroot(c(x_val-0.07*x_val^3/(2*0.05),-1,0.07*3*x_val/(2*0.05),-0.07/0.05)))
  N_roots <- Re(polyroot(c(-0.5*x_val-0.07*x_val^3/(4*0.05),-1,0.07*3*x_val/(4*0.05),-0.07/(2*0.05))))
  
  df_poly <- data.frame(
    b1 = bs,
    # poly_val1 = -x_val^3 + 3*bs^2*x_val - 3*bs^3 + alpha + bs, # N=0
    # poly_val2 = (3*bs^2-1)*x_tilde+bs-3*bs^3+y_val,
    # poly_val3 = (3*bs^2-1)*x_val+bs-3*bs^3+y_val, # N=F
    # poly_val4 = beta-gamma*(x_val-x_tilde-bs)
    poly_val5 = -3*y_val*bs^2+gamma*bs-(gamma-1)*y_val+beta
  )
  
  # Top k choices
  top_k <- b_results[order(b_results$error), ][1:k, ]
  
  # Build the plot
  p <- ggplot() +
    geom_path(data = ensure_xy(paths[[1]]),
              aes(x = x, y = y),
              color = "grey80", linewidth = 0.4) +
    
    # Nullclines
    geom_line(data = df_xnull, aes(x = x, y = y, linetype = which),
              color = "darkgrey",
              show.legend = FALSE) +
    geom_line(data = df_ynull, aes(x = x, y = y, linetype = which),
              color = "darkgrey",
              show.legend = FALSE) +
    
    # Polynomial overlay
    # geom_line(data = df_poly, aes(x = b1, y = poly_val2),
    #           color = "green", linewidth = 0.6) +
    # geom_line(data = df_poly, aes(x = b1, y = poly_val1),
    #           color = "purple", linewidth = 0.6) +
    # geom_line(data = df_poly, aes(x = b1, y = poly_val3),
    #           color = "red", linewidth = 0.6) +
    # geom_line(data = df_poly, aes(x = b1, y = poly_val4),
    #           color = "red", linewidth = 0.6) +
    # geom_line(data = df_poly, aes(x = b1, y = poly_val5),
    #           color = "red", linewidth = 0.6) +
    # geom_vline(xintercept = b_tilde,color = "red", linewidth = 0.6) +
    geom_vline(xintercept = b_trust[1],color = "red", linewidth = 0.6) +
    geom_hline(yintercept = b_trust[2],color = "red", linewidth = 0.6) +
    # geom_vline(xintercept = b_tilde[1],color = "blue", linewidth = 0.6) +
    # geom_hline(yintercept = b_tilde[2],color = "blue", linewidth = 0.6) +
    #geom_vline(xintercept = N_roots[1],color = "red", linewidth = 0.6) +
    #geom_vline(xintercept = N_roots[2],color = "red", linewidth = 0.6) +
    #geom_vline(xintercept = N_roots[3],color = "red", linewidth = 0.6) +

    # Top k points
    geom_point(data = top_k,
               aes(x = b1, y = b2, color = error),
               shape = 4, size = 4) +
    
    # Initial point
    geom_point(data = x0_df,
               aes(x = x, y = y, shape = which),
               color = "black", size = 2) +
    
    scale_color_gradientn(colors = c("yellow", "red"),
                          name = "Prediction Error") +
    scale_shape_manual(name = "", values = c("Initial point" = 19)) +
    scale_linetype_manual(name = "", values = c("Nullclines" = "dashed")) +
    
    labs(
      #title = paste0("Initial point x = (",round(x0[1], 2), ", ", round(x0[2], 2)),
      x = 'x',
      y = 'y'
    ) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme_bw()
  
  return(p)
}


plot_top_k_x_choices <- function(paths, x_results, b, par, h, k = 5,
                                 xlim = c(-2, 2), ylim = c(-0.2, 2.1)) {
  library(ggplot2)
  
  ensure_xy <- function(df) {
    df <- as.data.frame(df)
    if (all(c("x", "y") %in% names(df))) return(df)
    if (all(c("X1", "X2") %in% names(df))) {
      names(df)[names(df) == "X1"] <- "x"
      names(df)[names(df) == "X2"] <- "y"
      return(df)
    }
    names(df)[1:2] <- c("x", "y")
    return(df)
  }
  
  # Ensure all data are data.frames
  path_df <- as.data.frame(ensure_xy(paths[[1]]))
  x_results <- as.data.frame(x_results)
  top_k <- head(x_results[order(x_results$error), ], k)
  names(top_k)[1:2] <- c("x", "y")
  
  b_df <- data.frame(x = b[1], y = b[2], which = "Linearization point")
  
  # Nullclines
  xs <- seq(xlim[1], xlim[2], length.out = 400)
  df_xnull <- data.frame(x = xs, y = xs - xs^3 + par[2], which = "Nullclines")
  df_ynull <- data.frame(x = xs, y = par[3] * xs + par[4], which = "Nullclines")
  
  p <- ggplot() +
    geom_path(data = path_df, aes(x = x, y = y), color = "grey80", linewidth = 0.4) +
    geom_line(data = df_xnull, aes(x = x, y = y), color = "darkgrey", linetype = "dashed") +
    geom_line(data = df_ynull, aes(x = x, y = y), color = "darkgrey", linetype = "dashed") +
    geom_point(data = top_k, aes(x = x, y = y, color = error), shape = 4, size = 4) +
    geom_point(data = b_df, aes(x = x, y = y, shape = which), color = "black", size = 2) +
    scale_color_gradientn(colors = c("yellow", "red"), name = "Prediction Error") +
    scale_shape_manual(name = "", values = c("Linearization point" = 19)) +
    labs(x = "x", y = "y") +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme_bw()
  
  return(p)
}


# Kræver dataframe df med præcis 3 columns svarende til dimensioner
plot_3d_curve <- function(df,colorscheme = 'Jet'){
  fig <- plot_ly(
    x = df[,1], y = df[,2], z = df[,3],
    type = "scatter3d", mode = "lines+markers",
    line = list(width = 4,
                color = 1:nrow(df),       
                colorscale = colorscheme),     
    marker = list(size = 0.1,
                  color = 1:nrow(df),
                  colorscale = colorscheme,
                  showscale = TRUE)   
  )
  fig
}

plot_lorenz_without_background <- function(df) {
  
  t <- seq_len(nrow(df))  # time index
  
  fig <- plot_ly(
    x = df[,1],
    y = df[,2],
    z = df[,3],
    type = "scatter3d",
    mode = "lines",
    line = list(
      width = 4,
      color = t,
      colorscale = "Plasma",
      showscale = FALSE
    )
  ) %>%
    layout(
      scene = list(
        xaxis = list(visible = FALSE),
        yaxis = list(visible = FALSE),
        zaxis = list(visible = FALSE),
        bgcolor = "white"
      ),
      paper_bgcolor = "white",
      plot_bgcolor  = "white"
    )
  
  fig
}

plot_2d_curve_lorenz <- function(df,
                                 axisnames = c('X', 'Y'),
                                 line_size = 0.2,
                                 palette = "viridis",
                                 hide_legend=TRUE) {
  if (!all(c("X1", "X2") %in% names(df))) {
    names(df) <- c('X1', 'X2')
  }
  
  df$idx <- seq_len(nrow(df))
  cols <- viridis::viridis(256, option = palette)
  
  
  p <- ggplot(df, aes(x = X1, y = X2, color = idx)) +
    geom_path(linewidth = line_size, lineend = "round") +
    scale_color_gradientn(colors = cols, name = "step") +
    labs(
      x = axisnames[1],
      y = axisnames[2]
    ) +
    theme_bw(base_size = 17)
  
  if (hide_legend) p <- p + guides(color = "none") 
  return(p)
}
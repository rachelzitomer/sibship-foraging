plot_1d_likelihood_surface <- function(fit, simulations=NULL, scaling_factor=1, fontfamily="", sim_color="gray50",
                                       what=c("loglik", "floral_cover_on_capture_rate", "landscape_distance_on_capture_rate"))
{
  stopifnot(length(fit$mle) == 1)
  library(ggplot2)
  library(cowplot)
  library(dplyr)

  what <- match.arg(what)

  stopifnot(scaling_factor > 0)

  rescaling <- scales::trans_new(
    "rescaling", 
    function(x) x*scaling_factor, 
    function(x) x/scaling_factor
  )

  fitted_models <- data.frame(fit$fitted_models) %>%
    mutate(loglik=loglik - max(loglik)) %>%
    as.data.frame

  plt <- ggplot(fitted_models)  

  if (!is.null(simulations))
  {
    loglik_surfaces <- simulations$loglik_surfaces %>%
      group_by(seed) %>% 
      mutate(sim_loglik = sim_loglik - max(sim_loglik)) %>%
      as.data.frame

    plt <- plt + geom_line(
      data=loglik_surfaces, 
      aes_string(
        x=paste0("sim_", names(fit$mle)), 
        y=paste0("sim_", what), 
        group="seed"
      ), 
      color=sim_color, alpha=0.25, size=0.25
    )
  }

  plt <- plt +
    geom_vline(xintercept=fit$mle, color="black", lty=2) +
    annotate(geom="text", y=Inf, x=fit$mle, vjust=-0.1, color="black", 
             label="Estimate\n(MLE)", family=fontfamily, size=3) +
    geom_vline(xintercept=0, color=sim_color, lty=2) +
    annotate(geom="text", y=Inf, x=0.0, vjust=-0.1, color=sim_color, 
             label="Null model\n(no effect)", family=fontfamily, size=3) +
    geom_line(aes_string(x=names(fit$mle), y=what)) +
    geom_point(aes_string(x=names(fit$mle), y=what)) +
    theme_cowplot() +
    theme(text=element_text(family=fontfamily),
          plot.background = element_rect(fill="white"),
          plot.margin = grid::unit(c(0.10,0.05,0.05,0.05), "npc")) +
    coord_cartesian(clip="off") +
    scale_x_continuous(trans=rescaling) +
    ylab("Loglikelihood") +
    xlab(paste0("Effect of covariate on cell resistance (", names(fit$mle), ")"))

  return(plt)
}

plot_1d_sampling_distributions <- function(fit, null_simulations, parametric_bootstraps, scaling_factor=1, fontfamily="", num_bins=15, null_color="firebrick", boot_color="dodgerblue")
{
  stopifnot(length(fit$mle) == 1)
  library(ggplot2)
  library(cowplot)

  rescaling <- scales::trans_new(
    "rescaling", 
    function(x) x*scaling_factor, 
    function(x) x/scaling_factor
  )

  par_name <- names(fit$mle)
  null_quant <- quantile(null_simulations$MLEs[,par_name], c(0.025,0.975), na.rm=TRUE)
  boot_quant <- quantile(parametric_bootstraps$MLEs[,par_name], c(0.025,0.975), na.rm=TRUE)

  # Histograms
  # TODO: switch to density?
  ggplot() +
    geom_histogram(data=null_simulations$MLEs, aes_string(x=par_name), 
                   fill=null_color, group="null", alpha=0.25, bins=num_bins) +
    geom_histogram(data=parametric_bootstraps$MLEs, aes_string(x=par_name), 
                   fill=boot_color, group="boot", alpha=0.25, bins=num_bins) -> plt 

  # Get upper extent of histograms
  plt_data <- Reduce(rbind, ggplot_build(plt)$data)
  y_max <- max(plt_data$count)
  y_max_boot <- max(plt_data$count[plt_data$group=="boot"]) + y_max*0.05
  y_max_null <- max(plt_data$count[plt_data$group=="null"]) + y_max*0.05
  
  # Annotations, formatting
  plt <- plt +
    # bootstrap distribution
    annotate(geom="segment", y=y_max_boot, yend=y_max_boot, x=boot_quant[1], xend=boot_quant[2], 
             color=boot_color, arrow=grid::arrow(ends="both", angle=90, length=grid::unit(0.01,"npc"))) +
    annotate(geom="text", y=y_max_boot, x=mean(boot_quant), family=fontfamily, vjust=-0.5, 
             label="95% CI (bootstrap)", color=boot_color) +
    # MLE
    annotate(geom="segment", y=-Inf, yend=y_max_boot, x=fit$mle, 
             xend=fit$mle, color="black", lty=2) +
    annotate(geom="point", y=y_max_boot, x=fit$mle, color="black") +
    annotate(geom="text", y=-Inf, x=fit$mle, label="MLE", 
             vjust=1.5, color="black", family=fontfamily) +
    # null distribution
    annotate(geom="segment", y=y_max_null, yend=y_max_null, x=null_quant[1], xend=null_quant[2], 
             color=null_color, arrow=grid::arrow(ends="both", angle=90, length=grid::unit(0.01,"npc"))) +
    annotate(geom="text", y=y_max_null, x=0, family=fontfamily, vjust=-0.5, 
             label="Null 95% quantiles", color=null_color) +
    annotate(geom="segment", y=-Inf, yend=y_max_null, x=0, xend=0, color=null_color, lty=2) +
    annotate(geom="point", y=y_max_null, x=0, color=null_color) +
    # formatting
    theme_cowplot() +
    theme(text=element_text(family=fontfamily),
          plot.background = element_rect(fill="white"),
          plot.margin = grid::unit(c(0.10,0.05,0.05,0.05), "npc")) +
    coord_cartesian(clip="off") +
    scale_x_continuous(trans=rescaling) +
    #xlim(min(fit$parameter_grid), max(fit$parameter_grid)) +
    ylab("Number of simulations") +
    xlab(paste0("Effect of covariate on cell resistance (", names(fit$mle), ")"))

  return(plt)
}

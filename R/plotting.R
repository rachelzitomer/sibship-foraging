plot_1d_likelihood_surface <- function(fit, bootstraps=NULL, scaling_factor=1, fontfamily="Garamond", 
                                       what=c("loglik", "floral_cover_on_capture_rate", "landscape_distance_on_capture_rate"))
{
  stopifnot(length(fit$mle) == 1)
  library(ggplot2)
  library(cowplot)

  what <- match.arg(what)

  stopifnot(scaling_factor > 0)

  rescaling <- scales::trans_new(
    "rescaling", 
    function(x) x*scaling_factor, 
    function(x) x/scaling_factor
  )

  ggplot(fit$fitted_models) + 
    geom_vline(xintercept=fit$mle, color="black", lty=2) +
    annotate(geom="text", y=Inf, x=fitted_model_mle, vjust=-0.1, color="black", 
             label="Estimate\n(MLE)", family=fontfamily, size=3) +
    geom_vline(xintercept=0, color="firebrick", lty=2) +
    annotate(geom="text", y=Inf, x=0.0, vjust=-0.1, color="firebrick", 
             label="Null model\n(no effect)", family=fontfamily, size=3) +
    geom_line(aes_string(x=names(mle), y=what)) +
    geom_point(aes_string(x=names(mle), y=what)) +
    theme_cowplot() +
    theme(text=element_text(family=fontfamily),
          plot.background = element_rect(fill="white"),
          plot.margin = grid::unit(c(0.10,0.05,0.05,0.05), "npc")) +
    coord_cartesian(clip="off") +
    scale_x_continuous(trans=rescaling) +
    ylab("Loglikelihood") +
    xlab(paste0("Effect of covariate on cell resistance (", names(mle), ")")) -> plt

  if (!is.null(bootstraps))
  {
    plt <- plt + geom_line(
      data=bootstraps$loglik_surfaces, 
      aes_string(x=names(mle), y=what, group="seed"), 
      color="firebrick", alpha=0.25, size=0.25
    )
  }

  return(plt)
}

plot_1d_sampling_distributions <- function(fit, null_simulations, parametric_bootstraps, scaling_factor=1, fontfamily="Garamond", num_bins=15)
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
                   fill="firebrick", group="null", alpha=0.25, bins=num_bins) +
    geom_histogram(data=parametric_bootstraps$MLEs, aes_string(x=par_name), 
                   fill="dodgerblue", group="boot", alpha=0.25, bins=num_bins) -> plt 

  # Get upper extent of histograms
  plt_data <- ggplot_build(plt)$data[[1]]
  y_max <- max(plt_data$ncount)
  y_max_boot <- max(plt_data$ncount[plt_data$group=="boot"]) + y_max*0.05
  y_max_null <- max(plt_data$ncount[plt_data$group=="null"]) + y_max*0.05
  
  # Annotations, formatting
  plt <- plt +
    # bootstrap distribution
    annotate(geom="segment", y=y_max_boot, yend=y_max_boot, x=boot_quant[1], xend=boot_quant[2], 
             color="dodgerblue", arrow=grid::arrow(ends="both", angle=90, length=grid::unit(0.01,"npc"))) +
    annotate(geom="text", y=y_max_boot, x=mean(boot_quant), family=fontfamily, vjust=-0.5, 
             label="95% CI (bootstrap)", color="dodgerblue") +
    # MLE
    annotate(geom="segment", y=-Inf, yend=y_max_boot, x=fit$mle, 
             xend=fit$mle, color="black", lty=2) +
    annotate(geom="point", y=y_max_boot, x=fit$mle, color="black") +
    annotate(geom="text", y=-Inf, x=fit$mle, label="MLE", 
             vjust=1.5, color="black", family=fontfamily) +
    # null distribution
    annotate(geom="segment", y=y_max_null, yend=y_max_null, x=null_quant[1], xend=null_quant[2], 
             color="firebrick", arrow=grid::arrow(ends="both", angle=90, length=grid::unit(0.01,"npc"))) +
    annotate(geom="text", y=y_max_null, x=0, family=fontfamily, vjust=-0.5, 
             label="Null 95% quantiles", color="firebrick") +
    annotate(geom="segment", y=-Inf, yend=y_max_null, x=0, xend=0, color="firebrick", lty=2) +
    annotate(geom="point", y=y_max_null, x=0, color="firebrick") +
    # formatting
    theme_cowplot() +
    theme(text=element_text(family=fontfamily),
          plot.background = element_rect(fill="white"),
          plot.margin = grid::unit(c(0.10,0.05,0.05,0.05), "npc")) +
    coord_cartesian(clip="off") +
    scale_x_continuous(trans=rescaling) +
    #xlim(min(fit$parameter_grid), max(fit$parameter_grid)) +
    ylab("Number of simulations") +
    xlab(paste0("Effect of covariate on cell resistance (", names(mle), ")"))

  return(plt)
}

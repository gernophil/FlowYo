#' Load an FCS file
#'
#' This function loads an FCS file as dataframe and returns it.
#'
#' @param infile Path to the input FCS file
#' @return A dataframe of the FCS infile
#' @export
read_fcs <- function(fcs_path) {
  fcs_data_raw <- flowCore::read.FCSread.FCS(fcs_path,
                                             truncate_max_range = FALSE)

  fcs_data <- as.data.frame(fcs_data_raw@exprs)

  return(fcs_data)
}


#' Gate an FCS file
#'
#' This function takes and FCS dataframe from read_fcs() and a gate and returns
#' the original dataframe with information about the gate, a gated dataframe and
#' a plot of the gating.
#'
#' @param dataframe from read_fcs()
#' @param string plot title
#' @param string x column name
#' @param string y column name
#' @param gate which can be a dataframe, a list or a number
#' @param string for the gate type
#' @param string for the filter type
#' @param logical if x axis should be logarithmic
#' @param logical if y axis should be logarithmic
#' @param numeric for the max x value of
#' @param numeric for the max y value of
#' @param string for the plot type
#' @return A list containing the original dataframe with information about the
#' gate, a gated dataframe and a plot of the gating
#' @export
gate_fcs <- function(fcs_data,
                     name = "Flow Yo plot",
                     x = "FSC-A",
                     y = "SSC-A",
                     filter_gate = data.frame(x = c(30000, 30000, 180000, 262143, 262143, 30000),
                                              y = c(0, 40000, 262143, 262143, 0, 0)),
                     gate_type = "polygon",
                     gate_name = "gate",
                     filter_type = "inside",
                     # additional_gates = NULL,
                     log_x = FALSE,
                     log_y = FALSE,
                     max_x = 262143,
                     max_y = 262143,
                     plot_type = "hexagon") {

  if (gate_type == "polygon") {

    # implement check that first item of gate equals last item!
    fcs_data <-
      fcs_data |>
      dplyr::mutate(!!gate_name := ifelse(secr::pointsInPolygon(fcs_data[, c(x, y)], filter_gate, logical = TRUE),
                                          "inside",
                                          "outside"))

    if (filter_type == "inside") {

      fcs_data_gated <-
        fcs_data |>
        dplyr::filter(!!as.name(gate_name) == "inside")

    } else if (filter_type == "outside") {

      fcs_data_gated <-
        fcs_data |>
        dplyr::filter(!!as.name(gate_name) == "outside")

    } else {
      stop("Invalid filter type for gate type 'polygon'!")
    }

  } else if (gate_type == "threshold_y") {

    fcs_data <-
      fcs_data |>
      dplyr::mutate(!!gate_name := ifelse(!!as.name(y) > filter_gate,
                                          "higher",
                                          "lower_or_equal")) # if value equals threshold, where should it go to?

    if (filter_type == "higher") {

      fcs_data_gated <-
        fcs_data |>
        dplyr::filter(!!as.name(gate_name) == "higher")

    } else if (filter_type == "lower_or_equal") {

      fcs_data_gated <-
        fcs_data |>
        dplyr::filter(!!as.name(gate_name) == "lower_or_equal")

    } else {
      stop("Invalid filter type for gate type 'threshold_y'!")
    }

  } else if (gate_type == "threshold_x") {

    fcs_data <-
      fcs_data |>
      dplyr::mutate(!!gate_name := ifelse(!!as.name(x) > filter_gate,
                                          "bigger",
                                          "smaller_or_equal")) # if value equals threshold, where should it go to?

    if (filter_type == "bigger") {

      fcs_data_gated <-
        fcs_data |>
        dplyr::filter(!!as.name(gate_name) == "bigger")

    } else if (filter_type == "smaller_or_equal") {

      fcs_data_gated <-
        fcs_data |>
        dplyr::filter(!!as.name(gate_name) == "smaller_or_equal")

    } else {
      stop("Invalid filter type for gate type 'threshold_x'!")
    }

  } else if (gate_type == "quadrant") {

    fcs_data <-
      fcs_data |>
      dplyr::mutate(!!(paste0(gate_name, "_x")) := ifelse(!!as.name(x) > filter_gate$x,
                                                          "bigger",
                                                          "smaller_or_equal"),  # if value equals threshold, where should it go to?
                    !!(paste0(gate_name, "_y")) := ifelse(!!as.name(y) > filter_gate$y,
                                                          "higher",
                                                          "lower_or_equal")) # if value equals threshold, where should it go to?

    if (filter_type == "top_left") {

      fcs_data_gated <-
        fcs_data |>
        dplyr::filter(!!as.name(paste0(gate_name, "_x")) == "smaller_or_equal" &
                        !!as.name(paste0(gate_name, "_y")) == "higher")

    } else if (filter_type == "top_right") {

      fcs_data_gated <-
        fcs_data |>
        dplyr::filter(!!as.name(paste0(gate_name, "_x")) == "bigger" &
                        !!as.name(paste0(gate_name, "_y")) == "higher")

    } else if (filter_type == "bottom_left") {

      fcs_data_gated <-
        fcs_data |>
        dplyr::filter(!!as.name(paste0(gate_name, "_x")) == "smaller_or_equal" &
                        !!as.name(paste0(gate_name, "_y")) == "lower_or_equal")

    } else if (filter_type == "bottom_right") {

      fcs_data_gated <-
        fcs_data |>
        dplyr::filter(!!as.name(paste0(gate_name, "_x")) == "bigger" &
                        !!as.name(paste0(gate_name, "_y")) == "lower_or_equal")

    } else {
      stop("Invalid filter type for gate type 'quadrant'!")
    }

  } else {
    stop("Invalid gate type!")
  }

  fcs_plot <- ggplot2::ggplot() +
    ggplot2::xlab(x) +
    ggplot2::ylab(y) +
    ggplot2::ggtitle(name) +
    ostendorflab::theme_custom2()

  if (plot_type == "hexagon") {

    fcs_plot <-
      fcs_plot +
      ggplot2::geom_hex(aes(fcs_data[[x]], fcs_data[[y]]), bins = 256) +
      ggplot2::scale_fill_gradientn(colours = rev(brewer.pal(11, "Spectral")))

  } else if (plot_type == "density") {

    fcs_plot <-
      fcs_plot +
      ggplot2::geom_point(aes(fcs_data[[x]], fcs_data[[y]]),
                          size = .1,
                          color = "lightgrey") +
      ggplot2::geom_density_2d(aes(fcs_data[[x]], fcs_data[[y]]),
                               color = "black")

  } else {
    stop("Invalid plot type!")
  }

  if (log_x == TRUE) {
    fcs_plot <-
      fcs_plot +
      # scale_x_log10(limits = c(1, max(fcs_data[[x]])))
      ggplot2::scale_x_log10(limits = c(1, max_x))
  } else {
    fcs_plot <-
      fcs_plot +
      # scale_x_continuous(limits = c(0, max(fcs_data[[x]])))
      ggplot2::scale_x_continuous(limits = c(0, max_x))
  }

  if (log_y == TRUE) {
    fcs_plot <-
      fcs_plot +
      # scale_y_log10(limits = c(1, max(fcs_data[[y]])))
      scale_y_log10(limits = c(1, max_y))
  } else {
    fcs_plot <-
      fcs_plot +
      # scale_y_continuous(limits = c(0, max(fcs_data[[y]])))
      ggplot2::scale_y_continuous(limits = c(0, max_y))
  }

  if (gate_type == "polygon") {

    fcs_plot <-
      fcs_plot +
      ggplot2::geom_polygon(aes(filter_gate$x, filter_gate$y), fill = NA, color = "black") +
      gggplot2::eom_text(aes(max_x * 0.9,
                             max_y * 0.9,
                             label = paste0("inside: ",
                                            round(100 * sum(fcs_data[[gate_name]] == "inside")/nrow(fcs_data), 0), "%",
                                            "\n",
                                            "outside: ",
                                            round(100 * sum(fcs_data[[gate_name]] == "outside")/nrow(fcs_data), 0), "%")),
                         color = "black",
                         hjust = 1)

  } else if (gate_type == "threshold_y") {

    fcs_plot <-
      fcs_plot +
      ggplot2::geom_hline(yintercept = filter_gate) +
      ggplot2::geom_text(aes(max_x * 0.9,
                             max_y * 0.9,
                             label = paste0("higher: ",
                                            round(100 * sum(fcs_data[[gate_name]] == "higher")/nrow(fcs_data), 0), "%",
                                            "\n",
                                            "lower_or_equal: ",
                                            round(100 * sum(fcs_data[[gate_name]] == "lower_or_equal")/nrow(fcs_data), 0), "%")),
                         color = "black",
                         hjust = 1)

  } else if (gate_type == "threshold_x") {

    fcs_plot <-
      fcs_plot +
      ggplot2::geom_vline(xintercept = filter_gate) +
      ggplot2::geom_text(aes(max_x * 0.9,
                             max_y * 0.9,
                             label = paste0("bigger: ",
                                            round(100 * sum(fcs_data[[gate_name]] == "bigger")/nrow(fcs_data), 0), "%",
                                            "\n",
                                            "smaller_or_equal: ",
                                            round(100 * sum(fcs_data[[gate_name]] == "smaller_or_equal")/nrow(fcs_data), 0), "%")),
                         color = "black",
                         hjust = 1)

  } else if (gate_type == "quadrant") {

    fcs_plot <-
      fcs_plot +
      ggplot2::geom_hline(yintercept = filter_gate$y) +
      ggplot2::geom_vline(xintercept = filter_gate$x) +
      ggplot2::geom_text(aes(max_x * 0.9,
                             max_y * 0.7,
                             label = paste0("top_left: ",
                                            round(100 * sum(fcs_data[[paste0(gate_name, "_x")]] == "smaller_or_equal" &
                                                              fcs_data[[paste0(gate_name, "_y")]] == "higher")/nrow(fcs_data), 0), "%",
                                            "\n",
                                            "top_right: ",
                                            round(100 * sum(fcs_data[[paste0(gate_name, "_x")]] == "bigger" &
                                                              fcs_data[[paste0(gate_name, "_y")]] == "higher")/nrow(fcs_data), 0), "%",
                                            "\n",
                                            "bottom_left: ",
                                            round(100 * sum(fcs_data[[paste0(gate_name, "_x")]] == "smaller_or_equal" &
                                                              fcs_data[[paste0(gate_name, "_y")]] == "lower_or_equal")/nrow(fcs_data), 0), "%",
                                            "\n",
                                            "bottom_right: ",
                                            round(100 * sum(fcs_data[[paste0(gate_name, "_x")]] == "bigger" &
                                                              fcs_data[[paste0(gate_name, "_y")]] == "lower_or_equal")/nrow(fcs_data), 0), "%")),
                         color = "black",
                         hjust = 1)

  }

  print(fcs_plot)

  return(list("ungated" = fcs_data,
              "gated" = fcs_data_gated,
              "plot" = fcs_plot))
}


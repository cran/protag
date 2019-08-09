##### DEFINE FUNCTION FOR SPECTRA LIST PLOT  -<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>-----
#'
#'
#' tag.spectra.listplot
#'
#'
#' This function takes the output dataset from tag.search, and draw using ggplot2 the centroid mass spectra displayed in a listed manner. Peaks from the same "pair" (with designated m/z difference) are highlighted in differentiating colors, distinguished away from peaks of the "match" (with the same m/z) and the "mismatch" (neither of the prior two cases).
#'
#'
#' This function is designed for comparison of multiple mass spectra. In case of comparison of two mass spectra, it is recommended to use tag.spectra.butterflyplot for the highest annotation clarity.
#'
#'
#' @param search.output.list the output list from function tag.search
#' @param show.peak.pair if TRUE, show the paired peaks
#' @param show.peak.match if TRUE, show the matched peaks
#' @param show.peak.mismatch if TRUE, show the mismatched peaks
#' @param show.annotation.pair if TRUE, show the m/z annotations for the paired peaks
#' @param show.annotation.match if TRUE, show the m/z annotations for the mathced peaks
#' @param show.annotation.mismatch if TRUE, show the m/z annotations for the mismatched peaks
#'
#' @param size.peak.pair adjust the peak width of the paired peaks. All size.xxx arguments take a numeric value, same functionality as line width or text size control in ggplot2
#' @param size.peak.match adjust the peak width of the matched peaks
#' @param size.peak.mismatch adjust the peak width of the mismatched peaks
#' @param size.divider adjust divider width
#' @param size.annotation.pair adjust the m/z annotation text size for the paired peaks
#' @param size.annotation.match adjust the m/z annotation text size for the matched peaks
#' @param size.annotation.mismatch adjust the m/z annotation text size for the mismatched peaks
#' @param size.groupname adjust the text size for groupnames (e.g., "control", "experiment1", "experiment2", etc.).
#'
#' @param alpha.peak.pair adjust the transparency of the paired peaks. All alpha.xxx arguments take a numeric value [0,1]
#' @param alpha.peak.match adjust the transparency of the matched peaks
#' @param alpha.peak.mismatch adjust the transparency of the mismatched peaks
#' @param alpha.annotation.pair adjust the transparency of the m/z annotations for the paired peaks
#' @param alpha.annotation.match adjust the transparency of the m/z annotations for the matched peaks
#' @param alpha.annotation.mismatch adjust the transparency of the m/z annotations for the mismatched peaks
#'
#' @param color.pair control the color for the paired peaks and the associated m/z annotations. Each pair will be of the same color, and different pairs of differentiating colors. In case of multiple mass shifts being of interest within a pair, e.g., delta = c(14, 28, 56), then peaks with m/z difference of either 14, 28 or 56, all belonging to the same pair, will be of the same color.  Apart from the default color set, users could otherwise choose color from RColorBrewer palettes, e.g., color.pair = "Set1", or color.pair = "Blues".
#'
#' Colors for peaks (paired, matched, and mismatched) and the respectively associated annotations are designed to be of the same set of color for maximum clarity.
#'
#' @param color.match control the color for the matched peaks with the associated m/z annotations, with default in "black". Users may otherwise reset to different colors, e.g., color.match = "firebrick". As the matched peaks and mismatched peaks are usually of less research interest than paired peaks, the matched and mismatched peaks are respectively designed to be of monocolor.
#' @param color.mismatch control the color for the mismatched peaks with the associated m/z annotations, with default in "black".
#' @param color.groupname control the color for the groupnames, with default in "black".
#' @param color.divider control the color of the central divider
#' @param angle.annotation adjust the angle for the m/z annotations, taking a numeric value. This argument is useful to avoid annotation overlap, and is particularly handy when the plot is reoriented with coord_flip().
#' @param angle.groupname adjust the angle of the groupnames.
#' @param gap.groupname adjust the horizontal position of groupnames. A positive numeric value adjusts the distance between groupnames and the left bound of the mass spectra; negative values shifts the groupnames to the right side.
#' @param gap.annotation adjust the distance between m/z annotations and the top of the peak.
#' @param peak.height.shrink Taking a numeric value [0, 1], a small shrinking factor renders smaller peak height, and generates more space between peak and the central divider, leaving more space for annotations. This argument resolves overlap among annotations with upper-floor-residing peaks, a problem unique to listplot. Therefore, this argument is not used in the butterflyplot.
#' @return a ggplot2 plot.
#' @export
#' @examples
#' search.result <- tag.search(myoglobin, delta = c(14, 28), error.Da.pair = .3)
#' search.result
#' tag.spectra.listplot(search.result)

tag.spectra.listplot =
  function(search.output.list,
           show.peak.pair = TRUE, show.peak.match = TRUE, show.peak.mismatch = TRUE,
           show.annotation.pair = TRUE, show.annotation.match = TRUE, show.annotation.mismatch = FALSE,

           size.peak.pair = 2, size.peak.match = 1, size.peak.mismatch = .5,
           size.divider = .3,
           size.annotation.pair = NA, size.annotation.match = NA,
           size.annotation.mismatch = NA, size.groupname = NA,

           alpha.peak.pair = .8, alpha.peak.match = .5, alpha.peak.mismatch = .2,
           alpha.annotation.pair = .8, alpha.annotation.match = .5, alpha.annotation.mismatch = .2,

           color.pair = 1, color.match = "black", color.mismatch = "black",
           color.groupname = "black", color.divider = "black",
           angle.annotation = 90, angle.groupname = 90,

           gap.groupname = .02, gap.annotation = 0.15,

           peak.height.shrink = .7) {

    dataset = search.output.list[[1]]

    ## set up mass range
    mass.max = dataset$mass %>% max(); mass.max
    mass.min = dataset$mass %>% min(); mass.min
    mass.range = mass.max - mass.min; mass.range


    ## Set up color choice from RColorBrewer: palettes with maximum number of different values
    palettes.sequential = c("Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd")
    max.sequential = rep(9, times = palettes.sequential %>% length())

    palettes.qualitative = c("Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3")
    max.qualitative = c(8, 8, 12, 9, 8, 9, 8, 12)

    palettes.diverging = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral")
    max.diverging = rep(11, times = palettes.diverging %>% length())

    color.all = c(max.sequential, max.qualitative, max.diverging)
    names(color.all) = c(palettes.sequential, palettes.qualitative, palettes.diverging)


    ## Plot
    plot.spectra = dataset %>% ggplot(aes(x = mass, xend = mass,
                                          y = spectra.dividor,
                                          yend = spectra.dividor + intensity.scaled * peak.height.shrink)) +

      # central divider
      geom_segment(aes(x = mass.min - mass.range*.02, xend = mass.max + mass.range*.02,
                       y = spectra.dividor , yend = spectra.dividor),
                   size = size.divider, color = color.divider) +
      # aesthetic for peak and annotation unless otherwise specified; groups names and dividors apply their own aesthetic position
      # theme
      theme_minimal() +
      theme(panel.grid = element_blank(), legend.position = "None") +
      labs(x = "mass-charge ratio (m/z)", y = "relative intensity") +
      # groups names
      geom_text(aes(x = mass.min - mass.range * gap.groupname , y = spectra.dividor + 0.5, label = group),
                angle = angle.groupname, size = size.groupname, color = color.groupname)


    # PAIRED
    if (show.peak.pair == TRUE) { # peak
      plot.spectra = plot.spectra +
        geom_segment(data = dataset %>% filter(status == "pair"),
                     aes(color = factor(pair.tracker), alpha = alpha.peak.pair),
                     size = size.peak.pair)
    }
    if (show.annotation.pair == TRUE) { # m/z annotation
      plot.spectra = plot.spectra +
        geom_text(data = dataset %>% filter(status == "pair"),
                  aes(x = mass, y = spectra.dividor + intensity.scaled * peak.height.shrink + gap.annotation,
                      label = mass, color = factor(pair.tracker)),
                  angle = angle.annotation, size = size.annotation.pair, alpha = alpha.annotation.pair)
    }


    # MATCHED (plot them only when found)
    # Note: error will pop up when true match not found and the filtered match subset is empty!
    if (search.output.list[[2]] == "Found both paired peaks (mass differentiate by expected delta) and matched peaks (of the same mass)."){
      if (show.peak.match == TRUE) { # peak
        plot.spectra = plot.spectra +
          geom_segment(data = dataset %>% filter(status == "match"),
                       size = size.peak.match, alpha = alpha.peak.match, color = color.match)
      }

      if (show.annotation.match == TRUE) { # annotation
        plot.spectra = plot.spectra +
          geom_text(data = dataset %>% filter(status == "match"),
                    aes(x = mass, y = spectra.dividor + intensity.scaled * peak.height.shrink + gap.annotation,
                        label = mass),
                    angle = angle.annotation, size = size.annotation.match,
                    alpha = alpha.annotation.match, color = color.match)
      }
    }


    # MISMATCHED
    if ("mismatch" %in% dataset$status) { # draw only if there is true mismatches

      if (show.peak.mismatch == TRUE) {
        plot.spectra = plot.spectra +
          geom_segment(data = dataset %>% filter(status == "mismatch"),
                       size = size.peak.mismatch, alpha = alpha.peak.mismatch, color = color.mismatch)
        }

      if (show.annotation.mismatch == TRUE) {
        plot.spectra = plot.spectra +
          geom_text(data = dataset %>% filter(status == "mismatch"),
                    aes(x = mass, y = spectra.dividor + intensity.scaled * peak.height.shrink + gap.annotation,
                        label = mass),
                    angle = angle.annotation, size = size.annotation.mismatch,
                    alpha = alpha.annotation.mismatch, color = color.mismatch)
      }
    }


    ## Colors of the paired peak
    if (color.pair != 1) {
      func.color.pair = colorRampPalette( c(brewer.pal(color.all[color.pair], # max of different values in the given palette
                                                       color.pair)) ) # Divide given palette into gradients
      plot.spectra = plot.spectra +
        scale_color_manual(values = func.color.pair(n = dataset$pair.tracker %>% n_distinct()) %>% sample() ) # no. of gradient steps being no. of the label pairs
    }



    # Using the correct display format?
    all.groups.no = dataset$group %>% n_distinct()
    if (all.groups.no == 2){
      message("\nLooks good! \n\nWhen only two spectra are compared, a mirrored / butterfly plot may give higher clarity.",
              "\nDo you also want to try `tag.spectra.butterflyplot()`?\n")
      } else {
        message("\nGreat plot!\n")
        }

    return(plot.spectra)

  }
utils::globalVariables(c("search.output.list", "dataset", "spectra.dividor", "status", "pair.tracker", "colorRampPalette", "brewer.pal" ))



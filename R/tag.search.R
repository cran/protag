##### DEFINE FUNCTION FOR MASS LIST SEARCH  -<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>-----
#'
#' tag.search
#'
#'
#'This function takes a mass list dataset containing columns of "mass", "intensity" and "group" (contains the "control" observations) , and searches within specified error tolerance for "paired" peaks, "matched" peaks, and "mismatched" peaks. Mass spectra peaks with m/z diffrence being the designated variable "delta" (within error tolerance) are defined as a "pair", and peaks of the same m/z (within error tolerance) as "match"; otherwise defined as "mismatch".
#'
#'
#' @param dataset a tidy dataset containing the mass list. At least three numeric columns are required: "mass", "intensity" and "group". The "mass" refers to m/z values; "intensity" refers to peak height/area; "group" must contain the "control" observations.
#' @param delta a single numeric value, or a numeric vector when multiple m/z difference is of interest. The variable "delta" reflects the mass difference between the labelled proteins/peptides vs. the non-labelled (the control), caused by the chemically-labelling group.
#' @param error.Da.pair error tolerance for the paired peaks, in Dalton; default at 0.5.
#' @param error.Da.match error tolerance for the matched peaks, in Dalton; default at 0.5.
#'
#' @param error.ppm.pair error tolerance threshold for the paired peaks, in ppm. For paired peaks p and q, the tolerance threshold is defined as 0.5 * (p+q) * error.ppm.pair / 10^6. When the absolute difference between the measured vs. theoretical delta is lower than the error tolerance, then the associated two peaks are considered a pair. The default value of error.ppm.pair is Inf (positive infinite); that is, the error tolerance by default is controled by error.Da.pair. When error.ppm.pair is otherwise set, say at 100 (ppm), then the practical error tolerance value is the smallest of either the Dalton control or ppm control. When the ppm control is more desirable than Dalton control, consider setting error.Da.pair = Inf.
#' @param error.ppm.match error tolerance for the matched peaks, in ppm. Error tolerance control for matched peaks is similar to the case of paired peaks.
#' @param intens.log.transfrom default to FALSE. If set to TRUE, peak intensities will be logarithmically transformed. This is useful for displaying low-intensity peaks that would otherwise be overshadowed and less visible in the mass spectra.
#' @import dplyr
#' @import ggplot2
#' @return a tidy dataset, with the original input dataset augmented with additional columns. The content in the input dataset remain unchanged (though the display sequence may change).
#' @export
#' @examples
#' search.result <- tag.search(myoglobin, delta = c(14, 28), error.Da.pair = .3)
#' search.result
#' tag.spectra.listplot(search.result)


tag.search = function(dataset, delta = NA,
                      error.Da.pair = .5, error.Da.match = .5,
                      error.ppm.pair = Inf, error.ppm.match = Inf,
                      intens.log.transfrom = FALSE){

  ## Check column / level names
  if ( ! "mass" %in% colnames(dataset) ) { # have mass?
    stop("\nA variable/column named `mass` is required. Do you need to rename one of your columns?\n") }

  if ( ! "group" %in% colnames(dataset) ) { # have group?
    stop("\nA variable/column named `group` is required. Do you need to rename one of your columns?\n") }

  if ( ! "intensity" %in% colnames(dataset) ) { # have intensity?
    dataset = dataset %>% mutate(intensity = 100)
    cat("\n`Intensity` not found. Mass spectra peaks (if applicable) will be drawn as equal intensity.\n") }

  if ( ! "control" %in% dataset$group) { # have "control"?
    stop(("\n`control` is required in the `group` column.\nDo you need to rename some of the experimental levels?\n"))
  }


  ## Check delta input
  if ( sum(delta) %>% is.na() ){
    stop("\nMass shift input required.Use argument: delta = xxx (a numeric single value or vector) \n")}

  ## Remove duplicated rows in the 3 columns (to remove potential join errors)
  dataset = dataset[!duplicated(dataset %>% select(mass, group)), ]

  ## log transform intensity
  if (intens.log.transfrom == TRUE) {
    dataset = dataset %>% mutate(intensity = intensity %>% log10())
  }


  ## Names of experiments (NOT including control)
  label.names = dataset$group %>% unique();
  label.names = label.names[label.names != "control"]


  ## Set up control dataset
  d.control = dataset %>% filter(group == "control")


  ## Searching by pairwise m/z difference computation
  row.match = row.pair = 1 # row counter
  vctr.match.control = vctr.match.label = vctr.match.group =
    vctr.pair.control = vctr.pair.label = vctr.pair.group = c() # empty vectors recording matched and paired m/z


  for (i in 1:nrow(d.control)) { # iterate through m/z of the control

    for (x in 1:length(label.names)) { # iterate through different labelling experiments
      d.label.x = dataset %>% filter(group ==  label.names[x])


      for (j in 1:nrow(d.label.x)) { # iterate through m/z's of the labelled experiment

        diff = d.label.x$mass[j] - d.control$mass[i] # diff defined as label - control
        err.tolerance.match = mean( d.label.x$mass[j], d.control$mass[i]) * error.ppm.match / (10^6) # by ppm error
        err.tolerance.pair = mean( d.label.x$mass[j], d.control$mass[i] ) * error.ppm.pair / (10^6) # by ppm error

        # finding mathcing (equal)  m/z
        if ( ( abs(diff) <= err.tolerance.match)  & ( abs(diff) <= error.Da.match) ) { # diff: positive or negative for matched peaks; error by both Da and ppm control
          vctr.match.control[row.match] = d.control$mass[i] # control m/z
          vctr.match.label[row.match] = d.label.x$mass[j] # label m/z
          vctr.match.group[row.match] = label.names[x] # group names
          row.match = row.match + 1
        }


        for (p in 1:length(delta)){ # iterate though expected shift by tagged mass
          # finding paired (with expected m/z shift with labelled tag)
          diff.2 = abs( diff - delta[p] ) # diff as label - control, presumably positive for paired peaks

          if ( (diff.2 <= err.tolerance.pair) & (diff.2 <= error.Da.pair) ) {

            vctr.pair.control[row.pair] = d.control$mass[i] # control m/z
            vctr.pair.label[row.pair] = d.label.x$mass[j] # label m/z
            vctr.pair.group[row.pair] = label.names[x]
            row.pair = row.pair + 1
          }
        }
      }
    }
  }


  if (row.pair == 1) { # if no pair found, no continuous analysis and plot at all
    message("\nNo paired peaks (with expected delta) were found, and mass spectra cannot be drawn.",
    "\nCheck input mass shift (delta), or try increasing error tolerance for the pair search.\n")
    } else{

      ## Combine vectors into tidy dataset of the PAIRED m/z
      d.pair = tibble(pair.control = vctr.pair.control,
                      pair.label = vctr.pair.label,
                      group = vctr.pair.group)

      ## Combine vectors into dataset of the MATCHED m/z
      if (row.match == 1) {
        d.match = tibble(match.control = 1,
                         match.label = 1,
                         group = label.names[1])
      } else{
        d.match = tibble(match.control = vctr.match.control,
                         match.label = vctr.match.label,
                         group = vctr.match.group)
        }


      ## Set up MISMATCHED dataset
      found.control = c(d.match$match.control, d.pair$pair.control) # the matched and paired (collectively, "found")  in the control
      d.mismatch.control = d.control %>% filter(! mass %in% found.control) # mismatched from the control

      d.match.renamed = d.match %>% select(match.label, group) %>%
        rename(mass = match.label) %>% mutate(status = "match"); d.match.renamed
      d.pair.renamed = d.pair %>% select(pair.label, group) %>%
        rename(mass = pair.label) %>% mutate(status = "pair"); d.pair.renamed

      d.label.found = rbind(d.match.renamed, d.pair.renamed) # matched and paired from the labelled data
      d.mismatch.label = dataset %>% filter(group != "control") %>%
        anti_join(d.label.found, by = c("group", "mass")) # mismatched dataset for the labeled data

      d.mismatch.tidy = rbind(d.mismatch.control, d.mismatch.label) %>% # combine mismatch dataset from control & label
        mutate(status = "mismatch")



      ## Set up PAIRED dataset (augment with intensity and tidy up "d.pair")
      pair.counter = 1
      pair.tracker = c(1)

      # Sequence up pair numbers
      d.pair = d.pair %>% arrange(pair.control)

      if (nrow(d.pair) >=2 ){ # if =1, then pair.tracker = 1, as two lines above
        for (i in 2:nrow(d.pair)){
          if (d.pair$pair.control[i] == d.pair$pair.control[i-1]){
            pair.tracker[i] = pair.counter
          } else {
            pair.counter = pair.counter + 1
            pair.tracker[i] = pair.counter
          }
        }
      }


      d.pair = d.pair %>% mutate(pair.tracker = pair.tracker); d.pair

      # tidy up "d.pair", combine "pair.control" and "pair.label" into one column, with column name mass
      x1 = d.pair %>% select(-contains("label")) %>% rename(mass = pair.control) %>%
        mutate(group = "control") # change of table format requires to reset group into control

      x2 = d.pair %>% select(-contains("control")) %>% rename(mass = pair.label)
      d.pair.tidy = rbind(x1, x2); d.pair.tidy

      # augment with intensity
      d.pair.tidy = d.pair.tidy %>%
        left_join(dataset, by = c("mass", "group")) %>% mutate(status = "pair"); d.pair.tidy


      ## Set up MATCHED dataset (augment with intensity and tidy up "d.match")
      y1 = d.match %>% select(match.control, group) %>% rename(mass = match.control) %>%
        mutate(group = "control") # change of table format requires to reset group into control; like the case of x1 above
      y2 = d.match %>% select(match.label, group) %>% rename(mass = match.label); y2
      d.match.tidy = rbind(y1, y2) %>%
        mutate(status = "match") %>% left_join(dataset, by = c("mass", "group"))

      ## Combine all three tidy dataset
      d.tidy = rbind(d.mismatch.tidy, d.match.tidy) %>%
        mutate(pair.tracker = NA) %>% # match and mismatch has NA entries for pair.tracker column
        rbind(d.pair.tidy)


      ## Scale intensity
      d.tidy = d.tidy %>% group_by(group) %>%
        mutate(intensity.scaled = intensity / max(intensity, na.rm = TRUE)) %>% arrange(group)
        # na.rm =T: for faked matched peaks (when true match not found), the associated intensity is NA, and must be removed
        # This is particularly important to "sign" computation and control mirror pattern for butterfly plot

      ## Set up plot central dividor number

      xx = 1:(d.tidy$group %>% unique() %>% length())
      names(xx) = d.tidy$group %>% unique(); xx
      d.tidy = d.tidy %>% mutate(spectra.dividor = xx[group] - 1) # d.tidy being the dataset input into list plot function; minus 1 so that the scale start from zero


      ## Set up sign-differentiating column, reserved for butterfly plot
      d.tidy = d.tidy %>%
        mutate(`intensity.scaled.+.-` = ifelse(group == "control", -intensity.scaled, intensity.scaled),
               sign = `intensity.scaled.+.-`/abs(`intensity.scaled.+.-`)) %>%
        distinct() # remove duplicated rows (the paired peaks in the control; if not removed deeper color when overlapped multiple times)


      ## Remove faked match (their intensities are NA)
      d.tidy = d.tidy[!d.tidy$intensity %>% is.na(), ] # if not remove faked matched peaks (two rows, control + label), as intensity as NA, these two rows will be removed when plotting and return a console message syaing "two rows are removed with missing values". This does not hamper script run at all, but could cause confusion to users. Thus, remove the two rows to be confusion-free.


      return(list(d.tidy,
                  ifelse(row.match == 1,
                         "Found paired peaks (mass differentiate by expected delta) , but not matched peaks (of the same mass).",
                         "Found both paired peaks (mass differentiate by expected delta) and matched peaks (of the same mass).")))

    }
}
utils::globalVariables(c("intensity", "group", "mass", "match.label", "pair.label", "pair.control", "match.control", "intensity.scaled", "intensity.scaled.+.-"))



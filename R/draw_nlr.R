#' draw_nlr
#'
#' Create ggplot2 object with protein chains from feature database
#'
#' @param data Dataframe of one or more rows with the following column names: 'type', 'description', 'start', 'end'. Must contain a minimum of one "CHAIN" as data$type.
#' @param label_regions Add text to the canonical NLR domains, standard = TRUE
#' @param label_domains Add text to the integrated domains, standard = TRUE
#' @param label_motifs Label the NB-ARC motifs, standard = FALSE
#' @param label_repeats Label the LRRs, only available for the RefPlantNLR dataset, standard = FALSE
#' @param PDF Add text to the canonical NLR domains, standard = FALSE
#'
#' @return A ggplot2 object either in the plot window or as an object.
#' @export
#'
#' @examples
#' draw_nlr(RefPlantNLR, label_regions = TRUE, label_domains = TRUE, label_motifs = TRUE, label_repeats = FALSE, PDF = "RefPlantNLR.pdf")

draw_nlr <- function(data,
                     label_regions = TRUE,
                     label_domains = TRUE,
                     label_motifs = FALSE,
                     label_repeats = FALSE,
                     PDF = FALSE){

  # Add a column order to order the data in the plot
  data1 <- dplyr::distinct(data, seqname)
  data1 <- dplyr::arrange(data1, desc(seqname))
  data1 <- dplyr::mutate(data1, order = dplyr::row_number())
  data <- dplyr::left_join(data1, data, by = "seqname")

  # Rename NB-ARC motifs in order of occurrence
  data <- dplyr::mutate(data, description = dplyr::case_when(type == "MOTIF" & description == "P-loop" ~ "1",
                                                             type == "MOTIF" & description == "RNBS-A" ~ "2",
                                                             type == "MOTIF" & description == "Kin-2" ~ "3",
                                                             type == "MOTIF" & description == "RNBS-B" ~ "4",
                                                             type == "MOTIF" & description == "RNBS-C" ~ "5",
                                                             type == "MOTIF" & description == "GLPL" ~ "6",
                                                             type == "MOTIF" & description == "NB-ARC motif 12" ~ "7",
                                                             type == "MOTIF" & description == "RNBS-D" ~ "8",
                                                             type == "MOTIF" & description == "linker" ~ "9",
                                                             type == "MOTIF" & description == "MHD" ~ "10",
                                                             TRUE ~ description))

  points = ggplot2::scale_fill_manual(values = c("R1" = "#FBB4AE",
                                                 "CC" = "#B3CDE3",
                                                 "RPW8" = "#CCEBC5",
                                                 "TIR" = "#DECBE4",
                                                 "PLOOP" = "#FFFFCC",
                                                 "NBARC" = "#FED9A6",
                                                 "LRR" = "#E5D8BD",
                                                 "CJID" = "#FDDAEC",
                                                 "OTHER" = "#666666"))

  start=end=NULL

  # Make the canvas
  p <- ggplot2::ggplot() +
    ggplot2::ylim(0.5, max(data$order)+0.5) +
    ggplot2::xlim(-max(data$end, na.rm=TRUE)*0.2,
                  max(data$end, na.rm=TRUE) + max(data$end, na.rm=TRUE)*0.1) +
    ggplot2::scale_x_continuous(breaks = seq(from = 0, to = 100000, by = 250)) + # Arbitrarily large number
    ggplot2::labs(x = "Amino acid number") + # label x-axis
    ggplot2::labs(y = "") + # label y-axis
    ggplot2::theme_bw(base_size = 20) +  # white background and change text size
    ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),
                   panel.grid.major=ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   legend.key.size = ggplot2::unit(0.5, 'cm'),
                   legend.title = ggplot2::element_text(size=10),
                   legend.justification = "right",
                   legend.text = ggplot2::element_text(size=8))

  # Add chains
  p <- p + ggplot2::geom_rect(data = data[data$type == "CHAIN",],
                              mapping=ggplot2::aes(xmin=start,
                                                   xmax=end,
                                                   ymin=order-0.2,
                                                   ymax=order+0.2),
                              fill = "grey",
                              size = 0.5)

  # Add chain labels
  p <- p + ggplot2::annotate("text", x = -10,
                             y = data[data$type == "CHAIN",]$order,
                             label = data[data$type == "CHAIN",]$seqname,
                             hjust = 1,
                             size = 4)

  # Draw regions except for NBARC
  if(nrow(dplyr::filter(data, type == "REGION" & !description %in% c("NBARC", "PLOOP"))) != 0){
    p <- p + ggplot2::geom_rect(data= data[data$type == "REGION" & !data$description  %in% c("NBARC", "PLOOP"),],
                                mapping=ggplot2::aes(xmin=start,
                                                     xmax=end,
                                                     ymin=order-0.3,
                                                     ymax=order+0.3,
                                                     fill=description),
                                show.legend = TRUE) +
      points
  }

  # Draw regions: NBARC. This ensures the NBARC domain is plotted on top in case the boundaries overlap with other domains
  if(nrow(dplyr::filter(data, type == "REGION" & description %in% c("NBARC", "PLOOP"))) != 0){
    p <- p + ggplot2::geom_rect(data= data[data$type == "REGION" & data$description  %in% c("NBARC", "PLOOP"),],
                                mapping=ggplot2::aes(xmin=start,
                                                     xmax=end,
                                                     ymin=order-0.3,
                                                     ymax=order+0.3,
                                                     fill=description),
                                show.legend = TRUE) +
      points
  }

  # Add labels for regions, standard output = TRUE
  if(nrow(dplyr::filter(data, type == "REGION")) != 0 & label_regions == TRUE){
    p <- p + ggplot2::geom_text(data = data[data$type == "REGION", ],
                                ggplot2::aes(x = start + (end-start)/2,
                                             y = order,
                                             label = description),
                                colour = "white",
                                fontface = "bold",
                                size = 2)
  }

  # Add integrated domains
  if(nrow(dplyr::filter(data, type == "DOMAIN")) != 0){
    p <- p + ggplot2::geom_rect(data= data[data$type == "DOMAIN",],
                                mapping=ggplot2::aes(xmin=start,
                                                     xmax=end,
                                                     ymin=order-0.3,
                                                     ymax=order+0.3),
                                fill="#666666",
                                show.legend = FALSE)
  }

  # Add labels for integrated domains, standard output = TRUE
  if(nrow(dplyr::filter(data, type == "DOMAIN")) != 0 & label_domains == TRUE){
    p <- p + ggplot2::geom_label(data = data[data$type == "DOMAIN", ],
                                 ggplot2::aes(x = start + (end-start)/2,
                                              y = order,
                                              label = description),
                                 fill="#666666",
                                 colour = "white",
                                 fontface = "bold",
                                 size = 2)
  }

  # Add motifs, standard output = TRUE
  if(nrow(dplyr::filter(data, type == "MOTIF" & label_motifs == TRUE)) != 0){
    p <- p + ggplot2::geom_rect(data= data[data$type == "MOTIF",],
                                mapping=ggplot2::aes(xmin=start,
                                                     xmax=end,
                                                     ymin=order-0.25,
                                                     ymax=order-0.5),

                                fill="#c1c1c1")
  }

  # Add labels for motifs, standard output = TRUE
  if(nrow(dplyr::filter(data, type == "MOTIF")) != 0 & label_motifs == TRUE){
    p <- p + ggplot2::geom_text(data = data[data$type == "MOTIF", ],
                                ggplot2::aes(x = start + (end-start)/2,
                                             y = order,
                                             label = description),
                                size = 2) +
      ggplot2::labs(tag = "Motifs
      1 = P-loop
      2 = RNBS-A
      3 = Kin-2
      4 = RNBS-B
      5 = RNBS-C
      6 = GLPL
      7 = motif 12
      8 = RNBS-D
      9 = linker
      10 = MHD") + # label y-axis
      ggplot2::theme(plot.tag = ggplot2::element_text(size = 8),
                     plot.tag.position = "right")
  }

  # Add repeats, standard output = FALSE
  if(nrow(dplyr::filter(data, type == "REPEAT" & label_repeats == TRUE)) != 0){
    p <- p + ggplot2::geom_rect(data= data[data$type == "REPEAT",],
                                mapping=ggplot2::aes(xmin=start,
                                                     xmax=end,
                                                     ymin=order-0.25,
                                                     ymax=order-0.5),

                                fill="#c1c1c1",
                                show.legend = FALSE)
  }

  # Add labels for integrated domains, standard output = TRUE
  if(PDF != FALSE){
    ggplot2::ggsave(filename = PDF,p, width=max(data$end)/500, height=nrow(dplyr::distinct(data, seqname))/20, units="in", scale= 5, limitsize = FALSE)
  }
  return(p)
}

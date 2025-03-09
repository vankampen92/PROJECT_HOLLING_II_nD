mytheme <- theme_bw() + theme(
  legend.title  = element_text( size=17),
  #  legend.position = "bottom",
  #	legend.direction = "horizontal",
  legend.key = element_blank(),
  legend.text  = element_text( size=17),
  panel.background = element_blank(),
  panel.grid = element_blank(),
  text = element_text( family="Helvetica", size=19),
  strip.text.x = element_text(family = "Helvetica", size = 14),
  strip.text.y = element_text(family = "Helvetica", size = 14),
  panel.border = element_rect( colour = "black", size=1.3),
  axis.ticks = element_line(size = 1.3),
  strip.background = element_rect( fill = "transparent", size = 1.3, colour = "black"  ),
  strip.text = element_text(size = 19)
)

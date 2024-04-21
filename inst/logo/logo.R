# Logo creation
# From hexsticker.png to logo.png
# With thanks to: https://nanx.me/blog/post/rebranding-r-packages-with-hexagon-stickers/

hexSticker::sticker(
  subplot = "inst/logo/greenlogo.png",
  s_x = 1, s_y = 1, s_width = 0.8, s_height = 0.8,
  package = "", p_x = 1, p_y = 1, p_size = 8, h_size = 1.2, p_family = "pf",
  p_color = "#00857C", h_fill = "#FFFFFF", h_color = "#00857C",
  dpi = 640, filename = "man/figures/logo.png",
  white_around_sticker = TRUE
)



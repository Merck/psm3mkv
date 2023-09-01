# Logo creation
# From hexsticker.png to logo.png
# With thanks to: https://nanx.me/blog/post/rebranding-r-packages-with-hexagon-stickers/

hexSticker::sticker(subplot="inst/logo/hexsticker.png",
        s_x = 0.98, s_y = 0.95, s_width = 1.075, s_height = 1.075,
        package = "", p_x = 1, p_y = 1, p_size = 8, h_size = 1.2, p_family = "pf",
        p_color = "#7B1F58", h_fill = "#EFDBE5", h_color = "#7B1F58",
        dpi = 640, filename = "man/figures/logo.png",
        white_around_sticker = TRUE
        )

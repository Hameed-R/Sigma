#### Creating Geographic Maps with ggplot2

library(maps)
# Load the map data for the desired region
map_data <- map_data("world")
# Plot the map
ggplot() +
  geom_polygon(data = map_data, aes(x = long, y = lat, group = group),
               fill = "lightblue", color = "gray", size = 0.2) +
  coord_equal() +
  labs(title = "World Map") +
  theme_minimal()
####
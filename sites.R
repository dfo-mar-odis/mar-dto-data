require(arcpullr)
require(ggplot2)
require(ggspatial)
require(rnaturalearth)
require(rnaturalearthhires)
require(sf)
require(dplyr)

draft <- get_spatial_layer(
  "https://egisp.dfo-mpo.gc.ca/arcgis/rest/services/open_data_donnees_ouvertes/draft_conservation_network_sites/MapServer/1"
) |>
  dplyr::filter(Classification_E != "Existing site") |>
  select(Classification_E, SiteName_E)

oecm <- get_spatial_layer(
  "https://egisp.dfo-mpo.gc.ca/arcgis/rest/services/open_data_donnees_ouvertes/other_effective_area_based_conservation_measures/MapServer/0",
  where = "REGION_E='Maritimes'"
)


planning <- get_spatial_layer(
  "https://egisp.dfo-mpo.gc.ca/arcgis/rest/services/open_data_donnees_ouvertes/eastern_canada_marine_spatial_planning_areas/MapServer/0",
  where = "NAME_E='Scotian Shelf and Bay of Fundy'"
)

mpas <- get_spatial_layer(
  "https://egisp.dfo-mpo.gc.ca/arcgis/rest/services/open_data_donnees_ouvertes/oceans_act_marine_protected_areas_zones_de_protection_marines_de_la_loi_sur_les_oceans/MapServer/0"
) |>
  st_make_valid() |>
  dplyr::filter(NAME_E != "Laurentian Channel Marine Protected Area") |>
  st_filter(planning) |>
  group_by(NAME_E) |>
  summarize()

sites <- bind_rows(
  draft |>
    mutate(established = FALSE),
  oecm |>
    mutate(Classification_E = "Marine Refuge", established = TRUE),
  mpas |>
    mutate(Classification_E = "Marine Protected Area", established = TRUE)
) |>
  mutate(
    Classification = if_else(
      startsWith(Classification_E, "Tier"),
      "Draft Network Site",
      Classification_E
    ),

    Classification = factor(
      Classification,
      levels = c(
        "Marine Protected Area",
        "Marine Refuge",
        "Areas of Interest (AOI)",
        "Draft Network Site"
      )
    )
  )


cropbox <- planning |>
  st_bbox() |>
  st_as_sfc() |>
  st_buffer(200000) |>
  st_bbox() |>
  st_as_sfc()

plotbox <- planning |>
  st_buffer(100000) |>
  st_bbox()

canada <- ne_states(country = "canada", returnclass = "sf") |>
  st_crop(cropbox)
land <- ne_countries(returnclass = "sf", scale = "large") |>
  dplyr::filter(name != "Canada") |>
  st_make_valid() |>
  st_crop(cropbox)

targets::tar_load(data_glorys_bathy)
bathy_df <- data_glorys_bathy[planning] |>
  as.data.frame() |>
  mutate(deptho = as.numeric(deptho))


ggplot() +
  geom_sf(data = planning, fill = "transparent") +
  metR::geom_contour2(
    data = bathy_df,
    aes(x = longitude, y = latitude, z = deptho, label = after_stat(level)),
    breaks = c(100, 200, 1000, 2500, 3500),
    colour = "grey60",
    linewidth = 0.3
  ) +

  scale_fill_brewer(palette = "Blues", name = "Depth (m)") +
  geom_sf(data = canada) +
  geom_sf(data = land) +
  geom_sf(
    data = sites,
    aes(fill = Classification, linewidth = established, linetype = established)
  ) +
  geom_sf(
    data = sites[grepl("St. Anns Bank Marine Protected Area", sites$NAME_E), ],
    fill = "transparent",
    col = "red"
  ) +
  scale_fill_manual(
    values = c("#33a02c", "#b2df8a", "#1f78b4", "#a6cee3")
  ) +
  scale_linewidth_manual(values = c(0.3, 1), guide = "none") +
  scale_linetype_manual(values = c("dashed", "solid"), guide = "none") +
  annotation_scale(
    location = "br",
    width_hint = 0.2
  ) +
  annotation_north_arrow(
    location = "br",
    which_north = "true",
    pad_y = unit(0.5, "cm"),
    style = north_arrow_fancy_orienteering()
  ) +
  coord_sf(expand = FALSE, xlim = plotbox[c(1, 3)], ylim = plotbox[c(2, 4)]) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() +
  theme(
    legend.position = c(0.99, 0.99),
    legend.justification = c(1, 1),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 16)
  )

ggsave("MPAmap.png", units = "in", dpi = 300, width = 10, height = 10)

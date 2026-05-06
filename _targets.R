library(targets)
library(tarchetypes)
pkgs <- c(
  "marea",
  "sf",
  "arcpullr",
  "rnaturalearth",
  "mregions2",
  "reticulate",
  "stars",
  "ncmeta",
  "signal",
  "terra",
  "tidyterra",
  "patchwork",
  "RColorBrewer",
  "scales",
  "purrr",
  "ggplot2",
  "dplyr"
)
librarian::shelf(pkgs)

# Function to regularize a nearly regular stars grid
# Adapted from stars_to_terra by Remi Daigle, February 13, 2026, DTO project

regularize_stars <- function(s) {
  # Get dimension info
  dims <- st_dimensions(s)
  dim_names <- names(dims)

  # Find x, y, and z dimensions
  x_dim <- which(
    sapply(dims, function(d) {
      !is.na(d$refsys) && grepl("x", d$axis, ignore.case = TRUE)
    }) |
      grepl("^lon|^x$", dim_names, ignore.case = TRUE)
  )
  y_dim <- which(
    sapply(dims, function(d) {
      !is.na(d$refsys) && grepl("y", d$axis, ignore.case = TRUE)
    }) |
      grepl("^lat|^y$", dim_names, ignore.case = TRUE)
  )

  # Get coordinate values and compute regular spacing
  x_vals <- st_get_dimension_values(s, x_dim)
  y_vals <- st_get_dimension_values(s, y_dim)
  xres <- diff(range(x_vals)) / (length(x_vals) - 1)
  yres <- diff(range(y_vals)) / (length(y_vals) - 1)

  # Overwrite x/y dimension metadata with regular offset/delta
  dims[[x_dim]]$offset <- min(x_vals) - xres / 2
  dims[[x_dim]]$delta <- xres
  dims[[x_dim]]$values <- NULL # drop irregular point values if present

  dims[[y_dim]]$offset <- min(y_vals) - yres / 2
  dims[[y_dim]]$delta <- yres
  dims[[y_dim]]$values <- NULL

  st_dimensions(s) <- dims

  return(s)
}


raster_time_trend <- function(
  stars_obj,
  time_dim = "year",
  lat_dim = "latitude",
  lon_dim = "longitude"
) {
  time_numeric <- as.numeric(st_get_dimension_values(
    stars_obj,
    time_dim
  ))

  data_array <- stars_obj[[1]]

  result_array <- apply(data_array, c(1, 2), function(y) {
    valid <- !is.na(y)
    if (sum(valid) < 3) {
      return(c(NA_real_, NA_real_))
    }

    fit <- .lm.fit(cbind(1, time_numeric[valid]), y[valid])
    n <- sum(valid)
    residuals <- y[valid] - cbind(1, time_numeric[valid]) %*% fit$coefficients
    rss <- sum(residuals^2)

    se <- sqrt(
      (rss / (n - 2)) /
        sum((time_numeric[valid] - mean(time_numeric[valid]))^2)
    )
    t_stat <- fit$coefficients[2] / se
    pval <- 2 * pt(abs(t_stat), df = n - 2, lower.tail = FALSE)

    c(fit$coefficients[2], pval)
  })

  slope_array <- result_array[1, , ]
  pvalue_array <- result_array[2, , ]

  # convert back to stars
  result <- st_as_stars(
    list(slope = slope_array, pvalue = pvalue_array),
    dimensions = st_dimensions(stars_obj)[c(
      lon_dim,
      lat_dim
    )]
  )
}

# Set target options:
tar_option_set(
  packages = pkgs
)

list(
  tar_target(name = regions, command = {
    PA <- bind_rows(
      get_spatial_layer(
        "https://egisp.dfo-mpo.gc.ca/arcgis/rest/services/open_data_donnees_ouvertes/eastern_canada_marine_spatial_planning_areas/MapServer/0",
        where = "NAME_E='Newfoundland-Labrador Shelves'"
      ),
      get_spatial_layer(
        "https://egisp.dfo-mpo.gc.ca/arcgis/rest/services/open_data_donnees_ouvertes/eastern_canada_marine_spatial_planning_areas/MapServer/0",
        where = "NAME_E<>'Newfoundland-Labrador Shelves'"
      )
    )

    canada <- ne_states(country = "Canada", returnclass = "sf")

    A <- PA[grepl(PA$NAME_E, pattern = "Gulf"), ]
    B <- canada[canada$name_en %in% c("Quebec", "Newfoundland and Labrador"), ]
    C <- canada[
      canada$name_en %in%
        c("New Brunswick", "Nova Scotia", "Prince Edward Island"),
    ]

    # Create a fine resolution grid within A
    grid_points <- st_make_grid(A, cellsize = 0.025, what = "polygons") %>%
      st_as_sf() %>%
      st_filter(A) # Keep only points inside A

    # Determine which part is closer to which reference polygon
    regions <- grid_points |>
      st_as_sf() |>
      mutate(
        centroids = st_centroid(x),
        dist_to_B = as.numeric(st_distance(centroids, st_union(B))),
        dist_to_C = as.numeric(st_distance(centroids, st_union(C))),
        NAME_E = ifelse(dist_to_B < dist_to_C, "Quebec", "Gulf")
      ) |>
      group_by(NAME_E) |>
      summarise(geoms = st_union(x)) |>
      bind_rows(PA[!grepl(PA$NAME_E, pattern = "Gulf"), ]) |>
      mutate(
        NAME_E = if_else(
          NAME_E == "Scotian Shelf and Bay of Fundy",
          "Maritimes",
          NAME_E
        ),
        NAME_E = if_else(
          NAME_E == "Newfoundland-Labrador Shelves",
          "Newfoundland & Labrador",
          NAME_E
        )
      )
  }),

  tar_target(name = MPAs, command = {
    regions
    areas <- get_spatial_layer(
      "https://maps-cartes.ec.gc.ca/arcgis/rest/services/CWS_SCF/CPCAD/MapServer/0",
      where = "BIOME='M' AND MGMT_E='Fisheries And Oceans Canada'"
    ) |>
      group_by(NAME_E, NAME_F) |>
      summarise(geoms = st_make_valid(st_union(geoms)))

    centroids <- st_centroid(areas$geoms) |>
      st_as_sf() |>
      st_join(regions)

    areas$region <- centroids$NAME_E

    areas <- areas |>
      filter(!is.na(region))

    sf::sf_use_s2(FALSE)

    eez_buffer <- data_eez |>
      st_transform(3979) |>
      st_union() |>
      st_simplify(dTolerance = 50000) |>
      st_buffer(550000) |>
      st_make_valid() |>
      st_transform(st_crs(areas))

    # Get union of all MPAs
    all_mpa_union <- st_union(areas$geoms)

    # Compute difference (the "outside" area)
    outside_geom <- st_difference(eez_buffer, all_mpa_union) |> st_make_valid()

    sf::sf_use_s2(TRUE)

    # Step 6: Create the Outside row
    outside_row <- tibble(
      NAME_E = "Non_Conservation_Area",
      NAME_F = "Extérieur",
      geoms = outside_geom,
      region = NA
    )

    # Step 7: Bind with existing areas
    areas_full <- bind_rows(areas, outside_row)

    areas_full$date_of_establishment <- c(
      2019,
      2005,
      1981,
      2017,
      2017,
      2016,
      2007,
      2017,
      2022,
      2017,
      2017,
      2005,
      2005,
      2016,
      2017,
      2004,
      2017,
      2017,
      2017,
      2016,
      2019,
      2016,
      1998,
      2004,
      1996,
      2004,
      2004,
      1998,
      2002,
      2014,
      1986,
      2006,
      2017,
      2016,
      2017,
      2017,
      1985,
      2017,
      2017,
      2017,
      2017,
      2017,
      2017,
      2017,
      2000
    )
    areas_full
  }),

  tar_target(name = data_eez, command = {
    data_eez <- mregions2::gaz_geometry(8493)
  }),

  tar_target(name = data_glorys_bathy, command = {
    bbox <- st_bbox(regions[regions$NAME_E == "Maritimes", ])

    tmp <- tempfile(fileext = ".nc")

    file <- get_CMEMS_ncdf(
      username = Sys.getenv("CMEMS_USERNAME"),
      password = Sys.getenv("CMEMS_PASSWORD"),
      dataset_id = "cmems_mod_glo_phy_my_0.083deg_static",
      variables = list("deptho_lev", "deptho"),
      minimum_longitude = bbox[1],
      maximum_longitude = bbox[3],
      minimum_latitude = bbox[2],
      maximum_latitude = bbox[4],
      output_filename = tmp
    )

    read_ncdf(file, proxy = FALSE) |>
      regularize_stars()
  }),

  tar_target(name = data_glorys_sst, command = {
    bbox <- st_bbox(regions[regions$NAME_E == "Maritimes", ])

    tmp <- tempfile(fileext = ".nc")

    file <- get_CMEMS_ncdf(
      username = Sys.getenv("CMEMS_USERNAME"),
      password = Sys.getenv("CMEMS_PASSWORD"),
      dataset_id = "cmems_mod_glo_phy_my_0.083deg_P1D-m",
      variables = list("thetao"),
      minimum_longitude = bbox[1],
      maximum_longitude = bbox[3],
      minimum_latitude = bbox[2],
      maximum_latitude = bbox[4],
      output_filename = tmp,
      start_datetime = "1993-01-01T00:00:00",
      end_datetime = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    )

    read_ncdf(
      file,
      proxy = FALSE,
      ncsub = cbind(start = c(1, 1, 1, 1), count = c(NA, NA, 1, NA))
    ) |>
      regularize_stars() |>
      adrop()
  }),

  tar_target(name = data_glorys_sst_butter90, command = {
    butter_filter <- butter(n = 3, W = 2 / 90, type = "low")

    sst_filtered <- st_apply(
      data_glorys_sst,
      MARGIN = c(1, 2),
      FUN = function(x) {
        if (all(is.na(x) | is.nan(x))) {
          return(rep(NA_real_, length(x)))
        }
        filtfilt(butter_filter, x)
      },
      .fname = "time"
    ) |>
      aperm(c("longitude", "latitude", "time"))

    st_dimensions(sst_filtered)["time"] <- st_dimensions(data_glorys_sst)[
      "time"
    ]

    # Rename thetao
    names(sst_filtered) <- "sst_butter90"
    sst_filtered
  }),

  tar_target(name = data_glorys_bottomT, command = {
    bbox <- st_bbox(regions[regions$NAME_E == "Maritimes", ])

    tmp <- tempfile(fileext = ".nc")

    file <- get_CMEMS_ncdf(
      username = Sys.getenv("CMEMS_USERNAME"),
      password = Sys.getenv("CMEMS_PASSWORD"),
      dataset_id = "cmems_mod_glo_phy_my_0.083deg_P1D-m",
      variables = list("bottomT"),
      minimum_longitude = bbox[1],
      maximum_longitude = bbox[3],
      minimum_latitude = bbox[2],
      maximum_latitude = bbox[4],
      output_filename = tmp,
      start_datetime = "1993-01-01T00:00:00",
      end_datetime = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    )

    read_ncdf(
      file,
      proxy = FALSE
    ) |>
      regularize_stars()
  }),

  tar_target(name = data_onset_of_spring, command = {
    time_vals <- st_get_dimension_values(data_glorys_sst_butter90, "time")
    years <- as.integer(format(time_vals, "%Y"))

    minyear <- as.integer(min(format(time_vals, "%Y")))

    apriltemp <- data_glorys_sst_butter90[,,,] |>
      slice(
        "time",
        which(
          format(time_vals, "%m") == "04" &
            as.integer(format(time_vals, "%Y")) %in%
              minyear:(minyear +
                30)
        )
      ) |>
      st_apply(c("longitude", "latitude"), mean, na.rm = TRUE)

    result_list <- lapply(sort(unique(years)), function(yr) {
      idx <- which(years == yr)

      sst_yr <- st_apply(
        data_glorys_sst_butter90[,,, idx],
        c("longitude", "latitude"),
        function(x) {
          if (all(is.na(x))) {
            return(x)
          }
          min_idx <- which.min(x)
          max_idx <- which.max(x)
          x[seq_len(min_idx)] <- NA
          x[max_idx:length(x)] <- NA
          x
        },
        .fname = "time"
      ) |>
        aperm(c("longitude", "latitude", "time"))

      st_dimensions(sst_yr)[
        "time"
      ] <- st_dimensions(data_glorys_sst_butter90[,,, idx])[
        "time"
      ]

      spring <- (sst_yr > apriltemp) |>
        st_apply(c("longitude", "latitude"), function(x) {
          if (all(is.na(x))) {
            return(NA_real_)
          }
          min(which(x), na.rm = TRUE)
        })
      names(spring) <- "Onset of Spring (DOY)"
      spring
    })

    result_stars <- do.call(
      c,
      c(result_list, list(along = list(year = sort(unique(years)))))
    )

    result_stars$`Onset of Spring (DOY)`[
      result_stars$`Onset of Spring (DOY)` > 365
    ] <- NA

    anomaly <- result_stars |>
      st_apply(
        c("longitude", "latitude"),
        function(x) {
          if (all(is.na(x))) {
            rep(NA_real_, length(x))
          }
          anomalies <- mean(x[1:30], na.rm = TRUE) - x
          anomalies / 7
        },
        .fname = "year"
      ) |>
      aperm(c("longitude", "latitude", "year"))
    names(anomaly) <- "Onset of Spring Anomaly (Weeks)"
    st_dimensions(anomaly)["year"] <- st_dimensions(result_stars)["year"]

    c(anomaly, result_stars) |>
      slice(
        "year",
        which(sort(unique(years)) < as.integer(format(Sys.time(), "%Y")))
      )
  }),

  tar_target(name = data_onset_of_fall, command = {
    time_vals <- st_get_dimension_values(data_glorys_sst_butter90, "time")
    years <- as.integer(format(time_vals, "%Y"))

    minyear <- as.integer(min(format(time_vals, "%Y")))

    octtemp <- data_glorys_sst_butter90[,,,] |>
      slice(
        "time",
        which(
          format(time_vals, "%m") == "10" &
            as.integer(format(time_vals, "%Y")) %in%
              minyear:(minyear +
                30)
        )
      ) |>
      st_apply(c("longitude", "latitude"), mean, na.rm = TRUE)

    result_list <- lapply(sort(unique(years)), function(yr) {
      idx <- which(years == yr)

      sst_yr <- st_apply(
        data_glorys_sst_butter90[,,, idx],
        c("longitude", "latitude"),
        function(x) {
          if (all(is.na(x))) {
            return(x)
          }
          max_idx <- which.max(x)
          x[seq_len(max_idx)] <- NA # NA any temperature before max temp
          x
        },
        .fname = "time"
      ) |>
        aperm(c("longitude", "latitude", "time"))

      st_dimensions(sst_yr)[
        "time"
      ] <- st_dimensions(data_glorys_sst_butter90[,,, idx])[
        "time"
      ]

      fall <- (sst_yr < octtemp) |>
        st_apply(c("longitude", "latitude"), function(x) {
          if (all(is.na(x))) {
            return(NA_real_)
          }
          min(which(x), na.rm = TRUE)
        })
      names(fall) <- "Onset of Fall (DOY)"
      fall
    })

    result_stars <- do.call(
      c,
      c(result_list, list(along = list(year = sort(unique(years)))))
    )

    result_stars$`Onset of Fall (DOY)`[
      result_stars$`Onset of Fall (DOY)` > 365
    ] <- NA

    anomaly <- result_stars |>
      st_apply(
        c("longitude", "latitude"),
        function(x) {
          if (all(is.na(x))) {
            rep(NA_real_, length(x))
          }
          anomalies <- mean(x[1:30], na.rm = TRUE) - x
          anomalies / 7
        },
        .fname = "year"
      ) |>
      aperm(c("longitude", "latitude", "year"))
    names(anomaly) <- "Onset of Fall Anomaly (Weeks)"
    st_dimensions(anomaly)["year"] <- st_dimensions(result_stars)["year"]

    c(anomaly, result_stars) |>
      slice(
        "year",
        which(sort(unique(years)) < as.integer(format(Sys.time(), "%Y")))
      )
  }),

  tar_target(name = data_length_growing_season, command = {
    result_stars <- data_onset_of_fall["Onset of Fall (DOY)"] -
      data_onset_of_spring["Onset of Spring (DOY)"]
    names(result_stars) <- "Length of the Growing Season (days)"

    anomaly <- result_stars |>
      st_apply(
        c("longitude", "latitude"),
        function(x) {
          if (all(is.na(x))) {
            rep(NA_real_, length(x))
          }
          anomalies <- x - mean(x[1:30], na.rm = TRUE)
          anomalies / 7
        },
        .fname = "year"
      ) |>
      aperm(c("longitude", "latitude", "year"))
    names(anomaly) <- "Length of the Growing Season Anomaly (Weeks)"
    st_dimensions(anomaly)["year"] <- st_dimensions(result_stars)["year"]

    c(anomaly, result_stars)
  }),

  tar_target(name = ind_glorys_bottomT_trend, command = {
    annualmean <- aggregate(
      data_glorys_bottomT,
      by = "year",
      FUN = mean,
      na.rm = TRUE
    )

    annualmean <- annualmean |>
      st_set_dimensions(
        "time",
        names = "year",
        st_get_dimension_values(annualmean, "time") |>
          format("%Y") |>
          as.integer()
      ) |>
      aperm(c("longitude", "latitude", "year"))

    raster_time_trend(annualmean)
  }),

  tar_target(name = ind_glorys_sst_trend, command = {
    annualmean <- aggregate(
      data_glorys_sst,
      by = "year",
      FUN = mean,
      na.rm = TRUE
    )

    annualmean <- annualmean |>
      st_set_dimensions(
        "time",
        names = "year",
        st_get_dimension_values(annualmean, "time") |>
          format("%Y") |>
          as.integer()
      ) |>
      aperm(c("longitude", "latitude", "year"))

    raster_time_trend(annualmean)
  }),

  tar_target(name = ind_onset_of_spring_trend, command = {
    raster_time_trend(data_onset_of_spring["Onset of Spring Anomaly (Weeks)"])
  }),

  tar_target(name = ind_onset_of_fall_trend, command = {
    raster_time_trend(data_onset_of_fall["Onset of Fall Anomaly (Weeks)"])
  }),

  tar_target(name = ind_length_growing_season_trend, command = {
    raster_time_trend(data_length_growing_season[
      "Length of the Growing Season Anomaly (Weeks)"
    ])
  }),

  tar_target(name = plot_spring_fall_example, command = {
    example_years <- c(1993, 2009, 2025)

    bbox <- st_bbox(MPAs[
      MPAs$NAME_E == "Western and Emerald Banks Marine Refuge",
    ])
    ex_lon <- mean(c(bbox["xmin"], bbox["xmax"]))
    ex_lat <- mean(c(bbox["ymin"], bbox["ymax"]))

    pt <- st_point(c(ex_lon, ex_lat)) |> st_sfc(crs = 4326) |> st_as_sf()

    combined <- left_join(
      # normal sst
      st_extract(data_glorys_sst, pt) |>
        as.data.frame() |>
        mutate(sst = as.numeric(thetao)) |>
        select(-thetao),
      # filtered sst
      st_extract(data_glorys_sst_butter90, pt) |>
        as.data.frame() |>
        select(-x) |>
        rename(sst_filt = sst_butter90),
      by = "time"
    ) |>
      mutate(
        year = as.numeric(format(time, "%Y")),
        month = as.numeric(format(time, "%m"))
      )

    indicators <- st_extract(data_onset_of_spring, pt) |>
      as.data.frame() |>
      left_join(
        st_extract(data_onset_of_fall, pt) |>
          as.data.frame() |>
          select(-x),
        by = "year"
      )

    april_val <- combined |>
      filter(year <= (min(year) + 30), month == 4) |>
      summarise(mean(sst_filt, na.rm = TRUE)) |>
      as.numeric()

    oct_val <- combined |>
      filter(year <= (min(year) + 30), month == 10) |>
      summarise(mean(sst_filt, na.rm = TRUE)) |>
      as.numeric()

    # Build plots without legends
    plots <- lapply(example_years, function(yr) {
      df_yr <- combined |>
        dplyr::filter(year == yr) |>
        mutate(date = as.Date(time))

      spring_doy <- indicators$`Onset of Spring (DOY)`[indicators$year == yr]
      fall_doy <- indicators$`Onset of Fall (DOY)`[indicators$year == yr]
      season_ribbon <- df_yr[spring_doy:fall_doy, ]

      p <- ggplot(df_yr, aes(x = date)) +
        geom_line(
          aes(y = sst, colour = "Raw SST"),
          linewidth = 0.4,
          alpha = 0.6
        ) +
        geom_line(aes(y = sst_filt, colour = "Filtered SST"), linewidth = 0.9) +
        geom_hline(
          aes(yintercept = april_val, linetype = "April mean (clim.)"),
          colour = "#1b7837",
          linewidth = 0.7
        ) +
        geom_hline(
          aes(yintercept = oct_val, linetype = "October mean (clim.)"),
          colour = "#762a83",
          linewidth = 0.7
        ) +
        geom_ribbon(
          data = season_ribbon,
          aes(ymin = -Inf, ymax = Inf),
          fill = "#ffd700",
          alpha = 0.15
        ) +
        geom_vline(
          xintercept = df_yr$date[spring_doy],
          colour = "#1b7837",
          linetype = "dashed"
        ) +
        annotate(
          "text",
          x = df_yr$date[spring_doy],
          y = max(df_yr$sst, na.rm = TRUE),
          label = paste0("Spring\nDOY ", spring_doy),
          hjust = -0.1,
          size = 3,
          colour = "#1b7837"
        ) +
        geom_vline(
          xintercept = df_yr$date[fall_doy],
          colour = "#762a83",
          linetype = "dashed"
        ) +
        annotate(
          "text",
          x = df_yr$date[fall_doy],
          y = max(df_yr$sst, na.rm = TRUE),
          label = paste0("Fall\nDOY ", fall_doy),
          hjust = 1.1,
          size = 3,
          colour = "#762a83"
        ) +
        scale_colour_manual(
          name = NULL,
          values = c("Raw SST" = "grey50", "Filtered SST" = "#2166ac")
        ) +
        scale_linetype_manual(
          name = NULL,
          values = c(
            "April mean (clim.)" = "solid",
            "October mean (clim.)" = "solid"
          )
        ) +
        scale_x_date(expand = FALSE) +
        labs(
          title = paste("Year:", yr),
          subtitle = paste0(
            "Growing season: DOY ",
            spring_doy,
            " – ",
            fall_doy,
            " (",
            fall_doy - spring_doy,
            " days)"
          ),
          x = NULL,
          y = "SST (°C)"
        ) +
        theme_bw(base_size = 11) +
        theme(
          legend.position = "bottom"
        )
    })

    # Collect legends from one plot
    guide_plot <- plots[[1]] + theme(legend.position = "bottom")

    fig <- wrap_plots(plots, ncol = 1) +
      guide_area() +
      plot_layout(
        ncol = 1,
        guides = "collect", # merges all legends into one
        heights = c(1, 1, 1, 0.1) # 3 equal panels + small legend row
      ) +
      plot_annotation(
        title = "Onset of Spring / Fall methodology",
        subtitle = paste0(
          "Example pixel — lon: ",
          round(ex_lon, 2),
          ", lat: ",
          round(ex_lat, 2)
        ),
        theme = theme(
          plot.title = element_text(face = "bold"),
          legend.position = "bottom"
        )
      )
    fig
  }),

  tar_target(name = plot_spring_minimaps, command = {
    test <- st_set_dimensions(
      data_onset_of_spring["Onset of Spring Anomaly (Weeks)"],
      which = "year",
      values = as.character(st_get_dimension_values(
        data_onset_of_spring,
        "year"
      )),
      names = "band"
    ) |>
      rast()

    test |>
      plot()
  }),

  tar_target(name = plot_growing_season_minimaps, command = {
    test <- st_set_dimensions(
      data_length_growing_season[
        "Length of the Growing Season Anomaly (Weeks)"
      ],
      which = "year",
      values = as.character(st_get_dimension_values(
        data_length_growing_season,
        "year"
      )),
      names = "band"
    ) |>
      rast()

    test |>
      plot()
  }),

  tar_target(name = plot_spring_trend_map, command = {
    # Land / coastline for background
    land <- ne_countries(scale = "medium", returnclass = "sf") |>
      st_crop(
        xmin = -72,
        xmax = -52,
        ymin = 40,
        ymax = 52
      )

    r <- rast(ind_onset_of_spring_trend)

    slope <- r[["slope_slope"]]
    pvalue <- r[["pvalue_pvalue"]]

    # Mask slope where pvalue > 0.05 (set to NA)
    slope_masked <- mask(slope, pvalue > 0.05, maskvalues = 1)

    # Define color palette
    cols <- colorRampPalette(c(
      "#2166AC", # deep blue
      "#74ADD1", # mid blue
      "lightblue", # light blue
      "#FFFFFF", # white (zero)
      "#FEC981", # light orange
      "#F4820E", # orange
      "#C13B06" # deep red-orange
    ))(100)

    # Symmetric limits around zero
    lim <- max(abs(minmax(slope_masked)), na.rm = TRUE)

    plot(
      slope_masked,
      col = cols,
      range = c(-lim, lim),
      main = "Onset of Spring Trend (Slope)",
      plg = list(title = "slope\n(units/yr)"),
      colNA = "grey"
    )
    plot(st_geometry(land), add = TRUE, border = "black", col = "grey80")
  }),

  tar_target(name = plot_spring_anomaly_timeseries, command = {
    selected_mpas <- MPAs |>
      filter(
        NAME_E %in%
          c(
            "St. Anns Bank Marine Protected Area",
            "Western and Emerald Banks Marine Refuge",
            "Gully Marine Protected Area",
            "Northeast Channel Marine Refuge"
          )
      )
    #TODO: plot anomalies
  })
)

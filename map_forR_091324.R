# Load libraries
library(sf)
library(ggplot2)
library(rnaturalearth)

##################################JUST MAPS NO POINTS#########################33

# Set the path to the directory where the shapefile is located
# Adjust the path based on where you downloaded and extracted the shapefile
shapefile_path <- "C:/Users/Lindsay Alma/Downloads/cb_2018_us_state_500k(1)/cb_2018_us_state_500k.shp"

# Read the shapefile
us_states <- st_read("C:/Users/Lindsay Alma/Downloads/cb_2018_us_state_500k(1)/cb_2018_us_state_500k.shp")

# View the structure of the shapefile to find Washington
print(us_states)

# Filter to get only Washington State
washington<- us_states[us_states$NAME == "Washington", ]
ggplot(data = washington) +
  geom_sf(fill = "darkgray", color = "black") +
  theme_minimal()

california <- us_states[us_states$NAME == "California", ]
ggplot(data = washington) +
  geom_sf(fill = "darkgray", color = "black") +
  theme_minimal()

oregon <- us_states[us_states$NAME == "Oregon", ]
ggplot(data = oregon) +
  geom_sf(fill = "darkgray", color = "black") +
  theme_minimal()


alaska <- us_states[us_states$NAME == "Alaska", ]
ggplot(data = Alaska) +
  geom_sf(fill = "darkgray", color = "black") +
  theme_minimal()

#alaska fix lat and long
ggplot(data = Alaska) +
  geom_sf(fill = "darkgray", color = "black") +
  coord_sf(xlim = c(-179, -129),  # Approximate longitude limits for Alaska
           ylim = c(51, 72),      # Approximate latitude limits for Alaska
           expand = FALSE) +
    theme_minimal()



# ########Combine Washington, Oregon, California, alaska, and British Columbia into one sf object##############
################################################################################################

common_crs <- 4326

washington <- st_transform(washington, crs = common_crs)
oregon <- st_transform(oregon, crs = common_crs)
california <- st_transform(california, crs = common_crs)
british_columbia <- st_transform(british_columbia, crs = common_crs)
alaska <- st_transform(alaska, crs = common_crs)

colnames(washington)
colnames(oregon)
colnames(california)
colnames(british_columbia)
colnames(alaska)

names(washington)
names(oregon)
names(california)
names(british_columbia)
names(alaska)
common_columns <- Reduce(intersect, list(
  names(washington),
  names(oregon),
  names(california),
  names(alaska),
  names(british_columbia)
))

washington <- washington[, common_columns]
oregon <- oregon[, common_columns]
california <- california[, common_columns]
british_columbia <- british_columbia[, common_columns]
alaska <- alaska[, common_columns]

# Combine them after ensuring the CRS is the same
pacific_northwest <- rbind(washington, oregon, california, british_columbia, alaska)


# Define an orthographic CRS (for example, centered on the region)
orthographic_crs <- "+proj=ortho +lat_0=45 +lon_0=-120 +datum=WGS84 +units=m +no_defs"

# Transform to the orthographic projection
pacific_northwest_orthographic <- st_transform(pacific_northwest, crs = orthographic_crs)

# Calculate bounding box for Orthographic projection
bbox_orthographic <- st_bbox(pacific_northwest_orthographic)


# Create a data frame for the points
points_df <- data.frame(
  lat = c(55.4967115, 55.73943425, 55.70733725, 55.24908025, 55.597631, 55.2832995,
          38.318556, 38.333387, 38.20455198, 38.30995233, 38.10502768, 38.17850748,
          51.64505167, 51.7706365, 51.80798667, 51.92398, 51.89405167,
          43.34607183, 43.35929283, 44.61293517, 44.61944483, 44.62305733,
          32.62567667, 32.77164333, 32.97057017, 32.7141475, 32.78958667,
          48.4685835, 48.48291883, 48.69187517, 48.7046845, 48.5634),
  lon = c(-133.1694408, -133.3120588, -133.3422598, -132.8815268, -133.1627413, -133.3335205,
          -123.054278, -123.059461, -122.9268415, -123.0590258, -122.8463775, -122.9088285,
          -128.11946, -128.0599033, -128.2369017, -128.46993, -128.2348717,
          -124.3181705, -124.3098838, -124.0139907, -124.0325398, -124.0233798,
          -117.106445, -117.2196167, -117.2518788, -117.225227, -117.2311152,
          -123.0027512, -123.0741698, -122.9524467, -123.054054, -122.9355568)
)



# Convert to sf object, assuming the points are in WGS 84 (EPSG 4326)
points_sf <- st_as_sf(points_df, coords = c("lon", "lat"), crs = 4326)

# Transform the points to the orthographic projection (same as pacific_northwest_orthographic)
points_sf_ortho <- st_transform(points_sf, crs = orthographic_crs)

# Plot using Orthographic Projection
# Plot the map with the points
ggplot(data = pacific_northwest_orthographic) +
  geom_sf(fill = "darkgray", color = "black") +          # Plot the regions
  geom_sf(data = points_sf_ortho, color = "red", size = 3) +  # Add the points in red
  theme_minimal() +
  coord_sf(
    xlim = c(bbox_orthographic["xmin"], bbox_orthographic["xmax"]),
    ylim = c(bbox_orthographic["ymin"], bbox_orthographic["ymax"]),
    expand = FALSE
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "Pacific Northwest with Points in Orthographic Projection"
  )


#####################MAPS WITH POINTS#########################################

#Alaska
points_data <- data.frame(
    lon = c(-133.169655,	-133.169122,	-133.169696,	-133.16929,	-133.312096,	-133.312085,	-133.312032,	-133.312022,	-133.34257,	-133.341967,	-133.342562,	-133.34194,	-132.881748,	-132.881276,	-132.881797,	-132.881286,	-133.163002,	-133.162523,	-133.162951,	-133.162489,	-133.333546,	-133.333433,	-133.333594,	-133.333509),
  lat = c(55.496851,	55.496659,	55.496728,	55.496608,	55.739648,	55.739262,	55.739272,	55.739555,	55.707346,	55.70741,	55.70726,	55.707333,	55.249207,	55.248974,	55.249188,	55.248952,	55.59751,	55.597762,	55.597513,	55.597739,	55.28312,	55.283468,	55.283128,	55.283482)
)
points_sf <- st_as_sf(points_data, coords = c("lon", "lat"), crs = 4326)
Alaska <- us_states[us_states$NAME == "Alaska", ]
alaska_transformed <- st_transform(Alaska, crs = 3338)
points_sf_transformed <- st_transform(points_sf, crs = st_crs(alaska_transformed))
ggplot(data = Alaska) +
  geom_sf(fill = "darkgray", color = "black") + geom_sf(data = points_sf_transformed, color = "red", size = 1) +
  coord_sf(xlim = c(-179, -129),  # Approximate longitude limits for Alaska
           ylim = c(51, 72),      # Approximate latitude limits for Alaska
           expand = FALSE) +
  theme_minimal()


#WASHINGTON
points_data <- data.frame(
  lon = c(-123.002299,	-123.002602,	-123.002876,	-123.002664,	-123.002898,	-123.003168,	-123.074358,	-123.074129,	-123.073933,	-123.074413,	-123.074281,	-123.073905,	-122.952169,	-122.952466,	-122.952737,	-122.952152,	-122.952432,	-122.952724,	-123.05377,	-123.054053,	-123.054377,	-123.053744,	-123.054033,	-123.054347,	-122.935871,	-122.935541,	-122.93523,	-122.93589,	-122.935569,	-122.93524
),
  lat = c(48.468591,	48.468676,	48.468817,	48.468323,	48.46846,	48.468634,	48.482598,	48.482775,	48.482822,	48.482656,	48.482796,	48.483866,	48.692011,	48.69188,	48.691775,	48.691983,	48.691851,	48.691751,	48.704823,	48.70471,	48.704541,	48.704809,	48.704672,	48.704552,	48.563311,	48.563396,	48.563442,	48.563353,	48.563431,	48.563467
)
)
points_sf <- st_as_sf(points_data, coords = c("lon", "lat"), crs = 4326)
points_sf_transformed <- st_transform(points_sf, crs = st_crs(alaska_transformed))
# Filter to get only Washington State
washington <- us_states[us_states$NAME == "Washington", ]
# Plot Washington State
ggplot(data = washington) +
  geom_sf(fill = "darkgray", color = "black") + geom_sf(data = points_sf_transformed, color = "red", size = 1) +
  theme_minimal()

####WASHINGTON INSET
points_data <- data.frame(
  lon = c(-123.002299,	-123.002602,	-123.002876,	-123.002664,	-123.002898,	-123.003168,	-123.074358,	-123.074129,	-123.073933,	-123.074413,	-123.074281,	-123.073905,	-122.952169,	-122.952466,	-122.952737,	-122.952152,	-122.952432,	-122.952724,	-123.05377,	-123.054053,	-123.054377,	-123.053744,	-123.054033,	-123.054347,	-122.935871,	-122.935541,	-122.93523,	-122.93589,	-122.935569,	-122.93524
  ),
  lat = c(48.468591,	48.468676,	48.468817,	48.468323,	48.46846,	48.468634,	48.482598,	48.482775,	48.482822,	48.482656,	48.482796,	48.483866,	48.692011,	48.69188,	48.691775,	48.691983,	48.691851,	48.691751,	48.704823,	48.70471,	48.704541,	48.704809,	48.704672,	48.704552,	48.563311,	48.563396,	48.563442,	48.563353,	48.563431,	48.563467
  )
)
points_sf <- st_as_sf(points_data, coords = c("lon", "lat"), crs = 4326)
Alaska <- us_states[us_states$NAME == "Washington", ]
alaska_transformed <- st_transform(Alaska, crs = 3338)
points_sf_transformed <- st_transform(points_sf, crs = st_crs(alaska_transformed))
ggplot(data = Alaska) +
  geom_sf(fill = "darkgray", color = "black") + geom_sf(data = points_sf_transformed, color = "red", size = 3) +
  coord_sf(xlim = c(-122.75, -123.3),  # Approximate longitude limits for Alaska
           ylim = c(48.4, 48.73),      # Approximate latitude limits for Alaska
           expand = FALSE) +
  theme_minimal()

#OREGON
points_data <- data.frame(
  lon = c(-124.3181705,	-124.3098838,	-124.0139907,	-124.0325398,	-124.0233798
),
  lat = c(43.34607183,	43.35929283,	44.61293517,	44.61944483,	44.62305733
  )
)
points_sf <- st_as_sf(points_data, coords = c("lon", "lat"), crs = 4326)
points_sf_transformed <- st_transform(points_sf, crs = st_crs(alaska_transformed))
oregon <- us_states[us_states$NAME == "Oregon", ]
# Plot Washington State
ggplot(data = oregon) +
  geom_sf(fill = "darkgray", color = "black") + geom_sf(data = points_sf_transformed, color = "red", size = 2) +
  theme_minimal()

#California
points_data <- data.frame(
  lon = c(-123.053792,	-123.053538,	-123.053656,	-123.053948,	-123.053749,	-123.053713,	-123.059385,	-123.059405,	-123.05947,	-123.059494,	-123.059449,	-123.059563,	-122.926448,	-122.926358,	-122.926183,	-122.927469,	-122.927547,	-122.927044,	-123.059161,	-123.058908,	-123.058737,	-123.059039,	-123.05901,	-123.0593,	-122.846577,	-122.846131,	-122.846491,	-122.846248,	-122.846629,	-122.846189,	-122.909126,	-122.909032,	-122.908804,	-122.908865,	-122.908716,	-122.908428,	-117.10671,	-117.10646,	-117.10626,	-117.10664,	-117.10639,	-117.10621,	-117.21993,	-117.21967,	-117.21935,	-117.21988,	-117.21958,	-117.21929,	-117.25244,	-117.251969,	-117.251359,	-117.252467,	-117.251914,	-117.251124,	-117.22563,	-117.225248,	-117.224798,	-117.225623,	-117.225261,	-117.224802,	-117.23137,	-117.23112,	-117.2308,	-117.231421,	-117.23113,	-117.23085
),
  lat = c(38.31867,	38.318892,	38.319145,	38.31871,	38.318894,	38.319131,	38.33316104,	38.33338198,	38.333601,	38.33317101,	38.33336798,	38.33363897,	38.20456596,	38.20434996,	38.20434602,	38.20476796,	38.20465498,	38.20462698,	38.31034399,	38.30987602,	38.31004601,	38.30976698,	38.30978701,	38.30989396,	38.10507302,	38.10504402,	38.10517101,	38.10490103,	38.10499599,	38.10498099,	38.17865096,	38.17843596,	38.17826699,	38.17882799,	38.17861299,	38.17824997,	32.62597,	32.62565,	32.62539,	32.62603,	32.62567,	32.62535,	32.77188,	32.7716,	32.77132,	32.77195,	32.7717,	32.77141,	32.970515,	32.970618,	32.970645,	32.970573,	32.97068,	32.97039,	32.71385,	32.714114,	32.714321,	32.713975,	32.714207,	32.714418,	32.78932,	32.78956,	32.78982,	32.78935,	32.78961,	32.78986
)
)
points_sf <- st_as_sf(points_data, coords = c("lon", "lat"), crs = 4326)
points_sf_transformed <- st_transform(points_sf, crs = st_crs(alaska_transformed))
California<- us_states[us_states$NAME == "California", ]
# Plot Washington State
ggplot(data = California) +
  geom_sf(fill = "darkgray", color = "black") + geom_sf(data = points_sf_transformed, color = "red", size = 1) +
  theme_minimal()

#############CANADA###########################################

shapefile_path <- "C:/Users/Lindsay Alma/Downloads/lpr_000b21a_e/lpr_000b21a_e.shp"

# Read the shapefile
canada_provinces <- st_read(shapefile_path)

# View the structure of the shapefile to find British Columbia
print(canada_provinces)
head(canada_provinces)
str(canada_provinces)


points_data <- data.frame(
  lon = c(-128.119,	-128.11946,	-128.12004,	-128.11891,	-128.11926,	-128.12009,	-128.05959,	-128.05971,	-128.06019,	-128.05965,	-128.06001,	-128.06027,	-128.23622,	-128.23663,	-128.237,	-128.23691,	-128.23721,	-128.23744,	-128.46947,	-128.46956,	-128.46945,	-128.47015,	-128.47061,	-128.47034,	-128.23553,	-128.23503,	-128.23462,	-128.23503,	-128.23474,	-128.23428
),
  lat = c(  51.6454,	51.64555,	51.64532,	51.64495,	51.64474,	51.64435,	51.77027,	51.77032,	51.771249,	51.77021,	51.77053,	51.77124,	51.80806,	51.80814,	51.80824,	51.80761,	51.80782,	51.80805,	51.92354,	51.92381,	51.92413,	51.9239,	51.92413,	51.92437,	51.89411,	51.89399,	51.89383,	51.89445,	51.89413,	51.8938)

)
points_sf <- st_as_sf(points_data, coords = c("lon", "lat"), crs = 4326)
points_sf_transformed <- st_transform(points_sf, crs = st_crs(alaska_transformed))


# Filter to get only British Columbia
british_columbia <- canada_provinces[canada_provinces$PRENAME == "British Columbia", ]

colnames(canada_provinces)
british_columbia_transformed_albers <- st_transform(british_columbia, crs = 3005)
british_columbia_transformed_utm <- st_transform(british_columbia, crs = 32610)

british_columbia <- british_columbia %>%
  dplyr::select(PRNAME, LANDAREA, geometry) %>%
  dplyr::rename(NAME = PRNAME, ALAND = LANDAREA)


# Plot British Columbia
ggplot(data = british_columbia_transformed_utm) +
  geom_sf(fill = "darkgray", color = "black") + geom_sf(data = points_sf_transformed, color = "red", size = 1) +
    theme_minimal()




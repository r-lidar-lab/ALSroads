# Parameters to extract the point cloud, process it
# into small sections and computing profiles
extraction = list(
  road_buffer = 160,
  section_length = 10,
  profile_resolution = 0.5
)

# Parameters to detect embankments in the terrain
embankments = list(
  min_slope = 10,
  max_width = 20
)

# Parameters to detect road in the terrain (excluding embankment)
terrain = list(
  max_flat_slope = 5,
  max_sd_ground_points = 0.1,
  max_elevation_ground_points = 0.15
)

# Parameters to detect road in the vegetation
vegetation = list(
  max_percentage_point_above_50cm = 0.25,
  max_percentage_point_above_10cm = 0.25,
  max_stdz = 0.15
)

peak = list(
  lag = 30,
  threshold = 3,
  influence = 1
)

# Parameters to evaluate the state of a section
state = list(
  percentage_veg_thresholds = c(20,70),
  drivable_width_thresholds = c(1,5),
  score_thresholds = c(0.2, 0.5)
)

param = list(
  extraction = extraction,
  embankments = embankments,
  terrain = terrain,
  vegetation = vegetation,
  peak = peak,
  state = state
)

#' Parameters
#'
#' Parameters for function \link{measure_road}. Keeping 'as is' is recommended.
#'
#' **extraction**: parameters used to extract the point cloud, process it into small sections
#' and compute profiles
#'
#' - **road_buffer**: width of the point-cloud extracted around the road
#' - **section_length**: length of sections of road
#' - **profile_resolution**: resolution of the profiles and DTMs
#'
#' **embankments**: parameters to detect embankments in the terrain
#'
#' - **min_slope**: slope (degrees) greater may initiate or terminate the detection of embankments
#' - **max_width**: the distance between two embankments cannot be greater than that.
#' Otherwise it is likely to be only random terrain variations.
#'
#' **terrain**: parameters to detect road based on the terrain aspect (excluding embankment)
#'
#' - **max_flat_slope**: a road is flat, max slope (degree) to be considered flat
#' - **max_sd_ground_points**: point belonging on road have a low dispersion on Z.
#' maximum standard deviation to be considered as low dispersion
#' - **max_elevation_ground_points**: once normalized relatively to the potential road
#' all the points are expected to be at 0. This is the tolerance.
#'
#' **vegetation**: parameters to detect road in the vegetation gaps
#'
#' - **max_percentage_point_above_50cm**: maximum allowed ratio of point above 50 cm
#' - **max_percentage_point_above_10cm**: maximum allowed ratio of point above 10 cm
#' - **max_stdz**: maximum allowed standard deviation on Z
#'
#' **state**: parameters to evaluate the state of a section
#'
#' - **percentage_veg_thresholds**:
#' - **drivable_width_thresholds**:
#' @examples
#' str(mffproads_default_parameters)
#' @export
#' @md
mffproads_default_parameters = param

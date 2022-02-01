# Parameters to extract the point cloud, process it
# into small sections and computing profiles
extraction = list(
  road_max_len = 2000,
  road_max_width = 30,
  road_buffer = 80,
  section_length = 10,
  profile_resolution = 0.5
)

constraint = list(confidence = 0.1)

# Parameters to detect embankments in the terrain
embankments = list(
  min_slope = 10
  #max_width = 20
)

# Parameters to detect road in the terrain (excluding embankment)
terrain = list(
  max_sd_ground_points = 0.1,
  max_elevation_ground_points = 0.15
)

conductivity = list(
  s = c(5, 20),
  r = c(0.2, 0.4),
  e = c(15, 40),
  q = c(0.1, 0.25, 0.5),
  h = c(0.5, 1),
  d = c(0.25, 0.75, 0.95)
)

# Parameters to evaluate the state of a section
state = list(
  percentage_veg_thresholds = c(10,40),
  drivable_width_thresholds = c(1,5),
  conductivity_thresholds = c(0.25, 0.5),
  shoulder_thresholds =  c(50, 75)
)

param = list(
  extraction = extraction,
  constraint = constraint,
  conductivity = conductivity,
  embankments = embankments,
  terrain = terrain,
  state = state
)

#' Parameters
#'
#' Parameters for function \link{measure_road}. Keeping 'as is' is recommended.
#'
#' **extraction**: parameters used to extract the point cloud, process it into small sections
#' and compute profiles
#'
#' - **road_max_len**: maximum size of a processed road. If longer than that it will be split in
#' chunks of equal sizes such as no one is more than this size.
#' - **road_max_width**: maximum width of a road for metrics measurements
#' - **road_buffer**: width of a buffer around the road for point-cloud extraction
#' - **section_length**: length of sections of road for metrics measurements
#' - **profile_resolution**: resolution of the profiles and DTMs for metrics measurements
#'
#' **constraint**: parameters to contraint the road relocation step
#'  - **confidence** The confidence you have on the location of the reference roads. 1 means
#' that the road is 100\% a ground truth. This will skip the relocation step. High values mean high
#' confidence that the reference road is correct and this will help the algorithm. Low values leave
#' more freedom to the algorithm but it becomes also more prone to errors. However this parameter is
#' not very sensitive.
#'
#'
#' **embankments**: parameters to detect embankments in the terrain
#'
#' - **min_slope**: slope (degrees) greater may initiate or terminate the detection of embankments
#'
#' **terrain**: parameters to detect road based on the terrain aspect (excluding embankment)
#'
#' - **max_sd_ground_points**: point belonging on road have a low dispersion on Z.
#' maximum standard deviation to be considered as low dispersion
#' - **max_elevation_ground_points**: once normalized relatively to the potential road
#' all the points are expected to be at 0. This is the tolerance.
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

# plot_time_series_dygraph

## Interactive Time Series Visualization with Optional Dual Y-Axes and Event Highlighting

### Description

`plot_time_series_dygraph` is an R function designed to generate interactive time series plots using the `dygraphs` package. It supports multiple series on both primary and secondary y-axes, customizable colors, and integrated event highlighting. This function is ideal for visualizing complex datasets, such as tri-axial accelerometry data, where different metrics and events need to be clearly distinguished and analyzed over time.

### Features

- **Dual Y-Axes:** Plot multiple time series with different scales on primary and secondary y-axes.
- **Customizable Colors:** Assign specific colors to each series for clear differentiation.
- **Event Highlighting:** Shade and annotate specific events (e.g., behavioural bouts) within the time series for enhanced analysis.
- **Range Selector:** Zoom into specific time intervals with an interactive range selector.
- **Unzoom Plugin:** Easily reset the zoom to view the full data range.
- **Custom Titles and Labels:** Define main and axis titles for better context.
- **Flexible Time Handling:** Use existing timestamp columns or generate timestamps based on sampling rates.

### Required Packages

The function relies on the following R packages. If any are not installed, the function will automatically install them before proceeding.

- **data.table:** For efficient data manipulation.
- **dygraphs:** For creating interactive time series plots.
- **xts:** For handling extensible time series data.
- **htmlwidgets:** For embedding JavaScript widgets.
- **RColorBrewer:** For generating aesthetically pleasing color palettes.

### Installation

Ensure that you have R installed on your system. You can then source the script containing the `plot_time_series_dygraph` function or include it within your R project.

```r
# Example: Sourcing the script
source("path_to_your_script/plot_time_series_dygraph.R")
```

## Inputs

- **`data`**:  
  A `data.frame` or `data.table` containing the dataset to be plotted (the exact columns to be plotted specified below). 

- **`columns_y`**:  
  A character vector specifying the names of the columns to plot on the **primary y-axis** (e.g., for tri-axial accelerometery data, this might be: `c("ax", "ay", "az")`).

- **`columns_y2`** *(Optional)*:  
  A character vector specifying the names of the columns to plot on the **secondary y-axis** (e.g., `c("mx", "my", "mz")`). Use this parameter to plot additional series that have a different scale from those on the primary y-axis. **Note:** Ensure that there is no overlap between `columns_y` and `columns_y2`.

- **`time_col`** *(Optional)*:  
  A character string specifying the name of the timestamp column in your dataset. This column should be of type `POSIXct` or `POSIXlt`. If not provided, the function will require `sampling_rate` to create a timestamp.

- **`sampling_rate`** *(Optional)*:  
  A numeric value indicating the sampling rate in Hertz (samples per second). This parameter is **required** if `time_col` is not provided. It is used to generate a sequence of timestamps assuming a uniform sampling interval.

- **`colors_y`** *(Optional)*:  
  A character vector specifying colors for the series plotted on the **primary y-axis**. The length of this vector should match the number of columns specified in `columns_y`. If not provided, the function assigns default colors automatically.

- **`colors_y2`** *(Optional)*:  
  A character vector specifying colors for the series plotted on the **secondary y-axis**. The length of this vector should match the number of columns specified in `columns_y2`. If not provided, the function assigns default colors automatically.

- **`events_col`** *(Optional)*:  
  A character string specifying the name of the column that contains event markers in your dataset. This column should be numeric or integer, where different integer values represent different event types. Events will be shaded and annotated on the plot based on this column.

- **`event_labels`** *(Optional)*:  
  A character vector specifying labels for the events.  
  - If provided as a **named vector**, the names should correspond to event codes in `events_col`, and the values are the labels.
  - If provided as an **unnamed vector**, labels are assigned based on the sorted order of unique event codes.

- **`event_colors`** *(Optional)*:  
  A character vector specifying colors for the events.  
  - If provided as a **named vector**, the names should correspond to event codes in `events_col`, and the values are the colors.
  - If provided as an **unnamed vector**, colors are assigned based on the sorted order of unique event codes.

- **`event_alpha`** *(Optional)*:  
  A numeric value indicating the transparency level for event shading. Values range from 0 (completely transparent) to 1 (completely opaque). Default is `0.3`.

- **`y_range`** *(Optional)*:  
  A numeric vector of length 2 specifying the fixed range for the **primary y-axis** (e.g., `c(-10, 10)`). If not provided, the axis scales automatically based on the data.

- **`y2_range`** *(Optional)*:  
  A numeric vector of length 2 specifying the fixed range for the **secondary y-axis** (e.g., `c(-5, 15)`). If not provided, the axis scales automatically based on the data.

- **`main_title`** *(Optional)*:  
  A character string for the main title of the plot. Default is `"Time Series Plot"`.

- **`x_label`** *(Optional)*:  
  A character string for the x-axis label. Default is `"Time"`.

- **`y_label`** *(Optional)*:  
  A character string for the **primary y-axis** label. If not provided, the axis label is left blank.

- **`y2_label`** *(Optional)*:  
  A character string for the **secondary y-axis** label. If not provided, the axis label is left blank.

- **`verbose`**:  
  A logical flag indicating whether to print informative messages to the console during function execution. Default is `TRUE`.

  
  ### Example Usage

```R
# Timestamp column
df$timestamp <- ymd_hms(df$timestamp, tz= "UTC") #lubridate function

# Define columns to plot
columns_y <- c("aX", "aY", "aZ", "aX_mean_100", "aY_mean_100", "aZ_mean_100") # Primary y-axis (e.g., raw and smoothed tri-axial acceleration variables)
columns_y2 <- c("VeDBA.sm")  # Secondary y-axis

# Define colors for each axis
colors_y <- c("red", "green", "blue", "#FF9999", "#99FF99", "#9999FF")  # Colours for primary y-axis
colors_y2 <- c("purple") # Colors for secondary y-axis

# Define event labels and colors
event_labels_plain <- c("Low Activity", "High Activity")
event_colors_plain <- c("cyan", "yellow")

df$Marked_events = ifelse(df$VeDBA.sm < 0.1, 1, 2) # Marked event values based on a simple dynamic acceleration threshold

# Apply the plotting function with dual y-axes
dygraph_plot <- plot_time_series_dygraph(
  data = df,
  columns_y = columns_y,
  columns_y2 = columns_y2,
  time_col = "timestamp",  # Ensure 'timestamp' exists and is POSIXct
  sampling_rate = NULL,
  colors_y = colors_y,
  colors_y2 = colors_y2,
  main_title = "Time series Data with Dual Y-Axes",
  events_col = "Marked_events",
  event_labels = event_labels_plain,
  event_colors = event_colors_plain,
  event_alpha = 0.5,
  y_range = NULL,
  y2_range = c(0, 3),  # Example fixed range for secondary y-axis
  verbose = TRUE
)
```

## License

This project is licensed under the MIT License - see the `LICENSE` file for details.

## Contact
# R. Gunner
For bug reports, questions or suggestions, feel free to contact me at [rgunner@ab.mpg].

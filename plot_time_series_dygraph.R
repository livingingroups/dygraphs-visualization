#### plot_time_series_dygraph Function ####

# Purpose:

# Generate an interactive dygraph for multiple time series columns with options fo customizable colors, and integrated plugins for enhanced interactivity, and event highlighting.

# Inputs:

#data: A data.frame or data.table containing the dataset to be plotted (the exact columns to be plotted specified below). 
#columns_y: A character vector specifying the names of the columns to plot on the primary y-axis (e.g., for tri-axial accelerometery data, this might be: c("ax", "ay", "az")).
#columns_y2 (Optional): A character vector specifying the names of the columns to plot on the secondary y-axis (e.g., c("mx", "my", "mz")). Use this parameter to plot additional series that have a different scale from those on the primary y-axis. Note: Ensure that there is no overlap between columns_y and columns_y2.
#time_col (Optional): A character string specifying the name of the timestamp column in your dataset. This column should be of type POSIXct or POSIXlt. If not provided, the function will require sampling_rate to create a timestamp.
#sampling_rate (Optional): A numeric value indicating the sampling rate in Hertz (samples per second). This parameter is required if time_col is not provided. It is used to generate a sequence of timestamps assuming a uniform sampling interval.
#colors_y (Optional): A character vector specifying colors for the series plotted on the primary y-axis. The length of this vector should match the number of columns specified in columns_y. If not provided, the function assigns default colors automatically.
#colors_y2 (Optional): A character vector specifying colors for the series plotted on the secondary y-axis. The length of this vector should match the number of columns specified in columns_y2. If not provided, the function assigns default colors automatically.
#events_col (Optional): A character string specifying the name of the column that contains event markers in your dataset. This column should be numeric or integer, where different integer values represent different event types. Events will be shaded and annotated on the plot based on this column.
#event_labels (Optional): A character vector specifying labels for the events. If provided as a named vector, the names should correspond to event codes in events_col, and the values are the labels. If provided as an unnamed vector, labels are assigned based on the sorted order of unique event codes.
#event_colors (Optional): A character vector specifying colors for the events. If provided as a named vector, the names should correspond to event codes in events_col, and the values are the colors. If provided as an unnamed vector, colors are assigned based on the sorted order of unique event codes.
#event_alpha (Optional): A numeric value indicating the transparency level for event shading. Values range from 0 (completely transparent) to 1 (completely opaque). Default is 0.3.
#y_range (Optional): A numeric vector of length 2 specifying the fixed range for the primary y-axis (e.g., c(-10, 10)). If not provided, the axis scales automatically based on the data.
#y2_range (Optional): A numeric vector of length 2 specifying the fixed range for the secondary y-axis (e.g., c(-5, 15)). If not provided, the axis scales automatically based on the data.
#main_title (Optional): A character string for the main title of the plot. Default is "Time Series Plot".
#x_label (Optional): A character string for the x-axis label. Default is "Time".
#y_label (Optional): A character string for the primary y-axis label. If not provided, the axis label is left blank.
#y2_label (Optional): A character string for the secondary y-axis label. If not provided, the axis label is left blank.
#verbose: A logical flag indicating whether to print informative messages to the console during function execution. Default = TRUE.

# Outputs:

# An interactive dygraph object displaying the specified time series with integrated plugins and event demarcation.

# Define the plotting function
plot_time_series_dygraph <- function(data,
                                     columns_y,
                                     columns_y2 = NULL,
                                     time_col = NULL,
                                     sampling_rate = NULL,
                                     colors_y = NULL,
                                     colors_y2 = NULL,
                                     events_col = NULL,
                                     event_labels = NULL,
                                     event_colors = NULL,
                                     event_alpha = 0.3,  
                                     y_range = NULL,
                                     y2_range = NULL,
                                     main_title = "Time Series Plot",
                                     x_label = "Time",
                                     y_label = NULL,
                                     y2_label = NULL,
                                     verbose = TRUE) {
  
  # -------------------------
  # 0. Load Necessary Packages
  # -------------------------
  
  required_packages <- c("data.table", "dygraphs", "xts", "htmlwidgets", "RColorBrewer")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
  
  # -------------------------
  # 1. Define the Unzoom Plugin
  # -------------------------
  
  dyUnzoom <- function(dygraph_obj, plugin_path = NULL) {
    # Check if plugin_path is provided; if not, use default path
    if (is.null(plugin_path)) {
      plugin_path <- system.file("plugins/unzoom.js", package = "dygraphs")
      if (plugin_path == "") {
        stop("Unzoom plugin 'unzoom.js' not found in the dygraphs package.")
      }
    } else {
      if (!file.exists(plugin_path)) {
        stop(paste("Specified plugin_path does not exist:", plugin_path))
      }
    }
    
    dyPlugin(
      dygraph = dygraph_obj,
      name = "Unzoom",
      path = plugin_path,
      options = list()
    )
  }
  
  # -------------------------
  # 2. Define the Custom Value Formatter JavaScript Function
  # -------------------------
  
  CustomValueFormat <- '
  function (ms) {
    var d = new Date(ms);
    return Dygraph.zeropad(d.getDate()) + "/" + 
           Dygraph.zeropad(d.getMonth() + 1) + "/" +
           Dygraph.zeropad(d.getFullYear()) + " " + 
           Dygraph.zeropad(d.getHours()) + ":" + 
           Dygraph.zeropad(d.getMinutes()) + ":" + 
           Dygraph.zeropad(d.getSeconds()) + "." + 
           Dygraph.zeropad(d.getMilliseconds());
  }'
  
  # -------------------------
  # 3. Input Validation
  # -------------------------
  
  # Validate 'data' input
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or data.table.")
  }
  
  # Convert to data.table if not already
  if (!is.data.table(data)) {
    data <- as.data.table(data)
  }
  
  # Validate 'columns_y'
  if (!is.character(columns_y) | length(columns_y) < 1) {
    stop("'columns_y' must be a non-empty character vector specifying column names.")
  }
  
  if (!all(columns_y %in% names(data))) {
    missing_cols_y <- columns_y[!columns_y %in% names(data)]
    stop(paste("The following 'columns_y' are not present in 'data':", paste(missing_cols_y, collapse = ", ")))
  }
  
  # Validate 'columns_y2' if provided
  if (!is.null(columns_y2)) {
    if (!is.character(columns_y2) | length(columns_y2) < 1) {
      stop("'columns_y2' must be a non-empty character vector specifying column names.")
    }
    
    if (!all(columns_y2 %in% names(data))) {
      missing_cols_y2 <- columns_y2[!columns_y2 %in% names(data)]
      stop(paste("The following 'columns_y2' are not present in 'data':", paste(missing_cols_y2, collapse = ", ")))
    }
    
    # Ensure no overlap between columns_y and columns_y2
    overlap_cols <- intersect(columns_y, columns_y2)
    if (length(overlap_cols) > 0) {
      stop(paste("The following columns are specified in both 'columns_y' and 'columns_y2':", 
                 paste(overlap_cols, collapse = ", ")))
    }
  }
  
  # Validate 'time_col' if provided
  if (!is.null(time_col)) {
    if (!is.character(time_col) | length(time_col) != 1) {
      stop("'time_col' must be a single character string specifying the timestamp column name.")
    }
    
    if (!(time_col %in% names(data))) {
      stop(paste("The specified 'time_col' does not exist in 'data':", time_col))
    }
    
    # Check if 'time_col' is a proper timestamp object
    if (!inherits(data[[time_col]], c("POSIXct", "POSIXlt"))) {
      datetime_converted <- as.POSIXct(data[[time_col]], tz = "UTC")
      if (any(is.na(datetime_converted))) {
        stop("'time_col' must be a POSIXct/POSIXlt datetime object or convertible to one.", call. = FALSE)
      } else {
        data[[time_col]] <- datetime_converted
      }
    }
    
    # Ensure time is strictly increasing
    if (any(diff(data[[time_col]]) <= 0)) {
      stop("'time_col' must be strictly increasing.")
    }
  } else {
    # If 'time_col' is not provided, 'sampling_rate' must be provided
    if (is.null(sampling_rate)) {
      stop("Either 'time_col' must be provided to estimate the sampling rate or 'sampling_rate' must be specified.")
    }
    
    # Validate 'sampling_rate'
    if (!is.numeric(sampling_rate) | length(sampling_rate) != 1 | sampling_rate <= 0) {
      stop("'sampling_rate' must be a single positive numeric value.")
    }
  }
  
  # Validate 'colors_y' if provided
  num_series_y <- length(columns_y)
  if (!is.null(colors_y)) {
    if (!is.character(colors_y)) {
      stop("'colors_y' must be a character vector specifying colors.")
    }
    if (length(colors_y) != num_series_y) {
      stop("Length of 'colors_y' must match the number of 'columns_y' to plot.")
    }
  } else {
    # Assign default colors using RColorBrewer
    if (num_series_y <= 12) {
      colors_y <- brewer.pal(max(3, min(num_series_y, 12)), "Set3")
    } else {
      colors_y <- colorRampPalette(brewer.pal(12, "Set3"))(num_series_y)
    }
  }
  
  # Validate 'colors_y2' if provided
  if (!is.null(columns_y2)) {
    num_series_y2 <- length(columns_y2)
    if (!is.null(colors_y2)) {
      if (!is.character(colors_y2)) {
        stop("'colors_y2' must be a character vector specifying colors.")
      }
      if (length(colors_y2) != num_series_y2) {
        stop("Length of 'colors_y2' must match the number of 'columns_y2' to plot.")
      }
    } else {
      # Assign default colors using RColorBrewer
      if (num_series_y2 <= 12) {
        colors_y2 <- brewer.pal(max(3, min(num_series_y2, 12)), "Set1")
      } else {
        colors_y2 <- colorRampPalette(brewer.pal(12, "Set1"))(num_series_y2)
      }
    }
  }
  
  # Validate 'events_col' if provided
  if (!is.null(events_col)) {
    if (!is.character(events_col) | length(events_col) != 1) {
      stop("'events_col' must be a single character string specifying the events column name.")
    }
    
    if (!(events_col %in% names(data))) {
      stop(paste("The specified 'events_col' does not exist in 'data':", events_col))
    }
    
    if (!is.numeric(data[[events_col]]) & !is.integer(data[[events_col]])) {
      stop("'events_col' must be a numeric or integer vector representing event markers.")
    }
  }
  
  # Validate 'y_range' and 'y2_range' if provided
  if (!is.null(y_range)) {
    if (!is.numeric(y_range) | length(y_range) != 2) {
      stop("'y_range' must be a numeric vector of length 2 specifying the range for the primary y-axis.")
    }
  }
  
  if (!is.null(columns_y2) && !is.null(y2_range)) {
    if (!is.numeric(y2_range) | length(y2_range) != 2) {
      stop("'y2_range' must be a numeric vector of length 2 specifying the range for the secondary y-axis.")
    }
  }
  
  # -------------------------
  # 4. Handle Time Column and Sampling Rate
  # -------------------------
  
  if (!is.null(time_col)) {
    # Inform the user about sampling rate estimation
    if (verbose) {
      # Calculate sampling rate based on mean difference to handle irregularities
      time_diffs <- as.numeric(diff(data[[time_col]]))
      mean_diff <- mean(time_diffs)
      estimated_sampling_rate <- 1 / mean_diff
      message(sprintf("Estimated Sampling Rate from 'time_col': %.2f Hz", estimated_sampling_rate))
    }
  } else {
    # Create a timestamp column based on 'sampling_rate'
    n <- nrow(data)
    time_vector <- seq(0, by = 1 / sampling_rate, length.out = n)
    
    # Add as POSIXct assuming origin = "1970-01-01"
    time_vector <- as.POSIXct(time_vector, origin = "1970-01-01", tz = "UTC")
    data[, timestamp_auto := time_vector]
    time_col <- "timestamp_auto"
    
    # Inform the user about sampling rate usage
    if (verbose) {
      message(sprintf("Using provided Sampling Rate: %.2f Hz to create timestamp.", sampling_rate))
    }
  }
  
  # -------------------------
  # 5. Prepare Data for dygraph
  # -------------------------
  
  # Select relevant columns for primary y-axis
  plot_data_y <- data[, ..columns_y]
  
  # Select relevant columns for secondary y-axis if provided
  if (!is.null(columns_y2)) {
    plot_data_y2 <- data[, ..columns_y2]
  }
  
  # Combine all columns into one data.table
  if (!is.null(columns_y2)) {
    combined_data <- cbind(data[, ..time_col], plot_data_y, plot_data_y2)
  } else {
    combined_data <- cbind(data[, ..time_col], plot_data_y)
  }
  
  # Convert to xts object
  if (!is.null(columns_y2)) {
    plot_xts <- xts(combined_data[, -1, with = FALSE], order.by = combined_data[[time_col]])
  } else {
    plot_xts <- xts(combined_data[, -1, with = FALSE], order.by = combined_data[[time_col]])
  }
  
  # -------------------------
  # 6. Create dygraph
  # -------------------------
  
  dygraph_obj <- dygraph(plot_xts, main = main_title) %>%
    dyOptions(
      drawPoints = FALSE, 
      strokeWidth = 1.5, 
      useDataTimezone = TRUE,
      sigFigs = 3,
      axisLineWidth = 1.5, 
      rightGap = 90,
      drawGrid = FALSE
    ) %>%
    dyAxis("y", 
           label = y_label, 
           axisLabelWidth = 80, 
           labelWidth = 30, 
           axisLabelFontSize = 14, 
           drawGrid = FALSE,
           valueRange = y_range
    ) %>%
    dyAxis("x", 
           valueFormatter = JS(CustomValueFormat), 
           ticker = "Dygraph.dateTicker", 
           label = x_label, 
           axisLabelWidth = 80, 
           labelWidth = 30, 
           labelHeight = 30,
           axisLabelFontSize = 14, 
           drawGrid = FALSE
    ) %>%
    dyRangeSelector(retainDateWindow = TRUE) %>%
    dyLegend(show = "onmouseover",
             width = 200,
             hideOnMouseOut = TRUE, 
             labelsSeparateLines = TRUE) %>% 
    dyUnzoom()
  
  # Add dySeries for primary y-axis columns with specified colors
  for (i in seq_along(columns_y)) {
    col_name <- columns_y[i]
    dygraph_obj <- dygraph_obj %>%
      dySeries(col_name, label = col_name, color = colors_y[i])
  }
  
  # Add dySeries for secondary y-axis columns with specified colors
  if (!is.null(columns_y2)) {
    for (i in seq_along(columns_y2)) {
      col_name <- columns_y2[i]
      dygraph_obj <- dygraph_obj %>%
        dySeries(col_name, label = col_name, color = colors_y2[i], axis = "y2")
    }
    
    # Configure the secondary y-axis
    dygraph_obj <- dygraph_obj %>%
      dyAxis("y2",
             label = y2_label, 
             axisLabelWidth = 80, 
             labelWidth = 30, 
             axisLabelFontSize = 14, 
             drawGrid = FALSE,
             independentTicks = TRUE,
             valueRange = y2_range
      )
  }
  
  # -------------------------
  # 7. Handle Marked Events with dyShading and dyEvent
  # -------------------------
  
  if (!is.null(events_col)) {
    # Identify unique events excluding 0
    unique_events <- sort(unique(data[[events_col]][data[[events_col]] != 0]))
    
    if (length(unique_events) == 0) {
      if (verbose) {
        message("No events found in 'events_col'. Skipping event shading and labeling.")
      }
    } else {
      # Function to add alpha to colors
      addalpha <- function(col, alpha = 0.3) {
        rgb_val <- col2rgb(col) / 255
        rgb(rgb_val[1], rgb_val[2], rgb_val[3], alpha = alpha)
      }
      
      # Assign labels and colors based on user input or defaults
      if (!is.null(event_labels)) {
        if (!is.null(names(event_labels))) {
          # Named vector: map names to event codes
          event_labels_map <- event_labels
        } else {
          # Unnamed vector: assign labels based on sorted unique events
          if (length(event_labels) != length(unique_events)) {
            stop("Length of 'event_labels' must match the number of unique events.")
          }
          event_labels_map <- setNames(event_labels, unique_events)
        }
      } else {
        # Default labels as event codes
        event_labels_map <- setNames(as.character(unique_events), unique_events)
      }
      
      if (!is.null(event_colors)) {
        if (!is.null(names(event_colors))) {
          # Named vector: map names to event codes
          event_colors_map <- event_colors
        } else {
          # Unnamed vector: assign colors based on sorted unique events
          if (length(event_colors) != length(unique_events)) {
            stop("Length of 'event_colors' must match the number of unique events.")
          }
          event_colors_map <- setNames(event_colors, unique_events)
        }
        # Apply alpha transparency
        event_colors_map <- sapply(event_colors_map, function(col) addalpha(col, alpha = event_alpha), USE.NAMES = TRUE)
      } else {
        # Assign default semi-transparent colors using RColorBrewer
        if (length(unique_events) <= 12) {
          default_event_colors <- brewer.pal(max(3, min(length(unique_events), 12)), "Set1")
        } else {
          default_event_colors <- colorRampPalette(brewer.pal(12, "Set1"))(length(unique_events))
        }
        # Add alpha for semi-transparency
        default_event_colors <- sapply(default_event_colors, function(col) addalpha(col, alpha = event_alpha))
        # Map to event integers
        event_colors_map <- setNames(default_event_colors[1:length(unique_events)], unique_events)
      }
      
      # Iterate over each unique event
      for (event in unique_events) {
        # Get the indices where the event occurs
        event_indices <- which(data[[events_col]] == event)
        
        if (length(event_indices) == 0) next  # Skip if no occurrence
        
        # Identify contiguous blocks of events
        rle_event <- rle(diff(event_indices) == 1)
        lengths <- rle_event$lengths
        values <- rle_event$values
        
        # Initialize start index
        start_idx <- event_indices[1]
        shading_events <- list()
        
        for (i in seq_along(lengths)) {
          if (values[i]) {
            # Continuation of a contiguous block
            # No action needed; the end_idx will be updated in the next iteration
          } else {
            # End of a contiguous block
            end_idx <- event_indices[sum(lengths[1:i])]
            shading_events[[length(shading_events) + 1]] <- c(start_idx, end_idx)
            if (i < length(lengths)) {
              start_idx <- event_indices[sum(lengths[1:i]) + 1]
            }
          }
        }
        
        # Handle the last block if it's contiguous
        if (length(values) > 0 && values[length(values)]) {
          end_idx <- event_indices[length(event_indices)]
          shading_events[[length(shading_events) + 1]] <- c(start_idx, end_idx)
        }
        
        # Apply dyShading and dyEvent for each shading event
        for (shading in shading_events) {
          start_time <- data[[time_col]][shading[1]]
          end_time <- data[[time_col]][shading[2]]
          
          # Apply dyShading with the specified color and alpha
          dygraph_obj <- dygraph_obj %>%
            dyShading(from = start_time, to = end_time, color = event_colors_map[[as.character(event)]])
          
          # Add dyEvent at the start of the event with correct label
          dygraph_obj <- dygraph_obj %>%
            dyEvent(
              x = start_time, 
              label = as.character(event_labels_map[[as.character(event)]]),  # Ensure label is a string
              labelLoc = "bottom",
              color = "grey", 
              strokePattern = "dashed"
            )
        }
      }
    }
  }
  
  
  # -------------------------
  # 9. Return the Result
  # -------------------------
  print(dygraph_obj)
  return(dygraph_obj)
  
}
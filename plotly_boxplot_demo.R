# Get this figure: fig <- py$get_figure("etpinard", 199)
# Get this figure's data: data <- py$get_figure("etpinard", 199)$data
# Add data to this figure: py$plotly(list(x=c(4, 5), y=c(4, 5)), kwargs=list(filename="s4_grades-box-styled2", fileopt="extend"))
# Get y data of first trace: y1 <- py$get_figure("etpinard", 199)$data[[1]]$y

# Get figure documentation: https://plot.ly/r/get-requests/
# Add data documentation: https://plot.ly/r/file-options/

# You can reproduce this figure in R with the following code!

# Learn about API authentication here: https://plot.ly/r/getting-started
# Find your api_key here: https://plot.ly/settings/api

library(plotly)
py <- plotly(username='llsquare', key='90xamoida4')
trace1 <- list(
  y = c(90.0, 88.0, 55.0, 88.0, 72.0, 100.0, 88.0, 25.0, 92.0, 100.0, 82.0, 82.0, 90.0, 68.0, 85.0, "", 82.0, 40.0, 100.0, 92.0, 82.0, 55.0, 62.0, 85.0, 100.0, 75.0, 88.0, 78.0, 80.0, "", 92.0, 100.0, 88.0, 72.0, 95.0, 80.0, 90.0, 72.0, 100.0, "", 100.0, 75.0, 82.0, 60.0, 90.0, 85.0, 90.0, 38.0, 78.0, 82.0, 100.0, 90.0, 80.0, "", 80.0, 100.0, 70.0, 100.0, 82.0, 62.0, 92.0, "", 100.0, 80.0, "", 100.0, 88.0, 85.0), 
  boxmean = "sd", 
  boxpoints = "all", 
  jitter = 0.5, 
  line = list(color = "#1c9099"), 
  marker = list(
    color = "#feb24c", 
    line = list(
      color = "#FFFFFF", 
      width = 1
    ), 
    size = 10
  ), 
  name = "Homeworks", 
  pointpos = -2, 
  type = "box"
)
trace2 <- list(
  y = c(70.0, 65.0, 85.0, 75.0, 72.0, 75.0, 90.0, 88.0, 85.0, 80.0, 92.0, 85.0, 85.0, 75.0, 72.0, "", 80.0, 42.0, 80.0, 95.0, 90.0, 62.0, 65.0, 65.0, 82.0, 68.0, 48.0, 57.0, 95.0, 70.0, 100.0, 80.0, 95.0, 78.0, 80.0, 80.0, 85.0, 90.0, 100.0, 52.0, 85.0, 72.0, 70.0, 45.0, 75.0, 85.0, 95.0, 65.0, 70.0, 85.0, 70.0, 85.0, 35.0, "", 90.0, 95.0, 95.0, 65.0, 62.0, 48.0, 60.0, "", 85.0, 85.0, "", 90.0, 70.0, 68.0), 
  boxmean = "sd", 
  boxpoints = "all", 
  jitter = 0.5, 
  line = list(color = "#1c9099"), 
  marker = list(
    color = "#feb24c", 
    line = list(
      color = "#FFFFFF", 
      width = 1
    ), 
    size = 10
  ), 
  name = "Midterm Exam", 
  pointpos = -2, 
  type = "box"
)
trace3 <- list(
  y = c(95.0, 75.0, 70.0, 72.0, 52.0, 70.0, 82.0, 90.0, 95.0, 80.0, 68.0, 88.0, 82.0, 52.0, 80.0, "", 78.0, 57.0, 88.0, 88.0, 100.0, 50.0, 65.0, 78.0, 92.0, 65.0, 50.0, 60.0, 88.0, "", 100.0, 50.0, 90.0, 70.0, 60.0, 72.0, 75.0, 95.0, 100.0, 45.0, 68.0, 72.0, 45.0, 60.0, 78.0, 85.0, 92.0, 45.0, 68.0, 70.0, 85.0, 82.0, 62.0, "", 75.0, 100.0, 80.0, 65.0, 52.0, 48.0, 57.0, "", 100.0, 72.0, "", 100.0, 80.0, 65.0), 
  boxmean = "sd", 
  boxpoints = "all", 
  jitter = 0.5, 
  line = list(color = "#1c9099"), 
  marker = list(
    color = "#feb24c", 
    line = list(
      color = "#FFFFFF", 
      width = 1
    ), 
    size = 10
  ), 
  name = "Final Exam", 
  pointpos = -2, 
  type = "box"
)
trace4 <- list(
  y = c(86.0, 75.9, 70.0, 77.7, 64.0, 80.5, 86.2, 69.9, 91.1, 86.0, 79.4, 85.3, 85.3, 63.7, 79.1, "", 79.8, 47.4, 89.2, 91.3, 91.6, 55.1, 64.1, 76.2, 91.4, 68.9, 60.8, 64.5, 87.7, 21.0, 97.6, 74.0, 90.9, 73.0, 76.5, 76.8, 82.5, 86.6, 100.0, 33.6, 82.7, 72.9, 63.6, 55.5, 80.7, 85.0, 92.3, 48.9, 71.6, 78.1, 85.0, 85.3, 59.3, "", 81.0, 98.5, 81.5, 75.5, 64.0, 52.2, 68.4, "", 95.5, 78.3, "", 97.0, 79.4, 71.9), 
  boxmean = "sd", 
  boxpoints = "all", 
  jitter = 0.5, 
  line = list(color = "#1c9099"), 
  marker = list(
    color = "#feb24c", 
    line = list(
      color = "#FFFFFF", 
      width = 1
    ), 
    size = 10
  ), 
  name = "Course Grade", 
  pointpos = -2, 
  type = "box"
)
data <- list(trace1, trace2, trace3, trace4)
layout <- list(
  autosize = FALSE, 
  height = 500, 
  plot_bgcolor = "#EFECEA", 
  showlegend = FALSE, 
  title = "Fig 4.4c: Course Grade Distributions", 
  width = 650, 
  xaxis = list(
    gridcolor = "#FFFFFF", 
    showgrid = TRUE, 
    ticklen = 8, 
    ticks = "outside", 
    tickwidth = 1.5, 
    zeroline = FALSE
  ), 
  yaxis = list(
    gridcolor = "#FFFFFF", 
    showgrid = TRUE, 
    ticklen = 8, 
    ticks = "outside", 
    tickwidth = 1.5, 
    title = "Grade [%]", 
    zeroline = FALSE
  )
)
response <- py$plotly(data, kwargs=list(layout=layout))
url <- response$url
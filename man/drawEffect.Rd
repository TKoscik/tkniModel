\name{drawEffect}
\alias{drawEffect}
\title{Visualize Statistical Effects with Random Effect Adjustment}
\description{
  Generates publication-quality scatter plots or point-range plots for main effects
  and interactions (up to 4-way). For mixed-effects models (\code{lmer}), the function
  automatically "partials out" random effects to visualize fixed-effect relationships
  adjusted for the model's hierarchical structure.
}
\usage{
drawEffect(model, effect_name, axes = NULL, labels = NULL,
           plot_colors = get.color(), save_dir, file_name,
           img_format = "png", img_w = 10, img_unit = "cm",
           img_dpi = 600, save_plot = TRUE, return_plot = FALSE)
}
\arguments{
  \item{model}{A fitted model object (e.g., \code{lm}, \code{glm}, or \code{merModLmerTest}).}
  \item{effect_name}{Character. The term in the model to visualize (e.g., \code{"age"} or \code{"age:group"}).}
  \item{axes}{Character vector. Optional custom mapping of variables to plot aesthetics: \code{c(x_axis, y_axis, color, facet_row, facet_col)}.}
  \item{labels}{Character vector. Labels for the plot axes and legends.}
  \item{plot_colors}{Character vector. Colors for the plot levels.}
  \item{save_dir}{Character. Directory to save the output plot.}
  \item{file_name}{Character. Name for the saved image file.}
  \item{img_format}{Character. Image format (e.g., \code{"png"}, \code{"pdf"}).}
  \item{img_w}{Numeric. Width of the image. Height is calculated dynamically based on plot complexity.}
  \item{img_unit}{Character. Units for dimensions (\code{"cm"}, \code{"in"}, \code{"mm"}).}
  \item{img_dpi}{Integer. Resolution of the output.}
  \item{save_plot}{Logical. Whether to write the file to disk.}
  \item{return_plot}{Logical. Whether to return the ggplot object.}
}
\details{
  \bold{Partialling Out Random Effects:}
  When a mixed-effects model (\code{lmer}) is provided, the function adjusts the dependent
  variable by subtracting the Best Linear Unbiased Predictors (BLUPs) of the random effects.
  This visualizes the data points relative to the fixed effects, effectively showing
  the relationship after accounting for subject-level or group-level variance.

  \bold{Interaction Handling:}
  The function automatically handles up to 4-way interactions by mapping variables
  to the X-axis, colors, and facets (grid/wrap) based on the order specified in \code{effect_name}.
}
\value{
  Invisibly returns the path to the saved image file, or the \code{ggplot} object
  if \code{return_plot = TRUE}.
}
\examples{
\dontrun{
# For a simple linear model
drawEffect(model = fit, effect_name = "age", save_dir = "plots/")

# For a mixed-effects interaction
drawEffect(model = lmer_fit, effect_name = "age:group",
           plot_colors = c("red", "blue"), save_dir = "plots/")
}
}

# Assume ldata is loaded as in the file

load("/Users/jeremiahwala/Dropbox_HMS/2025-Amiryousefi-Wala/Analysis/Data/listdatacysift2_Phen_Ki_BT_TCF1_tDBSCAN_BTinteract.RData")
library(data.table)   
library(ggplot2)     
library(bitops)      

### also load directly from CSV to get the cflags
lf_prostate = list.files("~/Dropbox_HMS/2025-Amiryousefi-Wala/Analysis/Data/cyf_files2/csv/", full.names=TRUE, pattern = "csv$")

load_prostate <- function(x) {
  cat("loading", basename(x), "\n")
  dt <- fread(cmd = paste("cut -f1-6 -d','", x))
  attr(dt, "path") <- x
  attr(dt, "Slide_ID") <- sub(".*(LSP\\d+).*", "\\1", x)
  dt
}
dts_prostate <- parallel::mclapply(lf_prostate, load_prostate, mc.cores = 8)
names(dts_prostate) <- sub(".*(LSP\\d+).*", "\\1", basename(lf_prostate))

fcheck <- function(pf, x) {
  return (bitops::bitAnd(as.numeric(pf), as.numeric(x))==as.numeric(x))
}

# ggplot theme
theme_jw <- function() {
  theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black"),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA) ,
      text = element_text(family = "Helvetica"),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line()
    )
}

gleason_colors=c("Low"="#1b9e77","High"="#7570b3")

# ---------------------------------------------------------------------------
# roi_area()
#
# Load an ROI CSV (one polygon per row) and sum the area of every region whose
# label matches a keyword (case insensitive), returned in square millimeters.
#
# File format assumptions (verified against LSP12601_roi.csv):
#   - One ROI per row, with a header.
#   - Polygon vertices are in the 'all_points' column as space-separated
#     "X,Y" pairs, e.g. "22312.00,3060.00 22308.00,3058.00 ...".
#   - The human-readable label is in the 'Text' column (e.g. "Invasive tumor
#     4+4"); 'Name' is a fallback if 'Text' is missing.
#   - Coordinates are in pixels. Polygons are not explicitly closed (the last
#     vertex is not a repeat of the first); the shoelace formula closes them.
#
# Arguments:
#   path        path to the ROI csv
#   keyword     substring to match in the label, case insensitive (e.g. "tumor")
#   microns_per_pixel  scale factor applied to X and Y before area calc.
#                      e.g. 0.65 means each pixel is 0.65 um. Default 1.0
#                      (i.e. coordinates already in microns / no scaling).
#   label_cols  which columns to search for the keyword, in priority order.
#   types       polygon-like ROI types to include (others, e.g. ellipses or
#               points, are skipped because the shoelace area does not apply).
#   verbose     if TRUE, print each matched ROI and its individual area.
#
# Returns:
#   A single numeric: total area in mm^2 of all matched polygons. Returns 0 if
#   nothing matches.
# ---------------------------------------------------------------------------
roi_area <- function(path,
                     keyword,
                     microns_per_pixel = 1.0,
                     label_cols = c("Text", "Name"),
                     types = c("Polygon"),
                     verbose = FALSE) {
  
  # ---- argument checks ----
  if (!file.exists(path)) stop("File not found: ", path)
  if (!is.character(keyword) || length(keyword) != 1L || !nzchar(keyword))
    stop("'keyword' must be a single non-empty string.")
  if (!is.numeric(microns_per_pixel) || length(microns_per_pixel) != 1L ||
      microns_per_pixel <= 0)
    stop("'microns_per_pixel' must be a single positive number.")
  
  roi <- fread(path)
  
  # locate the label column to match against (first available from label_cols)
  label_col <- label_cols[label_cols %in% names(roi)][1]
  if (is.na(label_col))
    stop("None of the label columns (", paste(label_cols, collapse = ", "),
         ") are present in the file.")
  
  if (!"all_points" %in% names(roi))
    stop("Expected an 'all_points' column with the polygon vertices.")
  
  # ---- select matching rows ----
  # keyword match: case insensitive, fixed substring (not a regex), so a
  # keyword like "tumor 4+4" is treated literally rather than as a pattern.
  is_match <- grepl(tolower(keyword), tolower(roi[[label_col]]), fixed = TRUE)
  
  # restrict to polygon-like types if a 'type' column exists
  if ("type" %in% names(roi)) {
    is_match <- is_match & (tolower(roi$type) %in% tolower(types))
  }
  
  idx <- which(is_match)
  if (length(idx) == 0L) {
    if (verbose) cat("No ROIs matched keyword '", keyword, "'.\n", sep = "")
    return(0)
  }
  
  # um^2 per pixel^2: area scales by the square of the linear scale factor
  um2_per_px2 <- microns_per_pixel^2
  UM2_PER_MM2 <- 1e6   # (1000 um)^2 per mm^2
  
  # ---- shoelace area for one "X,Y X,Y ..." string, in pixel^2 ----
  polygon_area_px2 <- function(points_str) {
    toks <- strsplit(trimws(points_str), "\\s+")[[1]]
    toks <- toks[nzchar(toks)]
    if (length(toks) < 3L) return(0)   # need >= 3 vertices for an area
    xy <- vapply(toks, function(p) as.numeric(strsplit(p, ",")[[1]]),
                 numeric(2))            # 2 x N matrix: row1 = X, row2 = Y
    x <- xy[1, ]
    y <- xy[2, ]
    n <- length(x)
    j <- c(n, seq_len(n - 1L))          # previous-vertex index (wraps)
    abs(sum(x * y[c(seq_len(n - 1L) + 1L, 1L)] -
              x[c(seq_len(n - 1L) + 1L, 1L)] * y)) / 2
  }
  
  areas_mm2 <- numeric(length(idx))
  for (k in seq_along(idx)) {
    i <- idx[k]
    a_px2 <- polygon_area_px2(roi[["all_points"]][i])
    areas_mm2[k] <- a_px2 * um2_per_px2 / UM2_PER_MM2
    if (verbose)
      cat(sprintf("  match: %-30s area = %.4f mm^2\n",
                  roi[[label_col]][i], areas_mm2[k]))
  }
  
  total <- sum(areas_mm2)
  if (verbose)
    cat(sprintf("Total over %d ROI(s) matching '%s': %.4f mm^2\n",
                length(idx), keyword, total))
  total
}

# # function to estimate the area of a tumor from a collection of cells (as points) that have been 
# # labeled as containing tumor
cycif_disc_union_area <- function(
    dt,
    radius = 20,
    method = c("mc", "grid"),
    # MC controls
    samples = 200000,            # fewer = coarser, faster, less memory
    chunk_size = 50000,          # processed at a time to limit memory
    seed = NULL,
    # Grid controls (only used if method="grid" or NN pkgs unavailable)
    cell_size = NULL             # default: radius (very coarse, very light)
) {
  # ---- checks ----
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Package 'data.table' is required.")
  if (!data.table::is.data.table(dt)) stop("'dt' must be a data.table.")
  if (!all(c("x","y") %in% names(dt))) stop("Input must contain columns 'x' and 'y'.")
  if (!is.numeric(dt[["x"]]) || !is.numeric(dt[["y"]])) stop("'x' and 'y' must be numeric.")
  if (!is.numeric(radius) || length(radius) != 1L || radius <= 0) stop("'radius' must be a positive number.")

  method <- match.arg(method)
  dat <- data.table::copy(dt)[stats::complete.cases(x, y)]
  if (nrow(dat) == 0L) stop("No rows with non-missing 'x' and 'y'.")

  # Bounding box expanded so disks are fully contained
  xr <- range(dat$x); yr <- range(dat$y)
  xmin <- xr[1] - radius; xmax <- xr[2] + radius
  ymin <- yr[1] - radius; ymax <- yr[2] + radius
  box_area <- (xmax - xmin) * (ymax - ymin)

  # Helper: simple return wrapper
  ret <- function(area, method_used, extra = list()) {
    c(list(area = as.numeric(area),
           method_used = method_used,
           radius = radius,
           bounding_box = c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
           box_area = box_area),
      extra)
  }

  # ---- Monte Carlo estimator (preferred; very light) ----
  if (method == "mc") {
    engine <- if (requireNamespace("RANN", quietly = TRUE)) {
      "RANN"
    } else if (requireNamespace("FNN", quietly = TRUE)) {
      "FNN"
    } else {
      "grid-fallback"  # No NN package; fall back to grid
    }

    if (engine != "grid-fallback") {
      if (!is.null(seed)) set.seed(seed)
      if (!is.numeric(samples) || samples < 1) stop("'samples' must be >= 1.")
      if (!is.numeric(chunk_size) || chunk_size < 1) stop("'chunk_size' must be >= 1.")
      X <- as.matrix(dat[, .(x, y)])

      inside <- 0L
      done <- 0L
      while (done < samples) {
        m <- min(chunk_size, samples - done)
        Q <- cbind(stats::runif(m, xmin, xmax),
                   stats::runif(m, ymin, ymax))
        d <- if (engine == "RANN") {
          RANN::nn2(X, Q, k = 1)$nn.dists[, 1]
        } else {
          FNN::get.knnx(X, Q, k = 1)$nn.dist[, 1]
        }
        inside <- inside + sum(d <= radius)
        done <- done + m
      }
      p <- inside / samples
      area <- box_area * p

      # Rough 95% CI on area (binomial → proportion → scale by box_area)
      se_p <- sqrt(max(p * (1 - p), 1e-12) / samples)
      ci_p <- p + c(-1, 1) * 1.96 * se_p
      ci_area <- pmax(0, pmin(1, ci_p)) * box_area

      return(ret(area, paste0("mc/", engine),
                 list(coverage_fraction = p,
                      samples_used = samples,
                      ci95_area = ci_area)))
    }

    # fall through to grid method when no NN package is available
    method <- "grid"
  }

  # ---- Ultra-light grid dilation (approximates circles with squares) ----
  # Bias: tends to overestimate vs true circles (coarser = faster = more bias).
  if (is.null(cell_size)) cell_size <- radius # very coarse, minimal RAM
  if (!is.numeric(cell_size) || length(cell_size) != 1L || cell_size <= 0)
    stop("'cell_size' must be a positive number.")

  s <- ceiling(radius / cell_size)  # dilation radius in grid cells (Chebyshev metric)
  # Anchor grid to data min to keep indices small
  x0 <- min(dat$x); y0 <- min(dat$y)
  gx <- floor((dat$x - x0) / cell_size)
  gy <- floor((dat$y - y0) / cell_size)

  # Unique seed cells (points may share cells)
  seeds <- unique(data.table::data.table(gx = gx, gy = gy))

  # Record covered cells in a small hash map (environment) to avoid big matrices
  covered <- new.env(parent = emptyenv())
  dxs <- seq.int(-s, s)
  dys <- seq.int(-s, s)

  # Loop over unique seed cells (fast in practice; avoids expanding per-point)
  for (i in seq_len(nrow(seeds))) {
    gx0 <- seeds$gx[i]; gy0 <- seeds$gy[i]
    for (dx in dxs) {
      gxk <- gx0 + dx
      for (dy in dys) {
        covered[[paste0(gxk, ":", gy0 + dy)]] <- TRUE
      }
    }
  }
  covered_count <- length(ls(covered, all.names = TRUE))
  area_grid <- covered_count * (cell_size ^ 2)

  return(ret(area_grid, paste0("grid_squares(cell_size=", cell_size, ", s=", s, ")"),
             list(cells_covered = covered_count,
                  note = "Grid squares approximate circles; expect upward bias that shrinks as 'cell_size' decreases.")))
}

## function to compute cells / mm2
calc_rates <- function(dts, areas, flag_combo, scale = 1) {
  vapply(names(dts), function(sid) {
    dt <- dts[[sid]]
    sum(fcheck(dt$cflag, 8) & fcheck(dt$pflag, flag_combo), na.rm = TRUE) /
      areas[[sid]] * scale
  }, numeric(1))
}

#####
#### load the colorectal cancer data
#####

# load the cohort sample data from Google Sheets
gopath <- "https://docs.google.com/spreadsheets/d/1VVzcuDDi8_n9NXWFuuUaUR0f9ghpAJO84Ccb31LX-gI/edit?gid=0#gid=0"
crr_crc <- as.data.table(googlesheets4::read_sheet(gopath, sheet="Clinical_data"))
stopifnot(all(!duplicated(crr_crc$Slide_ID)))

# reformat 
crr_crc[, TIL_num := c("absent"=0,"mild"=1,"moderate"=2,"marked"=3)[TIL]]
crr_crc[, mucinous := as.factor(ifelse(grepl("ucinous", Histology), 1, 0))]
crr_crc[, side := {
  if (Location %in% c("Appendix","Cecum","Ascending","Transverse")) { "Right" } else if (Location %in% c("Sigmoid","Descending")) { "Left" }  else { "Rectal" }
}, by="Slide_ID"]
crr_crc$Location <- factor(crr_crc$Location, levels = c("Appendix", "Cecum", "Ascending", "Transverse", "Descending", "Sigmoid", "Rectosigmoid", "Rectum"))
crr_crc[, PFSC := ifelse(crr_crc$PFSCensor==0,1,0)]
crr_crc[, pMMR := factor(ifelse(tipMMR %in% c("tipMMR","tdpMMR"), "pMMR", "dMMR"))]
#crr_crc$Grade <- factor(crr$Grade, levels = c("Low", "High"))
crr_crc[, TIL_num := c("Absent"=0,"Mild"=1,"Moderate"=2,"Marked"=3)[TIL]]
crr_crc[, mucinous := as.factor(ifelse(grepl("ucinous", Histology), 1, 0))]

## load the CRC quantified data
lf_crc = list.files("~/Dropbox_HMS/projects/orion/orion_1_74/csv_chain_coy2", full.names=TRUE, pattern = "csv$")

CRCCD8FLAG  <- 256
CRCCD4FLAG  <- 64
CRCCD20FLAG <- 1024
CRCCD3FLAG  <- 4096

load_crc <- function(x) {
  cat("loading", basename(x), "\n")
  dt <- fread(cmd = paste("cut -f1-6 -d','", x))
  attr(dt, "path") <- x
  attr(dt, "Slide_ID") <- sub(".*(LSP\\d+).*", "\\1", x)
  dt
}
dts_crc <- parallel::mclapply(lf_crc, load_crc, mc.cores = 8)
names(dts_crc) <- sub(".*(LSP\\d+).*", "\\1", basename(lf_crc))

## loop area calculation from ROIs
files <- list.files("/Users/jeremiahwala/Dropbox_HMS/projects/orion/orion_1_74/rois", pattern = "roi\\.csv$",
                    full.names = TRUE, ignore.case = TRUE)
areas_crc <- sapply(files, function(f)
  roi_area(f, keyword = "tumor", microns_per_pixel = 0.325))
names(areas_crc) <- sub(".*(LSP\\d+).*", "\\1", basename(files))

### CD3CD8
cd8_cd3_rates <- calc_rates(dts_crc, areas_crc, CRCCD8FLAG + CRCCD3FLAG)
cd8_cd3_dt <- data.table(
  Slide_ID       = sub(".*(LSP\\d+).*", "\\1", names(cd8_cd3_rates)),
  cd3cd8_rate   = unlist(cd8_cd3_rates)
)
abb_crc = merge(cd8_cd3_dt, crr_crc, by="Slide_ID")

### CD3CD4
cd4_cd3_rates <- calc_rates(dts_crc, areas_crc, CRCCD4FLAG + CRCCD3FLAG)
cd4_cd3_dt <- data.table(
  Slide_ID       = sub(".*(LSP\\d+).*", "\\1", names(cd4_cd3_rates)),
  cd3cd4_rate   = unlist(cd4_cd3_rates)
)
abb_crc = merge(cd4_cd3_dt, abb_crc, by="Slide_ID")

### CD20 rates
cd20_rates <- calc_rates(dts_crc, areas_crc, CRCCD20FLAG)
cd20_dt <- data.table(
  Slide_ID       = sub(".*(LSP\\d+).*", "\\1", names(cd20_rates)),
  cd20_rate   = unlist(cd20_rates)
)
abb_crc = merge(cd20_dt, abb_crc, by="Slide_ID")

## waterfall plot (for internal use only) of CRC
setkeyv(abb_crc, "cd20_rate") 
tumor_colors_Fig1 = c("tipMMR"="red", "dMMR"="#7670B2", "tdpMMR"="#169B77")
clevels <- abb_crc$Slide_ID
abb_crc[, Slide_ID := factor(Slide_ID, levels=clevels)]
abb_crc[, samplecount := seq(.N)]
g <- ggplot(abb_crc, aes(x=Slide_ID, y=cd20_rate, fill=tipMMR)) + 
  geom_bar(stat="identity", color="black") + ylab("Cells / mm2") + xlab("Sample ID") +
  coord_flip() +  scale_fill_manual(values=tumor_colors_Fig1) + theme_jw() + 
  theme(axis.line = element_line(), axis.ticks = element_line()) + 
  theme(legend.position = c(0.65, 0.3), 
        legend.justification = c(-0.2, -0.3), # Coordinates relative to the plot area
        legend.box.margin = margin(6, 6, 6, 6), # Adjust spacing around the legend box if necessary
        legend.margin = margin(-10, -10, -10, -10),
        text = element_text(family = "Helvetica"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(size = 9, margin = margin(r = 5)),
        axis.title.x = element_text(size = 9, margin = margin(t = 2)),
        plot.title = element_text(size = 9),
        legend.title = element_blank(),
        #legend.spacing.y = unit(1, "cm"),
        legend.text = element_text(size = 9),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(colour = "black"),  # Add axis lines
        panel.background = element_rect(fill = "transparent", colour = NA),  # remove panel background
        plot.background = element_rect(fill = "transparent", colour = NA)
  ) +
  scale_y_continuous(expand=c(0,0)) +
  geom_hline(yintercept = c(abb_crc[pMMR=="pMMR", median(cd20_rate)], 
                            abb_crc[pMMR=="dMMR", median(cd20_rate)]), color = "red", linetype = "dashed")

#########
## calculations for prostate
########

crr <- as.data.table(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1pXzLJkaX2wkJI_uOdgXLaapg1CJlo459G84FRY1pJYI/edit?gid=0#gid=0", sheet="manifest"))

CD8FLAG=4096
CD4FLAG=1024
CD3FLAG=2048
CD20FLAG=64

## load the data from ldata into a more friendly format
ab = lapply(ldata, function(x) {
  dt = as.data.table(x)[,.(cflag, pflag, x, y, sid)]
  dt
})

## total number of tumor cells
tot_tumor = sapply(ab, function(x) x[,sum(fcheck(cflag, 8))])
ab <- ab[tot_tumor > 0] ## limit to only those with tumor regions
ab <- ab[!names(ab) %in% "LSP12629"]  ## excluding LSP12629 as in other analyses

## loop area calculation from ROIs
files <- list.files("/Users/jeremiahwala/Dropbox_HMS/2025-Amiryousefi-Wala/jeremiah/primary_data/roisgu", pattern = "roi\\.csv$",
                    full.names = TRUE, ignore.case = TRUE)
areas_prostate <- sapply(files, function(f)
  roi_area(f, keyword = "tumor", microns_per_pixel = 0.650))
names(areas_prostate) <- sub(".*(LSP\\d+).*", "\\1", basename(files))

### CD3CD8
cd8_cd3_rates <- calc_rates(ab, areas_prostate, CD8FLAG + CD3FLAG)
cd8_cd3_dt <- data.table(
  Slide_ID       = sub(".*(LSP\\d+).*", "\\1", names(cd8_cd3_rates)),
  cd3cd8_rate   = unlist(cd8_cd3_rates)
)
abb_prostate = merge(cd8_cd3_dt, crr, by="Slide_ID")

## CD3CD4
cd4_cd3_rates <- calc_rates(ab, areas_prostate, CD4FLAG + CD3FLAG)
cd4_cd3_dt <- data.table(
  Slide_ID       = sub(".*(LSP\\d+).*", "\\1", names(cd4_cd3_rates)),
  cd3cd4_rate   = unlist(cd4_cd3_rates)
)
abb_prostate = merge(cd4_cd3_dt, abb_prostate, by="Slide_ID")

## CD20
cd20_rates <- calc_rates(ab, areas_prostate, CD20FLAG)
cd20_dt <- data.table(
  Slide_ID       = sub(".*(LSP\\d+).*", "\\1", names(cd20_rates)),
  cd20_rate   = unlist(cd20_rates)
)
abb_prostate = merge(cd20_dt, abb_prostate, by="Slide_ID")

## waterfall plot CD20 (Figure 1g)
setkeyv(abb_prostate, "cd20_rate") 
tumor_colors_Fig1 = c("High"="#7670B2", Low="#169B77")
clevels <- abb_prostate$Slide_ID
abb_prostate[, Slide_ID := factor(Slide_ID, levels=clevels)]
abb_prostate[, samplecount := seq(.N)]
g <- ggplot(abb_prostate[Slide_ID != "LSP12629"], aes(x=Slide_ID, y=cd20_rate, fill=Gleason)) + 
  geom_bar(stat="identity", color="black") + ylab("Cells / mm2") + xlab("Sample ID") +
  coord_flip() +  scale_fill_manual(values=tumor_colors_Fig1) + theme_jw() + 
  theme(axis.line = element_line(), axis.ticks = element_line()) + 
  theme(legend.position = c(0.65, 0.3), 
        legend.justification = c(-0.2, -0.3), # Coordinates relative to the plot area
        legend.box.margin = margin(6, 6, 6, 6), # Adjust spacing around the legend box if necessary
        legend.margin = margin(-10, -10, -10, -10),
        text = element_text(family = "Helvetica"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(size = 9, margin = margin(r = 5)),
        axis.title.x = element_text(size = 9, margin = margin(t = 2)),
        plot.title = element_text(size = 9),
        legend.title = element_blank(),
        #legend.spacing.y = unit(1, "cm"),
        legend.text = element_text(size = 9),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(colour = "black"),  # Add axis lines
        panel.background = element_rect(fill = "transparent", colour = NA),  # remove panel background
        plot.background = element_rect(fill = "transparent", colour = NA)
  ) +
  scale_y_continuous(expand=c(0,0)) +
  geom_hline(yintercept = c(abb_crc[pMMR=="pMMR", median(cd20_rate)], 
                            abb_crc[pMMR=="dMMR", median(cd20_rate)]), color = "red", linetype = "dashed")
save_ppdf(g, "Fig1g.pdf", width=3.5, height=3.5)
wilcox.test(abb_prostate[Gleason=="High", cd20_rate],abb_prostate[Gleason=="Low", cd20_rate])

wilcox.test(abb_prostate[Gleason=="High", cd20_rate],abb_crc[pMMR=="dMMR", cd20_rate])
wilcox.test(abb_prostate[Gleason=="Low", cd20_rate],abb_crc[pMMR=="pMMR", cd20_rate])

abb_prostate[, min(cd20_rate)]
abb_prostate[, max(cd20_rate)]
abb_prostate[Gleason=="High", mean(cd20_rate)]
abb_prostate[Gleason=="Low", mean(cd20_rate)]
abb_prostate[, max(cd20_rate)] / abb_prostate[, min(cd20_rate)]
abb_crc[pMMR=="pMMR", mean(cd20_rate)]
abb_crc[pMMR=="dMMR", mean(cd20_rate)]

abb_crc[, mean(cd20_rate)] / abb_prostate[, mean(cd20_rate)]

## waterfall plot CD3CD8 (Figure 1h)
setkeyv(abb_prostate, "cd3cd8_rate") 
clevels <- abb_prostate$Slide_ID
abb_prostate[, Slide_ID := factor(Slide_ID, levels=clevels)]
abb_prostate[, samplecount := seq(.N)]
g <- ggplot(abb_prostate, aes(x=Slide_ID, y=cd3cd8_rate, fill=Gleason)) + 
  geom_bar(stat="identity", color="black") + ylab("Cells / mm2") + xlab("Sample ID") +
  coord_flip() +  scale_fill_manual(values=tumor_colors_Fig1) + theme_jw() + 
  theme(axis.line = element_line(), axis.ticks = element_line()) + 
  theme(legend.position = c(0.65, 0.3), 
        legend.justification = c(-0.2, -0.3), # Coordinates relative to the plot area
        legend.box.margin = margin(6, 6, 6, 6), # Adjust spacing around the legend box if necessary
        legend.margin = margin(-10, -10, -10, -10),
        text = element_text(family = "Helvetica"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(size = 9, margin = margin(r = 5)),
        axis.title.x = element_text(size = 9, margin = margin(t = 2)),
        plot.title = element_text(size = 9),
        legend.title = element_blank(),
        #legend.spacing.y = unit(1, "cm"),
        legend.text = element_text(size = 9),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(colour = "black"),  # Add axis lines
        panel.background = element_rect(fill = "transparent", colour = NA),  # remove panel background
        plot.background = element_rect(fill = "transparent", colour = NA)
  ) +
  scale_y_continuous(expand=c(0,0)) + 
  geom_hline(yintercept = c(abb_crc[pMMR=="pMMR", median(cd3cd8_rate)], 
                            abb_crc[pMMR=="dMMR", median(cd3cd8_rate)]), color = "red", linetype = "dashed")
save_ppdf(g, "Fig1h.pdf", width=3.5, height=3.5)
wilcox.test(abb_prostate[Gleason=="High", cd3cd8_rate],abb_prostate[Gleason=="Low", cd3cd8_rate])

abb_prostate[, min(cd3cd8_rate)]
abb_prostate[, max(cd3cd8_rate)]
abb_prostate[Gleason=="High", mean(cd3cd8_rate)]
abb_prostate[Gleason=="Low", mean(cd3cd8_rate)]
abb_prostate[, max(cd3cd8_rate)] / abb_prostate[, min(cd3cd8_rate)]
abb_crc[pMMR=="pMMR", mean(cd3cd8_rate)]
abb_crc[pMMR=="dMMR", mean(cd3cd8_rate)]

abb_crc[, mean(cd3cd8_rate)] / abb_prostate[, mean(cd3cd8_rate)]


## waterfall plot CD3CD4 (Figure 1i)
setkeyv(abb_prostate, "cd3cd4_rate") 
clevels <- abb_prostate$Slide_ID
abb_prostate[, Slide_ID := factor(Slide_ID, levels=clevels)]
abb_prostate[, samplecount := seq(.N)]
g <- ggplot(abb_prostate, aes(x=Slide_ID, y=cd3cd4_rate, fill=Gleason)) + 
  geom_bar(stat="identity", color="black") + ylab("Cells / mm2") + xlab("Sample ID") +
  coord_flip() +  scale_fill_manual(values=tumor_colors_Fig1) + theme_jw() + 
  theme(axis.line = element_line(), axis.ticks = element_line()) + 
  theme(legend.position = c(0.65, 0.3), 
        legend.justification = c(-0.2, -0.3), # Coordinates relative to the plot area
        legend.box.margin = margin(6, 6, 6, 6), # Adjust spacing around the legend box if necessary
        legend.margin = margin(-10, -10, -10, -10),
        text = element_text(family = "Helvetica"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(size = 9, margin = margin(r = 5)),
        axis.title.x = element_text(size = 9, margin = margin(t = 2)),
        plot.title = element_text(size = 9),
        legend.title = element_blank(),
        #legend.spacing.y = unit(1, "cm"),
        legend.text = element_text(size = 9),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(colour = "black"),  # Add axis lines
        panel.background = element_rect(fill = "transparent", colour = NA),  # remove panel background
        plot.background = element_rect(fill = "transparent", colour = NA)
  ) +
  scale_y_continuous(expand=c(0,0)) + 
  geom_hline(yintercept = c(abb_crc[pMMR=="pMMR", median(cd3cd4_rate)], 
                            abb_crc[pMMR=="dMMR", median(cd3cd4_rate)]), color = "red", linetype = "dashed")
save_ppdf(g, "Fig1i.pdf", width=3.5, height=3.5)
wilcox.test(abb_prostate[Gleason=="High", cd3cd4_rate],abb_prostate[Gleason=="Low", cd3cd4_rate])

abb_prostate[, min(cd3cd4_rate)]
abb_prostate[, max(cd3cd4_rate)]
abb_prostate[Gleason=="High", mean(cd3cd4_rate)]
abb_prostate[Gleason=="Low", mean(cd3cd4_rate)]
abb_prostate[, max(cd3cd4_rate)] / abb_prostate[, min(cd3cd4_rate)]
abb_crc[pMMR=="pMMR", mean(cd3cd4_rate)]
abb_crc[pMMR=="dMMR", mean(cd3cd4_rate)]

abb_crc[, mean(cd3cd4_rate)] / abb_prostate[, mean(cd3cd4_rate)]


#####
##### FIGURE 1F
#####
library(clinfun)

# >>> set these from bit_report() output <
grade_bits  <- c(G1 = 4096, G2 = 8192, G3 = 16384, G4 = 32768, G5 = 65536)
CD8FLAG <- 4096; CD4FLAG <- 1024; CD3FLAG <- 2048; CD20FLAG <- 64
# <<< 

groups      <- names(grade_bits)
CYCIF_SCALE <- 2
RATE_FACTOR <- 1e6

dd <- dts_prostate[!names(dts_prostate) %in% "LSP12629"]   # remove outlier
#dd <- dts_prostate[names(dts_prostate) %in% "LSP12657"]   # remove outlier

grade_n <- sapply(grade_bits, function(b)
  sum(vapply(dd, function(dt) sum(bitops::bitAnd(dt$cflag, b) > 0), numeric(1))))
print(grade_n)
if (all(grade_n == 0))
  stop("grade_bits don't match dts_prostate$cflag - set them from bit_report() output.")

subregion <- rbindlist(Map(function(dt, nm) {
  cf    <- dt$cflag
  grade <- fcase(
    bitops::bitAnd(cf, grade_bits["G1"]) > 0, "G1",
    bitops::bitAnd(cf, grade_bits["G2"]) > 0, "G2",
    bitops::bitAnd(cf, grade_bits["G3"]) > 0, "G3",
    bitops::bitAnd(cf, grade_bits["G4"]) > 0, "G4",
    bitops::bitAnd(cf, grade_bits["G5"]) > 0, "G5",
    default = "stroma"
  )
  id <- sub(".*(LSP\\d+).*", "\\1", nm)
  rbindlist(lapply(groups, function(GGG) {
    ix <- grade == GGG
    if (!any(ix)) return(NULL)
    xy   <- data.table(x = dt$x[ix] * CYCIF_SCALE, y = dt$y[ix] * CYCIF_SCALE)
    area <- cycif_disc_union_area(xy)$area
    data.table(
      Slide_ID  = id, region = GGG,
      cd8_rate  = sum(ix & fcheck(dt$pflag, CD8FLAG + CD3FLAG), na.rm = TRUE) / area * RATE_FACTOR,
      cd4_rate  = sum(ix & fcheck(dt$pflag, CD4FLAG + CD3FLAG), na.rm = TRUE) / area * RATE_FACTOR,
      cd20_rate = sum(ix & fcheck(dt$pflag, CD20FLAG),          na.rm = TRUE) / area * RATE_FACTOR
    )
  }))
}, dd, names(dd)), fill = TRUE)

subregion[, region := factor(region, levels = groups,
                             labels = c("3+3","3+4","4+3","4+4","Any 5"))]

## plot it
long <- melt(subregion, id.vars = c("Slide_ID","region"),
             measure.vars = c("cd8_rate","cd4_rate","cd20_rate"),
             variable.name = "marker", value.name = "rate")
long[, marker := factor(marker, levels = c("cd8_rate","cd4_rate","cd20_rate"),
                        labels = c("CD8","CD4","CD20"))]

# Jonckheere-Terpstra trend across ordered grades, per marker + overall
jt <- long[, .(p = jonckheere.test(rate, as.numeric(region),
                                   alternative = "increasing")$p.value), by = marker]
tot <- subregion[, .(rate = cd8_rate + cd4_rate + cd20_rate), by = .(Slide_ID, region)]
jt_total <- jonckheere.test(tot$rate, as.numeric(tot$region), alternative = "increasing")$p.value

means <- long[, .(mean_rate = mean(rate, na.rm = TRUE)), by = .(region, marker)]
means[, marker := factor(marker, levels = c("CD8","CD4","CD20"))]  # CD8 top, CD20 bottom
cols  <- c(CD8 = "#4DAF4A", CD4 = "#377EB8", CD20 = "#E41A1C")
labs  <- setNames(sprintf("%s  P = %.3f", jt$marker, jt$p), jt$marker)

g <- ggplot(means, aes(region, mean_rate, fill = marker)) +
  geom_col(width = 0.85, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = cols, breaks = c("CD8","CD4","CD20"),
                    labels = labs[c("CD8","CD4","CD20")], name = "Marker") +
  labs(x = NULL, y = expression("Cells / mm"^2)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  annotate("segment", x = 1.5, xend = 4.5,
           y = max(means$mean_rate)*0.4, yend = max(means$mean_rate)*0.9,
           arrow = arrow(length = unit(0.15,"cm"))) +
  annotate("text", x = 2.2, y = max(means$mean_rate)*0.75,
           label = sprintf("J-T test P = %.3f", jt_total), angle = 30, size = 3) +
  theme_jw() +
  theme(legend.position = c(0.04, 0.96), legend.justification = c(0,1),
        legend.title = element_blank(), legend.text = element_text(size = 8))
save_ppdf(g, "Fig1f.pdf", width = 3, height = 3)

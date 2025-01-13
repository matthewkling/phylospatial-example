
# analysis ====================================================================


# install the package from CRAN, load it, and view vignettes
# install.packages("phylospatial")
library(phylospatial)
# browseVignettes("phylospatial")


# load data and build phylospatial data set
sdm <- terra::rast("data/sdm.tif")
tree <- ape::read.tree("data/chronogram.tree")
d <- phylospatial(sdm, tree)


# calculate alpha diversity, significance, and CANAPE
div <- ps_diversity(d)
sig <- ps_rand(d, n_rand = 999, n_cores = 8)
cnp <- ps_canape(sig)


# calculate pairwise dissimilarity; use for ordination & clustering
d <- ps_add_dissim(d, method = "sorensen")
rgb <- ps_rgb(d, trans = rank)
reg <- ps_regions(d, k = 8, method = "kmeans")


# conservation prioritization
con <- ps_prioritize(d,
                     init = terra::rast("data/protection.tif"),
                     cost = terra::rast("data/popdens.tif"))



# figure ======================================================================

# load packages
library(terra)
library(tidyverse)
library(patchwork)

# run analyses
source("code/phylospatial.R")

# define custom theme
thm <- theme_void() +
      theme(legend.position = "inside",
            legend.position.inside = c(.48, .94),
            legend.justification = c(0, 1),
            legend.title = element_text(face = "bold"))

# PD
pdiv <- div %>%
      as.data.frame(xy = T) %>%
      ggplot(aes(x, y, fill = PD)) +
      geom_raster() +
      scale_fill_gradientn(colors = c("gray90", "limegreen", "darkgreen", "black")) +
      labs(fill = "(a) PD") +
      thm

# PD significance
psig <- sig %>%
      as.data.frame(xy = T) %>%
      mutate(qpd = pmax(1/1000, pmin(999/1000, qPD))) %>%
      ggplot(aes(x, y, fill = qpd)) +
      geom_raster() +
      scale_fill_gradientn(colors = c("darkorchid4", "orchid3", "gray90", "limegreen", "darkgreen"),
                           trans = "logit", breaks = c(.001, .05, .5, .95, .999)) +
      labs(fill = "(b) PD null quantile") +
      thm

# CANAPE
pcnp <- cnp %>%
      as.data.frame(xy = T) %>%
      ggplot(aes(x, y, fill = endemism)) +
      geom_raster() +
      scale_fill_manual(values = cats(cnp)[[1]]$color) +
      labs(fill = "(c) Endemism category") +
      thm

# RGB ordination
prgb <- rgb %>%
      as.data.frame(xy = T) %>%
      mutate(color = rgb(r, g, b)) %>%
      ggplot(aes(x, y)) +
      geom_tile(aes(color = r)) +
      geom_raster(aes(fill = color)) +
      scale_fill_identity() +
      labs(color = "(d) RGB ordination") +
      thm +
      guides(color = guide_colorbar(barwidth = 0)) +
      theme(legend.text = element_blank(),
            legend.key = element_blank())

# phylogenetic regions
preg <- reg %>%
      as.data.frame(xy = T) %>%
      ggplot(aes(x, y, fill = factor(phyloregion))) +
      geom_raster() +
      scale_fill_brewer(type = "qual", palette = 3) +
      guides(fill = guide_legend(ncol = 2)) +
      labs(fill = "(e) Phylogenetic region") +
      thm

# conservation priorities
pcon <- con %>%
      as.data.frame(xy = T) %>%
      ggplot(aes(x, y, fill = priority)) +
      geom_raster() +
      scale_fill_gradientn(colors = c(viridis::viridis_pal()(5)[c(1:5, 5)]),
                           trans = "log10") +
      guides(fill = guide_colorbar(reverse = T)) +
      labs(fill = "(f) Conservation priority") +
      thm

# combine and save
p <- pdiv + psig + pcnp + prgb + preg + pcon +
      plot_layout(nrow = 2)
ggsave("figures/maps.png", p, width = 12, height = 10, units = "in")

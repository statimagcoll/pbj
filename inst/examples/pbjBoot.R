# loading example data
set.seed(11)
pain <- pain21::pain21()
pain$data$group <- factor(sample(1:4, size = nrow(pain$data), replace=TRUE))
pain$data$x <- rnorm(nrow(pain$data))
# pain$data$Winv <- runif(nrow(pain$data))
pain$data$W <- 1/runif(nrow(pain$data))
# creates a fake ID variable
pain$data$ID = c(rep(1:10, each=2), 11)
imgs <- simplify2array(RNifti::readNifti(pain$data$images))
Winvs <- simplify2array(RNifti::readNifti(pain$data$varimages))
mask <- RNifti::readNifti(pain$mask) * c(apply(imgs!=0, 1:3, all))

# get one voxel
testvox <- which(mask==1, arr.ind = TRUE)[1:2,]
mask[,,] <- 0
mask[testvox] <- 1

pain$data$y <- imgs[testvox[1,1], testvox[1,2], testvox[1,3], ]
pain$data$Winv.img <- Winvs[testvox[1,1], testvox[1,2], testvox[1,3], ]

outdir <- tempfile("dir")
dir.create(outdir)
statMap <- lmPBJ(pain$data$images, form = ~ 1, formred = NULL, mask = mask,
                 id=pain$data$ID, template=pain$template, data = pain$data,
                 W = pain$data$W, robust = TRUE, zeros=TRUE, transform='none',
                 HC3 = TRUE, W_structure = "exchangeable", outdir = outdir)

sqrtSigma <- readRDS(statMap$sqrtSigma)
mask <- statMap$mask
dims <- dim(sqrtSigma$res)
rboot=function(n){ (2*stats::rbinom(n, size=1, prob=0.5)-1)}

# under the null
actual <- pbjBoot(sqrtSigma, method = 'wild')

# under the alternative
actual <- pbjBoot(sqrtSigma, method = 'wild', null = FALSE)

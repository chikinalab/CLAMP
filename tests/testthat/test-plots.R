test_that("CLAMPdotplot can use AUC or -log10(FDR) on the x-axis", {
    summ <- data.frame(
        pathway = c("Path A", "Path B", "Path C"),
        LV = rep("LV2", 3),
        AUC = c(0.99, 0.70, 0.80),
        FDR = c(0.05, 0.001, 0.01),
        stringsAsFactors = FALSE
    )

    p_auc <- CLAMPdotplot(
        list(summary = summ),
        lv = "LV2",
        top = 2,
        auc.cutoff = 0.6,
        fdr.cutoff = 0.1,
        x.axis = "AUC",
        order.by = "AUC"
    )

    expect_s3_class(p_auc, "ggplot")
    expect_equal(p_auc$theme$plot.title$hjust, 0.5)
    expect_equal(p_auc$data$x_value, p_auc$data$AUC)

    p_fdr <- CLAMPdotplot(
        list(summary = summ),
        lv = "LV2",
        top = 2,
        auc.cutoff = 0.6,
        fdr.cutoff = 0.1,
        x.axis = "-log10(FDR)",
        order.by = "-log10(FDR)"
    )

    expect_s3_class(p_fdr, "ggplot")
    expect_equal(p_fdr$data$x_value, p_fdr$data$log10_fdr)
    expect_equal(p_fdr$labels$x, "-log10(FDR)")
    expect_setequal(p_fdr$data$pathway_label, c("Path B", "Path C"))
})

test_that("CLAMPplotTopZ returns a loading plot for selected LVs", {
    set.seed(1)
    Z <- matrix(
        runif(20),
        nrow = 10,
        dimnames = list(paste0("Gene", seq_len(10)), c("LV1", "LV2"))
    )
    U <- matrix(1, nrow = 2, ncol = 2)
    clampRes <- list(Z = Z, U = U)

    p <- CLAMPplotTopZ(clampRes, top = 5, index = "LV2", label.top = 3)

    expect_s3_class(p, "ggplot")
    expect_equal(p$labels$title, "LV2")
    expect_equal(p$theme$plot.title$hjust, 0.5)
})

test_that("CLAMPplotU handles missing FDR values", {
    pathways <- c("Path1", "Path2")
    lvs <- "LV1"
    clampRes <- list(
        U = matrix(1, 2, 1, dimnames = list(pathways, lvs)),
        Uauc = matrix(0.5, 2, 1, dimnames = list(pathways, lvs)),
        Up = matrix(1, 2, 1, dimnames = list(pathways, lvs)),
        summary = data.frame(
            pathway = pathways,
            LV = rep(lvs, 2),
            AUC = 0.5,
            p_value = c(0.5, 0.6),
            FDR = c(NA, 1)
        )
    )

    expect_message(
        result <- CLAMPplotU(clampRes),
        "No entries pass the AUC/FDR thresholds"
    )
    expect_null(result)
})

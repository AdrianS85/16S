

corncob::bbdml()

corncob <- bbdml(formula = OTU.1 ~ 1,
                 phi.formula = ~ 1,
                 data = soil)

corncob_da <- bbdml(formula = OTU.1 ~ DayAmdmt,
                    phi.formula = ~ DayAmdmt,
                    data = soil)

lrtest(mod_null = corncob, mod = corncob_da)


da_analysis <- differentialTest(formula = ~ DayAmdmt,
                                phi.formula = ~ DayAmdmt,
                                formula_null = ~ 1,
                                phi.formula_null = ~ DayAmdmt,
                                test = "Wald", boot = FALSE,
                                data = soil,
                                fdr_cutoff = 0.05)


otu_to_taxonomy(OTU = da_analysis$significant_taxa, data = soil)